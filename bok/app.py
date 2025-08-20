

# File: activity_cliff_app/app.py
import streamlit as st
import pandas as pd
from qdrant_client import QdrantClient
import os
from pathlib import Path
from cliff_detector import (
    detect_activity_cliffs,
    detect_activity_cliffs_sali,
    detect_activity_cliffs_knn,
)
from utils import compute_pIC50, smiles_diff_to_images, highlight_canonical_smiles_diff
from prompts import generate_prompt, call_llm, generate_fewshot_prompt, generate_rag_prompt, translate_text # Import translate_text
from structure_diff import detect_diff_type
from rdkit import Chem # Import Chem for SMILES validation
from datetime import date # Import date for report generation date

st.set_page_config(page_title="Activity Cliff Analyzer", layout="wide")

# Use @st.cache_resource with a cleanup function (ttl is not strictly needed but good practice)
@st.cache_resource(ttl=3600)
def get_qdrant_client():
    """Initializes and caches the Qdrant client."""
    base_dir = Path(__file__).resolve().parent
    qdrant_path = base_dir / "rag" / "qdrant"
    return QdrantClient(path=str(qdrant_path))

st.title("🔬 Activity Cliff 자동 탐지 & 해석")

# Clear Cache Button
if st.sidebar.button("캐시 지우기"): # Clear Cache button
    st.cache_data.clear()
    st.cache_resource.clear()
    st.sidebar.success("캐시가 지워졌습니다!")

uploaded_file = st.file_uploader("데이터 업로드 (SMILES, IC50 포함)", type=["csv", "xlsx"])
if uploaded_file:
    # Determine file type and read accordingly
    if uploaded_file.name.endswith('.csv'):
        df = pd.read_csv(uploaded_file)
    elif uploaded_file.name.endswith('.xlsx'):
        df = pd.read_excel(uploaded_file)
    else:
        st.error("지원되지 않는 파일 형식입니다. CSV 또는 XLSX 파일을 업로드해주세요.")
        st.stop()

    st.write("업로드된 데이터 앞부분:", df.head(2))

    # Validate SMILES entries
    invalid_smiles = []
    if 'SMILES' in df.columns:
        for idx, smiles in df['SMILES'].items():
            if Chem.MolFromSmiles(str(smiles)) is None:
                invalid_smiles.append((idx, smiles))
    else:
        st.error("업로드된 파일에 'SMILES' 컬럼이 없습니다. 'SMILES' 컬럼이 포함된 파일을 업로드해주세요.")
        st.stop()

    if invalid_smiles:
        st.error("다음 SMILES 항목이 유효하지 않습니다. 파일을 수정하여 다시 업로드해주세요:")
        for idx, smiles in invalid_smiles:
            st.write(f"- 행 {idx}: {smiles}")
        st.stop() # Stop execution if invalid SMILES are found

    # Ensure pIC50 is present or computable
    if 'pIC50' not in df.columns:
        if 'IC50' in df.columns:
            try:
                df['pIC50'] = compute_pIC50(df['IC50'])
            except Exception:
                st.error("IC50 값을 pIC50로 변환하는 중 오류가 발생했습니다. 숫자 형식과 단위를 확인해주세요 (nM).")
                st.stop()
        else:
            st.error("업로드된 파일에 'pIC50' 또는 'IC50' 컬럼이 필요합니다. 둘 중 하나를 포함해 주세요.")
            st.stop()

    st.sidebar.header("🔧 탐지 설정")
    method = st.sidebar.selectbox(
        "탐지 방법 선택",
        ("기본 임계값", "SALI 순위", "k-NN 이웃 기반"),
        index=0,
    )

    # Common / method-specific parameters
    if method == "기본 임계값":
        sim_thres = st.sidebar.slider("구조 유사도 임계값", 0.7, 1.0, 0.85, 0.01)
        act_thres = st.sidebar.slider("pIC50 차이 임계값", 0.5, 3.0, 1.0, 0.1)
    elif method == "SALI 순위":
        top_n = st.sidebar.number_input("상위 N SALI 쌍", min_value=10, max_value=2000, value=200, step=10)
        sim_floor = st.sidebar.slider("최소 유사도 (SALI 필터)", 0.0, 1.0, 0.5, 0.01)
        radius = st.sidebar.select_slider("Morgan 반경", options=[2, 3], value=2)
        scaffold_constrained = st.sidebar.checkbox("동일 스캐폴드 내에서만 비교", value=False)
    else:  # k-NN
        k = st.sidebar.slider("이웃 수 k", 1, 50, 8)
        act_thres = st.sidebar.slider("pIC50 차이 임계값", 0.5, 3.0, 1.0, 0.1)
        min_sim = st.sidebar.slider("최소 유사도 (이웃 필터)", 0.0, 1.0, 0.0, 0.01)
        radius = st.sidebar.select_slider("Morgan 반경", options=[2, 3], value=2)
        scaffold_constrained = st.sidebar.checkbox("동일 스캐폴드 내에서만 비교", value=False)

    # Replaced checkboxes with radio button for prompt generation strategy
    prompt_strategy = st.sidebar.radio(
        "LLM 해석 요청 방식 선택",
        ("기본", "Few-shot 예시 기반", "RAG 기반 도메인 해석"),
        index=0 # Default to Basic
    )

    # LLM Output Language Selection
    llm_output_lang = st.sidebar.radio(
        "LLM 분석 결과 언어 선택",
        ("한국어", "English"),
        index=0 # Default to Korean
    )

    # Get the cached Qdrant client
    qdrant_client = get_qdrant_client()

    with st.spinner("Activity Cliff 탐지 중..."):
        if method == "기본 임계값":
            results = detect_activity_cliffs(df, sim_thres, act_thres)
        elif method == "SALI 순위":
            results = detect_activity_cliffs_sali(
                df,
                top_n=top_n,
                radius=radius,
                fp_bits=2048,
                sim_floor=sim_floor,
                scaffold_constrained=scaffold_constrained,
            )
        else:  # k-NN
            results = detect_activity_cliffs_knn(
                df,
                k=k,
                act_thres=act_thres,
                radius=radius,
                fp_bits=2048,
                min_sim=min_sim,
                scaffold_constrained=scaffold_constrained,
            )

    st.success(f"{len(results)}개의 Activity Cliff 쌍이 탐지되었습니다.")

    # Download Results Button
    if not results.empty:
        csv_data = results.to_csv(index=False).encode('utf-8')
        st.download_button(
            label="결과 다운로드 (CSV)",
            data=csv_data,
            file_name="activity_cliffs.csv",
            mime="text/csv",
        )

    for i, row in results.iterrows():
        col1, col2 = st.columns(2)

        # Generate images and highlighted SMILES text
        img1, img2 = smiles_diff_to_images(row['mol1_smiles'], row['mol2_smiles'])

        with col1:
            st.image(img1, caption=f"A: {row['mol1_activity']:.2f}")
            st.markdown(f"SMILES: `{row['mol1_smiles']}`")
            if 'mol1_canonical' in results.columns:
                if 'mol2_canonical' in results.columns:
                    a_html, _ = highlight_canonical_smiles_diff(row['mol1_canonical'], row['mol2_canonical'])
                    st.markdown("Canonical (A):", help="Red background indicates segments that differ from B")
                    st.markdown(a_html, unsafe_allow_html=True)
                else:
                    st.markdown(f"Canonical: `{row['mol1_canonical']}`")
        with col2:
            st.image(img2, caption=f"B: {row['mol2_activity']:.2f}")
            st.markdown(f"SMILES: `{row['mol2_smiles']}`")
            if 'mol2_canonical' in results.columns:
                if 'mol1_canonical' in results.columns:
                    _, b_html = highlight_canonical_smiles_diff(row['mol1_canonical'], row['mol2_canonical'])
                    st.markdown("Canonical (B):", help="Blue background indicates segments that differ from A")
                    st.markdown(b_html, unsafe_allow_html=True)
                else:
                    st.markdown(f"Canonical: `{row['mol2_canonical']}`")

        # Metrics line (include SALI if present)
        if 'sali' in results.columns:
            st.markdown(f"- 유사도: `{row['sim']:.2f}` / pIC50 차이: `{row['activity_diff']:.2f}` / SALI: `{row['sali']:.2f}`")
        else:
            st.markdown(f"- 유사도: `{row['sim']:.2f}` / pIC50 차이: `{row['activity_diff']:.2f}`")

        # Show structural difference type for transparency
        try:
            _mol_a = Chem.MolFromSmiles(row['mol1_smiles'])
            _mol_b = Chem.MolFromSmiles(row['mol2_smiles'])
            diff_type = detect_diff_type(_mol_a, _mol_b)
            st.markdown(f"- 구조 차이 유형: `{diff_type}`")
        except Exception:
            pass

        # Added a button to trigger prompt generation and LLM call
        if st.button(f"LLM 해석 생성 - 쌍 {i}", key=f"llm_gen_btn_{i}"):
            system_prompt = ""
            user_prompt = ""
            if prompt_strategy == "RAG 기반 도메인 해석":
                system_prompt, user_prompt = generate_rag_prompt(row, qdrant_client)
            elif prompt_strategy == "Few-shot 예시 기반":
                system_prompt, user_prompt = generate_fewshot_prompt(row)
            else:
                system_prompt, user_prompt = generate_prompt(row)

            # Call LLM with English prompts
            llm_output_raw = call_llm(system_prompt, user_prompt)
            
            # Translate if Korean is selected, otherwise use English
            if llm_output_lang == "한국어":
                llm_output_final = translate_text(llm_output_raw, target_language="Korean")
            else:
                llm_output_final = llm_output_raw

            # Fix line formatting for Markdown display
            llm_output_formatted = llm_output_final.replace("\n", "\n\n")

            # Construct the full report based on the example format
            # Display the two molecules and their info side by side in a single block
            st.markdown("""
< 자동 생성 SAR 요약 리포트 >

**분석 대상:** Janus kinase (JAK) 저해제 후보 화합물 512종
                        
**리포트 생성일:** {report_date}

**핵심 분석 1: 주요 활성 변화 요인 (Key Activity Cliff)**

**요약:** 분자의 특정 3차원 구조가 활성에 {activity_ratio:.0f}배 차이를 유발함.
""".format(
                report_date=date.today().strftime("%Y-%m-%d"),
                activity_ratio=max(row['mol1_activity'], row['mol2_activity']) / min(row['mol1_activity'], row['mol2_activity']),
            ))

            col_a, col_b = st.columns(2)
            with col_a:
                st.markdown(f"""
**화합물 ID:** Molecule {row['mol1_idx']} (활성: {row['mol1_activity']:.2f} nM)

**구조:**
""")
                st.image(img1, caption=f"Molecule {row['mol1_idx']}")
            with col_b:
                st.markdown(f"""
**화합물 ID:** Molecule {row['mol2_idx']} (활성: {row['mol2_activity']:.2f} nM)

**구조:**
""")
                st.image(img2, caption=f"Molecule {row['mol2_idx']}")


            st.markdown(f"""
**자동화된 해석 및 가설:**
{llm_output_formatted}
""")

        st.divider() # Add a separator after each pair

import streamlit as st
import pandas as pd
import datetime
import os
import json
import zipfile
import io
from modules.cheminformatics import find_activity_cliffs
from modules.visualization import visualize_structure_difference
from modules.llm_handler import generate_hypothesis
from modules.io_utils import load_smiles_activity_csv, save_hypothesis_to_md

# --- Helper Functions ---

def get_openai_api_key_from_file():
    try:
        with open("openAI_key.txt", "r") as f:
            return f.read().strip()
    except FileNotFoundError:
        return None

def format_hypothesis_for_markdown(data: dict) -> str:
    """주어진 가설 데이터(dict)를 가독성 좋은 마크다운 문자열로 변환합니다."""
    md_lines = []

    # 기본 정보
    md_lines.append(f"### 🏆 주요 가설: {data.get('primary_hypothesis', 'N/A')}")
    md_lines.append(f"- **더 활성이 높은 화합물:** `{data.get('more_active', 'N/A')}`")
    md_lines.append(f"- **활성도 변화 설명:** {data.get('delta_pAct_explained', 'N/A')}")
    md_lines.append(f"- **신뢰도:** {data.get('confidence', 0.0) * 100:.1f}%")
    md_lines.append("\n")

    # 기전 분석
    md_lines.append("### 🔬 기전 분석")
    rationale = data.get('mechanistic_rationale', {})
    for key, value in rationale.items():
        if value and value not in ["N/A", "선택"]:
            md_lines.append(f"- **{key.replace('_', ' ').title()}:** {value}")
    md_lines.append("\n")

    # 설계 제안
    md_lines.append("### 💡 검증을 위한 분자 설계 제안")
    suggestions = data.get('design_suggestions', [])
    if suggestions:
        df = pd.DataFrame(suggestions)
        md_lines.append(df.to_markdown(index=False))
    md_lines.append("\n")

    # 반대 가설
    md_lines.append("### 🤔 반대 가설")
    counter_hypotheses = data.get('counter_hypotheses', [])
    for i, counter in enumerate(counter_hypotheses, 1):
        md_lines.append(f"{i}. {counter}")
    md_lines.append("\n")

    # ADMET 위험성
    md_lines.append("### ⚠️ ADMET 위험성 예측")
    admet_flags = data.get('admet_flags', [])
    for flag in admet_flags:
        md_lines.append(f"- {flag}")
    md_lines.append("\n")
    
    # 가정 및 한계
    md_lines.append("### 📋 가정 및 한계")
    assumptions = data.get('assumptions_and_limits', [])
    for assumption in assumptions:
        md_lines.append(f"- {assumption}")

    return "\n".join(md_lines)


# --- Streamlit App ---

st.set_page_config(layout="wide")
st.title("🔬 SAR 분석 및 가설 생성 자동화 도구")
st.write("분자 구조와 활성 데이터 기반의 구조-활성 관계(SAR) 분석 및 가설 생성을 자동화합니다.")

# --- 1. 데이터 업로드 ---
st.header("1. 데이터 업로드")

uploaded_file = st.file_uploader("분자 구조(SMILES)와 활성 데이터가 포함된 CSV 파일을 업로드하세요.", type="csv")

if uploaded_file is not None:
    try:
        df, suggestion = load_smiles_activity_csv(uploaded_file)
        st.success("파일이 성공적으로 업로드되었습니다!")
    except Exception as e:
        st.error(f"CSV 로딩 중 오류가 발생했습니다: {e}")
        df = None

    if df is not None:
        st.subheader("업로드된 데이터 미리보기")
        st.dataframe(df.head())
        st.caption("자동 인식된 컬럼 제안값을 확인하세요. 필요 시 변경 가능합니다.")
        st.session_state['auto_suggestion'] = suggestion

    st.session_state['df'] = df

# --- 2. Activity Cliff 분석 설정 ---
st.header("2. Activity Cliff 분석")

if 'df' in st.session_state and st.session_state['df'] is not None:
    df = st.session_state['df']
    
    col1, col2 = st.columns(2)
    with col1:
        auto = st.session_state.get('auto_suggestion', {})
        smiles_col_default = auto.get('smiles_col') if auto.get('smiles_col') in df.columns else None
        activity_col_default = auto.get('activity_col') if auto.get('activity_col') in df.columns else None

        smiles_col = st.selectbox("SMILES 컬럼을 선택하세요:", df.columns, index=(list(df.columns).index(smiles_col_default) if smiles_col_default else 0))
        activity_col = st.selectbox("활성도 컬럼을 선택하세요:", df.columns, index=(list(df.columns).index(activity_col_default) if activity_col_default else (1 if len(df.columns) > 1 else 0)))

    with col2:
        similarity_threshold = st.slider("구조 유사도 임계값 (Tanimoto)", 0.7, 1.0, 0.85, 0.01)
        activity_diff_threshold = st.number_input("활성도 차이 임계값", min_value=0.0, value=1.0, step=0.1)

    if st.button("Activity Cliff 분석 실행"):
        with st.spinner("Activity Cliff를 분석 중입니다..."):
            work_df = df.copy()
            work_df = work_df.dropna(subset=[activity_col])
            work_df = work_df.reset_index(drop=True)
            
            cliff_df = find_activity_cliffs(
                work_df,
                smiles_col=smiles_col,
                activity_col=activity_col,
                similarity_threshold=similarity_threshold,
                activity_diff_threshold=activity_diff_threshold
            )
        
        st.success(f"{len(cliff_df)}개의 Activity Cliff 쌍을 찾았습니다!")
        st.dataframe(cliff_df)
        st.session_state['cliff_df'] = cliff_df

# --- 3. 결과 시각화 및 가설 생성 ---
st.header("3. 결과 시각화 및 가설 생성")
if 'cliff_df' in st.session_state and not st.session_state['cliff_df'].empty:
    cliff_df = st.session_state['cliff_df']
    
    selected_indices = st.multiselect("분석 및 시각화할 Activity Cliff 쌍을 선택하세요:", cliff_df.index)
    
    if selected_indices:
        openai_api_key = get_openai_api_key_from_file()
        if not openai_api_key:
            st.warning("LLM 가설 생성을 위해 openAI_key.txt 파일에 API Key를 입력해주세요.")
        else:
            st.session_state['openai_api_key'] = openai_api_key
            st.success("API 키가 openAI_key.txt 파일에서 로드되었습니다.")

            if st.button("선택된 쌍에 대한 가설 생성"):
                output_dir = "hypotheses"
                os.makedirs(output_dir, exist_ok=True)

                for i in selected_indices:
                    row = cliff_df.loc[i]
                    st.subheader(f"분석 쌍 #{i}")
                    
                    # 시각화
                    img = visualize_structure_difference(
                        smiles1=row['SMILES_1'],
                        smiles2=row['SMILES_2'],
                        legend1=f"SMILES: {row['SMILES_1']}\nActivity: {row['Activity_1']:.2f}",
                        legend2=f"SMILES: {row['SMILES_2']}\nActivity: {row['Activity_2']:.2f}"
                    )
                    st.image(img, caption=f"유사도: {row['Similarity']:.3f} | 활성도 차이: {row['Activity_Diff']:.2f}")

                    # 가설 생성
                    with st.spinner(f"쌍 #{i}에 대한 LLM 가설을 생성 중입니다..."):
                        if row['Activity_1'] > row['Activity_2']:
                            high_act_smiles, high_act_val = row['SMILES_1'], row['Activity_1']
                            low_act_smiles, low_act_val = row['SMILES_2'], row['Activity_2']
                        else:
                            high_act_smiles, high_act_val = row['SMILES_2'], row['Activity_2']
                            low_act_smiles, low_act_val = row['SMILES_1'], row['Activity_1']

                        json_response = generate_hypothesis(
                            api_key=openai_api_key,
                            smiles1=low_act_smiles,
                            activity1=low_act_val,
                            smiles2=high_act_smiles,
                            activity2=high_act_val,
                            structural_difference_description=f"화합물 1({low_act_smiles})과 화합물 2({high_act_smiles})의 구조적 차이점."
                        )
                        
                        try:
                            hypothesis_data = json.loads(json_response)
                            
                            display_md = format_hypothesis_for_markdown(hypothesis_data)
                            
                            file_header = f"""**분석 대상 분자:**\n- **화합물 1 (상대적 저활성):** `{low_act_smiles}` (활성도: {low_act_val:.2f})\n- **화합물 2 (상대적 고활성):** `{high_act_smiles}` (활성도: {high_act_val:.2f})\n\n---\n"""
                            file_md = file_header + display_md

                            st.markdown(file_md, unsafe_allow_html=True)

                            timestamp = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
                            filename = f"hypothesis_pair_{i}_{timestamp}.md"
                            filepath = os.path.join(output_dir, filename)
                            save_hypothesis_to_md(file_md, filepath)
                            st.success(f"가설이 '{filepath}' 파일로 저장되었습니다.")

                        except json.JSONDecodeError:
                            st.error("LLM 응답이 유효한 JSON 형식이 아닙니다. 원본 응답을 표시합니다:")
                            st.text(json_response)

# --- 4. 저장된 가설 관리 (보기/다운로드) ---
st.header("📜 저장된 가설 관리 (보기/다운로드)")

hypotheses_dir = "hypotheses"

if not os.path.isdir(hypotheses_dir):
    st.info("아직 저장된 가설이 없습니다. 가설을 먼저 생성해주세요.")
else:
    files = sorted([f for f in os.listdir(hypotheses_dir) if f.endswith(".md")], reverse=True)
    
    if not files:
        st.info("저장된 가설 파일을 찾을 수 없습니다.")
    else:
        selected_files_to_download = {}
        
        for f in files:
            with st.expander(f):
                # 체크박스를 expander 안에 배치
                selected_files_to_download[f] = st.checkbox("다운로드를 위해 선택", key=f"dl_{f}")
                
                # 파일 내용 읽고 표시
                try:
                    with open(os.path.join(hypotheses_dir, f), "r", encoding="utf-8") as file:
                        content = file.read()
                    st.markdown(content, unsafe_allow_html=True)
                except Exception as e:
                    st.error(f"파일을 읽는 중 오류가 발생했습니다: {e}")

        # 다운로드 버튼은 파일 목록 루프 바깥에 위치
        if st.button("선택된 파일 다운로드"):
            files_to_zip = [f for f, selected in selected_files_to_download.items() if selected]
            
            if not files_to_zip:
                st.warning("다운로드할 파일을 먼저 선택해주세요.")
            else:
                zip_buffer = io.BytesIO()
                with zipfile.ZipFile(zip_buffer, "a", zipfile.ZIP_DEFLATED, False) as zip_file:
                    for f in files_to_zip:
                        file_path = os.path.join(hypotheses_dir, f)
                        zip_file.write(file_path, f)
                
                st.download_button(
                    label="ZIP 파일 다운로드",
                    data=zip_buffer.getvalue(),
                    file_name="hypotheses.zip",
                    mime="application/zip"
                )
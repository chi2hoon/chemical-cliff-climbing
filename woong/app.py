import streamlit as st
import pandas as pd
import datetime
import os
import json
import zipfile
import io
from modules.cheminformatics import find_activity_cliffs
from modules.visualization import visualize_structure_difference, smiles_to_image_b64
from modules.llm_handler import generate_hypothesis, evaluate_hypothesis, revise_hypothesis
from modules.io_utils import load_smiles_activity_csv, save_hypothesis_to_md, parse_hypothesis_md

# --- Helper Functions ---

def get_openai_api_key_from_file():
    try:
        with open("openAI_key.txt", "r") as f:
            return f.read().strip()
    except FileNotFoundError:
        return None

def format_hypothesis_for_markdown(data: dict) -> str:
    """주어진 가설 데이터(dict)를 가독성 좋은 마크다운 및 HTML 문자열로 변환합니다."""
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

    # 설계 제안 (HTML 테이블로 변경)
    md_lines.append("### 💡 검증을 위한 분자 설계 제안")
    suggestions = data.get('design_suggestions', [])
    if suggestions:
        html_table = "<table><thead><tr><th>Design (SMILES)</th><th>Structure</th><th>Expected Effect</th><th>Rationale</th><th>Validation Metric</th></tr></thead><tbody>"
        for s in suggestions:
            smiles = s.get('design', '')
            b64_img = smiles_to_image_b64(smiles) if smiles else ''
            img_tag = f'<img src="data:image/png;base64,{b64_img}" width="200">' if b64_img else ''
            
            html_table += f"<tr>"
            html_table += f"<td>`{smiles}`</td>"
            html_table += f"<td>{img_tag}</td>"
            html_table += f"<td>{s.get('expected_effect', 'N/A')}</td>"
            html_table += f"<td>{s.get('rationale', 'N/A')}</td>"
            html_table += f"<td>{s.get('validation_metric', 'N/A')}</td>"
            html_table += f"</tr>"
        html_table += "</tbody></table>"
        md_lines.append(html_table)
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

def format_evaluation_for_markdown(data: dict) -> str:
    """Evaluation 결과를 마크다운으로 변환합니다."""
    md_lines = []
    summary = data.get('summary', {})
    md_lines.append(f"### 📝 평가 요약")
    md_lines.append(f"- **판정:** {summary.get('verdict', 'N/A')}")
    md_lines.append(f"- **심사 방법:** {summary.get('method_sketch', 'N/A')}")
    md_lines.append("\n")

    details = data.get('detailed_solution', {})
    md_lines.append("### 🔍 상세 평가")
    md_lines.append(f"- **기본 일치성:** {details.get('consistency_check', 'N/A')}")
    md_lines.append(f"- **관점별 검증:** {details.get('aspect_validation', 'N/A')}")
    md_lines.append(f"- **반대 가설 검토:** {details.get('counter_hypothesis_review', 'N/A')}")
    md_lines.append(f"- **설계 제안 검토:** {details.get('design_suggestion_review', 'N/A')}")
    md_lines.append(f"- **추가 필요 요소:** {details.get('additional_requirements', 'N/A')}")
    return "\n".join(md_lines)

# --- Streamlit App ---

st.set_page_config(layout="wide")
st.title("🔬 SAR 분석 및 가설 생성/평가/수정 자동화 도구")
st.write("분자 구조와 활성 데이터 기반의 구조-활성 관계(SAR) 분석 및 가설 생성, 평가, 수정을 자동화합니다.")

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

# --- 5. 가설 평가 및 수정 ---
st.header("🔄 가설 평가 및 수정")

hypotheses_dir = "hypotheses"
if not os.path.isdir(hypotheses_dir) or not os.listdir(hypotheses_dir):
    st.info("평가할 가설이 없습니다. 먼저 가설을 생성해주세요.")
else:
    # Initialize session state keys
    if 'selected_hypothesis_file' not in st.session_state:
        st.session_state.selected_hypothesis_file = None
    if 'evaluation_result' not in st.session_state:
        st.session_state.evaluation_result = None
    if 'revised_hypothesis' not in st.session_state:
        st.session_state.revised_hypothesis = None

    # Get hypothesis files
    hypothesis_files = sorted([f for f in os.listdir(hypotheses_dir) if f.endswith(".md")], reverse=True)
    
    # File selector
    selected_file = st.selectbox("평가/수정할 가설 파일을 선택하세요:", hypothesis_files, index=0)

    # If selection changes, reset the state
    if selected_file != st.session_state.selected_hypothesis_file:
        st.session_state.selected_hypothesis_file = selected_file
        st.session_state.evaluation_result = None
        st.session_state.revised_hypothesis = None

    if st.session_state.selected_hypothesis_file:
        filepath = os.path.join(hypotheses_dir, st.session_state.selected_hypothesis_file)
        
        try:
            with open(filepath, "r", encoding="utf-8") as f:
                content = f.read()
            
            parsed_data = parse_hypothesis_md(content)
            if not all(parsed_data.values()):
                st.error("선택된 파일에서 분자 정보를 파싱할 수 없습니다. 파일 형식을 확인해주세요.")
            else:
                st.subheader("선택된 원본 가설")
                st.markdown(parsed_data['hypothesis_body'], unsafe_allow_html=True)

                # --- Evaluation Step ---
                if st.button("가설 평가 실행"):
                    openai_api_key = get_openai_api_key_from_file()
                    if not openai_api_key:
                        st.warning("API 키를 openAI_key.txt 파일에서 로드해주세요.")
                    else:
                        with st.spinner("LLM이 가설을 평가 중입니다..."):
                            eval_response = evaluate_hypothesis(
                                api_key=openai_api_key,
                                hypothesis_text=parsed_data['hypothesis_body'],
                                smiles1=parsed_data['smiles1'],
                                activity1=parsed_data['activity1'],
                                smiles2=parsed_data['smiles2'],
                                activity2=parsed_data['activity2'],
                                structural_difference_description=""
                            )
                            st.session_state.evaluation_result = eval_response
                            st.session_state.revised_hypothesis = None # Reset revision

                # --- Display Evaluation and Trigger Revision ---
                if st.session_state.evaluation_result:
                    st.subheader("가설 평가 결과")
                    try:
                        eval_data = json.loads(st.session_state.evaluation_result)
                        st.markdown(format_evaluation_for_markdown(eval_data))

                        if st.button("평가 기반으로 가설 수정"):
                            openai_api_key = get_openai_api_key_from_file()
                            if not openai_api_key:
                                st.warning("API 키를 openAI_key.txt 파일에서 로드해주세요.")
                            else:
                                with st.spinner("LLM이 가설을 수정 중입니다..."):
                                    revise_response = revise_hypothesis(
                                        api_key=openai_api_key,
                                        original_hypothesis_text=parsed_data['hypothesis_body'],
                                        review_findings=st.session_state.evaluation_result,
                                        smiles1=parsed_data['smiles1'],
                                        activity1=parsed_data['activity1'],
                                        smiles2=parsed_data['smiles2'],
                                        activity2=parsed_data['activity2'],
                                        structural_difference_description=""
                                    )
                                    st.session_state.revised_hypothesis = revise_response

                    except json.JSONDecodeError:
                        st.error("LLM 평가 응답이 유효한 JSON이 아닙니다.")
                        st.text(st.session_state.evaluation_result)

                # --- Display Revision and Save ---
                if st.session_state.revised_hypothesis:
                    st.subheader("수정된 가설")
                    try:
                        revised_data = json.loads(st.session_state.revised_hypothesis)
                        
                        # Format the revised hypothesis for display
                        display_md = format_hypothesis_for_markdown(revised_data)
                        st.markdown(display_md, unsafe_allow_html=True)

                        # Prepare for saving
                        file_header = f"""**분석 대상 분자:**\n- **화합물 1 (상대적 저활성):** `{parsed_data['smiles1']}` (활성도: {parsed_data['activity1']:.2f})\n- **화합물 2 (상대적 고활성):** `{parsed_data['smiles2']}` (활성도: {parsed_data['activity2']:.2f})\n\n---\n"""
                        final_md_to_save = file_header + display_md

                        st.subheader("수정된 가설 저장")
                        new_filename = st.text_input("새 파일 이름:", f"revised_{st.session_state.selected_hypothesis_file}")
                        if st.button("수정된 가설 저장"):
                            if new_filename:
                                new_filepath = os.path.join(hypotheses_dir, new_filename)
                                save_hypothesis_to_md(final_md_to_save, new_filepath)
                                st.success(f"수정된 가설이 {new_filepath}에 저장되었습니다.")
                            else:
                                st.warning("파일 이름을 입력하세요.")

                    except json.JSONDecodeError:
                        st.error("LLM 수정 응답이 유효한 JSON이 아닙니다.")
                        st.text(st.session_state.revised_hypothesis)

        except Exception as e:
            st.error(f"파일 처리 중 오류 발생: {e}")

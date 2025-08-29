import streamlit as st
import pandas as pd
import datetime
import os
import json
import zipfile
import io
from modules.cheminformatics import find_activity_cliffs
from modules.visualization import visualize_structure_difference, smiles_to_image_b64
from modules.llm_handler import generate_hypothesis, evaluate_hypothesis, revise_hypothesis, create_activity_summary
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
            html_table += "</tr>"
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

st.set_page_config(layout="centered")
st.title("🔬 SAR 분석 및 가설 생성/평가/수정 자동화 도구")
st.write("분자 구조와 활성 데이터 기반의 구조-활성 관계(SAR) 분석 및 가설 생성, 평가, 수정을 자동화합니다.")

tab1, tab2, tab3, tab4, tab5 = st.tabs([
    "1. 데이터 업로드",
    "2. Activity Cliff 분석",
    "3. 가설 생성",
    "4. 가설 관리",
    "5. 가설 평가 및 수정"
])

with tab1:
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

        # 2017년 데이터 빠른 로딩 (권장)
    st.markdown("### 🚀 2017년 데이터 빠른 시작")
    col1, col2, col3 = st.columns([2, 2, 1])

    with col1:
        if st.button("📊 2017년 전체 데이터 로드", type="primary", use_container_width=True):
            try:
                with st.spinner("2017년 데이터를 불러오는 중..."):
                    # 2017년 데이터의 모든 그룹을 하나로 합치기
                    data_root = "hoon/data"
                    groups_2017 = list_hoon_groups(data_root, "2017")

                    if groups_2017.empty:
                        st.error("2017년 데이터가 없습니다.")
                    else:
                        all_dfs = []
                        for _, group in groups_2017.iterrows():
                            assay_id = group['assay_id'] if group['assay_id'] else None
                            cell_line = group['cell_line'] if group['cell_line'] else None

                            df_group = load_hoon_group_as_dataframe(
                                dataset_id="2017",
                                assay_id=assay_id,
                                cell_line=cell_line,
                                data_root=data_root
                            )
                            if not df_group.empty:
                                # 그룹 정보 추가
                                df_group['assay_id'] = assay_id or ""
                                df_group['cell_line'] = cell_line or ""
                                all_dfs.append(df_group)

                        if all_dfs:
                            df_combined = pd.concat(all_dfs, ignore_index=True)
                            st.session_state['df'] = df_combined
                            st.session_state['auto_suggestion'] = {"smiles_col": "SMILES", "activity_col": "Activity"}
                            st.success(f"2017년 전체 데이터 로드 완료! 총 {len(df_combined)}개 레코드")
                            st.dataframe(df_combined.head())
                        else:
                            st.warning("불러올 수 있는 데이터가 없습니다.")

            except Exception as e:
                st.error(f"데이터 로드 실패: {e}")

    with col2:
        st.markdown("**데이터 정보:**")
        st.info("2017년 cytotoxicity 데이터\n• 5개 세포주 (253J, 5637, J82, KU-19-19, MBT-2)\n• 각 그룹당 ~107개 화합물\n• 총 ~535개 데이터 포인트")

    with col3:
        if st.button("🔄 새로고침"):
            st.rerun()

    st.markdown("---")
    with st.expander("또는 세부 그룹 선택 (고급)"):
        data_root = "hoon/data"
        ds_list = list_hoon_datasets(data_root)
        if not ds_list:
            st.info("hoon/data/silver 하위에 데이터셋이 없습니다. 먼저 UDM 파이프라인을 실행하세요 (예: `python hoon/udm_cli.py silver --config hoon/configs/2017.yml --root hoon/data`).")
        else:
            colA, colB = st.columns([1, 2])
            with colA:
                ds = st.selectbox("데이터셋 선택", ds_list, index=ds_list.index("2017") if "2017" in ds_list else 0, key="hoon_ds")
                groups_df = list_hoon_groups(data_root, ds)
                # 안전한 기본값 처리
                assay_opts = sorted(groups_df['assay_id'].unique().tolist()) if 'assay_id' in groups_df.columns else []
                cell_opts = sorted(groups_df['cell_line'].unique().tolist()) if 'cell_line' in groups_df.columns else []
            with colB:
                assay_sel = st.selectbox("Assay ID", [""] + assay_opts, index=0, key="hoon_assay")
                cell_sel = st.selectbox("Cell line", [""] + cell_opts, index=0, key="hoon_cell")

            if st.button("선택한 그룹 불러오기 (SMILES/Activity)"):
                try:
                    df2 = load_hoon_group_as_dataframe(dataset_id=ds, assay_id=(assay_sel or None), cell_line=(cell_sel or None), data_root=data_root)
                    if df2 is None or df2.empty:
                        st.warning("선택한 조건에 해당하는 데이터가 없습니다.")
                    else:
                        st.success(f"불러오기 성공: {len(df2)}개 레코드")
                        st.dataframe(df2.head())
                        st.session_state['df'] = df2
                        st.session_state['auto_suggestion'] = {"smiles_col": "SMILES", "activity_col": "Activity"}
                except Exception as e:
                    st.error(f"불러오기 실패: {e}")

with tab2:
    st.header("2. Activity Cliff 분석")
    # 데이터 소스 선택
    source = st.radio("데이터 소스", ["업로드/불러온 데이터로 계산", "Hoon AC 쌍(사전계산) 불러오기"], index=0, horizontal=True)

    if source == "Hoon AC 쌍(사전계산) 불러오기":
        if st.button("Hoon 사전계산 AC 쌍 로드"):
            try:
                cliff_df = load_hoon_ac_pairs(data_root="hoon/data")
                if cliff_df is None or cliff_df.empty:
                    st.info("사전계산된 AC 쌍이 없습니다. hoon 파이프라인에서 `ac` 또는 `ac-all` 스테이지를 먼저 실행하세요.")
                else:
                    st.success(f"{len(cliff_df)}개의 Activity Cliff 쌍을 불러왔습니다.")
                    st.dataframe(cliff_df.head())
                    st.session_state['cliff_df'] = cliff_df
            except Exception as e:
                st.error(f"AC 쌍 로드 실패: {e}")
    elif 'df' in st.session_state and st.session_state['df'] is not None:
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

        activity_assumption = st.radio(
            "활성도 데이터의 의미를 선택해주세요:",
            ('값이 높을수록 활성도가 높음 (Higher is better)', '값이 낮을수록 활성도가 높음 (Lower is better)'),
            key='activity_assumption'
        )

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
                    activity_diff_threshold=activity_diff_threshold,
                    higher_is_better= (st.session_state.activity_assumption == '값이 높을수록 활성도가 높음 (Higher is better)')
                )
            
            st.success(f"{len(cliff_df)}개의 Activity Cliff 쌍을 찾았습니다!")
            st.dataframe(cliff_df)
            st.session_state['cliff_df'] = cliff_df

            # 활성도 지표에 대한 요약 정보 표시
            summary_text = create_activity_summary(activity_col, (st.session_state.activity_assumption == '값이 높을수록 활성도가 높음 (Higher is better)'))
            st.markdown("---")
            st.markdown(summary_text)
    else:
        st.info("1. 데이터 업로드 탭에서 데이터를 먼저 업로드해주세요.")

with tab3:
    st.header("3. 결과 시각화 및 가설 생성")
    if 'cliff_df' in st.session_state and not st.session_state['cliff_df'].empty:
        cliff_df = st.session_state['cliff_df']
        
        st.subheader("분석할 Activity Cliff 쌍 선택")
        st.dataframe(cliff_df)
        selected_indices = st.multiselect("분석 및 시각화할 Activity Cliff 쌍의 인덱스를 선택하세요:", cliff_df.index)
        
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
                        
                        img = visualize_structure_difference(
                            smiles1=row['SMILES_1'],
                            smiles2=row['SMILES_2'],
                            legend1=f"SMILES: {row['SMILES_1']}\nActivity: {row['Activity_1']:.2f}",
                            legend2=f"SMILES: {row['SMILES_2']}\nActivity: {row['Activity_2']:.2f}"
                        )
                        st.image(img, caption=f"유사도: {row['Similarity']:.3f} | 활성도 차이: {row['Activity_Diff']:.2f}")

                        with st.spinner(f"쌍 #{i}에 대한 LLM 가설을 생성 중입니다..."):
                            higher_is_better = st.session_state.get('activity_assumption') == '값이 높을수록 활성도가 높음 (Higher is better)'
                            
                            # 가정에 따라 고활성/저활성 분자 결정
                            if (higher_is_better and row['Activity_1'] > row['Activity_2']) or \
                               (not higher_is_better and row['Activity_1'] < row['Activity_2']):
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
                                structural_difference_description=f"화합물 1({low_act_smiles})과 화합물 2({high_act_smiles})의 구조적 차이점.",
                                similarity=row['Similarity']
                            )
                            
                            try:
                                hypothesis_data = json.loads(json_response)
                                display_md = format_hypothesis_for_markdown(hypothesis_data)
                                file_header = f"**분석 대상 분자:**\n- **화합물 1 (상대적 저활성):** `{low_act_smiles}` (활성도: {low_act_val:.2f})\n- **화합물 2 (상대적 고활성):** `{high_act_smiles}` (활성도: {high_act_val:.2f})\n\n---\n"
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
    else:
        st.info("2. Activity Cliff 분석 탭에서 분석을 먼저 실행해주세요.")

with tab4:
    st.header("📜 저장된 가설 관리 (보기/다운로드)")
    hypotheses_dir = "hypotheses"

    if st.button("🔄 목록 새로고침"):
        st.rerun()

    if not os.path.isdir(hypotheses_dir) or not os.listdir(hypotheses_dir):
        st.info("아직 저장된 가설이 없습니다. 3. 가설 생성 탭에서 가설을 먼저 생성해주세요.")
    else:
        files = sorted([f for f in os.listdir(hypotheses_dir) if f.endswith(".md")], reverse=True)
        
        if not files:
            st.info("저장된 가설 파일을 찾을 수 없습니다.")
        else:
            selected_file = st.selectbox("확인할 가설 파일을 선택하세요:", files)

            if selected_file:
                filepath = os.path.join(hypotheses_dir, selected_file)
                try:
                    with open(filepath, "r", encoding="utf-8") as file:
                        content = file.read()
                    
                    st.markdown("---")
                    st.subheader(f"📄 {selected_file}")
                    st.markdown(content, unsafe_allow_html=True)

                    st.download_button(
                        label=f"'{selected_file}' 다운로드",
                        data=content,
                        file_name=selected_file,
                        mime="text/markdown"
                    )

                except Exception as e:
                    st.error(f"파일을 읽는 중 오류가 발생했습니다: {e}")

with tab5:
    st.header("🔄 가설 평가 및 수정")
    hypotheses_dir = "hypotheses"

    if not os.path.isdir(hypotheses_dir) or not os.listdir(hypotheses_dir):
        st.info("평가할 가설이 없습니다. 3. 가설 생성 탭에서 가설을 먼저 생성해주세요.")
    else:
        hypothesis_files = sorted([f for f in os.listdir(hypotheses_dir) if f.endswith(".md")], reverse=True)
        
        # --- 1. 가설 파일 목록 표시 ---
        st.subheader("📋 가설 목록")
        selected_file = st.radio(
            "평가할 가설을 선택하세요:", 
            hypothesis_files, 
            key="selected_hypothesis_file",
            label_visibility="collapsed"
        )
        st.markdown("---")

        # --- 2. 선택된 가설 내용 및 평가/수정 ---
        if selected_file:
            st.subheader(f"📄 원본 가설: {selected_file}")
            filepath = os.path.join(hypotheses_dir, selected_file)

            try:
                with open(filepath, "r", encoding="utf-8") as file:
                    content = file.read()
                
                parsed_data = parse_hypothesis_md(content)
                if not all(parsed_data.values()):
                    st.error(f"파일({selected_file})에서 분자 정보를 파싱할 수 없습니다. 파일 형식을 확인해주세요.")
                    st.stop()

                # 원본 가설 표시
                st.markdown(content, unsafe_allow_html=True)
                st.markdown("---")

                # --- 가설 평가 ---
                eval_container = st.container(border=True)
                with eval_container:
                    st.subheader("🔬 가설 평가")
                    eval_key = f"eval_{selected_file}"
                    
                    if st.button("가설 평가 실행", key=f"eval_btn_{selected_file}"):
                        openai_api_key = get_openai_api_key_from_file()
                        if not openai_api_key:
                            st.warning("API 키를 openAI_key.txt 파일에서 로드해주세요.")
                        else:
                            with st.spinner("LLM이 가설을 평가 중입니다..."):
                                eval_response = evaluate_hypothesis(
                                    api_key=openai_api_key,
                                    hypothesis_text=parsed_data['hypothesis_body'],
                                    smiles1=parsed_data['smiles1'], activity1=parsed_data['activity1'],
                                    smiles2=parsed_data['smiles2'], activity2=parsed_data['activity2'],
                                    structural_difference_description=""
                                )
                                st.session_state[eval_key] = eval_response
                                # 새로운 평가가 시작되면 이전 수정 결과는 삭제
                                if f"revise_{selected_file}" in st.session_state:
                                    del st.session_state[f"revise_{selected_file}"]
                    
                    if eval_key in st.session_state:
                        st.markdown("##### 평가 결과")
                        try:
                            eval_data = json.loads(st.session_state[eval_key])
                            formatted_eval_md = format_evaluation_for_markdown(eval_data)
                            st.markdown(formatted_eval_md, unsafe_allow_html=True)

                            # 평가 결과 저장 섹션
                            st.markdown("---")
                            evaluations_dir = "evaluations"
                            os.makedirs(evaluations_dir, exist_ok=True)
                            
                            base_filename = selected_file.replace(".md", "")
                            eval_filename = f"{base_filename}_Eval.md"
                            
                            content_to_save = f"# 원본 가설: {selected_file}\n\n{content}\n\n---\n\n# 가설 평가 결과\n\n{formatted_eval_md}"
                            eval_filepath = os.path.join(evaluations_dir, eval_filename)
                            save_hypothesis_to_md(content_to_save, eval_filepath)
                            st.success(f"평가 결과가 '{eval_filepath}'에 저장되었습니다.")
                            st.toast("평가 저장 완료!")

                        except json.JSONDecodeError:
                            st.error("LLM 평가 응답이 유효한 JSON이 아닙니다.")
                            st.text(st.session_state[eval_key])

                # --- 가설 수정 및 저장 ---
                if eval_key in st.session_state:
                    revise_container = st.container(border=True)
                    with revise_container:
                        st.subheader("✍️ 평가 기반 가설 수정")
                        revise_key = f"revise_{selected_file}"

                        if st.button("평가 기반으로 가설 수정", key=f"revise_btn_{selected_file}"):
                            openai_api_key = get_openai_api_key_from_file()
                            if not openai_api_key:
                                st.warning("API 키를 openAI_key.txt 파일에서 로드해주세요.")
                            else:
                                with st.spinner("LLM이 가설을 수정 중입니다..."):
                                    revise_response = revise_hypothesis(
                                        api_key=openai_api_key,
                                        original_hypothesis_text=parsed_data['hypothesis_body'],
                                        review_findings=st.session_state[eval_key],
                                        smiles1=parsed_data['smiles1'], activity1=parsed_data['activity1'],
                                        smiles2=parsed_data['smiles2'], activity2=parsed_data['activity2'],
                                        structural_difference_description=""
                                    )
                                    st.session_state[revise_key] = revise_response
                        
                        if revise_key in st.session_state:
                            st.markdown("##### 수정된 가설")
                            try:
                                revised_data = json.loads(st.session_state[revise_key])
                                display_md = format_hypothesis_for_markdown(revised_data)
                                st.markdown(display_md, unsafe_allow_html=True)

                                # 수정된 가설 저장 섹션
                                st.markdown("---")
                                st.subheader("💾 수정된 가설 저장")
                                
                                file_header = f"**분석 대상 분자:**\n- **화합물 1 (상대적 저활성):** `{parsed_data['smiles1']}` (활성도: {parsed_data['activity1']:.2f})\n- **화합물 2 (상대적 고활성):** `{parsed_data['smiles2']}` (활성도: {parsed_data['activity2']:.2f})\n\n---"
                                final_md_to_save = file_header + display_md

                                new_filename = f"revised_{selected_file}"
                                new_filepath = os.path.join(hypotheses_dir, new_filename)
                                save_hypothesis_to_md(final_md_to_save, new_filepath)
                                st.success(f"수정된 가설이 '{new_filepath}'에 저장되었습니다.")
                                st.toast("저장 완료!")

                            except json.JSONDecodeError:
                                st.error("LLM 수정 응답이 유효한 JSON이 아닙니다.")
                                st.text(st.session_state[revise_key])

            except Exception as e:
                st.error(f"파일 처리 중 오류가 발생했습니다: {e}")


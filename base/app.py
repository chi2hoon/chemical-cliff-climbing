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
from modules.io_utils import (
    load_smiles_activity_csv,
    save_hypothesis_to_md,
    parse_hypothesis_md,
    load_gold_data,
    get_available_gold_years,
    get_cell_lines_for_panel,
    get_available_targets,
    get_available_panel_ids
)

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


tab1, tab2, tab3, tab4, tab5, tab6 = st.tabs([
    "1. 데이터셋 선택",
    "2. Activity Cliff 분석",
    "3. 가설 생성",
    "4. 가설 관리",
    "5. 가설 평가 및 수정",
    "6. 자동 가설 수정"
])

with tab1:
    st.header("1. 데이터 로드")
    st.markdown("표준화된 데이터셋을 로드하여 분석을 시작하세요.")

    # 데이터셋 선택 (앱 파일 위치를 기준으로 고정 경로 구성)
    base_dir = os.path.dirname(os.path.abspath(__file__))
    data_root = os.path.join(base_dir, "data")
    available_years_all = get_available_gold_years(data_root)

    if not available_years_all:
        st.warning("분석용 데이터가 없습니다. 터미널에서 `PYTHONPATH=base python -m pipeline.cli gold --years 2017 2018`로 생성하세요.")
    else:
        # 패널 이름 매핑
        panel_names_map = {
            "blca": "방광암세포주 패널",
            "prad": "전립선암세포주 패널", 
            "luad": "폐암세포주 패널",
            "brca": "유방암세포주 패널",
            "heme": "혈액암세포주 패널",
            "paad": "췌장암세포주 패널",
            "coad": "대장암세포주 패널",
            "misc12": "뇌암/기타 패널",
            "misc13": "기타 패널"
        }

        # 년도(왼쪽) - 패널/타겟(오른쪽)
        col_year, col_panel = st.columns([1, 2])

        with col_year:
            selected_year = st.selectbox("📅 데이터셋 년도", sorted(available_years_all), index=0)

        selected_panel = None
        selected_target = None

        with col_panel:
            # 연도의 데이터 특성을 감지: 패널 목록을 파생할 수 있으면 패널 기반으로 간주
            filtered_panel_ids = get_available_panel_ids(selected_year, data_root)
            if filtered_panel_ids:
                # 해당 연도는 패널 기반(2017은 target_id에서 파생)
                panel_display_options = ["전체 패널"]
                panel_id_to_display = {"전체 패널": None}
                for panel_id in sorted(filtered_panel_ids):
                    display_name = panel_names_map.get(panel_id, panel_id)
                    # 중복 라벨 방지: 같으면 한 번만 표기
                    if str(display_name).strip() == str(panel_id).strip():
                        display_option = f"{display_name}"
                    else:
                        display_option = f"{panel_id} ({display_name})"
                    panel_display_options.append(display_option)
                    panel_id_to_display[display_option] = panel_id
                selected_panel_display = st.selectbox("🧬 패널 선택", panel_display_options, index=0)
                selected_panel = panel_id_to_display[selected_panel_display]
            else:
                # 해당 연도는 타겟 기반
                targets = get_available_targets(selected_year, data_root)
                target_display_options = ["전체 타겟"] + targets
                selected_target_display = st.selectbox("🎯 타겟 선택", target_display_options, index=0)
                selected_target = None if selected_target_display == "전체 타겟" else selected_target_display

        # 세포주 셀렉터 (패널 선택 시)
        selected_cell_line = None
        if selected_panel:
            cell_lines = get_cell_lines_for_panel(selected_year, selected_panel, data_root)
            if cell_lines:
                selected_cell_line = st.selectbox("🧫 세포주 선택", ["전체 세포주"] + cell_lines, index=0)
                if selected_cell_line == "전체 세포주":
                    selected_cell_line = None

        # 로드 버튼 - selected_year가 있을 때만 표시
        if selected_year:
            st.markdown("### 🚀 데이터 로드")
            # 데이터 설명 토글(전체 너비)
            _desc = None
            if str(selected_year) == "2017":
                _desc = (
                    "2017: 이 데이터는 국립암센터와 한국화학연구원이 출원한 특허 패밀리(KR101920163B1, PCT/WO2018021849A1)의 부속 표를 정규화한 세트입니다.\n"
                    "c-Myc/Max/DNA 결합 억제 계열 화합물의 세포독성(IC₅₀, μM) 결과가 암종별 여러 세포주에서 보고되며, 표 일부에는 비교 화합물과 후보군의 독성 관련 보조 정보가 함께 제시됩니다."
                )
            elif str(selected_year) == "2018":
                _desc = (
                    "2018: 이 데이터는 PCT/EP2018/056824, WO 2018/172250 특허의 부속 표를 정규화한 세트입니다.\n"
                    "2-메틸 퀴나졸린 계열 화합물의 Ras–SOS1 상호작용 억제를 주 타깃으로 한 in vitro HTRF 기반 생화학 어세이(Assay 1~3)와 EGFR 키나아제 어세이 결과를 포함합니다.\n"
                    "측정값은 IC₅₀ 또는 20 µM 단일 농도에서의 % 억제율로 보고되며, 표에는 ‘n.d.’(not determined) 표기와 EGFR 선택성 관련 정보가 함께 제시됩니다."
                )
            with st.expander("데이터 설명", expanded=False):
                st.markdown(_desc or "해당 연도에 대한 설명은 준비 중입니다.")
            load_text = f"{selected_year}년 데이터 로드"
            if selected_panel:
                panel_name = panel_names_map.get(selected_panel, selected_panel)
                load_text += f" ({panel_name})"
            elif selected_target:
                load_text += f" (target: {selected_target})"

            if st.button(f"📊 {load_text}", type="primary", use_container_width=True):
                try:
                    with st.spinner(f"{selected_year}년 데이터를 불러오는 중..."):
                        df_gold = load_gold_data(
                            year=selected_year, 
                            data_root=data_root, 
                            panel_id=selected_panel,
                            cell_line=selected_cell_line,
                            target_id=selected_target
                        )

                        if df_gold.empty:
                            if selected_panel:
                                st.error(f"{selected_year}년 {selected_panel} 패널 데이터가 없습니다.")
                            elif selected_target:
                                st.error(f"{selected_year}년 target '{selected_target}' 데이터가 없습니다.")
                            else:
                                st.error(f"{selected_year}년 데이터가 없습니다.")
                        else:
                            st.session_state['df'] = df_gold
                            st.session_state['auto_suggestion'] = {"smiles_col": "SMILES", "activity_col": "Activity"}
                            
                            success_msg = f"{selected_year}년 데이터 로드 완료! 총 {len(df_gold)}개 레코드"
                            if selected_panel:
                                success_msg += f" ({panel_names_map.get(selected_panel, selected_panel)})"
                            elif selected_target:
                                success_msg += f" (target: {selected_target})"
                            
                            st.success(success_msg)
                            st.dataframe(df_gold.head())

                            # 데이터 스키마 정보 표시
                            st.info("**데이터 스키마:**\n"
                                   "• SMILES: 표준화된 캐노니컬 SMILES\n"
                                   "• Activity: 표준화된 활성도 값 (value_std)\n"
                                   "• 메타데이터: assay_id, panel_id/target_id, cell_line, inchikey 등")

                except Exception as e:
                    st.error(f"데이터 로드 실패: {e}")

    # (레거시) Gold 데이터 설명/파이프라인 정보 섹션 제거됨

with tab2:
    st.header("2. Activity Cliff 분석")
    if 'df' in st.session_state and st.session_state['df'] is not None:
        df = st.session_state['df']
        
        # 단위 분포 요약 표시
        if 'unit_std' in df.columns:
            try:
                counts = df['unit_std'].fillna('NaN').value_counts().to_dict()
                counts_str = ", ".join([f"{k}: {v}" for k, v in counts.items()])
                st.caption(f"단위 분포 요약: {counts_str}")
            except Exception:
                pass

        # 컬럼 자동 선택(고정): Gold 스키마 가정하에 자동 결정
        smiles_col = 'SMILES' if 'SMILES' in df.columns else df.columns[0]
        activity_col = 'Activity' if 'Activity' in df.columns else ( 'value_std' if 'value_std' in df.columns else (df.columns[1] if len(df.columns) > 1 else df.columns[0]))
        st.caption(f"자동 선택된 컬럼: SMILES='{smiles_col}', Activity='{activity_col}'")

        col1, col2 = st.columns(2)
        with col1:
            scale_choice = st.selectbox("활성도 스케일", ["원본(단위 유지)", "pAct (-log10[M])"], index=1)

        with col2:
            similarity_threshold = st.slider("구조 유사도 임계값 (Tanimoto)", 0.7, 1.0, 0.85, 0.01)
            # 스케일별 임계값 자동 제안
            default_diff = 1.0; step = 0.1; fmt = "%f"
            try:
                import numpy as _np
                if scale_choice == "pAct (-log10[M])":
                    default_diff = 0.5
                    step = 0.1
                    fmt = "%0.2f"
                else:
                    vals = pd.to_numeric(df[activity_col], errors="coerce").dropna()
                    if not vals.empty:
                        q10, q90 = _np.percentile(vals, [10, 90])
                        span = max(0.0, float(q90 - q10))
                        default_diff = max(1e-9, span * 0.25)
                        step = max(1e-9, default_diff / 10.0)
                        fmt = "%e"
            except Exception:
                pass
            activity_diff_threshold = st.number_input("활성도 차이 임계값", min_value=0.0, value=float(default_diff), step=float(step), format=fmt)
            if scale_choice == "pAct (-log10[M])":
                st.caption("기본 0.5 pAct ≈ 약 3.16배 potency 차이. 0.3(≈2배)도 자주 사용합니다.")
            else:
                st.caption("원본 단위 스케일에서는 분포에 따라 자동 제안됩니다. 필요 시 조정하세요.")

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
                # 스케일 변환 적용
                chosen_activity_col = activity_col
                if scale_choice == "pAct (-log10[M])":
                    unit_col = "unit_std" if "unit_std" in work_df.columns else None
                    def _to_m(val, unit):
                        try:
                            v = float(val)
                        except Exception:
                            return None
                        u = (str(unit) if unit is not None else "").strip()
                        if u == "M":
                            return v
                        if u == "mM":
                            return v * 1e-3
                        if u == "uM":
                            return v * 1e-6
                        if u == "nM":
                            return v * 1e-9
                        return v * 1e-6  # 폴백: uM 가정
                    def _to_pact(val_m):
                        try:
                            import math
                            if val_m is None or float(val_m) <= 0:
                                return None
                            return -math.log10(float(val_m))
                        except Exception:
                            return None
                    vals_m = [ _to_m(work_df.iloc[i][activity_col], (work_df.iloc[i][unit_col] if unit_col else "uM")) for i in range(len(work_df)) ]
                    pact = [ _to_pact(x) for x in vals_m ]
                    work_df["Activity_pAct"] = pact
                    work_df = work_df[pd.Series(work_df["Activity_pAct"]).notna()].reset_index(drop=True)
                    chosen_activity_col = "Activity_pAct"

                cliff_df = find_activity_cliffs(
                    work_df,
                    smiles_col=smiles_col,
                    activity_col=chosen_activity_col,
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
        st.info("1. 상단에서 Gold 데이터를 먼저 로드해주세요.")

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

with tab6:
    st.header("🤖 자동 가설 수정")
    hypotheses_dir = "hypotheses"

    if not os.path.isdir(hypotheses_dir) or not os.listdir(hypotheses_dir):
        st.info("자동 수정할 가설이 없습니다. 3. 가설 생성 탭에서 가설을 먼저 생성해주세요.")
    else:
        hypothesis_files = sorted([f for f in os.listdir(hypotheses_dir) if f.endswith(".md")], reverse=True)
        
        selected_file_auto = st.selectbox(
            "자동 수정을 시작할 가설을 선택하세요:", 
            hypothesis_files, 
            key="selected_hypothesis_file_auto"
        )
        
        col1, col2 = st.columns(2)
        with col1:
            min_iterations = st.number_input("최소 반복 횟수:", min_value=1, max_value=10, value=1, step=1)
        with col2:
            max_iterations = st.number_input("최대 반복 횟수:", min_value=1, max_value=10, value=3, step=1)

        if st.button("🤖 자동 수정 시작", key="auto_revise_start"):
            openai_api_key = get_openai_api_key_from_file()
            if not openai_api_key:
                st.warning("API 키를 openAI_key.txt 파일에서 로드해주세요.")
                st.stop()

            filepath = os.path.join(hypotheses_dir, selected_file_auto)
            try:
                with open(filepath, "r", encoding="utf-8") as file:
                    content = file.read()
                
                parsed_data = parse_hypothesis_md(content)
                if not all(parsed_data.values()):
                    st.error(f"파일({selected_file_auto})에서 분자 정보를 파싱할 수 없습니다. 파일 형식을 확인해주세요.")
                    st.stop()

            except Exception as e:
                st.error(f"원본 가설 파일을 읽는 중 오류가 발생했습니다: {e}")
                st.stop()

            current_hypothesis_body = parsed_data['hypothesis_body']
            final_hypothesis_body = ""

            with st.status(f"'{selected_file_auto}'에 대한 자동 수정을 시작합니다...", expanded=True) as status:
                for i in range(max_iterations):
                    st.write(f"---")
                    st.write(f"**🚀 반복 {i+1}/{max_iterations}**")
                    
                    # 1. 평가
                    st.write("1️⃣ 가설을 평가합니다...")
                    try:
                        eval_response = evaluate_hypothesis(
                            api_key=openai_api_key,
                            hypothesis_text=current_hypothesis_body,
                            smiles1=parsed_data['smiles1'], activity1=parsed_data['activity1'],
                            smiles2=parsed_data['smiles2'], activity2=parsed_data['activity2'],
                            structural_difference_description=""
                        )
                        eval_data = json.loads(eval_response)
                        verdict = eval_data.get('summary', {}).get('verdict', 'Unknown').upper()
                        
                        with st.expander("평가 결과 보기"):
                            st.markdown(format_evaluation_for_markdown(eval_data), unsafe_allow_html=True)

                    except Exception as e:
                        st.error(f"반복 {i+1}에서 평가 중 오류 발생: {e}")
                        status.update(label="오류로 인해 중단됨", state="error")
                        st.stop()

                    # 2. 판정에 따른 분기
                    st.write(f"2️⃣ 평가 판정: **{verdict}**")
                    if ("GOOD" in verdict or "SOUND" in verdict) and (i + 1) >= min_iterations:
                        st.success(f"✅ 가설이 'Good' 또는 'Unsound'로 판정되고 최소 반복 횟수({min_iterations})에 도달하여 프로세스를 종료합니다.")
                        final_hypothesis_body = current_hypothesis_body
                        status.update(label="자동 수정 완료!", state="complete")
                        break
                    
                    elif "WEAK" in verdict:
                        st.info("🤔 가설이 'Weak'로 판정되어 수정을 진행합니다.")
                        # 3. 수정
                        st.write("3️⃣ 평가 기반으로 가설을 수정합니다...")
                        try:
                            revise_response = revise_hypothesis(
                                api_key=openai_api_key,
                                original_hypothesis_text=current_hypothesis_body,
                                review_findings=eval_response,
                                smiles1=parsed_data['smiles1'], activity1=parsed_data['activity1'],
                                smiles2=parsed_data['smiles2'], activity2=parsed_data['activity2'],
                                structural_difference_description=""
                            )
                            revised_data = json.loads(revise_response)
                            current_hypothesis_body = format_hypothesis_for_markdown(revised_data)
                            
                            with st.expander("수정된 가설 내용 보기"):
                                st.markdown(current_hypothesis_body, unsafe_allow_html=True)

                        except Exception as e:
                            st.error(f"반복 {i+1}에서 수정 중 오류 발생: {e}")
                            status.update(label="오류로 인해 중단됨", state="error")
                            st.stop()
                    else:
                        st.warning(f"⚠️ 알 수 없는 판정('{verdict}')으로 인해 프로세스를 중단합니다.")
                        final_hypothesis_body = current_hypothesis_body
                        status.update(label="알 수 없는 판정으로 중단됨", state="error")
                        break
                
                else: # for-else loop: break 없이 끝났을 경우
                    st.warning(f"🔔 최대 반복 횟수({max_iterations})에 도달했습니다.")
                    final_hypothesis_body = current_hypothesis_body
                    status.update(label="최대 반복 후 완료", state="complete")

            # 최종 결과 저장
            if final_hypothesis_body:
                st.markdown("---")
                st.subheader("💾 최종 결과 저장")
                
                file_header = f"**분석 대상 분자:**\n- **화합물 1 (상대적 저활성):** `{parsed_data['smiles1']}` (활성도: {parsed_data['activity1']:.2f})\n- **화합물 2 (상대적 고활성):** `{parsed_data['smiles2']}` (활성도: {parsed_data['activity2']:.2f})\n\n---"
                final_md_to_save = file_header + final_hypothesis_body

                # 새 파일명 생성
                base_name = selected_file_auto.replace('.md', '')
                if base_name.startswith('auto_'): # 기존 접두사 제거
                    base_name = base_name[5:]
                
                new_filename = f"auto_{base_name}.md"
                new_filepath = os.path.join(hypotheses_dir, new_filename)
                
                save_hypothesis_to_md(final_md_to_save, new_filepath)
                st.success(f"최종 수정된 가설이 '{new_filepath}'에 저장되었습니다.")
                
                st.markdown("##### 최종 가설 내용:")
                st.markdown(final_md_to_save, unsafe_allow_html=True)

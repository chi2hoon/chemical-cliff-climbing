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
from modules.context_builder import build_pair_context
from modules.io_utils import (
    load_smiles_activity_csv,
    save_hypothesis_to_md,
    parse_hypothesis_md,
    load_gold_data,
    get_available_gold_years,
    get_available_targets,
    get_available_panel_ids,
    get_all_available_panels_and_years,
    get_cell_lines_for_panel
)

# --- Helper Functions ---

def get_openai_api_key_from_file():
    # 1) 환경변수 우선
    env_key = os.getenv("OPENAI_API_KEY")
    if env_key and env_key.strip():
        return env_key.strip()
    # 2) 현재 작업 디렉토리(base)에서 탐색
    try:
        with open("openAI_key.txt", "r") as f:
            return f.read().strip()
    except FileNotFoundError:
        pass
    # 3) 레포 루트(../)에서도 탐색
    try:
        here = os.path.dirname(os.path.abspath(__file__))
        root_key = os.path.abspath(os.path.join(here, os.pardir, "openAI_key.txt"))
        if os.path.exists(root_key):
            with open(root_key, "r") as f:
                return f.read().strip()
    except Exception:
        pass
    return None


# --- Token usage helpers ---
def _init_token_usage_state():
    if 'token_usage' not in st.session_state:
        st.session_state['token_usage'] = {
            'calls': [],
            'totals': {'prompt_tokens': 0, 'completion_tokens': 0, 'total_tokens': 0},
        }


def _reset_token_usage():
    st.session_state['token_usage'] = {
        'calls': [],
        'totals': {'prompt_tokens': 0, 'completion_tokens': 0, 'total_tokens': 0},
    }


def _add_token_usage(phase: str, model: str, usage: dict):
    _init_token_usage_state()
    usage = usage or {}
    pt = int(usage.get('prompt_tokens', 0) or 0)
    ct = int(usage.get('completion_tokens', 0) or 0)
    tt = int(usage.get('total_tokens', pt + ct) or (pt + ct))
    st.session_state['token_usage']['calls'].append({
        'phase': phase, 'model': model, 'prompt_tokens': pt, 'completion_tokens': ct, 'total_tokens': tt,
    })
    st.session_state['token_usage']['totals']['prompt_tokens'] += pt
    st.session_state['token_usage']['totals']['completion_tokens'] += ct
    st.session_state['token_usage']['totals']['total_tokens'] += tt


def _show_last_and_total_tokens():
    tu = st.session_state.get('token_usage')
    if not tu or not tu.get('calls'):
        return
    last = tu['calls'][-1]
    totals = tu['totals']
    st.info(
        f"토큰 사용량 — 이번 호출({last['phase']}, {last['model']}): "
        f"prompt {last['prompt_tokens']}, completion {last['completion_tokens']}, total {last['total_tokens']} | "
        f"누적: prompt {totals['prompt_tokens']}, completion {totals['completion_tokens']}, total {totals['total_tokens']}"
    )

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

        # 년도(왼쪽) - 보기축/필터(오른쪽)
        col_year, col_right = st.columns([1, 2])


        with col_year:
            selected_year = st.selectbox("📅 데이터셋 년도", sorted(available_years_all), index=0)


        available_panels = get_available_panel_ids(selected_year, data_root)
        available_targets = get_available_targets(selected_year, data_root)
        view_options = []
        if available_panels:
            view_options.append("세포 패널 보기")
        if available_targets:
            view_options.append("타깃 보기")
        if not view_options:
            view_options = ["타깃 보기"]

        selected_view = st.radio("보기 축", view_options, index=0, horizontal=True)

        selected_panel = None
        selected_target = None

        with col_right:
            if selected_view == "세포 패널 보기" and available_panels:
                panel_display_options = ["전체 패널"]
                panel_id_to_display = {"전체 패널": None}
                for panel_id in sorted(available_panels):
                    display_name = panel_names_map.get(panel_id, panel_id)
                    display_option = f"{panel_id} ({display_name})" if str(display_name).strip() != str(panel_id).strip() else f"{display_name}"
                    panel_id_to_display[display_option] = panel_id
                selected_panel_display = st.selectbox("🧬 패널 선택", panel_display_options, index=0)
                selected_panel = panel_id_to_display[selected_panel_display]
            else:
                targets = available_targets
                target_display_options = ["전체 타깃"] + targets
                selected_target_display = st.selectbox("🎯 타깃 선택", target_display_options, index=0)
                selected_target = None if selected_target_display == "전체 타깃" else selected_target_display


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
                    "2017: 국립암센터·한국화학연 특허(KR101920163B1, PCT/WO2018021849A1) 부속 표를 정리한 세트입니다. 암종별 세포주 패널에서 세포독성 IC₅₀(μM)이 보고되며, 일부 표엔 도스·치사율·심박변화 같은 보조 정보가 포함됩니다.\n\n"
                    "활용 가이드\n"
                    "- 세포 패널 보기에서 원하는 패널/세포주를 선택해 동일 맥락에서 구조–활성 비교를 진행하세요. pAct 모드(-log10[M])를 권장합니다.\n"
                    "- Activity Cliff: 유사도 임계값(≥0.85)과 ΔpAct 임계값을 조정하며 선택성 높은 쌍을 선별하세요.\n"
                    "- 보조 정보(도스·치사율·심박변화)는 후보 우선순위 판단 시 안전성/독성 리스크 관점에서 함께 참고하세요."
                )
            elif str(selected_year) == "2018":
                _desc = (
                    "2018: PCT/EP2018/056824, WO 2018/172250 부속 표를 정리한 세트입니다. Ras–SOS1 상호작용(HTRF, Assay 1~3)과 EGFR kinase 결과가 포함되며, 값은 IC₅₀ 또는 20 µM 단일 농도의 % 억제로 보고됩니다.\n\n"
                    "활용 가이드\n"
                    "- Ras–SOS1 관련 지표를 선택해 계열 내 SAR을 확인하세요. pAct 모드(-log10[M])를 권장합니다.\n"
                    "- EGFR 선택성: EGFR 값을 함께 확인해 오프타깃 민감 후보를 빠르게 가려내세요.\n"
                    "- 20 µM 단일농도 % 억제는 보조 지표로만 참고하고, IC₅₀와 함께 해석하세요."
                )
            elif str(selected_year) == "2020":
                _desc = (
                    "2020: 특허 WO2020132269 부속 표를 정리한 USP1 저해제 세트입니다. Examples(구조)·Intermediates/BB(부품)와 함께, 효소 저해 IC₅₀가 등급 기호(+, ++, +++, ++++)로 보고됩니다. 앱에서 등급→수치로 변환되어 비교할 수 있습니다.\n"
                    "· 민감도(Table 5): <316 nM 기준 Yes/No 라벨(패널/유전형 수준).\n"
                    "· ADME(Table 6~7): 용해도(pH 2.0/7.4)·간 미소체 안정성(HLM/RLM t₁/₂)을 +/<임계·++/≥임계(10 μM, 25 min) 등급으로 요약.\n"
                    "· Asterisk 예외: 이미지–SMILES 불일치 교정표의 NEW SMILES가 적용됩니다.\n\n"
                    "활용 가이드\n"
                    "- 타깃을 ‘usp1’로 두고 Activity Cliff 분석에서 pAct 모드로 비교하세요(등급→수치 변환 반영).\n"
                    "- 후보 우선순위: 상위 등급(저농도) 중심으로 ΔpAct를 확인해 구조–활성 차이를 해석하세요.\n"
                    "- 민감도/ADME 요약은 품질/적합성의 보조 기준으로 참고하되, 직접 비교에 과신하지 마세요."
                )
            elif str(selected_year) == "2021":
                _desc = (
                    "2021: 특허 WO2021163344 부속 표를 정리한 PRMT5 저해제 세트입니다. Table 19에 생화학 Ki(nM, Method A/B)와 세포 증식 IC₅₀(μM, HCT116 MTAP-null/WT)이 수치로 보고되며 일부는 ‘>10’ 같은 검열 표기가 있습니다. 값과 단위는 앱에서 일관 표기로 표시됩니다.\n\n"
                    "활용 가이드\n"
                    "- 효소–세포 상관: 타깃 보기(prmt5, Ki A/B)와 세포 패널 보기(HCT116: MTAP-null/WT)를 각각 선택해 상관/격차를 확인하세요.\n"
                    "- 선택성 지표: MTAP-null vs WT의 ΔpIC₅₀ 또는 (WT/null) 비율로 선택성을 계산해 우선 후보를 추립니다.\n"
                    "- Method 우선순위: A를 우선 사용(결측 시 B). ‘>10’ 검열값은 하한으로만 해석해 비교 시 주의하세요."
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
                            panel_id=(selected_panel if selected_view == "세포 패널 보기" else None),
                            cell_line=(selected_cell_line if selected_view == "세포 패널 보기" else None),
                            target_id=(selected_target if selected_view == "타깃 보기" else None)
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


with tab2:
    st.header("2. Activity Cliff 분석")
    if 'df' in st.session_state and st.session_state['df'] is not None:
        df = st.session_state['df']
        
        # (표시 제거) 단위 분포 요약

        # 컬럼 자동 선택(고정): Gold 스키마 가정하에 자동 결정
        smiles_col = 'SMILES' if 'SMILES' in df.columns else df.columns[0]
        activity_col = 'Activity' if 'Activity' in df.columns else ( 'value_std' if 'value_std' in df.columns else (df.columns[1] if len(df.columns) > 1 else df.columns[0]))
        # (표시 제거) 자동 선택된 컬럼 안내

        col1, col2 = st.columns(2)
        with col1:
            scale_choice = st.selectbox("활성도 스케일", ["원본(단위 유지)", "pAct (-log10[M])"], index=1)
            # 스케일에 따라 활성도 의미 자동 설정(UI 비노출)
            if scale_choice == "pAct (-log10[M])":
                st.session_state['activity_assumption'] = '값이 높을수록 활성도가 높음 (Higher is better)'
                st.caption("활성도 의미: pAct 스케일 → 값이 높을수록 활성도가 높음")
            else:
                st.session_state['activity_assumption'] = '값이 낮을수록 활성도가 높음 (Lower is better)'
                st.caption("활성도 의미: 원본(IC50 등) 스케일 → 값이 낮을수록 활성도가 높음")

        with col2:
            similarity_threshold = st.slider("구조 유사도 임계값 (Tanimoto)", 0.7, 1.0, 0.85, 0.01)
            # 스케일별 Δ 기본값: pAct는 관례적으로 2.0 권장
            default_diff = 2.0 if scale_choice == "pAct (-log10[M])" else 1.0
            step = 0.1 if scale_choice == "pAct (-log10[M])" else 0.1
            fmt = "%0.2f" if scale_choice == "pAct (-log10[M])" else "%f"

        # (UI 숨김) 활성도 의미는 위 스케일 선택에 의해 자동 설정됨

        # --- 프리뷰: 유사도만 적용한 분포 시각화 및 Δ 임계값 선택 ---
        # 스케일 변환 적용(프리뷰와 분석 모두에서 재사용)
        work_df = df.copy()
        work_df = work_df.dropna(subset=[activity_col]).reset_index(drop=True)
        chosen_activity_col = activity_col
        if scale_choice == "pAct (-log10[M])":
            unit_col = "unit_std" if "unit_std" in work_df.columns else None
            def _to_m(val, unit):
                try:
                    v = float(val)
                except Exception:
                    return None
                # 단위 정규화(소문자, µ→u 치환)
                u = (str(unit) if unit is not None else "").strip()
                u_norm = u.replace("µ", "u").strip().lower()
                if u_norm == "m":
                    return v
                if u_norm == "mm":
                    return v * 1e-3
                if u_norm == "um":
                    return v * 1e-6
                if u_norm == "nm":
                    return v * 1e-9
                if u_norm == "pm":
                    return v * 1e-12
                # % 등 비농도 단위는 pAct 계산 대상에서 제외
                if u_norm in {"%", "percent", "pct"}:
                    return None
                # 알 수 없는 경우 보수적으로 uM 가정(호환성 유지)
                return v * 1e-6
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

        # 프리뷰 쌍(Δ 필터 미적용) — 동일 파라미터에서 재사용되도록 세션 캐시
        _cache_key = (
            "pairs_preview",
            int(similarity_threshold * 1000),
            chosen_activity_col,
            len(work_df),
        )
        prev_cache = st.session_state.get("_preview_cache", {})
        preview_df = None
        if prev_cache.get("key") == _cache_key:
            preview_df = prev_cache.get("df")
        if preview_df is None:
            preview_df = find_activity_cliffs(
                work_df,
                smiles_col=smiles_col,
                activity_col=chosen_activity_col,
                similarity_threshold=similarity_threshold,
                activity_diff_threshold=0.0,
                higher_is_better=(st.session_state.activity_assumption == '값이 높을수록 활성도가 높음 (Higher is better)')
            )
            st.session_state["_preview_cache"] = {"key": _cache_key, "df": preview_df}

        # Δ 임계값 입력(히트맵 보기 전에 기본값 산정)
        try:
            import numpy as _np
            q95 = float(preview_df["Activity_Diff"].quantile(0.95)) if len(preview_df) else default_diff
            max_for_slider = max(0.1, q95)
        except Exception:
            max_for_slider = default_diff
        activity_diff_threshold = st.slider(
            "활성도 차이 임계값",
            min_value=0.0,
            max_value=float(max_for_slider),
            value=float(min(default_diff, max_for_slider)),
            step=float(step)
        )

        # --- 단일 히트맵 영역(프리뷰/결과 토글) ---
        try:
            import altair as alt
            import pandas as _pd
            # 결과가 이미 생성되어 있으면 토글 제공
            has_result = ('cliff_df' in st.session_state) and (st.session_state['cliff_df'] is not None) and (len(st.session_state['cliff_df']) > 0)
            options = ["프리뷰(유사도만 적용)"] + (["결과(임계값 모두 적용)"] if has_result else [])
            src_choice = st.radio("시각화 소스", options, horizontal=True, index=0, key="viz_source_choice_main")
            viz_df = preview_df if src_choice.startswith("프리뷰") else st.session_state['cliff_df']

            total_pairs = len(viz_df)
            mask_sel = (viz_df["Similarity"] >= similarity_threshold) & (viz_df["Activity_Diff"] >= activity_diff_threshold)
            selected_pairs = int(mask_sel.sum())
            ratio = (selected_pairs / total_pairs * 100.0) if total_pairs else 0.0
            m1, m2, m3 = st.columns(3)
            m1.metric("전체 쌍 수", f"{total_pairs:,}")
            m2.metric("선택 영역 쌍 수", f"{selected_pairs:,}")
            m3.metric("비율(%)", f"{ratio:0.1f}")

            base = alt.Chart(viz_df)
            heat = base.mark_rect().encode(
                alt.X("Similarity:Q", bin=alt.Bin(maxbins=30), scale=alt.Scale(domain=[0.7, 1.0])),
                alt.Y("Activity_Diff:Q", bin=alt.Bin(maxbins=30)),
                alt.Color("count():Q", scale=alt.Scale(scheme="magma")),
                tooltip=[alt.Tooltip("count():Q", title="Count")]
            ).properties(height=320)
            v_rule = alt.Chart(_pd.DataFrame({"x": [similarity_threshold]})).mark_rule(color="#00FFFF", strokeDash=[6,4]).encode(x="x:Q")
            h_rule = alt.Chart(_pd.DataFrame({"y": [activity_diff_threshold]})).mark_rule(color="#00FFFF", strokeDash=[6,4]).encode(y="y:Q")
            st.altair_chart((heat + v_rule + h_rule).resolve_scale(color="independent"), use_container_width=True)
        except Exception:
            import numpy as _np
            st.info("시각화 엔진(Altair) 사용이 어려워 간단한 산점도로 대체합니다.")
            src = preview_df
            if ('cliff_df' in st.session_state) and len(st.session_state['cliff_df']) > 0:
                src = st.session_state['cliff_df']
            sample = src.sample(min(len(src), 5000), random_state=42) if len(src) > 5000 else src
            st.scatter_chart(sample[["Similarity", "Activity_Diff"]])

        if st.button("Activity Cliff 분석 실행"):
            with st.spinner("Activity Cliff를 분석 중입니다..."):
                # work_df / chosen_activity_col는 프리뷰 단계에서 이미 계산됨

                cliff_df = find_activity_cliffs(
                    work_df,
                    smiles_col=smiles_col,
                    activity_col=chosen_activity_col,
                    similarity_threshold=similarity_threshold,
                    activity_diff_threshold=activity_diff_threshold,
                    higher_is_better= (st.session_state.activity_assumption == '값이 높을수록 활성도가 높음 (Higher is better)')
                )
            
            st.success(f"{len(cliff_df)}개의 Activity Cliff 쌍을 찾았습니다!")
            # IUPAC 매핑(표시용)
            try:
                base_df = st.session_state.get('df', pd.DataFrame())
                name_map = {}
                if len(base_df) > 0 and 'SMILES' in base_df.columns:
                    if 'iupac_name' in base_df.columns:
                        tmp = base_df[['SMILES','iupac_name']].drop_duplicates()
                        name_map = {str(r['SMILES']): (str(r['iupac_name']).strip() if str(r['iupac_name']).strip() else None) for _, r in tmp.iterrows()}
                def _name_or_smiles(s):
                    s = str(s)
                    n = name_map.get(s)
                    return n if n and n.lower() != 'nan' else s
                display_df = cliff_df.copy()
                display_df['IUPAC_1'] = display_df['SMILES_1'].map(_name_or_smiles)
                display_df['IUPAC_2'] = display_df['SMILES_2'].map(_name_or_smiles)
                # 표시 컬럼 재배치 및 정렬(Δ, sim 내림차순)
                cols = ['IUPAC_1','Activity_1','IUPAC_2','Activity_2','Similarity','Activity_Diff','SMILES_1','SMILES_2']
                cols = [c for c in cols if c in display_df.columns]
                display_df = display_df[cols].sort_values(by=['Activity_Diff','Similarity'], ascending=[False, False]).reset_index(drop=True)
                # 1-based 인덱스
                display_df.index = range(1, len(display_df) + 1)
                display_df.index.name = 'Pair #'
                st.dataframe(display_df, use_container_width=True)
                try:
                    csv_bytes = display_df.to_csv(index=True).encode('utf-8')
                    st.download_button(
                        label="결과 CSV 다운로드",
                        data=csv_bytes,
                        file_name="activity_cliffs.csv",
                        mime="text/csv",
                        use_container_width=True,
                    )
                except Exception:
                    pass
            except Exception:
                # 문제가 있으면 원본 표라도 제공
                df0 = cliff_df.copy().sort_values(by=['Activity_Diff','Similarity'], ascending=[False, False]).reset_index(drop=True)
                df0.index = range(1, len(df0) + 1)
                df0.index.name = 'Pair #'
                st.dataframe(df0, use_container_width=True)

            st.session_state['cliff_df'] = cliff_df

            # 활성도 지표에 대한 요약 정보 표시
            summary_text = create_activity_summary(activity_col, (st.session_state.activity_assumption == '값이 높을수록 활성도가 높음 (Higher is better)'))
            st.markdown("---")
            st.markdown(summary_text)

            # (히트맵은 상단의 단일 영역에서 이미 제공)
    else:
        st.info("1. 상단에서 Gold 데이터를 먼저 로드해주세요.")

with tab3:
    st.header("3. 결과 시각화 및 가설 생성")
    if 'cliff_df' in st.session_state and not st.session_state['cliff_df'].empty:
        cliff_df = st.session_state['cliff_df']

        # IUPAC 매핑 및 정렬(표시용)
        base_df = st.session_state.get('df', pd.DataFrame())
        name_map = {}
        if len(base_df) > 0 and 'SMILES' in base_df.columns and 'iupac_name' in base_df.columns:
            tmp = base_df[['SMILES','iupac_name']].drop_duplicates()
            name_map = {str(r['SMILES']): (str(r['iupac_name']).strip() if str(r['iupac_name']).strip() else None) for _, r in tmp.iterrows()}
        def _name_or_smiles(s):
            s = str(s)
            n = name_map.get(s)
            return n if n and n.lower() != 'nan' else s
        disp = cliff_df.copy()
        disp['IUPAC_1'] = disp['SMILES_1'].map(_name_or_smiles)
        disp['IUPAC_2'] = disp['SMILES_2'].map(_name_or_smiles)
        sort_idx = disp.sort_values(by=['Activity_Diff','Similarity'], ascending=[False, False]).index.tolist()
        st.session_state['cliff_sorted_idx'] = sort_idx
        disp = disp.loc[sort_idx]
        cols = ['IUPAC_1','Activity_1','IUPAC_2','Activity_2','Similarity','Activity_Diff','SMILES_1','SMILES_2']
        cols = [c for c in cols if c in disp.columns]
        disp = disp[cols].reset_index(drop=True)
        disp.index = range(1, len(disp)+1)
        disp.index.name = 'Pair #'

        st.subheader("분석할 Activity Cliff 쌍 선택")
        # 표에서 직접 체크박스로 선택하도록 지원
        try:
            # 표시용 테이블에 선택 컬럼 추가
            disp_show = disp.copy()
            disp_show.insert(0, 'Pair #', list(range(1, len(disp_show)+1)))
            default_selected = set(st.session_state.get('selected_display_ids', []))
            disp_show.insert(0, '선택', [ (i in default_selected) for i in disp_show['Pair #'] ])
            edited = st.data_editor(
                disp_show,
                hide_index=True,
                use_container_width=True,
                disabled=[c for c in disp_show.columns if c not in ['선택']],
                key='cliff_table_editor'
            )
            selected_display_ids = [ int(v) for v in edited[edited['선택'] == True]['Pair #'].tolist() ]
            st.session_state['selected_display_ids'] = selected_display_ids
        except Exception:
            # 폴백: 멀티셀렉트 유지
            st.dataframe(disp, use_container_width=True)
            options_1based = list(range(1, len(disp)+1))
            selected_display_ids = st.multiselect("분석 및 시각화할 Pair 번호(1-based)를 선택하세요:", options_1based)
        # 1-based → 원본 인덱스 매핑
        selected_indices = [ st.session_state['cliff_sorted_idx'][i-1] for i in selected_display_ids ]
        
        if selected_indices:
            # 선택 쌍 미리보기(가설 생성 전 즉시 표시)
            with st.container(border=True):
                max_preview = 6
                pairs = list(zip(selected_display_ids, selected_indices))
                if len(pairs) > max_preview:
                    st.info(f"미리보기는 최대 {max_preview}쌍까지만 표시합니다. (총 선택 {len(pairs)}쌍)")
                for disp_id, i in pairs[:max_preview]:
                    row = cliff_df.loc[i]
                    nm1 = name_map.get(str(row['SMILES_1'])) if name_map else None
                    nm2 = name_map.get(str(row['SMILES_2'])) if name_map else None
                    l1 = f"IUPAC: {nm1}" if nm1 else f"SMILES: {row['SMILES_1']}"
                    l2 = f"IUPAC: {nm2}" if nm2 else f"SMILES: {row['SMILES_2']}"
                    st.markdown(f"### 분석 쌍 #{disp_id}")
                    img = visualize_structure_difference(
                        smiles1=row['SMILES_1'],
                        smiles2=row['SMILES_2'],
                        legend1=f"{l1}\nActivity: {row['Activity_1']:.2f}",
                        legend2=f"{l2}\nActivity: {row['Activity_2']:.2f}"
                    )
                    st.image(img, caption=f"유사도: {row['Similarity']:.3f} | 활성도 차이: {row['Activity_Diff']:.2f}")

            openai_api_key = get_openai_api_key_from_file()
            if not openai_api_key:
                st.warning("LLM 가설 생성을 위해 openAI_key.txt 파일에 API Key를 입력해주세요.")
            else:
                st.session_state['openai_api_key'] = openai_api_key

                if st.button("선택된 쌍에 대한 가설 생성"):
                    output_dir = "hypotheses"
                    os.makedirs(output_dir, exist_ok=True)

                    for i in selected_indices:
                        row = cliff_df.loc[i]
                        st.subheader(f"가설 생성 중 · 쌍 #{i}")

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

                            # 컨텍스트 JSON 구성
                            try:
                                selected_axis = 'cell' if selected_view == '세포 패널 보기' else 'target'
                                scale_used = 'pAct' if (st.session_state.get('activity_assumption') == '값이 높을수록 활성도가 높음 (Higher is better)') else 'raw'
                                ctx = build_pair_context(
                                    year=str(selected_year),
                                    pair_row=row,
                                    df_base=st.session_state['df'],
                                    data_root=data_root,
                                    selected_axis=selected_axis,
                                    scale_used=scale_used,
                                )
                            except Exception:
                                ctx = None

                            gen_result = generate_hypothesis(
                                api_key=openai_api_key,
                                smiles1=low_act_smiles,
                                activity1=low_act_val,
                                smiles2=high_act_smiles,
                                activity2=high_act_val,
                                structural_difference_description=f"화합물 1({low_act_smiles})과 화합물 2({high_act_smiles})의 구조적 차이점.",
                                similarity=row['Similarity'],
                                context_json=ctx
                            )
                            _add_token_usage('generation', gen_result.get('model', 'unknown'), gen_result.get('usage', {}))
                            _show_last_and_total_tokens()
                            content = gen_result.get('content', '')

                            try:
                                hypothesis_data = json.loads(content)
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
                                # LLM 오류 메시지 처리
                                msg = str(content or "").strip()
                                rate_hit = ("rate" in msg.lower()) or ("429" in msg) or ("insufficient_quota" in msg) or ("quota" in msg.lower())
                                if rate_hit:
                                    st.error("LLM 호출이 제한되었습니다(리밋/한도). Billing/Usage를 확인하세요.")
                                else:
                                    st.error("LLM 응답이 유효한 JSON 형식이 아닙니다. 잠시 후 다시 시도하세요.")
                                # 가설 저장/표시는 생략(생성하지 않음)
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

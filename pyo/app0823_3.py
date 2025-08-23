# -*- coding: utf-8 -*-
from __future__ import annotations

# =========================
# MSAR Streamlit UI (app.py)
# =========================
import os, sys, json, glob, time, subprocess, tempfile, textwrap
from typing import List, Optional, Dict, Tuple
import pathlib
import yaml
import pandas as pd
import numpy as np
import streamlit as st

# ------ Optional deps (RDKit / 3D) ------
RDKit_OK, RDKit_ERR = True, ""
try:
    from rdkit import Chem
    from rdkit.Chem import Draw, AllChem, Crippen, Lipinski, Descriptors, rdMolDescriptors
except Exception as e:
    RDKit_OK, RDKit_ERR = False, str(e)

P3D_OK = True
try:
    import py3Dmol  # conda-forge: py3Dmol
except Exception:
    P3D_OK = False

# ------ Paths ------
ROOT = pathlib.Path(__file__).resolve().parents[1]
DATA = ROOT / "data"
DATA_PROCESSED = DATA / "processed"
REPORTS = ROOT / "reports"
LOGS = REPORTS / "logs"
CONF_FILE = ROOT / "configs" / "ac_thresholds.yaml"
HYP_DIR = REPORTS / "explanations"
TH_HIST = REPORTS / "thresholds_history.csv"

# ------ Look & feel (mild/eye-friendly buttons) ------
st.set_page_config(page_title="MSAR : Molecular Structure Activity Reasoning", layout="wide")
st.markdown("""
<style>
/* softer primary color for buttons */
div.stButton > button[kind="primary"] {
  background: #0ea5a5; border: 0; color: white; font-weight: 600;
}
div.stButton > button[kind="secondary"] {
  background: #334155; color: #e2e8f0; border: 0; font-weight: 600;
}
div.stButton > button {
  border-radius: 8px;
}
</style>
""", unsafe_allow_html=True)

# ------ small helpers ------
def _py() -> str:
    return sys.executable

def _ts() -> str:
    return time.strftime("%Y-%m-%d %H:%M:%S")

def load_yaml(path: pathlib.Path) -> dict:
    if not path.exists():
        return {}
    with open(path, "r", encoding="utf-8") as f:
        return yaml.safe_load(f) or {}

def save_yaml(path: pathlib.Path, data: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", encoding="utf-8") as f:
        yaml.safe_dump(data, f, sort_keys=False, allow_unicode=True)

def list_prefixes() -> List[str]:
    # candidates like data/processed/<prefix>_clean.csv
    cands = []
    for p in DATA_PROCESSED.glob("*_clean.csv"):
        cands.append(p.name.split("_")[0])
    return sorted(set(cands))

def run_step(cmd: List[str], log_area, title: str) -> bool:
    """Run a python step and stream last ~200 lines into the UI."""
    try:
        proc = subprocess.Popen(cmd, cwd=str(ROOT), stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
        lines: List[str] = []
        while True:
            line = proc.stdout.readline()
            if not line and proc.poll() is not None:
                break
            if line:
                lines.append(line.rstrip())
                log_area.code("\n".join(lines[-200:]))
        rc = proc.wait()
        if rc != 0:
            st.error(f"❌ {title} 실패 (exit={rc})")
            return False
        st.success(f"✅ {title} 완료")
        return True
    except Exception as e:
        st.exception(e)
        return False

def hyp_file(prefix: str) -> pathlib.Path:
    return HYP_DIR / f"{prefix}_hypotheses.json"

def load_hyp(prefix: str) -> list[dict]:
    p = hyp_file(prefix)
    if not p.exists():
        return []
    try:
        return json.loads(p.read_text(encoding="utf-8"))
    except Exception:
        return []

def save_hyp(prefix: str, items: list[dict]) -> None:
    p = hyp_file(prefix)
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text(json.dumps(items, ensure_ascii=False, indent=2), encoding="utf-8")

# ------ RDKit/3D helpers ------
def mol_from_smiles(smi: str):
    if not RDKit_OK:
        return None
    try:
        m = Chem.MolFromSmiles(smi)
        return m
    except Exception:
        return None

def mol2d_png(mi, size=(250, 250)):
    if not RDKit_OK or mi is None:
        return None
    try:
        return Draw.MolToImage(mi, size=size)
    except Exception:
        return None

def mol3d_view(mi):
    """Embed simple 3D conformer & render in py3Dmol (if possible)."""
    if not (RDKit_OK and P3D_OK and mi is not None):
        return None
    try:
        m3d = Chem.AddHs(Chem.Mol(mi))
        AllChem.EmbedMolecule(m3d, AllChem.ETKDG())
        AllChem.UFFOptimizeMolecule(m3d, maxIters=200)
        mb = Chem.MolToMolBlock(m3d)
        v = py3Dmol.view(width=350, height=300)
        v.addModel(mb, 'mol')
        v.setStyle({'stick': {}})
        v.zoomTo()
        return v
    except Exception:
        return None

# ------ side: prefix + thresholds ------
def sidebar_controls() -> dict:
    st.sidebar.header("데이터 선택")
    prefixes = list_prefixes()
    if not prefixes:
        st.sidebar.warning("data/processed 에 *_clean.csv 가 없습니다.")
    prefix = st.sidebar.selectbox("Prefix (data/processed의 접두어)", prefixes, index=0 if prefixes else None)

    st.sidebar.markdown("### 임계값 (configs/ac_thresholds.yaml)")
    conf = load_yaml(CONF_FILE) or {}
    sim_default = float(conf.get("sim_threshold", 0.80))
    fold_default = float(conf.get("fold_threshold", 10.0))
    topk_default = int(conf.get("top_k", 200))

    sim = st.sidebar.slider("Tanimoto (sim_threshold)", 0.50, 1.00, value=sim_default, step=0.01)
    fold = st.sidebar.slider("Fold change (fold_threshold)", 1.0, 200.0, value=fold_default, step=0.5)
    top_k = st.sidebar.number_input("Top K pairs (top_k)", min_value=10, max_value=2000, value=topk_default, step=10)

    manual = st.sidebar.checkbox("직접 입력(정확한 수치 넣어쓰기)")
    if manual:
        sim = st.sidebar.number_input("Tanimoto (+수동)", min_value=0.0, max_value=1.0, value=float(sim), step=0.01, key="sim_manual")
        fold = st.sidebar.number_input("Fold change (+수동)", min_value=0.0, max_value=1e6, value=float(fold), step=0.5, key="fold_manual")

    # save thresholds
    save_col = st.sidebar.container()
    with save_col:
        if st.button("💾 임계값 저장", type="secondary", use_container_width=True):
            new_conf = dict(sim_threshold=float(sim), fold_threshold=float(fold), top_k=int(top_k))
            # keep optional grids if they exist
            for k in ("tanimoto_grid", "fold_grid"):
                if k in conf:
                    new_conf[k] = conf[k]
            save_yaml(CONF_FILE, new_conf)
            # history
            TH_HIST.parent.mkdir(parents=True, exist_ok=True)
            pd.DataFrame([{
                "ts": _ts(), "sim_threshold": float(sim), "fold_threshold": float(fold), "top_k": int(top_k)
            }]).to_csv(TH_HIST, index=False, mode=("a" if TH_HIST.exists() else "w"))
            st.sidebar.success(f"saved: {CONF_FILE}")

    st.sidebar.markdown("### 실행 단계")
    do_p1 = st.sidebar.checkbox("Phase 1 (EDA)", value=True)
    do_p2 = st.sidebar.checkbox("Phase 2 (Activity Cliffs)", value=True)

    # Phase 3 block + indent
    st.sidebar.markdown("**Phase 3 (LLM 해석)**")
    do_p3 = st.sidebar.checkbox("  LLM 실행 포함", value=False)
    st.sidebar.caption("  Phase3는 `OPENAI_API_KEY` 환경 변수 필요.")
    top_n = st.sidebar.slider("  LLM 상위 N개 (Phase3)", 5, 300, 30, 5)

    do_p4 = st.sidebar.checkbox("Phase 4 (리포트 빌드)", value=True)

    with st.sidebar.expander("🧭 도움말", expanded=False):
        st.markdown("""
- **임계값 저장**을 누르면 YAML과 로그 CSV(변경 이력)에 기록됩니다.  
- Phase 2를 여러 번 사용해도 단계별 로그는 보고서와 Streamlit에 표시됩니다.  
- Phase 3: `OPENAI_API_KEY` 환경변수가 없으면 실행되지 않습니다.  
  - PowerShell(현재 세션): `$env:OPENAI_API_KEY="sk-xxxx"`  
  - 사용자 고정: `[Environment]::SetEnvironmentVariable("OPENAI_API_KEY","sk-xxxx","User")`  
""")

    # run buttons
    run = st.sidebar.button("▶ 실행", type="primary", use_container_width=True)
    return dict(prefix=prefix, sim=sim, fold=fold, top_k=top_k, do_p1=do_p1, do_p2=do_p2, do_p3=do_p3, do_p4=do_p4, top_n=top_n, run=run)

# ------ data/summary tab ------
def ui_data_summary(prefix: str):
    st.subheader("📊 데이터/요약")
    sj = REPORTS / "phase1" / f"summary_{prefix}.json"
    if not sj.exists():
        st.info("요약 파일이 아직 없습니다. Phase1을 먼저 실행하세요.")
        return
    s = json.loads(sj.read_text(encoding="utf-8"))

    # headline metrics
    c1, c2, c3 = st.columns(3)
    c1.metric("N compounds", s.get("n_compounds", "N/A"))
    c2.metric("SOS1_best_nM (median)", f"{s.get('activity_stats',{}).get('SOS1_best_nM',{}).get('p50','N/A')}")
    c3.metric("selectivity_lb (median)", f"{s.get('activity_stats',{}).get('selectivity_lb',{}).get('p50','N/A')}")

    # activity_stats
    act = s.get("activity_stats", {})
    if act:
        rows = []
        for k in ["SOS1_best_nM", "EGFR kinase inhibition IC50_nM", "selectivity_lb"]:
            v = act.get(k, {})
            rows.append(dict(metric=k, count=v.get("count"), mean=v.get("mean"), std=v.get("std"), min=v.get("min"),
                             q25=v.get("p25"), q50=v.get("p50"), q75=v.get("p75"), max=v.get("max")))
        st.markdown("#### activity_stats")
        st.dataframe(pd.DataFrame(rows), use_container_width=True)
    # physchem
    pc = s.get("physchem_stats", {})
    if pc:
        rows = []
        order = ["MW", "TPSA", "HBD", "HBA", "RotB", "Rings", "AromaticRings", "FractionCSP3"]
        for k in order:
            v = pc.get(k, {})
            if v:
                rows.append(dict(metric=k, count=v.get("count"), mean=v.get("mean"), std=v.get("std"), min=v.get("min"),
                                 q25=v.get("p25"), q50=v.get("p50"), q75=v.get("p75"), max=v.get("max")))
        st.markdown("#### physchem_stats (샘플)")
        st.dataframe(pd.DataFrame(rows), use_container_width=True)

    with st.expander("📑 리포트/다운로드", expanded=False):
        pdf = REPORTS / "_report_all.pdf"
        md = REPORTS / f"auto_report_{prefix}.md"
        if pdf.exists():
            with open(pdf, "rb") as f:
                st.download_button("📎 _report_all.pdf 다운로드", f, file_name="_report_all.pdf", mime="application/pdf")
        if md.exists():
            st.download_button("📎 auto_report.md 다운로드", md.read_bytes(), file_name=md.name, mime="text/markdown")

    with st.expander("🖼️ 그림 미리보기", expanded=False):
        # 단순히 파일이 있으면 PNG 4개 정도 보여주기
        figs = list((REPORTS / "figures").glob("*.png"))
        if not figs:
            st.info("예시 그래프를 표시할 데이터가 없습니다.")
        else:
            for p in figs[:4]:
                st.image(str(p), caption=p.name, use_column_width=True)

# ------ pair inspector tab ------
def ui_pair_inspector(prefix: str):
    st.subheader("🔬 Pair inspector")
    path_exp = DATA_PROCESSED / f"{prefix}_cliffs_explained.csv"
    path_raw = DATA_PROCESSED / f"{prefix}_cliffs.csv"
    if not path_exp.exists() and not path_raw.exists():
        st.info("cliffs 파일이 없습니다. Phase 2를 먼저 실행하세요.")
        return
    path = path_exp if path_exp.exists() else path_raw
    try:
        df = pd.read_csv(path)
    except Exception as e:
        st.warning(f"파일 로드 실패: {e}")
        return

    show_cols = [c for c in ["i","j","Num_i","Num_j","smiles_i","smiles_j","tanimoto","fold_change","explanation"] if c in df.columns]
    st.dataframe(df[show_cols].head(200), use_container_width=True)
    if df.empty:
        st.info("표시할 페어가 없습니다. 임계값을 낮추거나 Phase 2를 다시 실행하세요.")
        return

    # index slider (safe bounds)
    max_idx = max(0, len(df)-1)
    idx = st.number_input("행 번호 선택", min_value=0, max_value=max_idx, value=0, step=1)
    row = df.iloc[int(idx)]

    # text diff / key values
    st.json({k: row.get(k) for k in ["Num_i","Num_j","tanimoto","fold_change"] if k in row})

    # 2D & 3D drawings (best-effort)
    col2d, col3d = st.columns(2)
    smi_i, smi_j = str(row.get("smiles_i","")), str(row.get("smiles_j",""))

    if RDKit_OK:
        mi, mj = mol_from_smiles(smi_i), mol_from_smiles(smi_j)
        img = Draw.MolsToGridImage([mi, mj], legends=["A", "B"], useSVG=False, molsPerRow=2, subImgSize=(300,300))
        col2d.image(img, caption="2D structures (A vs B)")
        # 간단 하이라이트 설명
        col2d.caption("Note: 구조적 차이를 좌우하는 치환/작용기 변화를 중심으로 LLM 해석(Phase3)과 함께 보세요.")
        if P3D_OK:
            vA, vB = mol3d_view(mi), mol3d_view(mj)
            if vA is not None:
                col3d.components.v1.html(vA._make_html(), height=320)
            if vB is not None:
                col3d.components.v1.html(vB._make_html(), height=320)
    else:
        st.warning(f"RDKit 사용 불가: {RDKit_ERR}\n\nNumPy 2.x로 올라갔다면 아래를 권장합니다.\n"
                   "```conda install \"numpy<2\" rdkit=2024.03 -c conda-forge```")

# ------ hypotheses CRUD tab ------
def ui_hypotheses(prefix: str):
    st.subheader("🧠 가설")
    hyps = load_hyp(prefix)
    st.write(f"현재 저장된 가설: **{len(hyps)}** 건")
    with st.form("new_hyp"):
        pair_idx = st.number_input("pair index", 0, 10**9, 0, step=1)
        rationale = st.text_area("주요 가설(핵심 근거 요약)")
        mechanism = st.text_area("기전 가설")
        counter = st.text_area("반대가설/위험요인")
        admet = st.text_area("ADMET 위험성")
        design = st.text_area("검증을 위한 분자설계 아이디어")
        priority = st.slider("우선순위 (1~5)", 1, 5, 3)
        status = st.selectbox("상태", ["draft", "accepted", "reject"])
        submitted = st.form_submit_button("➕ 가설 추가")
    if submitted:
        hyps.append({"id": int(time.time()), "pair_idx": int(pair_idx), "rationale": rationale,
                     "mechanism": mechanism, "counter": counter, "admet": admet,
                     "design": design, "priority": int(priority), "status": status, "score": None})
        save_hyp(prefix, hyps)
        st.success("추가 저장 완료")

    if hyps:
        # edit / delete
        sel = st.selectbox("선택 가설 (id)", [h["id"] for h in hyps])
        cur = next(h for h in hyps if h["id"] == sel)
        nstat = st.selectbox("상태 변경", ["draft","accepted","reject"], index=["draft","accepted","reject"].index(cur["status"]))
        npri = st.slider("우선순위 변경", 1, 5, cur["priority"])
        if st.button("💾 수정 저장"):
            cur["status"] = nstat; cur["priority"] = int(npri)
            save_hyp(prefix, hyps); st.success("업데이트 완료")
        if st.button("🗑️ 선택 가설 삭제"):
            hyps = [h for h in hyps if h["id"] != sel]
            save_hyp(prefix, hyps); st.warning("삭제 완료")

        st.download_button("📥 가설(전체) JSON 다운로드",
                           data=json.dumps(hyps, ensure_ascii=False, indent=2),
                           file_name=f"{prefix}_hypotheses.json",
                           mime="application/json")
        st.table(pd.DataFrame(hyps)[["id","pair_idx","priority","status","score","rationale"]])

# ------ hypothesis scoring tab ------
def ui_hypothesis_scoring(prefix: str):
    st.subheader("🧪 가설 평가(간단 · 0~100)")
    hyps = load_hyp(prefix)
    if not hyps:
        st.info("가설이 없습니다."); return

    path = DATA_PROCESSED / f"{prefix}_cliffs_explained.csv"
    if path.exists():
        df = pd.read_csv(path)
    else:
        df = None

    for h in hyps:
        sc = 0.0
        if df is not None and isinstance(h.get("pair_idx"), int) and 0 <= h["pair_idx"] < len(df):
            row = df.iloc[int(h["pair_idx"])]
            fc = float(row.get("fold_change", 1.0))
            tan = float(row.get("tanimoto", 0.0))
            sc = min(100.0, 20.0*min(fc, 10.0) + 80.0*max(0.0, tan-0.5))
        h["score"] = round(sc, 1)
    save_hyp(prefix, hyps)
    st.dataframe(pd.DataFrame(hyps)[["id","pair_idx","priority","status","score","rationale"]], use_container_width=True)

# ------ ADMET quick check tab ------
def ui_admet(prefix: str):
    st.subheader("🩺 ADMET 빠른 체크 (Lipinski/Veber)")
    clean = DATA_PROCESSED / f"{prefix}_clean.csv"
    if not clean.exists():
        st.info("clean.csv가 없습니다. Phase1을 먼저 실행하세요.")
        return
    df = pd.read_csv(clean)
    if "SMILES" in df.columns:
        smiles_col = "SMILES"
    elif "smiles" in df.columns:
        smiles_col = "smiles"
    else:
        st.info("SMILES 컬럼을 찾을 수 없습니다."); return

    if not RDKit_OK:
        st.warning("RDKit이 없어 ADMET 계산을 건너뜁니다."); return

    rows = []
    for smi in df[smiles_col].dropna().head(500).astype(str):
        m = mol_from_smiles(smi)
        if m is None: continue
        mw = Descriptors.MolWt(m)
        logp = Crippen.MolLogP(m)
        hbd = Lipinski.NumHDonors(m)
        hba = Lipinski.NumHAcceptors(m)
        tpsa = rdMolDescriptors.CalcTPSA(m)
        rotb = Lipinski.NumRotatableBonds(m)
        lipinski_viol = int(mw>500) + int(logp>5) + int(hbd>5) + int(hba>10)
        veber_ok = int((rotb<=10) and (tpsa<=140))
        rows.append([mw, logp, hbd, hba, tpsa, rotb, lipinski_viol, veber_ok])
    out = pd.DataFrame(rows, columns=["MW","logP","HBD","HBA","TPSA","RotB","LipinskiViol","VeberOK"])
    st.write("요약")
    st.dataframe(out.describe().T, use_container_width=True)
    st.write("룰 위반 카운트 (Lipinski)")
    st.bar_chart(out["LipinskiViol"].value_counts().sort_index())
    st.write("Veber 적합 비율")
    st.progress(out["VeberOK"].mean())

# ------ Threshold history tab ------
def ui_threshold_history():
    st.subheader("🧾 임계값 변경 이력")
    if not TH_HIST.exists():
        st.info("이력 파일이 없습니다.")
        return
    df = pd.read_csv(TH_HIST)
    st.dataframe(df.tail(50), use_container_width=True)

# ------ Sensitivity sweep tab ------
def ui_sweep(prefix: str):
    st.subheader("📐 민감도 스윕")
    conf = load_yaml(CONF_FILE)
    tg = conf.get("tanimoto_grid", [0.75, 0.80, 0.85])
    fg = conf.get("fold_grid", [5.0, 10.0, 15.0])
    st.caption(f"grid: tanimoto_grid={tg}, fold_grid={fg}")
    if st.button("▶ 민감도 스윕 실행", type="primary"):
        res = []
        log = st.empty()
        for sim in tg:
            for fold in fg:
                # 임시 YAML 만들고 phase2만 호출
                with tempfile.NamedTemporaryFile("w", suffix=".yaml", delete=False) as fp:
                    yaml.safe_dump({"sim_threshold": float(sim), "fold_threshold": float(fold),
                                    "top_k": int(conf.get("top_k", 200))}, fp)
                    tmp_conf = fp.name
                ok = run_step([_py(), "src/phase2_cliffs.py", "--prefix", prefix, "--config", tmp_conf],
                              log, f"Phase2(sim={sim}, fold={fold})")
                # pairs=0/파일없음도 수용
                path = DATA_PROCESSED / f"{prefix}_cliffs.csv"
                n = 0
                if path.exists():
                    try:
                        n = max(0, len(pd.read_csv(path)))
                    except Exception:
                        n = 0
                res.append({"sim": sim, "fold": fold, "pairs": n})
        st.dataframe(pd.DataFrame(res).pivot(index="sim", columns="fold", values="pairs").fillna(0).astype(int),
                     use_container_width=True)

# ------ Print controls (do not hide tabs) ------
def ui_print_controls():
    st.toggle("🖨️ 인쇄 모드", key="print_mode", help="선택 섹션만 아래에 인쇄용으로 재표시합니다.")
    sel = st.multiselect("인쇄할 섹션 선택", ["데이터/요약","Pair inspector","가설"], key="print_sel", help="탭은 그대로 유지됩니다.")
    return st.session_state.get("print_mode", False), sel

# =========================
# Main
# =========================
ctl = sidebar_controls()

# Title row + print controls on the right
tcol1, tcol2 = st.columns([0.72, 0.28])
with tcol1:
    st.title("MSAR : Molecular Structure Activity Reasoning")
with tcol2:
    print_mode, print_sel = ui_print_controls()

# RDKit import 문제가 있으면 상단 안내
if not RDKit_OK:
    with st.expander("⚠ RDKit 불러오기 실패 안내 (클릭)", expanded=True):
        st.error(f"{RDKit_ERR}")
        st.write("**권장 해결**")
        st.code('conda install "numpy<2" rdkit=2024.03 -c conda-forge', language="bash")

# Tabs (항상 보이게)
tab1, tab2, tab3, tab4, tab5, tab6, tab7 = st.tabs(["📊 데이터/요약", "🔬 Pair inspector", "🧠 가설", "🧪 가설 평가", "🩺 ADMET", "🧾 임계값 이력", "📐 민감도 스윕"])

with tab1: ui_data_summary(ctl["prefix"])
with tab2: ui_pair_inspector(ctl["prefix"])
with tab3: ui_hypotheses(ctl["prefix"])
with tab4: ui_hypothesis_scoring(ctl["prefix"])
with tab5: ui_admet(ctl["prefix"])
with tab6: ui_threshold_history()
with tab7: ui_sweep(ctl["prefix"])

# Print-only view (선택한 섹션만 아래에 재렌더)
if print_mode and print_sel:
    st.divider()
    st.subheader("🖨 인쇄용 뷰")
    if "데이터/요약" in print_sel: ui_data_summary(ctl["prefix"])
    if "Pair inspector" in print_sel: ui_pair_inspector(ctl["prefix"])
    if "가설" in print_sel: ui_hypotheses(ctl["prefix"])

# Run pipeline
if ctl["run"]:
    log_area = st.empty()
    try:
        if ctl["do_p1"]:
            ok = run_step([_py(), "src/phase1_eda.py", "--prefix", ctl["prefix"]], log_area, "Phase 1");  st.stop() if not ok else None
        if ctl["do_p2"]:
            ok = run_step([_py(), "src/phase2_cliffs.py", "--prefix", ctl["prefix"], "--config", str(CONF_FILE)], log_area, "Phase 2"); st.stop() if not ok else None
        if ctl["do_p3"]:
            if not os.getenv("OPENAI_API_KEY"):
                st.warning("OPENAI_API_KEY 환경변수가 설정되어 있지 않아 Phase3를 건너뜁니다.")
            else:
                ok = run_step([_py(), "src/phase3_llm.py", "--prefix", ctl["prefix"], "--top", str(int(ctl["top_n"]))], log_area, "Phase 3"); st.stop() if not ok else None
        if ctl["do_p4"]:
            ok = run_step([_py(), "src/build_report.py", "--prefix", ctl["prefix"]], log_area, "Phase 4"); st.stop() if not ok else None
        st.toast("완료!", icon="✅")
    except Exception as e:
        st.exception(e); st.stop()

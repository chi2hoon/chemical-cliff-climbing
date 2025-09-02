import os
import pandas as pd

from pipeline.validate.hooks import stable_sort, write_manifest, write_quarantine
from pipeline.validate.gold_checks import validate_gold_assays, validate_gold_compounds
from pipeline.transforms.chem_ids import derive_chem_ids


def _ensure_dir(p):
    d = os.path.dirname(os.path.abspath(p))
    if d and not os.path.exists(d):
        os.makedirs(d, exist_ok=True)


def build_gold(years):
    """Args: years(list[str|int]) -> dict

    Silver 산출을 모아 고정 골드 스키마 CSV를 생성한다.
    - gold/assay_readings.csv
    - gold/compounds.csv
    메타 테이블은 후속 단계에서 추가 가능.
    """
    years = [str(y) for y in years]
    # assay_readings 집계
    assay_rows = []
    for y in years:
        path = os.path.join("data", "refined", y, "assay_readings_silver.csv")
        if not os.path.exists(path):
            # 해당 연도에 측정치가 없으면 건너뜀
            continue
        df = pd.read_csv(path, dtype=str)
        # dtype 정리
        def _to_float(x):
            try:
                return float(x)
            except Exception:
                return None
        df["value_std"] = df["value_std"].map(_to_float)
        df["year"] = y
        df["qc_pass"] = df.get("qc_pass", True)
        keep = ["compound_id","target_id","assay_id","qualifier","value_std","unit_std","year","qc_pass","provenance_file","provenance_sheet","provenance_row"]
        assay_rows.append(df[keep])
    assays = pd.concat(assay_rows, ignore_index=True) if assay_rows else pd.DataFrame(columns=["compound_id","target_id","assay_id","qualifier","value_std","unit_std","year","qc_pass","provenance_file","provenance_sheet","provenance_row"])
    assays = stable_sort(assays, ["compound_id","assay_id","provenance_row"])

    # compounds 집계: 구조가 있으면 캐노니컬/인치키 생성
    comp_rows = []
    for y in years:
        path = os.path.join("data", "refined", y, "compounds_silver.csv")
        if not os.path.exists(path):
            continue
        df = pd.read_csv(path, dtype=str)
        rows = []
        for _, r in df.iterrows():
            ids = derive_chem_ids(r.get("smiles_raw"))
            rows.append({
                "compound_key": ids.get("compound_key", ""),
                "smiles_canonical": ids.get("smiles_canonical", ""),
                "has_structure": bool(ids.get("has_structure")),
                "iupac_name": r.get("iupac_name"),
                "inchikey14": ids.get("inchikey14", ""),
            })
        comp_rows.append(pd.DataFrame(rows))
    comps = pd.concat(comp_rows, ignore_index=True) if comp_rows else pd.DataFrame(columns=["compound_key","smiles_canonical","has_structure","iupac_name","inchikey14"])
    comps = stable_sort(comps, ["compound_key","smiles_canonical"])

    out_dir = os.path.join("data", "gold")
    os.makedirs(out_dir, exist_ok=True)
    assays_path = os.path.join(out_dir, "assay_readings.csv")
    comps_path = os.path.join(out_dir, "compounds.csv")
    # 검증 및 격리
    assays_ok, assays_bad = validate_gold_assays(assays)
    comps_ok, comps_bad = validate_gold_compounds(comps)
    if len(assays_bad) > 0:
        write_quarantine("gold", "all", "assay_rule_violation", assays_bad)
    if len(comps_bad) > 0:
        write_quarantine("gold", "all", "compounds_rule_violation", comps_bad)
    assays_ok.to_csv(assays_path, index=False, encoding="utf-8")
    comps_ok.to_csv(comps_path, index=False, encoding="utf-8")
    return {"assay_readings": assays_path, "compounds": comps_path}

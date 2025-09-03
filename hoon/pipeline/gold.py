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

    Silver 산출을 연도별로 Gold 고정 스키마로 저장한다.
    - data/gold/{year}/assay_readings.csv
    - data/gold/{year}/compounds.csv
    """
    years = [str(y) for y in years]
    results = {}
    for y in years:
        # assays per year
        assay_path_in = os.path.join("data", "refined", y, "assay_readings_silver.csv")
        assays = pd.DataFrame(columns=["compound_id","target_id","assay_id","qualifier","value_std","unit_std","year","qc_pass","provenance_file","provenance_sheet","provenance_row"])
        if os.path.exists(assay_path_in):
            df = pd.read_csv(assay_path_in, dtype=str)
            def _to_float(x):
                try:
                    return float(x)
                except Exception:
                    return None
            df["value_std"] = df["value_std"].map(_to_float)
            df["year"] = y
            df["qc_pass"] = df.get("qc_pass", True)
            keep = ["compound_id","target_id","assay_id","qualifier","value_std","unit_std","year","qc_pass","provenance_file","provenance_sheet","provenance_row"]
            assays = df[keep].copy()
            assays = stable_sort(assays, ["compound_id","assay_id","provenance_row"])

        # compounds per year
        comp_path_in = os.path.join("data", "refined", y, "compounds_silver.csv")
        comps = pd.DataFrame(columns=["compound_key","smiles_canonical","has_structure","iupac_name","inchikey14"])
        if os.path.exists(comp_path_in):
            dfc = pd.read_csv(comp_path_in, dtype=str)
            rows = []
            for _, r in dfc.iterrows():
                ids = derive_chem_ids(r.get("smiles_raw"))
                rows.append({
                    "compound_key": ids.get("compound_key", ""),
                    "smiles_canonical": ids.get("smiles_canonical", ""),
                    "has_structure": bool(ids.get("has_structure")),
                    "iupac_name": r.get("iupac_name"),
                    "inchikey14": ids.get("inchikey14", ""),
                })
            comps = pd.DataFrame(rows, columns=["compound_key","smiles_canonical","has_structure","iupac_name","inchikey14"]).fillna("")
            comps = stable_sort(comps, ["compound_key","smiles_canonical"])

        # 검증 및 격리(연도별)
        assays_ok, assays_bad = validate_gold_assays(assays)
        comps_ok, comps_bad = validate_gold_compounds(comps)
        if len(assays_bad) > 0:
            write_quarantine("gold", y, "assay_rule_violation", assays_bad)
        if len(comps_bad) > 0:
            write_quarantine("gold", y, "compounds_rule_violation", comps_bad)

        out_dir = os.path.join("data", "gold", y)
        os.makedirs(out_dir, exist_ok=True)
        assays_path = os.path.join(out_dir, "assay_readings.csv")
        comps_path = os.path.join(out_dir, "compounds.csv")
        assays_ok.to_csv(assays_path, index=False, encoding="utf-8")
        comps_ok.to_csv(comps_path, index=False, encoding="utf-8")
        results[y] = {"assay_readings": assays_path, "compounds": comps_path}
    return results

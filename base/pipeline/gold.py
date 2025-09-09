import os
import pandas as pd

from pipeline.validate.hooks import stable_sort, write_manifest, write_quarantine
from pipeline.validate.gold_checks import validate_gold_assays, validate_gold_compounds
from pipeline.transforms.chem_ids import derive_chem_ids
from pipeline.transforms.normalize import convert_unit, ascii_units


def _ensure_dir(p):
    d = os.path.dirname(os.path.abspath(p))
    if d and not os.path.exists(d):
        os.makedirs(d, exist_ok=True)


def build_gold(years):
    """Args: years(list[str|int]) -> dict

    Silver 산출을 연도별로 Gold 고정 스키마로 저장한다.
    - base/data/gold/{year}/assay_readings.csv
    - base/data/gold/{year}/compounds.csv
    """
    years = [str(y) for y in years]
    results = {}
    for y in years:
        assay_path_in = os.path.join("base", "data", "silver", y, "assay_readings_silver.csv")
        assays = pd.DataFrame(columns=["compound_id","target_id","assay_id","qualifier","value_std","unit_std","year","qc_pass","provenance_file","provenance_sheet","provenance_row"])
        if os.path.exists(assay_path_in):
            df = pd.read_csv(assay_path_in, dtype=str)
            # 연도/출처 간 단위 통일(기본 uM). 퍼센트는 assay_readings에 존재하지 않는다는 전제.
            def _conv(v, u):
                v2, u2 = convert_unit(v, u, 'uM')
                return v2, ascii_units(u2)
            if 'value_std' in df.columns and 'unit_std' in df.columns:
                vv, uu = [], []
                for _, r in df.iterrows():
                    v2, u2 = _conv(r.get('value_std'), r.get('unit_std'))
                    vv.append(v2)
                    uu.append(u2)
                df['value_std'] = vv
                df['unit_std'] = uu
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

        comp_path_in = os.path.join("base", "data", "silver", y, "compounds_silver.csv")
        comps = pd.DataFrame(columns=["compound_key","smiles_canonical","has_structure","iupac_name","inchikey14"])
        props_rows = []
        if os.path.exists(comp_path_in):
            dfc = pd.read_csv(comp_path_in, dtype=str)
            rows = []
            for _, r in dfc.iterrows():
                ids = derive_chem_ids(r.get("smiles_raw"), r.get("iupac_name"))
                rows.append({
                    "compound_key": ids.get("compound_key", ""),
                    "smiles_canonical": ids.get("smiles_canonical", ""),
                    "has_structure": bool(ids.get("has_structure")),
                    "iupac_name": r.get("iupac_name"),
                    "inchikey14": ids.get("inchikey14", ""),
                })
                props_rows.append({
                    "compound_key": ids.get("compound_key", ""),
                    "compound_id": r.get("compound_id"),
                    "iupac_name": r.get("iupac_name"),
                    "mw": r.get("mw"),
                    "lcms_text": r.get("lcms_text"),
                    "nmr_1h_text": r.get("nmr_1h_text"),
                    "flag_asterisk": r.get("flag_asterisk"),
                    "flag_imaging_conflict": r.get("flag_imaging_conflict"),
                    "flag_smiles_o3_changed": r.get("flag_smiles_o3_changed"),
                    "year": y,
                    "provenance_file": r.get("provenance_file"),
                    "provenance_sheet": r.get("provenance_sheet"),
                    "provenance_row": r.get("provenance_row"),
                })
            comps = pd.DataFrame(rows, columns=["compound_key","smiles_canonical","has_structure","iupac_name","inchikey14"]).fillna("")
            comps = stable_sort(comps, ["compound_key","smiles_canonical"])

        assays_ok, assays_bad = validate_gold_assays(assays)
        comps_ok, comps_bad = validate_gold_compounds(comps)
        if len(assays_bad) > 0:
            write_quarantine("gold", y, "assay_rule_violation", assays_bad)
        if len(comps_bad) > 0:
            write_quarantine("gold", y, "compounds_rule_violation", comps_bad)

        out_dir = os.path.join("base", "data", "gold", y)
        os.makedirs(out_dir, exist_ok=True)
        assays_path = os.path.join(out_dir, "assay_readings.csv")
        comps_path = os.path.join(out_dir, "compounds.csv")
        assays_ok.to_csv(assays_path, index=False, encoding="utf-8")
        comps_ok.to_csv(comps_path, index=False, encoding="utf-8")
        if props_rows:
            props_df = pd.DataFrame(props_rows, columns=["compound_key","compound_id","iupac_name","mw","lcms_text","nmr_1h_text","flag_asterisk","flag_imaging_conflict","flag_smiles_o3_changed","year","provenance_file","provenance_sheet","provenance_row"]).fillna("")
            props_df = stable_sort(props_df, ["compound_key","compound_id","provenance_row"]) if "compound_key" in props_df.columns else props_df
            props_out = os.path.join(out_dir, "compound_props.csv")
            props_df.to_csv(props_out, index=False, encoding="utf-8")
        meta_in = os.path.join("base","data","silver", y, "assay_context_silver.csv")
        if os.path.exists(meta_in):
            meta_df = pd.read_csv(meta_in, dtype=str)
            meta_df = stable_sort(meta_df, [c for c in ["compound_id","provenance_row"] if c in meta_df.columns])
            meta_out = os.path.join(out_dir, "assay_context.csv")
            meta_df.to_csv(meta_out, index=False, encoding="utf-8")
        results[y] = {"assay_readings": assays_path, "compounds": comps_path}
    return results

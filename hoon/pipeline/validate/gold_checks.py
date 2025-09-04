import pandas as pd


def validate_gold_assays(df):
    """Args: df(DataFrame) -> (ok_df, bad_df)

    골드 assay_readings 최소 규칙 검사.
    """
    req = ["compound_id","target_id","assay_id","qualifier","value_std","unit_std","year","provenance_file","provenance_sheet","provenance_row"]
    missing = [c for c in req if c not in df.columns]
    if missing:
        bad = df.copy()
        bad["reason"] = "missing_columns:" + ",".join(missing)
        return df.iloc[0:0].copy(), bad
    allowed_units = {"uM","nM","%"}
    unit_ok = df["unit_std"].isin(list(allowed_units)) | df["unit_std"].isna()
    qual_ok = df["qualifier"].isin(["=", ">", "<"]) | df["qualifier"].isna()
    def _to_float(x):
        try:
            return float(x)
        except Exception:
            return None
    val = df["value_std"].map(_to_float)
    val_ok = val.isna() | (val > 0)
    ok = df[unit_ok & qual_ok & val_ok].copy()
    bad = df[~(unit_ok & qual_ok & val_ok)].copy()
    if len(bad) > 0:
        bad["reason"] = "rule_violation"
    return ok, bad


def validate_gold_compounds(df):
    """Args: df(DataFrame) -> (ok_df, bad_df)

    골드 compounds 최소 규칙 검사.
    has_structure True면 inchikey가 있어야 한다는 규칙을 확인한다.
    """
    df2 = df.copy()
    has_struct = df2.get("has_structure")
    inchikey = df2.get("compound_key")
    if has_struct is None or inchikey is None:
        return df2, df2.iloc[0:0].copy()
    mask_bad = (has_struct == True) & ((inchikey.astype(str).str.strip() == "") | inchikey.isna())
    bad = df2[mask_bad].copy()
    if len(bad) > 0:
        bad["reason"] = "has_structure_without_inchikey"
    ok = df2[~mask_bad].copy()
    return ok, bad


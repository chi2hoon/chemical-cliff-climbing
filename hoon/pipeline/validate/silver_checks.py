import pandas as pd


def _required_cols(df, cols):
    missing = [c for c in cols if c not in df.columns]
    if missing:
        return False, missing
    return True, []


def validate_silver(compounds_df, assays_df):
    """Args: compounds_df(DataFrame), assays_df(DataFrame) -> (ok:dict, quarantine:dict)

    Silver 레이어 최소 검증 및 격리 대상 분리.
    """
    ok = {"compounds": compounds_df.copy(), "assays": assays_df.copy()}
    quarantine = {}

    # 필수 컬럼 검사
    req_comp = ["compound_id", "provenance_file", "provenance_sheet", "provenance_row"]
    present, missing = _required_cols(ok["compounds"], req_comp)
    if not present:
        quarantine["compounds_missing_columns"] = ok["compounds"].copy()
        ok["compounds"] = ok["compounds"].iloc[0:0].copy()
    else:
        mask = None
        for c in req_comp:
            m = ok["compounds"][c].isna() | (ok["compounds"][c].astype(str).str.strip() == "")
            mask = m if mask is None else (mask | m)
        bad = ok["compounds"][mask].copy()
        if len(bad) > 0:
            bad["reason"] = "missing_required_value"
            quarantine["compounds_missing_required_value"] = bad
            ok["compounds"] = ok["compounds"][~mask].copy()

    req_assay = ["compound_id","target_id","assay_id","qualifier","value_std","unit_std","year","provenance_file","provenance_sheet","provenance_row"]
    present, missing = _required_cols(ok["assays"], req_assay)
    if not present:
        quarantine["assays_missing_columns"] = ok["assays"].copy()
        ok["assays"] = ok["assays"].iloc[0:0].copy()
    else:
        # 규칙
        df = ok["assays"].copy()
        # 단위
        allowed_units = {"uM","nM","%"}
        df["unit_ok"] = df["unit_std"].isin(list(allowed_units)) | df["unit_std"].isna()
        # qualifier
        df["qual_ok"] = df["qualifier"].isin(["=", ">", "<"]) | df["qualifier"].isna()
        # value_std > 0 when present
        def _to_float(x):
            try:
                return float(x)
            except Exception:
                return None
        df["value_num"] = df["value_std"].map(_to_float)
        df["val_ok"] = df["value_num"].isna() | (df["value_num"] > 0)
        bad = df.loc[~(df["unit_ok"] & df["qual_ok"] & df["val_ok"])].copy()
        if len(bad) > 0:
            bad["reason"] = "rule_violation"
            quarantine["assays_rule_violation"] = bad[req_assay + ["reason"]]
            df = df.loc[df["unit_ok"] & df["qual_ok"] & df["val_ok"]]
        # 색상 플래그가 있으면 격리
        flag_cols = [c for c in ["flag_asterisk","flag_imaging_conflict"] if c in df.columns]
        if flag_cols:
            flagged_mask = False
            for fc in flag_cols:
                flagged_mask = (df[fc] == True) if flagged_mask is False else (flagged_mask | (df[fc] == True))
            flagged = df[flagged_mask].copy()
            if len(flagged) > 0:
                flagged["reason"] = "flagged_by_color"
                quarantine["assays_flagged"] = flagged[req_assay + ["reason"]]
                df = df[~flagged_mask]
        ok["assays"] = df[req_assay].copy()

    return ok, quarantine

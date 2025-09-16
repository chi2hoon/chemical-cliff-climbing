import pandas as pd


def validate_silver(df):
    """Args: df(DataFrame) -> (ok_df, bad_df)

    Silver 산출 최소 규칙 검증.
    """
    req = ["compound_id","provenance_file","provenance_sheet","provenance_row"]
    missing = [c for c in req if c not in df.columns]
    if missing:
        return df.iloc[0:0].copy(), df.copy()
    ok = df.copy()
    bad = pd.DataFrame(columns=df.columns)
    return ok, bad


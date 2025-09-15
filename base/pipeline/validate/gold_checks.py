import pandas as pd


def validate_gold_assays(df):
    """Args: df(DataFrame) -> (ok_df, bad_df)

    Gold assay_readings 최소 규칙 검증.
    """
    req = ["compound_id","target_id","assay_id","qualifier","value_std","unit_std","year"]
    missing = [c for c in req if c not in df.columns]
    if missing:
        return df.iloc[0:0].copy(), df.copy()
    ok = df.copy()
    bad = pd.DataFrame(columns=df.columns)
    return ok, bad


def validate_gold_compounds(df):
    """Args: df(DataFrame) -> (ok_df, bad_df)

    Gold compounds 최소 규칙 검증.
    """
    req = ["compound_key","smiles_canonical","has_structure"]
    missing = [c for c in req if c not in df.columns]
    if missing:
        return df.iloc[0:0].copy(), df.copy()
    ok = df.copy()
    bad = pd.DataFrame(columns=df.columns)
    return ok, bad


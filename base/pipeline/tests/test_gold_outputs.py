import os
import pandas as pd


def test_gold_outputs_columns_2017():
    ar = os.path.join("base", "data", "gold", "2017", "assay_readings.csv")
    cp = os.path.join("base", "data", "gold", "2017", "compounds.csv")
    assert os.path.exists(ar)
    assert os.path.exists(cp)
    df_ar = pd.read_csv(ar)
    req_ar = [
        "compound_id","target_id","assay_id","qualifier","value_std","unit_std","year",
        "provenance_file","provenance_sheet","provenance_row"
    ]
    for c in req_ar:
        assert c in df_ar.columns
    df_cp = pd.read_csv(cp)
    req_cp = ["compound_key","smiles_canonical","has_structure"]
    for c in req_cp:
        assert c in df_cp.columns


import os
import pandas as pd


def test_silver_outputs_columns_2017():
    path = os.path.join("base", "data", "silver", "2017", "assay_readings_silver.csv")
    assert os.path.exists(path)
    df = pd.read_csv(path)
    required = [
        "compound_id","target_id","assay_id","qualifier","value_std","unit_std","year",
        "provenance_file","provenance_sheet","provenance_row"
    ]
    for c in required:
        assert c in df.columns


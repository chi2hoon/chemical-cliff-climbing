import os
import pandas as pd


def test_silver_assay_columns_2018():
    path = os.path.join("data", "silver", "2018", "assay_readings_silver.csv")
    assert os.path.exists(path)
    df = pd.read_csv(path)
    required = [
        "compound_id","target_id","assay_id","qualifier","value_std","unit_std","year","provenance_file","provenance_sheet","provenance_row"
    ]
    for c in required:
        assert c in df.columns


def test_manifest_2018_exists():
    path = os.path.join("logs", "manifest", "2018.json")
    assert os.path.exists(path)
    import json
    m = json.load(open(path, encoding="utf-8"))
    # 입력/스키마 해시 중 하나 이상 존재해야 함
    assert ("input_sha256" in m) or ("yaml_sha256" in m)
    assert "rows_out" in m


def test_silver_compounds_columns_2020_2021():
    for y in ("2020", "2021"):
        path = os.path.join("data", "silver", y, "compounds_silver.csv")
        assert os.path.exists(path)
        df = pd.read_csv(path)
        for c in ["compound_id", "smiles_raw", "provenance_file", "provenance_sheet", "provenance_row"]:
            assert c in df.columns


def test_silver_2017_bridge_outputs_exist():
    for fname in ("compounds_silver.csv", "assay_readings_silver.csv"):
        path = os.path.join("data", "silver", "2017", fname)
        assert os.path.exists(path)

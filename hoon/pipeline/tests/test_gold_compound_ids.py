import os
import pandas as pd


def test_gold_compound_ids_when_structure():
    path = os.path.join("data", "gold", "compounds.csv")
    assert os.path.exists(path)
    df = pd.read_csv(path)
    if "has_structure" in df.columns and df["has_structure"].any():
        sub = df[df["has_structure"] == True]
        assert (sub["compound_key"].astype(str).str.strip() != "").all()


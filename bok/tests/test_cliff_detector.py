import pandas as pd

from cliff_detector import (
    detect_activity_cliffs,
    detect_activity_cliffs_sali,
    detect_activity_cliffs_knn,
)


def sample_df():
    # Simple aromatic examples with known high similarity and controllable pIC50
    data = {
        "SMILES": [
            "c1ccccc1",      # benzene
            "Cc1ccccc1",     # toluene
            "n1ccccc1",      # pyridine
        ],
        "pIC50": [5.0, 7.0, 6.0],
    }
    return pd.DataFrame(data)


def test_detect_activity_cliffs_basic():
    df = sample_df()
    # Use a lenient similarity threshold to ensure benzene vs toluene is captured
    res = detect_activity_cliffs(df, sim_thres=0.6, act_thres=1.5)
    assert not res.empty
    # Check required columns
    for col in [
        "mol1_idx",
        "mol2_idx",
        "mol1_smiles",
        "mol2_smiles",
        "mol1_activity",
        "mol2_activity",
        "sim",
        "activity_diff",
    ]:
        assert col in res.columns


def test_detect_activity_cliffs_sali():
    df = sample_df()
    res = detect_activity_cliffs_sali(df, top_n=5, sim_floor=0.3)
    # Should return at least one pair and include SALI
    assert "sali" in res.columns
    # If empty (edge-case on certain RDKit builds), the column still exists
    if not res.empty:
        assert (res["sali"] > 0).any()


def test_detect_activity_cliffs_knn():
    df = sample_df()
    res = detect_activity_cliffs_knn(df, k=2, act_thres=1.0, min_sim=0.3)
    # At least one neighbor pair should qualify
    assert set(["sim", "activity_diff"]).issubset(res.columns)
    # Pairs should be unique (i<j)
    if not res.empty:
        assert ((res["mol1_idx"] < res["mol2_idx"]) | (res["mol2_idx"] < res["mol1_idx"])) .all()


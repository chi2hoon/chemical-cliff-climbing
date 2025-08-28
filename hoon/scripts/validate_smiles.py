#!/usr/bin/env python3
from __future__ import annotations

import argparse
import os
from typing import Dict

import pandas as pd
import yaml

from src.udm.io_csv import read_csv_safe, write_csv_safe
from src.udm.silver_smiles import _try_import_rdkit


def load_cfg(cfg_path: str) -> Dict:
    with open(cfg_path, "r", encoding="utf-8") as f:
        return yaml.safe_load(f)


def compute_canonical(smiles_raw: str) -> Dict[str, str]:
    Chem, inchi = _try_import_rdkit()
    can = ""
    ik = ""
    if smiles_raw:
        mol = Chem.MolFromSmiles(str(smiles_raw))
        if mol is not None:
            can = Chem.MolToSmiles(mol, canonical=True)
            try:
                ik = inchi.MolToInchiKey(mol)
            except Exception:
                ik = ""
    return {"smiles_canonical_rdkit": can, "inchikey_rdkit": ik}


def main() -> int:
    parser = argparse.ArgumentParser(description="Validate canonical SMILES by recomputing with RDKit and comparing to silver output")
    parser.add_argument("--config", default="configs/2017.yml", help="Path to YAML config")
    parser.add_argument("--root", default=os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "data"), help="Data root directory (defaults to ./data)")
    parser.add_argument("--outdir", default=os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "data", "silver", "all"), help="Directory to write validation CSVs")
    args = parser.parse_args()

    cfg = load_cfg(args.config)
    root_dir = os.path.abspath(args.root)
    dataset_id = str(cfg.get("dataset_id", ""))

    bronze_dir = os.path.join(root_dir, cfg["paths"]["bronze_dir"])
    comp_path = os.path.join(bronze_dir, "compounds.csv")
    silver_path = os.path.join(root_dir, "silver", dataset_id, "compounds_canonical.csv")

    if not os.path.exists(comp_path):
        raise FileNotFoundError(comp_path)
    if not os.path.exists(silver_path):
        raise FileNotFoundError(silver_path)

    df_bronze = read_csv_safe(comp_path)
    df_silver = read_csv_safe(silver_path)

    # bronze smiles_raw에서 표준 형태 재계산
    calc_rows = []
    for _, r in df_bronze.iterrows():
        compound_id = str(r.get("compound_id", ""))
        smiles_raw = str(r.get("smiles_raw", ""))
        rec = {"compound_id": compound_id, "smiles_raw": smiles_raw}
        rec.update(compute_canonical(smiles_raw))
        calc_rows.append(rec)
    df_calc = pd.DataFrame(calc_rows)

    # silver 표준 형태와 결합
    cols_keep = [
        "compound_id",
        "dataset_id",
        "smiles_raw",
        "smiles_canonical",
        "inchikey",
    ]
    df_silver_small = df_silver[[c for c in cols_keep if c in df_silver.columns]].copy()
    merged = df_calc.merge(df_silver_small, on=["compound_id", "smiles_raw"], how="left", suffixes=("_calc", "_silver"))

    # 동등성 비교
    merged["match_smiles"] = (merged.get("smiles_canonical_rdkit", "").astype(str) == merged.get("smiles_canonical", "").astype(str))
    merged["match_inchikey"] = (merged.get("inchikey_rdkit", "").astype(str) == merged.get("inchikey", "").astype(str))

    os.makedirs(args.outdir, exist_ok=True)
    out_all = os.path.join(args.outdir, f"{dataset_id}_smiles_validation.csv")
    out_bad = os.path.join(args.outdir, f"{dataset_id}_smiles_mismatches.csv")
    write_csv_safe(merged, out_all)
    write_csv_safe(merged.loc[~(merged["match_smiles"] & merged["match_inchikey"])], out_bad)

    print(f"Wrote: {out_all}")
    print(f"Wrote mismatches: {out_bad}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())





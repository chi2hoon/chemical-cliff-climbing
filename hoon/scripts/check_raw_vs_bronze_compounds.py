#!/usr/bin/env python3
from __future__ import annotations

import argparse
import os
from typing import Dict, List

import pandas as pd
import yaml

from src.udm.io_csv import read_csv_safe, write_csv_safe


def load_cfg(cfg_path: str) -> Dict:
    with open(cfg_path, "r", encoding="utf-8") as f:
        return yaml.safe_load(f)


def read_raw_excel(raw_path: str, sheet: str, id_col: str, smiles_col: str, header: int | None) -> pd.DataFrame:
    # Bronze 단계 규칙: 타입을 추측하지 말고 문자열로 읽어서 공백 유지
    df = pd.read_excel(
        raw_path,
        sheet_name=sheet,
        header=header,
        dtype=str,
        keep_default_na=False,
    )
    # 열 이름을 문자열로 정규화 (정확한 매치 예상)
    cols = [str(c) for c in df.columns]
    if id_col not in cols or smiles_col not in cols:
        raise ValueError(f"Column not found. Available: {cols}. Need id_col='{id_col}', smiles_col='{smiles_col}'")
    out = df[[id_col, smiles_col]].copy()
    out.columns = ["compound_id", "smiles_raw_rawexcel"]
    # 문자열 보장
    out["compound_id"] = out["compound_id"].astype(str)
    out["smiles_raw_rawexcel"] = out["smiles_raw_rawexcel"].astype(str)
    return out


def main() -> int:
    parser = argparse.ArgumentParser(description="Check raw Excel vs bronze compounds consistency by compound_id and smiles")
    parser.add_argument("--config", default="configs/2017.yml", help="Path to YAML config")
    parser.add_argument("--root", default=os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "data"), help="Data root directory (defaults to ./data)")
    parser.add_argument("--sheet", required=True, help="Raw Excel sheet name containing compound_id and SMILES")
    parser.add_argument("--id-col", required=True, help="Column name in raw sheet for compound_id")
    parser.add_argument("--smiles-col", required=True, help="Column name in raw sheet for SMILES")
    parser.add_argument("--header", type=int, default=0, help="Row index to use as header in raw sheet (use -1 for no header)")
    parser.add_argument("--outdir", default=os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "data", "bronze", "2017"), help="Directory to write comparison outputs")
    args = parser.parse_args()

    cfg = load_cfg(args.config)
    root_dir = os.path.abspath(args.root)
    raw_path = os.path.join(root_dir, cfg["paths"]["raw_file"])  # absolute raw file path
    bronze_comp_path = os.path.join(root_dir, cfg["paths"]["bronze_files"]["compounds"])  # bronze compounds.csv

    if not os.path.exists(raw_path):
        raise FileNotFoundError(raw_path)
    if not os.path.exists(bronze_comp_path):
        raise FileNotFoundError(bronze_comp_path)

    header = None if args.header < 0 else args.header
    df_raw = read_raw_excel(raw_path, args.sheet, args.id_col, args.smiles_col, header=header)
    df_bronze = read_csv_safe(bronze_comp_path)

    # 문자열 타입 강제 적용
    df_bronze["compound_id"] = df_bronze["compound_id"].astype(str)
    if "smiles_raw" in df_bronze.columns:
        df_bronze["smiles_raw"] = df_bronze["smiles_raw"].astype(str)
    else:
        df_bronze["smiles_raw"] = ""

    # compound_id로 raw->bronze 왼쪽 결합
    merged = df_raw.merge(df_bronze[["compound_id", "smiles_raw"]], on="compound_id", how="left")
    merged["id_found_in_bronze"] = merged["smiles_raw"].notna()
    merged["smiles_equal"] = merged.apply(lambda r: str(r.get("smiles_raw_rawexcel", "")) == str(r.get("smiles_raw", "")), axis=1)

    # 출력
    os.makedirs(args.outdir, exist_ok=True)
    out_all = os.path.join(args.outdir, "raw_vs_bronze_compounds_comparison.csv")
    out_missing = os.path.join(args.outdir, "raw_ids_missing_in_bronze.csv")
    out_mismatch = os.path.join(args.outdir, "raw_bronze_smiles_mismatch.csv")

    write_csv_safe(merged, out_all)
    write_csv_safe(merged.loc[~merged["id_found_in_bronze"]], out_missing)
    write_csv_safe(merged.loc[merged["id_found_in_bronze"] & ~merged["smiles_equal"]], out_mismatch)

    # 간단한 요약 출력
    total = len(merged)
    num_missing = int((~merged["id_found_in_bronze"]).sum())
    num_mismatch = int((merged["id_found_in_bronze"] & ~merged["smiles_equal"]).sum())
    num_match = total - num_missing - num_mismatch
    print(f"total={total} match={num_match} missing_ids={num_missing} smiles_mismatch={num_mismatch}")
    print(f"Wrote: {out_all}\nMissing IDs: {out_missing}\nSMILES mismatches: {out_mismatch}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())





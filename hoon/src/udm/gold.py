from __future__ import annotations

import os
from typing import Dict, List, Optional

import pandas as pd

from .io_csv import read_csv_safe, write_csv_safe


def _to_float(x: object) -> Optional[float]:
	try:
		return float(x)
	except Exception:
		return None


def _load_measurements_std(root_dir: str, cfg: Dict) -> pd.DataFrame:
	silver_dir = os.path.join(root_dir, "silver", cfg.get("dataset_id", ""))
	path = os.path.join(silver_dir, "measurements_std.csv")
	if not os.path.exists(path):
		raise FileNotFoundError(path)
	df = read_csv_safe(path)
	df["dataset_id"] = cfg.get("dataset_id", "")
	return df


def _load_compounds_canonical(root_dir: str, cfg: Dict) -> pd.DataFrame:
	silver_dir = os.path.join(root_dir, "silver", cfg.get("dataset_id", ""))
	path = os.path.join(silver_dir, "compounds_canonical.csv")
	if not os.path.exists(path):
		raise FileNotFoundError(path)
	df = read_csv_safe(path)
	# expected columns: compound_id, dataset_id, smiles_raw, smiles_canonical, inchikey
	keep = [c for c in ["compound_id", "dataset_id", "smiles_canonical", "inchikey"] if c in df.columns]
	return df[keep]


def build_measurements_gold(root_dir: str, cfg: Dict) -> str:
	"""
	Create gold dataset by joining standardized measurements with canonical SMILES and InChIKey.

	Inputs:
	- silver/{dataset_id}/measurements_std.csv
	- silver/{dataset_id}/compounds_canonical.csv

	Output:
	- gold/{dataset_id}/measurements_gold.csv

	Columns (when available):
	- compound_id, dataset_id, assay_id, panel_id, cell_line
	- smiles_canonical, inchikey
	- value_std (float), unit_std, censor
	- readout, matrix
	- provenance_file, provenance_sheet, provenance_row
	- std_rule_id, std_confidence

	Rows are filtered to ensure non-empty smiles_canonical and numeric value_std.
	"""
	meas = _load_measurements_std(root_dir, cfg)
	comp = _load_compounds_canonical(root_dir, cfg)

	# join on compound_id within dataset
	if "dataset_id" in comp.columns and "dataset_id" in meas.columns:
		merged = meas.merge(comp, on=["compound_id"], how="left", suffixes=("", "_comp"))
	else:
		merged = meas.merge(comp, on=["compound_id"], how="left", suffixes=("", "_comp"))

	# numeric cast and filtering
	merged["value_std_num"] = merged.get("value_std").apply(_to_float)
	merged = merged[(merged["value_std_num"].notna()) & (merged.get("smiles_canonical", "").astype(str).str.len() > 0)]
	# ensure value_std is float for downstream usage before selecting columns
	merged["value_std"] = merged["value_std_num"].astype(float)

	# select and order columns
	front_cols: List[str] = [
		"compound_id", "dataset_id", "assay_id", "panel_id", "cell_line",
		"smiles_canonical", "inchikey",
		"value_std", "unit_std", "censor",
		"readout", "matrix",
		"provenance_file", "provenance_sheet", "provenance_row",
		"std_rule_id", "std_confidence",
	]
	present_front = [c for c in front_cols if c in merged.columns]
	other_cols = [c for c in merged.columns if c not in present_front and c != "value_std_num"]
	out_df = merged[present_front + other_cols].copy()
	# drop helper if it slipped through
	if "value_std_num" in out_df.columns:
		out_df = out_df.drop(columns=["value_std_num"])

	gold_dir = os.path.join(root_dir, "gold", cfg.get("dataset_id", ""))
	os.makedirs(gold_dir, exist_ok=True)
	out_path = os.path.join(gold_dir, "measurements_gold.csv")
	write_csv_safe(out_df, out_path)
	return out_path



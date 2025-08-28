from __future__ import annotations

import os
import re
from typing import Dict, List, Tuple

import pandas as pd

from .io_csv import read_csv_safe, write_csv_safe



def _save_csv(df: pd.DataFrame, root_dir: str, cfg: Dict, filename: str) -> str:
	out_dir = os.path.join(root_dir, cfg["paths"]["bronze_dir"])  # unify under bronze/{year}
	os.makedirs(out_dir, exist_ok=True)
	path = os.path.join(out_dir, filename)
	write_csv_safe(df, path)
	return path


def qc_reports(root_dir: str, cfg: Dict) -> Dict[str, str]:
	paths = cfg["paths"]["bronze_files"]
	measurements_path = os.path.join(root_dir, paths["measurements"])  # bronze source
	assays_path = os.path.join(root_dir, paths["assays"])  # for assay_id
	panel_meta_path = os.path.join(root_dir, cfg["paths"]["bronze_files"]["panel_meta"]) if "panel_meta" in cfg["paths"]["bronze_files"] else None

	df_m = read_csv_safe(measurements_path)
	df_a = read_csv_safe(assays_path)

	outputs: Dict[str, str] = {}

	# qc_bad_compound_id.csv
	bad_re = re.compile(cfg["qc"]["bad_id_regex"]) if cfg["qc"].get("bad_id_regex") else None
	if bad_re is not None and "compound_id" in df_m.columns:
		mask_bad = df_m["compound_id"].astype(str).str.match(bad_re) == False
		outputs["qc_bad_compound_id.csv"] = _save_csv(df_m.loc[mask_bad], root_dir, cfg, "qc_bad_compound_id.csv")

	# qc_dup_pair.csv (일부 열이 누락된 경우에도 정상 작동)
	pair_cols_cfg = cfg["qc"]["check_duplicates"]["pair"]
	pair_cols = [c for c in pair_cols_cfg if c in df_m.columns]
	if len(pair_cols) >= 1:
		dup_pair = df_m[df_m.duplicated(subset=pair_cols, keep=False)].sort_values(pair_cols)
	else:
		dup_pair = df_m.iloc[0:0]
	outputs["qc_dup_pair.csv"] = _save_csv(dup_pair, root_dir, cfg, "qc_dup_pair.csv")

	# qc_dup_triplet.csv (일부 열이 누락된 경우에도 정상 작동)
	triplet_cols_cfg = cfg["qc"]["check_duplicates"]["triplet"]
	triplet_cols = [c for c in triplet_cols_cfg if c in df_m.columns]
	if len(triplet_cols) >= 1:
		dup_triplet = df_m[df_m.duplicated(subset=triplet_cols, keep=False)].sort_values(triplet_cols)
	else:
		dup_triplet = df_m.iloc[0:0]
	outputs["qc_dup_triplet.csv"] = _save_csv(dup_triplet, root_dir, cfg, "qc_dup_triplet.csv")

	# qc_missing_assay_id.csv
	if "assay_id" in df_m.columns:
		missing_assay = df_m[df_m["assay_id"].astype(str).str.len() == 0]
	else:
		missing_assay = df_m.iloc[0:0]
	outputs["qc_missing_assay_id.csv"] = _save_csv(missing_assay, root_dir, cfg, "qc_missing_assay_id.csv")

	# qc_nonstandard_cell_label.csv: 매핑 키로 나타나는 레이블 플래그 (비표준 형태)
	cell_map = cfg.get("parsing", {}).get("cell_line_map", {})
	if cell_map and "cell_line" in df_m.columns:
		nonstd_keys = set(cell_map.keys())
		mask_nonstd = df_m["cell_line"].isin(list(nonstd_keys))
		out_df = df_m.loc[mask_nonstd, ["cell_line"]].drop_duplicates().sort_values("cell_line")
	else:
		out_df = pd.DataFrame(columns=["cell_line"])  # empty placeholder
	outputs["qc_nonstandard_cell_label.csv"] = _save_csv(out_df, root_dir, cfg, "qc_nonstandard_cell_label.csv")

	# qc_overall.md (bronze/{year}에 작성)
	rows_total = len(df_m)
	missing_pct = (df_m.isna() | (df_m.astype(str) == "")).mean().round(3)
	# 검열 분포
	censor_counts = df_m["censor"].value_counts(dropna=False).to_dict() if "censor" in df_m.columns else {}
	md_lines: List[str] = []
	md_lines.append(f"# QC Overall - {cfg.get('dataset_id','')}\n")
	md_lines.append(f"Total rows: {rows_total}\n")
	md_lines.append("\n## Missing percent by column\n")
	md_lines.append("column,missing_pct\n")
	for col, pct in missing_pct.items():
		md_lines.append(f"{col},{pct}\n")
	md_lines.append("\n## Censor distribution\n")
	md_lines.append("censor,count\n")
	for k, v in censor_counts.items():
		md_lines.append(f"{k},{v}\n")
	bronze_dir = os.path.join(root_dir, cfg["paths"]["bronze_dir"])
	os.makedirs(bronze_dir, exist_ok=True)
	md_path = os.path.join(bronze_dir, "qc_overall.md")
	with open(md_path, "w", encoding="utf-8") as f:
		f.writelines(md_lines)
	outputs["qc_overall.md"] = md_path

	return outputs



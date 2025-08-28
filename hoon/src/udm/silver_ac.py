from __future__ import annotations

import os
from typing import Dict, List, Tuple

import pandas as pd

from .io_csv import read_csv_safe, write_csv_safe


def _try_import_rdkit():
	try:
		from rdkit import Chem  # type: ignore
		from rdkit.Chem import AllChem, DataStructs  # type: ignore
		return Chem, AllChem, DataStructs
	except Exception:
		raise RuntimeError("RDKit is required for activity cliff calculation. Use an environment with rdkit.")


def _load_measurements_std(data_root: str, cfg: Dict) -> pd.DataFrame:
	silver_dir = os.path.join(data_root, "silver", cfg.get("dataset_id", ""))
	path = os.path.join(silver_dir, "measurements_std.csv")
	if not os.path.exists(path):
		return pd.DataFrame()
	df = read_csv_safe(path)
	df["dataset_id"] = cfg.get("dataset_id", "")
	return df


def _load_canonical_compounds(data_root: str, cfg: Dict) -> pd.DataFrame:
	silver_dir = os.path.join(data_root, "silver", cfg.get("dataset_id", ""))
	path = os.path.join(silver_dir, "compounds_canonical.csv")
	if not os.path.exists(path):
		return pd.DataFrame(columns=["compound_id", "smiles_canonical"])
	try:
		df = read_csv_safe(path)
		cols = list(df.columns)
		if (len(cols) == 0) or ("compound_id" not in cols) or ("smiles_canonical" not in cols):
			return pd.DataFrame(columns=["compound_id", "smiles_canonical"])
		return df[["compound_id", "smiles_canonical"]]
	except Exception:
		# 빈/유효하지 않은 표준 파일을 우아하게 처리
		return pd.DataFrame(columns=["compound_id", "smiles_canonical"])


def build_activity_cliffs(
	data_root: str,
	cfgs: List[Dict],
	similarity_threshold: float = 0.85,
	delta_threshold: float = 1.0,
) -> str:
	"""
	Compute activity cliffs across multiple datasets and write CSV to data/silver/all/ac_pairs.csv.
	Cliff definition: Tanimoto(Morgan2) >= similarity_threshold AND |value_std_i - value_std_j| >= delta_threshold
	Grouping: within same assay_id and (if present) same cell_line.
	"""
	Chem, AllChem, DataStructs = _try_import_rdkit()

	# 측정값 로드 및 결합
	meas_list: List[pd.DataFrame] = []
	comp_map_list: List[pd.DataFrame] = []
	for cfg in cfgs:
		m = _load_measurements_std(data_root, cfg)
		if not m.empty:
			meas_list.append(m)
		c = _load_canonical_compounds(data_root, cfg)
		if not c.empty:
			c["dataset_id"] = cfg.get("dataset_id", "")
			comp_map_list.append(c)
	if not meas_list:
		raise RuntimeError("No measurements_std.csv available for the provided configs.")
	meas = pd.concat(meas_list, ignore_index=True)
	comp_map = pd.concat(comp_map_list, ignore_index=True) if comp_map_list else pd.DataFrame(columns=["compound_id", "smiles_canonical", "dataset_id"])

	# 표준 SMILES 결합
	meas = meas.merge(comp_map, on=["compound_id"], how="left")

	# 지문 캐시 준비
	fp_cache: Dict[str, object] = {}
	def get_fp(smiles: str):
		if not isinstance(smiles, str) or smiles == "":
			return None
		if smiles in fp_cache:
			return fp_cache[smiles]
		mol = Chem.MolFromSmiles(smiles)
		if mol is None:
			fp_cache[smiles] = None
			return None
		fp = AllChem.GetMorganFingerprintAsBitVect(mol, radius=2, nBits=2048)
		fp_cache[smiles] = fp
		return fp

	# 숫자형 보장
	def to_float(x):
		try:
			return float(x)
		except Exception:
			return None
	meas["value_std_num"] = meas["value_std"].apply(to_float)

	# 그룹화 열 선택
	group_cols = [c for c in ["assay_id", "cell_line"] if c in meas.columns]
	if not group_cols:
		group_cols = ["assay_id"] if "assay_id" in meas.columns else []

	rows: List[Dict] = []
	for group_key, g in meas.groupby(group_cols):
		g = g.dropna(subset=["compound_id"])  # keep valid compounds
		g = g[(g["value_std_num"].notna()) & (g["smiles_canonical"].astype(str).str.len() > 0)]
		if len(g) < 2:
			continue
		# 인덱스로 쌍 구축
		if "dataset_id" not in g.columns:
			g = g.copy()
			g["dataset_id"] = ""
		records = g[["compound_id", "dataset_id", "smiles_canonical", "value_std_num"]].reset_index(drop=True)
		n = len(records)
		for i in range(n):
			ci = records.loc[i]
			fp_i = get_fp(ci["smiles_canonical"])
			if fp_i is None:
				continue
			for j in range(i + 1, n):
				cj = records.loc[j]
				fp_j = get_fp(cj["smiles_canonical"])
				if fp_j is None:
					continue
				sim = DataStructs.TanimotoSimilarity(fp_i, fp_j)
				delta = abs(ci["value_std_num"] - cj["value_std_num"])
				if sim >= similarity_threshold and delta >= delta_threshold:
					row: Dict[str, object] = {
						"assay_id": g["assay_id"].iloc[0] if "assay_id" in g.columns else "",
						"cell_line": g["cell_line"].iloc[0] if "cell_line" in g.columns else "",
						"compound_id_1": ci["compound_id"],
						"compound_id_2": cj["compound_id"],
						"dataset_id_1": ci["dataset_id"],
						"dataset_id_2": cj["dataset_id"],
						"smiles_1": ci["smiles_canonical"],
						"smiles_2": cj["smiles_canonical"],
						"similarity": sim,
						"value_std_1": ci["value_std_num"],
						"value_std_2": cj["value_std_num"],
						"delta_value_std": delta,
					}
					# 그룹 키 추가
					if isinstance(group_key, tuple):
						pass
					rows.append(row)

	all_dir = os.path.join(data_root, "silver", "all")
	os.makedirs(all_dir, exist_ok=True)
	out_path = os.path.join(all_dir, "ac_pairs.csv")
	write_csv_safe(pd.DataFrame(rows), out_path)
	return out_path



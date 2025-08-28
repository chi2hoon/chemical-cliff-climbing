from __future__ import annotations

import os
from typing import Dict

import pandas as pd

from .io_csv import read_csv_safe, write_csv_safe


def _try_import_rdkit():
	try:
		from rdkit import Chem  # type: ignore
		from rdkit.Chem import AllChem  # noqa: F401
		from rdkit.Chem import inchi  # type: ignore
		return Chem, inchi
	except Exception as e:
		raise RuntimeError("RDKit is required for SMILES canonicalization. Use an environment with rdkit.")


def build_canonical_compounds(root_dir: str, cfg: Dict) -> str:
	"""
	data/silver/{dataset_id}/compounds_canonical.csv에 데이터셋별 표준화된 화합물 테이블 생성
	열: compound_id, dataset_id, smiles_raw, smiles_canonical, inchikey
	"""
	Chem, inchi = _try_import_rdkit()
	bronze_dir = os.path.join(root_dir, cfg["paths"]["bronze_dir"])
	comp_path = os.path.join(bronze_dir, "compounds.csv")
	if not os.path.exists(comp_path):
		raise FileNotFoundError(comp_path)
	df = read_csv_safe(comp_path)

	rows = []
	for _, r in df.iterrows():
		compound_id = str(r.get("compound_id", ""))
		smiles_raw = str(r.get("smiles_raw", ""))
		dataset_id = str(r.get("dataset_id", cfg.get("dataset_id", "")))
		smiles_canonical = ""
		inchikey = ""
		if smiles_raw:
			mol = Chem.MolFromSmiles(smiles_raw)
			if mol is not None:
				smiles_canonical = Chem.MolToSmiles(mol, canonical=True)
				try:
					inchikey = inchi.MolToInchiKey(mol)
				except Exception:
					inchikey = ""
		rows.append({
			"compound_id": compound_id,
			"dataset_id": dataset_id,
			"smiles_raw": smiles_raw,
			"smiles_canonical": smiles_canonical,
			"inchikey": inchikey,
		})

	silver_dir = os.path.join(root_dir, "silver", cfg.get("dataset_id", ""))
	os.makedirs(silver_dir, exist_ok=True)
	out_path = os.path.join(silver_dir, "compounds_canonical.csv")
	write_csv_safe(pd.DataFrame(rows), out_path)
	return out_path



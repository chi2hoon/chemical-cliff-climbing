import os

import pandas as pd

from .io_csv import read_csv_safe, write_csv_safe
from .parsers import normalize_unit_label


def _convert_unit(value: float | None, unit: str | None, target: str | None) -> float | None:
	if value is None:
		return None
	if unit is None or target is None or unit == target:
		return value
	# M, μM, nM 간 일반적인 변환
	conv = {
		("M", "μM"): 1e6,
		("M", "nM"): 1e9,
		("μM", "M"): 1e-6,
		("μM", "nM"): 1e3,
		("nM", "μM"): 1e-3,
		("nM", "M"): 1e-9,
	}
	f = conv.get((unit, target))
	if f is None:
		return value
	return value * f


def build_measurements_std(root_dir: str, cfg: dict) -> str:
	"""
	Create standardized measurements based on YAML config.
	- Preserve censor
	- unit_std decided from readout/matrix target units
	- value_std computed with unit conversion if needed
	Output: data/silver/{dataset_id}/measurements_std.csv
	"""
	paths = cfg["paths"]
	measurements_path = os.path.join(root_dir, paths["bronze_files"]["measurements"])  # using bronze as source
	if not os.path.exists(measurements_path):
		raise FileNotFoundError(measurements_path)
	df = read_csv_safe(measurements_path)

	# 소스에서 단위 정규화 보장
	df["unit"] = df.get("unit", "").map(normalize_unit_label)

	# 어세이 메타데이터 결합 (readout/matrix, 사용 가능한 경우)
	assays_path = os.path.join(root_dir, cfg["paths"]["bronze_files"]["assays"]) if "assays" in cfg["paths"]["bronze_files"] else None
	if assays_path and os.path.exists(assays_path):
		assays_df = read_csv_safe(assays_path)
		join_cols = [c for c in ["assay_id", "readout", "matrix", "label_raw"] if c in assays_df.columns]
		if "assay_id" in join_cols:
			df = df.merge(assays_df[join_cols], on="assay_id", how="left", suffixes=("", "_assay"))

	# 적합한 행: 숫자 판독값이 있는 효소/세포; 백분율 행 건너뛰기
	mask_numeric = df.get("unit", "").astype(str) != "%"
	mask_matrix = (df.get("matrix", "").isin(["cell", "enzyme"])) | df.get("assay_id", "").astype(str).str.startswith("assay.cell.")
	df_elig = df.loc[mask_numeric & mask_matrix].copy()

	# YAML readout/matrix 대상 단위를 기반으로 표준화
	target_units = cfg["rules"]["silver"].get("target_units", {})
	identity_for_cell_um = cfg["rules"]["silver"].get("identity_for_cell_um", False)
	identity_rules = cfg["rules"]["silver"].get("identity_rules", [])

	def decide_target_unit(row: pd.Series) -> str | None:
		readout = row.get("readout")
		matrix = row.get("matrix")
		if isinstance(readout, str) and readout in target_units:
			return target_units[readout]
		if matrix == "cell":
			return target_units.get("cell")
		if matrix == "enzyme":
			# 판독값이 제공되지 않는 한 일반적인 효소 기본값 없음
			return target_units.get(readout) if isinstance(readout, str) else None
		return None

	df_elig["unit_std"] = df_elig.apply(decide_target_unit, axis=1)

	# 가능한 단위 변환으로 value_std 계산; 검열 값 보존
	def to_float_safe(x: object) -> float | None:
		try:
			return float(x)
		except Exception:
			return None

	df_elig["value_std"] = [
		_convert_unit(to_float_safe(v), u, t)
		for v, u, t in zip(df_elig.get("value_num"), df_elig.get("unit"), df_elig.get("unit_std"))
	]
	df_elig["std_rule_id"] = cfg["rules"]["silver"]["transform"]["scale_rule_id"]
	df_elig["std_confidence"] = 1.0

	# 원본 필드 보존
	cols_keep = [
		"compound_id", "assay_id", "panel_id", "cell_line", "value_raw", "value_num", "unit", "censor",
		"provenance_file", "provenance_sheet", "provenance_row",
		"readout", "matrix",
		"unit_std", "value_std", "std_rule_id", "std_confidence",
	]
	cols_present = [c for c in cols_keep if c in df_elig.columns]
	# 새로 추가된 열이 존재하는지 확인
	for c in ["unit_std", "value_std", "std_rule_id", "std_confidence"]:
		if c not in cols_present:
			cols_present.append(c)

	silver_dir = os.path.join(root_dir, "silver", cfg.get("dataset_id", ""))
	os.makedirs(silver_dir, exist_ok=True)
	out_path = os.path.join(silver_dir, "measurements_std.csv")
	write_csv_safe(df_elig[cols_present], out_path)
	return out_path



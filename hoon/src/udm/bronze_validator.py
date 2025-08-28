from __future__ import annotations

import os
from typing import Dict, List, Tuple

import pandas as pd

from .io_csv import read_csv_safe, write_csv_safe
from .parsers import normalize_unit_label, split_censor_and_value, parse_scientific_notation

def _get_bronze_builder(cfg: Dict):
    dataset_id = str(cfg.get("dataset_id", ""))
    if dataset_id == "2017":
        from .bronze_build_2017 import build_bronze_from_raw as builder
        return builder
    else:
        # Lazy import to avoid importing 2017-incompatible builders
        from .bronze_build import build_bronze_from_raw as builder
        return builder


class BronzeValidator:
	def __init__(self, cfg: Dict):
		self.cfg = cfg
		self.paths = cfg["paths"]
		self.schema = cfg["schema"]
		self.rules = cfg["rules"]["bronze"]

	def _required_ok(self, df: pd.DataFrame, required_cols: List[str]) -> Tuple[bool, List[str]]:
		missing = [c for c in required_cols if c not in df.columns]
		return len(missing) == 0, missing

	def validate_csv(self, path: str, entity: str) -> Tuple[bool, List[str]]:
		df = read_csv_safe(path)
		required = self.schema[entity]["required"]
		ok, missing = self._required_ok(df, required)
		errors: List[str] = []
		if not ok:
			errors.append(f"Missing required columns in {entity}: {missing}")
		return len(errors) == 0, errors

	def _fix_measurements(self, df: pd.DataFrame, root_dir: str) -> pd.DataFrame:
		# 'unit'에서 단위 레이블 정규화
		if "unit" in df.columns:
			df["unit"] = df["unit"].map(normalize_unit_label)
		# 요청되고 value_num이 누락/유효하지 않은 경우에만 value_raw에서 censor/value 분리
		if self.rules.get("split_censor", True) and "value_raw" in df.columns:
			new_censor: List[str] = []
			new_value_num: List[str] = []
			for raw in df["value_raw"].tolist():
				censor, number_text, unit_token = split_censor_and_value(raw)
				if number_text is not None:
					val, audit = parse_scientific_notation(number_text)
					new_value_num.append(None if val is None else f"{val}")
				else:
					new_value_num.append(None)
				new_censor.append(censor)
			# 원본 bronze를 보존하기 위해 빈 곳만 채움
			if "censor" in df.columns:
				df["censor"] = df["censor"].where(df["censor"].astype(str).str.len() > 0, pd.Series(new_censor, index=df.index))
			else:
				df["censor"] = new_censor
			if "value_num" in df.columns:
				mask_empty = df["value_num"].astype(str).str.len() == 0
				df.loc[mask_empty, "value_num"] = pd.Series(new_value_num, index=df.index)[mask_empty]
			else:
				df["value_num"] = new_value_num
		# 절대 경로 포함을 피하기 위해 provenance_file을 데이터 루트에 상대적으로 만들기
		if "provenance_file" in df.columns and isinstance(root_dir, str) and len(root_dir) > 0:
			def _to_rel(p: object) -> object:
				try:
					sp = str(p)
					return os.path.relpath(sp, start=root_dir) if os.path.isabs(sp) else sp
				except Exception:
					return p
			df["provenance_file"] = df["provenance_file"].map(_to_rel)
		return df

	def build_or_validate(self, root_dir: str) -> Dict[str, List[str]]:
		"""Validate existing bronze CSVs under bronze/{year}; fix only deviations in measurements."""
		bronze_dir = os.path.join(root_dir, self.paths["bronze_dir"])
		files = self.paths["bronze_files"]
		results: Dict[str, List[str]] = {"errors": [], "fixed": []}

		# YAML에 assay_patterns가 존재하는 경우 새 규칙을 적용하기 위해 raw에서 재구축
		assay_patterns = self.cfg.get("parsing", {}).get("assay_patterns", [])
		if assay_patterns:
			_get_bronze_builder(self.cfg)(root_dir, self.cfg)
			results["fixed"].append(os.path.join(bronze_dir, "assays.csv"))
			results["fixed"].append(os.path.join(bronze_dir, "measurements.csv"))

		for entity, relpath in files.items():
			path = os.path.join(root_dir, relpath)
			if not os.path.exists(path):
				# 아티팩트가 누락된 경우 raw에서 한 번 구축 시도
				built = _get_bronze_builder(self.cfg)(root_dir, self.cfg)
				# 여전히 누락된 경우 오류 보고
				if not os.path.exists(path):
					results["errors"].append(f"Missing bronze file: {entity} at {path}")
					continue
			# 검증 전 최소 자동 수정: 누락된 경우 compounds에 dataset_id 추가
			if entity == "compounds":
				df_comp = read_csv_safe(path)
				if "dataset_id" not in df_comp.columns:
					df_comp["dataset_id"] = str(self.cfg.get("dataset_id", ""))
					write_csv_safe(df_comp, path)
					results["fixed"].append(path)
			ok, errs = self.validate_csv(path, entity)
			if not ok:
				results["errors"].extend(errs)
				continue
			# 측정값에 대해서만 규칙 준수를 위한 최소 수정 적용
			if entity == "measurements":
				df = read_csv_safe(path)
				df_fixed = self._fix_measurements(df.copy(), root_dir)
				if not df_fixed.equals(df):
					write_csv_safe(df_fixed, path)
					results["fixed"].append(path)
		return results



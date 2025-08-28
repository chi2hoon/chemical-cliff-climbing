import os
from typing import Any, Dict, Optional

import pandas as pd


def ensure_parent_dir(path: str) -> None:
	"""부모 디렉토리가 존재하지 않으면 생성합니다."""
	parent = os.path.dirname(os.path.abspath(path))
	if parent and not os.path.exists(parent):
		os.makedirs(parent, exist_ok=True)


def read_csv_safe(path: str, **kwargs: Any) -> pd.DataFrame:
	"""
	스키마 기반 수집을 위한 안전한 기본값으로 CSV를 읽습니다.

	기본값:
	- dtype=str
	- keep_default_na=False
	- na_filter=False
	"""
	read_kwargs: Dict[str, Any] = {
		"dtype": str,
		"keep_default_na": False,
		"na_filter": False,
	}
	read_kwargs.update(kwargs)
	return pd.read_csv(path, **read_kwargs)


def write_csv_safe(df: pd.DataFrame, path: str, **kwargs: Any) -> None:
	"""
	안전한 기본값으로 CSV를 작성합니다.

	기본값:
	- encoding="utf-8"
	- index=False
	- float_format="%.8g"
	"""
	ensure_parent_dir(path)
	write_kwargs: Dict[str, Any] = {
		"encoding": "utf-8",
		"index": False,
		"float_format": "%.8g",
	}
	write_kwargs.update(kwargs)
	df.to_csv(path, **write_kwargs)



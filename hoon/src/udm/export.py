from __future__ import annotations

import os
import shutil
from typing import Dict

import pandas as pd

from .io_csv import read_csv_safe, write_csv_safe


def _to_relative(path_value: str, root_dir: str) -> str:
    try:
        abs_root = os.path.abspath(root_dir)
        abs_val = os.path.abspath(path_value)
        if abs_val.startswith(abs_root):
            rel = os.path.relpath(abs_val, abs_root)
            return rel
    except Exception:
        pass
    return path_value


def export_bronze_release(root_dir: str, cfg: Dict, mode: str = "relative") -> str:
    """
    Create a sanitized export of bronze artifacts under data/bronze_release/{year}/.
    - For CSVs containing provenance_file, convert to root-relative paths (mode='relative').
    - Original bronze files remain untouched.
    Returns the release directory path.
    """
    bronze_dir = os.path.join(root_dir, cfg["paths"]["bronze_dir"])
    release_dir = os.path.join(root_dir, "bronze_release", cfg.get("dataset_id", ""))
    os.makedirs(release_dir, exist_ok=True)

    if not os.path.isdir(bronze_dir):
        raise FileNotFoundError(bronze_dir)

    for name in os.listdir(bronze_dir):
        src_path = os.path.join(bronze_dir, name)
        dst_path = os.path.join(release_dir, name)
        if not os.path.isfile(src_path):
            continue
        # provenance_file을 정리하기 위해 CSV 처리
        if name.lower().endswith('.csv'):
            df = read_csv_safe(src_path)
            if mode == "relative" and "provenance_file" in df.columns:
                df["provenance_file"] = df["provenance_file"].apply(lambda v: _to_relative(str(v), root_dir))
            write_csv_safe(df, dst_path)
        else:
            # 비-CSV 파일 (예: qc_overall.md)은 그대로 복사
            shutil.copy2(src_path, dst_path)

    return release_dir



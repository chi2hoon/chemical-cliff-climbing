from __future__ import annotations

from typing import Dict, Callable


def build_bronze_from_raw(root_dir: str, cfg: Dict) -> Dict[str, str]:
    """Dispatch to the appropriate year-specific bronze builder.

    - 2017: dedicated builder tuned for panel-aware parsing of tables 5–14.
    - Others: generic builder with year-specific regex and YAML pattern support.
    """
    dataset_id = str(cfg.get("dataset_id", ""))
    if dataset_id == "2017":
        from .bronze_build_2017 import build_bronze_from_raw as builder  # local import to avoid side-effects
    else:
        from .bronze_build import build_bronze_from_raw as builder  # type: ignore
    return builder(root_dir, cfg)


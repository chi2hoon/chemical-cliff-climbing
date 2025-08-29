from __future__ import annotations

from typing import Dict, Callable


def build_bronze_from_raw(root_dir: str, cfg: Dict) -> Dict[str, str]:
    """적절한 연도별 브론즈 빌더로 디스패치합니다.

    - 2017: 표 5-14의 패널 인식 파싱을 위해 조정된 전용 빌더.
    - 기타: 연도별 정규식과 YAML 패턴 지원이 포함된 일반 빌더.
    """
    dataset_id = str(cfg.get("dataset_id", ""))
    if dataset_id == "2017":
        from .bronze_build_2017 import build_bronze_from_raw as builder  # 부작용을 피하기 위한 로컬 임포트
    else:
        from .bronze_build import build_bronze_from_raw as builder  # type: ignore
    return builder(root_dir, cfg)


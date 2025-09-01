from importlib import import_module
import logging

# 연도별 특수 빌더 모듈 매핑
# 값은 "상대 모듈 경로"(현재 패키지 기준)로 둡니다.
_BUILDER_MODULE_BY_YEAR: dict[str, str] = {
    "2017": ".bronze_build_2017",
}

_DEFAULT_BUILDER_MODULE = ".bronze_build"  # 공통 빌더


def _import_builder(module_ref: str):
    """모듈에서 build_bronze_from_raw를 가져와 반환합니다.

    module_ref가 '.'로 시작하면 현재 패키지 기준의 상대 경로로, 아니면 절대 경로로 임포트합니다.
    """
    try:
        if module_ref.startswith("."):
            mod = import_module(module_ref, package=__package__)
        else:
            mod = import_module(module_ref)
    except Exception as e:
        raise ImportError(f"브론즈 빌더 모듈 임포트 실패: {module_ref}") from e

    try:
        return getattr(mod, "build_bronze_from_raw")
    except AttributeError as e:
        raise ImportError(
            f"모듈에 build_bronze_from_raw가 없습니다: {module_ref}"
        ) from e


def build_bronze_from_raw(root_dir: str, cfg: dict) -> dict[str, str]:
    """연도/설정에 맞는 브론즈 빌더로 디스패치합니다.

    우선순위:
      1) cfg["builder"]가 있으면 그 모듈을 사용
      2) 연도별 매핑(_BUILDER_MODULE_BY_YEAR)
      3) 기본 공통 빌더(_DEFAULT_BUILDER_MODULE)
    """
    dataset_id = str(cfg.get("dataset_id", "")).strip()
    if not dataset_id:
        raise ValueError("cfg['dataset_id']가 비어있습니다.")

    # YAML 등에서 직접 빌더를 지정할 수 있게 허용 (실험/오버라이드 용)
    override = str(cfg.get("builder", "")).strip()
    if override:
        module_ref = override
    else:
        module_ref = _BUILDER_MODULE_BY_YEAR.get(dataset_id, _DEFAULT_BUILDER_MODULE)

    builder = _import_builder(module_ref)

    logging.debug(
        "bronze builder selected",
        extra={
            "dataset_id": dataset_id,
            "module": getattr(builder, "__module__", str(module_ref)),
        },
    )

    return builder(root_dir, cfg)


import importlib

# 연도별 어댑터 레지스트리
# 값은 "모듈경로:함수" 형태로, 각 어댑터는 run(yaml_path, out_dir) 시그니처를 갖습니다.
ADAPTERS = {
    "2017": "pipeline.adapters.year_2017:run",
    "2018": "pipeline.adapters.year_2018:run",
    "2020": "pipeline.adapters.year_2020:run",
    "2021": "pipeline.adapters.year_2021:run",
}


def _resolve_callable(spec):
    """Args: spec(str) -> callable

    "모듈:함수" 문자열을 임포트하여 호출가능 객체를 반환한다.
    """
    if ":" not in spec:
        raise ValueError("어댑터 스펙은 'module.path:callable' 형식이어야 합니다.")
    mod_name, func_name = spec.split(":", 1)
    mod = importlib.import_module(mod_name)
    func = getattr(mod, func_name)
    return func


def build_bronze_from_raw(dataset_id, yaml_path, out_dir):
    """Args: dataset_id(str), yaml_path(str), out_dir(str) -> dict

    연도별 ADAPTERS 레지스트리를 사용해 해당 어댑터의 run을 호출한다.
    반환값은 산출물 경로 맵(예: {"compounds": path, ...}).
    """
    key = str(dataset_id)
    if key not in ADAPTERS:
        raise KeyError("지원하지 않는 연도: %s" % key)
    func = _resolve_callable(ADAPTERS[key])
    return func(yaml_path, out_dir)


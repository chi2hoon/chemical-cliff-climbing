from pipeline.adapters._yaml_ingest import run_yaml_bronze_ingest


def run(yaml_path, out_dir):
    """Args: yaml_path(str), out_dir(str) -> dict

    2017 브론즈를 공통 YAML 인제스트 엔진으로 수행한다.
    반환: 산출물 경로 맵(dict)
    """
    return run_yaml_bronze_ingest("2017", yaml_path, out_dir)

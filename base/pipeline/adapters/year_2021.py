def run(yaml_path, out_dir):
    """Args: yaml_path(str), out_dir(str) -> dict

    2021 브론즈는 공통 YAML 인제스트를 사용한다.
    """
    from pipeline.adapters._yaml_ingest import run_yaml_bronze_ingest
    return run_yaml_bronze_ingest(2021, yaml_path, out_dir)


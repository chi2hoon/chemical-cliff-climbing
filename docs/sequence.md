### Execution Flow: CLI → Bronze → Silver → Gold

#### 개요

이 문서는 연도별 특허 Excel에서 Gold 테이블까지 도달하는 실행 플로우를 시퀀스 다이어그램으로 정리합니다.

#### Medallion 파이프라인 시퀀스

```mermaid
sequenceDiagram
    participant User as User (CLI)
    participant CLI as pipeline.cli
    participant Adapter as YearAdapter (pipeline.adapters.year_*)
    participant Engine as BronzeEngine (_yaml_ingest)
    participant Silver as SilverBuilder (pipeline.silver)
    participant Gold as GoldBuilder (pipeline.gold)
    participant FS as Filesystem

    Note over User,Gold: 예시: 2017년 데이터를 Bronze→Silver→Gold로 재생성

    User->>CLI: python -m pipeline.cli bronze --year 2017
    CLI->>Adapter: build_bronze_from_raw(\"2017\", yaml, out_dir)
    Adapter->>Engine: run_yaml_bronze_ingest(\"2017\", yaml_path, out_dir)
    Engine->>FS: Read raw Excel (2017_raw.xlsx)
    Engine->>FS: Write data/bronze/2017/tables/*.csv + matrix/*_long.csv
    Engine->>FS: Write logs/manifest/2017.json + data/quarantine/bronze/*

    User->>CLI: python -m pipeline.cli silver --year 2017
    CLI->>Silver: build_silver(\"2017\", yaml)
    Silver->>FS: Read Bronze tables + YAML silver 설정
    Silver->>FS: Write data/silver/2017/compounds_silver.csv
    Silver->>FS: Write data/silver/2017/assay_readings_silver.csv
    Silver->>FS: Write data/silver/2017/assay_context_silver.csv
    Silver->>FS: Update logs/manifest/2017.json + data/quarantine/silver/*

    User->>CLI: python -m pipeline.cli gold --years 2017
    CLI->>Gold: build_gold([\"2017\"])
    Gold->>FS: Read data/silver/2017/*.csv
    Gold->>FS: Write data/gold/2017/assay_readings.csv
    Gold->>FS: Write data/gold/2017/compounds.csv
    Gold->>FS: Write data/gold/2017/compound_props.csv
    Gold->>FS: Write data/gold/2017/assay_context.csv
    Gold->>FS: Write data/quarantine/gold/*
```

#### 앱 레벨 데이터 플로우 (요약)

- `modules/io_utils.load_gold_data` 는 `data/gold/{year}/...`를 읽어 **SMILES + Activity** 스키마로 변환합니다.
- `app.py` 의 탭 1–2는 이 Gold 데이터를 사용해:
  - 연도/패널/타깃/세포주 선택
  - 구조-활성(Activity Cliff) 분석 및 히트맵 시각화를 수행합니다.
- 탭 3–6은 선택된 쌍과 Gold 메타데이터를 바탕으로 LLM API를 호출해:
  - 가설 생성 → 평가 → 수동/자동 수정 워크플로우를 제공합니다.



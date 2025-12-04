### Architecture: Bronze–Silver–Gold 파이프라인

이 문서는 특허 기반 데이터셋을 Medallion 스타일( **Bronze → Silver → Gold** )로 정리하는 파이프라인 구조를 정리합니다.

#### 개요

- **입력**: 연도별 특허 Excel (`data/bronze/artifacts/{year}_raw.xlsx`)
- **중간 산출물**
  - Bronze: 연도·시트별 슬라이스/롱 테이블 (`data/bronze/{year}/tables/*.csv`, `*_long.csv`)
  - Silver: 공통 스키마 기반 정제 테이블 (`data/silver/{year}/compounds_silver.csv`, `assay_readings_silver.csv`, `assay_context_silver.csv`)
- **최종 산출물(Gold)**  
  - `data/gold/{year}/assay_readings.csv`
  - `data/gold/{year}/compounds.csv`
  - `data/gold/{year}/compound_props.csv`
  - `data/gold/{year}/assay_context.csv`
- **품질 관리**
  - 격리(quarantine): `data/quarantine/{stage}/...`
  - 실행 메타데이터(manifest): `logs/manifest/{year}.json`

#### 계층별 역할

- **Bronze (원본 정리 / 슬라이싱)**
  - 연도별 어댑터: `pipeline/adapters/year_2017.py` 등
  - YAML 스키마: `schemas/silver/{year}.yaml`
  - 공통 인제스트 엔진: `pipeline/adapters/_yaml_ingest.py`
  - 역할:
    - 특허 Excel 시트/셀 범위를 YAML에 정의된 규칙으로 슬라이스
    - 헤더 보정 및 컬럼 리네이밍
    - 색상 기반 QC 플래그(`flag_asterisk`, `flag_imaging_conflict` 등) 주입
    - `provenance_file/sheet/row` 메타데이터 부여
    - 필요 시 매트릭스형 표를 롱 테이블로 전개 (`*_long.csv`)

- **Silver (표준화 / 단위 통일 / 맥락 분리)**
  - 엔트리: `pipeline/silver.py` 의 `build_silver`
  - 입력: Bronze 테이블 + YAML silver 섹션 (`schemas/silver/{year}.yaml`)
  - 역할:
    - 단위·기호 통일 (예: `>, >=` 기호 분리, `µM` → `uM` 정규화)
    - 정량 값(`value_std`, `unit_std`, `qualifier`)과 맥락(context)을 분리
    - 2018년 % 억제율, 2020년 asterisk 교정 등 연도별 예외 처리
    - 안정 정렬(`stable_sort`)을 통해 실행 시점과 무관한 재현성 확보
    - 은행 역할:
      - `compounds_silver.csv`: 구조/식별자/플래그 중심
      - `assay_readings_silver.csv`: 분석용 정량 측정값
      - `assay_context_silver.csv`: 도스, 치사율, % 억제율 등 보조 맥락

- **Gold (고정 스키마 / 모델 입력용)**  
  - 엔트리: `pipeline/gold.py` 의 `build_gold`
  - 입력: Silver 산출물
  - 역할:
    - `value_std`를 공통 단위(uM 등)로 통일
    - `derive_chem_ids`를 통해 `compound_key`, canonical SMILES, InChIKey14 등 생성
    - 최소 스키마 검증(`pipeline/validate/gold_checks.py`) 후 규칙 위반 레코드 격리
    - 4개 Gold 테이블 생성:
      - `assay_readings.csv`
      - `compounds.csv`
      - `compound_props.csv`
      - `assay_context.csv`

#### 품질 관리 / 격리 전략

- **격리(quarantine)**: 규칙 위반 또는 메타데이터 누락 레코드는 `data/quarantine/{stage}/...csv` 로 이동
  - Bronze: 프로비넌스 정보가 비어 있는 행
  - Gold: 필수 컬럼 누락, 규칙 위반 등
- **Manifest 로그**: `logs/manifest/{year}.json`
  - 사용된 YAML 경로 및 SHA256
  - 원본 Excel 경로 및 SHA256
  - 시트별 산출 행 수
  - 격리 파일 경로

이 구조를 통해, 연도별 특허 데이터의 난해한 시트 구조와 단위·표기 차이를 숨기고, **LLM/분석 앱이 바로 사용할 수 있는 안정된 Gold 스키마**로 수렴시키는 것이 목표입니다.



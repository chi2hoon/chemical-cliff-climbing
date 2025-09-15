## 실행 방법

아래 순서대로 그대로 따라 하면 됩니다. 

1) 가상환경 만들기 및 라이브러리 설치

"base/README.md"을 읽고, 가상환경을 생성한 이후를 가정합니다. 

2) 원본 엑셀 파일 복사(현재 2017만 사용)

아래 파일명을 정확히 맞춰 `base/data/bronze/artifacts/` 폴더에 두세요. 폴더가 없으면 만들어도 됩니다.

- 2017: `base/data/bronze/artifacts/2017_raw.xlsx`

주의: data/ 폴더는 Git에 올리지 않습니다(대용량/민감 데이터 제외 규칙).

3) 브론즈 → 실버 → 골드 실행(2017 예시)

```bash
# 모듈 경로 인식: PYTHONPATH=hoon 를 항상 붙입니다.

# 2017 Bronze: 시트→테이블/매트릭스 롱 (원본 보존)
PYTHONPATH=hoon python -m pipeline.cli bronze --year 2017 --cfg hoon/schemas/silver/2017.yaml

# 2017 Silver: 정규화/표준화(부등호/단위/도스/치사율/컨텍스트)
PYTHONPATH=hoon python -m pipeline.cli silver --year 2017 --cfg hoon/schemas/silver/2017.yaml

# 2017 Gold: 고정 스키마 + 메타(compound_props/assay_context)
PYTHONPATH=hoon python -m pipeline.cli gold --years 2017
```

4) (참고) 다른 연도는 추후 통합 예정

2018/2020/2021 스키마는 포함되어 있으나, 현재는 2017 데이터만 공식 지원/검증되었습니다.

5) 간단 검증(널바이트 스캔)

```bash
PYTHONPATH=hoon python -m pipeline.cli validate --stage nulls
# 출력이 NO_NULLBYTES 면 정상
```

산출물 위치(중요):
- Bronze: `base/data/bronze/{year}/tables/*.csv`, `base/data/bronze/{year}/matrix/*.csv`
- Silver: `base/data/silver/{year}/compounds_silver.csv`, `assay_readings_silver.csv`, `assay_context_silver.csv`
- Gold: `base/data/gold/{year}/compounds.csv`, `assay_readings.csv`, `compound_props.csv`, `assay_context.csv`

실패 시 가장 흔한 원인
- 원본 파일 경로/이름 오타(예: 2017_raw.xlsx)
- RDKit 설치 문제(pip 환경). conda-forge 설치를 권장합니다.
- Excel 엔진 누락: `openpyxl`이 설치되어야 합니다(base/requirements.txt에 포함됨).

# 파이프라인 안내(메달리온·YAML 설정 주도)

본 디렉토리는 4개 연도(2017/2018/2020/2021) 엑셀 원본을 메달리온 아키텍처(Bronze→Silver→Gold)로 처리하는 설정 주도형 파이프라인을 제공합니다. 모든 산출물은 CSV이며, 결정적 정렬과 출처(provenance)를 보존합니다.

## 디렉토리 구조
- pipeline/: 실행 가능한 CLI와 어댑터/파서/변환/검증 모듈
- schemas/silver/: 연도별 YAML 스키마(시트/범위/매트릭스/규칙/색상 플래그 등)
- data/
  - bronze/
    - artifacts/: 엑셀 원본(브론즈 원천)
    - {year}/: 시트→테이블/매트릭스 롱 CSV + manifest.json
  - silver/: Silver 표준화(compounds_silver.csv, assay_readings_silver.csv, assay_context_silver.csv)
  - gold/: Gold 고정 스키마(compounds.csv, assay_readings.csv, compound_props.csv, assay_context.csv)
  - quarantine/: 검증 실패/플래그 격리 CSV
- logs/manifest/: 실행 매니페스트(JSON)

## 원본 데이터 준비(참고)
- 파일 배치: `base/data/bronze/artifacts/{year}_raw.xlsx` (위 “실행 방법” 2) 단계 참고)
- YAML 설정: `hoon/schemas/silver/{year}.yaml` — 기본 제공 스키마의 `file:` 키가 해당 경로를 가리키도록 되어 있습니다.

## Gold 스키마(고정 계약)
- gold/{year}/compounds.csv
  - compound_key(=inchikey), smiles_canonical, has_structure, (옵션)iupac_name, inchikey14
- gold/{year}/assay_readings.csv
  - compound_id, target_id, assay_id, qualifier(=,>,<), value_std, unit_std(ASCII: uM,nM,%), year, qc_pass, provenance_file/sheet/row
 - 메타 CSV(분리 보존)
   - gold/{year}/assay_context.csv: 도스/성별/종, 치사율(n/k, % 등)
   - gold/{year}/compound_props.csv: 화합물 속성(iupac_name, MW, LC/MS, ¹H-NMR 등)

## Silver 표준화 규칙(요약)
- 부등호/값/단위 분리: qualifier/value_std/unit_std
- ASCII 단위 통일: µ/μ→u (uM/nM/% 등)
- 저장 전 안정 정렬: 자연 정렬 적용(예: 1,2,3, …, 10,11) — 키: compound_id/assay_id/provenance_row
- provenance 필수: 모든 행에 file/sheet/row 유지
- 색상 플래그: YAML flags_from_color로 읽어 flag_asterisk/flag_imaging_conflict 주입, 검증에서 격리
 - 날짜 오인 보정: `mortality_text`가 엑셀 직렬값/ISO 날짜로 들어와도 `1/10`, `10/10` 형태로 복원

## YAML 스키마 키(발췌)
- file, read: 엑셀 경로/옵션(dtype: str, engine: openpyxl)
- sheets[].tables[]: {id, range:A1, header_row_offset, rename, melt{id_cols,value_cols}}
- sheets[].matrix: {id, range, panel_row_offset, cellline_row_offset, data_row_start, id_col, panels[], value_parser{unit_default,rules[]}}
- flags_from_color: {col_index, map{RGB→flag_col}}
- silver.compounds: {from: data/ingest/{year}/tables/xxx.csv}
- silver.assays: {from_melt: name | from_matrix: {...}, variable_rules{label→{target_id, assay_id, parse{default_unit,rules[]}}}}

규칙 예시:
```
">=?\s*([\d\.]+)\s*µ?M -> qualifier:'>', value:\1, unit:uM"
"([\d\.]+)\s*µ?M -> qualifier:'=', value:\1, unit:uM"
```

## RDKit 경고 안내
- 예시 메시지: `Explicit valence for atom # 11 N, 4, is greater than permitted`
- 의미: 비정상 SMILES(과발렌스 등)로 RDKit sanitize 단계에서 경고가 발생
- 영향: 치명적이지 않음(파이프라인은 계속 진행). 해당 레코드는 `has_structure=False` 처리되거나 InChIKey 생성이 생략될 수 있음
- 대응: 원본 SMILES의 전하/결합 표기 등을 보정하면 경고가 사라짐. 필요 시 로거 레벨을 낮춰 경고 억제 가능

## 검증/격리
- silver_checks
  - compounds_silver: compound_id, provenance_* 필수
  - assay_readings_silver: compound_id, target_id, assay_id, qualifier, value_std, unit_std, year, provenance_* 필수
  - 규칙: unit_std∈{uM,nM,%}, qualifier∈{=,>,<}, value_std>0(있을 때)
  - 색상 플래그 True → data/quarantine/silver/{year}_assays_flagged.csv
- gold_checks
  - assay_readings: 단위/부등호/양수 확인
  - compounds: has_structure=True ⇒ inchikey(compound_key) 필수

모든 격리는 CSV로 저장하고, logs/manifest/{year}.json에 카운트를 기록합니다.

## 결정성/재현성
- 입력/스키마 해시를 manifest에 기록(input_sha256, yaml_sha256)
- 저장 전 항상 안정 정렬(동률 시 입력 순서 유지)

## 테스트
- 실행: `PYTHONPATH=hoon python -m pytest -q`
- 포함 테스트
  - normalize: parse_qual_value_unit, convert_unit
  - engine: matrix_to_long 스모크
  - silver: 2018/2020/2021 컬럼 준수, 2017 브릿지 산출 존재
  - gold: 컬럼 준수, has_structure=True ⇒ inchikey 존재

## 레거시 경고
- `hoon/archived/` 하위 레거시 코드는 사용하지 않습니다.
- 모든 신규 작업은 YAML 설정과 `pipeline/` 하위 모듈 + `PYTHONPATH=hoon python -m pipeline.cli ...`를 사용하십시오.

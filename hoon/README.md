# 파이프라인 안내(메달리온·YAML 설정 주도)

본 디렉토리는 4개 연도(2017/2018/2020/2021) 엑셀 원본을 메달리온 아키텍처(Bronze→Silver→Gold)로 처리하는 설정 주도형 파이프라인을 제공합니다. 모든 산출물은 CSV이며, 결정적 정렬과 출처(provenance)를 보존합니다.

## 디렉토리 구조
- pipeline/: 실행 가능한 CLI와 어댑터/파서/변환/검증 모듈
- schemas/silver/: 연도별 YAML 스키마(시트/범위/매트릭스/규칙/색상 플래그 등)
- data/
  - raw/: 엑셀 원본(Bronze 원천)
  - ingest/: Bronze 인제스트(시트→테이블/매트릭스 롱) + manifest.json
  - refined/: Silver 표준화(compounds_silver.csv, assay_readings_silver.csv)
  - gold/: Gold 고정 스키마(compounds.csv, assay_readings.csv)
  - quarantine/: 검증 실패/플래그 격리 CSV
- logs/manifest/: 실행 매니페스트(JSON)

## 실행 방법
```bash
# Bronze: 시트→테이블/매트릭스 롱(원본 보존)
python -m pipeline.cli bronze --year 2018 --cfg schemas/silver/2018.yaml

# Silver: 단위/부등호/텍스트 정규화(ASCII 단위, 결정적 정렬)
python -m pipeline.cli silver --year 2018 --cfg schemas/silver/2018.yaml

# Gold: 연도 합산 고정 스키마 생성
python -m pipeline.cli gold --years 2017 2018 2020 2021

# 유효성/환경 검사
python -m pipeline.cli validate --stage nulls   # 널바이트 스캔(NO_NULLBYTES 기대)
```

## Gold 스키마(고정 계약)
- gold/compounds.csv
  - compound_key(=inchikey), smiles_canonical, has_structure, (옵션)iupac_name, inchikey14
- gold/assay_readings.csv
  - compound_id, target_id, assay_id, qualifier(=,>,<), value_std, unit_std(ASCII: uM,nM,%), year, qc_pass, provenance_file/sheet/row

실험 맥락(도스/성별/종, LC/MS, ¹H-NMR 등)은 Gold 메타 CSV(추가 예정)로 분리합니다.

## Silver 표준화 규칙(요약)
- 부등호/값/단위 분리: qualifier/value_std/unit_std
- ASCII 단위 통일: µ/μ→u (uM/nM/% 등)
- 저장 전 안정 정렬: sort_values([compound_id, assay_id, provenance_row], kind="stable")
- provenance 필수: 모든 행에 file/sheet/row 유지
- 색상 플래그: YAML flags_from_color로 읽어 flag_asterisk/flag_imaging_conflict 주입, 검증에서 격리

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
- 실행: `python -m pytest -q`
- 포함 테스트
  - normalize: parse_qual_value_unit, convert_unit
  - engine: matrix_to_long 스모크
  - silver: 2018/2020/2021 컬럼 준수, 2017 브릿지 산출 존재
  - gold: 컬럼 준수, has_structure=True ⇒ inchikey 존재

## 레거시 경고
- `hoon/udm_cli.py`는 곧 폐기 예정입니다. 새 파이프라인 CLI(`python -m pipeline.cli ...`) 사용을 권장합니다.
  - 하위 호환을 위해 기존 경로는 유지하지만, 모든 신규 작업은 YAML 설정과 pipeline/ 하위 모듈을 사용하십시오.


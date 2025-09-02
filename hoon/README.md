# UDM Silver 데이터 안내

## 소개
- 본 레포지토리는 분석과 가설 생성을 위한 Silver 및 Gold 데이터를 제공합니다.
- 모든 수치는 표준 단위로 정규화되어 있으며, 원문 정보(provenance)와 센서(censor)는 보존되어 있습니다.

## Base 앱 연동 팁
- base에서 유사도/활성도차 계산 시 매핑은 smiles_col=smiles_canonical, activity_col=value_std로 사용하시면 됩니다.
- 의미 있는 비교를 위해 동일 `assay_id`(필요 시 `cell_line`) 그룹 내에서 분석하시길 권장드립니다.
- `censor`가 있는 경우, 절대차 계산에서 제외하거나 별도 처리 정책을 적용해 주십시오.

## 디렉토리 구조(제공물)
- `data/silver/{year}/`
  - `measurements_std.csv`
  - `compounds_canonical.csv`
- `data/gold/{year}/`
  - `measurements_gold.csv`
- 선택 제공: `data/silver/all/ac_pairs.csv`(활동 절벽 쌍)
- 주의: 2020의 `measurements_std.csv`는 등급형 지표 제외 정책으로 비어 있을 수 있습니다(헤더만 존재).

## Silver 데이터란?
- 의미: 원문을 손상하지 않고 단위·수치를 표준화한, 바로 분석 가능한 데이터 계층입니다.
- 핵심 컬럼
  - unit_std: 표준 단위(예: Ki는 nM, IC50은 μM)
  - value_std: 표준 단위로 변환된 수치
  - censor: 검출 한계/등호 정보(gt, lt, eq)
  - readout/matrix: 측정 종류/실험 매트릭스(Ki/enzyme, IC50/cell 등)
  - provenance_file/sheet/row: 원문 출처(상대경로)
- 구조 표준화: `compounds_canonical.csv`에는 `compound_id`, `smiles_canonical`, `inchikey`가 포함됩니다.

## 파일별 컬럼 개요
- measurements_std.csv
  - 필수: `compound_id`, `assay_id`, `value_raw`, `unit`, `censor`, `provenance_file`, `provenance_sheet`, `provenance_row`
  - 표준화: `readout`, `matrix`, `unit_std`, `value_std`, `std_rule_id`, `std_confidence`
  - 선택적 맥락: `cell_line`(해당 시)
- compounds_canonical.csv
  - `compound_id`, `dataset_id`, `smiles_raw`, `smiles_canonical`, `inchikey`

## Silver 데이터 활용: 가설 생성 가이드
- 최소 입력
  - `data/silver/{year}/measurements_std.csv`
  - `data/silver/{year}/compounds_canonical.csv`
- 권장 조인
  - `compound_id`로 조인하여 각 측정값에 캐노니컬 SMILES를 연결해 주십시오.
- 분석 팁
  - 동일 맥락 비교: `assay_id`(필요 시 `cell_line`) 그룹 내 비교를 권장드립니다.
  - 단위 확인: `unit_std` 기준으로 `value_std`를 비교하시면 안전합니다.
  - 센서 처리: `gt/lt`는 검열값으로 해석하시고, 회귀·순위 학습에서 적절히 반영해 주십시오.
- 모델 아이디어(예시)
  - 유사 구조 대비 성능 차이를 학습하는 Pairwise 랭킹
  - SMILES + 실험 메타(`readout`, `matrix`, `cell_line`) 기반 다과제 회귀/분류

## Activity Cliff 기준(유사도·활성 차이)
- 정의: 동일 `assay_id`(및 필요 시 `cell_line`) 그룹 내에서
  - 화학 유사도: Tanimoto(Morgan fingerprint, radius=2, 2048 bit) ≥ 0.85
  - 활성 차이: |value_std_i − value_std_j| ≥ 1.0
- 단위 해석: 활성 차이는 표준 단위 기준입니다. 예를 들어 Ki 세트는 nM, IC50(세포) 세트는 μM 단위에서 1.0 차이를 의미합니다.
- 선정 이유
  - 유사도 0.85는 활동 절벽 탐지에 널리 쓰이는 보수적 기준입니다.
  - 활성 차이 1.0은 측정 잡음을 넘어서는 실질적인 효능 차이를 보장합니다.
- 참고: 탐색 범위를 넓히고자 하시면(가설 확대) 0.80/0.5 등으로 완화하는 것도 가능합니다.

## 데이터 원칙(요약)
- 추정 금지: No guessing/imputation
- 센서 분리 유지: `censor`는 항상 분리 보존합니다(>, <, =).
- 출처 보존: 모든 레코드에 `provenance_file/sheet/row`를 유지합니다(상대경로).
- 형식: CSV만 사용하며, 스키마/규칙은 YAML로 관리합니다.

## 협업/배포
- 본 레포지토리는 Silver 데이터만 제공합니다.
- LLM 가설 생성·분석 파이프라인의 입력으로 `measurements_std.csv`와 `compounds_canonical.csv`를 사용해 주십시오.

## 파이프라인 개요
- 스테이지
  - bronze: 원천 데이터 수집/검증
  - silver: 단위 표준화(`value_std`, `unit_std`, `censor` 유지)
  - smiles: 구조 표준화(`smiles_canonical`, `inchikey` 생성)
  - gold: silver 측정값과 캐노니컬 구조를 조인해 분석 친화 테이블 생성
  - ac / ac-all: 사전 계산된 Activity Cliff 쌍 산출(`data/silver/all/ac_pairs.csv`)
  - export: bronze 릴리스 산출물 내보내기

### 실행 순서(예: 2017)
```bash
# 1) Silver: 측정값 표준화
python hoon/udm_cli.py silver --config hoon/configs/2017.yml --root hoon/data

# 2) SMILES: 캐노니컬 구조/인치키
python hoon/udm_cli.py smiles --config hoon/configs/2017.yml --root hoon/data

# 3) Gold: 분석 친화 테이블(조인 결과)
python hoon/udm_cli.py gold --config hoon/configs/2017.yml --root hoon/data

# (선택) Activity Cliff 사전 계산
python hoon/udm_cli.py ac --config hoon/configs/2017.yml --root hoon/data
```

## Gold 데이터란?
- 경로: `data/gold/{year}/measurements_gold.csv`
- 목적: 유사도/활성도 차이 계산과 모델링에 바로 사용할 수 있도록 표준 측정값과 캐노니컬 구조를 한 테이블로 제공합니다.
- 주요 컬럼(가능 시)
  - `smiles_canonical`, `inchikey`
  - `compound_id`, `dataset_id`
  - `assay_id`, `panel_id`, `cell_line`
  - `value_std`(float), `unit_std`, `censor`
  - `readout`, `matrix`
  - `provenance_file`, `provenance_sheet`, `provenance_row`
  - `std_rule_id`, `std_confidence`
- 필터: `smiles_canonical`이 비어있지 않고, `value_std`가 숫자로 변환 가능한 행만 포함합니다.


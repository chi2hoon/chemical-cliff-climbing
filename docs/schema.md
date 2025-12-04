### Gold 스키마 개요

이 문서는 Medallion 파이프라인의 최종 산출물인 **Gold 레이어**의 테이블 구조를 정리합니다.  
모든 파일은 `data/gold/{year}/` 하위에 생성됩니다.

---

### 1. `assay_readings.csv`

단일 측정값(예: IC₅₀, Ki 등)을 한 행으로 나타내는 **정량 측정 테이블**입니다.

- **주요 컬럼**
  - `compound_id` : 특허 표에서 사용된 화합물 ID (예: `73`, `Example 4`)
  - `target_id` : 타깃 또는 세포주/패널 ID (예: `cell:방광암세포주.KU-19-19`, `prmt5.enzyme`)
  - `assay_id` : 실험 종류 (예: `cell.cytotoxicity.ic50`, `prmt5.enzyme.ki.A`)
  - `qualifier` : 부등호/검열 기호 (`=`, `>`, `>=` 등)
  - `value_std` : 숫자 값 (float, 공통 단위 기준)
  - `unit_std` : 표준화된 단위 (예: `uM`)
  - `year` : 특허/데이터 연도 (`2017`, `2018`, `2020`, `2021`)
  - `qc_pass` : 품질 플래그 (기본 True, 규칙 위반 시 quarantine로 이동)
  - `provenance_file`, `provenance_sheet`, `provenance_row` : 원본 Excel 위치
  - `compound_key` : 구조 기반 고유 ID (InChIKey14를 이용해 `pipeline.transforms.chem_ids`에서 생성)

- **개념적 키**
  - ( `compound_key`, `target_id`, `assay_id`, `year`, `provenance_row` ) 조합이 한 측정값을 유일하게 식별

---

### 2. `compounds.csv`

구조 ID와 표준화된 구조 정보를 담는 **화합물 마스터 테이블**입니다.

- **컬럼**
  - `compound_key` : 구조 기반 고유 키 (PK)
  - `smiles_canonical` : canonical SMILES
  - `has_structure` : 구조가 유효하게 파싱되었는지 (`True/False`)
  - `iupac_name` : IUPAC 이름 (있을 경우)
  - `inchikey14` : InChIKey의 앞 14자

- **관계**
  - `assay_readings.compound_key` → `compounds.compound_key` (FK)
  - `compound_props.compound_key` → `compounds.compound_key` (FK)

---

### 3. `compound_props.csv`

유효 구조에 대한 **실험/문헌 기반 속성 및 QC 플래그**를 묶은 테이블입니다.

- **주요 컬럼**
  - `compound_key` : 구조 ID (FK → `compounds.compound_key`)
  - `compound_id` : 특허 상 화합물 번호
  - `iupac_name` : IUPAC 이름
  - `mw` : 분자량
  - `lcms_text` : LC/MS 텍스트
  - `nmr_1h_text` : ¹H NMR 텍스트
  - `flag_asterisk` : 이미지–SMILES 교정(a* 예외)이 적용되었는지 여부
  - `flag_imaging_conflict` : 이미지/텍스트 간 충돌 여부
  - `flag_smiles_o3_changed` : SMILES 변경 플래그
  - `year` : 연도
  - `provenance_*` : 원본 위치

- **개념적 키**
  - ( `compound_key`, `compound_id`, `year`, `provenance_row` )

---

### 4. `assay_context.csv`

정량 값 주변의 **보조 맥락 정보**를 연도·특허별로 다르게 수용하는 테이블입니다.

- 예시 컬럼 (연도별로 부분 집합)
  - 2017
    - `species`, `sex`
    - `dose_route`, `dose_value`, `dose_unit`
    - `mortality_pct`, `mortality_n`, `mortality_k`, `mortality_ratio`
    - `hr_change_pct`
  - 2018
    - `percent_at_20uM`
  - 2020
    - asterisk 교정/ADME 등급 연결에 필요한 메타 필드
  - 2021
    - `sample_id`, PRMT5 Ki A/B, 세포 증식 IC₅₀와 매핑되는 등급/부가 정보

- **키/관계**
  - 공통적으로 `compound_id` 또는 `sample_id`를 통해 `assay_readings` 와 느슨하게 연결됩니다.
  - 강한 FK 제약 대신, 분석 시점에 조인하는 **맥락 테이블**로 사용하는 것을 전제로 설계되었습니다.

---

### 5. SQL 매핑 & 예시 쿼리

Gold는 그대로 **데이터 웨어하우스 테이블**로 적재하기 쉽게 설계되었습니다.  
아래는 PostgreSQL/BigQuery 스타일을 가정한 예시입니다.

```sql
-- 1) 타깃별 상위 활성 화합물 랭킹 (pAct 관점)
SELECT
  ar.target_id,
  ar.compound_key,
  c.smiles_canonical AS smiles,
  MIN(ar.value_std)      AS best_ic50_uM,
  -LOG(10, MIN(ar.value_std * 1e-6)) AS best_pAct  -- 단순 예시
FROM assay_readings AS ar
JOIN compounds      AS c
  ON ar.compound_key = c.compound_key
WHERE ar.assay_id = 'cell.cytotoxicity.ic50'
  AND ar.unit_std = 'uM'
GROUP BY ar.target_id, ar.compound_key, c.smiles_canonical
ORDER BY ar.target_id, best_pAct DESC
LIMIT 100;
```

```sql
-- 2) 구조–활성 Cliff 후보 추출(같은 타깃, 서로 다른 compound_key 쌍)
-- (실제 앱에서는 RDKit/Tanimoto를 사용하지만, 여기서는 값 차이만 예시로 제시)
WITH per_compound AS (
  SELECT
    compound_key,
    target_id,
    MIN(value_std) AS best_ic50_uM
  FROM assay_readings
  WHERE assay_id = 'cell.cytotoxicity.ic50'
    AND unit_std = 'uM'
  GROUP BY compound_key, target_id
)
SELECT
  a.target_id,
  a.compound_key AS compound_key_1,
  b.compound_key AS compound_key_2,
  a.best_ic50_uM,
  b.best_ic50_uM,
  ABS(LOG(10, a.best_ic50_uM) - LOG(10, b.best_ic50_uM)) AS delta_pAct
FROM per_compound AS a
JOIN per_compound AS b
  ON a.target_id = b.target_id
 AND a.compound_key < b.compound_key
WHERE ABS(LOG(10, a.best_ic50_uM) - LOG(10, b.best_ic50_uM)) >= 2.0;
```

이처럼 Gold 레이어는:

- `compounds` 로 **구조 마스터**를 정의하고,
- `assay_readings` 로 **정량 활성값**을 표준 단위로 모으며,
- `compound_props` 와 `assay_context` 로 **QC·독성·ADME 맥락**을 추가해,

LLM/모델·대시보드·웨어하우스 적재 모두에서 재사용 가능한 형태로 설계되어 있습니다.



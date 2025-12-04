## TL;DR

- **무엇**: 특허 기반 신약 후보 데이터(4개 연도)를 Medallion 스타일 **Bronze–Silver–Gold 파이프라인**으로 정리하고, Gold 레이어를 기반으로 **Activity Cliff 분석 + LLM 가설 생성**까지 이어지는 프로토타입입니다.
- **내 역할**: `/pipeline` 전반(어댑터 레지스트리, YAML 기반 Bronze 인제스트, Silver/Gold 표준화, QC/격리/로그)을 설계·구현하고, 앱이 바로 쓸 수 있는 `data/gold/*` 스키마를 정의했습니다.
- **결과**: 연도·특허별로 제각각이던 표/단위/기호를 **고정 스키마(`assay_readings`, `compounds`, `compound_props`, `assay_context`)**로 수렴시키고, 격리/manifest로 **재현성과 추적 가능성**을 확보했습니다.
- **검증 범위**
  - **Level 1 (키 없이)**: Bronze→Silver→Gold 파이프라인 재실행 및 Gold 스키마/테스트 검증
  - **Level 2 (OpenAI API 키 필요)**: OpenAI API Key 설정 후 Activity Cliff 기반 **가설 생성·평가·수정 LLM 플로우**
- **Links**
  - YouTube 데모: `https://www.youtube.com/watch?v=q2CVhYdBhU4`(2:28 ~ 8:50 까지 제 파트입니다)
  - 현재 레포: `https://github.com/chi2hoon/chemical-cliff-climbing`
  - 파이프라인 엔트리: [`pipeline/cli.py`](pipeline/cli.py), [`pipeline/silver.py`](pipeline/silver.py), [`pipeline/gold.py`](pipeline/gold.py)
---

## 프로젝트 한 줄 소개 & 목적

- **한 줄 소개**: 특허 PDF에서 추출한 난해한 Excel들을 정제·표준화하여, 연구원이 **구조–활성 관계(SAR)를 탐색하고 LLM 기반 가설을 빠르게 생성**할 수 있게 돕는 데이터/AI 파이프라인입니다.
- **목적**
  - 신약 개발 전체 5단계 중 **3단계(리드 최적화)** 구간에서, 구조 변경에 따른 활성/독성 가설을 더 빨리 만들고 검증 후보를 줄이는 것에 초점을 맞췄습니다.
  - 이 레포에서는 그중에서도 **데이터 엔지니어링 관점**(Medallion 파이프라인, QC, 재현성)을 중심으로 정리해 두었습니다.

---

## 데이터 문제 (원본 특허 → Excel)

원본 데이터는 Tabula 등으로 특허 PDF를 Excel로 추출한 연도별 4개 데이터셋입니다. 주요 문제는 다음 네 가지입니다.

- **단위 / 표기 불일치**
  - 동일 지표라도 `µM`, `nM`, `%`, `+, ++, +++` 등 **단위·등급·부등호 표기**가 제각각입니다.
  - 예: `>= 2 µM`, `++++`, `45%` 등은 그대로는 하나의 수치 축으로 비교가 어렵습니다.

- **시트 구조 난해**
  - 하나의 워크북에 **다중 테이블**, 좌우 병렬 배치, 헤더 누락/오탈자가 섞여 있습니다.
  - 시작 행을 수동으로 지정할 수 없고, 연도별로 시트 구조가 달라 **일관된 파서가 필요**했습니다.

- **정량 / 정성 혼재**
  - 동일 표 안에 IC₅₀ 수치, “0/10”, “전 개체 사망”, 심박수 변화 같은 정성/로그가 함께 들어 있습니다.
  - 모델/분석은 **정량 값만 정제**해서 쓰고, 나머지는 **맥락(context) 메타데이터**로 분리해야 했습니다.

- **변환 오류 / QC 흔적**
  - PDF→Excel 과정에서 이미지가 문자로 오인되거나, 수동 QC 히스토리(색상 플래그, 주석)가 남아 있습니다.
  - 이 흔적을 그대로 쓰면 나중에 구조–활성 해석에서 **혼선과 편향**을 유발할 수 있어, 플래그로 분리·격리했습니다.

---

## 제가 한 일

이 레포는 팀 프로젝트 레포를 포크한 것이고, 저는 주로 **데이터 아키텍처와 파이프라인 구현**을 맡았습니다.  
아래 파일들이 제가 실제로 작업한 부분들입니다.

- **파이프라인 엔트리 & CLI**
  - `pipeline/cli.py`: 연도별 Bronze/Silver 실행과 Gold/validate 서브커맨드를 제공하는 **Medallion CLI 진입점**입니다.
  - `pipeline/bronze_dispatch.py`: 연도별 어댑터 레지스트리(`ADAPTERS`)와 동적 import 로직을 통해, 새 연도 추가 시 **레지스트리 한 줄로 확장** 가능하게 설계했습니다.

- **YAML 기반 Bronze 인제스트**
  - `pipeline/adapters/_yaml_ingest.py`: YAML 스키마(`schemas/silver/{year}.yaml`)를 읽어
    - 시트/범위 슬라이스
    - 헤더 행 자동 적용
    - 색상 기반 QC 플래그(`flag_asterisk`, `flag_imaging_conflict`) 추출
    - `provenance_*` 메타데이터 주입
    를 공통 엔진으로 처리합니다.
  - `schemas/silver/2017.yaml` 등: 연도별 **시트 구조/예외를 코드가 아닌 설정(YAML)에 캡슐화**했습니다.

- **Silver 표준화 레이어**
  - `pipeline/silver.py`: 단위/기호 통일, 퍼센트/등급 파싱, 정량(readings) vs 맥락(context) 분리를 담당합니다.
    - `parse_qual_value_unit`, `convert_unit` 등을 이용해 `value_std`, `unit_std`, `qualifier` 스키마로 통일
    - 2018 % 억제율, 2020 asterisk 교정, 도스/치사율 파싱 등 **연도별 예외 규칙을 Silver 단계에서 흡수**

- **Gold 스키마 & QC**
  - `pipeline/gold.py`: Silver 출력을 읽어
    - `assay_readings.csv`, `compounds.csv`, `compound_props.csv`, `assay_context.csv` 4개 Gold 테이블을 생성
    - `pipeline.transforms.chem_ids.derive_chem_ids` 를 통해 `compound_key`, canonical SMILES, InChIKey14를 유도
    - `pipeline.validate.gold_checks` 와 `write_quarantine` 로 규칙 위반 레코드를 **data/quarantine/gold/** 로 분리

- **정렬/격리/manifest 유틸**
  - `pipeline/validate/hooks.py`: 안정 정렬(`stable_sort`), 격리 CSV 기록(`write_quarantine`), 실행 매니페스트(`write_manifest`)를 구현해, **재실행 시에도 동일한 정렬/경로를 보장**합니다.

- **앱과 Gold를 잇는 인터페이스**
  - `modules/io_utils.py`: `data/gold/{year}/*`를 읽어 **SMILES + Activity 스키마**로 변환하고, 패널/타깃/세포주 별로 필터링하는 로직을 구현했습니다.
  - `modules/context_builder.py`: Gold/ Silver/ context를 조합해 LLM 가설 생성에 필요한 **구조·활성·독성/ADME 메타 JSON**을 만들어 줍니다.

---

## Architecture (Bronze–Silver–Gold)

상세 설명은 `docs/architecture.md`와 `docs/sequence.md`를 참고하세요. 아래는 요약 다이어그램입니다.

```mermaid
flowchart LR
    Raw[Patent Excel (2017_raw.xlsx 등)] -->|YAML 스키마: schemas/silver/{year}.yaml| Bronze[Bronze: YAML Ingest Engine]
    Bronze -->|tables/*.csv, *_long.csv| Silver[Silver: 정제/표준화]
    Silver -->|compounds_silver, assay_readings_silver, assay_context_silver| Gold[Gold: 고정 스키마]

    Bronze -->|manifest.json| Logs[Logs: logs/manifest/{year}.json]
    Silver -->|manifest.json 업데이트| Logs
    Bronze -->|이상값/메타 누락| QBronze[Quarantine: data/quarantine/bronze/*]
    Gold -->|규칙 위반 레코드| QGold[Quarantine: data/quarantine/gold/*]
```

- **Bronze**
  - 연도별 YAML + 공통 엔진으로 특허 Excel을 슬라이스하고, 시트 구조/색상 플래그를 처리합니다.
- **Silver**
  - 정량 값/맥락 분리, 단위/부등호 통일, 연도별 예외 보정(asterisk, % 억제율 등)을 수행합니다.
- **Gold**
  - 구조 ID(`compound_key`)와 고정 컬럼 셋을 갖는 테이블로 수렴시키고, QC 규칙 위반 레코드를 격리합니다.

---

## Results (관측 가능한 결과)

- **Gold 산출물**
  - 연도별 공통 구조:
    - `data/gold/{year}/assay_readings.csv`
    - `data/gold/{year}/compounds.csv`
    - `data/gold/{year}/compound_props.csv`
    - `data/gold/{year}/assay_context.csv`
  - 각 파일의 스키마 및 관계는 `docs/schema.md` 에 정리되어 있습니다.

- **정규화/표준화 효과**
  - IC₅₀/Ki 등 농도 값은 **단일 단위(uM)** 기준으로 `value_std` + `unit_std` 컬럼에 정규화됩니다.
  - `> 2 µM`, `++++` 같은 표기는 `qualifier` (`>`, `=` 등)와 수치/등급 컬럼으로 분해되어, **수치 비교와 검열값 처리를 분리**할 수 있습니다.
  - 퍼센트 억제율, 도스/치사율 등은 `assay_context_silver.csv` 및 Gold `assay_context.csv`로 흡수되어, 정량 테이블과 분리됩니다.

- **재현성/품질 검증**
  - 동일한 YAML/원본 Excel에 대해 `python -m pipeline.cli ...`를 여러 번 실행하면, `stable_sort` 와 manifest를 통해 **동일한 CSV/로그를 재생성**할 수 있습니다.
  - 필수 메타데이터(provenance) 누락이나 Gold 규칙 위반 레코드는 `data/quarantine/*`로 격리되어, 후속 분석에서 **눈에 잘 띄는 형태로 분리**됩니다.

---

## Reproducibility (재현 방법)

### 공통 환경

- OS: macOS (zsh 기준)
- Python: 3.10 이상 권장
- 필수 패키지: `requirements.txt` 참고

```bash
# 레포 클론 후(또는 이 브랜치 체크아웃 후) 루트에서 실행
python3 -m venv .venv
source .venv/bin/activate  # macOS / zsh
pip install -r requirements.txt
```

### Windows (PowerShell 예시)

```powershell
cd C:\path\to\chemical-cliff-climbing
python -m venv .venv
.\.venv\Scripts\Activate.ps1
pip install -r requirements.txt
```

### Level 1 — 키 없이 검증 가능한 범위

**목표**: OpenAI 키 없이도 **파이프라인과 Gold 스키마**를 검증할 수 있는 최소 경로입니다.

```bash
# 1) 2017년 데이터에 대해 Bronze → Silver → Gold 재생성
python -m pipeline.cli bronze --year 2017
python -m pipeline.cli silver --year 2017
python -m pipeline.cli gold --years 2017

# 2) Silver/Gold 스키마 검증 테스트 실행 (선택)
pytest pipeline/tests/test_silver_outputs.py pipeline/tests/test_gold_outputs.py
```

- 위 명령만으로:
  - `data/bronze/2017`, `data/silver/2017`, `data/gold/2017` 가 재생성됩니다.
  - 테스트는 2017년 Silver/Gold 산출물의 **필수 컬럼이 모두 존재하는지**를 확인합니다.

### Level 2 — OpenAI API Key가 있을 때

**목표**: Gold 데이터를 사용한 Activity Cliff 분석과, LLM 기반 **가설 생성·평가·수정 탭**까지 end-to-end로 확인합니다.

```bash
# 1) (권장) 환경변수로 API Key 설정 — 키 값은 절대 커밋하지 마세요.
export OPENAI_API_KEY="sk-..."   # 사용자가 직접 발급한 키를 입력

# 2) Streamlit 앱 실행
source .venv/bin/activate
streamlit run app.py
```

- **탭 1–2 (키 없이도 동작)**: Gold 디렉토리에서 연도/패널/타깃을 선택하고, Activity Cliff 히트맵과 구조–활성 분포를 탐색합니다.
- **탭 3–6 (키 필요)**: 선택된 cliff 쌍에 대해
  - LLM이 구조–활성 가설을 JSON 스키마로 생성
  - 평가/수정 루프를 통해 가설을 강화
  - 결과를 Markdown 파일로 저장 (`hypotheses/`, `evaluations/`)

---

## Known Issues / Limitations & Next Steps

- **한계 / 제약**
  - 원본 특허 데이터는 4개 연도·특허 세트에 한정되어 있으며, **새로운 특허 구조에는 YAML 스키마 추가 작업**이 필요합니다.
  - 단위/등급 파싱은 현재 규칙 집합에 맞춰져 있어, 완전히 새로운 표기법이 들어오면 수동 보정이 필요할 수 있습니다.
  - Gold 검증(`pipeline/validate/gold_checks.py`)은 **최소 스키마 규칙**만을 적용하며, 도메인별 복잡한 통계적 QC까지는 다루지 않습니다.
  - LLM 기반 가설 생성/평가는 OpenAI API에 의존하므로, **과금·Rate limit·모델 업데이트**에 따라 동작이 달라질 수 있습니다.

- **다음 개선 아이디어**
  - 연도/특허 추가 시 공통으로 재사용 가능한 **YAML 스키마 템플릿/체커**를 추가해, 설정 오류를 사전에 방지.
  - Silver/Gold 단계에 **통계적 QC(분포/아웃라이어 탐지)** 및 시각화 리포트 자동 생성.
  - LLM 플로우를 벤더 중립 인터페이스로 추상화해, **오픈소스 LLM·사내 모델**로도 쉽게 바꿔 끼울 수 있도록 확장.

### Security / 공개 안전성

- 이 레포에는 **API Key, 토큰, 내부 URL, 개인 데이터 경로가 포함되어 있지 않습니다.**
  - `.gitignore` 에서 `.env`, `openAI_key.txt`, `base/data` 원본, 대용량 산출물 등이 제외되어 있습니다.
  - OpenAI 키는 **환경변수(`OPENAI_API_KEY`) 또는 로컬 `openAI_key.txt`** 에서만 읽으며, 레포에 커밋하지 않습니다.
- 사용자는 반드시 **본인 계정에서 직접 발급한 키**를 사용해야 하며, 키 값은 Pull Request나 이슈 등에 공유하지 않는 것을 전제로 합니다.

---

## SQL 매핑 (Gold → Warehouse)

Gold 레이어는 그대로 데이터 웨어하우스(예: BigQuery, Snowflake) 테이블로 적재하기 쉽게 설계되어 있습니다.  
각 테이블과 관계는 `docs/schema.md` 에 정리되어 있으며, 핵심 관계는 다음과 같습니다.

- `compounds(compound_key)` : 구조 마스터 (PK)
- `assay_readings(compound_key, target_id, assay_id, year, provenance_row, ...)` : 정량 측정 (FK → `compounds`)
- `compound_props(compound_key, compound_id, ...)` : 구조별 속성 (FK → `compounds`)
- `assay_context` : 도스/독성/등급 등 맥락 정보 (느슨한 키 기반 조인)

간단한 예시 쿼리는 다음과 같습니다 (자세한 예시는 `docs/schema.md` 참고).

```sql
-- 타깃별 상위 활성 화합물 (개념 예시)
SELECT
  ar.target_id,
  ar.compound_key,
  c.smiles_canonical AS smiles,
  MIN(ar.value_std) AS best_ic50_uM
FROM assay_readings AS ar
JOIN compounds      AS c
  ON ar.compound_key = c.compound_key
WHERE ar.assay_id = 'cell.cytotoxicity.ic50'
  AND ar.unit_std = 'uM'
GROUP BY ar.target_id, ar.compound_key, c.smiles_canonical
ORDER BY ar.target_id, best_ic50_uM ASC
LIMIT 100;
```

---

## Attribution (팀 프로젝트 & 포크 출처)

- 이 레포는 팀 프로젝트 레포를 포크한 것으로, **원본 레포의 전체 아이디어/LLM 프롬프트/앱 UI는 팀 작업 결과**입니다.
- 원본 팀 레포: `https://github.com/nazirite96/chemical-cliff-climbing`
- 이 포트폴리오 브랜치에서는 다음에 초점을 맞췄습니다.
  - `/pipeline` 및 관련 YAML/validate 모듈 정리
  - Medallion 아키텍처 및 Gold 스키마 문서화 (`docs/architecture.md`, `docs/sequence.md`, `docs/schema.md`)
  - 포트폴리오 관점에서 **데이터 파이프라인/재현성/QC** 를 이해하기 쉽게 보여주는 README 개편


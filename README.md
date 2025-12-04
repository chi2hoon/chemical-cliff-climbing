# Chemical Cliff Climbing

> #### 제각각인 실험 Excel → 깨끗한 표준 테이블 → AI가 구조–활성 가설 자동 생성

---

## TL;DR — 이 프로젝트를 한 줄로 말하면

**제각각인 실험 데이터를, 모델·분석이 바로 쓸 수 있는 표준 테이블로 바꾸는 데이터 파이프라인입니다.**  
도메인은 신약이지만, 제가 보여드리고 싶은 건 **데이터/ML 엔지니어링 역량**입니다.

### 무엇을 만들었나요

특허 PDF에서 뽑아낸 연도별 실험표를, 규칙 기반 파이프라인(Bronze→Silver→Gold)으로 정제해서  
**깨끗한 CSV 테이블 4개(Gold 레이어)** 로 만들었습니다.

### 제가 한 일

`/pipeline` 폴더 전반을 담당했습니다.

- 연도별 어댑터 설계
- YAML 기반 데이터 추출 엔진(Bronze)
- 단위/표기 통일 로직(Silver)
- 최종 스키마 생성 + QC 검증(Gold)
- 격리/로그/테스트/문서화

### 결과

- 연도마다 형식이 다른 원본을 → **공통 스키마 4개 테이블**로 통합
- 문제 있는 데이터는 → **자동으로 격리 폴더에 분리**
- 언제 다시 돌려도 → **같은 결과가 나오게 재현성 확보**

### 검증 방법

- **Level 1 (키 없이)**: CLI로 파이프라인 재실행 + pytest로 스키마 검증
- **Level 2 (OpenAI 키 있으면)**: 웹 앱에서 Activity Cliff 분석 + LLM 가설 생성까지 확인

### 링크

- [YouTube 데모](https://www.youtube.com/watch?v=q2CVhYdBhU4) (2:28~8:50이 제 파트)
- [현재 레포](https://github.com/chi2hoon/chemical-cliff-climbing)
- 핵심 파일: [`pipeline/cli.py`](pipeline/cli.py), [`pipeline/silver.py`](pipeline/silver.py), [`pipeline/gold.py`](pipeline/gold.py)

---

## 프로젝트 배경

### 한 줄 소개

특허 PDF에서 추출한 엉망인 Excel들을 정제·표준화해서,  
연구원이 **구조–활성 관계(어떤 구조 변화가 효과를 높이는지)**를 빠르게 탐색하고 가설을 만들 수 있게 돕는 파이프라인입니다.

### 목적

신약 개발은 보통 5단계를 거칩니다.

1. 타깃 발굴
2. 초기 후보 물질 찾기
3. **리드 최적화** ← 이 프로젝트 범위
4. 전임상(동물) 시험
5. 임상 시험

이 프로젝트는 **3단계(리드 최적화)** 에 초점을 맞췄습니다.  
"어떤 구조 변경이 효과를 높이고 독성을 낮출까?"라는 가설을 빨리 만들고 검증 후보를 줄이는 게 목표였습니다.

다만 이 레포에서는 그중에서도 **데이터 엔지니어링 관점**(난해한 원본 데이터 정리, QC, 재현성)에 집중해 정리해 두었습니다.

---

## 데이터가 왜 문제였나

원본은 Tabula 같은 도구로 특허 PDF를 Excel로 추출한 파일들인데, 그대로는 쓸 수가 없었습니다.

### 1) 단위/표기가 뒤죽박죽

- 같은 "효과 측정값"인데 표기가 제각각: `µM`, `nM`, `%`, `++++` 등
- 예: `>= 2 µM`, `++++`, `45%` 는 그대로는 하나의 숫자 축으로 비교 불가

### 2) 시트 구조가 난해

- 한 워크북에 테이블 여러 개가 좌우로 붙어 있음
- 헤더가 없거나 오탈자 투성이
- 연도별로 시트 구조가 달라서 일관된 파서가 필요했음

### 3) 숫자와 텍스트가 섞여 있음

- 같은 표 안에 측정값 `1.23`, 실험 메모 `"0/10"`, `"전 개체 사망"` 같은 정성 데이터가 혼재
- 모델은 **숫자만** 써야 하니까, 나머지는 별도 테이블(context)로 분리해야 했음

### 4) 변환 오류/QC 흔적

- PDF → Excel 변환 과정에서 이미지가 글자로 잘못 인식됨
- 수동 검수 히스토리(색상 표시, 주석)가 그대로 남아 있음
- 이걸 그대로 쓰면 나중에 분석에서 혼선 발생

---

## 제가 한 일

이 레포는 팀 프로젝트를 포크한 것이고, 저는 **데이터 파이프라인 설계·구현**을 맡았습니다.

### 핵심 파일들

| 파일 | 역할 |
|------|------|
| `pipeline/cli.py` | Bronze/Silver/Gold 명령을 제공하는 CLI 진입점 |
| `pipeline/bronze_dispatch.py` | 연도별 어댑터 레지스트리 — 새 연도 추가 시 여기에 한 줄만 쓰면 됨 |
| `pipeline/adapters/_yaml_ingest.py` | YAML 파일로 시트 구조를 정의하고, 공통 엔진으로 추출하는 로직 |
| `schemas/silver/*.yaml` | 연도별 시트 구조/예외를 코드가 아닌 설정으로 분리 |
| `pipeline/silver.py` | 단위 통일, 기호 파싱, 숫자/텍스트 분리 |
| `pipeline/gold.py` | 최종 스키마 생성 + QC 검증 + 이상값 격리 |
| `pipeline/validate/hooks.py` | 안정 정렬, 격리 CSV 기록, 실행 로그 저장 |
| `modules/io_utils.py` | Gold 테이블을 읽어서 앱이 쓸 수 있는 형태로 변환 |
| `modules/context_builder.py` | LLM 가설 생성에 필요한 메타 정보 JSON 구성 |

### 간단히 요약하면

- **Bronze**: YAML 설정으로 원본 Excel을 슬라이스하고, QC 플래그/출처 정보를 붙임
- **Silver**: 단위 통일하고, 기호/등급 파싱하고, 숫자 vs 맥락(텍스트)를 분리
- **Gold**: 고정 스키마 4개 테이블로 수렴 + 이상값은 격리 폴더로

---

## 파이프라인 구조

```mermaid
flowchart LR
    Raw["원본 Excel<br/>(2017_raw.xlsx 등)"] -->|YAML 설정<br/>schemas/silver/{year}.yaml| Bronze["Bronze<br/>추출 엔진"]
    Bronze -->|tables/*.csv<br/>*_long.csv| Silver["Silver<br/>정제/표준화"]
    Silver -->|compounds_silver<br/>assay_readings_silver<br/>assay_context_silver| Gold["Gold<br/>최종 스키마"]

    Bronze -->|실행 로그| Logs["Logs<br/>logs/manifest/{year}.json"]
    Silver -->|로그 업데이트| Logs
    Bronze -->|문제 데이터| QBronze["Quarantine<br/>data/quarantine/bronze/*"]
    Gold -->|규칙 위반| QGold["Quarantine<br/>data/quarantine/gold/*"]
```

자세한 설명: [`docs/architecture.md`](docs/architecture.md), [`docs/sequence.md`](docs/sequence.md)

---

## 결과물

### Gold 테이블 4개

연도별로 아래 4개 CSV가 만들어집니다.

```
data/gold/{year}/
├── assay_readings.csv    # 효과 측정값 (어떤 화합물이 어떤 타깃에서 얼마나 효과가 있는지)
├── compounds.csv         # 화합물 마스터 (구조 정보)
├── compound_props.csv    # 화합물 속성 (분자량, NMR, QC 플래그 등)
└── assay_context.csv     # 맥락 정보 (실험 메모, 독성 데이터 등)
```

스키마 상세: [`docs/schema.md`](docs/schema.md)

### 정제 효과

**Before (원본)**
- 단위: `µM`, `nM`, `%`, `++++` 제각각
- 표기: `>= 2 µM`, `0/10`, `전 개체 사망` 섞여 있음
- 시트: 테이블 여러 개 붙어 있고 헤더 누락

**After (Gold)**
- 단위: `uM` 하나로 통일 (`value_std`, `unit_std`)
- 기호: `>`, `=` 따로 분리 (`qualifier`)
- 숫자/텍스트: 측정값은 `assay_readings`, 메모는 `assay_context`로 분리
- 재현성: 언제 다시 돌려도 같은 결과

---

## 직접 돌려보기

### 환경 세팅 (macOS 기준)

```bash
cd /path/to/this/repo
python3 -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
```

### Windows는 이렇게

```powershell
cd C:\path\to\chemical-cliff-climbing
python -m venv .venv
.\.venv\Scripts\Activate.ps1
pip install -r requirements.txt
```

### Level 1: 키 없이 파이프라인 검증

```bash
# 2017년 데이터 기준으로 Bronze → Silver → Gold 재생성
python -m pipeline.cli bronze --year 2017
python -m pipeline.cli silver --year 2017
python -m pipeline.cli gold --years 2017

# (선택) 스키마 테스트
pytest pipeline/tests/test_silver_outputs.py pipeline/tests/test_gold_outputs.py
```

실행하면:
- `data/bronze/2017`, `data/silver/2017`, `data/gold/2017` 폴더가 다시 만들어집니다
- pytest는 필수 컬럼이 다 있는지 확인합니다

### Level 2: LLM 가설 생성까지 (OpenAI 키 필요)

```bash
export OPENAI_API_KEY="sk-..."   # 본인이 발급받은 키
streamlit run app.py
```

브라우저에서 앱이 열리면:
- **탭 1~2 (키 없어도 됨)**: Gold 데이터 불러와서 Activity Cliff 히트맵으로 탐색
- **탭 3~6 (키 필요)**: 선택한 화합물 쌍에 대해 LLM이 구조–활성 가설 생성·평가·수정

---

## 한계와 개선 아이디어

### 지금 한계

- 4개 연도 데이터에만 맞춰져 있어서, 새 특허가 오면 YAML 추가 작업 필요
- 단위/등급 파싱 규칙이 고정되어 있어서, 완전히 새로운 표기는 수동 보정 필요
- Gold 검증은 최소 스키마만 체크 (통계적 아웃라이어 탐지는 없음)
- LLM 부분은 OpenAI API에 의존 (과금/Rate limit 영향)

### 다음에 하고 싶은 것

- YAML 스키마 템플릿/체커 만들어서 설정 오류 사전 방지
- Silver/Gold에 분포 기반 아웃라이어 탐지 붙이기
- LLM 부분을 벤더 중립으로 추상화 (오픈소스 LLM 대응)

---

## SQL로 쓰려면

Gold 테이블은 그대로 웨어하우스(BigQuery, Snowflake 등)에 올릴 수 있게 설계했습니다.

### 테이블 관계

- `compounds(compound_key)` ← 화합물 마스터 (PK)
- `assay_readings(..., compound_key, ...)` ← 측정값 (FK)
- `compound_props(..., compound_key, ...)` ← 속성 (FK)
- `assay_context` ← 실험 메모/독성 정보 (느슨한 조인)

### 예시 쿼리

```sql
-- 타깃별 상위 활성 화합물
SELECT
  ar.target_id,
  ar.compound_key,
  c.smiles_canonical AS structure,
  MIN(ar.value_std) AS best_value
FROM assay_readings ar
JOIN compounds c ON ar.compound_key = c.compound_key
WHERE ar.assay_id = 'cell.cytotoxicity.ic50'
  AND ar.unit_std = 'uM'
GROUP BY ar.target_id, ar.compound_key, c.smiles_canonical
ORDER BY ar.target_id, best_value ASC
LIMIT 100;
```

더 많은 예시: [`docs/schema.md`](docs/schema.md)

---

## 보안

- 이 레포에 **API 키, 토큰, 내부 URL은 없습니다**
- `.gitignore`에서 `.env`, `openAI_key.txt`, 원본 데이터 등 제외됨
- 키는 환경변수나 로컬 파일로만 읽고, 커밋하지 않습니다

---

## 팀 프로젝트 출처

이 레포는 [원본 팀 레포](https://github.com/nazirite96/chemical-cliff-climbing)를 포크한 것입니다.

- 전체 아이디어, LLM 프롬프트, 앱 UI는 **팀 공동 작업**
- 저는 `/pipeline` **데이터 파이프라인과 Gold 스키마 설계·구현**을 담당
- 이 포트폴리오 브랜치에서는:
  - `/pipeline` 정리
  - Medallion 아키텍처 문서화 ([`docs/`](docs/))
  - README를 **데이터 파이프라인/재현성/QC** 중심으로 개편

---

## 상세 문서

- [`docs/architecture.md`](docs/architecture.md) — Bronze/Silver/Gold 각 단계 역할
- [`docs/sequence.md`](docs/sequence.md) — CLI 실행 플로우 시퀀스 다이어그램
- [`docs/schema.md`](docs/schema.md) — Gold 테이블 스키마/관계/SQL 예시

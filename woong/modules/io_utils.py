import io
import re
from typing import Dict, Optional, Tuple

import pandas as pd
try:
    from rdkit import Chem
except Exception:  # RDKit가 없거나 로드 실패 시 이름 기반 휴리스틱만 사용
    Chem = None


def _normalize_column_name(name: str) -> str:
    if name is None:
        return ""
    # 소문자, 앞뒤 공백 제거, 영숫자만 남김
    lowered = str(name).strip().lower()
    return re.sub(r"[^a-z0-9]", "", lowered)


def _smiles_parse_ratio(series: pd.Series, sample_size: int = 300) -> float:
    if Chem is None:
        return 0.0
    values = series.dropna().astype(str).str.strip()
    if values.empty:
        return 0.0
    sample = values.sample(min(sample_size, len(values)), random_state=42)
    total = 0
    ok = 0
    for v in sample:
        if not v:
            continue
        total += 1
        mol = Chem.MolFromSmiles(v)
        if mol is not None:
            ok += 1
    return (ok / total) if total > 0 else 0.0


def _numeric_ratio(series: pd.Series, sample_size: int = 500) -> float:
    values = series.dropna().astype(str).str.strip()
    if values.empty:
        return 0.0
    sample = values.sample(min(sample_size, len(values)), random_state=42)
    nums = pd.to_numeric(sample, errors="coerce")
    return float(nums.notna().mean())


def _detect_header_row(df_without_header: pd.DataFrame) -> Optional[int]:
    # 어느 행에 실제 헤더(예: SMILES)가 있는지 탐지
    for row_idx in range(len(df_without_header)):
        row = df_without_header.iloc[row_idx].astype(str).fillna("").tolist()
        normalized_cells = [_normalize_column_name(cell) for cell in row]
        if any(("smile" in cell) for cell in normalized_cells):
            return row_idx
    return None


def _choose_smiles_and_activity_columns(df: pd.DataFrame) -> Tuple[str, Optional[str]]:
    normalized_map = {col: _normalize_column_name(col) for col in df.columns}

    # 1) SMILES 후보: 이름에 'smile' 포함, 'structure', 'mol' 등
    name_based_candidates = [
        col for col, norm in normalized_map.items()
        if ("smile" in norm) or ("structure" in norm) or (norm in {"mol", "molecule"})
    ]

    # 2) RDKit로 실제 SMILES 파싱률 기반 스코어링 (가능할 때)
    scored: Dict[str, float] = {}
    candidates = name_based_candidates if name_based_candidates else list(df.columns)
    for col in candidates:
        try:
            ratio = _smiles_parse_ratio(df[col])
        except Exception:
            ratio = 0.0
        scored[col] = ratio

    # 최적 컬럼 선택: 파싱률 0.2 이상 중 최고, 없으면 이름 기반 첫 번째, 그것도 없으면 첫 컬럼
    smiles_col = None
    if scored:
        best_col = max(scored, key=lambda c: scored[c])
        if scored[best_col] >= 0.2:
            smiles_col = best_col
    if smiles_col is None:
        smiles_col = name_based_candidates[0] if name_based_candidates else df.columns[0]

    # 3) Activity 후보: 이름 키워드 가중치 + 숫자 비율
    activity_keywords = (
        "pic50", "ic50", "ec50", "ki", "kd", "potency", "activity", "response", "value"
    )
    act_scores: Dict[str, float] = {}
    for col, norm in normalized_map.items():
        if col == smiles_col:
            continue
        keyword_bonus = 1.0 if any(k in norm for k in activity_keywords) else 0.0
        try:
            num_ratio = _numeric_ratio(df[col])
        except Exception:
            num_ratio = 0.0
        # 기본 점수: 숫자 비율 + 키워드 보너스
        act_scores[col] = num_ratio + keyword_bonus

    activity_col = None
    if act_scores:
        activity_col = max(act_scores, key=lambda c: act_scores[c])

    # 폴백: 그래도 없으면 두 번째(또는 세 번째) 컬럼
    if activity_col is None and len(df.columns) >= 2:
        second = df.columns[1]
        activity_col = second if second != smiles_col else (df.columns[2] if len(df.columns) >= 3 else None)

    return smiles_col, activity_col


def _convert_activity_to_numeric(series: pd.Series) -> pd.Series:
    # 1) 숫자 변환 시도
    numeric = pd.to_numeric(series, errors="coerce")
    if numeric.notna().sum() >= max(1, int(0.5 * series.notna().sum())):
        return numeric

    # 2) + 기호 등급 변환: + -> 1, ++ -> 2, ...
    def plus_to_num(x):
        if pd.isna(x):
            return pd.NA
        s = str(x).strip()
        # +만 추출
        plus_only = re.sub(r"[^+]", "", s)
        if plus_only:
            return len(plus_only)
        # 3) 문자열 내 숫자 추출 (예: "45 nM")
        m = re.search(r"-?\d+\.?\d*", s)
        if m:
            try:
                return float(m.group())
            except Exception:
                return pd.NA
        return pd.NA

    converted = series.apply(plus_to_num)
    return pd.to_numeric(converted, errors="coerce")


def load_smiles_activity_csv(file_like) -> Tuple[pd.DataFrame, Dict[str, Optional[str]]]:
    """
    다양한 변형 CSV(헤더 행 위치 불일치, 등급형 활동값 등)를 유연하게 읽고
    SMILES/활성도 컬럼을 자동 추정합니다.

    Returns:
        df: 정리된 DataFrame (활성도 컬럼은 가능한 경우 숫자형으로 변환)
        suggestion: {"smiles_col": str, "activity_col": Optional[str]}
    """
    # 업로드 파일은 스트림이므로 여러 번 읽기 위해 바이트로 고정
    data_bytes = file_like.getvalue() if hasattr(file_like, "getvalue") else file_like.read()

    # 헤더 없는 상태로 일단 로드해서 헤더 행 탐지
    tmp = pd.read_csv(io.BytesIO(data_bytes), header=None, dtype=str, engine="python")
    header_row = _detect_header_row(tmp)

    if header_row is not None:
        df = pd.read_csv(io.BytesIO(data_bytes), header=header_row, dtype=str, engine="python")
    else:
        # 폴백: 기본 헤더 가정
        df = pd.read_csv(io.BytesIO(data_bytes), header=0, dtype=str, engine="python")

    # 공백/NA 정리
    df = df.dropna(how="all")
    df.columns = [str(c).strip() for c in df.columns]

    # 컬럼 자동 선택
    smiles_col, activity_col = _choose_smiles_and_activity_columns(df)

    # SMILES 비어있는 행 제거
    df = df[df[smiles_col].notna() & (df[smiles_col].astype(str).str.strip() != "")]

    # 활동값 숫자화 시도 (활동 컬럼이 존재하는 경우)
    if activity_col is not None and activity_col in df.columns:
        df[activity_col] = _convert_activity_to_numeric(df[activity_col])
        # 숫자 결측이 너무 많으면 그대로 두되, 연산 전 필터링할 수 있도록 둔다

    # 인덱스 리셋
    df = df.reset_index(drop=True)

    return df, {"smiles_col": smiles_col, "activity_col": activity_col}


def save_hypothesis_to_md(hypothesis: str, filename: str) -> None:
    """
    주어진 가설 텍스트를 마크다운 파일로 저장합니다.

    Args:
        hypothesis (str): 저장할 가설 텍스트.
        filename (str): 저장할 파일 이름 (e.g., "hypothesis_1.md").
    """
    try:
        with open(filename, 'w', encoding='utf-8') as f:
            f.write(hypothesis)
    except IOError as e:
        # 실제 앱에서는 st.error() 등을 사용하여 사용자에게 알릴 수 있습니다.
        print(f"Error saving file {filename}: {e}")

def parse_hypothesis_md(file_content: str) -> Dict:
    """
    저장된 가설 마크다운 파일의 내용을 파싱하여 초기 데이터와 가설 본문을 분리합니다.

    Args:
        file_content (str): 마크다운 파일의 전체 내용.

    Returns:
        Dict: 파싱된 데이터를 담은 딕셔너리. 실패 시 일부 필드가 None일 수 있음.
              예: {
                  'smiles1': str, 'activity1': float, 
                  'smiles2': str, 'activity2': float, 
                  'hypothesis_body': str
              }
    """
    data = {
        'smiles1': None, 'activity1': None,
        'smiles2': None, 'activity2': None,
        'hypothesis_body': None
    }

    # 정규식을 사용하여 화합물 정보 추출
    # 화합물 1 (저활성)
    match1 = re.search(r"- \*\*화합물 1 \(상대적 저활성\):\*\* `([^`]+)` \(활성도: ([\d\.]+)\)", file_content)
    if match1:
        data['smiles1'] = match1.group(1)
        data['activity1'] = float(match1.group(2))

    # 화합물 2 (고활성)
    match2 = re.search(r"- \*\*화합물 2 \(상대적 고활성\):\*\* `([^`]+)` \(활성도: ([\d\.]+)\)", file_content)
    if match2:
        data['smiles2'] = match2.group(1)
        data['activity2'] = float(match2.group(2))

    # 헤더와 본문 분리
    parts = re.split(r"\n---\n", file_content, maxsplit=1)
    if len(parts) > 1:
        data['hypothesis_body'] = parts[1].strip()
    else:
        # 헤더가 없는 경우 전체를 본문으로 간주
        data['hypothesis_body'] = file_content.strip()

    return data

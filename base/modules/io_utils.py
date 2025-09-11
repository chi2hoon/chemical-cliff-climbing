import io
import re
from typing import Dict, Optional, Tuple, List
import yaml

import pandas as pd
import os


def _normalize_column_name(name: str) -> str:
    if name is None:
        return ""
    # 소문자, 앞뒤 공백 제거, 영숫자만 남김
    lowered = str(name).strip().lower()
    return re.sub(r"[^a-z0-9]", "", lowered)


def _detect_header_row(df_without_header: pd.DataFrame) -> Optional[int]:
    # 어느 행에 실제 헤더(예: SMILES)가 있는지 탐지
    for row_idx in range(len(df_without_header)):
        row = df_without_header.iloc[row_idx].astype(str).fillna("").tolist()
        normalized_cells = [_normalize_column_name(cell) for cell in row]
        if any(cell == "smiles" for cell in normalized_cells):
            return row_idx
    return None


def _choose_smiles_and_activity_columns(df: pd.DataFrame) -> Tuple[str, Optional[str]]:
    normalized_map = {col: _normalize_column_name(col) for col in df.columns}

    # SMILES 컬럼 선택: normalize == "smiles" 우선
    smiles_candidates = [col for col, norm in normalized_map.items() if norm == "smiles"]
    smiles_col = smiles_candidates[0] if smiles_candidates else df.columns[0]

    # Activity 컬럼 선택 우선순위: ic50 포함 > activity 포함 > 두 번째 컬럼
    activity_candidates = [
        col for col, norm in normalized_map.items() if ("ic50" in norm or "activity" in norm)
    ]
    activity_col = None
    for col in activity_candidates:
        if col != smiles_col:
            activity_col = col
            break
    if activity_col is None and len(df.columns) >= 2:
        # 안전한 폴백: 두 번째 컬럼(단, SMILES와 다를 것)
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


# (레거시) hoon 관련 헬퍼 제거됨

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


def load_gold_data(year: str = "2017", data_root: str = "base/data", panel_id: Optional[str] = None, cell_line: Optional[str] = None, target_id: Optional[str] = None) -> pd.DataFrame:
    """base/data/gold/{year}에서 assay_readings/compounds/compound_props를 조인해
    SMILES/Activity 스키마를 반환한다.
    """
    gold_dir = os.path.join(data_root, "gold", str(year))
    assay_path = os.path.join(gold_dir, "assay_readings.csv")
    comps_path = os.path.join(gold_dir, "compounds.csv")
    props_path = os.path.join(gold_dir, "compound_props.csv")
    if not os.path.exists(assay_path):
        return pd.DataFrame(columns=["SMILES", "Activity"])
    try:
        df_assay = pd.read_csv(assay_path, dtype=str, keep_default_na=False, na_filter=False)
        df_comps = pd.read_csv(comps_path, dtype=str, keep_default_na=False, na_filter=False) if os.path.exists(comps_path) else pd.DataFrame()
        df_props = pd.read_csv(props_path, dtype=str, keep_default_na=False, na_filter=False) if os.path.exists(props_path) else pd.DataFrame()
        def to_float(x):
            try:
                return float(x)
            except Exception:
                return None
        df_assay["Activity"] = df_assay.get("value_std", "").apply(to_float)
        if panel_id:
            if "panel_id" in df_assay.columns:
                df_assay = df_assay[df_assay["panel_id"] == panel_id]
            elif "target_id" in df_assay.columns:
                def _match_panel(t):
                    s = str(t)
                    if not s.startswith("cell:") or "." not in s:
                        return False
                    body = s.split(":", 1)[1]
                    return body.split(".", 1)[0] == panel_id
                df_assay = df_assay[df_assay["target_id"].map(_match_panel)]
        if target_id:
            if "target_id" in df_assay.columns:
                df_assay = df_assay[df_assay["target_id"] == target_id]
        if cell_line:
            if "cell_line" in df_assay.columns:
                df_assay = df_assay[df_assay["cell_line"] == cell_line]
            elif "target_id" in df_assay.columns:
                def _match_cell(t):
                    s = str(t)
                    if not s.startswith("cell:") or "." not in s:
                        return False
                    return s.split(":", 1)[1].split(".", 1)[1] == cell_line
                df_assay = df_assay[df_assay["target_id"].map(_match_cell)]
        df = df_assay.copy()
        # 1) compound_key 보강: props에서 compound_id→compound_key 주입(빈 키만)
        if ("compound_id" in df.columns) and (len(df_props) > 0) and ("compound_id" in df_props.columns) and ("compound_key" in df_props.columns):
            props_map = df_props[["compound_id","compound_key"]].drop_duplicates()
            df = df.merge(props_map, on="compound_id", how="left", suffixes=("", "_props"))
            # 우선순위: 기존 키 > props 키
            if 'compound_key_props' in df.columns:
                df['compound_key'] = df['compound_key'].where(df['compound_key'].astype(str).str.strip()!="", df['compound_key_props'])
                df = df.drop(columns=['compound_key_props'])
        # 2) compounds.csv로부터 canonical SMILES 조인
        if ("compound_key" in df.columns) and (len(df_comps) > 0) and ("compound_key" in df_comps.columns) and ("smiles_canonical" in df_comps.columns):
            df = df.merge(df_comps[["compound_key","smiles_canonical"]].drop_duplicates(), on="compound_key", how="left")
        # 3) silver compounds_silver의 smiles_raw로 폴백(표시용)
        comp_silver_path = os.path.join(data_root, "silver", str(year), "compounds_silver.csv")
        comp_sil = pd.read_csv(comp_silver_path, dtype=str) if os.path.exists(comp_silver_path) else pd.DataFrame()
        if (len(comp_sil)>0) and ("compound_id" in df.columns) and ("compound_id" in comp_sil.columns) and ("smiles_raw" in comp_sil.columns):
            df = df.merge(comp_sil[["compound_id","smiles_raw"]].drop_duplicates(), on="compound_id", how="left")
        # 4) 최종 SMILES 선택: canonical 우선, 없으면 raw
        def _pick(a,b):
            a = str(a).strip() if a is not None else ""
            b = str(b).strip() if b is not None else ""
            return a if a else b
        df["SMILES"] = [ _pick(a,b) for a,b in zip(df.get("smiles_canonical",""), df.get("smiles_raw","")) ]
        cols = ["SMILES", "Activity"]
        for meta_col in ["compound_id", "assay_id", "target_id", "unit_std", "qualifier"]:
            if meta_col in df.columns:
                cols.append(meta_col)
        out = df[cols].copy()
        # 유효 행만 유지
        out = out.dropna(subset=["SMILES", "Activity"]).copy()
        smi_str = out["SMILES"].astype(str).str.strip()
        out = out[(smi_str != "") & (smi_str.str.lower() != "nan")]
        return out.reset_index(drop=True)
    except Exception as e:
        print(f"Error loading gold data: {e}")
        return pd.DataFrame(columns=["SMILES", "Activity"])


def get_available_gold_years(data_root: str = "base/data") -> List[str]:
    """Get list of available years in gold data directory."""
    gold_root = os.path.join(data_root, "gold")
    if not os.path.exists(gold_root):
        return []
    
    years = []
    for item in os.listdir(gold_root):
        if os.path.isdir(os.path.join(gold_root, item)):
            gold_file = os.path.join(gold_root, item, "assay_readings.csv")
            if os.path.exists(gold_file):
                years.append(item)
    
    return sorted(years)


def get_available_panel_ids(year: str = "2017", data_root: str = "base/data") -> List[str]:
    """Get list of available panel IDs for a given year."""
    gold_path = os.path.join(data_root, "gold", str(year), "assay_readings.csv")
    
    if not os.path.exists(gold_path):
        return []
    
    try:
        df = pd.read_csv(gold_path, dtype=str, keep_default_na=False, na_filter=False)
        if "panel_id" in df.columns:
            return sorted([p for p in df["panel_id"].unique().tolist() if str(p).strip() != ""])
        if "target_id" in df.columns:
            panels = []
            for t in df["target_id"].astype(str).tolist():
                if t.startswith("cell:") and "." in t:
                    body = t.split(":", 1)[1]
                    panel = body.split(".", 1)[0]
                    panels.append(panel)
            return sorted(sorted(set([p for p in panels if p])))
        return []
    except Exception:
        return []


def get_all_available_panels_and_years(data_root: str = "base/data") -> Dict[str, List[str]]:
    """Get all available panels across all years.
    
    Returns:
        Dict[str, List[str]]: Dictionary mapping panel_id to list of years where it's available
    """
    panel_years = {}
    
    # Get all available years
    years = get_available_gold_years(data_root)
    
    # For each year, get its panels
    for year in years:
        panels = get_available_panel_ids(year, data_root)
        for panel in panels:
            if panel not in panel_years:
                panel_years[panel] = []
            panel_years[panel].append(year)
    
    return panel_years


def get_cell_lines_for_panel(year: str, panel_id: str, data_root: str = "base/data") -> List[str]:
    """assay_readings.csv의 target_id에서 해당 패널의 cell_line 목록을 파생한다."""
    gold_path = os.path.join(data_root, "gold", str(year), "assay_readings.csv")
    if not os.path.exists(gold_path):
        return []
    try:
        df = pd.read_csv(gold_path, dtype=str, keep_default_na=False, na_filter=False)
        clines = []
        if "cell_line" in df.columns and "panel_id" in df.columns:
            sub = df[df["panel_id"] == panel_id]
            if len(sub) == 0:
                return []
            return sorted([c for c in sub["cell_line"].unique().tolist() if str(c).strip() != ""])
        for t in df.get("target_id", []):
            s = str(t)
            if not s.startswith("cell:") or "." not in s:
                continue
            body = s.split(":", 1)[1]
            p, rest = body.split(".", 1)
            if p != panel_id:
                continue
            clines.append(rest)
        return sorted(sorted(set([c for c in clines if c])))
    except Exception:
        return []


# ---- Panel/Target helpers for UI (strict behaviors) ----

def has_panel_column(year: str = "2017", data_root: str = "base/data") -> bool:
    """해당 연도의 gold/assay_readings.csv에 panel_id 컬럼이 존재하는지 반환한다."""
    gold_path = os.path.join(data_root, "gold", str(year), "assay_readings.csv")
    if not os.path.exists(gold_path):
        return False
    try:
        df = pd.read_csv(gold_path, nrows=1, dtype=str, keep_default_na=False, na_filter=False)
        return "panel_id" in df.columns
    except Exception:
        return False


def get_available_panels_strict(year: str = "2017", data_root: str = "base/data") -> List[str]:
    """panel_id 컬럼이 있는 경우에만 panel_id 목록을 반환한다. 없으면 빈 리스트."""
    gold_path = os.path.join(data_root, "gold", str(year), "assay_readings.csv")
    if not os.path.exists(gold_path):
        return []
    try:
        df = pd.read_csv(gold_path, dtype=str, keep_default_na=False, na_filter=False)
        if "panel_id" not in df.columns:
            return []
        vals = [str(x).strip() for x in df["panel_id"].dropna().tolist()]
        return sorted([v for v in set(vals) if v])
    except Exception:
        return []


def get_available_targets(year: str = "2017", data_root: str = "base/data") -> List[str]:
    """해당 연도의 target_id 고유값 목록을 반환한다. 없으면 빈 리스트."""
    gold_path = os.path.join(data_root, "gold", str(year), "assay_readings.csv")
    if not os.path.exists(gold_path):
        return []
    try:
        df = pd.read_csv(gold_path, dtype=str, keep_default_na=False, na_filter=False)
        if "target_id" not in df.columns:
            return []
        vals = [str(x).strip() for x in df["target_id"].dropna().tolist()]
        return sorted([v for v in set(vals) if v])
    except Exception:
        return []

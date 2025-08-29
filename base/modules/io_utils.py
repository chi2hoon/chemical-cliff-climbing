import io
import re
from typing import Dict, Optional, Tuple, List

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


# ==== Hoon pipeline integration helpers ====

def _hoon_silver_dir(data_root: str, dataset_id: str) -> str:
    return os.path.join(data_root, "silver", str(dataset_id))


def list_hoon_datasets(data_root: str = "hoon/data") -> List[str]:
    """Return dataset IDs that have silver outputs under hoon/data/silver/{id}."""
    silver_root = os.path.join(data_root, "silver")
    try:
        return sorted([d for d in os.listdir(silver_root) if os.path.isdir(os.path.join(silver_root, d)) and d != "all"])
    except Exception:
        return []


def list_hoon_groups(data_root: str = "hoon/data", dataset_id: Optional[str] = None) -> pd.DataFrame:
    """
    List available groups (assay_id, cell_line) from measurements_std.csv.
    Returns a DataFrame with columns: dataset_id, assay_id, cell_line, n_rows.
    """
    rows = []
    ds_list = [dataset_id] if dataset_id else list_hoon_datasets(data_root)
    for ds in ds_list:
        silver_dir = _hoon_silver_dir(data_root, ds)
        path = os.path.join(silver_dir, "measurements_std.csv")
        if not os.path.exists(path):
            continue
        try:
            df = pd.read_csv(path, dtype=str, keep_default_na=False, na_filter=False)
            # use available columns robustly
            assay_col = "assay_id" if "assay_id" in df.columns else None
            cell_col = "cell_line" if "cell_line" in df.columns else None
            if assay_col is None and cell_col is None:
                continue
            gcols = [c for c in [assay_col, cell_col] if c]
            grouped = df.groupby(gcols).size().reset_index(name="n_rows")
            for _, r in grouped.iterrows():
                rows.append({
                    "dataset_id": ds,
                    "assay_id": r.get(assay_col, "") if assay_col else "",
                    "cell_line": r.get(cell_col, "") if cell_col else "",
                    "n_rows": int(r["n_rows"]) if "n_rows" in r else 0,
                })
        except Exception:
            continue
    return pd.DataFrame(rows, columns=["dataset_id", "assay_id", "cell_line", "n_rows"]).sort_values(["dataset_id", "assay_id", "cell_line"]).reset_index(drop=True)


def load_hoon_group_as_dataframe(
    dataset_id: str,
    assay_id: Optional[str] = None,
    cell_line: Optional[str] = None,
    data_root: str = "hoon/data",
) -> pd.DataFrame:
    """
    Load hoon silver measurements + canonical smiles for a group and return a simple DataFrame
    with columns: SMILES, Activity (float). Additional metadata columns may be included.
    """
    silver_dir = _hoon_silver_dir(data_root, dataset_id)
    meas_path = os.path.join(silver_dir, "measurements_std.csv")
    comp_path = os.path.join(silver_dir, "compounds_canonical.csv")
    if not (os.path.exists(meas_path) and os.path.exists(comp_path)):
        return pd.DataFrame(columns=["SMILES", "Activity"])  # graceful empty

    m = pd.read_csv(meas_path, dtype=str, keep_default_na=False, na_filter=False)
    c = pd.read_csv(comp_path, dtype=str, keep_default_na=False, na_filter=False)
    # robust float cast
    def to_float(x):
        try:
            return float(x)
        except Exception:
            return None
    # filter group
    if assay_id and "assay_id" in m.columns:
        m = m[m["assay_id"] == assay_id]
    if cell_line and "cell_line" in m.columns:
        m = m[m["cell_line"] == cell_line]
    # join canonical smiles
    if "compound_id" in m.columns and "compound_id" in c.columns:
        df = m.merge(c[["compound_id", "smiles_canonical"]], on="compound_id", how="left")
    else:
        df = m.copy()
        df["smiles_canonical"] = ""
    # build simple schema
    df["Activity"] = [to_float(x) for x in df.get("value_std", [])]
    df.rename(columns={"smiles_canonical": "SMILES"}, inplace=True)
    df = df.dropna(subset=["SMILES", "Activity"])  # ensure usable rows
    # keep a minimal set upfront; but retain metadata for reference
    cols_front = ["SMILES", "Activity"]
    other_cols = [c for c in df.columns if c not in cols_front]
    return df[cols_front + other_cols]


def load_hoon_ac_pairs(data_root: str = "hoon/data") -> pd.DataFrame:
    """Load precomputed activity cliff pairs from hoon/data/silver/all/ac_pairs.csv and map to base schema."""
    path = os.path.join(data_root, "silver", "all", "ac_pairs.csv")
    if not os.path.exists(path):
        return pd.DataFrame(columns=["SMILES_1", "Activity_1", "SMILES_2", "Activity_2", "Similarity", "Activity_Diff"])
    ac = pd.read_csv(path, dtype=str, keep_default_na=False, na_filter=False)
    # numeric cast helpers
    def tf(x):
        try:
            return float(x)
        except Exception:
            return None
    out = pd.DataFrame({
        "SMILES_1": ac.get("smiles_1"),
        "Activity_1": [tf(x) for x in ac.get("value_std_1", [])],
        "SMILES_2": ac.get("smiles_2"),
        "Activity_2": [tf(x) for x in ac.get("value_std_2", [])],
        "Similarity": [tf(x) for x in ac.get("similarity", [])],
        "Activity_Diff": [tf(x) for x in ac.get("delta_value_std", [])],
    })
    # drop invalid rows
    out = out.dropna(subset=["SMILES_1", "SMILES_2", "Activity_1", "Activity_2", "Similarity", "Activity_Diff"])  
    return out.reset_index(drop=True)

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


def load_hoon_gold_data(year: str = "2017", data_root: str = "hoon/data", panel_id: Optional[str] = None) -> pd.DataFrame:
    """
    Load gold dataset from hoon/data/gold/{year}/measurements_gold.csv
    and return a DataFrame ready for base app analysis.

    Args:
        year (str): Dataset year (default: "2017")
        data_root (str): Root directory for hoon data
        panel_id (Optional[str]): Panel ID to filter by (e.g., "blca", "brca", etc.)

    Returns:
        pd.DataFrame: Gold dataset with columns ready for SMILES and activity analysis
    """
    gold_dir = os.path.join(data_root, "gold", str(year))
    gold_path = os.path.join(gold_dir, "measurements_gold.csv")

    if not os.path.exists(gold_path):
        return pd.DataFrame(columns=["SMILES", "Activity"])

    try:
        df = pd.read_csv(gold_path, dtype=str, keep_default_na=False, na_filter=False)

        # Filter by panel_id if specified
        if panel_id and "panel_id" in df.columns:
            df = df[df["panel_id"] == panel_id]

        # Convert value_std to float
        def to_float(x):
            try:
                return float(x)
            except Exception:
                return None

        df["value_std_float"] = df.get("value_std", "").apply(to_float)

        # Map gold schema to base schema
        # smiles_canonical -> SMILES, value_std -> Activity
        df_mapped = pd.DataFrame({
            "SMILES": df.get("smiles_canonical", ""),
            "Activity": df["value_std_float"],
        })

        # Filter out rows with missing SMILES or Activity
        df_mapped = df_mapped.dropna(subset=["SMILES", "Activity"])
        df_mapped = df_mapped[df_mapped["SMILES"].astype(str).str.len() > 0]

        # Add metadata columns for reference
        metadata_cols = ["compound_id", "dataset_id", "assay_id", "panel_id", "cell_line",
                        "inchikey", "unit_std", "censor", "readout", "matrix",
                        "provenance_file", "provenance_sheet", "provenance_row",
                        "std_rule_id", "std_confidence"]

        for col in metadata_cols:
            if col in df.columns:
                df_mapped[col] = df[col]

        return df_mapped.reset_index(drop=True)

    except Exception as e:
        print(f"Error loading gold data: {e}")
        return pd.DataFrame(columns=["SMILES", "Activity"])


def get_available_gold_years(data_root: str = "hoon/data") -> List[str]:
    """Get list of available years in gold data directory."""
    gold_root = os.path.join(data_root, "gold")
    if not os.path.exists(gold_root):
        return []
    
    years = []
    for item in os.listdir(gold_root):
        if os.path.isdir(os.path.join(gold_root, item)):
            gold_file = os.path.join(gold_root, item, "measurements_gold.csv")
            if os.path.exists(gold_file):
                years.append(item)
    
    return sorted(years)


def get_available_panel_ids(year: str = "2017", data_root: str = "hoon/data") -> List[str]:
    """Get list of available panel IDs for a given year."""
    gold_path = os.path.join(data_root, "gold", str(year), "measurements_gold.csv")
    
    if not os.path.exists(gold_path):
        return []
    
    try:
        df = pd.read_csv(gold_path, dtype=str, keep_default_na=False, na_filter=False)
        if "panel_id" in df.columns:
            return sorted(df["panel_id"].unique().tolist())
        return []
    except Exception:
        return []

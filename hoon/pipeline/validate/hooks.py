import os
import json
import pandas as pd


def stable_sort(df, preferred_keys):
    """Args: df(DataFrame), preferred_keys(list[str]) -> DataFrame

    주어진 키 중 존재하는 컬럼으로 안정 정렬. 동률은 원래 순서 유지.
    특수 컬럼(compound_id, row_id, provenance_row)은 값이 모두 숫자형 문자열이면
    자연 정렬(숫자 기준)로 정렬한다.
    """
    present = [k for k in preferred_keys if k in df.columns]
    if not present:
        return df
    import re
    idx = df.reset_index().rename(columns={"index": "__ord__"})
    sort_cols = []
    tmp_cols = []
    numeric_candidates = {"compound_id", "row_id", "provenance_row"}
    for k in present:
        col = idx[k]
        use_col = k
        if k in numeric_candidates:
            try:
                # 모두 숫자형 문자열인지 검사(None/빈값 제외)
                ser = col.dropna().astype(str).str.strip()
                if len(ser) > 0 and ser.str.fullmatch(r"\d+").all():
                    tmp = f"__sort_{k}__"
                    idx[tmp] = ser.replace("", None).astype(float)
                    # NaN 보존: 원 컬럼 None일 때 NaN, 아닌 경우 숫자
                    # 정렬 키로 임시 컬럼 사용
                    use_col = tmp
                    tmp_cols.append(tmp)
            except Exception:
                pass
        sort_cols.append(use_col)
    out = idx.sort_values(by=sort_cols + ["__ord__"], kind="mergesort").drop(columns=["__ord__"] + tmp_cols, errors="ignore").reset_index(drop=True)
    return out


def bronze_checks(df, required_cols):
    """Args: df(DataFrame), required_cols(list[str]) -> (ok_df, fails_df)

    필수 컬럼 존재/결측 검사. 실패 레코드와 통과 레코드를 분리.
    """
    missing_cols = [c for c in required_cols if c not in df.columns]
    if missing_cols:
        # 컬럼 자체가 없으면 전체 실패로 간주
        fails = df.copy()
        fails["reason"] = "missing_columns:" + ",".join(missing_cols)
        return df.iloc[0:0].copy(), fails
    mask = None
    for c in required_cols:
        m = df[c].isna() | (df[c].astype(str).str.strip() == "")
        mask = m if mask is None else (mask | m)
    fails = df[mask].copy()
    if len(fails) > 0 and "reason" not in fails.columns:
        fails["reason"] = "missing_required_value"
    ok = df[~mask].copy()
    return ok, fails


def write_quarantine(stage, year, name, df):
    """Args: stage(str), year(str|int), name(str), df(DataFrame) -> str

    격리 CSV 경로로 저장하고 경로를 반환한다.
    """
    year = str(year)
    out_dir = os.path.join("data", "quarantine", stage)
    os.makedirs(out_dir, exist_ok=True)
    out = os.path.join(out_dir, f"{year}_{name}.csv")
    df.to_csv(out, index=False, encoding="utf-8")
    return out


def write_manifest(year, manifest):
    """Args: year(str|int), manifest(dict) -> str

    실행 매니페스트를 JSON으로 기록.
    """
    year = str(year)
    out_dir = os.path.join("logs", "manifest")
    os.makedirs(out_dir, exist_ok=True)
    out = os.path.join(out_dir, f"{year}.json")
    with open(out, "w", encoding="utf-8") as f:
        json.dump(manifest, f, ensure_ascii=False, indent=2)
    return out

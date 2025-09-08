import os
import json
import pandas as pd


def stable_sort(df, preferred_keys):
    """Args: df(DataFrame), preferred_keys(list[str]) -> DataFrame

    주어진 키 중 존재하는 컬럼으로 안정 정렬. 동률은 원래 순서 유지.
    특수 컬럼(compound_id, row_id, provenance_row)은 숫자 기준으로 정렬한다.
    - 가능한 경우 정수로 캐스팅하고, 실패 시 큰 숫자(말미)로 취급하여 안정적으로 정렬.
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
                # 숫자 정렬: 전처리 후 가능한 값을 숫자로 변환, 실패는 말미로 보냄
                ser_all = col.astype(str).str.strip()
                # 순수 숫자(또는 1.0 형태)인 경우: 직접 변환
                is_pure_num = ser_all.str.match(r"^\d+(\.\d+)?$")
                # 그 외에는 숫자 부분 추출(예: ' 1 ' → 1, '01' → 1)
                extracted = ser_all.where(is_pure_num, other=ser_all.str.extract(r"(\d+)", expand=False))
                tmp = f"__sort_{k}__"
                num = pd.to_numeric(extracted, errors="coerce")
                idx[tmp] = num.fillna(1e18)
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

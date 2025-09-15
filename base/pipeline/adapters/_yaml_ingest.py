import os
import yaml
import pandas as pd

from pipeline.parsers.engine import read_excel, extract_block, apply_header, melt_table, sanitize_strings, file_sha256, read_cell_colors
from pipeline.validate.hooks import stable_sort, bronze_checks, write_quarantine, write_manifest


def _repo_root():
    """Args: None -> str

    레포 루트(최상위)를 반환.
    """
    here = os.path.dirname(os.path.abspath(__file__))
    # .../base/pipeline/adapters → 세 단계 상위가 repo 루트
    return os.path.abspath(os.path.join(here, os.pardir, os.pardir, os.pardir))


def _apply_rename(df, mapping):
    if not mapping:
        return df
    cols = {str(k): str(v) for k, v in mapping.items()}
    return df.rename(columns=cols)


def run_yaml_bronze_ingest(year, yaml_path, out_dir):
    """Args: year(str|int), yaml_path(str), out_dir(str|None) -> dict

    YAML 스키마 우선으로 브론즈(ingest) 산출을 생성한다.
    산출물: tables/{id}.csv, *_long.csv, manifest.json
    """
    root = _repo_root()
    year = str(year)
    out_dir = out_dir or os.path.join("base", "data", "bronze", year)
    os.makedirs(out_dir, exist_ok=True)

    with open(yaml_path, encoding="utf-8") as f:
        cfg = yaml.safe_load(f)

    xls_path = cfg.get("file")
    if not xls_path:
        raise ValueError("YAML에 file 경로가 없습니다.")
    if not os.path.isabs(xls_path):
<<<<<<< HEAD:base/pipeline/adapters/_yaml_ingest.py
        # base 상대 또는 레포 루트 상대(base/...) 모두 지원
=======
        # hoon/ 상대 또는 레포 루트 상대(base/...) 모두 지원
>>>>>>> origin/main:hoon/pipeline/adapters/_yaml_ingest.py
        cand1 = os.path.join(root, xls_path)
        cand2 = os.path.join(os.path.abspath(os.path.join(root, os.pardir)), xls_path)
        xls_path = cand1 if os.path.exists(cand1) else cand2
    if not os.path.exists(xls_path):
        # 호환: base/data → hoon/data 폴백
        if "/base/data/" in ("/" + xls_path):
            alt = xls_path.replace("/base/data/", "/hoon/data/")
            if os.path.exists(alt):
                xls_path = alt
        if not os.path.exists(xls_path):
            raise FileNotFoundError(xls_path)

    yaml_abs = yaml_path if os.path.isabs(yaml_path) else os.path.join(root, yaml_path)

    manifest = {
        "input_file": os.path.relpath(xls_path, root),
        "input_sha256": file_sha256(xls_path),
        "yaml_path": os.path.relpath(yaml_abs, root),
        "yaml_sha256": file_sha256(yaml_abs) if os.path.exists(yaml_abs) else None,
        "sheets": [],
        "rows_out": {},
        "quarantine": {},
        "decisions": {"header_strategy": "yaml_first"},
    }

    result = {}

    def _a1_start(a1):
        """Args: a1(str|None) -> (row0, col0)

        A1 범위 시작 좌표의 0-기반 (행,열) 인덱스를 반환한다. 범위가 없으면 (0,0).
        """
        import re
        if not a1:
            return 0, 0
        m = re.match(r"([A-Za-z]+)(\d+):([A-Za-z]+)(\d+)", str(a1).strip())
        if not m:
            return 0, 0
        def col_to_idx(col):
            col = col.upper()
            acc = 0
            for c in col:
                acc = acc * 26 + (ord(c) - ord('A') + 1)
            return acc - 1
        c1, r1, c2, r2 = m.groups()
        return int(r1) - 1, col_to_idx(c1)
    sheets = cfg.get("sheets", []) or []
    for s in sheets:
        sheet_name = s.get("name")
        if sheet_name is None:
            continue
        manifest["sheets"].append(sheet_name)

        for t in s.get("tables", []) or []:
            tid = t.get("id", "table")
            a1 = t.get("range")
            header_off = int(t.get("header_row_offset", 0))
            df_all = read_excel(xls_path, sheet_name)
            blk = extract_block(df_all, a1) if a1 else df_all.copy()
            data = apply_header(blk, header_off, nrows=1)
            data = sanitize_strings(data)
            data = _apply_rename(data, t.get("rename", {}))
            # coerce/clean
            def _apply_clean(df, cfg):
                if not cfg:
                    return df
                out = df.copy()
                # 강제 문자열 캐스팅
                for c in (cfg.get("coerce_string") or []):
                    if c in out.columns:
                        out[c] = out[c].map(lambda x: None if x is None else str(x))
                # 값 치환
                repl = cfg.get("replace") or {}
                if repl:
                    out = out.replace(repl)
                return out
            data = _apply_clean(data, t.get("clean") or {})

            # 정규화 규칙(normalize_rules)
            def _apply_normalize_rules(df, rules_cfg):
                if not rules_cfg:
                    return df
                out = df.copy()
                import re
                import datetime as _dt
                def _excel_serial_to_ratio(val):
                    try:
                        s = str(val).strip()
                        if not re.fullmatch(r"\d{4,6}", s):
                            return None
                        days = int(s)
                        base = _dt.datetime(1899, 12, 30)
                        dt = base + _dt.timedelta(days=days)
                        if dt.day == 10:
                            return f"{int(dt.month)}/10"
                    except Exception:
                        return None
                    return None
                def _iso_date_to_ratio(val):
                    try:
                        s = str(val).strip()
                        m = re.match(r"^\s*(\d{4})[-/.](\d{1,2})[-/.](\d{1,2})(?:\s+.*)?$", s)
                        if m and m.group(3) == "10":
                            return f"{int(m.group(2))}/10"
                    except Exception:
                        return None
                    return None
                # rules_cfg: {col: rule_or_list}
                for col, spec in (rules_cfg or {}).items():
                    if col not in out.columns:
                        continue
                    if not isinstance(spec, (list, tuple)):
                        spec = [spec]
                    def _normalize_one(v):
                        if v is None:
                            return v
                        cur = str(v)
                        for rule in spec:
                            if rule is None:
                                continue
                            # 내장 규칙
                            if rule == "excel_serial_day10_to_ratio":
                                r = _excel_serial_to_ratio(cur)
                                if r is not None:
                                    cur = r
                                continue
                            if rule == "iso_date_day10_to_ratio":
                                r = _iso_date_to_ratio(cur)
                                if r is not None:
                                    cur = r
                                continue
                            if rule == "dot_to_slash_10":
                                m = re.match(r"^\s*(\d+)\.(\d+)\s*$", cur)
                                if m:
                                    cur = f"{m.group(1)}/10"
                                continue
                            # 패턴 → 템플릿 규칙: "regex -> template" (re.sub 사용, 전역 치환)
                            if isinstance(rule, str) and "->" in rule:
                                try:
                                    pat, rhs = rule.split("->", 1)
                                    pat = pat.strip()
                                    rhs_raw = rhs
                                    rhs = rhs.strip()
                                    if len(rhs) >= 2 and ((rhs[0] == rhs[-1]) and rhs[0] in ("'", '"')):
                                        rhs_use = rhs[1:-1]
                                    else:
                                        rhs_use = rhs
                                    new_s = re.sub(pat, rhs_use, cur)
                                    cur = new_s
                                except Exception:
                                    pass
                        return cur
                    out[col] = out[col].map(_normalize_one)
                return out

            data = _apply_normalize_rules(data, t.get("normalize_rules") or {})
            # 색상 플래그 주입 (테이블 우선 → 전역)
            flags_cfg_any = t.get("flags_from_color") or cfg.get("flags_from_color") or {}
            # A1 시작 좌표 계산(항상 산출하여 provenance_row에도 사용)
            start_r, start_c = _a1_start(a1)
            if flags_cfg_any:
                # 단일/다중 컬럼 매핑 모두 허용
                entries = []
                if isinstance(flags_cfg_any, dict):
                    entries = [flags_cfg_any]
                elif isinstance(flags_cfg_any, list):
                    entries = flags_cfg_any
                if entries:
                    colors = read_cell_colors(xls_path, sheet_name)
                    data = data.reset_index(drop=True)
                    for ent in entries:
                        if not isinstance(ent, dict):
                            continue
                        col_idx = int(ent.get("col_index", -1))
                        fmap = ent.get("map", {}) or {}
                        if col_idx < 0 or not fmap:
                            continue
                        # 플래그 컬럼 생성
                        for hexval, colname in fmap.items():
                            if colname not in data.columns:
                                data[colname] = False
                        # 데이터 행 루프
                        for i in range(len(data)):
                            abs_row = start_r + header_off + 1 + i
                            rgb = colors.get((abs_row, col_idx))
                            if rgb and rgb in fmap:
                                colname = fmap[rgb]
                                try:
                                    data.at[i, colname] = True
                                except Exception:
                                    pass

            data = data.copy()
            data["provenance_file"] = os.path.relpath(xls_path, root)
            data["provenance_sheet"] = sheet_name
            idx = data.reset_index().index.astype(int)
            data["provenance_row"] = (start_r + header_off + 1) + 1 + idx

            data_sorted = stable_sort(data, ["compound_id", "assay_id", "provenance_row"])

            ok, fails = bronze_checks(data_sorted, ["provenance_file", "provenance_sheet", "provenance_row"])
            if len(fails) > 0:
                qpath = write_quarantine("bronze", year, f"{tid}_missing_provenance", fails)
                manifest["quarantine"][tid] = qpath

            tdir = os.path.join(out_dir, "tables")
            os.makedirs(tdir, exist_ok=True)
            tpath = os.path.join(tdir, f"{tid}.csv")
            ok.to_csv(tpath, index=False, encoding="utf-8")
            result[f"table_{tid}"] = tpath
            manifest["rows_out"][f"table_{tid}"] = int(len(ok))

            melt_cfg = t.get("melt") or {}
            if melt_cfg:
                id_cols = list(melt_cfg.get("id_cols", []))
                value_cols = list(melt_cfg.get("value_cols", []))
                flag_cols = [c for c in data.columns if str(c).startswith("flag_")]
                for fc in flag_cols:
                    if fc not in id_cols:
                        id_cols.append(fc)
                if id_cols and value_cols:
                    long_df = melt_table(ok, id_cols, value_cols)
                    long_df = stable_sort(long_df, id_cols + ["provenance_row"])
                    lpath = os.path.join(out_dir, f"{tid}_long.csv")
                    long_df.to_csv(lpath, index=False, encoding="utf-8")
                    result[f"{tid}_long"] = lpath
                    manifest["rows_out"][f"{tid}_long"] = int(len(long_df))

        m = s.get("matrix")
        if m:
            a1 = m.get("range")
            df_all = read_excel(xls_path, sheet_name)
            blocks = []
            matrix_cfg = m
            detect = bool(matrix_cfg.get("detect_blocks", False))
            if detect:
                try:
                    id_col = int(matrix_cfg.get("id_col", 0))
                    first_col = df_all.iloc[:, id_col].map(lambda x: str(x).strip().lower() if x is not None else "")
                    starts = [i for i, v in enumerate(first_col.tolist()) if v.startswith("table ")]
                    if starts:
                        for i, st in enumerate(starts):
                            en = starts[i+1] if i+1 < len(starts) else len(df_all)
                            blocks.append(df_all.iloc[st:en, :])
                except Exception:
                    blocks = []
            if not blocks:
                blk = extract_block(df_all, a1) if a1 else df_all
                blocks = [blk]

            from pipeline.parsers.engine import matrix_to_long
            longs = []
            for blk in blocks:
                ldf = matrix_to_long(blk, matrix_cfg)
                if ldf is None or len(ldf) == 0:
                    continue
                longs.append(ldf)
            import pandas as _pd
            long_df = _pd.concat(longs, ignore_index=True) if longs else _pd.DataFrame(columns=["row_id","panel","cell_line","value_raw","qualifier","value","unit"]) 
            long_df["provenance_file"] = os.path.relpath(xls_path, root)
            long_df["provenance_sheet"] = sheet_name
            long_df["provenance_row"] = None
            long_df = stable_sort(long_df, ["row_id", "cell_line"]) 
            mdir = os.path.join(out_dir, "matrix")
            os.makedirs(mdir, exist_ok=True)
            mpath = os.path.join(mdir, f"{m.get('id','matrix')}_long.csv")
            long_df.to_csv(mpath, index=False, encoding="utf-8")
            result[f"matrix_{m.get('id','matrix')}_long"] = mpath
            manifest["rows_out"][f"matrix_{m.get('id','matrix')}_long"] = int(len(long_df))

    write_manifest(year, manifest)
    try:
        if str(year) == "2017":
            tdir = os.path.join(out_dir, "tables")
            for old, new in [("table1.csv","table_1.csv"),("table2.csv","table_2.csv")]:
                op = os.path.join(tdir, old)
                np = os.path.join(tdir, new)
                if os.path.exists(op) and os.path.exists(np):
                    os.remove(op)
            t3 = os.path.join(tdir, "table_3.csv")
            if os.path.exists(t3):
                try:
                    import pandas as _pd
                    dfc = _pd.read_csv(t3, dtype=str)
                    dfc = stable_sort(dfc, ["compound_id","provenance_row"]) 
                    comp_out = os.path.join(tdir, "compounds.csv")
                    dfc.to_csv(comp_out, index=False, encoding="utf-8")
                    result["table_compounds_alias"] = comp_out
                except Exception:
                    pass
    except Exception:
        pass
    return result

import os
import yaml
import pandas as pd

from pipeline.parsers.engine import read_excel, extract_block, apply_header, melt_table, sanitize_strings, file_sha256, read_cell_colors
from pipeline.validate.hooks import stable_sort, bronze_checks, write_quarantine, write_manifest


def _repo_root():
    """Args: None -> str

    레포 루트(hoon)를 반환.
    """
    here = os.path.dirname(os.path.abspath(__file__))
    return os.path.abspath(os.path.join(here, os.pardir, os.pardir))


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
    out_dir = out_dir or os.path.join("data", "bronze", year)
    os.makedirs(out_dir, exist_ok=True)

    with open(yaml_path, encoding="utf-8") as f:
        cfg = yaml.safe_load(f)

    xls_path = cfg.get("file")
    if not xls_path:
        raise ValueError("YAML에 file 경로가 없습니다.")
    if not os.path.isabs(xls_path):
        xls_path = os.path.join(root, xls_path)
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
            # 색상 플래그 주입
            flags_cfg = cfg.get("flags_from_color") or {}
            if flags_cfg:
                # A1 시작 좌표 계산
                def _a1_start(a1):
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
                start_r, start_c = _a1_start(a1)
                col_idx = int(flags_cfg.get("col_index", -1))
                fmap = flags_cfg.get("map", {})
                if col_idx >= 0 and isinstance(fmap, dict) and len(fmap) > 0:
                    colors = read_cell_colors(xls_path, sheet_name)
                    # 플래그 컬럼 생성
                    for hexval, colname in fmap.items():
                        if colname not in data.columns:
                            data[colname] = False
                    # 데이터 행 루프
                    data = data.reset_index(drop=True)
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
            data["provenance_row"] = (header_off + 2) + data.reset_index().index.astype(int)

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
                # 플래그가 있으면 id_cols에 포함시켜 propagate
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
            blk = extract_block(df_all, a1) if a1 else df_all
            from pipeline.parsers.engine import matrix_to_long
            long_df = matrix_to_long(blk, m)
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
    return result

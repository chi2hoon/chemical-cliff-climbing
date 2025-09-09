import os
import json
import pandas as pd

from pipeline.parsers.engine import sanitize_strings, file_sha256
from pipeline.validate.hooks import stable_sort, write_manifest, write_quarantine
from pipeline.transforms.normalize import parse_qual_value_unit, convert_unit, ascii_units, strip_all, dash_to_nan


def _repo_root():
    here = os.path.dirname(os.path.abspath(__file__))
    # repo 루트: .../base/pipeline → 두 단계 상위
    return os.path.abspath(os.path.join(here, os.pardir, os.pardir))


def _load_yaml(yaml_path):
    import yaml
    with open(yaml_path, encoding="utf-8") as f:
        return yaml.safe_load(f)


def _ensure_dir(p):
    d = os.path.dirname(os.path.abspath(p))
    if d and not os.path.exists(d):
        os.makedirs(d, exist_ok=True)


def _year_refined_dir(year):
    return os.path.join("base", "data", "silver", str(year))


def build_silver(year, yaml_path=None):
    """Args: year(str|int), yaml_path(str|None) -> dict

    ingest 산출을 공통 스키마로 표준화하여 silver 산출을 생성.
    출력: compounds_silver.csv, assay_readings_silver.csv
    """
    year = str(year)
    yaml_path = yaml_path or os.path.join("schemas", "silver", f"{year}.yaml")
    cfg = _load_yaml(yaml_path)
    root = _repo_root()

    yaml_abs = yaml_path if os.path.isabs(yaml_path) else os.path.join(root, yaml_path)
    src_file = (cfg.get("file") or "")
    src_abs = src_file
    if src_file and not os.path.isabs(src_file):
        cand1 = os.path.join(root, src_file)
        cand2 = os.path.join(os.path.abspath(os.path.join(root, os.pardir)), src_file)
        src_abs = cand1 if os.path.exists(cand1) else cand2
        # 호환: base/data 경로가 없고 hoon/data에 존재하면 폴백
        if not os.path.exists(src_abs) and "/base/data/" in ("/" + src_file):
            alt = os.path.join(root, src_file.replace("base/data/", "hoon/data/"))
            if os.path.exists(alt):
                src_abs = alt
    manifest = {
        "stage": "silver",
        "year": year,
        "yaml_path": os.path.relpath(yaml_abs, root) if os.path.exists(yaml_abs) else yaml_path,
        "yaml_sha256": file_sha256(yaml_abs) if os.path.exists(yaml_abs) else None,
        "input_file": os.path.relpath(src_abs, root) if os.path.exists(src_abs) else src_file,
        "input_sha256": file_sha256(src_abs) if os.path.exists(src_abs) else None,
        "rows_out": {},
        "quarantine": {},
    }

    outdir = _year_refined_dir(year)
    os.makedirs(outdir, exist_ok=True)
    comp_out = os.path.join(outdir, "compounds_silver.csv")
    assays_out = os.path.join(outdir, "assay_readings_silver.csv")

    silver_cfg = cfg.get("silver", {})
    comp_cfg = (silver_cfg.get("compounds") or {})

    ingest_dir = os.path.join("base", "data", "bronze", year, "tables")
    comp_src = comp_cfg.get("from") or None
    if not comp_src:
        files = []
        if os.path.exists(ingest_dir):
            for name in os.listdir(ingest_dir):
                if name.endswith(".csv"):
                    files.append(os.path.join(ingest_dir, name))
        comp_src = files[0] if files else None
    comp_df = pd.read_csv(comp_src, dtype=str) if comp_src and os.path.exists(comp_src) else pd.DataFrame()
    comp_df = sanitize_strings(comp_df)

    for c in comp_df.columns:
        comp_df[c] = comp_df[c].map(strip_all).map(dash_to_nan)

    props_cols = ["mw", "lcms_text", "nmr_1h_text"]
    props_df = None
    try:
        if os.path.exists(ingest_dir):
            for name in sorted(os.listdir(ingest_dir)):
                if not name.endswith('.csv'):
                    continue
                path = os.path.join(ingest_dir, name)
                try:
                    tdf = pd.read_csv(path, dtype=str)
                except Exception:
                    continue
                tdf = sanitize_strings(tdf)
                for c in tdf.columns:
                    tdf[c] = tdf[c].map(strip_all).map(dash_to_nan)
                if 'compound_id' in tdf.columns and any(col in tdf.columns for col in props_cols):
                    keep = ['compound_id'] + [c for c in props_cols if c in tdf.columns]
                    props_df = tdf[keep].copy()
                    break
    except Exception:
        props_df = None
    if props_df is not None and 'compound_id' in comp_df.columns:
        comp_df = comp_df.merge(props_df, on='compound_id', how='left', suffixes=("", ""))

    flag_cols = [c for c in comp_df.columns if str(c).startswith("flag_")]
    keep_cols = [c for c in ["compound_id", "smiles_raw", "iupac_name", "mw", "lcms_text", "nmr_1h_text", "provenance_file", "provenance_sheet", "provenance_row"] if c in comp_df.columns] + flag_cols
    comp_silver = comp_df[keep_cols].copy() if len(keep_cols) else comp_df.copy()
    comp_silver = stable_sort(comp_silver, ["compound_id", "provenance_row"])
    comp_silver.to_csv(comp_out, index=False, encoding="utf-8")
    manifest["rows_out"]["compounds_silver"] = int(len(comp_silver))

    assay_cfg = (silver_cfg.get("assays") or {})
    rows = []
    def _parse_with_rules(text, parse_cfg):
        from pipeline.transforms.normalize import parse_qual_value_unit
        s = None if text is None else str(text).strip()
        if not parse_cfg:
            return parse_qual_value_unit(s, None)
        rules = (parse_cfg or {}).get("rules") or []
        default_unit = (parse_cfg or {}).get("default_unit")
        for rule in rules:
            if "->" not in rule:
                continue
            pat, rhs = rule.split("->", 1)
            import re
            m = re.search(pat.strip(), s or "")
            if not m:
                continue
            q, v, u = parse_qual_value_unit(s, default_unit)
            kvs = {}
            for part in rhs.split(','):
                part = part.strip()
                if not part or ":" not in part:
                    continue
                k, vexpr = part.split(":", 1)
                k = k.strip(); vexpr = vexpr.strip()
                if vexpr.startswith("'") and vexpr.endswith("'"):
                    kvs[k] = vexpr[1:-1]
                elif vexpr.startswith("\"") and vexpr.endswith("\""):
                    kvs[k] = vexpr[1:-1]
                elif vexpr.startswith("\\") and vexpr[1:].isdigit():
                    idx = int(vexpr[1:]); kvs[k] = m.group(idx)
                else:
                    kvs[k] = vexpr
            return kvs.get("qualifier", q), kvs.get("value", v), kvs.get("unit", u)
        return parse_qual_value_unit(s, default_unit)

    if assay_cfg.get("from_melt"):
        melt_name = assay_cfg["from_melt"]
        melt_path = os.path.join("base", "data", "bronze", year, f"{melt_name}.csv")
        df = pd.read_csv(melt_path, dtype=str)
        df = sanitize_strings(df)
        for c in df.columns:
            df[c] = df[c].map(strip_all).map(dash_to_nan)
        rules = assay_cfg.get("variable_rules", {})
        target_id_default = assay_cfg.get("target_id")
        for _, r in df.iterrows():
            var = str(r.get("variable", ""))
            rule = rules.get(var, {})
            assay_id = rule.get("assay_id") or assay_cfg.get("assay_id")
            target_id = rule.get("target_id") or target_id_default
            parse_cfg = rule.get("parse") or assay_cfg.get("parse") or {}
            q, v_raw, u_in = _parse_with_rules(r.get("value"), parse_cfg)
            unit_std = ascii_units((parse_cfg or {}).get("default_unit") or u_in)
            v_std, u_std = convert_unit(v_raw, u_in, unit_std)
            row = {
                "compound_id": r.get("compound_id"),
                "target_id": target_id,
                "assay_id": assay_id,
                "qualifier": q,
                "value_std": v_std,
                "unit_std": u_std,
                "year": year,
                "qc_pass": True,
                "provenance_file": r.get("provenance_file"),
                "provenance_sheet": r.get("provenance_sheet"),
                "provenance_row": r.get("provenance_row"),
            }
            for fc in ["flag_asterisk", "flag_imaging_conflict"]:
                if fc in r.index:
                    row[fc] = r.get(fc)
            rows.append(row)

    assay_cols = ["compound_id","target_id","assay_id","qualifier","value_std","unit_std","year","qc_pass","provenance_file","provenance_sheet","provenance_row"]
    assay_df = pd.DataFrame(rows, columns=assay_cols) if rows else pd.DataFrame(columns=assay_cols)
    assay_df = stable_sort(assay_df, ["compound_id","assay_id","provenance_row"])

    meta_blocks = (silver_cfg.get("meta") or [])
    meta_rows = []
    def _coerce_na(x):
        return None if x is None or str(x).strip().lower() in {"nan", ""} else x
    if meta_blocks:
        for mb in meta_blocks:
            src = mb.get("from")
            if not src:
                continue
            path = src
            if not os.path.isabs(path):
                path = os.path.join(_repo_root(), src)
            if not os.path.exists(path):
                alt = os.path.join("base","data","bronze",year,"tables", os.path.basename(src))
                path = alt if os.path.exists(alt) else path
            try:
                mdf = pd.read_csv(path, dtype=str)
            except Exception:
                continue
            mdf = sanitize_strings(mdf)
            inj = mb.get("inject") or {}
            parse = mb.get("parse") or {}
            for c in mdf.columns:
                mdf[c] = mdf[c].map(_coerce_na)
            for k, v in inj.items():
                mdf[k] = v
            from pipeline.transforms.normalize import parse_dose, parse_ratio, parse_percent
            dose_cfg = parse.get("dose") or {}
            if dose_cfg:
                col = dose_cfg.get("col")
                into = dose_cfg.get("into") or {}
                if col in mdf.columns:
                    res = mdf[col].map(lambda s: parse_dose(s))
                    mdf[into.get("dose_route","dose_route")] = res.map(lambda x: x[0] if x else None)
                    mdf[into.get("dose_value","dose_value")] = res.map(lambda x: x[1] if x else None)
                    mdf[into.get("dose_unit","dose_unit")] = res.map(lambda x: x[2] if x else None)
            ratio_cfg = parse.get("ratio") or {}
            if ratio_cfg:
                col = ratio_cfg.get("col")
                into = ratio_cfg.get("into") or {}
                if col in mdf.columns:
                    res = mdf[col].map(lambda s: parse_ratio(s))
                    mdf[into.get("mortality_pct","mortality_pct")] = res.map(lambda x: x[0] if x else None)
                    mdf[into.get("mortality_n","mortality_n")] = res.map(lambda x: x[1] if x else None)
                    mdf[into.get("mortality_k","mortality_k")] = res.map(lambda x: x[2] if x else None)
                    def _mk_ratio(t):
                        if not t:
                            return None
                        n, k = t[1], t[2]
                        return (f"{n}/{k}" if (n and k) else None)
                    ratio_col = into.get("mortality_ratio","mortality_ratio")
                    mdf[ratio_col] = res.map(_mk_ratio)
                    def _needs_rewrite(s):
                        if s is None:
                            return False
                        ss = str(s)
                        return ('/' not in ss) and (ss.strip() != '')
                    mask = mdf[col].map(_needs_rewrite)
                    if mask.any():
                        mdf.loc[mask, col] = mdf.loc[mask, ratio_col]
            pct_cfg = parse.get("percent") or {}
            if pct_cfg:
                col = pct_cfg.get("col")
                into = pct_cfg.get("into") or {}
                if col in mdf.columns:
                    mdf[into.get("hr_change_pct","hr_change_pct")] = mdf[col].map(lambda s: parse_percent(s))
            cid_cfg = parse.get("compound_id_from_text") or {}
            if cid_cfg:
                col = cid_cfg.get("col")
                into = cid_cfg.get("into") or {}
                keep_as = cid_cfg.get("keep_original_as")
                regex = cid_cfg.get("regex") or r"(\d+)"
                import re
                if col in mdf.columns:
                    if keep_as:
                        mdf[keep_as] = mdf[col]
                    mdf[into.get("compound_id","compound_id")] = mdf[col].map(lambda s: (re.search(regex, s).group(1) if (s and re.search(regex, s)) else None))
            meta_rows.append(mdf)

    meta_out_path = None
    if meta_rows:
        meta_df = pd.concat(meta_rows, ignore_index=True)
        meta_df = stable_sort(meta_df, ["compound_id","provenance_row"]) if "compound_id" in meta_df.columns else meta_df
        meta_out_path = os.path.join(outdir, "assay_context_silver.csv")
        meta_df.to_csv(meta_out_path, index=False, encoding="utf-8")

    if assay_cfg.get("from_matrix"):
        fm = assay_cfg.get("from_matrix")
        if isinstance(fm, dict):
            fm = [fm]
        raw_file = cfg.get("file")
        if raw_file and not os.path.isabs(raw_file):
            raw_file = os.path.join(_repo_root(), raw_file)
        ases = (cfg.get("silver") or {}).get("assays") or {}
        target_map = ases.get("target_map", {})
        cell_line_map = ases.get("cell_line_map", {})
        for spec in fm:
            sheet = spec.get("sheet") or spec.get("name")
            if not sheet:
                continue
            from pipeline.parsers.engine import read_excel, extract_block, matrix_to_long
            df_all = read_excel(raw_file, sheet)
            a1 = (spec.get("matrix") or {}).get("range")
            matrix_cfg = spec.get("matrix") or {}
            detect = bool(matrix_cfg.get("detect_blocks", False))
            blocks = []
            if detect:
                id_col = int(matrix_cfg.get("id_col", 0))
                first_col = df_all.iloc[:, id_col].map(lambda x: str(x).strip().lower() if x is not None else "")
                starts = [i for i, v in enumerate(first_col.tolist()) if v.startswith("table ")]
                if starts:
                    for i, st in enumerate(starts):
                        en = starts[i+1] if i+1 < len(starts) else len(df_all)
                        blocks.append(df_all.iloc[st:en, :])
            if not blocks:
                blk = extract_block(df_all, a1) if a1 else df_all
                blocks = [blk]
            longs = []
            for blk in blocks:
                ldf = matrix_to_long(blk, matrix_cfg)
                if ldf is None or len(ldf) == 0:
                    continue
                longs.append(ldf)
            import pandas as _pd
            long_df = _pd.concat(longs, ignore_index=True) if longs else _pd.DataFrame(columns=["row_id","panel","cell_line","value_raw","qualifier","value","unit"]) 
            if not long_df.empty:
                def _norm_cell(x):
                    x = str(x)
                    return cell_line_map.get(x, x)
                long_df["cell_line"] = long_df["cell_line"].map(_norm_cell)
                def _derive_target(r):
                    cl = str(r.get("cell_line") or "").strip()
                    pid = str(r.get("panel") or "").strip()
                    m = target_map.get(cl)
                    if m:
                        return m
                    cln = cl.replace(" ", "").replace("/", "-")
                    pidn = pid.replace(" ", "").replace("/", "-")
                    return f"cell:{pidn}.{cln}" if cl else "cell:unknown"
                long_df["target_id"] = long_df.apply(_derive_target, axis=1)
                def _parse_val(text):
                    q, v, u = parse_qual_value_unit(text, default_unit=(spec.get("value_parser") or {}).get("unit_default"))
                    v2, u2 = convert_unit(v, u, ascii_units((spec.get("value_parser") or {}).get("unit_default")))
                    return q, v2, u2
                long_df[["qualifier","value_std","unit_std"]] = long_df["value_raw"].map(lambda s: _parse_val(s)).apply(pd.Series)
                long_df["assay_id"] = ases.get("assay_id") or "cell.cytotoxicity.ic50"
                long_df["year"] = year
                long_df["qc_pass"] = True
                # row_id가 곧 compound_id이므로 보존 후 매핑
                long_df["compound_id"] = long_df.get("row_id")
                for c in ["row_id", "panel", "cell_line", "value_raw"]:
                    if c in long_df.columns:
                        del long_df[c]
                longs_df = long_df[["target_id","qualifier","value_std","unit_std","assay_id","year","compound_id"]].copy()
                longs_df["provenance_file"] = os.path.relpath(raw_file, _repo_root())
                longs_df["provenance_sheet"] = sheet
                longs_df["provenance_row"] = None
                rows.extend(longs_df.to_dict("records"))

    assay_df2 = pd.DataFrame(rows, columns=assay_cols) if rows else pd.DataFrame(columns=assay_cols)
    assay_df2 = stable_sort(assay_df2, ["compound_id","assay_id","provenance_row"])

    assay_df = pd.concat([assay_df, assay_df2], ignore_index=True) if not assay_df2.empty else assay_df

    assay_df.to_csv(assays_out, index=False, encoding="utf-8")
    write_manifest(year, manifest)

    return {"compounds_silver": comp_out, "assay_readings_silver": assays_out}

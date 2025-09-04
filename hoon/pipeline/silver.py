import os
import json
import pandas as pd

from pipeline.parsers.engine import sanitize_strings, file_sha256
from pipeline.validate.hooks import stable_sort, write_manifest, write_quarantine
from pipeline.transforms.normalize import parse_qual_value_unit, convert_unit, ascii_units, strip_all, dash_to_nan


def _repo_root():
    here = os.path.dirname(os.path.abspath(__file__))
    return os.path.abspath(os.path.join(here, os.pardir))


def _load_yaml(yaml_path):
    import yaml
    with open(yaml_path, encoding="utf-8") as f:
        return yaml.safe_load(f)


def _ensure_dir(p):
    d = os.path.dirname(os.path.abspath(p))
    if d and not os.path.exists(d):
        os.makedirs(d, exist_ok=True)


def _year_refined_dir(year):
    return os.path.join("data", "silver", str(year))


def build_silver(year, yaml_path=None):
    """Args: year(str|int), yaml_path(str|None) -> dict

    ingest 산출을 공통 스키마로 표준화하여 silver 산출을 생성.
    출력: compounds_silver.csv, assay_readings_silver.csv
    """
    year = str(year)
    yaml_path = yaml_path or os.path.join("schemas", "silver", f"{year}.yaml")
    cfg = _load_yaml(yaml_path)
    root = _repo_root()

    # 입력/스키마 경로/해시 수집
    yaml_abs = yaml_path if os.path.isabs(yaml_path) else os.path.join(root, yaml_path)
    src_file = (cfg.get("file") or "")
    src_abs = src_file
    if src_file and not os.path.isabs(src_file):
        src_abs = os.path.join(root, src_file)
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

    # 2017도 YAML 기반 규칙으로 수행 (브릿지 제거)

    # compounds
    outdir = _year_refined_dir(year)
    os.makedirs(outdir, exist_ok=True)
    comp_out = os.path.join(outdir, "compounds_silver.csv")
    assays_out = os.path.join(outdir, "assay_readings_silver.csv")

    silver_cfg = cfg.get("silver", {})
    comp_cfg = (silver_cfg.get("compounds") or {})

    # 기본: ingest/tables/* 를 받아 최소 컬럼만 정리
    ingest_dir = os.path.join("data", "bronze", year, "tables")
    # 우선순위 파일: comp_cfg.file 또는 첫 csv
    comp_src = comp_cfg.get("from") or None
    if not comp_src:
        # pick first csv in ingest/tables
        files = []
        if os.path.exists(ingest_dir):
            for name in os.listdir(ingest_dir):
                if name.endswith(".csv"):
                    files.append(os.path.join(ingest_dir, name))
        comp_src = files[0] if files else None
    comp_df = pd.read_csv(comp_src, dtype=str) if comp_src and os.path.exists(comp_src) else pd.DataFrame()
    comp_df = sanitize_strings(comp_df)

    # 정리: 공백 제거, '-' → None
    for c in comp_df.columns:
        comp_df[c] = comp_df[c].map(strip_all).map(dash_to_nan)

    # 추가 성질 병합: ingest/tables 내에 mw/lcms_text/nmr_1h_text가 있는 테이블과 병합
    # - 우선 YAML로 지정된 테이블이 있다면 사용하고, 없으면 휴리스틱으로 탐색
    props_cols = ["mw", "lcms_text", "nmr_1h_text"]
    props_df = None
    try:
        # 후보 파일 탐색
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
                # 공백/대시 처리
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

    # 출력 컬럼 최소화 (원본 보존 + provenance)
    keep_cols = [c for c in ["compound_id", "smiles_raw", "iupac_name", "mw", "lcms_text", "nmr_1h_text", "provenance_file", "provenance_sheet", "provenance_row"] if c in comp_df.columns]
    comp_silver = comp_df[keep_cols].copy() if len(keep_cols) else comp_df.copy()
    comp_silver = stable_sort(comp_silver, ["compound_id", "provenance_row"])
    comp_silver.to_csv(comp_out, index=False, encoding="utf-8")
    manifest["rows_out"]["compounds_silver"] = int(len(comp_silver))

    # assays
    assay_cfg = (silver_cfg.get("assays") or {})
    rows = []
    def _parse_with_rules(text, parse_cfg):
        # 규칙 기반 파서: 규칙이 없으면 기본 파서로 처리
        from pipeline.transforms.normalize import parse_qual_value_unit
        s = None if text is None else str(text).strip()
        if not parse_cfg:
            return parse_qual_value_unit(s, None)
        rules = (parse_cfg or {}).get("rules") or []
        default_unit = (parse_cfg or {}).get("default_unit")
        # 정규식 우선
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
        melt_path = os.path.join("data", "bronze", year, f"{melt_name}.csv")
        df = pd.read_csv(melt_path, dtype=str)
        df = sanitize_strings(df)
        # 공백/대시 처리
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
            # 색상 플래그 propagate
            for fc in ["flag_asterisk", "flag_imaging_conflict"]:
                if fc in r.index:
                    row[fc] = r.get(fc)
            rows.append(row)

    assay_cols = ["compound_id","target_id","assay_id","qualifier","value_std","unit_std","year","qc_pass","provenance_file","provenance_sheet","provenance_row"]
    assay_df = pd.DataFrame(rows, columns=assay_cols) if rows else pd.DataFrame(columns=assay_cols)
    assay_df = stable_sort(assay_df, ["compound_id","assay_id","provenance_row"])

    # meta 처리: assay_context 등
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
                # 상대경로가 아닌 경우 bronze/{year}/tables 기준으로도 시도
                alt = os.path.join("data","bronze",year,"tables", os.path.basename(src))
                path = alt if os.path.exists(alt) else path
            try:
                mdf = pd.read_csv(path, dtype=str)
            except Exception:
                continue
            mdf = sanitize_strings(mdf)
            inj = mb.get("inject") or {}
            parse = mb.get("parse") or {}
            # 공백/대시 처리
            for c in mdf.columns:
                mdf[c] = mdf[c].map(_coerce_na)
            # 주입 상수 적용
            for k, v in inj.items():
                mdf[k] = v
            # 파싱 규칙
            from pipeline.transforms.normalize import parse_dose, parse_ratio, parse_percent
            # dose
            dose_cfg = parse.get("dose") or {}
            if dose_cfg:
                col = dose_cfg.get("col")
                into = dose_cfg.get("into") or {}
                if col in mdf.columns:
                    res = mdf[col].map(lambda s: parse_dose(s))
                    mdf[into.get("dose_route","dose_route")] = res.map(lambda x: x[0] if x else None)
                    mdf[into.get("dose_value","dose_value")] = res.map(lambda x: x[1] if x else None)
                    mdf[into.get("dose_unit","dose_unit")] = res.map(lambda x: x[2] if x else None)
            # ratio
            ratio_cfg = parse.get("ratio") or {}
            if ratio_cfg:
                col = ratio_cfg.get("col")
                into = ratio_cfg.get("into") or {}
                if col in mdf.columns:
                    res = mdf[col].map(lambda s: parse_ratio(s))
                    mdf[into.get("mortality_pct","mortality_pct")] = res.map(lambda x: x[0] if x else None)
                    mdf[into.get("mortality_n","mortality_n")] = res.map(lambda x: x[1] if x else None)
                    mdf[into.get("mortality_k","mortality_k")] = res.map(lambda x: x[2] if x else None)
                    # 사람이 읽기 쉬운 문자열 형태도 제공 (예: 10/10)
                    def _mk_ratio(t):
                        if not t:
                            return None
                        n, k = t[1], t[2]
                        return (f"{n}/{k}" if (n and k) else None)
                    ratio_col = into.get("mortality_ratio","mortality_ratio")
                    mdf[ratio_col] = res.map(_mk_ratio)
                    # 만약 원본 텍스트가 엑셀 날짜 직렬값/ISO 날짜/소수점 표기 등으로 들어왔고
                    # 현재 텍스트에 '/'가 없다면, 사람이 읽기 쉬운 형태로 교체한다.
                    def _needs_rewrite(s):
                        if s is None:
                            return False
                        ss = str(s)
                        return ('/' not in ss) and (ss.strip() != '')
                    mask = mdf[col].map(_needs_rewrite)
                    if mask.any():
                        mdf.loc[mask, col] = mdf.loc[mask, ratio_col]
            # percent
            pct_cfg = parse.get("percent") or {}
            if pct_cfg:
                col = pct_cfg.get("col")
                into = pct_cfg.get("into") or {}
                if col in mdf.columns:
                    mdf[into.get("hr_change_pct","hr_change_pct")] = mdf[col].map(lambda s: parse_percent(s))
            # compound_id_from_text
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
            # text_pass (그대로 보존)
            # 나머지는 그대로 둠
            # provenance 컬럼은 브론즈에서 이미 주입됨
            meta_rows.append(mdf)

    meta_out_path = None
    if meta_rows:
        meta_df = pd.concat(meta_rows, ignore_index=True)
        meta_df = stable_sort(meta_df, ["compound_id","provenance_row"]) if "compound_id" in meta_df.columns else meta_df
        meta_out_path = os.path.join(outdir, "assay_context_silver.csv")
        meta_df.to_csv(meta_out_path, index=False, encoding="utf-8")
    # matrix 기반 정의 처리
    if assay_cfg.get("from_matrix"):
        fm = assay_cfg.get("from_matrix")
        if isinstance(fm, dict):
            fm = [fm]
        # raw 파일 경로: YAML 최상단 file 사용
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
                # 라벨 행("table X")으로 블록을 자동 분할
                id_col = int(matrix_cfg.get("id_col", 0))
                first_col = df_all.iloc[:, id_col].map(lambda x: str(x).strip().lower() if x is not None else "")
                starts = [i for i, v in enumerate(first_col.tolist()) if v.startswith("table ")]
                if starts:
                    for i, st in enumerate(starts):
                        en = starts[i+1] if i+1 < len(starts) else len(df_all)
                        blocks.append(df_all.iloc[st:en, :])
                else:
                    blocks.append(extract_block(df_all, a1) if a1 else df_all)
            else:
                blocks.append(extract_block(df_all, a1) if a1 else df_all)

            for blk in blocks:
                long_df = matrix_to_long(blk, matrix_cfg)
                # provenance
                long_df["provenance_file"] = os.path.relpath(raw_file, _repo_root())
                long_df["provenance_sheet"] = sheet
                long_df["provenance_row"] = None
                # 규칙: qualifier/value/unit → value_std/unit_std
                for _, r in long_df.iterrows():
                    q = r.get("qualifier"); v = r.get("value"); u = r.get("unit")
                    v_std, u_std = convert_unit(v, u, ascii_units(u or "uM"))
                    # target_id 매핑(있으면 우선)
                    cell_line = r.get("cell_line")
                    if cell_line in cell_line_map:
                        cell_line = cell_line_map[cell_line]
                    if target_map and cell_line in target_map:
                        tgt = target_map[cell_line]
                    else:
                        pnl = r.get("panel") or "panel"
                        cln = (cell_line or "").replace(" ", "").replace("/", "-")
                        tgt = f"cell:{pnl}.{cln}" if cln else f"cell:{pnl}"
                    row = {
                        "compound_id": r.get("row_id"),
                        "target_id": spec.get("target_id") or tgt,
                        "assay_id": spec.get("assay_id") or "cell.cytotoxicity.ic50",
                        "qualifier": q,
                        "value_std": v_std,
                        "unit_std": u_std,
                        "year": year,
                        "qc_pass": True,
                        "provenance_file": long_df.get("provenance_file").iloc[0],
                        "provenance_sheet": sheet,
                        "provenance_row": None,
                    }
                    for fc in ["flag_asterisk", "flag_imaging_conflict"]:
                        if fc in long_df.columns:
                            row[fc] = r.get(fc)
                    rows.append(row)
        assay_df = pd.DataFrame(rows, columns=assay_cols) if rows else pd.DataFrame(columns=assay_cols)
        assay_df = stable_sort(assay_df, ["compound_id","assay_id","provenance_row"])        

    # silver 검증 및 격리
    from pipeline.validate.silver_checks import validate_silver
    ok, quarantine = validate_silver(comp_silver, assay_df)
    qcounts = {}
    for name, qdf in quarantine.items():
        if len(qdf) > 0:
            qpath = write_quarantine("silver", year, name, qdf)
            manifest.setdefault("quarantine", {})[name] = qpath
            qcounts[name] = int(len(qdf))

    ok["compounds"].to_csv(comp_out, index=False, encoding="utf-8")
    ok["assays"].to_csv(assays_out, index=False, encoding="utf-8")
    manifest["rows_out"]["assay_readings_silver"] = int(len(ok["assays"]))
    if qcounts:
        manifest["quarantine_counts"] = qcounts

    # manifest 기록
    write_manifest(year, manifest)
    return {"compounds_silver": comp_out, "assay_readings_silver": assays_out}


def _build_silver_2017_bridge():
    """Args: None -> dict

    레거시 hoon 은행의 silver 산출물을 읽어 data/silver/2017 스키마로 매핑한다.
    - measurements_std.csv → assay_readings_silver.csv
    - compounds_canonical.csv(+bronze/compounds.csv) → compounds_silver.csv
    """
    import yaml
    import pandas as pd

    here = os.path.dirname(os.path.abspath(__file__))
    hoon_dir = os.path.abspath(os.path.join(here, os.pardir))
    cfg_path = os.path.join(hoon_dir, "configs", "2017.yml")
    with open(cfg_path, encoding="utf-8") as f:
        cfg = yaml.safe_load(f)

    data_root = os.path.join(hoon_dir, "data")
    # 레거시 브론즈/실버 보장
    try:
        from src.udm.bronze_build_dispatch import build_bronze_from_raw as legacy_bronze
        legacy_bronze(data_root, cfg)
    except Exception:
        pass
    try:
        from src.udm.silver import build_measurements_std as legacy_silver
        from src.udm.silver_smiles import build_canonical_compounds as legacy_smiles
        legacy_silver(data_root, cfg)
        legacy_smiles(data_root, cfg)
    except Exception:
        pass

    # 경로
    meas_src = os.path.join(data_root, "silver", "2017", "measurements_std.csv")
    comp_canon = os.path.join(data_root, "silver", "2017", "compounds_canonical.csv")
    comp_bronze = os.path.join(data_root, "bronze", "2017", "compounds.csv")

    # 로드
    meas = pd.read_csv(meas_src, dtype=str) if os.path.exists(meas_src) else pd.DataFrame()
    canon = pd.read_csv(comp_canon, dtype=str) if os.path.exists(comp_canon) else pd.DataFrame()
    bronze = pd.read_csv(comp_bronze, dtype=str) if os.path.exists(comp_bronze) else pd.DataFrame()

    # ASCII 단위, qualifier
    from pipeline.transforms.normalize import ascii_units
    def _map_qual(censor):
        m = {"gt": ">", "ge": ">", "lt": "<", "le": "<", "eq": "="}
        c = str(censor) if censor is not None else ""
        return m.get(c, "=")

    assay_cols = ["compound_id","target_id","assay_id","qualifier","value_std","unit_std","year","qc_pass","provenance_file","provenance_sheet","provenance_row"]
    # 2017 target_id 매핑 YAML 로드(선택)
    target_map = {}
    syl_yaml = os.path.join("schemas", "silver", "2017.yaml")
    try:
        import yaml as _y
        if os.path.exists(syl_yaml):
            with open(syl_yaml, encoding="utf-8") as f:
                ycfg = _y.safe_load(f) or {}
                target_map = ((ycfg.get("silver") or {}).get("assays") or {}).get("target_map", {}) or {}
    except Exception:
        target_map = {}

    def _derive_target_id(row):
        cl = str(row.get("cell_line") or "").strip()
        if cl in target_map:
            return target_map[cl]
        pid = str(row.get("panel_id") or "").strip()
        cln = cl.replace(" ", "").replace("/", "-")
        return "cell:" + (pid if pid else "panel") + "." + cln if cl else "cell:unknown"

    rows = []
    if len(meas) > 0:
        for _, r in meas.iterrows():
            unit_std = ascii_units(r.get("unit_std") or r.get("unit"))
            rows.append({
                "compound_id": r.get("compound_id"),
                "target_id": _derive_target_id(r),
                "assay_id": "cell.cytotoxicity.ic50",
                "qualifier": _map_qual(r.get("censor")),
                "value_std": r.get("value_std"),
                "unit_std": unit_std,
                "year": "2017",
                "qc_pass": True,
                "provenance_file": r.get("provenance_file"),
                "provenance_sheet": r.get("provenance_sheet"),
                "provenance_row": r.get("provenance_row"),
            })
    assay_df = pd.DataFrame(rows, columns=assay_cols) if rows else pd.DataFrame(columns=assay_cols)
    assay_df = stable_sort(assay_df, ["compound_id","assay_id","provenance_row"])

    # compounds: canonical + provenance from bronze
    comp = canon.copy()
    if len(bronze) > 0:
        merge_cols = [c for c in ["compound_id","provenance_file","provenance_sheet","provenance_row"] if c in bronze.columns]
        if "compound_id" in merge_cols:
            comp = comp.merge(bronze[merge_cols], on="compound_id", how="left")
    keep = [c for c in ["compound_id","smiles_raw","iupac_name","provenance_file","provenance_sheet","provenance_row"] if c in comp.columns]
    comp_df = comp[keep].copy() if keep else comp.copy()
    comp_df = stable_sort(comp_df, ["compound_id","provenance_row"])

    outdir = os.path.join("data", "silver", "2017")
    os.makedirs(outdir, exist_ok=True)
    comp_out = os.path.join(outdir, "compounds_silver.csv")
    assays_out = os.path.join(outdir, "assay_readings_silver.csv")
    comp_df.to_csv(comp_out, index=False, encoding="utf-8")
    assay_df.to_csv(assays_out, index=False, encoding="utf-8")

    # 매니페스트 기록
    from pipeline.validate.hooks import write_manifest
    from pipeline.parsers.engine import file_sha256
    manifest = {
        "stage": "silver",
        "year": "2017",
        "input_file": os.path.relpath(os.path.join(hoon_dir, cfg["paths"]["raw_file"]), hoon_dir),
        "input_sha256": file_sha256(os.path.join(hoon_dir, cfg["paths"]["raw_file"])) if os.path.exists(os.path.join(hoon_dir, cfg["paths"]["raw_file"])) else None,
        "yaml_path": os.path.join("schemas", "silver", "2017.yaml"),
        "yaml_sha256": file_sha256(os.path.join("schemas", "silver", "2017.yaml")) if os.path.exists(os.path.join("schemas", "silver", "2017.yaml")) else None,
        "rows_out": {
            "compounds_silver": int(len(comp_df)),
            "assay_readings_silver": int(len(assay_df)),
        },
        "quarantine": {},
    }
    write_manifest("2017", manifest)

    return {"compounds_silver": comp_out, "assay_readings_silver": assays_out}

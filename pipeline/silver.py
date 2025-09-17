import os
import json
import pandas as pd

from pipeline.parsers.engine import sanitize_strings, file_sha256
from pipeline.validate.hooks import stable_sort, write_manifest, write_quarantine
from pipeline.transforms.normalize import parse_qual_value_unit, convert_unit, ascii_units, strip_all, dash_to_nan


def _repo_root():
    here = os.path.dirname(os.path.abspath(__file__))
    # repo 루트: .../fin_repo/pipeline → 한 단계 상위
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

    yaml_abs = yaml_path if os.path.isabs(yaml_path) else os.path.join(root, yaml_path)
    src_file = (cfg.get("file") or "")
    src_abs = src_file
    if src_file:
        # 절대/상대 모두에 대해 base→data/hoon 폴백 적용
        candidates = []
        if os.path.isabs(src_file):
            candidates.append(src_file)
            if "/base/data/" in src_file:
                candidates.append(src_file.replace("/base/data/", "/data/"))
                candidates.append(src_file.replace("/base/data/", "/hoon/data/"))
        else:
            cand1 = os.path.join(root, src_file)
            candidates.append(cand1)
            # root 상위 폴더까지도 시도(과거 구조 호환)
            candidates.append(os.path.join(os.path.abspath(os.path.join(root, os.pardir)), src_file))
            if "base/data/" in src_file:
                candidates.append(os.path.join(root, src_file.replace("base/data/", "data/")))
                candidates.append(os.path.join(root, src_file.replace("base/data/", "hoon/data/")))
        for c in candidates:
            if os.path.exists(c):
                src_abs = c
                break
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

    ingest_dir = os.path.join("data", "bronze", year, "tables")
    if not os.path.exists(ingest_dir):
        legacy = os.path.join("base", "data", "bronze", year, "tables")
        ingest_dir = legacy if os.path.exists(legacy) else ingest_dir
    comp_src = comp_cfg.get("from") or None
    if comp_src:
        cand = None
        cands = []
        if os.path.isabs(comp_src):
            cands.append(comp_src)
            if "/base/data/" in comp_src or "base/data/" in comp_src:
                cands.append(comp_src.replace("/base/data/", "/data/").replace("base/data/", "data/"))
        else:
            # repo 루트 상대, base→data 폴백, ingest_dir/basename 탐색
            cands.append(os.path.join(root, comp_src))
            if "base/data/" in comp_src:
                cands.append(os.path.join(root, comp_src.replace("base/data/", "data/")))
            cands.append(os.path.join(ingest_dir, os.path.basename(comp_src)))
        for p in cands:
            if p and os.path.exists(p):
                cand = p
                break
        comp_src = cand or comp_src
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
        # 이미 comp_df에 존재하는 컬럼은 중복 병합하지 않도록 필터링
        add_cols = [c for c in props_df.columns if (c == 'compound_id') or (c not in comp_df.columns)]
        comp_df = comp_df.merge(props_df[add_cols], on='compound_id', how='left')

    flag_cols = [c for c in comp_df.columns if str(c).startswith("flag_")]
    keep_cols = [c for c in ["compound_id", "smiles_raw", "iupac_name", "mw", "lcms_text", "nmr_1h_text", "provenance_file", "provenance_sheet", "provenance_row"] if c in comp_df.columns] + flag_cols
    comp_silver = comp_df[keep_cols].copy() if len(keep_cols) else comp_df.copy()
    # 보강: compound_id 누락 시 가능한 소스에서 보충
    if 'compound_id' not in comp_silver.columns:
        if 'compound_label' in comp_df.columns:
            comp_silver.insert(0, 'compound_id', comp_df['compound_label'])
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
        melt_path = os.path.join("data", "bronze", year, f"{melt_name}.csv")
        if not os.path.exists(melt_path):
            legacy = os.path.join("base", "data", "bronze", year, f"{melt_name}.csv")
            melt_path = legacy if os.path.exists(legacy) else melt_path
        df = pd.read_csv(melt_path, dtype=str)
        df = sanitize_strings(df)
        for c in df.columns:
            df[c] = df[c].map(strip_all).map(dash_to_nan)
        rules = assay_cfg.get("variable_rules", {})
        target_id_default = assay_cfg.get("target_id")
        import re
        meta_extra = []
        for _, r in df.iterrows():
            var = str(r.get("variable", ""))
            rule = rules.get(var, {})
            assay_id = rule.get("assay_id") or assay_cfg.get("assay_id")
            target_id = rule.get("target_id") or target_id_default
            parse_cfg = rule.get("parse") or assay_cfg.get("parse") or {}
            raw_val = r.get("value")
            s = str(raw_val) if raw_val is not None else ""
            # 퍼센트 소수(예: 0.44) 또는 44% 표기 → 메타로만 보존하고 assay에서는 제외
            is_pct_symbol = bool(re.match(r"^\s*[\d\.]+\s*%\s*$", s))
            is_pct_decimal = ("e" not in s.lower()) and bool(re.match(r"^\s*0?\.\d+\s*$", s))
            if is_pct_symbol or is_pct_decimal:
                from pipeline.transforms.normalize import parse_percent
                try:
                    pnum = parse_percent(s)
                    if pnum is not None:
                        try:
                            fv = float(pnum)
                            if 0 <= fv <= 1:
                                fv = fv * 100.0
                            pnum = str(fv)
                        except Exception:
                            pass
                    meta_extra.append({
                        "compound_id": r.get("compound_id"),
                        "target_id": target_id,
                        "assay_id": assay_id,
                        "percent_at_20uM": pnum,
                        "concentration_uM": "20",
                        "context_label": "percent_at_20uM",
                        "provenance_file": r.get("provenance_file"),
                        "provenance_sheet": r.get("provenance_sheet"),
                        "provenance_row": r.get("provenance_row"),
                    })
                except Exception:
                    pass
                continue
            q, v_raw, u_in = _parse_with_rules(s, parse_cfg)
            # 파싱된 단위를 우선하고, 없을 때만 default_unit 사용
            unit_std = ascii_units(u_in or (parse_cfg or {}).get("default_unit"))
            v_std, u_std = convert_unit(v_raw, u_in, unit_std)
            # compound_id 보강: ex_no가 있고 compound_id가 비었으면 'Example {ex_no}'로 설정
            comp_id = r.get("compound_id")
            if (comp_id is None or str(comp_id).strip() == "") and ("ex_no" in r.index):
                eno = r.get("ex_no")
                comp_id = (f"Example {eno}" if eno is not None and str(eno).strip() != "" else None)
            row = {
                "compound_id": comp_id,
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
            # 원본 롱에 smiles_raw가 있으면 보존(후속 골드에서 compound_key 유도에 활용)
            if "smiles_raw" in r.index:
                row["smiles_raw"] = r.get("smiles_raw")
            for fc in ["flag_asterisk", "flag_imaging_conflict"]:
                if fc in r.index:
                    row[fc] = r.get(fc)
            rows.append(row)
        # 퍼센트 메타를 파일에 병합(보강: df에서 직접 재검출)
    # 2018 특허(20 uM 단일 농도 % 억제)에서만 퍼센트 맥락을 유도한다.
    try:
        if year == "2018":
            meta_path = os.path.join(outdir, "assay_context_silver.csv")
            # df 기반 재검출
            df2 = df.copy()
            df2["__is_pct_sym__"] = df2["value"].astype(str).str.match(r"^\s*[\d\.]+\s*%\s*$", na=False)
            df2["__is_pct_dec__"] = (~df2["value"].astype(str).str.lower().str.contains("e")) & df2["value"].astype(str).str.match(r"^\s*0?\.\d+\s*$", na=False)
            dfp = df2[df2["__is_pct_sym__"] | df2["__is_pct_dec__"]].copy()
            if not dfp.empty or meta_extra:
                from pipeline.transforms.normalize import parse_percent
                if not dfp.empty:
                    # 변수명→타깃/어세이 맵 적용
                    var_to = {}
                    for k, rule in rules.items():
                        var_to[k] = {
                            "target_id": rule.get("target_id") or target_id_default,
                            "assay_id": rule.get("assay_id") or assay_cfg.get("assay_id"),
                        }
                    recs = []
                    for _, rr in dfp.iterrows():
                        s2 = str(rr.get("value") or "")
                        pstr = parse_percent(s2)
                        try:
                            fv = float(pstr) if pstr is not None else None
                            if fv is not None and 0 <= fv <= 1:
                                fv = fv * 100.0
                        except Exception:
                            fv = pstr
                        to = var_to.get(str(rr.get("variable")) , {})
                        recs.append({
                            "compound_id": rr.get("compound_id"),
                            "target_id": to.get("target_id"),
                            "assay_id": to.get("assay_id"),
                            "percent_at_20uM": str(fv) if fv is not None else None,
                            "concentration_uM": "20",
                            "context_label": "percent_at_20uM",
                            "provenance_file": rr.get("provenance_file"),
                            "provenance_sheet": rr.get("provenance_sheet"),
                            "provenance_row": rr.get("provenance_row"),
                        })
                    meta_extra.extend(recs)
                prev = pd.read_csv(meta_path, dtype=str) if os.path.exists(meta_path) else pd.DataFrame()
                meta_all = pd.concat([prev, pd.DataFrame(meta_extra)], ignore_index=True)
                sort_keys = [c for c in ["compound_id","provenance_row"] if c in meta_all.columns]
                if sort_keys:
                    meta_all = stable_sort(meta_all, sort_keys)
                meta_all.to_csv(meta_path, index=False, encoding="utf-8")
    except Exception:
        pass

    assay_cols = ["compound_id","target_id","assay_id","qualifier","value_std","unit_std","year","qc_pass","provenance_file","provenance_sheet","provenance_row","smiles_raw"]
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
                alt = os.path.join("data","bronze",year,"tables", os.path.basename(src))
                if not os.path.exists(alt):
                    legacy = os.path.join("base","data","bronze",year,"tables", os.path.basename(src))
                    alt = legacy if os.path.exists(legacy) else alt
                path = alt if os.path.exists(alt) else path
            try:
                mdf = pd.read_csv(path, dtype=str)
            except Exception:
                continue
            mdf = sanitize_strings(mdf)
            # 필터링 옵션 적용(percent_only, regex)
            fcfg = mb.get("filter") or {}
            if isinstance(fcfg, dict) and len(mdf) > 0:
                # percent_only: 해당 컬럼에 '%'
                po = fcfg.get("percent_only")
                if isinstance(po, str) and po in mdf.columns:
                    def _has_pct(x):
                        if x is None:
                            return False
                        s = str(x)
                        return ('%' in s)
                    mdf = mdf[mdf[po].map(_has_pct)]
                # 일반 정규식 필터
                rx = fcfg.get("regex") or fcfg.get("where_regex")
                if isinstance(rx, dict):
                    col = rx.get("col")
                    pat = rx.get("pattern") or rx.get("regex")
                    neg = bool(rx.get("negate", False))
                    if col in mdf.columns and pat:
                        mask = mdf[col].astype(str).str.contains(pat, regex=True, na=False)
                        mdf = mdf[~mask] if neg else mdf[mask]
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

    # 2020: Asterisk 예외 교정표를 compounds_silver에 반영(NEW SMILES 우선 교체)
    try:
        if year == "2020":
            comp_path = os.path.join(outdir, "compounds_silver.csv")
            if os.path.exists(comp_path):
                comp_fix = pd.read_csv(comp_path, dtype=str)
                # asterisk 소스: bronze/tables/asterisk_exceptions.csv 또는 메타 파일에서 유도
                asterisk_path = os.path.join("data","bronze",year,"tables","asterisk_exceptions.csv")
                if not os.path.exists(asterisk_path):
                    legacy = os.path.join("base","data","bronze",year,"tables","asterisk_exceptions.csv")
                    asterisk_path = legacy if os.path.exists(legacy) else asterisk_path
                aster = None
                if os.path.exists(asterisk_path):
                    aster = pd.read_csv(asterisk_path, dtype=str)
                elif meta_out_path and os.path.exists(meta_out_path):
                    tmp = pd.read_csv(meta_out_path, dtype=str)
                    # asterisk 관련 컬럼이 있는 행만 필터링(NEW/OLD/Yes/No 중 하나라도 존재)
                    cols = [c for c in tmp.columns]
                    has_any = [c for c in cols if ("SMILES" in str(c)) or ("Yes/No" in str(c))]
                    aster = tmp[has_any].copy() if has_any else None
                if aster is not None and len(aster) > 0:
                    # 컬럼 정규화
                    def _pick_col(cands):
                        for c in aster.columns:
                            cl = str(c).strip().lower()
                            for tok in cands:
                                if tok in cl:
                                    return c
                        return None
                    col_ex = None
                    # Example 라벨 컬럼 탐지("Example 123" 패턴 비율로 판단)
                    for c in aster.columns:
                        s = aster[c].astype(str)
                        try:
                            ratio = (s.str.match(r"^\s*Example\s+\d+", na=False)).mean()
                        except Exception:
                            ratio = 0.0
                        if ratio and ratio > 0.2:
                            col_ex = c; break
                    col_yes = _pick_col(["yes/no"])
                    # 'NEW SMILES' 또는 'smiles_new' 모두 허용
                    col_new = _pick_col(["new smiles"]) or ("smiles_new" if any(str(c).strip().lower()=="smiles_new" for c in aster.columns) else None)
                    if col_ex is None and "nan" in [str(c) for c in aster.columns]:
                        col_ex = "nan"
                    # 유효 케이스만 선별
                    cols_keep = [c for c in [col_ex, col_yes, col_new] if c is not None]
                    if cols_keep and col_ex is not None:
                        af = aster[cols_keep].copy()
                        af.columns = ["example_label"] + (["asterisk_yesno"] if col_yes else []) + (["smiles_new"] if col_new else [])
                        # 공백/NA 정리
                        for c in af.columns:
                            af[c] = af[c].astype(str).str.strip().replace({"nan": None, "None": None, "": None})
                        # 매핑 준비
                        comp_fix["flag_asterisk"] = False
                        if "flag_smiles_o3_changed" not in comp_fix.columns:
                            comp_fix["flag_smiles_o3_changed"] = None
                        # join on compound_id == example_label
                        j = comp_fix.merge(af, left_on="compound_id", right_on="example_label", how="left", suffixes=("","_ast"))
                        # NEW SMILES 적용
                        mask_new = j.get("smiles_new").notna() if "smiles_new" in j.columns else False
                        if isinstance(mask_new, bool):
                            # no-op when column missing
                            mask_new = j.index == -1  # all False
                        j.loc[mask_new, "smiles_raw"] = j.loc[mask_new, "smiles_new"]
                        j.loc[mask_new, "flag_asterisk"] = True
                        j.loc[mask_new, "flag_smiles_o3_changed"] = True
                        # Yes/No만 있고 NEW 없음: Yes면 검수완료로 flag_asterisk만 표시
                        if "asterisk_yesno" in j.columns:
                            yn = j["asterisk_yesno"].astype(str).str.strip().str.upper()
                            mask_yes = (yn == "Y") | (yn == "YES")
                            mask_no = (yn == "N") | (yn == "NO")
                            j.loc[mask_yes & ~mask_new, "flag_asterisk"] = True
                            # No인데 NEW도 없으면 교정 필요(플래그만 남기고 SMILES는 유지)
                            j.loc[mask_no & ~mask_new, "flag_asterisk"] = True
                            # 변경 여부는 명시 불가 → None 유지
                        # 정리 및 저장
                        keep = [c for c in comp_fix.columns]
                        comp_fix2 = j[keep].copy()
                        comp_fix2 = stable_sort(comp_fix2, ["compound_id","provenance_row"]) if "compound_id" in comp_fix2.columns else comp_fix2
                        comp_fix2.to_csv(comp_path, index=False, encoding="utf-8")
    except Exception:
        pass

    rows_matrix = []
    if assay_cfg.get("from_matrix"):
        fm = assay_cfg.get("from_matrix")
        if isinstance(fm, dict):
            fm = [fm]
        raw_file = cfg.get("file")
        if raw_file:
            rf_candidates = []
            if os.path.isabs(raw_file):
                rf_candidates.append(raw_file)
                if "/base/data/" in raw_file:
                    rf_candidates.append(raw_file.replace("/base/data/", "/data/"))
                    rf_candidates.append(raw_file.replace("/base/data/", "/hoon/data/"))
            else:
                rf_candidates.append(os.path.join(_repo_root(), raw_file))
                if "base/data/" in raw_file:
                    rf_candidates.append(os.path.join(_repo_root(), raw_file.replace("base/data/", "data/")))
                    rf_candidates.append(os.path.join(_repo_root(), raw_file.replace("base/data/", "hoon/data/")))
            for c in rf_candidates:
                if os.path.exists(c):
                    raw_file = c
                    break
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
                rows_matrix.extend(longs_df.to_dict("records"))

    if rows_matrix:
        assay_df2 = pd.DataFrame(rows_matrix, columns=assay_cols)
        assay_df2 = stable_sort(assay_df2, ["compound_id","assay_id","provenance_row"])
        assay_df = pd.concat([assay_df, assay_df2], ignore_index=True)

    # 퍼센트 단위('%')는 assay_readings에서 제외하고 메타로 이동(20 uM 단일농도 억제율)
    try:
        percent_mask = (assay_df.get("unit_std") == "%") if "unit_std" in assay_df.columns else False
        if isinstance(percent_mask, pd.Series) and percent_mask.any():
            pct_df = assay_df[percent_mask].copy()
            # percent_at_20uM: 0~1 사이 값이면 x100 스케일, 그 외는 그대로 사용
            def _pct_scale(x):
                try:
                    v = float(x)
                    if 0 <= v <= 1:
                        return v * 100.0
                    return v
                except Exception:
                    return None
            pct_df["percent_at_20uM"] = pct_df["value_std"].map(_pct_scale)
            pct_df["concentration_uM"] = "20"
            pct_df["context_label"] = "percent_at_20uM"
            # 메타 파일 병합 저장
            meta_path = os.path.join(outdir, "assay_context_silver.csv")
            if os.path.exists(meta_path):
                try:
                    prev = pd.read_csv(meta_path, dtype=str)
                except Exception:
                    prev = pd.DataFrame()
                meta_all = pd.concat([prev, pct_df], ignore_index=True)
            else:
                meta_all = pct_df
            # 안정 정렬(가능한 경우)
            sort_keys = [c for c in ["compound_id", "provenance_row"] if c in meta_all.columns]
            if sort_keys:
                meta_all = stable_sort(meta_all, sort_keys)
            meta_all.to_csv(meta_path, index=False, encoding="utf-8")
            # assay에서 제거
            assay_df = assay_df[~percent_mask].copy()
    except Exception:
        pass

    # 보강: SMILES 기반으로 compound_id(Example N) 매핑 시도(직접 키가 없는 연도용)
    try:
        if "compound_id" in assay_df.columns and assay_df["compound_id"].isna().any() and "smiles_raw" in assay_df.columns:
            comp_path = os.path.join(outdir, "compounds_silver.csv")
            if os.path.exists(comp_path):
                cdf = pd.read_csv(comp_path, dtype=str)
                if "smiles_raw" in cdf.columns and "compound_id" in cdf.columns:
                    cmap = cdf.dropna(subset=["smiles_raw","compound_id"]).drop_duplicates(subset=["smiles_raw"])[["smiles_raw","compound_id"]]
                    assay_df = assay_df.merge(cmap, on="smiles_raw", how="left", suffixes=("","_by_smiles"))
                    mask_fill = assay_df["compound_id"].isna() & assay_df["compound_id_by_smiles"].notna()
                    if mask_fill.any():
                        assay_df.loc[mask_fill, "compound_id"] = assay_df.loc[mask_fill, "compound_id_by_smiles"]
                    if "compound_id_by_smiles" in assay_df.columns:
                        del assay_df["compound_id_by_smiles"]
    except Exception:
        pass

    assay_df.to_csv(assays_out, index=False, encoding="utf-8")

    # 보강: (2018 전용) from_melt 원천에서 20uM 퍼센트 맥락을 재계산하여 메타 파일로 생성/갱신
    try:
        if year == "2018" and assay_cfg.get("from_melt"):
            melt_name = assay_cfg["from_melt"]
            melt_path = os.path.join("data", "bronze", year, f"{melt_name}.csv")
            if not os.path.exists(melt_path):
                legacy = os.path.join("base", "data", "bronze", year, f"{melt_name}.csv")
                melt_path = legacy if os.path.exists(legacy) else melt_path
            if os.path.exists(melt_path):
                mdf = pd.read_csv(melt_path, dtype=str)
                mdf = sanitize_strings(mdf)
                # 규칙 매핑
                rules = assay_cfg.get("variable_rules", {})
                vmap = {}
                for k, rule in rules.items():
                    vmap[k] = {
                        "target_id": rule.get("target_id") or assay_cfg.get("target_id"),
                        "assay_id": rule.get("assay_id") or assay_cfg.get("assay_id"),
                    }
                # 퍼센트 판정: 기호 또는 지수 미포함 소수
                s = mdf["value"].astype(str)
                mask = s.str.match(r"^\s*[\d\.]+\s*%\s*$", na=False) | ((~s.str.lower().str.contains("e")) & s.str.match(r"^\s*0?\.\d+\s*$", na=False))
                pf = mdf[mask].copy()
                if len(pf) > 0:
                    from pipeline.transforms.normalize import parse_percent
                    def _to_pct(x):
                        p = parse_percent(str(x))
                        try:
                            fv = float(p) if p is not None else None
                            if fv is not None and 0 <= fv <= 1:
                                fv = fv * 100.0
                            return str(fv) if fv is not None else None
                        except Exception:
                            return p
                    pf["percent_at_20uM"] = pf["value"].map(_to_pct)
                    pf["concentration_uM"] = "20"
                    pf["context_label"] = "percent_at_20uM"
                    pf["target_id"] = pf["variable"].map(lambda v: (vmap.get(v) or {}).get("target_id"))
                    pf["assay_id"] = pf["variable"].map(lambda v: (vmap.get(v) or {}).get("assay_id"))
                    keep = [c for c in ["compound_id","target_id","assay_id","percent_at_20uM","concentration_uM","context_label","provenance_file","provenance_sheet","provenance_row"] if c in pf.columns or c in {"percent_at_20uM","concentration_uM","context_label","target_id","assay_id"}]
                    meta_all = pf[keep].copy()
                    meta_all = stable_sort(meta_all, ["compound_id","provenance_row"]) if "compound_id" in meta_all.columns else meta_all
                    meta_path = os.path.join(outdir, "assay_context_silver.csv")
                    meta_all.to_csv(meta_path, index=False, encoding="utf-8")
    except Exception:
        pass

    write_manifest(year, manifest)

    return {"compounds_silver": comp_out, "assay_readings_silver": assays_out}

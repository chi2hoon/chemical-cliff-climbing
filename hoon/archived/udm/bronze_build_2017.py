import os
import re

import pandas as pd

from .io_csv import write_csv_safe
from .parsers import split_censor_and_value, normalize_unit_label, parse_scientific_notation


# --- helpers (localized to 2017 builder) ---
NUMERIC_VALUE_RE = re.compile(r'^(?:\s*(?:>=|<=|[<>＝＝＜＞=])\s*)?[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?(?:\s*(?:μM|uM|µM|nM|mM|mg/kg|%))?\s*$')


def _is_numeric_value(text: object) -> bool:
    if text is None:
        return False
    s = str(text).strip()
    if s == "":
        return False
    return bool(NUMERIC_VALUE_RE.match(s))


def _is_plausible_smiles(text: str | None) -> bool:
    if text is None:
        return False
    val = str(text).strip()
    if val == "":
        return False
    if val.lower() in {"nan", "none", "null", "na", "n/a"}:
        return False
    if (" " in val) or ("," in val):
        return False
    if not re.fullmatch(r"[A-Za-z0-9Hh@=+\-#\\/\\\\\[\]()\.]+", val):
        return False
    return any(ch in val for ch in "BCNOSPFIclonps")


def _detect_header_row(df: pd.DataFrame, max_scan: int = 60) -> int:
    limit = min(max_scan, len(df))
    for i in range(limit):
        vals = [str(x).strip().lower() for x in df.iloc[i].tolist() if str(x).strip() != ""]
        has_comp = any(v in {"compound", "compound #", "compound no", "compound no."} or v.startswith("compound ") for v in vals)
        has_smiles = any(("smiles" in v) or ("structure" in v) or ("구조" in v) for v in vals)
        if has_comp and has_smiles:
            return i
    for i in range(limit):
        row = " ".join(str(x) for x in df.iloc[i].tolist()).lower()
        if "smiles" in row or "compound" in row or "num" in row:
            return i
    return 0


def _find_col(cols: list[str], keywords: list[str]) -> str | None:
    # match allowing whitespace variants (e.g., '화 합물' vs '화합물')
    cols_l = [str(c).strip().lower() for c in cols]
    for kw in keywords:
        kw_l = kw.strip().lower().replace(" ", "")
        for c, lc in zip(cols, cols_l):
            lc_norm = lc.replace(" ", "")
            if kw_l in lc_norm:
                return c
    return None


def _norm_cell_for_id(s: str) -> str:
    s = str(s or "").strip().lower()
    s = re.sub(r"[^a-z0-9]+", "-", s)
    return s.strip("-")


def build_bronze_from_raw(root_dir: str, cfg: dict) -> dict[str, str]:
    outputs: dict[str, str] = {}
    raw_path = os.path.join(root_dir, cfg["paths"]["raw_file"])
    bronze_dir = os.path.join(root_dir, cfg["paths"]["bronze_dir"])
    os.makedirs(bronze_dir, exist_ok=True)

    # accumulators
    comp_rows: dict[int, dict[str, str]] = {}
    assay_rows: dict[str, dict[str, str]] = {}
    meas_rows: list[dict[str, object]] = []
    panel_meta_acc: dict[str, dict[str, object]] = {}

    parsing = cfg.get("parsing", {})
    allowed_compound_sheets = parsing.get("compounds_from_sheets", None)
    whitelist_re = parsing.get("assay_sheets_whitelist_regex")
    cell_map = {**parsing.get("cell_line_map", {}), **parsing.get("cell_synonyms", {})}
    panel_defs = parsing.get("panel_definitions", {})

    # assay patterns exact mapping and allowed cell columns
    col_to_assign: dict[str, dict[str, object]] = {}
    norm_to_assign: dict[str, dict[str, object]] = {}
    allowed_labels: set[str] = set()
    for p in parsing.get("assay_patterns", []):
        assign = dict(p.get("assign", {}))
        for col in p.get("columns", []):
            c = str(col)
            col_to_assign[c] = assign
            norm_to_assign[str(c).strip().lower()] = assign
            allowed_labels.add(c)
    # also allow labels from panel_definitions cell_lines
    for tdef in panel_defs.values():
        for cl in tdef.get("cell_lines", []) or []:
            allowed_labels.add(str(cl))
    allowed_labels_norm = {str(x).strip().lower() for x in allowed_labels}

    xl = pd.ExcelFile(raw_path)

    for sheet_name in xl.sheet_names:
        # compounds collection (표3)
        df_all = xl.parse(sheet_name, header=None)
        # detect all header rows (rows containing '화합물') within the sheet
        header_rows: list[int] = []
        for i in range(len(df_all)):
            vals = [str(x) for x in df_all.iloc[i].tolist()]
            def _norm_txt(s: str) -> str:
                return s.replace(" ", "")
            if any(("화합물" in _norm_txt(v)) or ("compound" in v.lower()) for v in vals):
                header_rows.append(i)
        if not header_rows:
            # fallback to single detection
            header_rows = [_detect_header_row(df_all)]
        # use single header variant for 2017 표5~표14

        # collect compounds only from whitelisted sheets (if provided)
        if (allowed_compound_sheets is None) or (sheet_name in allowed_compound_sheets):
            # use the first detected header for compound sheet
            df_comp = xl.parse(sheet_name, header=header_rows[0])
            df_comp.columns = [str(c).strip() for c in df_comp.columns]
            comp_col = _find_col(list(df_comp.columns), [
                "compound #","compound no","compound","cmpd","num","no.","no","entry","id","화합물","번호"
            ])
            smiles_col = _find_col(list(df_comp.columns), ["smiles","structure","구조"])
            if comp_col and smiles_col:
                for idx, row in df_comp[[comp_col, smiles_col]].dropna(how="all").iterrows():
                    cid_raw = str(row[comp_col]).strip()
                    smi = str(row[smiles_col]).strip()
                    if cid_raw == "" and smi == "":
                        continue
                    try:
                        cid_int = int(cid_raw)
                    except Exception:
                        continue
                    if not _is_plausible_smiles(smi):
                        continue
                    rec = {
                        "compound_id": str(cid_int),
                        "smiles_raw": smi,
                        "dataset_id": str(cfg.get("dataset_id", "")),
                        "provenance_file": cfg["paths"]["raw_file"],
                        "provenance_sheet": sheet_name,
                        "provenance_row": str(int(idx) + 1),
                    }
                    if cid_int not in comp_rows:
                        comp_rows[cid_int] = rec

        # assay ingestion only for whitelisted tables (표5~14)
        if whitelist_re and not re.search(whitelist_re, sheet_name, flags=re.IGNORECASE):
            continue

        # find columns and melt per assay-like column based on YAML patterns or numeric content
        def _process_assay_column(cur_df: pd.DataFrame, comp_col2: str, orig_label: str):
            norm_label = str(orig_label).strip().lower()
            if (orig_label not in col_to_assign) and (norm_label not in norm_to_assign) and (norm_label not in allowed_labels_norm):
                return
            assign = dict(col_to_assign.get(orig_label) or norm_to_assign.get(norm_label) or {})
            if not assign:
                return
            readout = str(assign.get("readout", orig_label))
            matrix = str(assign.get("matrix", "cell"))
            unit_default = assign.get("unit", "μM")
            extras = dict(assign.get("extras", {}))
            table_id = str(extras.get("table_id", "")).strip().lower()
            cell_line = extras.get("cell_line") or orig_label
            if cell_line in cell_map:
                cell_line = cell_map[cell_line]
            panel_info = panel_defs.get(table_id, {}) if table_id else {}
            panel_id = str(panel_info.get("panel_id", "")).strip()
            disease_area = str(panel_info.get("disease_area", "")).strip()
            panel_label = str(panel_info.get("panel_label", "")).strip()
            assay_id = f"assay.cell.cytotoxicity.ic50.um.{panel_id}.{_norm_cell_for_id(cell_line)}" if panel_id else f"assay.cell.cytotoxicity.ic50.um.unknown.{_norm_cell_for_id(cell_line)}"
            if assay_id not in assay_rows:
                assay_rows[assay_id] = {
                    "assay_id": assay_id,
                    "readout": readout,
                    "matrix": matrix,
                    "label_raw": orig_label,
                    "unit": str(unit_default),
                    "panel_id": panel_id,
                    "cell_line": str(cell_line),
                    "disease_area": disease_area,
                    "meas_kind": "IC50",
                }
            if str(orig_label).strip().lower().startswith("unnamed"):
                return
            iter_frame = cur_df[[comp_col2, orig_label]].dropna(how="all")
            for ridx, r in iter_frame.iterrows():
                comp_raw = str(r[comp_col2]).strip()
                val_raw = str(r[orig_label]).strip()
                if comp_raw == "" and val_raw == "":
                    continue
                try:
                    comp_id = str(int(float(comp_raw)))
                except Exception:
                    m = re.search(r"(\d+)", comp_raw)
                    comp_id = str(int(m.group(1))) if m else comp_raw
                if not _is_numeric_value(val_raw):
                    continue
                is_percent = "%" in val_raw
                censor, number_text, unit_tok = split_censor_and_value(val_raw)
                unit_from_cell = normalize_unit_label(unit_tok) if unit_tok else None
                unit_final = "%" if is_percent else (unit_from_cell or unit_default or "")
                val_num = None
                if number_text is not None and not is_percent:
                    val_num, _ = parse_scientific_notation(number_text)
                mrow: dict[str, object] = {
                    "compound_id": comp_id,
                    "assay_id": assay_id,
                    "panel_id": panel_id,
                    "cell_line": str(cell_line),
                    "value_raw": val_raw,
                    "value_num": val_num if val_num is not None else "",
                    "unit": unit_final,
                    "censor": censor or "",
                    "assay_label": orig_label,
                    "provenance_file": cfg["paths"]["raw_file"],
                    "provenance_sheet": sheet_name,
                    "provenance_row": str(int(ridx) + 1),
                    "meas_kind": "IC50" if not is_percent else "percent_inhibition_at_20uM",
                    "table_id": table_id,
                }
                if not re.fullmatch(r"\d+", str(comp_id)):
                    mrow["compound_name_raw"] = comp_raw
                if panel_id:
                    pm = panel_meta_acc.setdefault(panel_id, {
                        "panel_label": panel_label,
                        "disease_area": disease_area,
                        "provenance_sheet": sheet_name,
                        "cell_lines": set(),
                    })
                    try:
                        cast: set[str] = pm["cell_lines"]  # type: ignore
                        cast.add(str(cell_line))
                    except Exception:
                        pass
                meas_rows.append(mrow)

        handled_any = False
        # process every detected table block within the sheet
        for bi, header_row0 in enumerate(header_rows):
            # build a DataFrame slice for block: header at header_row0, rows until next header
            next_row = header_rows[bi + 1] if bi + 1 < len(header_rows) else len(df_all)
            raw_vals = [str(x).strip() for x in df_all.iloc[header_row0].tolist()]
            header_vals = [v for v in raw_vals if v and v.lower() != 'nan']
            data_block = df_all.iloc[header_row0 + 1: next_row, :len(header_vals)].copy()
            # construct block dataframe with cleaned headers only
            cur_df = pd.DataFrame(data_block.values, columns=header_vals)
            cur_df.columns = [str(c).strip() for c in cur_df.columns]
            comp_col2 = _find_col(list(cur_df.columns), [
                "compound #","compound no","compound","cmpd","num","no.","no","entry","id","화합물","번호"
            ])
            if not comp_col2:
                continue
            handled_any = True
            for c in cur_df.columns:
                if c == comp_col2:
                    continue
                _process_assay_column(cur_df, comp_col2, str(c))

        # if no variant handled, skip sheet
        if whitelist_re and not handled_any:
            continue

    # build dataframes
    comp_cols = ["compound_id", "smiles_raw", "dataset_id", "provenance_file", "provenance_sheet", "provenance_row"]
    assay_cols = ["assay_id", "readout", "matrix", "label_raw", "unit", "panel_id", "cell_line", "disease_area", "meas_kind"]
    meas_cols = ["compound_id", "assay_id", "panel_id", "cell_line", "value_raw", "value_num", "unit", "censor",
                 "assay_label", "provenance_file", "provenance_sheet", "provenance_row", "meas_kind", "table_id", "compound_name_raw"]

    comp_df = pd.DataFrame(sorted(comp_rows.values(), key=lambda r: int(r["compound_id"])) if comp_rows else [], columns=comp_cols).drop_duplicates()
    assay_df = pd.DataFrame(list(assay_rows.values()), columns=assay_cols).drop_duplicates()
    meas_df = pd.DataFrame(meas_rows, columns=meas_cols)
    # Deduplicate repeated page segments: keep first occurrence per (compound_id, assay_id)
    if len(meas_df) > 0:
        meas_df = meas_df.drop_duplicates(subset=["compound_id", "assay_id"], keep="first").reset_index(drop=True)

    # write outputs
    comp_path = os.path.join(bronze_dir, "compounds.csv")
    assay_path = os.path.join(bronze_dir, "assays.csv")
    meas_path = os.path.join(bronze_dir, "measurements.csv")
    panel_meta_path = os.path.join(bronze_dir, "panel_block_meta.csv")

    write_csv_safe(comp_df, comp_path)
    write_csv_safe(assay_df, assay_path)
    write_csv_safe(meas_df, meas_path)

    # 2017: panel_block_meta by panel_id
    rows_panel = []
    for pid, meta in panel_meta_acc.items():
        cells = sorted(list(meta.get("cell_lines", []))) if isinstance(meta.get("cell_lines"), (set, list)) else []
        rows_panel.append({
            "panel_id": pid,
            "panel_label": str(meta.get("panel_label", "")),
            "disease_area": str(meta.get("disease_area", "")),
            "provenance_sheet": str(meta.get("provenance_sheet", "")),
            "cell_lines": ";".join(cells),
        })
    df_panel = pd.DataFrame(rows_panel, columns=["panel_id", "panel_label", "disease_area", "provenance_sheet", "cell_lines"]).drop_duplicates()
    write_csv_safe(df_panel, panel_meta_path)

    outputs.update({
        "compounds": comp_path,
        "assays": assay_path,
        "measurements": meas_path,
        "panel_meta": panel_meta_path,
    })
    return outputs


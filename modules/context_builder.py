import os
import json
import math
import pandas as pd

from rdkit import Chem
from rdkit.Chem import Descriptors, Crippen, Lipinski, rdMolDescriptors


def _to_m(val, unit):
    """Args: val(str|float|None), unit(str|None) -> float|None

    단위를 고려해 몰(M)로 환산한다. 알 수 없으면 None.
    """
    if val is None:
        return None
    try:
        v = float(val)
    except Exception:
        return None
    u = (str(unit) if unit is not None else "").strip().replace("µ", "u").lower()
    if u == "m":
        return v
    if u == "mm":
        return v * 1e-3
    if u == "um":
        return v * 1e-6
    if u == "nm":
        return v * 1e-9
    if u == "pm":
        return v * 1e-12
    # 알 수 없으면 uM 가정
    return v * 1e-6


def _to_pact(value, unit):
    """Args: value(str|float|None), unit(str|None) -> float|None

    값과 단위를 받아 pAct(-log10[M])로 변환.
    """
    m = _to_m(value, unit)
    try:
        if m is None or float(m) <= 0:
            return None
        return -math.log10(float(m))
    except Exception:
        return None


def _calc_descriptors(smiles):
    """Args: smiles(str) -> dict

    RDKit 기본 분자 지표 계산.
    """
    if not smiles:
        return {}
    try:
        mol = Chem.MolFromSmiles(str(smiles))
        if mol is None:
            return {}
        return {
            "logP": float(Crippen.MolLogP(mol)),
            "TPSA": float(rdMolDescriptors.CalcTPSA(mol)),
            "HBD": float(Lipinski.NumHDonors(mol)),
            "HBA": float(Lipinski.NumHAcceptors(mol)),
            "rotB": float(Lipinski.NumRotatableBonds(mol)),
            "aromatic_rings": float(rdMolDescriptors.CalcNumAromaticRings(mol)),
        }
    except Exception:
        return {}


def _first_or_none(values):
    for v in values:
        if v is not None and str(v).strip() != "" and str(v).strip().lower() != "nan":
            return v
    return None


def _read_csv_safe(path):
    try:
        if os.path.exists(path):
            return pd.read_csv(path, dtype=str)
    except Exception:
        pass
    return pd.DataFrame()


def _extract_2017_meta(meta_df, cid):
    sub = meta_df[meta_df.get("compound_id", "").astype(str) == str(cid)] if "compound_id" in meta_df.columns else pd.DataFrame()
    if len(sub) == 0:
        return None
    r = sub.iloc[0]
    return {
        "species": r.get("species"),
        "sex": r.get("sex"),
        "dose": {"route": r.get("dose_route"), "value": r.get("dose_value"), "unit": r.get("dose_unit")},
        "mortality": {"pct": r.get("mortality_pct"), "n": r.get("mortality_n"), "k": r.get("mortality_k"), "ratio": r.get("mortality_ratio")},
        "hr_change_pct": r.get("hr_change_pct"),
    }


def _extract_2018_meta(meta_df, cid):
    sub = meta_df[meta_df.get("compound_id", "").astype(str) == str(cid)] if "compound_id" in meta_df.columns else pd.DataFrame()
    if len(sub) == 0:
        return None
    # percent_at_20uM가 있는 첫 행을 선택
    sub2 = sub[sub.get("percent_at_20uM", "").astype(str).str.strip() != ""] if "percent_at_20uM" in sub.columns else sub.iloc[0:0]
    if len(sub2) == 0:
        return None
    r = sub2.iloc[0]
    return {"percent_at_20uM": r.get("percent_at_20uM")}


def _extract_2020_flags(comp_silver_df, cid):
    sub = comp_silver_df[comp_silver_df.get("compound_id", "").astype(str) == str(cid)] if "compound_id" in comp_silver_df.columns else pd.DataFrame()
    if len(sub) == 0:
        return {"applied": False}
    r = sub.iloc[0]
    fa = str(r.get("flag_asterisk") or "").strip().lower()
    fchg = str(r.get("flag_smiles_o3_changed") or "").strip().lower()
    return {
        "applied": (fa in {"true", "y", "yes"}) or (fchg in {"true", "y", "yes"})
    }


def _extract_2021_allvalues(assay_df, cid):
    out = {"enzyme": {}, "cell": {}}
    sub = assay_df[assay_df.get("compound_id", "").astype(str) == str(cid)] if "compound_id" in assay_df.columns else pd.DataFrame()
    if len(sub) == 0:
        return out
    # enzyme Ki A/B (nM in silver → uM in gold)
    def pick(assay_id):
        s = sub[sub.get("assay_id", "") == assay_id]
        if len(s) == 0:
            return None
        r = s.iloc[0]
        return {
            "value": r.get("value_std"),
            "qualifier": r.get("qualifier"),
            "unit": r.get("unit_std"),
        }
    out["enzyme"]["ki_A"] = pick("prmt5.enzyme.ki.A")
    out["enzyme"]["ki_B"] = pick("prmt5.enzyme.ki.B")
    out["cell"]["ic50_null_uM"] = pick("cell.proliferation.ic50")  # target과 무관 저장, 해석 시 target_id 참고
    out["cell"]["ic50_wt_uM"] = None
    # target 기반으로 WT/null 분리
    # null
    s_null = sub[(sub.get("target_id", "").astype(str) == "cell:HCT116.MTAP-null") & (sub.get("assay_id", "").astype(str) == "cell.proliferation.ic50")]
    if len(s_null) > 0:
        r = s_null.iloc[0]
        out["cell"]["ic50_null_uM"] = {"value": r.get("value_std"), "qualifier": r.get("qualifier"), "unit": r.get("unit_std")}
    s_wt = sub[(sub.get("target_id", "").astype(str) == "cell:HCT116.WT") & (sub.get("assay_id", "").astype(str) == "cell.proliferation.ic50")]
    if len(s_wt) > 0:
        r = s_wt.iloc[0]
        out["cell"]["ic50_wt_uM"] = {"value": r.get("value_std"), "qualifier": r.get("qualifier"), "unit": r.get("unit_std")}
    return out


def build_pair_context(year, pair_row, df_base, data_root, selected_axis, scale_used):
    """Args: year(str|int), pair_row(Series), df_base(DataFrame), data_root(str), selected_axis(str), scale_used(str) -> dict

    가설 생성용 컨텍스트 JSON을 구성한다.
    """
    year = str(year)
    # 매칭: SMILES 기반으로 df_base에서 행 추출
    s1 = str(pair_row.get("SMILES_1"))
    s2 = str(pair_row.get("SMILES_2"))
    a1 = pair_row.get("Activity_1")
    a2 = pair_row.get("Activity_2")
    sim = pair_row.get("Similarity")

    base1 = df_base[df_base.get("SMILES", "").astype(str) == s1]
    base2 = df_base[df_base.get("SMILES", "").astype(str) == s2]
    r1 = base1.iloc[0] if len(base1) > 0 else pd.Series(dtype=object)
    r2 = base2.iloc[0] if len(base2) > 0 else pd.Series(dtype=object)

    def _mk_comp(r, smiles, act):
        cid = r.get("compound_id")
        ck = r.get("compound_key")
        unit = r.get("unit_std")
        qual = r.get("qualifier")
        val = r.get("value_std") if r.get("value_std") is not None else r.get("Activity")
        pact = _to_pact(val, unit)
        return {
            "smiles": smiles,
            "activity": {
                "value": val,
                "unit": unit,
                "qualifier": qual,
                "pAct": pact,
                "shown": act,
            },
            "compound_id": cid,
            "compound_key": ck,
            "provenance": {
                "file": r.get("provenance_file"),
                "sheet": r.get("provenance_sheet"),
                "row": r.get("provenance_row"),
            },
        }

    comp_low = _mk_comp(r1, s1, a1)
    comp_high = _mk_comp(r2, s2, a2)

    # assay/target 요약(가능 시)
    assays = list({v for v in [r1.get("assay_id"), r2.get("assay_id")] if v})
    targets = list({v for v in [r1.get("target_id"), r2.get("target_id")] if v})

    # 구조 차이 간단 지표
    d1 = _calc_descriptors(s1)
    d2 = _calc_descriptors(s2)
    dd = {}
    for k in sorted(set(list(d1.keys()) + list(d2.keys()))):
        try:
            dd[k] = None if (k not in d1 or k not in d2) else float(d2[k] - d1[k])
        except Exception:
            dd[k] = None

    ctx = {
        "dataset": {
            "year": year,
            "selected_axis": selected_axis,
            "scale_used": ("pAct" if str(scale_used).lower().startswith("pact") else "raw"),
            "assay_id": assays if assays else None,
            "target_id": targets if targets else None,
        },
        "pair": {
            "similarity": sim,
            "low": comp_low,
            "high": comp_high,
        },
        "structural_diff": {
            "descriptors_low": d1,
            "descriptors_high": d2,
            "descriptors_delta": dd,
            "summary": None,
        },
        "meta": {},
    }

    # 연도별 메타
    gold_meta = _read_csv_safe(os.path.join("data", "gold", year, "assay_context.csv"))
    if len(gold_meta) == 0:
        gold_meta = _read_csv_safe(os.path.join("base", "data", "gold", year, "assay_context.csv"))

    if year == "2017":
        ctx["meta"]["low_in_vivo"] = _extract_2017_meta(gold_meta, comp_low.get("compound_id"))
        ctx["meta"]["high_in_vivo"] = _extract_2017_meta(gold_meta, comp_high.get("compound_id"))

    if year == "2018":
        ctx["meta"]["low_percent_at_20uM"] = _extract_2018_meta(gold_meta, comp_low.get("compound_id"))
        ctx["meta"]["high_percent_at_20uM"] = _extract_2018_meta(gold_meta, comp_high.get("compound_id"))

    if year == "2020":
        comp_silver = _read_csv_safe(os.path.join("data", "silver", year, "compounds_silver.csv"))
        if len(comp_silver) == 0:
            comp_silver = _read_csv_safe(os.path.join("base", "data", "silver", year, "compounds_silver.csv"))
        ctx["meta"]["asterisk_low"] = _extract_2020_flags(comp_silver, comp_low.get("compound_id"))
        ctx["meta"]["asterisk_high"] = _extract_2020_flags(comp_silver, comp_high.get("compound_id"))
        # ADME 매핑 파일(선택적)
        adme_map = _read_csv_safe(os.path.join("data", "mappings", "2020_adme_links.csv"))
        if len(adme_map) == 0:
            adme_map = _read_csv_safe(os.path.join("base", "data", "mappings", "2020_adme_links.csv"))
        ctx["meta"]["adme"] = {"linked": False}
        try:
            if len(adme_map) > 0 and "compound_id" in adme_map.columns and "sample_id" in adme_map.columns:
                # meta에서 해당 sample_id의 등급을 찾아 반환
                def _adme_for(cid):
                    m = adme_map[adme_map["compound_id"].astype(str) == str(cid)]
                    if len(m) == 0:
                        return None
                    sid = _first_or_none(m.get("sample_id", []))
                    if sid is None:
                        return None
                    mm = gold_meta[gold_meta.get("sample_id", "").astype(str) == str(sid)] if "sample_id" in gold_meta.columns else pd.DataFrame()
                    if len(mm) == 0:
                        return None
                    rr = mm.iloc[0]
                    return {
                        "sample_id": sid,
                        "solubility": {"pH2": rr.get("solub_pH2_grade"), "pH7.4": rr.get("solub_pH7_grade")},
                        "stability": {"HLM_t12": rr.get("hlm_t12_grade"), "RLM_t12": rr.get("rlm_t12_grade")},
                    }
                low_adme = _adme_for(comp_low.get("compound_id"))
                high_adme = _adme_for(comp_high.get("compound_id"))
                if low_adme or high_adme:
                    ctx["meta"]["adme"] = {"linked": True, "low": low_adme, "high": high_adme}
        except Exception:
            pass

    if year == "2021":
        # 두 축 모두 요약을 제공
        gold_assay = _read_csv_safe(os.path.join("data", "gold", year, "assay_readings.csv"))
        if len(gold_assay) == 0:
            gold_assay = _read_csv_safe(os.path.join("base", "data", "gold", year, "assay_readings.csv"))
        ctx["meta"]["low_all_values"] = _extract_2021_allvalues(gold_assay, comp_low.get("compound_id"))
        ctx["meta"]["high_all_values"] = _extract_2021_allvalues(gold_assay, comp_high.get("compound_id"))

    return ctx


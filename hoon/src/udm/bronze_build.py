from __future__ import annotations

import os
import re
from typing import Dict, List, Optional, Tuple, Set

import pandas as pd

from .io_csv import read_csv_safe, write_csv_safe, ensure_parent_dir
from .parsers import split_censor_and_value, normalize_unit_label, parse_scientific_notation

# --- numeric token whitelist for measurement cells ---
NUMERIC_VALUE_RE = re.compile(r'^(?:\s*(?:>=|<=|[<>＝＝＜＞=])\s*)?[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?(?:\s*(?:μM|uM|µM|nM|mM|mg/kg))?\s*$')

def _is_numeric_value(text: object) -> bool:
    if text is None:
        return False
    s = str(text).strip()
    if s == "":
        return False
    return bool(NUMERIC_VALUE_RE.match(s))

# --- 2017 specific helpers (panel mapping, cell-line normalization) ---
# Panel code by table id for 2017 dataset
PANEL_ID_BY_TABLE_2017: Dict[str, Dict[str, str]] = {
    # table_id: {panel_id, panel_label, disease_area}
    "table5":  {"panel_id": "blca", "panel_label": "Bladder cancer panel",   "disease_area": "bladder"},
    "table6":  {"panel_id": "prad", "panel_label": "Prostate cancer panel",  "disease_area": "prostate"},
    "table7":  {"panel_id": "luad", "panel_label": "Lung cancer panel",      "disease_area": "lung"},
    "table8":  {"panel_id": "brca", "panel_label": "Breast cancer panel",    "disease_area": "breast"},
    "table9":  {"panel_id": "heme", "panel_label": "Hematologic panel",      "disease_area": "hematologic"},
    "table10": {"panel_id": "paad", "panel_label": "Pancreatic cancer panel","disease_area": "pancreas"},
    # Tables 11–13 are mixed/other in the source; keep generic labels
    "table11": {"panel_id": "gi",   "panel_label": "GI/Other panel",         "disease_area": "gi"},
    "table12": {"panel_id": "misc", "panel_label": "Misc tumor panel",       "disease_area": "misc"},
    "table13": {"panel_id": "misc", "panel_label": "Misc tumor panel",       "disease_area": "misc"},
    "table14": {"panel_id": "blca", "panel_label": "Bladder cancer panel",   "disease_area": "bladder"},
}

def _norm_cell_for_id(s: str) -> str:
    s = str(s or "").strip().lower()
    s = re.sub(r"[^a-z0-9]+", "-", s)
    return s.strip("-")

# --- compile YAML regex patterns (header_regex + capture_to_extras) ---
def _compile_yaml_patterns(cfg: Dict) -> List[Tuple[re.Pattern, Dict[str, object], Dict[str, str]]]:
    patterns: List[Tuple[re.Pattern, Dict[str, object], Dict[str, str]]] = []
    for p in cfg.get("parsing", {}).get("assay_patterns", []):
        header_regex = p.get("header_regex")
        assign = dict(p.get("assign", {}))
        capmap = dict(p.get("capture_to_extras", {}))
        if header_regex:
            try:
                patterns.append((re.compile(header_regex, re.IGNORECASE), assign, capmap))
            except re.error:
                continue
    return patterns


def _slugify(text: str) -> str:
	text = re.sub(r"[^0-9a-zA-Z]+", "-", text).strip("-").lower()
	return text or "unnamed"


def _is_plausible_smiles(text: Optional[str]) -> bool:
	"""보수적인 SMILES 타당성 검사를 통해 주석 문자열이 SMILES로 수집되는 것을 방지합니다.

	- 일반적인 SMILES 문자만 허용
	- 최소 하나의 원자 문자가 있어야 함
	"""
	if text is None:
		return False
	val = str(text).strip()
	if val == "":
		return False
	# common null markers should not be treated as SMILES
	if val.lower() in {"nan", "none", "null", "na", "n/a"}:
		return False
	# 공백/콤마 포함 문자열 배제 (문장형/화학명 제거)
	if (" " in val) or ("," in val):
		return False
	# 허용된 문자셋만 포함 (엄격)
	# 허용 집합: 원자기호/결합/고리/입체/숫자
	if not re.fullmatch(r"[A-Za-z0-9Hh@=+\-#\\/\\\\\[\]()\.]+", val):
		return False
	# 최소 하나의 원자 문자가 포함되어야 함 (방향족 소문자 포함)
	return any(ch in val for ch in "BCNOSPFIclonps")


def normalize_label(s: str) -> str:
	"""정규식 규칙과 일치시키기 위한 강력한 헤더 정규화.

	- 소문자로 변환, 공백 제거
	- 유니코드 대시를 '-'로 정규화
	- NBSP를 포함한 공백을 단일 공백으로 정규화
	- ASCII 괄호, en/em 대시, 마이크로 변형을 일치시키기 위해 교체
	- 일반적인 오타 통합 (icso/ic_so/lc50 -> ic50, hctl16/hctu6/hcti16 -> hct116)
	- 구두점 주변 공백 축소 (),-/
	- 여러 공백 축소
	"""
	s = (s or "").strip().lower()
	# 유니코드 대시 → '-'
	s = re.sub(r"[\u2010-\u2015]", "-", s)
	# 공백 정규화 (NBSP 포함)
	s = re.sub(r"[\u00a0\s]+", " ", s)
	# 밑줄을 구분자로 처리
	s = re.sub(r"_+", " ", s)
	# ASCII 괄호, 잘못된 쉼표 제거, 대시 정규화
	s = s.replace("（", "(").replace("）", ")").replace("–", "-").replace("—", "-")
	# 일치시키기 위해 마이크로 단위 토큰 정규화 (μM, μΜ, uM, μμ) → 'um'
	s = re.sub(r"[μu][mΜμ]", "um", s)
	# IC_50/IC 50/IC-50 → IC50으로 정규화
	s = re.sub(r"ic[\s_\-]*50", "ic50", s)
	# 일반적인 오타 수정 'prolil' → 'prolif'
	s = s.replace("prolil", "prolif")
	# 일반적인 오타 통합
	s = s.replace("icso", "ic50").replace("ic_so", "ic50").replace("ic_50", "ic50").replace("lc50", "ic50")
	# 공백이 있는 오타를 포함한 HCT116 변형
	s = re.sub(r"hct[l|i]\s*16", "hct116", s)
	s = s.replace("hctl16", "hct116").replace("hctu6", "hct116").replace("hcti16", "hct116")
	# 'prolif'와 'hct116' 사이 누락된 공백 삽입 보장
	s = re.sub(r"(prolif)hct116", r"\1 hct116", s)
	# 구두점 주변 공백 축소 (패턴 안전을 위해 괄호 간격 유지)
	s = re.sub(r"\s*([\-/,])\s*", r"\1", s)
	# 여러 공백 축소
	s = re.sub(r"\s{2,}", " ", s)
	return s


def _year_rules(dataset_id: str):
	"""주어진 데이터셋 연도에 대한 컴파일된 정규식 규칙과 후처리기를 반환합니다.

	각 규칙: (pattern, assign_dict, postprocess_func)
	- assign_dict는 target/readout/matrix와 unit(또는 unit_base) 및 선택적 extras를 제공
	- postprocess_func(match) -> (assay_id: str, extras_overrides: Dict)
	"""
	year = str(dataset_id)
	rules: List[Tuple[re.Pattern, Dict[str, object], object]] = []
	if year == "2021":
		# 2021 PRMT5 (표 19)
		# 효소 Ki (방법 A/B)
		pat_ki = re.compile(r"prmt5.*mtase\s*glo.*ki\(n?m,meth\.([ab])\)")
		assign_ki = {"target": "PRMT5", "readout": "Ki", "matrix": "enzyme", "unit": "nM"}
		def post_ki(m: re.Match) -> Tuple[str, Dict[str, str]]:
			method = m.group(1).upper()
			return f"prmt5_ki_meth{method}", {"method": method}
		# 세포 IC50 (HCT116 MTAP-null / WT)
		# 사양에 따른 세포 IC50 매핑
		pat_cell = re.compile(r"prolif\s+hct116[-–]?(mtap\s+null|wt)\s+ic50\s*\((?:um|μμ)\)")
		assign_cell = {"target": "PRMT5", "readout": "IC50", "matrix": "cell", "unit": "μM", "extras": {"cell_line": "HCT116"}}
		def post_cell(m: re.Match) -> Tuple[str, Dict[str, str]]:
			geno_cap = m.group(1).strip()
			geno_norm = "mtap null" if geno_cap.lower().startswith("mtap") else "wt"
			assay_suffix = "mtapnull" if geno_norm == "mtap null" else "wt"
			return f"hct116_{assay_suffix}_ic50", {"genotype": geno_norm}
		rules = [
			(pat_ki, assign_ki, post_ki),
			(pat_cell, assign_cell, post_cell),
		]
	elif year == "2018":
		# 2018 KRAS/EGFR
		pat_bind = re.compile(r"kras.*sos1.*interaction")
		assign_bind = {"target": "KRAS", "readout": "binding", "matrix": "biochemical", "unit_base": "M"}
		def post_bind(m: re.Match) -> Tuple[str, Dict[str, str]]:
			return "kras_sos1_interaction", {}
		pat_act_hi = re.compile(r"kras.*activation.*highgtp")
		assign_act_hi = {"target": "KRAS", "readout": "activation_highGTP", "matrix": "biochemical", "unit_base": "M"}
		def post_act_hi(m: re.Match) -> Tuple[str, Dict[str, str]]:
			return "kras_activation_highgtp", {}
		pat_act_no = re.compile(r"kras.*activation.*nogtp")
		assign_act_no = {"target": "KRAS", "readout": "activation_noGTP", "matrix": "biochemical", "unit_base": "M"}
		def post_act_no(m: re.Match) -> Tuple[str, Dict[str, str]]:
			return "kras_activation_nogtp", {}
		pat_egfr = re.compile(r"egfr.*kinase.*ic50")
		assign_egfr = {"target": "EGFR", "readout": "IC50", "matrix": "enzyme", "unit_base": "M", "right_censor_limit": 20e-6}
		def post_egfr(m: re.Match) -> Tuple[str, Dict[str, str]]:
			return "egfr_ic50", {}
		rules = [
			(pat_bind, assign_bind, post_bind),
			(pat_act_hi, assign_act_hi, post_act_hi),
			(pat_act_no, assign_act_no, post_act_no),
			(pat_egfr, assign_egfr, post_egfr),
		]
	return rules


def _detect_header_row(df: pd.DataFrame, max_scan: int = 30) -> int:
	"""진정한 열 헤더를 선호하는 휴리스틱 헤더 탐지기.

	첫 max_scan 행 내 선호 순서:
	1) 'compound' 헤더 토큰 (예: 'compound #', 'compound no.')과 'smiles' 헤더 토큰을 모두 포함하는 행
	2) {'smiles','compound','num'} 중 하나를 포함하는 첫 행으로 폴백
	"""
	limit = min(max_scan, len(df))
	# 1단계: 엄격한 헤더 토큰
	for i in range(limit):
		vals = [str(x).strip() for x in df.iloc[i].tolist() if str(x).strip() != ""]
		vals_l = [v.lower() for v in vals]
		def is_compound_token(v: str) -> bool:
			return bool(re.fullmatch(r"compound(\s*#|\s*no\.?|)?", v)) or v in {"compound", "compound #", "compound no", "compound no."}
		def is_smiles_token(v: str) -> bool:
			# 정확한 'smiles' 또는 구분자 뒤에 'smiles'로 끝남
			return v == "smiles" or v.endswith(" smiles") or v.endswith("- smiles")
		if any(is_compound_token(v) for v in vals_l) and any(is_smiles_token(v) for v in vals_l):
			return i
	# 2단계: 레거시 느슨한 휴리스틱
	for i in range(limit):
		row = " ".join(str(x) for x in df.iloc[i].tolist()).lower()
		if "smiles" in row or "compound" in row or "num" in row:
			return i
	return 0


def _find_col(cols: List[str], keywords: List[str]) -> Optional[str]:
	cols_l = [c.lower() for c in cols]
	for kw in keywords:
		for c, lc in zip(cols, cols_l):
			if kw in lc:
				return c
	return None


def build_bronze_from_raw(root_dir: str, cfg: Dict) -> Dict[str, str]:
	"""
	일반 빌더: 최소 브론즈 아티팩트를 채우기 위해 원본 Excel을 스캔합니다.
	- compounds.csv: compound_id, smiles_raw + provenance + dataset_id
	- assays.csv: 탐지된 숫자 어세이 열당 한 행 (레이블 그대로)
	- measurements.csv: 어세이 열에서 녹인 긴 형식; 검열 분리; value_raw 보존
	- panel_block_meta.csv: 표(panel) 단위 메타데이터 집계(table_id, cell_lines, max_test_conc)
	"""
	outputs: Dict[str, str] = {}
	raw_path = os.path.join(root_dir, cfg["paths"]["raw_file"])
	bronze_dir = os.path.join(root_dir, cfg["paths"]["bronze_dir"])
	os.makedirs(bronze_dir, exist_ok=True)

	# 홀더 준비
	# compounds는 중복키를 방지하기 위해 dict로 수집 (키: 정수 compound_id)
	comp_rows: Dict[int, Dict[str, str]] = {}
	assay_rows: Dict[str, Dict[str, str]] = {}
	meas_rows: List[Dict[str, str]] = []
	# panel 메타 집계 (table_id 단위)
	panel_acc: Dict[str, Dict[str, object]] = {}
	# 원본 표의 열 순서 보존을 위한 (sheet_name, original_column_label) → order 인덱스 맵
	sheet_col_order: Dict[Tuple[str, str], int] = {}

	# 2017: 패널 메타(패널 단위) 누적용
	panel_meta_acc: Dict[str, Dict[str, object]] = {}

	# YAML에서 어세이 패턴 맵 준비 (선택사항, 정확한 헤더 매치 폴백)
	patterns = cfg.get("parsing", {}).get("assay_patterns", [])
	col_to_assign: Dict[str, Dict[str, str]] = {}
	norm_to_assign: Dict[str, Dict[str, str]] = {}
	for p in patterns:
		assign_map = p.get("assign", {})
		for col in p.get("columns", []):
			col_s = str(col)
			col_to_assign[col_s] = assign_map
			# 정규화된 라벨로도 매칭 가능하도록 보조 맵 구축
			norm_label_key = normalize_label(col_s)
			norm_to_assign[norm_label_key] = assign_map

	# 데이터셋 연도를 기반으로 정규식 규칙 컴파일
	regex_rules = _year_rules(cfg.get("dataset_id", ""))
	# YAML header_regex 패턴 컴파일
	regex_from_yaml = _compile_yaml_patterns(cfg)
	# 패널명 매핑(선택): panel_block_meta에 기록
	panel_name_map: Dict[str, str] = cfg.get("parsing", {}).get("panel_names", {})

	# Excel 로드
	xl = pd.ExcelFile(raw_path)
	# config에서 compounds를 수집할 시트를 선택적으로 제한
	allowed_compound_sheets = cfg.get("parsing", {}).get("compounds_from_sheets", None)
	# 2017: assay 시트 화이트리스트 정규식 (표5~14만 인제스트)
	assay_whitelist_re: Optional[str] = cfg.get("parsing", {}).get("assay_sheets_whitelist_regex")
	# 2017: 세포주 동의어 병합 맵
	cell_syn = cfg.get("parsing", {}).get("cell_synonyms", {})
	cell_map_yaml = cfg.get("parsing", {}).get("cell_line_map", {})
	merged_cell_map: Dict[str, str] = {**cell_map_yaml, **cell_syn}
	# 2017: 패널 정의(YAML 우선)
	panel_defs_yaml: Dict[str, Dict[str, object]] = cfg.get("parsing", {}).get("panel_definitions", {})
	for sheet_name in xl.sheet_names:
		is_2017 = str(cfg.get("dataset_id", "")) == "2017"
		assay_ingest_allowed = True
		if is_2017 and assay_whitelist_re:
			try:
				assay_ingest_allowed = bool(re.search(assay_whitelist_re, sheet_name, flags=re.IGNORECASE))
			except re.error:
				assay_ingest_allowed = True
		# 1차 통과: 헤더 행 탐지
		head_df = xl.parse(sheet_name, header=None, nrows=120)
		header_row = _detect_header_row(head_df, max_scan=100)
		df = xl.parse(sheet_name, header=header_row)
		# 열 정규화 및 기준 정규화 헤더 맵 구축 (단일 행 헤더)
		df.columns = [str(c).strip() for c in df.columns]
		norm_map: Dict[str, str] = {orig: normalize_label(str(orig)) for orig in df.columns}

		def _to_flat_columns(df_multi: pd.DataFrame) -> List[str]:
			cols: List[str] = []
			for tup in df_multi.columns.to_flat_index():
				parts = [str(x).strip() for x in tup if str(x).strip() != "" and str(x).lower() != "nan"]
				cols.append(" ".join(parts))
			# 2021 표 19에서 흔히 나타나는 노이즈 있는 선행 접두사 연결 및 제거
			flat = [re.sub(r"\s+", " ", c).strip() for c in cols]
			flat = [re.sub(r"^red\s+highlight\s+is\s+just\s+for\s+page\s+reference\s*", "", c) for c in flat]
			return flat

		# 2행 헤더 시도 (header_row, header_row+1) 및 더 나은 매치 점수 선택
		best_score = 0
		# _match_score가 정의된 후 나중에 base_score를 계산하지만, 순서를 유지하기 위해 지금 근사
		# _match_score 정의 후 정확한 점수 매기기 연기; 임시로 후보 저장
		single_df, single_norm = df, norm_map
		try:
			df2_try = xl.parse(sheet_name, header=[header_row, header_row + 1])
			df2_try.columns = _to_flat_columns(df2_try)
			n2_try = {orig: normalize_label(str(orig)) for orig in df2_try.columns}
			double_candidate = (df2_try, n2_try)
		except Exception:
			double_candidate = None

		# 폴백: 최대 정규식 매치로 대안 헤더 행 검색
		def _match_score(cols_map: Dict[str, str]) -> int:
			"""YAML 정확 매치 또는 정규식 규칙과 일치하는 열 수를 계산합니다."""
			score = 0
			for orig, norm in cols_map.items():
				if orig in col_to_assign:
					score += 1
					continue
				for pat, _, __ in regex_rules:
					if pat.search(norm):
						score += 1
						break
			return score

		# 이제 단일 행 vs 2행 변형 점수 매기기 (가능한 경우)
		base_score = _match_score(norm_map)
		best = (base_score, header_row, df, norm_map)
		if double_candidate is not None:
			s2 = _match_score(double_candidate[1])
			# 분할 헤더를 캡처하기 위해 동점 시 2행 변형 선호
			if s2 >= best[0]:
				best = (s2, header_row, double_candidate[0], double_candidate[1])

		# 휴리스틱: 단일 행 헤더에는 'compound'와 'smiles/structure'가 모두 없지만,
		# 2행 헤더에는 둘 다 존재하면 2행 헤더를 채택
		def _has_comp_and_smiles(cols: List[str]) -> bool:
			lc = [str(c).strip().lower() for c in cols]
			comp_kws = [
				"compound #", "compound no", "compound", "cmpd", "entry", "id", "num",
				"화합물 번호", "화합물", "번호",
			]
			smiles_kws = ["smiles", "structure", "구조"]
			comp_ok = any(any(kw in c for kw in comp_kws) for c in lc)
			smi_ok = any(any(kw in c for kw in smiles_kws) for c in lc)
			return comp_ok and smi_ok

		single_has = _has_comp_and_smiles(list(df.columns))
		double_has = False
		if double_candidate is not None:
			double_has = _has_comp_and_smiles(list(double_candidate[0].columns))
		if double_candidate is not None and double_has and not single_has:
			df = double_candidate[0]
			norm_map = double_candidate[1]
		for cand in range(0, min(120, len(head_df))):
			try:
				df2 = xl.parse(sheet_name, header=cand)
				df2.columns = [str(c).strip() for c in df2.columns]
				n2 = {orig: normalize_label(str(orig)) for orig in df2.columns}
				s = _match_score(n2)
				if s >= best[0]:
					best = (s, cand, df2, n2)
				# 이 후보에 대해 2행 변형 시도
				try:
					df2b = xl.parse(sheet_name, header=[cand, cand + 1])
					df2b.columns = _to_flat_columns(df2b)
					n2b = {orig: normalize_label(str(orig)) for orig in df2b.columns}
					sb = _match_score(n2b)
					# 동점 또는 개선 시 2행 선호
					if sb >= s and sb >= best[0]:
						best = (sb, cand, df2b, n2b)
				except Exception:
					pass
			except Exception:
				continue
		# 개선된 경우 최고 점수 헤더 채택
		if best[0] > base_score and best[0] > 0:
			df = best[2]
			norm_map = best[3]
		# 2021 디버그 헤더 로깅: PRMT5/세포 IC50과 관련된 연결된 헤더 후보 출력
		if str(cfg.get("dataset_id", "")) == "2021":
			for orig, norm in norm_map.items():
				if any(k in norm for k in ["prolif", "hct", "ic50", "prmt5"]):
					print(f"[HDR2021] sheet={sheet_name} orig={orig} norm={norm}")

		# 여러 헤더 변형 후보를 수집하여(2017 복수 테이블 대응) 순차 처리
		header_variants: List[Tuple[pd.DataFrame, Dict[str, str]]] = []
		header_variants.append((df, norm_map))
		seen_fps: Set[str] = set()
		seen_fps.add("|".join(sorted([f"{k}:{v}" for k, v in norm_map.items()])))
		# scan additional header rows to catch other tables in the same sheet
		limit_scan = min(120, len(head_df))
		for cand in range(0, limit_scan):
			try:
				df2 = xl.parse(sheet_name, header=cand)
				df2.columns = [str(c).strip() for c in df2.columns]
				n2 = {orig: normalize_label(str(orig)) for orig in df2.columns}
				s = _match_score(n2)
				fp = "|".join(sorted([f"{k}:{v}" for k, v in n2.items()]))
				if s > 0 and fp not in seen_fps:
					header_variants.append((df2, n2))
					seen_fps.add(fp)
					# also try two-row header for this candidate
					try:
						df2b = xl.parse(sheet_name, header=[cand, cand + 1])
						df2b.columns = _to_flat_columns(df2b)
						n2b = {orig: normalize_label(str(orig)) for orig in df2b.columns}
						fpb = "|".join(sorted([f"{k}:{v}" for k, v in n2b.items()]))
						sb = _match_score(n2b)
						if sb > 0 and fpb not in seen_fps:
							header_variants.append((df2b, n2b))
							seen_fps.add(fpb)
					except Exception:
						pass
			except Exception:
				continue

		for cur_df, cur_norm in header_variants:
			# 핵심 열 식별
			compound_col = _find_col(
				list(cur_df.columns),
				[
					"compound #",
					"compound no",
					"compound",
					"cmpd",
					"num",
					"no.",
					"no",
					"example",
					"ex.",
					"ex ",
					"entry",
					"id",
					"화합물 번호",
					"화합물",
					"번호",
				]
			)
			smiles_col = _find_col(list(cur_df.columns), ["smiles", "structure", "구조"])
			# 콘텐츠 기반 보정: 열 이름이 애매할 때 상단 몇 행을 검사하여 SMILES/번호를 추정
			if not compound_col or not smiles_col:
				try:
					cur_head = cur_df.head(50)
					for c in cur_df.columns:
						lc = str(c).strip().lower()
						col_str = cur_head[c].astype(str)
						is_smiles_like = col_str.str.fullmatch(r"[A-Za-z0-9@=+\-#\\/\\\\\[\]()\.]+").fillna(False).any() and \
							col_str.str.contains(r"[BCNOSPFIclonps]").fillna(False).any()
						is_compound_like = col_str.str.fullmatch(r"\s*\d+\s*").fillna(False).any() or ("compound" in lc)
						if not smiles_col and ("smiles" in lc or "structure" in lc or "구조" in lc or is_smiles_like):
							smiles_col = c
						if not compound_col and ("compound" in lc or "화합물" in lc or "번호" in lc or is_compound_like):
							compound_col = c
				except Exception:
					pass

			# compounds 수집: 허용된 시트가 명시된 경우 해당 시트에서만 수집
			collect_compounds_here = compound_col and smiles_col and (
				(allowed_compound_sheets is None) or (sheet_name in allowed_compound_sheets)
			)
			if collect_compounds_here:
				for idx, row in cur_df[[compound_col, smiles_col]].dropna(how="all").iterrows():
					compound_id_raw = str(row[compound_col]).strip()
					smiles_raw = str(row[smiles_col]).strip()
					if compound_id_raw == "" and smiles_raw == "":
						continue
					try:
						cid = int(compound_id_raw)
					except Exception:
						continue
					if not _is_plausible_smiles(smiles_raw):
						continue
					row_rec: Dict[str, str] = {
						"compound_id": str(cid),
						"smiles_raw": smiles_raw,
						"dataset_id": str(cfg.get("dataset_id", "")),
						"provenance_file": cfg["paths"]["raw_file"],
						"provenance_sheet": sheet_name,
						"provenance_row": str(int(idx) + 1),
					}
					if cid not in comp_rows:
						comp_rows[cid] = row_rec
					# else: keep first

			# 타당한 어세이 열 탐지: 숫자형 또는 assay_patterns에 나열된 열
			numeric_like_cols: List[str] = []
			if assay_ingest_allowed:
				for c in cur_df.columns:
					if c == compound_col or c == smiles_col:
						continue
					series = cur_df[c].astype(str)
					norm_c = cur_norm.get(c, normalize_label(str(c)))
					is_pattern_listed = (c in col_to_assign) or (norm_c in norm_to_assign)
					is_numericish = series.str.contains(r"^[\s]*(?:(?:>=|<=|[<>＝=＜＞])?)[\s]*(?:[0-9.]+(?:[eE][+-]?[0-9]+)?)", regex=True).any() or series.str.contains("%", regex=False).any()
					if is_pattern_listed or is_numericish:
						numeric_like_cols.append(c)
			# 현재 시트에서 감지된 측정 컬럼의 좌→우 순서 기록
			for ord_idx, col_name in enumerate(numeric_like_cols):
				sheet_col_order[(sheet_name, str(col_name))] = ord_idx
			# 어세이 등록
			for col in numeric_like_cols:
				# 방어적: 현재 변형에 해당 열이 없으면 스킵
				if str(col) not in [str(x) for x in cur_df.columns]:
					continue
				orig_label = str(col)
				norm_label = cur_norm.get(orig_label, normalize_label(orig_label))
				assign = col_to_assign.get(orig_label, {})  # YAML exact-match fallback
				if not assign:
					assign = norm_to_assign.get(norm_label, {})
				is_pattern_listed_here = (orig_label in col_to_assign) or (norm_label in norm_to_assign)
				assay_id: str
				post_extras: Dict[str, str] = {}
				matched = False

				# 1) try year-specific regex rules first
				for pat, assign_dict, post_fn in regex_rules:
					m = pat.search(norm_label)
					if m:
						matched = True
						assign = {**assign_dict}
						assay_id, extra = post_fn(m)
						post_extras = extra or {}
						break

				# 2) if not matched, try YAML header_regex patterns
				if not matched:
					for pat, assign_dict, capmap in regex_from_yaml:
						m = pat.search(norm_label)
						if m:
							matched = True
							assign = {**assign_dict}
							assay_id = assign.get("assay_id") or f"assay.unknown.{_slugify(orig_label)}"
							extras = dict(assign.get("extras", {}))
							for grp, key in capmap.items():
								try:
									extras[key] = m.group(int(grp))
								except Exception:
									pass
							clmap = cfg.get("parsing", {}).get("cell_line_map", {})
							if "cell_line" in extras and extras["cell_line"] in clmap:
								extras["cell_line"] = clmap[extras["cell_line"]]
							assign["extras"] = extras
							break

				# 3) if still not matched but exact columns listed, accept it
				if not matched and is_pattern_listed_here:
					matched = True
					assay_id = assign.get("assay_id") or f"assay.unknown.{_slugify(orig_label)}"

				# 4) if still not matched, skip this column entirely
				if not matched:
					continue

				# 매핑 필드 읽기
				readout = str(assign.get("readout", str(orig_label)))
				matrix = str(assign.get("matrix", ""))
				unit_default = assign.get("unit")  # may be "%" or base unit like "M"/"μM"/"nM"
				unit_base = assign.get("unit_base", unit_default)
				extras = dict(assign.get("extras", {}))
				extras.update(post_extras)
				# 2017: 패널/세포주 구분이 가능한 assay_id로 확장 및 패널 정보 주입
				if is_2017:
					# 세포주 동의어 정규화
					if "cell_line" in extras and extras["cell_line"] in merged_cell_map:
						extras["cell_line"] = merged_cell_map[extras["cell_line"]]
					table_id_here = str(extras.get("table_id", "")).strip().lower()
					panel_info = None
					if extras.get("panel_id"):
						panel_info = {
							"panel_id": str(extras.get("panel_id")),
							"panel_label": str(extras.get("panel_label", "")),
							"disease_area": str(extras.get("disease_area", "")),
						}
					elif table_id_here in panel_defs_yaml:
						panel_def = panel_defs_yaml.get(table_id_here, {})
						panel_info = {
							"panel_id": str(panel_def.get("panel_id", "")),
							"panel_label": str(panel_def.get("panel_label", "")),
							"disease_area": str(panel_def.get("disease_area", "")),
						}
					elif table_id_here in PANEL_ID_BY_TABLE_2017:
						panel_info = PANEL_ID_BY_TABLE_2017[table_id_here]
					if panel_info:
						cl_for_id = _norm_cell_for_id(extras.get("cell_line", ""))
						assay_id = f"assay.cell.cytotoxicity.ic50.um.{panel_info['panel_id']}.{cl_for_id}"
						extras.setdefault("panel_id", panel_info["panel_id"])
						extras.setdefault("panel_label", panel_info.get("panel_label", ""))
						extras.setdefault("disease_area", panel_info.get("disease_area", ""))
					# IC50로 강제
					assign["meas_kind_numeric"] = "IC50"
				else:
					# 2017 이외: 동일 assay 라벨이 여러 패널에 나타나면 table_id suffix로 구분
					if "table_id" in extras and isinstance(assay_id, str) and len(assay_id) > 0:
						tbl = str(extras.get("table_id", "")).strip()
						if tbl and not assay_id.endswith(f".{tbl}"):
							assay_id = f"{assay_id}.{tbl}"
				# 어세이 메타데이터 등록 (추적성을 위한 선택적 target/unit 필드 포함)
				if assay_id not in assay_rows:
					row_assay: Dict[str, str] = {
						"assay_id": assay_id,
						"readout": readout,
						"matrix": matrix,
						"label_raw": orig_label,
					}
					if "target" in assign:
						row_assay["target"] = str(assign.get("target", ""))
					if unit_default:
						row_assay["unit"] = str(unit_default)
					if unit_base and unit_base != unit_default:
						row_assay["unit_base"] = str(unit_base)
					# 2017 확장 필드
					if is_2017:
						row_assay["panel_id"] = str(extras.get("panel_id", ""))
						row_assay["cell_line"] = str(extras.get("cell_line", ""))
						row_assay["disease_area"] = str(extras.get("disease_area", ""))
						row_assay["meas_kind"] = "IC50"
					assay_rows[assay_id] = row_assay
					# 측정값으로 녹이기
					if compound_col and (compound_col in cur_df.columns):
						try:
							iter_frame = cur_df[[compound_col, col]].dropna(how="all")
						except Exception:
							continue
						for idx, cell in iter_frame.iterrows():
								compound_id = str(cell[compound_col]).strip()
								value_raw = str(cell[col]).strip()
								if compound_id == "" and value_raw == "":
									continue
								# 2017: 빈 compound_id를 만나면 해당 열 파싱 종료(footer 방지)
								if is_2017 and compound_id == "":
									break
								# 2017 table14는 자유 텍스트 허용, 그 외는 숫자 강제
								is_table14_here = is_2017 and str(extras.get("table_id", "")).strip().lower() == "table14"
								if not is_table14_here:
									try:
										_cid_int = int(float(compound_id))
									except Exception:
										continue
								# filter out header/label strings like RT4, LoVo, 암종명 등
								if not _is_numeric_value(value_raw):
									continue
								# use normalized numeric id (table14는 숫자 추출 가능 시만)
								if not is_table14_here:
									compound_id = str(_cid_int)
								else:
									m = re.search(r"(\d+)", compound_id)
									if m:
										compound_id = str(int(m.group(1)))
								# 값 모드: 백분율 vs 숫자
								is_percent = "%" in value_raw
								censor, number_text, unit_tok = split_censor_and_value(value_raw)
								unit_from_cell = normalize_unit_label(unit_tok) if unit_tok else None
								unit_final = "%" if is_percent else (unit_from_cell or unit_base or unit_default or "")
								value_num: Optional[float] = None
								if number_text is not None and not is_percent:
									val, _ = parse_scientific_notation(number_text)
									value_num = val
								row_meas: Dict[str, object] = {
									"compound_id": compound_id,
									"assay_id": assay_id,
									"value_raw": value_raw,
									"value_num": value_num if value_num is not None else "",
									"unit": unit_final,
									"censor": censor or "",
									"assay_label": orig_label,
									"provenance_file": cfg["paths"]["raw_file"],
									"provenance_sheet": sheet_name,
									"provenance_row": str(int(idx) + 1),
								}
								if is_percent:
									row_meas["meas_kind"] = "percent_inhibition_at_20uM"
								else:
									row_meas["meas_kind"] = assign.get("meas_kind_numeric", "IC50")
								if censor == "gt" and assign.get("right_censor_limit"):
									row_meas["limit"] = str(assign.get("right_censor_limit"))
								for k, v in extras.items():
									row_meas[k] = v
								if is_2017 and extras.get("panel_id"):
									row_meas["panel_id"] = extras.get("panel_id")
								if is_table14_here and re.search(r"[A-Za-z가-힣]", str(cell[compound_col])):
									row_meas["compound_name_raw"] = str(cell[compound_col]).strip()
								if "table_id" not in row_meas:
									row_meas["table_id"] = extras.get("table_id", "")
								tbl_here = str(extras.get("table_id", "")).strip()
								if tbl_here:
									acc = panel_acc.setdefault(tbl_here, {"cell_lines": set(), "assay_label": "", "thresholds": []})
									cl = str(extras.get("cell_line", "")).strip()
									if cl:
										acc["cell_lines"].add(cl)
									assay_lbl = (readout.capitalize() + (f" ({unit_base or unit_default})" if (unit_base or unit_default) else "")).strip()
									if assay_lbl and not acc["assay_label"]:
										acc["assay_label"] = assay_lbl
									if censor == "gt" and (value_num is not None):
										try:
											acc["thresholds"].append(float(value_num))
										except Exception:
											pass
								# 2017: 패널 메타(패널 id 단위) 누적
								if is_2017 and extras.get("panel_id"):
									pid = str(extras.get("panel_id"))
									pmeta = panel_meta_acc.setdefault(pid, {
										"panel_label": str(extras.get("panel_label", "")),
										"disease_area": str(extras.get("disease_area", "")),
										"provenance_sheet": sheet_name,
										"cell_lines": set(),
									})
									try:
										cast_set: Set[str] = pmeta["cell_lines"]  # type: ignore
										cast_set.add(str(extras.get("cell_line", "")))
									except Exception:
										pass
							meas_rows.append(row_meas)
				continue
				orig_label = str(col)
			norm_label = norm_map.get(orig_label, normalize_label(orig_label))
			assign = col_to_assign.get(orig_label, {})  # YAML exact-match fallback
			if not assign:
				assign = norm_to_assign.get(norm_label, {})
			is_pattern_listed_here = (orig_label in col_to_assign) or (norm_label in norm_to_assign)
			assay_id: str
			post_extras: Dict[str, str] = {}
			matched = False

			# 1) try year-specific regex rules first
			for pat, assign_dict, post_fn in regex_rules:
				m = pat.search(norm_label)
				if m:
					matched = True
					assign = {**assign_dict}
					assay_id, extra = post_fn(m)
					post_extras = extra or {}
					break

			# 2) if not matched, try YAML header_regex patterns
			if not matched:
				for pat, assign_dict, capmap in regex_from_yaml:
					m = pat.search(norm_label)
					if m:
						matched = True
						assign = {**assign_dict}
						assay_id = assign.get("assay_id") or f"assay.unknown.{_slugify(orig_label)}"
						extras = dict(assign.get("extras", {}))
						for grp, key in capmap.items():
							try:
								extras[key] = m.group(int(grp))
							except Exception:
								pass
						clmap = cfg.get("parsing", {}).get("cell_line_map", {})
						if "cell_line" in extras and extras["cell_line"] in clmap:
							extras["cell_line"] = clmap[extras["cell_line"]]
						assign["extras"] = extras
						break

			# 3) if still not matched but exact columns listed, accept it
			if not matched and is_pattern_listed_here:
				matched = True
				assay_id = assign.get("assay_id") or f"assay.unknown.{_slugify(orig_label)}"

			# 4) if still not matched, skip this column entirely
			if not matched:
				continue

			# 매핑 필드 읽기
			readout = str(assign.get("readout", str(orig_label)))
			matrix = str(assign.get("matrix", ""))
			unit_default = assign.get("unit")  # may be "%" or base unit like "M"/"μM"/"nM"
			unit_base = assign.get("unit_base", unit_default)
			extras = dict(assign.get("extras", {}))
			extras.update(post_extras)
			# 2017: 패널/세포주 구분이 가능한 assay_id로 확장 및 패널 정보 주입
			if is_2017:
				# 세포주 동의어 정규화
				if "cell_line" in extras and extras["cell_line"] in merged_cell_map:
					extras["cell_line"] = merged_cell_map[extras["cell_line"]]
				table_id_here = str(extras.get("table_id", "")).strip().lower()
				panel_info = None
				if extras.get("panel_id"):
					panel_info = {
						"panel_id": str(extras.get("panel_id")),
						"panel_label": str(extras.get("panel_label", "")),
						"disease_area": str(extras.get("disease_area", "")),
					}
				elif table_id_here in panel_defs_yaml:
					panel_def = panel_defs_yaml.get(table_id_here, {})
					panel_info = {
						"panel_id": str(panel_def.get("panel_id", "")),
						"panel_label": str(panel_def.get("panel_label", "")),
						"disease_area": str(panel_def.get("disease_area", "")),
					}
				elif table_id_here in PANEL_ID_BY_TABLE_2017:
					panel_info = PANEL_ID_BY_TABLE_2017[table_id_here]
				if panel_info:
					cl_for_id = _norm_cell_for_id(extras.get("cell_line", ""))
					assay_id = f"assay.cell.cytotoxicity.ic50.um.{panel_info['panel_id']}.{cl_for_id}"
					extras.setdefault("panel_id", panel_info["panel_id"])
					extras.setdefault("panel_label", panel_info.get("panel_label", ""))
					extras.setdefault("disease_area", panel_info.get("disease_area", ""))
				# IC50로 강제
				assign["meas_kind_numeric"] = "IC50"
			else:
				# 2017 이외: 동일 assay 라벨이 여러 패널에 나타나면 table_id suffix로 구분
				if "table_id" in extras and isinstance(assay_id, str) and len(assay_id) > 0:
					tbl = str(extras.get("table_id", "")).strip()
					if tbl and not assay_id.endswith(f".{tbl}"):
						assay_id = f"{assay_id}.{tbl}"
			# 어세이 메타데이터 등록 (추적성을 위한 선택적 target/unit 필드 포함)
			if assay_id not in assay_rows:
				row_assay: Dict[str, str] = {
					"assay_id": assay_id,
					"readout": readout,
					"matrix": matrix,
					"label_raw": orig_label,
				}
				if "target" in assign:
					row_assay["target"] = str(assign.get("target", ""))
				if unit_default:
					row_assay["unit"] = str(unit_default)
				if unit_base and unit_base != unit_default:
					row_assay["unit_base"] = str(unit_base)
				# 2017: assay 메타에 panel/cell_line/disease/meas_kind 추가
				if is_2017:
					row_assay["panel_id"] = str(extras.get("panel_id", ""))
					row_assay["cell_line"] = str(extras.get("cell_line", ""))
					row_assay["disease_area"] = str(extras.get("disease_area", ""))
					row_assay["meas_kind"] = "IC50"
				assay_rows[assay_id] = row_assay
				# 측정값으로 녹이기
					if compound_col and (compound_col in cur_df.columns):
						try:
							iter_frame = cur_df[[compound_col, col]].dropna(how="all")
						except Exception:
							continue
						for idx, cell in iter_frame.iterrows():
						compound_id = str(cell[compound_col]).strip()
						value_raw = str(cell[col]).strip()
						if compound_id == "" and value_raw == "":
							continue
						# 2017: 빈 compound_id를 만나면 해당 열 파싱 종료(footer 방지)
						if is_2017 and compound_id == "":
							break
						# 2017 table14는 자유 텍스트 허용, 그 외는 숫자 강제
						is_table14_here = is_2017 and str(extras.get("table_id", "")).strip().lower() == "table14"
						if not is_table14_here:
							# enforce numeric compound ids to avoid header leakage
							try:
								_cid_int = int(float(compound_id))
							except Exception:
								continue
						# filter out header/label strings like RT4, LoVo, 암종명 등
						if not _is_numeric_value(value_raw):
							continue
						# use normalized numeric id (table14는 숫자 추출 가능 시만)
						if not is_table14_here:
							compound_id = str(_cid_int)
						else:
							m = re.search(r"(\d+)", compound_id)
							if m:
								compound_id = str(int(m.group(1)))
						# 값 모드: 백분율 vs 숫자
						is_percent = "%" in value_raw
						censor, number_text, unit_tok = split_censor_and_value(value_raw)
						# 셀에 있는 텍스트 단위 선호; 그렇지 않으면 매핑 단위(또는 기본값)로 폴백
						unit_from_cell = normalize_unit_label(unit_tok) if unit_tok else None
						unit_final = "%" if is_percent else (unit_from_cell or unit_base or unit_default or "")
						# 사용 가능한 경우 숫자 문자열을 value_num의 float로 변환
						value_num: Optional[float] = None
						if number_text is not None and not is_percent:
							val, _ = parse_scientific_notation(number_text)
							value_num = val
						row_meas: Dict[str, object] = {
							"compound_id": compound_id,
							"assay_id": assay_id,
							"value_raw": value_raw,
						"value_num": value_num if value_num is not None else "",
						"unit": unit_final,
						"censor": censor or "",
						"assay_label": orig_label,
						"provenance_file": cfg["paths"]["raw_file"],
						"provenance_sheet": sheet_name,
						"provenance_row": str(int(idx) + 1),
					}
					# 특별 규칙에 따라 meas_kind 및 제한 추가
					if is_percent:
						# 작업에서 명시적으로 요청된 문자열
						row_meas["meas_kind"] = "percent_inhibition_at_20uM"
					else:
						row_meas["meas_kind"] = assign.get("meas_kind_numeric", "IC50")
					if censor == "gt" and assign.get("right_censor_limit"):
						row_meas["limit"] = str(assign.get("right_censor_limit"))
					# 문맥적 추가 정보 추가 (예: cell_line, genotype)
					for k, v in extras.items():
						row_meas[k] = v
					# 2017: panel_id/meas_kind 추가, table14 자유 텍스트 보존
					if is_2017 and extras.get("panel_id"):
						row_meas["panel_id"] = extras.get("panel_id")
					if is_table14_here and re.search(r"[A-Za-z가-힣]", str(cell[compound_col])):
						row_meas["compound_name_raw"] = str(cell[compound_col]).strip()
					# table_id 열 보존 (패널 식별)
					if "table_id" not in row_meas:
						row_meas["table_id"] = extras.get("table_id", "")
					# panel 메타 집계: table_id가 있는 경우 cell_line 수집 및 우측 검출 상한 수집
					tbl_here = str(extras.get("table_id", "")).strip()
					if tbl_here:
						acc = panel_acc.setdefault(tbl_here, {"cell_lines": set(), "assay_label": "", "thresholds": []})
						cl = str(extras.get("cell_line", "")).strip()
						if cl:
							acc["cell_lines"].add(cl)
						# panel 공통 assay 라벨: Readout(Unit) 형식으로 보수적으로 구성
						assay_lbl = (readout.capitalize() + (f" ({unit_base or unit_default})" if (unit_base or unit_default) else "")).strip()
						if assay_lbl and not acc["assay_label"]:
							acc["assay_label"] = assay_lbl
						if censor == "gt" and (value_num is not None):
							try:
								acc["thresholds"].append(float(value_num))
							except Exception:
								pass
						# 2017: 패널 메타(패널 id 단위) 누적
						if is_2017 and extras.get("panel_id"):
							pid = str(extras.get("panel_id"))
							pmeta = panel_meta_acc.setdefault(pid, {
								"panel_label": str(extras.get("panel_label", "")),
								"disease_area": str(extras.get("disease_area", "")),
								"provenance_sheet": sheet_name,
								"cell_lines": set(),
							})
							try:
								cast_set: Set[str] = pmeta["cell_lines"]  # type: ignore
								cast_set.add(str(extras.get("cell_line", "")))
							except Exception:
								pass
					meas_rows.append(row_meas)

	# 출력 작성
	comp_cols = ["compound_id", "smiles_raw", "dataset_id", "provenance_file", "provenance_sheet", "provenance_row"]
	assay_cols = ["assay_id", "readout", "matrix", "label_raw", "target", "unit", "unit_base",
				 "panel_id", "cell_line", "disease_area", "meas_kind"]
	meas_cols = [
		"compound_id", "assay_id", "panel_id", "cell_line", "genotype", "value_raw", "value_num", "unit", "censor",
		"assay_label", "provenance_file", "provenance_sheet", "provenance_row", "meas_kind", "limit",
		"table_id", "max_test_conc", "compound_name_raw",
	]

	# dict → dataframe (compound_id 오름차순)
	comp_df = pd.DataFrame(sorted(comp_rows.values(), key=lambda r: int(r["compound_id"])), columns=comp_cols).drop_duplicates()
	assay_df = pd.DataFrame(list(assay_rows.values()), columns=assay_cols).drop_duplicates()
	meas_df = pd.DataFrame(meas_rows, columns=meas_cols)
	# 정렬: 같은 시트 내에서 원본 행 순서(provenance_row)와 열 순서(sheet_col_order)를 보존
	if len(meas_df) > 0:
		meas_df["__prov_row_num"] = pd.to_numeric(meas_df["provenance_row"], errors="coerce")
		def _col_ord(r):
			key = (str(r["provenance_sheet"]), str(r["assay_label"]))
			return sheet_col_order.get(key, 10**9)
		meas_df["__col_ord"] = meas_df.apply(_col_ord, axis=1)
		meas_df = meas_df.sort_values(["provenance_sheet", "__prov_row_num", "__col_ord"])\
			.reset_index(drop=True)
		meas_df = meas_df.drop(columns=["__prov_row_num", "__col_ord"], errors="ignore")

	comp_path = os.path.join(bronze_dir, "compounds.csv")
	assay_path = os.path.join(bronze_dir, "assays.csv")
	meas_path = os.path.join(bronze_dir, "measurements.csv")
	panel_meta_path = os.path.join(bronze_dir, "panel_block_meta.csv")

	# panel 메타 집계 마무리: table_id별 max_test_conc 산출 및 측정치에 주입
	from collections import Counter
	max_conc_map: Dict[str, str] = {}
	for tbl, acc in panel_acc.items():
		thrs = acc.get("thresholds", []) or []
		thr_val: str = ""
		if len(thrs) > 0:
			cnt = Counter([f"{float(x):.8g}" for x in thrs])
			thr_val = max(cnt.items(), key=lambda kv: kv[1])[0]
		max_conc_map[tbl] = thr_val
	if len(meas_df) > 0:
		if "table_id" in meas_df.columns:
			meas_df["max_test_conc"] = meas_df["table_id"].map(lambda t: max_conc_map.get(str(t), ""))

	# panel_block_meta.csv 생성
	if str(cfg.get("dataset_id", "")) == "2017" and len(panel_meta_acc) > 0:
		panel_rows: List[Dict[str, str]] = []
		for pid, meta in sorted(panel_meta_acc.items(), key=lambda kv: kv[0]):
			cells = sorted(list(meta.get("cell_lines", set()))) if isinstance(meta.get("cell_lines"), (set, list)) else []
			panel_rows.append({
				"panel_id": str(pid),
				"panel_label": str(meta.get("panel_label", "")),
				"disease_area": str(meta.get("disease_area", "")),
				"provenance_sheet": str(meta.get("provenance_sheet", "")),
				"cell_lines": ";".join(cells),
			})
		panel_meta_df = pd.DataFrame(panel_rows, columns=["panel_id", "panel_label", "disease_area", "provenance_sheet", "cell_lines"]).drop_duplicates()
	else:
		panel_rows: List[Dict[str, str]] = []
		for tbl, acc in sorted(panel_acc.items(), key=lambda kv: kv[0]):
			cell_lines_sorted = ",".join(sorted(list(acc.get("cell_lines", set()))))
			panel_rows.append({
				"table_id": str(tbl),
				"assay_label": str(acc.get("assay_label", "")),
				"cell_lines": cell_lines_sorted,
				"max_test_conc": max_conc_map.get(str(tbl), ""),
				"panel_name": str(panel_name_map.get(str(tbl), panel_name_map.get(tbl, ""))),
			})
		panel_meta_df = pd.DataFrame(panel_rows, columns=["table_id", "assay_label", "cell_lines", "max_test_conc", "panel_name"]).drop_duplicates()

	write_csv_safe(comp_df, comp_path)
	write_csv_safe(assay_df, assay_path)
	write_csv_safe(meas_df, meas_path)
	write_csv_safe(panel_meta_df, panel_meta_path)

	outputs.update({
		"compounds": comp_path,
		"assays": assay_path,
		"measurements": meas_path,
		"panel_meta": panel_meta_path,
	})
	return outputs

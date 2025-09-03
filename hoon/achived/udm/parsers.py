import re
from dataclasses import dataclass


MU = "μ"


@dataclass
class AuditNote:
	field: str
	original: str
	corrected: str
	message: str


def normalize_unit_label(unit: str | None) -> str | None:
	"""단위 레이블 변형 정규화: μΜ/uM->μM, None/empty는 그대로 유지."""
	if unit is None:
		return None
	orig = unit
	val = unit.strip()
	# 변종 통일: 'μΜ'(그리스 대문자 Mu), 'µ'(micro sign), 'uM', 'μ M'
	val = (
		val
		.replace("μΜ", "μM")
		.replace("µM", f"{MU}M")
		.replace("uM", f"{MU}M")
		.replace("μ M", f"{MU}M")
	)
	return val


_CENSOR_RE = re.compile(r"^\s*((?:>=|<=|[<>＝=＜＞])?)[\s\t]*([+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?)\s*(\S+)?\s*$")

# '2-.9' 또는 '1-.7'과 같은 드문 잘못된 형식의 숫자 토큰에 대한 대체 처리
_MALFORMED_NUM_RE = re.compile(r"^(\d+)\s*[-\u2212]\s*\.?(\d+)\s*$")


def split_censor_and_value(value_raw: str | None) -> tuple[str | None, str | None, str | None]:
	"""
	원본 값 문자열에서 검열 기호와 숫자 부분을 분리합니다.

	반환값: (censor, number, unit_token)
	- censor: {'gt','lt','eq', None} 중 하나
	- number: 숫자 리터럴 문자열, 파싱되지 않은 경우 None
	- unit_token: 후행 공백이 아닌 토큰 (예: 'μM'), 없는 경우 None
	"""
	if value_raw is None:
		return None, None, None
	text = str(value_raw).strip()
	if text == "":
		return None, None, None
	m = _CENSOR_RE.match(text)
	if not m:
		# 공백으로 분리하여 단위와 잘못된 형식의 숫자를 탐지
		parts = text.split()
		if len(parts) >= 1:
			cand = parts[0]
			mm = _MALFORMED_NUM_RE.match(cand)
			if mm:
				# 예: '2-.9' -> '2.9'
				number = f"{mm.group(1)}.{mm.group(2)}"
				unit_token = parts[1] if len(parts) > 1 else None
				return None, number, unit_token
		return None, None, None
	sign, number, unit_token = m.groups()
	censor = None
	if sign == ">":
		censor = "gt"
	elif sign == "<":
		censor = "lt"
	elif sign == ">=":
		censor = "ge"
	elif sign == "<=":
		censor = "le"
	elif sign in {"=", "＝", "＜", "＞"}:
		censor = "eq"
	return censor, number, unit_token


_SCI_RE = re.compile(r"^([+-]?(?:\d+(?:\.\d*)?|\.\d+))\s*([tTeE])\s*([+-]?\d+)$")


def parse_scientific_notation(num_text: str) -> tuple[float | None, AuditNote | None]:
	"""
	8.8E-07 또는 8.8TE-07과 같은 숫자를 파싱합니다 (잘못된 'T' 포함).
	'T'가 발견되면 'E'로 수정하고 AuditNote를 반환합니다.
	"""
	if num_text is None:
		return None, None
	text = num_text.strip()
	m = _SCI_RE.match(text)
	if not m:
		# 일반적인 float 변환 시도
		try:
			return float(text), None
		except Exception:
			return None, None
	base, exponent_letter, exponent_power = m.groups()
	audit: AuditNote | None = None
	if exponent_letter in {"t", "T"}:
		corrected = f"{base}E{exponent_power}"
		audit = AuditNote(field="value_num", original=text, corrected=corrected, message="과학적 표기법에서 'TE'를 'E'로 수정")
		text = corrected
	try:
		return float(text), audit
	except Exception:
		return None, audit



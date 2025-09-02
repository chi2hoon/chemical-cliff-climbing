import re

from pipeline.parsers.measures import parse_measure


_ASCII_UNIT_MAP = {
    "µM": "uM",
    "μM": "uM",
    "uM": "uM",
    "nM": "nM",
    "M": "M",
    "%": "%",
    "mM": "mM",
}


def ascii_units(u):
    """Args: u(str|None) -> str|None

    µ/μ를 u로 통일(ASCII). 알 수 없는 단위는 원문 반환.
    """
    if u is None:
        return None
    return _ASCII_UNIT_MAP.get(str(u), str(u).replace("µ", "u").replace("μ", "u"))


def strip_all(s):
    """Args: s(str|None) -> str|None

    문자열 양끝 공백 제거. None은 그대로.
    """
    if s is None:
        return None
    return str(s).strip()


def dash_to_nan(s):
    """Args: s(str|None) -> str|None

    "-" 또는 빈 문자열을 None으로 변환.
    """
    if s is None:
        return None
    s2 = str(s).strip()
    if s2 in {"", "-", "--", "—"}:
        return None
    return s


def parse_qual_value_unit(s, default_unit=None):
    """Args: s(str|None), default_unit(str|None) -> (qualifier,str|None,str|None)

    측정 텍스트에서 qualifier/value/unit을 분리하고 ASCII 단위로 정규화한다.
    """
    q, v, u = parse_measure(s, default_unit=default_unit)
    return q, v, ascii_units(u)


def _to_float(v):
    try:
        return float(v)
    except Exception:
        return None


def convert_unit(value_str, unit_in, unit_out):
    """Args: value_str(str|None), unit_in(str|None), unit_out(str|None) -> (float|None,str|None)

    값 문자열을 unit_out(예: uM,nM,%) 기준으로 변환. 불가능하면 원 단위 유지.
    반환: (value_out, unit_out_ascii)
    """
    if value_str is None:
        return None, ascii_units(unit_out or unit_in)
    vin = _to_float(value_str)
    if vin is None:
        return None, ascii_units(unit_out or unit_in)
    ui = ascii_units(unit_in)
    uo = ascii_units(unit_out or unit_in)
    if ui == uo or uo is None:
        return vin, ui
    # 농도 단위 변환(M/mM/uM/nM)
    factors = {
        ("M", "uM"): 1e6,
        ("M", "nM"): 1e9,
        ("mM", "uM"): 1e3,
        ("mM", "nM"): 1e6,
        ("uM", "nM"): 1e3,
        ("nM", "uM"): 1e-3,
        ("uM", "M"): 1e-6,
        ("nM", "M"): 1e-9,
        ("uM", "mM"): 1e-3,
        ("nM", "mM"): 1e-6,
    }
    if ui == "%" and uo == "%":
        return vin, "%"
    f = factors.get((ui, uo))
    if f is None:
        # 변환 불가 → 원래 단위 유지
        return vin, ui
    return vin * f, uo


def normalize_mortality(text):
    """Args: text(str|None) -> str|None

    치사율 문자열 내 a.b를 a/10과 같은 규칙으로 정규화(옵션 기능). 미구현 시 원본 반환.
    """
    if text is None:
        return None
    s = str(text).strip()
    if re.fullmatch(r"\d+\.\d+", s):
        a, b = s.split(".", 1)
        if len(b) == 1:
            return f"{a}/{b}0"
    return text


def parse_dose(text):
    """Args: text(str|None) -> (route,str,unit)

    도스 문자열 파서: "IV 40 mg/kg" → ("IV","40","mg/kg")
    """
    if text is None:
        return None, None, None
    m = re.match(r"^([A-Za-z/]+)\s+([\d\.]+)\s*(\S+)?$", str(text).strip())
    if not m:
        return None, None, None
    route, val, unit = m.groups()
    return route, val, unit


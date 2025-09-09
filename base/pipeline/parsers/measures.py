import re


_UNIT_MAP = {
    "µM": "uM",
    "μM": "uM",
    "uM": "uM",
    "nM": "nM",
    "%": "%",
}


def parse_measure(s, default_unit=None):
    """Args: s(str|None), default_unit(str|None) -> (qualifier,str,str|None)

    측정 문자열에서 부등호/값/단위를 분리한다. 단위는 ASCII로 정규화(uM).
    예시: ">= 2 µM" -> ('>', '2', 'uM')
    """
    if s is None:
        return None, None, default_unit
    text = str(s).strip()
    if text == "" or text.lower() == "nan":
        return None, None, default_unit

    m = re.match(r"^(>=|>|<=|<|=)?\s*([+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?)\s*([a-zA-Z%µμ]+)?", text)
    if not m:
        return None, text, default_unit
    op, val, unit = m.groups()

    if op in (">=", ">"):
        q = ">"
    elif op in ("<=", "<"):
        q = "<"
    elif op == "=":
        q = "="
    else:
        q = "=" if val is not None else None

    unit_ascii = None
    if unit:
        unit_ascii = _UNIT_MAP.get(unit, unit)
        unit_ascii = unit_ascii.replace("µ", "u").replace("μ", "u")
    else:
        unit_ascii = default_unit

    return q, val, unit_ascii


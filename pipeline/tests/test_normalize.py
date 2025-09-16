from pipeline.transforms.normalize import parse_qual_value_unit, convert_unit


def test_parse_qual_value_unit_ascii():
    q, v, u = parse_qual_value_unit("> 2 µM", default_unit=None)
    assert q == ">"
    assert v == "2"
    assert u == "uM"


def test_convert_unit_m_to_uM():
    v, u = convert_unit("1e-6", "M", "uM")
    assert abs(v - 1.0) < 1e-9 and u == "uM"


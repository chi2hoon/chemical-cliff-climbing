from pipeline.parsers.measures import parse_measure


def test_parse_measure_basic():
    q, v, u = parse_measure(">= 2 µM", default_unit=None)
    assert q == ">"
    assert v == "2"
    assert u == "uM"


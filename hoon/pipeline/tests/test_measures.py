import pytest

from pipeline.parsers.measures import parse_measure


@pytest.mark.parametrize(
    "s, default, expected",
    [
        ("> 2 µM", None, (">", "2", "uM")),
        ("=2 uM", None, ("=", "2", "uM")),
        ("<= 10 nM", None, ("<", "10", "nM")),
        ("15%", None, ("=", "15", "%")),
        ("2", "uM", ("=", "2", "uM")),
        (None, "uM", (None, None, "uM")),
        ("", None, (None, None, None)),
    ],
)
def test_parse_measure(s, default, expected):
    assert parse_measure(s, default) == expected


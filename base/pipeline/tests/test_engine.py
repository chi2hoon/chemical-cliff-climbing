from pipeline.parsers.engine import matrix_to_long


def test_matrix_to_long_smoke():
    import pandas as pd
    # 간단한 2행 헤더와 1개 데이터 행
    data = [
        ["panel", "A", "A"],
        ["cell", "c1", "c2"],
        ["row1", "1 uM", "> 2 uM"],
    ]
    df = pd.DataFrame(data)
    cfg = {
        "panel_row_offset": 0,
        "cellline_row_offset": 1,
        "data_row_start": 2,
        "id_col": 0,
        "panels": [{"name_from_header": True, "cols": [1, 2]}],
        "value_parser": {"unit_default": "uM"},
    }
    out = matrix_to_long(df, cfg)
    assert len(out) == 2
    assert set(out["cell_line"]) == {"c1", "c2"}


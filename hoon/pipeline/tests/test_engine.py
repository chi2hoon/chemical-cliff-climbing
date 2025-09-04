import pandas as pd
from pipeline.parsers.engine import matrix_to_long


def test_matrix_to_long_smoke():
    # 간단한 2단 헤더 매트릭스 구성
    data = [
        ["", "PanelA", "PanelA", "PanelB"],
        ["ID", "CL1", "CL2", "CL3"],
        ["row1", "> 2 uM", "1 nM", "15%"],
        ["row2", "3", "4", "5"],
    ]
    df = pd.DataFrame(data)
    cfg = {
        "panel_row_offset": 0,
        "cellline_row_offset": 1,
        "data_row_start": 2,
        "id_col": 0,
        "panels": [
            {"name": "A", "cols": [1, 2]},
            {"name": "B", "cols": [3]},
        ],
        "value_parser": {"unit_default": "uM"},
    }
    out = matrix_to_long(df, cfg)
    assert set(["row_id", "panel", "cell_line", "value_raw", "qualifier", "value", "unit"]).issubset(out.columns)


from load_data import load_csv_to_db
from src.analysis import get_cell_frequency_data
from src.statistics import compare_responders


def setup_module() -> None:
    load_csv_to_db()


def test_cell_frequency_columns() -> None:
    df = get_cell_frequency_data()
    expected = {"sample_id", "cell_type", "count", "total_count", "percentage"}
    assert expected.issubset(set(df.columns))
    assert len(df) > 0


def test_compare_responders_columns() -> None:
    stats_df, filtered_df = compare_responders()
    expected = {
        "cell_type",
        "p_value",
        "stat_score",
        "significant",
        "avg_responder",
        "avg_non_responder",
    }
    assert expected.issubset(set(stats_df.columns))
    assert len(filtered_df) > 0

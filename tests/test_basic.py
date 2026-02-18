from typing import cast

import pandas as pd
import plotly.express as px

from load_data import load_csv_to_db
from src.analysis import get_cell_frequency_data
from src.queries import get_subset_stats
from src.reporting import build_html_report, build_pdf_report
from src.statistics import compare_responders


def setup_module() -> None:
    load_csv_to_db()


def test_cell_frequency_columns() -> None:
    df = get_cell_frequency_data()
    expected = {"sample_id", "cell_type", "count", "total_count", "percentage"}
    assert expected.issubset(set(df.columns))
    assert len(df) > 0


def test_compare_responders_columns() -> None:
    stats_df, filtered_df, summary = compare_responders()
    expected = {
        "cell_type",
        "n_yes",
        "n_no",
        "p_value",
        "q_value",
        "stat_score",
        "effect",
        "cliffs_delta",
        "significant",
        "avg_responder",
        "avg_non_responder",
    }
    assert expected.issubset(set(stats_df.columns))
    assert len(filtered_df) > 0
    assert summary["correction_label"] == "BH-FDR"


def test_part4_subset_counts_are_consistent() -> None:
    stats = get_subset_stats(
        condition="melanoma",
        treatment="miraclib",
        sample_type="PBMC",
        time_filter="baseline_only",
    )

    by_project_samples = cast(pd.Series, stats["by_project_samples"])
    by_project_subjects = cast(pd.Series, stats["by_project_subjects"])
    n_samples = cast(int, stats["n_samples"])
    n_subjects = cast(int, stats["n_subjects"])
    assert n_samples == int(by_project_samples.sum())
    assert n_subjects == int(by_project_subjects.sum())


def test_compare_responders_clr_transform() -> None:
    stats_df, clr_plot_df, summary = compare_responders(transform="clr")
    assert len(stats_df) > 0
    assert len(clr_plot_df) > 0
    assert summary["transform_label"] == "CLR"

    unit_means = cast(pd.Series, clr_plot_df.groupby("unit_id")["metric_value"].mean())
    max_abs_mean = float(unit_means.abs().max())
    assert max_abs_mean < 1e-9


def test_report_builders_return_valid_bytes() -> None:
    stats_df = pd.DataFrame(
        [
            {
                "cell_type": "b_cell",
                "n_yes": 3,
                "n_no": 4,
                "effect": 0.2,
                "p_value": 0.03,
                "q_value": 0.05,
                "significant": False,
            }
        ]
    )
    flow_df = pd.DataFrame(
        [
            {"step": "All samples", "n_samples": 10, "n_subjects": 7},
            {"step": "Condition=melanoma", "n_samples": 6, "n_subjects": 5},
        ]
    )
    plot_df = pd.DataFrame(
        [
            {"cell_type": "b_cell", "metric_value": 1.0, "response": "yes"},
            {"cell_type": "b_cell", "metric_value": 0.8, "response": "no"},
        ]
    )
    fig = px.box(plot_df, x="cell_type", y="metric_value", color="response")

    html_bytes = build_html_report(
        filters_text="Indication=melanoma | Treatment=miraclib | SampleType=PBMC",
        cohort_counts={"n_samples": 6, "n_subjects": 5},
        summary={
            "test_label": "Mann-Whitney U",
            "correction_label": "BH-FDR",
            "unit": "sample",
            "metric": "percentage",
            "transform_label": "Raw",
        },
        stats_df=stats_df,
        flow_df=flow_df,
        fig=fig,
    )
    assert b"Clinical Trial Analysis Report" in html_bytes

    pdf_bytes = build_pdf_report(
        filters_text="Indication=melanoma | Treatment=miraclib | SampleType=PBMC",
        cohort_counts={"n_samples": 6, "n_subjects": 5},
        summary={
            "test_label": "Mann-Whitney U",
            "correction_label": "BH-FDR",
            "unit": "sample",
            "metric": "percentage",
            "transform_label": "Raw",
        },
        stats_df=stats_df,
        flow_df=flow_df,
    )
    assert pdf_bytes.startswith(b"%PDF")

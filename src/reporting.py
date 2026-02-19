from datetime import UTC, datetime
from io import BytesIO
from typing import Any

import pandas as pd
import plotly.io as pio


def build_html_report(
    filters_text: str,
    cohort_counts: dict[str, int],
    summary: dict[str, str],
    stats_df: pd.DataFrame,
    flow_df: pd.DataFrame,
    fig: Any,
) -> bytes:
    generated_at = datetime.now(UTC).isoformat()
    stats_html = stats_df.to_html(index=False)
    flow_html = flow_df.to_html(index=False)
    fig_html = pio.to_html(fig, full_html=False, include_plotlyjs=True)

    html = f"""
<html>
  <head>
    <meta charset=\"utf-8\" />
    <title>Clinical Trial Analysis Report</title>
    <style>
      body {{ font-family: Arial, sans-serif; margin: 24px; }}
      h1, h2 {{ margin-bottom: 8px; }}
      .meta {{ margin-bottom: 20px; }}
      table {{ border-collapse: collapse; width: 100%; margin-bottom: 16px; }}
      th, td {{ border: 1px solid #ddd; padding: 6px; font-size: 12px; }}
      th {{ background: #f5f5f5; }}
    </style>
  </head>
  <body>
    <h1>Loblaw Bio Clinical Trial Report</h1>
    <div class=\"meta\">Generated at: {generated_at}</div>
    <h2>Active Filters</h2>
    <div>{filters_text}</div>
    <h2>Cohort Counts</h2>
    <div>n_samples={cohort_counts['n_samples']} | n_subjects={cohort_counts['n_subjects']}</div>
    <h2>Method Summary</h2>
    <div>Test={summary['test_label']} | Correction={summary['correction_label']} | Unit={summary['unit']} | Metric={summary['metric']} | Transform={summary['transform_label']}</div>
    <h2>Statistical Results</h2>
    {stats_html}
    <h2>Distribution Plot</h2>
    {fig_html}
    <h2>Cohort Flow</h2>
    {flow_html}
  </body>
</html>
"""
    return html.encode("utf-8")


def build_pdf_report(
    filters_text: str,
    cohort_counts: dict[str, int],
    summary: dict[str, str],
    stats_df: pd.DataFrame,
    flow_df: pd.DataFrame,
) -> bytes:
    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_pdf import PdfPages

    buffer = BytesIO()
    with PdfPages(buffer) as pdf:
        fig1, ax1 = plt.subplots(figsize=(8.27, 11.69))
        ax1.axis("off")

        lines = [
            "Loblaw Bio Clinical Trial Report",
            f"Generated at: {datetime.now(UTC).isoformat()}",
            "",
            f"Filters: {filters_text}",
            f"Cohort: n_samples={cohort_counts['n_samples']} | n_subjects={cohort_counts['n_subjects']}",
            (
                f"Method: {summary['test_label']} | Correction={summary['correction_label']} "
                f"| Unit={summary['unit']} | Metric={summary['metric']} | Transform={summary['transform_label']}"
            ),
        ]
        ax1.text(0.02, 0.98, "\n".join(lines), va="top", family="monospace", fontsize=10)
        pdf.savefig(fig1, bbox_inches="tight")
        plt.close(fig1)

        stats_cols = [
            "cell_type",
            "n_yes",
            "n_no",
            "median_diff",
            "direction",
            "ci_95_low",
            "ci_95_high",
            "effect",
            "p_value",
            "q_value",
            "significant",
        ]
        stats_table = stats_df.reindex(columns=stats_cols).copy()
        if len(stats_table) == 0:
            stats_table = pd.DataFrame({col: [] for col in stats_cols})
        stats_table = stats_table.fillna("")

        fig2, ax2 = plt.subplots(figsize=(11.69, 8.27))
        ax2.axis("off")
        table = ax2.table(
            cellText=stats_table.astype(str).values.tolist(),
            colLabels=stats_table.columns.tolist(),
            loc="center",
        )
        table.auto_set_font_size(False)
        table.set_fontsize(8)
        table.scale(1, 1.2)
        ax2.set_title("Part 3 Statistical Results", fontsize=12)
        pdf.savefig(fig2, bbox_inches="tight")
        plt.close(fig2)

        flow_table = flow_df.fillna("")
        fig3, ax3 = plt.subplots(figsize=(11.69, 8.27))
        ax3.axis("off")
        flow = ax3.table(
            cellText=flow_table.astype(str).values.tolist(),
            colLabels=flow_table.columns.tolist(),
            loc="center",
        )
        flow.auto_set_font_size(False)
        flow.set_fontsize(9)
        flow.scale(1, 1.3)
        ax3.set_title("Cohort Flow", fontsize=12)
        pdf.savefig(fig3, bbox_inches="tight")
        plt.close(fig3)

    buffer.seek(0)
    return buffer.getvalue()

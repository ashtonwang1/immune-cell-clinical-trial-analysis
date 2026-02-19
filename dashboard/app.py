import os
from typing import cast

import pandas as pd
import plotly.express as px
import streamlit as st

from src.config import DB_PATH

# Auto-generate database on first run (e.g. Streamlit Cloud)
if not os.path.exists(DB_PATH):
    from load_data import load_csv_to_db

    load_csv_to_db()

from src.analysis import get_cohort_counts, get_filter_options, get_filtered_data, get_part2_frequency_table
from src.config import CELL_TYPES
from src.queries import build_cohort_flow, get_subset_stats
from src.reporting import build_html_report, build_pdf_report
from src.statistics import compare_responders


@st.cache_data(show_spinner=False)
def cached_filter_options() -> dict[str, list[str]]:
    return get_filter_options()


@st.cache_data(show_spinner=False)
def cached_filtered_data(
    condition: str,
    treatment: str,
    sample_type: str,
    time_filter: str,
) -> pd.DataFrame:
    return get_filtered_data(
        condition=condition,
        treatment=treatment,
        sample_type=sample_type,
        time_filter=time_filter,
    )


@st.cache_data(show_spinner=False)
def cached_part2_table() -> pd.DataFrame:
    return get_part2_frequency_table()


@st.cache_data(show_spinner=False)
def cached_compare_responders(
    condition: str,
    treatment: str,
    sample_type: str,
    time_filter: str,
    unit: str,
    metric: str,
    transform: str,
    test: str,
    correction: str,
) -> tuple[pd.DataFrame, pd.DataFrame, dict[str, str]]:
    return compare_responders(
        condition=condition,
        treatment=treatment,
        sample_type=sample_type,
        time_filter=time_filter,
        unit=unit,
        metric=metric,
        transform=transform,
        test=test,
        correction=correction,
    )


@st.cache_data(show_spinner=False)
def cached_subset_stats(
    condition: str,
    treatment: str,
    sample_type: str,
    time_filter: str,
) -> dict[str, pd.Series | pd.DataFrame | int | float | None]:
    return get_subset_stats(
        condition=condition,
        treatment=treatment,
        sample_type=sample_type,
        time_filter=time_filter,
    )


@st.cache_data(show_spinner=False)
def cached_cohort_flow(
    condition: str,
    treatment: str,
    sample_type: str,
    time_filter: str,
) -> pd.DataFrame:
    return build_cohort_flow(
        condition=condition,
        treatment=treatment,
        sample_type=sample_type,
        time_filter=time_filter,
    )


def to_csv_bytes(df: pd.DataFrame) -> bytes:
    return df.to_csv(index=False).encode("utf-8")


def format_sensitivity_status(is_significant: bool, q_value: float | None) -> str:
    if q_value is None or pd.isna(q_value):
        return "insufficient data"
    state = "significant" if is_significant else "not significant"
    return f"{state} (q={q_value:.3f})"


st.set_page_config(page_title="Loblaw Bio Analysis", layout="wide")

st.title("Loblaw Bio: Clinical Trial Analysis")
st.markdown("Industrial-grade cohort analytics for immune-cell populations in Miraclib clinical trial data.")

options = cached_filter_options()

default_condition = "melanoma" if "melanoma" in [v.lower() for v in options["conditions"]] else options["conditions"][0]
default_treatment = "miraclib" if "miraclib" in [v.lower() for v in options["treatments"]] else options["treatments"][0]
default_sample_type = "PBMC" if "PBMC" in options["sample_types"] else options["sample_types"][0]

condition = cast(
    str,
    st.sidebar.selectbox("Indication", options["conditions"], index=options["conditions"].index(default_condition)),
)
treatment = cast(
    str,
    st.sidebar.selectbox("Treatment", options["treatments"], index=options["treatments"].index(default_treatment)),
)
sample_type = cast(
    str,
    st.sidebar.selectbox("Sample Type", options["sample_types"], index=options["sample_types"].index(default_sample_type)),
)
time_label = st.sidebar.selectbox("Time From Treatment Start", ["All", "Baseline only"], index=1)
time_filter = "baseline_only" if time_label == "Baseline only" else "all"

unit_label = st.sidebar.selectbox("Unit of Analysis", ["Sample", "Subject"], index=1)
unit = "sample" if unit_label == "Sample" else "subject"

metric_label = st.sidebar.selectbox("Metric", ["Percentage", "Count"], index=0)
metric = "percentage" if metric_label == "Percentage" else "count"

clr_transform_enabled = st.sidebar.toggle("Advanced: CLR transform (compositional)", value=False)
transform = "clr" if (clr_transform_enabled and metric == "percentage") else "none"

if clr_transform_enabled and metric != "percentage":
    st.sidebar.info("CLR transform is applied only when Metric=Percentage.")

show_all_points = st.sidebar.toggle("Show all points in boxplot", value=False)
point_mode = "all" if show_all_points else "outliers"

filtered_df = cached_filtered_data(condition, treatment, sample_type, time_filter)
cohort_counts = get_cohort_counts(filtered_df)

active_filters_text = (
    f"Indication={condition} | Treatment={treatment} | SampleType={sample_type} | "
    f"Time={'Baseline only' if time_filter == 'baseline_only' else 'All'} | Unit={unit_label} | "
    f"Metric={metric_label} | Transform={'CLR' if transform == 'clr' else 'Raw'}"
)
st.markdown(f"**Active Filters:** {active_filters_text}")

top_col_1, top_col_2 = st.columns(2)
top_col_1.metric("n_samples", cohort_counts["n_samples"])
top_col_2.metric("n_subjects", cohort_counts["n_subjects"])

tab_part2, tab_part3, tab_part4, tab_sensitivity, tab_methods = st.tabs(
    [
        "Data Overview (Part 2)",
        "Statistical Analysis (Part 3)",
        "Subset Analysis (Part 4)",
        "Sensitivity / Robustness",
        "Methods & Definitions",
    ]
)

with tab_part2:
    st.header("Frequency Table by Sample")

    part2_all_df = cached_part2_table()
    if len(filtered_df) == 0:
        part2_view = part2_all_df.iloc[0:0].copy()
    else:
        sample_ids = filtered_df["sample_id"].astype(str).unique().tolist()
        part2_view = part2_all_df.loc[part2_all_df["sample"].astype(str).isin(sample_ids)].copy()

    st.caption("Relative frequency (%) of each immune cell population per sample.")
    st.caption(f"Rows under active filters: {len(part2_view)}")

    if len(part2_view) == 0:
        st.warning("No Part 2 rows available for the current filter set.")
    else:
        st.dataframe(part2_view, use_container_width=True, hide_index=True)

    p2_col_1, p2_col_2 = st.columns(2)
    p2_col_1.download_button(
        "Download Part 2 table (active filters)",
        data=to_csv_bytes(part2_view),
        file_name="part2_frequency_table_filtered.csv",
        mime="text/csv",
    )
    p2_col_2.download_button(
        "Download Part 2 table (all samples)",
        data=to_csv_bytes(part2_all_df),
        file_name="part2_frequency_table.csv",
        mime="text/csv",
    )

with tab_part3:
    st.header("Responder vs Non-Responder")

    stats_df, plot_df, summary = cached_compare_responders(
        condition=condition,
        treatment=treatment,
        sample_type=sample_type,
        time_filter=time_filter,
        unit=unit,
        metric=metric,
        transform=transform,
        test="mannwhitney",
        correction="bh_fdr",
    )

    st.caption(
        f"Test: {summary['test_label']} | Multiple testing: {summary['correction_label']} | "
        f"Unit: {summary['unit']} | Metric: {summary['metric']} | {summary['bootstrap_ci']}"
    )

    if len(stats_df) == 0:
        st.warning("No records available for the current filter set.")
    else:
        table_df = stats_df.loc[
            :,
            [
                "cell_type",
                "n_yes",
                "n_no",
                "median_yes",
                "median_no",
                "median_diff",
                "direction",
                "ci_95_low",
                "ci_95_high",
                "effect",
                "cliffs_delta",
                "p_value",
                "q_value",
                "significant",
            ],
        ].copy()

        table_df["significant"] = table_df["significant"].map(
            lambda x: "significant (q<0.05)" if bool(x) else "not significant"
        )

        st.dataframe(table_df, use_container_width=True, hide_index=True)

        dcol1, dcol2 = st.columns(2)
        dcol1.download_button(
            "Download stats table (CSV)",
            data=to_csv_bytes(table_df),
            file_name="part3_stats.csv",
            mime="text/csv",
        )
        dcol2.download_button(
            "Download filtered analysis data (CSV)",
            data=to_csv_bytes(plot_df),
            file_name="part3_filtered_data.csv",
            mime="text/csv",
        )

        fig = px.box(
            plot_df,
            x="cell_type",
            y="metric_value",
            color="response",
            points=point_mode,
            category_orders={"cell_type": CELL_TYPES, "response": ["no", "yes"]},
            color_discrete_map={"yes": "#1f9d55", "no": "#d64545"},
            labels={"metric_value": metric_label, "response": "Response", "cell_type": "Cell Type"},
            title="Distribution by Response Group",
        )

        max_by_cell = plot_df.groupby("cell_type")["metric_value"].max().to_dict()
        max_global = float(plot_df["metric_value"].max()) if len(plot_df) > 0 else 0.0
        offset = max(1.0, max_global * 0.08)

        for _, stat_row in stats_df.iterrows():
            cell = str(stat_row["cell_type"])
            q_obj = stat_row["q_value"]
            q_val = float(q_obj) if isinstance(q_obj, (int, float)) else None
            text = f"q={q_val:.3g}" if q_val is not None else "q=NA"
            y = float(max_by_cell.get(cell, max_global)) + offset
            fig.add_annotation(x=cell, y=y, text=text, showarrow=False, font={"size": 11})

        if max_global > 0:
            fig.update_yaxes(range=[0, max_global + 2 * offset])

        st.plotly_chart(fig, use_container_width=True)

        flow_for_report = cached_cohort_flow(condition, treatment, sample_type, time_filter)
        html_report = build_html_report(
            filters_text=active_filters_text,
            cohort_counts=cohort_counts,
            summary=summary,
            stats_df=stats_df,
            flow_df=flow_for_report,
            fig=fig,
        )
        pdf_report = build_pdf_report(
            filters_text=active_filters_text,
            cohort_counts=cohort_counts,
            summary=summary,
            stats_df=stats_df,
            flow_df=flow_for_report,
        )

        r1, r2 = st.columns(2)
        r1.download_button(
            "Generate report (HTML)",
            data=html_report,
            file_name="analysis_report.html",
            mime="text/html",
        )
        r2.download_button(
            "Generate report (PDF)",
            data=pdf_report,
            file_name="analysis_report.pdf",
            mime="application/pdf",
        )

with tab_part4:
    st.header("Baseline / Cohort Characterization")

    subset_stats = cached_subset_stats(condition, treatment, sample_type, time_filter)
    flow_df = cached_cohort_flow(condition, treatment, sample_type, time_filter)

    project_count = cast(int, subset_stats["n_projects"])
    total_samples = cast(int, subset_stats["n_samples"])
    total_subjects = cast(int, subset_stats["n_subjects"])
    avg_b = subset_stats["avg_b_cell_male_responders"]
    avg_b_text = f"{avg_b:.2f}" if isinstance(avg_b, float) else "N/A"

    k1, k2, k3, k4 = st.columns(4)
    k1.metric("Projects", project_count)
    k2.metric("Total Samples (cohort)", total_samples)
    k3.metric("Total Subjects (cohort)", total_subjects)
    k4.metric("Avg B-cell Count (Male Responders)", avg_b_text)
    st.caption("Average B-cell count is computed at subject level (mean within subject, then cohort mean).")

    st.subheader("Cohort Flow")
    st.dataframe(flow_df, use_container_width=True, hide_index=True)
    st.bar_chart(flow_df.set_index("step")[["n_samples", "n_subjects"]])

    p_samples = cast(pd.Series, subset_stats["by_project_samples"]).rename("n_samples").reset_index()
    p_samples.columns = ["project_id", "n_samples"]
    p_subjects = cast(pd.Series, subset_stats["by_project_subjects"]).rename("n_subjects").reset_index()
    p_subjects.columns = ["project_id", "n_subjects"]

    c1, c2 = st.columns(2)
    c1.subheader("Samples by Project")
    c1.dataframe(p_samples, use_container_width=True, hide_index=True)
    c2.subheader("Subjects by Project")
    c2.dataframe(p_subjects, use_container_width=True, hide_index=True)

    r_col, s_col = st.columns(2)
    r_col.subheader("Subjects by Response")
    r_col.bar_chart(cast(pd.Series, subset_stats["by_response"]))
    s_col.subheader("Subjects by Sex")
    s_col.bar_chart(cast(pd.Series, subset_stats["by_sex"]))

    raw_subset = cast(pd.DataFrame, subset_stats["df_raw"])

    e1, e2 = st.columns(2)
    e1.download_button(
        "Download subset raw data (CSV)",
        data=to_csv_bytes(raw_subset),
        file_name="part4_subset_raw.csv",
        mime="text/csv",
    )
    e2.download_button(
        "Download cohort flow (CSV)",
        data=to_csv_bytes(flow_df),
        file_name="cohort_flow.csv",
        mime="text/csv",
    )

    with st.expander("Show Query Logic"):
        st.code(
            """
WHERE LOWER(sub.condition)=<condition>
  AND LOWER(sub.treatment)=<treatment>
  AND LOWER(s.sample_type)=<sample_type>
  AND (s.visit_time=0 when Time=Baseline only)
            """.strip(),
            language="sql",
        )

    with st.expander("View Raw Subset Data"):
        st.dataframe(raw_subset, use_container_width=True)

with tab_sensitivity:
    st.header("Sensitivity / Robustness")

    scenario_configs = [
        ("Baseline | MW | BH-FDR", "baseline_only", "mannwhitney", "bh_fdr"),
        ("All Time | MW | BH-FDR", "all", "mannwhitney", "bh_fdr"),
        ("Baseline | Welch t | BH-FDR", "baseline_only", "welch_t", "bh_fdr"),
        ("Baseline | MW | None", "baseline_only", "mannwhitney", "none"),
    ]

    scenario_sig: dict[str, dict[str, bool]] = {}
    scenario_q: dict[str, dict[str, float | None]] = {}
    for label, sc_time, sc_test, sc_corr in scenario_configs:
        s_df, _, _ = cached_compare_responders(
            condition=condition,
            treatment=treatment,
            sample_type=sample_type,
            time_filter=sc_time,
            unit=unit,
            metric=metric,
            transform=transform,
            test=sc_test,
            correction=sc_corr,
        )
        scenario_sig[label] = {
            str(row["cell_type"]): bool(row["significant"]) for _, row in s_df.iterrows()
        }
        scenario_q[label] = {}
        for _, row in s_df.iterrows():
            q_obj = row["q_value"]
            q_val: float | None = None
            if isinstance(q_obj, (int, float)):
                q_float = float(q_obj)
                if not pd.isna(q_float):
                    q_val = q_float
            scenario_q[label][str(row["cell_type"])] = q_val

    all_cells = sorted(set().union(*[set(v.keys()) for v in scenario_sig.values()]))
    reference_label = scenario_configs[0][0]

    rows = []
    for cell in all_cells:
        ref_value = scenario_sig[reference_label].get(cell, False)
        cell_row: dict[str, str] = {"cell_type": str(cell)}
        stable = True
        for label, _, _, _ in scenario_configs:
            current = scenario_sig[label].get(cell, False)
            q_val = scenario_q[label].get(cell)
            cell_row[label] = format_sensitivity_status(current, q_val)
            stable = stable and (current == ref_value)
        cell_row["robustness"] = "stable" if stable else "sensitive"
        rows.append(cell_row)

    robustness_df = pd.DataFrame(rows)
    if len(robustness_df) == 0:
        st.info("No testable cell types under the current filters. Relax filters to run sensitivity checks.")
    else:
        scenario_summary = []
        for label, _, _, _ in scenario_configs:
            tested = len(scenario_sig[label])
            significant = int(sum(scenario_sig[label].values()))
            scenario_summary.append(f"{label}: {significant}/{tested} significant")
        st.caption("Scenario signal counts: " + " | ".join(scenario_summary))
        st.dataframe(robustness_df, use_container_width=True, hide_index=True)

with tab_methods:
    st.header("Methods & Definitions")
    st.markdown(
        """
- Response groups use `response=yes` vs `response=no`.
- Main inference uses Mann-Whitney U, with BH-FDR correction across tested cell types.
- Significance flag is based on `q_value < 0.05`.
- Unit-of-analysis can be sample-level or subject-level median aggregation.
- Percentage metric is compositional and should be interpreted with domain caution.
- Cohort flow reports both unique samples and unique subjects at each filtering step.
        """
    )

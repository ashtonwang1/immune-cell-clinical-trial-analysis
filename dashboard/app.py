import os
import sys

import plotly.express as px
import streamlit as st


sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

from src.queries import get_baseline_melanoma_stats, get_male_responder_b_cell_avg
from src.statistics import compare_responders


st.set_page_config(page_title="Loblaw Bio Analysis", layout="wide")

st.title("Loblaw Bio: Clinical Trial Analysis")
st.markdown("Analysis of immune cell populations in Melanoma patients treated with Miraclib.")

tab1, tab2 = st.tabs(["Statistical Analysis (Part 3)", "Subset Analysis (Part 4)"])


with tab1:
    st.header("Responder vs. Non-Responder Analysis")
    st.markdown("Comparison of cell population relative frequencies (PBMC samples).")

    stats_df, plot_data = compare_responders()

    st.subheader("Statistical Significance (Mann-Whitney U)")

    def highlight_significant(row) -> list[str]:
        is_sig = bool(row["significant"])
        return ["background-color: #d1e7dd" if is_sig else "" for _ in range(len(row))]

    styled = (
        stats_df[["cell_type", "p_value", "significant", "avg_responder", "avg_non_responder"]]
        .style.apply(highlight_significant, axis=1)
        .format(
            {
                "p_value": "{:.4f}",
                "avg_responder": "{:.2f}%",
                "avg_non_responder": "{:.2f}%",
            }
        )
    )
    st.dataframe(styled, use_container_width=True)

    st.subheader("Population Distributions")
    fig = px.box(
        plot_data,
        x="cell_type",
        y="percentage",
        color="response",
        points="all",
        title="Cell Frequency Distribution by Response Group",
        labels={"percentage": "Relative Frequency (%)", "response": "Response"},
        color_discrete_map={"yes": "#28a745", "no": "#dc3545"},
    )
    st.plotly_chart(fig, use_container_width=True)


with tab2:
    st.header("Baseline Characterization (Part 4)")
    st.markdown("Cohort analysis for Melanoma patients at Time 0 (Baseline).")

    stats = get_baseline_melanoma_stats()
    avg_b_cell = get_male_responder_b_cell_avg()

    col1, col2, col3 = st.columns(3)
    with col1:
        st.metric("Samples per Project", f"{len(stats['by_project'])} Projects")
    with col2:
        st.metric("Total Subjects", f"{int(stats['by_project'].sum())}")
    with col3:
        st.metric("Avg B-Cells (Male Responders)", f"{avg_b_cell:.2f}")

    st.divider()

    c1, c2, c3 = st.columns(3)
    with c1:
        st.subheader("By Project")
        st.dataframe(stats["by_project"], use_container_width=True)
    with c2:
        st.subheader("By Response")
        st.bar_chart(stats["by_response"])
    with c3:
        st.subheader("By Gender")
        st.bar_chart(stats["by_sex"])

    with st.expander("View Raw Subset Data"):
        st.dataframe(stats["df_raw"], use_container_width=True)

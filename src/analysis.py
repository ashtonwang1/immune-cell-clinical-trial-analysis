from typing import cast

import numpy as np
import pandas as pd

from src.database import get_db_connection


def get_cell_frequency_data() -> pd.DataFrame:
    conn = get_db_connection()
    query = """
    SELECT
        s.sample_id,
        sub.subject_pk,
        sub.subject_id,
        sub.project_id,
        sub.treatment,
        sub.response,
        sub.condition,
        sub.sex,
        s.sample_type,
        s.visit_time,
        c.cell_type,
        c.count
    FROM samples s
    JOIN subjects sub ON s.subject_pk = sub.subject_pk
    JOIN cell_counts c ON s.sample_id = c.sample_id
    """

    df = cast(pd.DataFrame, pd.read_sql_query(query, conn))
    conn.close()

    df["total_count"] = df.groupby("sample_id")["count"].transform("sum")
    df["percentage"] = (df["count"] / df["total_count"]) * 100

    df["condition"] = df["condition"].astype(str)
    df["treatment"] = df["treatment"].astype(str)
    df["sample_type"] = df["sample_type"].astype(str)
    df["response"] = df["response"].astype(str).str.lower()
    df["sex"] = df["sex"].astype(str)

    return df


def get_part2_frequency_table() -> pd.DataFrame:
    df = get_cell_frequency_data().copy()
    out = df.loc[:, ["sample_id", "total_count", "cell_type", "count", "percentage"]].copy()
    out = out.rename(columns={"sample_id": "sample", "cell_type": "population"})
    return cast(pd.DataFrame, out.loc[:, ["sample", "total_count", "population", "count", "percentage"]])


def get_filtered_data(
    condition: str = "melanoma",
    treatment: str = "miraclib",
    sample_type: str = "PBMC",
    time_filter: str = "all",
) -> pd.DataFrame:
    df = get_cell_frequency_data()
    mask = pd.Series(True, index=df.index)

    if condition and condition != "all":
        mask &= df["condition"].str.lower() == condition.lower()
    if treatment and treatment != "all":
        mask &= df["treatment"].str.lower() == treatment.lower()
    if sample_type and sample_type != "all":
        mask &= df["sample_type"].str.lower() == sample_type.lower()
    if time_filter == "baseline_only":
        mask &= df["visit_time"] == 0

    return cast(pd.DataFrame, df.loc[mask].copy())


def get_filter_options() -> dict[str, list[str]]:
    df = get_cell_frequency_data()
    return {
        "conditions": sorted(df["condition"].dropna().astype(str).unique().tolist()),
        "treatments": sorted(df["treatment"].dropna().astype(str).unique().tolist()),
        "sample_types": sorted(df["sample_type"].dropna().astype(str).unique().tolist()),
    }


def get_cohort_counts(df: pd.DataFrame) -> dict[str, int]:
    subject_col = "subject_pk" if "subject_pk" in df.columns else "subject_id"
    return {
        "n_samples": int(df["sample_id"].nunique()),
        "n_subjects": int(df[subject_col].nunique()),
    }


def prepare_unit_level_data(
    df: pd.DataFrame,
    unit: str = "sample",
    metric: str = "percentage",
) -> pd.DataFrame:
    if metric not in {"percentage", "count"}:
        raise ValueError("metric must be 'percentage' or 'count'")

    value_col = metric

    if unit == "sample":
        output = df.loc[:, ["cell_type", "response", "sample_id", "subject_pk", value_col]].copy()
        output = output.rename(columns={value_col: "metric_value", "sample_id": "unit_id"})
        return cast(pd.DataFrame, output)

    if unit == "subject":
        grouped = df.groupby(["subject_pk", "response", "cell_type"], as_index=False)[value_col].median()
        grouped.columns = ["unit_id", "response", "cell_type", "metric_value"]
        return cast(pd.DataFrame, grouped)

    raise ValueError("unit must be 'sample' or 'subject'")


def apply_clr_transform(unit_df: pd.DataFrame, pseudocount: float = 1e-6) -> pd.DataFrame:
    pivot = cast(
        pd.DataFrame,
        unit_df.pivot_table(index="unit_id", columns="cell_type", values="metric_value", aggfunc="mean"),
    )
    adjusted = cast(pd.DataFrame, pivot.fillna(0.0) + pseudocount)

    log_vals = cast(pd.DataFrame, adjusted.apply(np.log))
    clr_vals = cast(pd.DataFrame, log_vals.sub(log_vals.mean(axis=1), axis=0))

    melted = cast(
        pd.DataFrame,
        clr_vals.reset_index().melt(id_vars=["unit_id"], var_name="cell_type", value_name="metric_value"),
    )
    responses = cast(pd.DataFrame, unit_df.loc[:, ["unit_id", "response"]].drop_duplicates())
    merged = cast(pd.DataFrame, melted.merge(responses, on="unit_id", how="left"))
    merged = merged.loc[:, ["cell_type", "response", "unit_id", "metric_value"]].copy()
    return cast(pd.DataFrame, merged)

from typing import cast

import pandas as pd

from src.analysis import get_cell_frequency_data
from src.database import get_db_connection


def _subset_where_clause(time_filter: str) -> str:
    base = """
    WHERE LOWER(sub.condition) = ?
      AND LOWER(sub.treatment) = ?
      AND LOWER(s.sample_type) = ?
    """
    if time_filter == "baseline_only":
        base += " AND s.visit_time = ?"
    return base


def _subset_params(
    condition: str,
    treatment: str,
    sample_type: str,
    time_filter: str,
) -> list[str | float]:
    params: list[str | float] = [condition.lower(), treatment.lower(), sample_type.lower()]
    if time_filter == "baseline_only":
        params.append(0.0)
    return params


def count_samples_by_project(
    condition: str,
    treatment: str,
    sample_type: str,
    time_filter: str,
) -> pd.Series:
    conn = get_db_connection()
    query = f"""
    SELECT sub.project_id, COUNT(DISTINCT s.sample_id) AS n_samples
    FROM samples s
    JOIN subjects sub ON s.subject_pk = sub.subject_pk
    {_subset_where_clause(time_filter)}
    GROUP BY sub.project_id
    """
    params = _subset_params(condition, treatment, sample_type, time_filter)
    df = cast(pd.DataFrame, pd.read_sql_query(query, conn, params=params))
    conn.close()
    if len(df) == 0:
        return pd.Series(dtype="int64")
    return cast(pd.Series, df.set_index("project_id")["n_samples"])


def count_subjects_by_project(
    condition: str,
    treatment: str,
    sample_type: str,
    time_filter: str,
) -> pd.Series:
    conn = get_db_connection()
    query = f"""
    SELECT sub.project_id, COUNT(DISTINCT sub.subject_pk) AS n_subjects
    FROM samples s
    JOIN subjects sub ON s.subject_pk = sub.subject_pk
    {_subset_where_clause(time_filter)}
    GROUP BY sub.project_id
    """
    params = _subset_params(condition, treatment, sample_type, time_filter)
    df = cast(pd.DataFrame, pd.read_sql_query(query, conn, params=params))
    conn.close()
    if len(df) == 0:
        return pd.Series(dtype="int64")
    return cast(pd.Series, df.set_index("project_id")["n_subjects"])


def count_subjects_by_response_and_sex(
    condition: str,
    treatment: str,
    sample_type: str,
    time_filter: str,
) -> pd.DataFrame:
    conn = get_db_connection()
    query = f"""
    SELECT response, sex, COUNT(*) AS n_subjects
    FROM (
        SELECT DISTINCT sub.subject_pk, LOWER(sub.response) AS response, sub.sex AS sex
        FROM samples s
        JOIN subjects sub ON s.subject_pk = sub.subject_pk
        {_subset_where_clause(time_filter)}
    ) dedup
    GROUP BY response, sex
    """
    params = _subset_params(condition, treatment, sample_type, time_filter)
    df = cast(pd.DataFrame, pd.read_sql_query(query, conn, params=params))
    conn.close()
    return df


def avg_b_cell_male_responders_baseline(
    condition: str,
    treatment: str,
    sample_type: str,
    time_filter: str,
) -> float | None:
    conn = get_db_connection()
    query = f"""
    WITH per_subject AS (
        SELECT sub.subject_pk, AVG(c.count) AS subject_mean_b
        FROM samples s
        JOIN subjects sub ON s.subject_pk = sub.subject_pk
        JOIN cell_counts c ON s.sample_id = c.sample_id
        {_subset_where_clause(time_filter)}
          AND sub.sex = 'M'
          AND LOWER(sub.response) = 'yes'
          AND c.cell_type = 'b_cell'
        GROUP BY sub.subject_pk
    )
    SELECT AVG(subject_mean_b) AS avg_b FROM per_subject
    """
    params = _subset_params(condition, treatment, sample_type, time_filter)
    row = cast(pd.DataFrame, pd.read_sql_query(query, conn, params=params))
    conn.close()
    if len(row) == 0 or pd.isna(row.loc[0, "avg_b"]):
        return None
    return float(row.loc[0, "avg_b"])


def _fetch_subset(
    condition: str,
    treatment: str,
    sample_type: str,
    time_filter: str,
) -> pd.DataFrame:
    conn = get_db_connection()

    query = f"""
    SELECT
        sub.project_id,
        sub.subject_pk,
        s.sample_id,
        sub.subject_id,
        LOWER(sub.response) AS response,
        sub.sex,
        s.visit_time,
        c.cell_type,
        c.count
    FROM samples s
    JOIN subjects sub ON s.subject_pk = sub.subject_pk
    JOIN cell_counts c ON s.sample_id = c.sample_id
    {_subset_where_clause(time_filter)}
    """

    params = _subset_params(condition, treatment, sample_type, time_filter)

    df = cast(pd.DataFrame, pd.read_sql_query(query, conn, params=params))
    conn.close()
    return df


def get_subset_stats(
    condition: str,
    treatment: str,
    sample_type: str,
    time_filter: str,
) -> dict[str, pd.Series | pd.DataFrame | int | float | None]:
    df = _fetch_subset(condition, treatment, sample_type, time_filter)

    by_project_samples = count_samples_by_project(condition, treatment, sample_type, time_filter)
    by_project_subjects = count_subjects_by_project(condition, treatment, sample_type, time_filter)
    by_response_sex = count_subjects_by_response_and_sex(condition, treatment, sample_type, time_filter)

    if len(by_response_sex) == 0:
        by_response = pd.Series(dtype="int64")
        by_sex = pd.Series(dtype="int64")
        n_subjects = 0
    else:
        by_response = cast(pd.Series, by_response_sex.groupby("response", as_index=True)["n_subjects"].sum())
        by_sex = cast(pd.Series, by_response_sex.groupby("sex", as_index=True)["n_subjects"].sum())
        n_subjects = int(by_response.sum())

    avg_b_cell = avg_b_cell_male_responders_baseline(condition, treatment, sample_type, time_filter)

    return {
        "df_raw": df,
        "by_project_samples": by_project_samples,
        "by_project_subjects": by_project_subjects,
        "by_response": by_response,
        "by_sex": by_sex,
        "n_projects": int(by_project_samples.index.nunique()),
        "n_samples": int(by_project_samples.sum()) if len(by_project_samples) > 0 else 0,
        "n_subjects": n_subjects,
        "avg_b_cell_male_responders": avg_b_cell,
    }


def get_baseline_melanoma_stats() -> dict[str, pd.Series | pd.DataFrame | int | float | None]:
    return get_subset_stats(
        condition="melanoma",
        treatment="miraclib",
        sample_type="PBMC",
        time_filter="baseline_only",
    )


def get_male_responder_b_cell_avg() -> float:
    stats = get_baseline_melanoma_stats()
    value = stats["avg_b_cell_male_responders"]
    if isinstance(value, float):
        return value
    return float("nan")


def build_cohort_flow(
    condition: str,
    treatment: str,
    sample_type: str,
    time_filter: str,
    sex: str = "all",
    response: str = "all",
) -> pd.DataFrame:
    frame = get_cell_frequency_data()

    rows: list[dict[str, int | str]] = []

    def add_step(label: str, current: pd.DataFrame) -> None:
        subject_col = "subject_pk" if "subject_pk" in current.columns else "subject_id"
        rows.append(
            {
                "step": label,
                "n_samples": int(current["sample_id"].nunique()),
                "n_subjects": int(current[subject_col].nunique()),
            }
        )

    add_step("All samples", frame)

    frame = cast(pd.DataFrame, frame.loc[frame["condition"].str.lower() == condition.lower()].copy())
    add_step(f"Condition={condition}", frame)

    frame = cast(pd.DataFrame, frame.loc[frame["sample_type"].str.lower() == sample_type.lower()].copy())
    add_step(f"SampleType={sample_type}", frame)

    frame = cast(pd.DataFrame, frame.loc[frame["treatment"].str.lower() == treatment.lower()].copy())
    add_step(f"Treatment={treatment}", frame)

    if time_filter == "baseline_only":
        frame = cast(pd.DataFrame, frame.loc[frame["visit_time"] == 0].copy())
        add_step("Time=Baseline", frame)
    else:
        add_step("Time=All", frame)

    if sex != "all":
        frame = cast(pd.DataFrame, frame.loc[frame["sex"] == sex].copy())
        add_step(f"Sex={sex}", frame)

    if response != "all":
        frame = cast(pd.DataFrame, frame.loc[frame["response"] == response.lower()].copy())
        add_step(f"Response={response}", frame)

    return pd.DataFrame(rows)

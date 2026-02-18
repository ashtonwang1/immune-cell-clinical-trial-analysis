from typing import cast

import pandas as pd

from src.analysis import get_cell_frequency_data
from src.database import get_db_connection


def _fetch_subset(
    condition: str,
    treatment: str,
    sample_type: str,
    time_filter: str,
) -> pd.DataFrame:
    conn = get_db_connection()

    query = """
    SELECT
        sub.project_id,
        s.sample_id,
        sub.subject_id,
        LOWER(sub.response) AS response,
        sub.sex,
        s.visit_time,
        c.cell_type,
        c.count
    FROM samples s
    JOIN subjects sub ON s.subject_id = sub.subject_id
    JOIN cell_counts c ON s.sample_id = c.sample_id
    WHERE LOWER(sub.condition) = ?
      AND LOWER(sub.treatment) = ?
      AND LOWER(s.sample_type) = ?
    """

    params: list[str | float] = [condition.lower(), treatment.lower(), sample_type.lower()]
    if time_filter == "baseline_only":
        query += " AND s.visit_time = ?"
        params.append(0.0)

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

    sample_unique = cast(pd.DataFrame, df.loc[:, ["project_id", "sample_id", "subject_id"]].drop_duplicates())
    subject_unique = cast(
        pd.DataFrame,
        df.loc[:, ["project_id", "subject_id", "response", "sex"]].drop_duplicates(),
    )

    by_project_samples = pd.Series(sample_unique["project_id"].tolist(), dtype="string").value_counts()
    by_project_subjects = pd.Series(subject_unique["project_id"].tolist(), dtype="string").value_counts()
    by_response = pd.Series(subject_unique["response"].tolist(), dtype="string").value_counts()
    by_sex = pd.Series(subject_unique["sex"].tolist(), dtype="string").value_counts()

    b_cell = cast(
        pd.DataFrame,
        df.loc[
            (df["sex"] == "M") & (df["response"] == "yes") & (df["cell_type"] == "b_cell"),
            ["count"],
        ].copy(),
    )
    avg_b_cell: float | None
    if len(b_cell) == 0:
        avg_b_cell = None
    else:
        avg_b_cell = float(cast(pd.Series, b_cell["count"]).mean())

    return {
        "df_raw": df,
        "by_project_samples": by_project_samples,
        "by_project_subjects": by_project_subjects,
        "by_response": by_response,
        "by_sex": by_sex,
        "n_projects": int(sample_unique["project_id"].nunique()),
        "n_samples": int(sample_unique["sample_id"].nunique()),
        "n_subjects": int(subject_unique["subject_id"].nunique()),
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
        rows.append(
            {
                "step": label,
                "n_samples": int(current["sample_id"].nunique()),
                "n_subjects": int(current["subject_id"].nunique()),
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

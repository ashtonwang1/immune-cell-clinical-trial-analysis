import pandas as pd
from typing import cast

from src.database import get_db_connection


def get_baseline_melanoma_stats() -> dict[str, pd.Series | pd.DataFrame]:
    conn = get_db_connection()
    query = """
    SELECT
        sub.project_id,
        sub.subject_id,
        sub.response,
        sub.sex,
        c.cell_type,
        c.count
    FROM samples s
    JOIN subjects sub ON s.subject_id = sub.subject_id
    JOIN cell_counts c ON s.sample_id = c.sample_id
    WHERE s.visit_time = 0
      AND sub.condition = 'melanoma'
      AND sub.treatment = 'miraclib'
      AND s.sample_type = 'PBMC'
    """

    df = cast(pd.DataFrame, pd.read_sql_query(query, conn))
    conn.close()

    unique_subjects = cast(
        pd.DataFrame, df.loc[:, ["project_id", "subject_id", "response", "sex"]].drop_duplicates()
    )
    counts_by_project = pd.Series(unique_subjects["project_id"].tolist(), dtype="string").value_counts()
    counts_by_response = pd.Series(unique_subjects["response"].tolist(), dtype="string").value_counts()
    counts_by_sex = pd.Series(unique_subjects["sex"].tolist(), dtype="string").value_counts()

    return {
        "df_raw": df,
        "by_project": counts_by_project,
        "by_response": counts_by_response,
        "by_sex": counts_by_sex,
    }


def get_male_responder_b_cell_avg() -> float:
    data = get_baseline_melanoma_stats()
    df = cast(pd.DataFrame, data["df_raw"])
    mask = (df["sex"] == "M") & (df["response"] == "yes") & (df["cell_type"] == "b_cell")
    target_cells = cast(pd.DataFrame, df.loc[mask].copy())
    avg_count = cast(pd.Series, target_cells["count"]).mean()
    return float(avg_count)

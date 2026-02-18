import sqlite3
from typing import cast

import pandas as pd

from src.config import DB_PATH


def get_db_connection() -> sqlite3.Connection:
    return sqlite3.connect(DB_PATH)


def get_cell_frequency_data() -> pd.DataFrame:
    conn = get_db_connection()

    query = """
    SELECT
        s.sample_id,
        sub.subject_id,
        sub.treatment,
        sub.response,
        sub.condition,
        s.sample_type,
        s.visit_time,
        c.cell_type,
        c.count
    FROM samples s
    JOIN subjects sub ON s.subject_id = sub.subject_id
    JOIN cell_counts c ON s.sample_id = c.sample_id
    """

    df = pd.read_sql_query(query, conn)
    conn.close()

    df["total_count"] = df.groupby("sample_id")["count"].transform("sum")
    df["percentage"] = (df["count"] / df["total_count"]) * 100

    return cast(pd.DataFrame, df)


def get_filtered_data(
    condition: str = "melanoma", treatment: str = "miraclib", sample_type: str = "PBMC"
) -> pd.DataFrame:
    df = get_cell_frequency_data()
    mask = (
        (df["condition"] == condition)
        & (df["treatment"] == treatment)
        & (df["sample_type"] == sample_type)
    )
    return cast(pd.DataFrame, df.loc[mask].copy())

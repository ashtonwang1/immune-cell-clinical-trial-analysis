import os

import pandas as pd

from src.config import CELL_TYPES, CSV_FILE
from src.database import get_db_connection, init_db


def load_csv_to_db() -> None:
    if not os.path.exists(CSV_FILE):
        print(f"Error: {CSV_FILE} not found.")
        return

    print(f"Reading data from {CSV_FILE}...")
    df = pd.read_csv(CSV_FILE)

    init_db()
    conn = get_db_connection()

    try:
        print("Processing Subjects...")
        subject_cols = {
            "subject": "subject_id",
            "project": "project_id",
            "condition": "condition",
            "age": "age",
            "sex": "sex",
            "treatment": "treatment",
            "response": "response",
        }
        subjects_df = df.loc[:, list(subject_cols.keys())].copy()
        subjects_df.columns = [subject_cols[column] for column in subjects_df.columns]
        subjects_df = subjects_df.drop_duplicates(subset=["project_id", "subject_id"])

        conn.execute("DELETE FROM cell_counts")
        conn.execute("DELETE FROM samples")
        conn.execute("DELETE FROM subjects")

        subjects_df.to_sql("subjects", conn, if_exists="append", index=False)
        print(f"-> Loaded {len(subjects_df)} subjects.")

        subject_keys = pd.read_sql_query(
            """
            SELECT subject_pk, project_id, subject_id
            FROM subjects
            """,
            conn,
        )

        print("Processing Samples...")
        sample_cols = {
            "sample": "sample_id",
            "project": "project_id",
            "subject": "subject_id",
            "time_from_treatment_start": "visit_time",
            "sample_type": "sample_type",
        }
        samples_df = df.loc[:, list(sample_cols.keys())].copy()
        samples_df.columns = [sample_cols[column] for column in samples_df.columns]
        samples_df = samples_df.drop_duplicates(subset=["sample_id"])

        samples_df = samples_df.merge(
            subject_keys,
            on=["project_id", "subject_id"],
            how="left",
            validate="many_to_one",
        )
        if samples_df["subject_pk"].isna().any():
            missing = int(samples_df["subject_pk"].isna().sum())
            raise ValueError(f"Missing subject mapping for {missing} sample rows")

        samples_df = samples_df.loc[:, ["sample_id", "subject_pk", "visit_time", "sample_type"]].copy()
        samples_df["subject_pk"] = samples_df["subject_pk"].astype(int)

        samples_df.to_sql("samples", conn, if_exists="append", index=False)
        print(f"-> Loaded {len(samples_df)} samples.")

        print("Processing Cell Counts...")
        counts_df = df.melt(
            id_vars=["sample"],
            value_vars=CELL_TYPES,
            var_name="cell_type",
            value_name="count",
        ).rename(columns={"sample": "sample_id"})

        counts_df.to_sql("cell_counts", conn, if_exists="append", index=False)
        print(f"-> Loaded {len(counts_df)} cell count records.")

        conn.commit()
        print("Data ingestion complete successfully.")

    except Exception as e:
        print(f"An error occurred: {e}")
        conn.rollback()
    finally:
        conn.close()


if __name__ == "__main__":
    load_csv_to_db()

import sqlite3

from src.config import DB_PATH


def get_db_connection() -> sqlite3.Connection:
    conn = sqlite3.connect(DB_PATH)
    conn.row_factory = sqlite3.Row
    conn.execute("PRAGMA foreign_keys = ON;")
    return conn


def init_db() -> None:
    conn = get_db_connection()
    cursor = conn.cursor()

    _ = cursor.executescript(
        """
    DROP TABLE IF EXISTS cell_counts;
    DROP TABLE IF EXISTS samples;
    DROP TABLE IF EXISTS subjects;

    CREATE TABLE subjects (
        subject_id TEXT PRIMARY KEY,
        project_id TEXT NOT NULL,
        condition TEXT NOT NULL,
        age INTEGER,
        sex TEXT CHECK (sex IN ('M', 'F')),
        treatment TEXT NOT NULL,
        response TEXT CHECK (response IN ('yes', 'no'))
    );

    CREATE TABLE samples (
        sample_id TEXT PRIMARY KEY,
        subject_id TEXT NOT NULL,
        visit_time REAL NOT NULL,
        sample_type TEXT NOT NULL,
        FOREIGN KEY (subject_id) REFERENCES subjects (subject_id) ON DELETE CASCADE
    );

    CREATE TABLE cell_counts (
        id INTEGER PRIMARY KEY AUTOINCREMENT,
        sample_id TEXT NOT NULL,
        cell_type TEXT NOT NULL,
        count INTEGER NOT NULL CHECK (count >= 0),
        FOREIGN KEY (sample_id) REFERENCES samples (sample_id) ON DELETE CASCADE,
        UNIQUE (sample_id, cell_type)
    );

    CREATE INDEX idx_subjects_project ON subjects(project_id);
    CREATE INDEX idx_subjects_condition_treatment ON subjects(condition, treatment);
    CREATE INDEX idx_subjects_response_sex ON subjects(response, sex);
    CREATE INDEX idx_samples_subject ON samples(subject_id);
    CREATE INDEX idx_samples_type_time ON samples(sample_type, visit_time);
    CREATE INDEX idx_cell_counts_sample ON cell_counts(sample_id);
    CREATE INDEX idx_cell_counts_type ON cell_counts(cell_type);
    """
    )

    conn.commit()
    conn.close()
    print(f"Database initialized at {DB_PATH}")

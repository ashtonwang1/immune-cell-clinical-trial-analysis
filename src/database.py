import sqlite3

from src.config import DB_PATH


def get_db_connection():
    conn = sqlite3.connect(DB_PATH)
    conn.row_factory = sqlite3.Row
    return conn


def init_db():
    conn = get_db_connection()
    cursor = conn.cursor()

    _ = cursor.execute(
        """
    CREATE TABLE IF NOT EXISTS subjects (
        subject_id TEXT PRIMARY KEY,
        project_id TEXT,
        condition TEXT,
        age INTEGER,
        sex TEXT,
        treatment TEXT,
        response TEXT
    );
    """
    )

    _ = cursor.execute(
        """
    CREATE TABLE IF NOT EXISTS samples (
        sample_id TEXT PRIMARY KEY,
        subject_id TEXT,
        visit_time REAL,
        sample_type TEXT,
        FOREIGN KEY (subject_id) REFERENCES subjects (subject_id)
    );
    """
    )

    _ = cursor.execute(
        """
    CREATE TABLE IF NOT EXISTS cell_counts (
        id INTEGER PRIMARY KEY AUTOINCREMENT,
        sample_id TEXT,
        cell_type TEXT,
        count INTEGER,
        FOREIGN KEY (sample_id) REFERENCES samples (sample_id)
    );
    """
    )

    conn.commit()
    conn.close()
    print(f"Database initialized at {DB_PATH}")

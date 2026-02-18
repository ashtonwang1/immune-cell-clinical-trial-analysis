# Immune Cell Clinical Trial Analysis

Production-style analytical pipeline for an immune-cell clinical trial assignment.

This project is structured as a small, maintainable analytics system rather than a single notebook/script. It separates ingestion, domain analysis, statistical testing, and presentation so each layer can evolve independently.

## Business Questions Covered

1. What is the relative frequency of each immune-cell population per sample?
2. For Melanoma + Miraclib + PBMC cohort, do responders and non-responders differ by cell population?
3. For baseline Melanoma/Miraclib/PBMC subset:
   - How many subjects are in each project?
   - What is the response/sex distribution?
   - What is the average B-cell count for male responders?

## Technical Architecture

### Layered design

- `load_data.py` is the ETL entrypoint (CSV -> normalized SQLite schema).
- `src/database.py` owns connection and schema lifecycle.
- `src/analysis.py` owns reusable feature engineering (counts -> percentages).
- `src/statistics.py` owns inferential logic (Mann-Whitney U test).
- `src/queries.py` owns targeted business queries for Part 4.
- `dashboard/app.py` is the UI layer only (no heavy business logic embedded).

This layout keeps concerns isolated, reduces coupling, and supports future extension (e.g., alternate statistical method or alternate front-end).

## Repository Structure

```text
.
├── cell-count.csv
├── load_data.py
├── run_analysis.py
├── immune_cells.db
├── requirements.txt
├── README.md
├── dashboard/
│   └── app.py
├── src/
│   ├── __init__.py
│   ├── analysis.py
│   ├── config.py
│   ├── database.py
│   ├── queries.py
│   └── statistics.py
└── tests/
    ├── __init__.py
    └── test_basic.py
```

## Data Model

SQLite database: `immune_cells.db`

- `subjects`
  - `subject_id` (PK), `project_id`, `condition`, `age`, `sex`, `treatment`, `response`
- `samples`
  - `sample_id` (PK), `subject_id` (FK), `visit_time`, `sample_type`
- `cell_counts`
  - `id` (PK), `sample_id` (FK), `cell_type`, `count`

## Why Mann-Whitney U

`src/statistics.py` uses `scipy.stats.mannwhitneyu` (two-sided) instead of a t-test because biological count/frequency data is often non-normal and can be robustly compared with a non-parametric test.

## Setup

### 1) Install dependencies

```bash
python3 -m pip install -r requirements.txt
```

### 2) Build the database from CSV

```bash
python3 load_data.py
```

Expected behavior:

- Creates/updates `immune_cells.db`
- Initializes schema
- Loads subjects, samples, and melted cell-count rows

### 3) Run command-line analysis report

```bash
python3 run_analysis.py
```

This script prints:

- Part 2 sample frequency preview
- Part 3 responder vs non-responder statistical summary
- Significant findings interpretation

### 4) Launch interactive dashboard

```bash
streamlit run dashboard/app.py
```

Dashboard tabs:

- `Statistical Analysis (Part 3)`
  - Mann-Whitney results table
  - Interactive Plotly boxplot by responder status
- `Subset Analysis (Part 4)`
  - KPI cards (projects, subjects, avg male-responder B-cell)
  - Project/response/sex distributions
  - Raw subset table explorer

## Verification

Run tests:

```bash
pytest -q
```

Current baseline from this implementation:

- Tests: `2 passed`
- ETL row counts after load:
  - `subjects`: 3500
  - `samples`: 10500
  - `cell_counts`: 52500
- Part 3 default cohort found a significant signal for `cd4_t_cell` (`p ~= 0.0133`)

## Engineering Notes

- Configuration values (paths/cell types) are centralized in `src/config.py`.
- SQL and transformation logic are explicit and reviewable.
- Core domain logic is reusable outside Streamlit (used by both CLI and UI).
- Tests validate the public analysis/statistics interfaces.

## Potential Next Improvements

1. Add parameterized SQL for all query interfaces and stricter input validation.
2. Add uniqueness guard/constraint for `(sample_id, cell_type)` to prevent accidental duplication.
3. Add subject-level aggregation mode to avoid pseudoreplication in repeated-measure settings.
4. Cache UI computation stages (`st.cache_data`) keyed by filter tuple for scale.
5. Expand test suite with integration checks for Part 4 query semantics.

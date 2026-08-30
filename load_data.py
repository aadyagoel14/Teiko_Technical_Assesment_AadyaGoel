"""
load_data.py

Initializes a SQLite database (cell-count.db) with a normalized schema and
loads all rows from cell-count.csv in the repository root.

Usage:
    python load_data.py

No command-line arguments required. Re-running this script drops and
rebuilds the database from the current contents of cell-count.csv, so it is
safe to re-run after replacing the CSV with a new export.
"""

import csv
import os
import sqlite3
import sys

CSV_PATH = os.path.join(os.path.dirname(
    os.path.abspath(__file__)), "cell-count.csv")
DB_PATH = os.path.join(os.path.dirname(
    os.path.abspath(__file__)), "cell-count.db")

# The five immune cell population columns present as wide columns in the CSV.
POPULATIONS = ["b_cell", "cd8_t_cell", "cd4_t_cell", "nk_cell", "monocyte"]

SCHEMA = """
CREATE TABLE projects (
    project_id      TEXT PRIMARY KEY
);

CREATE TABLE subjects (
    subject_id      TEXT PRIMARY KEY,
    project_id      TEXT NOT NULL REFERENCES projects(project_id),
    condition       TEXT,
    age             INTEGER,
    sex             TEXT,
    treatment       TEXT,
    response        TEXT   -- 'yes' / 'no' / NULL (e.g. untreated / healthy)
);

CREATE TABLE samples (
    sample_id                  TEXT PRIMARY KEY,
    subject_id                 TEXT NOT NULL REFERENCES subjects(subject_id),
    sample_type                TEXT,
    time_from_treatment_start  INTEGER
);

CREATE TABLE cell_counts (
    sample_id       TEXT NOT NULL REFERENCES samples(sample_id),
    population      TEXT NOT NULL,
    count           INTEGER NOT NULL,
    PRIMARY KEY (sample_id, population)
);

CREATE INDEX idx_subjects_project ON subjects(project_id);
CREATE INDEX idx_samples_subject ON samples(subject_id);
CREATE INDEX idx_cellcounts_sample ON cell_counts(sample_id);
"""

# Design notes:
# - `condition`, `age`, `sex`, `treatment`, and `response` describe the
#   subject's trial arm / outcome, not an individual blood draw, so they live
#   on `subjects` rather than being repeated on every sample row.
# - `sample_type` and `time_from_treatment_start` vary per blood draw, so
#   they live on `samples`.
# - Cell counts are stored long-format (one row per sample/population) in
#   `cell_counts` rather than as five wide columns, so adding a new
#   population later doesn't require a schema/table change, and aggregate
#   queries (totals, percentages) are simple joins/GROUP BYs.


def build_database(csv_path: str, db_path: str) -> None:
    if not os.path.exists(csv_path):
        sys.exit(f"ERROR: could not find {csv_path}")

    if os.path.exists(db_path):
        os.remove(db_path)

    conn = sqlite3.connect(db_path)
    conn.execute("PRAGMA foreign_keys = ON;")
    conn.executescript(SCHEMA)

    seen_projects = set()
    seen_subjects = set()
    n_samples = 0

    with open(csv_path, newline="") as f:
        reader = csv.DictReader(f)
        for row in reader:
            project_id = row["project"].strip()
            subject_id = row["subject"].strip()
            sample_id = row["sample"].strip()

            if project_id not in seen_projects:
                conn.execute(
                    "INSERT INTO projects (project_id) VALUES (?)", (project_id,)
                )
                seen_projects.add(project_id)

            if subject_id not in seen_subjects:
                response = row["response"].strip() or None
                conn.execute(
                    """INSERT INTO subjects
                       (subject_id, project_id, condition, age, sex, treatment, response)
                       VALUES (?, ?, ?, ?, ?, ?, ?)""",
                    (
                        subject_id,
                        project_id,
                        row["condition"].strip(),
                        int(row["age"]) if row["age"].strip() else None,
                        row["sex"].strip(),
                        row["treatment"].strip(),
                        response,
                    ),
                )
                seen_subjects.add(subject_id)

            conn.execute(
                """INSERT INTO samples
                   (sample_id, subject_id, sample_type, time_from_treatment_start)
                   VALUES (?, ?, ?, ?)""",
                (
                    sample_id,
                    subject_id,
                    row["sample_type"].strip(),
                    int(row["time_from_treatment_start"])
                    if row["time_from_treatment_start"].strip() != ""
                    else None,
                ),
            )
            n_samples += 1

            for pop in POPULATIONS:
                conn.execute(
                    "INSERT INTO cell_counts (sample_id, population, count) VALUES (?, ?, ?)",
                    (sample_id, pop, int(row[pop])),
                )

    conn.commit()
    conn.close()

    print(f"Loaded {len(seen_projects)} project(s), {len(seen_subjects)} subject(s), "
          f"{n_samples} sample(s) into {db_path}")


if __name__ == "__main__":
    build_database(CSV_PATH, DB_PATH)

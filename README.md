# Loblaw Bio : Cell Count Analysis

Analyzes immune cell population data from Bob Loblaw's miraclib clinical
trial: a relational SQLite schema, a loading script, statistical comparisons
of responders vs. non-responders, and an interactive dashboard.

**Live dashboard:** _[add your deployed Streamlit URL here — see "Dashboard" section below]_

> **Note on the provided data:** the version of `cell-count.csv` committed
> here is a partial transcription (33 rows / 11 subjects, project `prj1`
> only) used to build and validate the pipeline end-to-end before the full
> dataset was available. Drop the full dataset in as `cell-count.csv` in the
> repo root and re-run `make pipeline` — nothing else needs to change; the
> schema and every query are written generically over projects/subjects/
> samples, not hardcoded to this subset. `output/` and `cell-count.db` in
> this repo currently reflect that partial subset for the same reason.

---

## Running this project (GitHub Codespaces)

```bash
make setup       # creates a virtualenv (.venv) and installs dependencies
make pipeline    # builds cell-count.db (Part 1) and writes output/ (Parts 2-4)
make dashboard   # starts the Streamlit dashboard on port 8501
```

`make dashboard` binds to `0.0.0.0:8501`; Codespaces will detect the open
port and offer to forward/open it in the browser (or use the **Ports** tab).

No manual steps are required between commands, and `make pipeline` runs
`load_data.py` then `run_analysis.py` in sequence — one command builds
everything from a clean checkout.

Running the pieces individually (equivalent to `make setup`/`make pipeline`
without the venv) also works if you prefer:

```bash
pip install -r requirements.txt
python load_data.py          # builds cell-count.db from cell-count.csv
python run_analysis.py       # writes tables + boxplot to output/
streamlit run dashboard.py   # interactive dashboard
```

---

## Part 1: Schema design

Four normalized tables instead of one wide table:

```
projects (project_id PK)
subjects (subject_id PK, project_id FK, condition, age, sex, treatment, response)
samples  (sample_id PK, subject_id FK, sample_type, time_from_treatment_start)
cell_counts (sample_id FK, population, count)   -- PK (sample_id, population)
```

**Rationale:**

- `condition`, `age`, `sex`, `treatment`, and `response` describe the
  _subject's_ trial arm and outcome — they don't change between that
  subject's blood draws — so they're stored once per subject rather than
  repeated on every sample row (the raw CSV repeats them per row).
- `sample_type` and `time_from_treatment_start` vary per blood draw, so they
  live on `samples`, one row per draw.
- Cell counts are stored **long-format** (one row per sample × population)
  rather than as five wide columns. This means a query never needs to know
  the population names in advance — `GROUP BY population` works whether
  there are 5 populations or 50.

`load_data.py` creates `cell-count.db` with this schema and loads every row
of `cell-count.csv`. It's a plain script (`python load_data.py`, no args,
no `-m`) that lives in the repo root as required, and is idempotent — it
drops and rebuilds the database each run, so re-running after replacing the
CSV is safe.

### How this scales

**Hundreds of projects, thousands of samples:** the schema already handles
this without changes — `projects`/`subjects`/`samples` are natural
one-to-many relationships regardless of volume, and `cell_counts` grows
linearly as one row per sample × population (a full clinical dataset with
5 populations and 100,000 samples is 500,000 `cell_counts` rows — trivial
for SQLite, which handles tens of millions of rows fine). The indexes
already defined in `load_data.py` (`idx_subjects_project`,
`idx_samples_subject`, `idx_cellcounts_sample`) keep the joins used by
`compute_frequencies`/`compare_responders`/`summarize_baseline_subset`
efficient at that scale. If concurrent writes ever became a requirement
(e.g. multiple pipelines loading data simultaneously), SQLite's
single-writer model would become the limiting factor before row count did —
at that point this schema would port directly to Postgres with no
redesign, just a connection-string change.

**New/varied analytics:** the long-format `cell_counts` table is the key
design choice here. Because population is a _value_, not a column name,
adding a sixth cell population, or a completely different assay type,
is a data change (`INSERT`), not a schema migration (`ALTER TABLE`). Two
extensions that would come up in practice, and how they fit without
breaking anything already built:

- **A new assay type or measurement kind** (e.g. cytokine levels alongside
  cell counts) — add a table analogous to `cell_counts` (e.g.
  `cytokine_levels(sample_id, cytokine, level)`), keyed against the same
  `samples` table. Nothing about `subjects`/`samples`/`projects` needs to
  change, and existing queries are unaffected.
- **Analyses that need more subject-level detail** (e.g. treatment start
  date, dosage arm, prior therapies) — add columns to `subjects`, or, if a
  subject can have several distinct treatment courses over time, split
  `treatment`/`response` out into their own `treatment_courses` table
  keyed by `subject_id` with a start/end date range, and have `samples`
  reference the course active at draw time rather than the subject
  directly. That's the one place the current schema simplifies (assuming
  one treatment course per subject) that a longer-running, multi-arm trial
  would eventually need to relax.

## Part 2: Frequency table

`src/analysis.compute_frequencies()` joins `samples` and `cell_counts`,
sums counts per sample for the denominator, and returns one row per
sample × population with `sample, total_count, population, count,
percentage`. Surfaced in `run_analysis.py` (→ `output/frequencies.csv`) and
in the dashboard's first tab.

## Part 3: Responders vs. non-responders

Restricted to melanoma subjects treated with miraclib, PBMC samples only
(`src/analysis.get_responder_comparison_data`). For each of the five
populations:

- **Boxplot** (responders vs. non-responders), rendered as a static PNG via
  `run_analysis.py` and interactively (with individual points) in the
  dashboard.
- **Mann-Whitney U test per population**, run on **one aggregated value per
  subject** (mean % across that subject's samples/timepoints) — chosen over
  a t-test because trial arms here are small-n and percentage data isn't
  guaranteed to be normal, and aggregated to the subject level because each
  subject contributes multiple samples (baseline/day 7/day 14) that are not
  independent observations; testing at the raw-sample level would inflate
  the effective n and distort the p-values.
- **Benjamini-Hochberg FDR correction** across the 5 simultaneous tests, so
  the reported significance isn't inflated by multiple comparisons.

`output/responder_stats.csv` (or the dashboard's Part 3 tab) shows which
populations pass the FDR < 0.05 threshold — that's the evidence to bring to
Yah.

## Part 4: Baseline subset

`src/analysis.summarize_baseline_subset()` filters to melanoma, PBMC,
miraclib, `time_from_treatment_start = 0`, then breaks the resulting
subjects/samples down by project, response, and sex. Available via
`run_analysis.py` (`output/baseline_*.csv`) and the dashboard's Part 4 tab.

---

## Code structure

```
load_data.py         # Part 1 — builds cell-count.db from cell-count.csv
run_analysis.py       # CLI — Parts 2-4, writes output/
dashboard.py           # Streamlit app — Parts 2-4, interactive, same logic as above
src/analysis.py        # shared query/stats functions imported by both entry points
cell-count.csv          # input data
cell-count.db            # generated by load_data.py (committed for convenience/grading)
output/                   # generated by run_analysis.py (committed for convenience/grading)
Makefile
requirements.txt
.streamlit/config.toml     # light theme for the dashboard
```

**Why split this way:** `src/analysis.py` holds every query and statistical
function, with no I/O side effects beyond reading from the DB connection
it's given. `load_data.py` is the only file that writes to the database.
`run_analysis.py` and `dashboard.py` are two different _views_ onto the
same analysis code — a non-interactive CLI that writes files for
grading/scripting, and an interactive UI for exploration — rather than two
separate implementations. That way a change to, say, how the Mann-Whitney
test is computed only needs to happen in one place, and the CLI output and
the dashboard can never silently disagree with each other.

## Dashboard

`make dashboard` (or `streamlit run dashboard.py`) — three tabs (Parts
2-4), each reading live from `cell-count.db`:

- **Part 2**: filterable frequency table + stacked composition chart.
- **Part 3**: interactive boxplot, the Mann-Whitney/BH results table, and a
  plain-language summary of which populations are significant.
- **Part 4**: the baseline subset plus the three breakdown tables.

To get a persistent shareable link (rather than a Codespaces port-forward
URL that only works while the Codespace is running), deploy to
[Streamlit Community Cloud](https://share.streamlit.io) — connect this
repo, set the main file to `dashboard.py`, and it will build from
`requirements.txt` automatically. Add the resulting URL to the top of this
README.

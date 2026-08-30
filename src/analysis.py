"""
analysis.py

Query and analysis functions used for Parts 2-4 of the assignment.
Kept separate from load_data.py and the dashboard so both a CLI script and
the Streamlit app can import the same logic.
"""

import sqlite3
from itertools import combinations

import pandas as pd
from scipy import stats


def get_connection(db_path: str) -> sqlite3.Connection:
    # check_same_thread=False: Streamlit's @st.cache_resource can hand this
    # connection to a script rerun happening on a different internal
    # thread. This app only ever reads (never writes concurrently), so
    # sharing one connection across threads is safe here.
    return sqlite3.connect(db_path, check_same_thread=False)


# Part 2: relative frequency table

def compute_frequencies(conn: sqlite3.Connection) -> pd.DataFrame:
    """
    Returns one row per (sample, population) with columns:
    sample, total_count, population, count, percentage
    """
    df = pd.read_sql_query(
        """
        SELECT s.sample_id AS sample, cc.population, cc.count
        FROM cell_counts cc
        JOIN samples s ON s.sample_id = cc.sample_id
        """,
        conn,
    )
    totals = df.groupby("sample")["count"].transform("sum")
    df["total_count"] = totals
    df["percentage"] = (df["count"] / df["total_count"] * 100).round(4)
    return df[["sample", "total_count", "population", "count", "percentage"]].sort_values(
        ["sample", "population"]
    ).reset_index(drop=True)


# Part 3: responders vs non-responders (melanoma, miraclib, PBMC)

def get_responder_comparison_data(conn: sqlite3.Connection) -> pd.DataFrame:
    """
    Relative-frequency table restricted to melanoma patients treated with
    miraclib, PBMC samples only, joined with each subject's response.
    """
    freq = compute_frequencies(conn)

    meta = pd.read_sql_query(
        """
        SELECT s.sample_id AS sample, s.sample_type, sub.subject_id,
               sub.condition, sub.treatment, sub.response
        FROM samples s
        JOIN subjects sub ON sub.subject_id = s.subject_id
        WHERE sub.condition = 'melanoma'
          AND sub.treatment = 'miraclib'
          AND s.sample_type = 'PBMC'
        """,
        conn,
    )

    merged = freq.merge(meta, on="sample", how="inner")
    return merged


def compare_responders(conn: sqlite3.Connection) -> pd.DataFrame:
    """
    For each population, runs a Mann-Whitney U test comparing relative
    frequency between responders and non-responders. Returns a table with
    raw and Benjamini-Hochberg (FDR) adjusted p-values, sorted by p-value.

    Mann-Whitney U (rather than a t-test) is used because clinical trial
    arms like this are typically small-n and percentage data is not
    guaranteed to be normally distributed; Mann-Whitney makes no such
    assumption. Because five populations are tested simultaneously, a
    Benjamini-Hochberg correction is applied to control the false discovery
    rate rather than relying on the raw p-values alone.

    IMPORTANT: subjects contribute multiple samples (one per
    time_from_treatment_start), and those repeated measurements are not
    independent observations. To avoid pseudoreplication inflating the
    apparent sample size (and distorting the test), each subject's samples
    are first collapsed to a single mean relative-frequency value per
    population before the test is run — so n reflects subjects, not samples.
    """
    data = get_responder_comparison_data(conn)
    populations = sorted(data["population"].unique())

    # Collapse repeated measures per subject to one value per
    # subject x population (mean across their samples), so each subject
    # contributes exactly one independent observation to the test.
    subject_level = (
        data.groupby(["subject_id", "population", "response"])["percentage"]
        .mean()
        .reset_index()
    )

    results = []
    for pop in populations:
        sub = subject_level[subject_level["population"] == pop]
        resp = sub.loc[sub["response"] == "yes", "percentage"]
        nonresp = sub.loc[sub["response"] == "no", "percentage"]
        if len(resp) == 0 or len(nonresp) == 0:
            continue
        u_stat, p_value = stats.mannwhitneyu(
            resp, nonresp, alternative="two-sided")
        results.append(
            {
                "population": pop,
                "n_responders": len(resp),
                "n_non_responders": len(nonresp),
                "median_responders_pct": round(resp.median(), 3),
                "median_non_responders_pct": round(nonresp.median(), 3),
                "u_statistic": u_stat,
                "p_value": p_value,
            }
        )

    result_df = pd.DataFrame(results).sort_values(
        "p_value").reset_index(drop=True)

    # Benjamini-Hochberg FDR correction
    m = len(result_df)
    result_df["p_value_rank"] = result_df["p_value"].rank(method="first")
    result_df["p_adj_bh"] = (
        result_df["p_value"] * m / result_df["p_value_rank"]
    ).clip(upper=1.0)
    # Enforce monotonicity of BH-adjusted p-values
    result_df["p_adj_bh"] = result_df["p_adj_bh"][::-1].cummin()[::-1]
    result_df["significant_fdr_0.05"] = result_df["p_adj_bh"] < 0.05
    result_df = result_df.drop(columns=["p_value_rank"])

    return result_df


# Part 4: baseline melanoma / miraclib / PBMC subset

def get_baseline_subset(conn: sqlite3.Connection) -> pd.DataFrame:
    """
    All melanoma, PBMC samples at baseline (time_from_treatment_start = 0)
    from subjects treated with miraclib. One row per sample/subject.
    """
    return pd.read_sql_query(
        """
        SELECT s.sample_id AS sample, s.time_from_treatment_start,
               sub.subject_id, sub.project_id, sub.condition, sub.treatment,
               sub.response, sub.sex, sub.age
        FROM samples s
        JOIN subjects sub ON sub.subject_id = s.subject_id
        WHERE sub.condition = 'melanoma'
          AND sub.treatment = 'miraclib'
          AND s.sample_type = 'PBMC'
          AND s.time_from_treatment_start = 0
        """,
        conn,
    )


def summarize_baseline_subset(conn: sqlite3.Connection) -> dict:
    subset = get_baseline_subset(conn)

    by_project = (
        subset.groupby("project_id")["sample"].nunique()
        .rename("n_samples").reset_index()
    )
    by_response = (
        subset.drop_duplicates("subject_id").groupby("response")["subject_id"]
        .nunique().rename("n_subjects").reset_index()
    )
    by_sex = (
        subset.drop_duplicates("subject_id").groupby("sex")["subject_id"]
        .nunique().rename("n_subjects").reset_index()
    )

    return {
        "subset": subset,
        "by_project": by_project,
        "by_response": by_response,
        "by_sex": by_sex,
    }

"""
run_analysis.py

Runs Parts 2-4 of the analysis against cell-count.db and writes results to
the output/ directory:

  output/frequencies.csv          Part 2: per-sample population frequencies
  output/responder_stats.csv      Part 3: responder vs non-responder stats
  output/responder_boxplot.png    Part 3: boxplot of frequencies by response
  output/baseline_subset.csv      Part 4: baseline melanoma/miraclib/PBMC subset
  output/baseline_by_project.csv  Part 4: sample counts by project
  output/baseline_by_response.csv Part 4: subject counts by response
  output/baseline_by_sex.csv      Part 4: subject counts by sex

Run `python load_data.py` first to (re)build cell-count.db.

Usage:
    python run_analysis.py
"""

import seaborn as sns
import matplotlib.pyplot as plt
import os
import sys

import matplotlib
matplotlib.use("Agg")

sys.path.insert(0, os.path.join(
    os.path.dirname(os.path.abspath(__file__)), "src"))
from analysis import (  # noqa: E402
    get_connection,
    compute_frequencies,
    get_responder_comparison_data,
    compare_responders,
    summarize_baseline_subset,
)

ROOT = os.path.dirname(os.path.abspath(__file__))
DB_PATH = os.path.join(ROOT, "cell-count.db")
OUT_DIR = os.path.join(ROOT, "output")


def main():
    if not os.path.exists(DB_PATH):
        sys.exit("ERROR: cell-count.db not found. Run `python load_data.py` first.")

    os.makedirs(OUT_DIR, exist_ok=True)
    conn = get_connection(DB_PATH)

    # Part 2
    freq = compute_frequencies(conn)
    freq.to_csv(os.path.join(OUT_DIR, "frequencies.csv"), index=False)
    print(f"[Part 2] wrote {len(freq)} rows -> output/frequencies.csv")

    # Part 3
    comparison_data = get_responder_comparison_data(conn)
    stats_df = compare_responders(conn)
    stats_df.to_csv(os.path.join(OUT_DIR, "responder_stats.csv"), index=False)
    print("[Part 3] responder vs non-responder statistics:")
    print(stats_df.to_string(index=False))

    if not comparison_data.empty:
        # Aggregate to one value per subject x population (mean across that
        # subject's samples) so the plot matches what compare_responders()
        # actually tests -- one independent point per subject, not per
        # sample/timepoint.
        subject_level = (
            comparison_data.groupby(["subject_id", "population", "response"])[
                "percentage"]
            .mean()
            .reset_index()
        )

        plt.figure(figsize=(10, 6))
        order = sorted(subject_level["population"].unique())
        sns.boxplot(
            data=subject_level,
            x="population",
            y="percentage",
            hue="response",
            order=order,
            hue_order=["yes", "no"],
        )
        sns.stripplot(
            data=subject_level,
            x="population",
            y="percentage",
            hue="response",
            order=order,
            hue_order=["yes", "no"],
            dodge=True,
            palette=["black", "black"],
            size=4,
            alpha=0.6,
            legend=False,
        )
        plt.title("Relative frequency by population: responders vs non-responders\n"
                  "(melanoma, miraclib, PBMC — one point per subject, averaged across timepoints)")
        plt.ylabel("Relative frequency (%)")
        plt.xlabel("Cell population")
        plt.legend(title="Response")
        plt.tight_layout()
        plot_path = os.path.join(OUT_DIR, "responder_boxplot.png")
        plt.savefig(plot_path, dpi=150)
        print(f"[Part 3] wrote boxplot -> output/responder_boxplot.png")
    else:
        print(
            "[Part 3] no matching melanoma/miraclib/PBMC samples found; skipping boxplot")

    # Part 4
    baseline = summarize_baseline_subset(conn)
    baseline["subset"].to_csv(os.path.join(
        OUT_DIR, "baseline_subset.csv"), index=False)
    baseline["by_project"].to_csv(os.path.join(
        OUT_DIR, "baseline_by_project.csv"), index=False)
    baseline["by_response"].to_csv(os.path.join(
        OUT_DIR, "baseline_by_response.csv"), index=False)
    baseline["by_sex"].to_csv(os.path.join(
        OUT_DIR, "baseline_by_sex.csv"), index=False)
    print(
        f"[Part 4] {len(baseline['subset'])} baseline samples -> output/baseline_*.csv")

    conn.close()


if __name__ == "__main__":
    main()

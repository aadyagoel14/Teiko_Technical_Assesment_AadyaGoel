"""
dashboard.py

Interactive Streamlit dashboard for Bob's cell-count analysis.

Usage:
    python load_data.py      # build/refresh cell-count.db (run once, or after updating the CSV)
    streamlit run dashboard.py
"""

import importlib.util
import os

import pandas as pd
import plotly.express as px
import streamlit as st

ROOT = os.path.dirname(os.path.abspath(__file__))
DB_PATH = os.path.join(ROOT, "cell-count.db")

# Load src/analysis.py by file path rather than `sys.path.insert` + a normal
# `from analysis import ...`. A normal import is a static import statement
# that editor/linter "organize imports" features will happily move above
# the sys.path setup, breaking it. Loading by path this way isn't a
# reorderable import statement, so it survives auto-formatting.
_analysis_spec = importlib.util.spec_from_file_location(
    "analysis", os.path.join(ROOT, "src", "analysis.py")
)
analysis = importlib.util.module_from_spec(_analysis_spec)
_analysis_spec.loader.exec_module(analysis)

get_connection = analysis.get_connection
compute_frequencies = analysis.compute_frequencies
get_responder_comparison_data = analysis.get_responder_comparison_data
compare_responders = analysis.compare_responders
summarize_baseline_subset = analysis.summarize_baseline_subset

st.set_page_config(
    page_title="Loblaw Bio - Cell Count Dashboard", layout="wide")


@st.cache_resource
def _connect():
    return get_connection(DB_PATH)


def main():
    st.title("Immune Cell Population Dashboard")
    st.caption("Loblaw Bio : miraclib clinical trial cell-count analysis")

    if not os.path.exists(DB_PATH):
        st.error(
            "cell-count.db not found. Run `python load_data.py` in the repo root first, "
            "then reload this page."
        )
        st.stop()

    conn = _connect()

    tab2, tab3, tab4 = st.tabs(
        ["Frequencies", "Responders vs Non-responders",
            "Baseline Subset"]
    )

    # Part 2
    with tab2:
        st.subheader("Relative frequency of each cell population, per sample")
        freq = compute_frequencies(conn)

        samples = sorted(freq["sample"].unique())
        selected = st.multiselect(
            "Filter by sample (leave empty to show all)", samples)
        view = freq[freq["sample"].isin(selected)] if selected else freq

        st.dataframe(view, use_container_width=True, hide_index=True)
        st.download_button(
            "Download full table as CSV",
            freq.to_csv(index=False),
            file_name="frequencies.csv",
            mime="text/csv",
        )

        st.markdown("**Population composition by sample**")
        fig = px.bar(
            view, x="sample", y="percentage", color="population",
            title="Relative frequency (%) by sample",
        )
        fig.update_layout(xaxis_tickangle=-45)
        st.plotly_chart(fig, use_container_width=True)

    # Part 3
    with tab3:
        st.subheader(
            "Melanoma patients on miraclib (PBMC samples): responders vs non-responders")

        comparison_data = get_responder_comparison_data(conn)

        if comparison_data.empty:
            st.warning(
                "No samples match melanoma + miraclib + PBMC in the current database. "
                "Load the full dataset to see this analysis."
            )
        else:
            subject_level = (
                comparison_data.groupby(["subject_id", "population", "response"])[
                    "percentage"]
                .mean()
                .reset_index()
            )
            fig = px.box(
                subject_level,
                x="population",
                y="percentage",
                color="response",
                points="all",
                category_orders={"response": ["yes", "no"]},
                labels={
                    "percentage": "Relative frequency (%)", "population": "Cell population"},
                title="Relative frequency by population: responders (yes) vs non-responders (no)<br>"
                      "<sup>one point per subject, averaged across their timepoints</sup>",
            )
            st.plotly_chart(fig, use_container_width=True)

            st.markdown(
                "**Statistical comparison (Mann–Whitney U, BH-adjusted)**")
            stats_df = compare_responders(conn)
            st.dataframe(stats_df, use_container_width=True, hide_index=True)

            sig = stats_df[stats_df["significant_fdr_0.05"]]
            if not sig.empty:
                pops = ", ".join(sig["population"])
                st.success(
                    f"Statistically significant difference (FDR < 0.05) between responders and "
                    f"non-responders in: **{pops}**."
                )
            else:
                st.info(
                    "No population reached statistical significance at FDR < 0.05 in the current data.")

            st.caption(
                "Each subject contributes one point per population (mean % across their "
                "samples/timepoints), avoiding pseudoreplication from repeated measures. "
                "Mann-Whitney U test used per population (no normality assumption, appropriate "
                "for small clinical n); Benjamini-Hochberg correction applied across the 5 "
                "populations tested simultaneously to control the false discovery rate."
            )

    # Part 4
    with tab4:
        st.subheader("Baseline (day 0) melanoma, PBMC, miraclib samples")

        baseline = summarize_baseline_subset(conn)
        subset = baseline["subset"]

        if subset.empty:
            st.warning(
                "No baseline (time_from_treatment_start = 0) samples match the filters in the current data.")
        else:
            st.dataframe(subset, use_container_width=True, hide_index=True)

            c1, c2, c3 = st.columns(3)
            with c1:
                st.markdown("**Samples by project**")
                st.dataframe(baseline["by_project"], hide_index=True)
            with c2:
                st.markdown("**Subjects by response**")
                st.dataframe(baseline["by_response"], hide_index=True)
            with c3:
                st.markdown("**Subjects by sex**")
                st.dataframe(baseline["by_sex"], hide_index=True)


if __name__ == "__main__":
    main()

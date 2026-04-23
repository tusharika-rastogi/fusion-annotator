"""
app.py -- EML4-ALK Fusion Annotator Streamlit UI.

All annotation logic is in annotator.py. This file handles UI only.
"""

import io

import pandas as pd
import streamlit as st

from annotator import annotate_df, load_exon_reference

st.set_page_config(page_title="Fusion Annotator", layout="centered")
st.title("Fusion Annotator")

# ── Genome build selector ─────────────────────────────────────────────────────
build = st.radio(
    "Genome build",
    options=["hg19", "hg38"],
    index=0,
    horizontal=True,
    help="Select the genome build that matches your breakpoint coordinates.",
)


@st.cache_data
def get_exon_df(genome_build: str) -> pd.DataFrame:
    """Load and cache the exon reference DataFrame for the given build."""
    return load_exon_reference(genome_build)


try:
    exon_df = get_exon_df(build)
except FileNotFoundError as e:
    st.error(str(e))
    st.stop()


def show_results(annotated: pd.DataFrame, warnings: list[str], filename: str) -> None:
    """Render warnings, annotated table, frequency table, and download button."""
    if warnings:
        st.warning("\n".join(warnings))

    st.subheader("Annotated results")
    st.dataframe(annotated, use_container_width=True)

    st.subheader("Variant type counts")
    freq = annotated["EML4-ALK_VariantType"].value_counts().reset_index()
    if "count" in freq.columns:
        freq = freq.rename(columns={"EML4-ALK_VariantType": "Variant type", "count": "Count"})
    else:
        freq = freq.rename(columns={"index": "Variant type", "EML4-ALK_VariantType": "Count"})
    st.dataframe(freq, use_container_width=True)

    csv_bytes = annotated.to_csv(index=False).encode("utf-8")
    st.download_button(
        label="Download annotated CSV",
        data=csv_bytes,
        file_name=filename,
        mime="text/csv",
    )


tab_single, tab_paste, tab_upload = st.tabs(["Single entry", "Paste rows", "Upload CSV"])

# ── Tab A: Single fusion manual entry ────────────────────────────────────────
with tab_single:
    st.markdown("Enter one fusion manually.")
    sample_id = st.text_input("Sample ID (optional)", key="single_sample")
    fusion_name = st.text_input("Fusion name", value="EML4-ALK", key="single_name")
    bp_a = st.text_input("Breakpoint A (e.g. 2_29446394)", key="single_bpa")
    bp_b = st.text_input("Breakpoint B (e.g. 2_42522656)", key="single_bpb")

    if st.button("Annotate", key="btn_single"):
        if not fusion_name or not bp_a or not bp_b:
            st.error("Fusion name, Breakpoint A, and Breakpoint B are required.")
        else:
            try:
                df_single = pd.DataFrame([{
                    "sample_id": sample_id or "unknown",
                    "fusion_name": fusion_name,
                    "bp_a": bp_a,
                    "bp_b": bp_b,
                }])
                col_map = {
                    "sample_id": "sample_id",
                    "fusion_name": "fusion_name",
                    "bp_a": "bp_a",
                    "bp_b": "bp_b",
                }
                annotated, warnings = annotate_df(df_single, col_map, exon_df)
                show_results(annotated, warnings, "fusions_annotated.csv")
            except Exception as exc:
                st.error(f"Annotation failed: {exc}")

# ── Tab B: Paste multiple fusions ────────────────────────────────────────────
with tab_paste:
    st.markdown(
        "Paste tab- or comma-separated rows including a header row. "
        "Column order: **Sample ID, Fusion name, Breakpoint A, Breakpoint B**."
    )
    pasted = st.text_area(
        "Paste data here",
        height=200,
        placeholder=(
            "Sample\tFusion\tBreakpointA\tBreakpointB\n"
            "PT001\tEML4-ALK\t2_29446394\t2_42522656"
        ),
        key="paste_area",
    )

    if pasted.strip():
        try:
            sep = "\t" if "\t" in pasted else ","
            df_paste = pd.read_csv(io.StringIO(pasted), sep=sep, dtype=str)
            if df_paste.shape[1] < 4:
                st.error("Need at least 4 columns. Check your separator (tab or comma).")
            else:
                cols = df_paste.columns.tolist()
                col_map_paste = {
                    "sample_id":   cols[0],
                    "fusion_name": cols[1],
                    "bp_a":        cols[2],
                    "bp_b":        cols[3],
                }
                st.markdown("**Preview (first 5 rows)**")
                st.dataframe(df_paste.head(5), use_container_width=True)

                if st.button("Annotate", key="btn_paste"):
                    try:
                        annotated, warnings = annotate_df(df_paste, col_map_paste, exon_df)
                        show_results(annotated, warnings, "fusions_annotated.csv")
                    except Exception as exc:
                        st.error(f"Annotation failed: {exc}")
        except Exception as exc:
            st.error(f"Could not parse pasted data: {exc}")

# ── Tab C: Upload CSV ─────────────────────────────────────────────────────────
with tab_upload:
    st.markdown(
        "Upload a CSV file. Column order (names may vary): "
        "**1st = Sample ID, 2nd = Fusion name, 3rd = Breakpoint A, 4th = Breakpoint B**. "
        "A 5-row preview is shown before annotating."
    )
    uploaded = st.file_uploader("Choose a CSV file", type=["csv"], key="upload_file")

    if uploaded is not None:
        try:
            df_upload = pd.read_csv(uploaded, dtype=str)
            if df_upload.shape[1] < 4:
                st.error("File must have at least 4 columns.")
            else:
                cols = df_upload.columns.tolist()
                col_map_upload = {
                    "sample_id":   cols[0],
                    "fusion_name": cols[1],
                    "bp_a":        cols[2],
                    "bp_b":        cols[3],
                }

                st.markdown("**Preview (first 5 rows)**")
                st.dataframe(df_upload.head(5), use_container_width=True)

                out_filename = uploaded.name.replace(".csv", "_annotated.csv")

                if st.button("Annotate", key="btn_upload"):
                    try:
                        annotated, warnings = annotate_df(df_upload, col_map_upload, exon_df)
                        show_results(annotated, warnings, out_filename)
                    except Exception as exc:
                        st.error(f"Annotation failed: {exc}")
        except Exception as exc:
            st.error(f"Could not read file: {exc}")

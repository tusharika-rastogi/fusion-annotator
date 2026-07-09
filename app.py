"""
app.py -- EML4-ALK Fusion Annotator Streamlit UI.

All annotation logic is in annotator.py. This file handles UI only.
"""

import io

import pandas as pd
import streamlit as st

from annotator import annotate_df, load_exon_reference

APP_VERSION = "v1.4.0"

_VARIANT_CARD: dict[str, tuple[str, str, str]] = {
    "V1":                       ("13",                       "20",   "Most common (~33%); V6 (E13;ins69A20) also maps here"),
    "V2":                       ("20",                       "20",   "~9%"),
    "V3a/b":                    ("6",                        "20",   "~29%; RNA or long-read required to distinguish a/b"),
    "V4":                       ("15",                       "20",   "Rare (~2%); deletions at junction"),
    "V4'":                      ("14",                       "20",   "Rare (<1%); ALK bp ~49bp into intron 19 (ins11del49 at junction)"),
    "V7":                       ("14",                       "20",   "Rare (<1%); ALK bp ~12bp into intron 19 (del12 at junction)"),
    "V4'/V7":                   ("14",                       "20",   "Rare (<1%); ambiguous — ALK bp distance falls between V7 and V4' thresholds"),
    "V5a/b":                    ("2",                        "20",   "Rare (~2%); sub-variant requires junction sequencing"),
    "V5'":                      ("18",                       "20",   "Rare (~2%)"),
    "V8a/b":                    ("17",                       "20",   "Very rare (<1%); two sub-variants differ in functional ALK domain"),
    "EML4-ALK_NonCanonicalALK": ("any",                      "≠ 20", "Rare; verify caller output"),
    "EML4-ALK_OtherVariant":    ("≠ 2,6,13,14,15,17,18,20", "20",   "Novel or uncharacterized EML4 exon"),
}

st.set_page_config(page_title="Fusion Annotator", layout="wide")

st.markdown("""
<style>
.left-header {
    background-color: #0D7377;
    padding: 1.2rem 1.5rem;
    border-radius: 0.5rem;
    margin-bottom: 1rem;
}
.left-header h1 { color: white; font-size: 1.6rem; margin: 0; font-weight: 700; }
.left-header p  { color: #A8D5D7; font-size: 0.85rem; margin: 0.25rem 0 0 0; }
</style>
""", unsafe_allow_html=True)

if "results" not in st.session_state:
    st.session_state.results = None
if "last_build" not in st.session_state:
    st.session_state.last_build = "hg19"


@st.cache_data
def get_exon_df(genome_build: str) -> pd.DataFrame:
    return load_exon_reference(genome_build)


def _variant_ref_df() -> pd.DataFrame:
    return pd.DataFrame([
        {"Variant": k, "5′ Exon": v[0], "3′ Exon": v[1], "Notes": v[2]}
        for k, v in _VARIANT_CARD.items()
    ])


def show_results(annotated: pd.DataFrame, warnings: list[str], filename: str) -> None:
    if warnings:
        st.warning("\n".join(warnings))

    st.subheader("Annotated results")
    st.dataframe(annotated, use_container_width=True, hide_index=True)

    st.subheader("Variant type counts")
    freq = annotated["EML4-ALK_VariantType"].value_counts().reset_index()
    if "count" in freq.columns:
        freq = freq.rename(columns={"EML4-ALK_VariantType": "Variant type", "count": "Count"})
    else:
        freq = freq.rename(columns={"index": "Variant type", "EML4-ALK_VariantType": "Count"})
    st.dataframe(freq, use_container_width=True, hide_index=True)

    csv_bytes = annotated.to_csv(index=False).encode("utf-8")
    st.download_button(
        label="Download annotated CSV",
        data=csv_bytes,
        file_name=filename,
        mime="text/csv",
    )


# ── Two-panel layout ──────────────────────────────────────────────────────────
col_left, col_right = st.columns([5, 6], gap="large")

with col_left:
    st.markdown(
        f'<div class="left-header"><h1>Fusion Annotator</h1>'
        f'<p>EML4-ALK variant classifier &nbsp;|&nbsp; {APP_VERSION}</p></div>',
        unsafe_allow_html=True,
    )

    build = st.radio(
        "Genome build",
        options=["hg19", "hg38"],
        index=0,
        horizontal=True,
        key="build_radio",
        help="Select the genome build that matches your breakpoint coordinates.",
    )

    if build != st.session_state.last_build:
        st.session_state.results = None
        st.session_state.last_build = build

    try:
        exon_df = get_exon_df(build)
    except FileNotFoundError as e:
        st.error(str(e))
        st.stop()

    tab_single, tab_paste, tab_upload = st.tabs(["Single entry", "Paste rows", "Upload CSV"])

    # ── Tab A: Single fusion manual entry ─────────────────────────────────────
    with tab_single:
        # hg19 exon-midpoint coords for each canonical variant
        _EXAMPLES = {
            "V1 (exon 13–20)": ("EXAMPLE", "EML4-ALK", "2_42522588", "2_29446301"),
            "V2 (exon 20–20)": ("EXAMPLE", "EML4-ALK", "2_42552650", "2_29446301"),
            "V3a/b (exon 6–20)": ("EXAMPLE", "EML4-ALK", "2_42491858", "2_29446301"),
        }
        with st.expander("Try an example (hg19)"):
            ex_cols = st.columns(len(_EXAMPLES))
            for col, (label, vals) in zip(ex_cols, _EXAMPLES.items()):
                if col.button(label, use_container_width=True):
                    (st.session_state["single_sample"],
                     st.session_state["single_name"],
                     st.session_state["single_bpa"],
                     st.session_state["single_bpb"]) = vals

        sample_id  = st.text_input("Sample ID (optional)", key="single_sample")
        fusion_name = st.text_input("Fusion name", value="EML4-ALK", key="single_name")
        bp_a       = st.text_input("Breakpoint A (e.g. 2_29446394)", key="single_bpa")
        bp_b       = st.text_input("Breakpoint B (e.g. 2_42522656)", key="single_bpb")

        if st.button("Annotate Fusion", key="btn_single", type="primary", use_container_width=True):
            if not fusion_name or not bp_a or not bp_b:
                st.error("Fusion name, Breakpoint A, and Breakpoint B are required.")
            else:
                try:
                    df_single = pd.DataFrame([{
                        "sample_id":   sample_id or "unknown",
                        "fusion_name": fusion_name,
                        "bp_a":        bp_a,
                        "bp_b":        bp_b,
                    }])
                    col_map = {
                        "sample_id":   "sample_id",
                        "fusion_name": "fusion_name",
                        "bp_a":        "bp_a",
                        "bp_b":        "bp_b",
                    }
                    annotated, warnings = annotate_df(df_single, col_map, exon_df)
                    st.session_state.results = (annotated, warnings, f"fusions_annotated_{build}.csv")
                except Exception as exc:
                    st.error(f"Annotation failed: {exc}")

    # ── Tab B: Paste multiple fusions ─────────────────────────────────────────
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

                    if st.button("Annotate Fusion", key="btn_paste", type="primary", use_container_width=True):
                        try:
                            annotated, warnings = annotate_df(df_paste, col_map_paste, exon_df)
                            st.session_state.results = (annotated, warnings, f"fusions_annotated_{build}.csv")
                        except Exception as exc:
                            st.error(f"Annotation failed: {exc}")
            except Exception as exc:
                st.error(f"Could not parse pasted data: {exc}")

    # ── Tab C: Upload CSV ──────────────────────────────────────────────────────
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

                    base = uploaded.name.removesuffix(".csv")
                    out_filename = f"{base}_annotated_{build}.csv"

                    if st.button("Annotate Fusion", key="btn_upload", type="primary", use_container_width=True):
                        try:
                            annotated, warnings = annotate_df(df_upload, col_map_upload, exon_df)
                            st.session_state.results = (annotated, warnings, out_filename)
                        except Exception as exc:
                            st.error(f"Annotation failed: {exc}")
            except Exception as exc:
                st.error(f"Could not read file: {exc}")

with col_right:
    st.subheader("EML4-ALK Variant Reference")
    st.dataframe(_variant_ref_df(), use_container_width=True, hide_index=True)

    st.divider()

    if st.session_state.results is not None:
        annotated, warnings, filename = st.session_state.results
        show_results(annotated, warnings, filename)
    else:
        st.caption("After annotation, results appear here.")

    with st.expander("Documentation"):
        st.markdown(f"**Fusion Annotator {APP_VERSION}**")
        st.markdown("Full documentation, source code, and issue tracking are on GitHub.")
        st.link_button(
            "Open documentation on GitHub",
            "https://github.com/tusharika-rastogi/fusion-annotator",
        )

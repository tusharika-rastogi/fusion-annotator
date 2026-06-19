"""
annotator.py -- EML4-ALK fusion variant classification logic.

All annotation logic lives here. No Streamlit imports allowed.
"""

from pathlib import Path

import pandas as pd

REFERENCE_DIR = Path(__file__).parent / "reference"

REFERENCE_FILES: dict[str, Path] = {
    "hg38": REFERENCE_DIR / "eml4_alk_exons_hg38.csv",
    "hg19": REFERENCE_DIR / "eml4_alk_exons_hg19.csv",
}

# Variant rules: (eml4_exon, alk_exon) -> label
# Junction-level insertions/deletions (V4 del60, V5b ins117, V6 ins69, V8a/b ins30/ins61)
# are not detectable from coarse breakpoint coordinates; those sub-variants are collapsed.
# V4'/V7 are both E14:A20 but differ in how far the ALK breakpoint sits into intron 19
# (V7 ~12bp, V4' ~49bp); sub-classified at runtime via _classify_v4prime_v7().
_VARIANT_RULES: dict[tuple[int, int], str] = {
    (13, 20): "V1",       # V6 (E13;ins69A20) also maps here; distinguish by junction seq
    (20, 20): "V2",
    (6,  20): "V3a/b",
    (15, 20): "V4",
    (14, 20): "V4'/V7",  # refined to V4' or V7 by intron distance at runtime
    (2,  20): "V5a/b",
    (18, 20): "V5'",
    (17, 20): "V8a/b",
}

# Distance thresholds (bp from ALK exon 20 boundary) for V4'/V7 sub-classification.
# V7: E14;del12A20 → breakpoint ~12bp into intron 19
# V4': E14;ins11del49A20 → breakpoint ~49bp into intron 19
_V4P_V7_V7_MAX = 20    # ≤ this → V7
_V4P_V7_V4P_MIN = 40   # ≥ this → V4'


def _classify_v4prime_v7(
    alk_chrom: str, alk_pos: int, exon_df: pd.DataFrame
) -> tuple[str, str]:
    """Sub-classify an E14:A20 fusion as V4', V7, or V4'/V7 (ambiguous).

    Uses the distance from the ALK breakpoint to the nearest boundary of ALK exon 20.
    V7 (del12) places the breakpoint ~12bp into intron 19; V4' (del49) ~49bp.

    Args:
        alk_chrom: Chromosome of the ALK breakpoint (no 'chr' prefix).
        alk_pos: Genomic position of the ALK breakpoint.
        exon_df: DataFrame from load_exon_reference().

    Returns:
        Tuple of (label, warning). Warning is empty string when classification is unambiguous.
    """
    exon20 = exon_df[
        (exon_df["gene"] == "ALK") &
        (exon_df["chrom"] == alk_chrom) &
        (exon_df["exon_number"] == 20)
    ]
    if exon20.empty:
        return "V4'/V7", "Could not locate ALK exon 20 in reference; E14:A20 fusion reported as V4'/V7."

    row = exon20.iloc[0]
    dist = int(min(abs(alk_pos - row["exon_start"]), abs(alk_pos - row["exon_end"])))

    if dist <= _V4P_V7_V7_MAX:
        return "V7", ""
    if dist >= _V4P_V7_V4P_MIN:
        return "V4'", ""
    return "V4'/V7", (
        f"ALK breakpoint is {dist}bp from ALK exon 20 boundary "
        f"(V7 ≤{_V4P_V7_V7_MAX}bp, V4' ≥{_V4P_V7_V4P_MIN}bp); "
        "reported as V4'/V7 (ambiguous)."
    )


def load_exon_reference(build: str = "hg38") -> pd.DataFrame:
    """Load the exon reference CSV for the given genome build.

    Args:
        build: Genome build string, either 'hg38' or 'hg19'.

    Returns:
        DataFrame with columns: gene, exon_number, chrom, exon_start, exon_end, strand.

    Raises:
        ValueError: If build is not 'hg38' or 'hg19'.
        FileNotFoundError: If the reference file does not exist.
    """
    if build not in REFERENCE_FILES:
        raise ValueError(f"Unsupported genome build '{build}'. Choose 'hg38' or 'hg19'.")
    path = REFERENCE_FILES[build]
    if not path.exists():
        raise FileNotFoundError(
            f"Exon reference not found at {path}. "
            "Run scripts/build_exon_reference.py to generate it."
        )
    df = pd.read_csv(path)
    df["exon_number"] = df["exon_number"].astype(int)
    df["exon_start"] = df["exon_start"].astype(int)
    df["exon_end"] = df["exon_end"].astype(int)
    df["chrom"] = df["chrom"].astype(str)
    return df


def parse_breakpoint(bp_string: str) -> tuple[str, int]:
    """Parse a breakpoint string into (chrom, position).

    Accepts formats:
      - '2_42522656'
      - 'chr2_42522656'

    Args:
        bp_string: Raw breakpoint string from user input.

    Returns:
        Tuple of (chrom_without_prefix, position_as_int).

    Raises:
        ValueError: If the string cannot be parsed.
    """
    bp_string = str(bp_string).strip()
    if "_" not in bp_string:
        raise ValueError(
            f"Cannot parse breakpoint '{bp_string}'. "
            "Expected format: '2_42522656' or 'chr2_42522656'."
        )
    chrom_part, pos_part = bp_string.rsplit("_", 1)
    chrom = chrom_part.lstrip("chr").lstrip("Chr")
    try:
        pos = int(pos_part)
    except ValueError:
        raise ValueError(
            f"Position '{pos_part}' in breakpoint '{bp_string}' is not an integer."
        )
    return chrom, pos


def lookup_exon(
    gene: str, chrom: str, pos: int, exon_df: pd.DataFrame
) -> tuple[int, str]:
    """Map a genomic position to an exon number for the given gene.

    Checks whether the position falls within any exon interval [exon_start, exon_end].
    If yes, returns (exon_number, "exon").
    If no, returns the nearest exon number by minimum boundary distance and "intron".

    Args:
        gene: 'EML4' or 'ALK'.
        chrom: Chromosome without 'chr' prefix (e.g. '2').
        pos: 1-based genomic position (integer).
        exon_df: DataFrame from load_exon_reference().

    Returns:
        Tuple of (exon_number, feature_type) where feature_type is 'exon' or 'intron'.
        Returns (-1, 'unknown') if no exon data found for the gene/chrom combination.
    """
    subset = exon_df[(exon_df["gene"] == gene) & (exon_df["chrom"] == chrom)]
    if subset.empty:
        return -1, "unknown"

    in_exon = subset[(subset["exon_start"] <= pos) & (pos <= subset["exon_end"])]
    if not in_exon.empty:
        return int(in_exon.iloc[0]["exon_number"]), "exon"

    # Intronic: return the exon whose end is the largest value ≤ pos.
    # Correct for both strands without needing to inspect the strand column:
    #   + strand (EML4): selects the last exon before the breakpoint (last retained).
    #   - strand (ALK):  selects the transcript-downstream exon (first retained),
    #                    because for minus-strand genes, downstream-in-transcript
    #                    equals lower genomic coordinates, so its exon_end < pos.
    upstream = subset[subset["exon_end"] <= pos]
    if not upstream.empty:
        nearest_idx = upstream["exon_end"].idxmax()
        return int(subset.loc[nearest_idx, "exon_number"]), "intron"
    # pos is upstream of all exons; return the first exon in the reference
    nearest_idx = subset["exon_start"].idxmin()
    return int(subset.loc[nearest_idx, "exon_number"]), "intron"


def _min_dist_to_gene(chrom: str, pos: int, gene: str, exon_df: pd.DataFrame) -> float:
    """Return the minimum distance from pos to any exon boundary of gene.

    Used to auto-detect which breakpoint belongs to which gene when the fusion
    name column order does not match the breakpoint column order.

    Args:
        chrom: Chromosome without 'chr' prefix.
        pos: Genomic position (integer).
        gene: 'EML4' or 'ALK'.
        exon_df: DataFrame from load_exon_reference().

    Returns:
        Minimum distance in base pairs, or float('inf') if no matching rows.
    """
    subset = exon_df[(exon_df["gene"] == gene) & (exon_df["chrom"] == chrom)]
    if subset.empty:
        return float("inf")
    distances = subset.apply(
        lambda r: min(abs(pos - r["exon_start"]), abs(pos - r["exon_end"])), axis=1
    )
    return float(distances.min())


def classify_fusion(
    row: pd.Series, exon_df: pd.DataFrame
) -> tuple[str, list[str], int, int]:
    """Classify a single fusion row into an EML4-ALK variant label.

    Implements the rules from Section 5 of fusions_app_spec.md:
      1. Skip non-EML4-ALK fusions.
      2. Auto-assign each breakpoint to EML4 or ALK by proximity to the reference.
      3. Look up exon number and feature type (exon/intron) for each breakpoint.
      4. Apply variant rules (V1, V2, V3a/b, NonCanonicalALK, OtherVariant).
      5. Append '_intron_junction' if either breakpoint is intronic.

    Args:
        row: Pandas Series with keys 'fusion_name', 'bp_a', 'bp_b'.
        exon_df: DataFrame from load_exon_reference().

    Returns:
        Tuple of (variant_label, warnings_list, eml4_exon_number, alk_exon_number).
        Exon numbers are -1 when the breakpoint could not be resolved.
    """
    warnings: list[str] = []
    fusion_name = str(row.get("fusion_name", "")).strip()

    parts = fusion_name.replace(" ", "").split("-")
    if len(parts) < 2:
        return "Not_EML4-ALK", warnings, -1, -1

    gene_a, gene_b = parts[0].upper(), parts[1].upper()
    if "EML4" not in (gene_a, gene_b) or "ALK" not in (gene_a, gene_b):
        return "Not_EML4-ALK", warnings, -1, -1

    bp_a_raw = row.get("bp_a")
    bp_b_raw = row.get("bp_b")

    if pd.isna(bp_a_raw) or pd.isna(bp_b_raw):
        warnings.append(f"Row has NA breakpoint value (fusion: {fusion_name}).")
        return "Unclassified_EML4-ALK", warnings, -1, -1

    try:
        chrom_a, pos_a = parse_breakpoint(str(bp_a_raw))
        chrom_b, pos_b = parse_breakpoint(str(bp_b_raw))
    except ValueError as exc:
        warnings.append(str(exc))
        return "Unclassified_EML4-ALK", warnings, -1, -1

    # Auto-detect gene assignment by proximity: try both orientations and pick the
    # one where each breakpoint is closer to its assigned gene's exons.
    dist_natural = (
        _min_dist_to_gene(chrom_a, pos_a, "EML4", exon_df)
        + _min_dist_to_gene(chrom_b, pos_b, "ALK", exon_df)
    )
    dist_swapped = (
        _min_dist_to_gene(chrom_a, pos_a, "ALK", exon_df)
        + _min_dist_to_gene(chrom_b, pos_b, "EML4", exon_df)
    )
    if dist_swapped < dist_natural:
        eml4_chrom, eml4_pos, eml4_bp_raw = chrom_b, pos_b, bp_b_raw
        alk_chrom, alk_pos, alk_bp_raw = chrom_a, pos_a, bp_a_raw
    else:
        eml4_chrom, eml4_pos, eml4_bp_raw = chrom_a, pos_a, bp_a_raw
        alk_chrom, alk_pos, alk_bp_raw = chrom_b, pos_b, bp_b_raw

    eml4_exon, eml4_feature = lookup_exon("EML4", eml4_chrom, eml4_pos, exon_df)
    if eml4_exon == -1:
        warnings.append(
            f"EML4 breakpoint {eml4_bp_raw} not found in reference (fusion: {fusion_name})."
        )
        return "Unclassified_EML4-ALK", warnings, -1, -1

    alk_exon, alk_feature = lookup_exon("ALK", alk_chrom, alk_pos, exon_df)
    if alk_exon == -1:
        warnings.append(
            f"ALK breakpoint {alk_bp_raw} not found in reference (fusion: {fusion_name})."
        )
        return "Unclassified_EML4-ALK", warnings, eml4_exon, -1

    rule_key = (eml4_exon, alk_exon)
    if rule_key in _VARIANT_RULES:
        label = _VARIANT_RULES[rule_key]
        if label == "V4'/V7":
            label, sub_warn = _classify_v4prime_v7(alk_chrom, alk_pos, exon_df)
            if sub_warn:
                warnings.append(sub_warn)
    elif alk_exon != 20:
        label = "EML4-ALK_NonCanonicalALK"
    else:
        label = "EML4-ALK_OtherVariant"
        warnings.append(
            f"Unexpected EML4 exon {eml4_exon} found in fusion '{fusion_name}' "
            "(known exons: 2, 6, 13, 14, 15, 17, 18, 20). "
            "Labeled as EML4-ALK_OtherVariant."
        )

    if eml4_feature == "intron" or alk_feature == "intron":
        label = f"{label}_intron_junction"

    return label, warnings, eml4_exon, alk_exon


def annotate_df(
    df: pd.DataFrame,
    col_map: dict[str, str],
    exon_df: pd.DataFrame,
) -> tuple[pd.DataFrame, list[str]]:
    """Annotate a DataFrame of fusions with EML4-ALK variant labels.

    Args:
        df: Input DataFrame with at least columns for fusion_name, bp_a, bp_b.
        col_map: Dict mapping logical names to actual column names:
                 {"fusion_name": ..., "bp_a": ..., "bp_b": ...}
                 "sample_id" is optional.
        exon_df: DataFrame from load_exon_reference().

    Returns:
        Tuple of (annotated_df, all_warnings) where annotated_df has the original
        columns plus 'EML4-ALK_VariantType', and all_warnings is a deduplicated
        list of warning strings.
    """
    result = df.copy()
    all_warnings: list[str] = []
    labels: list[str] = []
    gene_a_exons: list = []
    gene_b_exons: list = []

    for _, row in df.iterrows():
        work_row = pd.Series({
            "fusion_name": row[col_map["fusion_name"]],
            "bp_a": row[col_map["bp_a"]],
            "bp_b": row[col_map["bp_b"]],
        })
        label, row_warnings, eml4_exon, alk_exon = classify_fusion(work_row, exon_df)
        labels.append(label)
        gene_a_exons.append(eml4_exon if eml4_exon != -1 else pd.NA)
        gene_b_exons.append(alk_exon if alk_exon != -1 else pd.NA)
        all_warnings.extend(row_warnings)

    result["EML4-ALK_VariantType"] = labels
    result["Gene_A_Exon"] = gene_a_exons
    result["Gene_B_Exon"] = gene_b_exons

    seen: set[str] = set()
    unique_warnings: list[str] = []
    for w in all_warnings:
        if w not in seen:
            seen.add(w)
            unique_warnings.append(w)

    return result, unique_warnings

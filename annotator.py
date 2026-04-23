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
_VARIANT_RULES: dict[tuple[int, int], str] = {
    (13, 20): "V1",
    (20, 20): "V2",
    (6, 20): "V3a/b",
}


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

    # Intronic: return nearest exon by minimum distance to any boundary
    distances = subset.apply(
        lambda row: min(abs(pos - row["exon_start"]), abs(pos - row["exon_end"])),
        axis=1,
    )
    nearest_idx = distances.idxmin()
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


def classify_fusion(row: pd.Series, exon_df: pd.DataFrame) -> tuple[str, list[str]]:
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
        Tuple of (variant_label, warnings_list).
    """
    warnings: list[str] = []
    fusion_name = str(row.get("fusion_name", "")).strip()

    parts = fusion_name.replace(" ", "").split("-")
    if len(parts) < 2:
        return "Not_EML4-ALK", warnings

    gene_a, gene_b = parts[0].upper(), parts[1].upper()
    if "EML4" not in (gene_a, gene_b) or "ALK" not in (gene_a, gene_b):
        return "Not_EML4-ALK", warnings

    bp_a_raw = row.get("bp_a")
    bp_b_raw = row.get("bp_b")

    if pd.isna(bp_a_raw) or pd.isna(bp_b_raw):
        warnings.append(f"Row has NA breakpoint value (fusion: {fusion_name}).")
        return "Unclassified_EML4-ALK", warnings

    try:
        chrom_a, pos_a = parse_breakpoint(str(bp_a_raw))
        chrom_b, pos_b = parse_breakpoint(str(bp_b_raw))
    except ValueError as exc:
        warnings.append(str(exc))
        return "Unclassified_EML4-ALK", warnings

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
        return "Unclassified_EML4-ALK", warnings

    alk_exon, alk_feature = lookup_exon("ALK", alk_chrom, alk_pos, exon_df)
    if alk_exon == -1:
        warnings.append(
            f"ALK breakpoint {alk_bp_raw} not found in reference (fusion: {fusion_name})."
        )
        return "Unclassified_EML4-ALK", warnings

    rule_key = (eml4_exon, alk_exon)
    if rule_key in _VARIANT_RULES:
        label = _VARIANT_RULES[rule_key]
    elif alk_exon != 20:
        label = "EML4-ALK_NonCanonicalALK"
    else:
        label = "EML4-ALK_OtherVariant"
        warnings.append(
            f"Unexpected EML4 exon {eml4_exon} found in fusion '{fusion_name}'. "
            "Labeled as EML4-ALK_OtherVariant."
        )

    if eml4_feature == "intron" or alk_feature == "intron":
        label = f"{label}_intron_junction"

    return label, warnings


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

    for _, row in df.iterrows():
        work_row = pd.Series({
            "fusion_name": row[col_map["fusion_name"]],
            "bp_a": row[col_map["bp_a"]],
            "bp_b": row[col_map["bp_b"]],
        })
        label, row_warnings = classify_fusion(work_row, exon_df)
        labels.append(label)
        all_warnings.extend(row_warnings)

    result["EML4-ALK_VariantType"] = labels

    seen: set[str] = set()
    unique_warnings: list[str] = []
    for w in all_warnings:
        if w not in seen:
            seen.add(w)
            unique_warnings.append(w)

    return result, unique_warnings

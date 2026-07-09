"""
LEGACY / RETIRED: builds the old single-transcript Ensembl BioMart CSV schema
(gene, exon_number, chrom, exon_start, exon_end, strand) that annotator.py no
longer reads. annotator.py now loads reference/eml4_alk_exons_{build}.tsv,
whole-genome UCSC RefSeq (refGene) dumps, filtered at load time to
NM_019063.5 (EML4) and NM_004304.5 (ALK) -- see load_exon_reference() and
_CANONICAL_TRANSCRIPTS in annotator.py.

To regenerate the current .tsv reference files, pull the UCSC refGene table
for the build in question (e.g. via the UCSC Table Browser or `hgsql`/
`genome-mysql.soe.ucsc.edu`) with columns
name, chrom, strand, exonStarts, exonEnds, name2 and save as
reference/eml4_alk_exons_{build}.tsv. No filtering to EML4/ALK is required at
export time -- annotator.py filters on load.

This script is kept for provenance/history only and is not part of the
current reference build path.

---

Original docstring (BioMart-based EML4/ALK exon fetch):

Queries canonical transcripts:
  - EML4: ENST00000318522
  - ALK:  ENST00000389048

Writes:
  reference/eml4_alk_exons_hg38.csv  (GRCh38, grch38.ensembl.org)
  reference/eml4_alk_exons_hg19.csv  (GRCh37, grch37.ensembl.org)

Schema: gene, exon_number, chrom, exon_start, exon_end, strand

Usage:
  python3 scripts/build_exon_reference.py           # builds both
  python3 scripts/build_exon_reference.py --hg19    # hg19 only
  python3 scripts/build_exon_reference.py --hg38    # hg38 only
"""

import argparse
import csv
import sys
from pathlib import Path

import requests

BUILDS: dict[str, str] = {
    "hg38": "https://www.ensembl.org/biomart/martservice",
    "hg19": "https://grch37.ensembl.org/biomart/martservice",
}

TRANSCRIPTS: dict[str, str] = {
    "ENST00000318522": "EML4",
    "ENST00000389048": "ALK",
}

REFERENCE_DIR = Path(__file__).parent.parent / "reference"


def build_biomart_query(transcript_ids: list[str]) -> str:
    """Return a BioMart XML query string for the given transcript IDs.

    Requests exon rank, coordinates, strand, and chromosome for each transcript.
    """
    ids_csv = ",".join(transcript_ids)
    return f"""<?xml version="1.0" encoding="UTF-8"?>
<!DOCTYPE Query>
<Query virtualSchemaName="default" formatter="TSV" header="1" uniqueRows="1"
       count="" datasetConfigVersion="0.6">
  <Dataset name="hsapiens_gene_ensembl" interface="default">
    <Filter name="ensembl_transcript_id" value="{ids_csv}"/>
    <Attribute name="ensembl_transcript_id"/>
    <Attribute name="external_gene_name"/>
    <Attribute name="rank"/>
    <Attribute name="chromosome_name"/>
    <Attribute name="exon_chrom_start"/>
    <Attribute name="exon_chrom_end"/>
    <Attribute name="strand"/>
  </Dataset>
</Query>"""


def fetch_exon_data(biomart_url: str, transcript_ids: list[str]) -> list[dict]:
    """Query Ensembl BioMart at the given URL and return raw exon records.

    Args:
        biomart_url: BioMart service endpoint (hg38 or hg19).
        transcript_ids: List of Ensembl transcript IDs to query.

    Returns:
        List of dicts with keys matching the BioMart TSV header.

    Raises:
        RuntimeError: If the request fails or returns an unexpected response.
    """
    query_xml = build_biomart_query(transcript_ids)
    response = requests.get(biomart_url, params={"query": query_xml}, timeout=60)
    response.raise_for_status()

    text = response.text.strip()
    if not text:
        raise RuntimeError("BioMart returned an empty response.")
    if text.startswith("ERROR"):
        raise RuntimeError(f"BioMart error: {text}")

    lines = text.splitlines()
    if len(lines) < 2:
        raise RuntimeError(f"BioMart returned fewer lines than expected:\n{text[:500]}")

    header = lines[0].split("\t")
    records = []
    for line in lines[1:]:
        if not line.strip():
            continue
        fields = line.split("\t")
        records.append(dict(zip(header, fields)))
    return records


def parse_records(records: list[dict], transcript_map: dict[str, str]) -> list[dict]:
    """Convert raw BioMart records into the output schema rows.

    Args:
        records: Raw dicts from fetch_exon_data.
        transcript_map: Mapping of Ensembl transcript ID -> gene name.

    Returns:
        List of dicts with keys: gene, exon_number, chrom, exon_start, exon_end, strand.
        Sorted by gene then exon_number.
    """
    rows = []
    for rec in records:
        transcript_id = rec.get("Transcript stable ID", "")
        gene = transcript_map.get(transcript_id)
        if gene is None:
            continue

        chrom = rec.get("Chromosome/scaffold name", "").replace("chr", "")
        exon_number = int(rec.get("Exon rank in transcript", 0))
        exon_start = int(rec.get("Exon region start (bp)", 0))
        exon_end = int(rec.get("Exon region end (bp)", 0))
        strand_raw = rec.get("Strand", "1")
        strand = "+" if str(strand_raw) == "1" else "-"

        rows.append({
            "gene": gene,
            "exon_number": exon_number,
            "chrom": chrom,
            "exon_start": exon_start,
            "exon_end": exon_end,
            "strand": strand,
        })

    rows.sort(key=lambda r: (r["gene"], r["exon_number"]))
    return rows


def write_csv(rows: list[dict], output_path: Path) -> None:
    """Write exon reference rows to a CSV file.

    Args:
        rows: List of dicts from parse_records.
        output_path: Destination path; parent directories are created if needed.
    """
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = ["gene", "exon_number", "chrom", "exon_start", "exon_end", "strand"]
    with open(output_path, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def build_one(build: str) -> None:
    """Fetch and write the exon reference for a single genome build.

    Args:
        build: 'hg38' or 'hg19'.
    """
    url = BUILDS[build]
    output_path = REFERENCE_DIR / f"eml4_alk_exons_{build}.csv"
    transcript_ids = list(TRANSCRIPTS.keys())

    print(f"[{build}] Querying {url} ...")
    records = fetch_exon_data(url, transcript_ids)
    print(f"[{build}]   Received {len(records)} raw records.")

    rows = parse_records(records, TRANSCRIPTS)
    print(f"[{build}]   Parsed {len(rows)} exon rows.")

    if not rows:
        print(f"[{build}] ERROR: No rows parsed. Check BioMart response field names.")
        sys.exit(1)

    write_csv(rows, output_path)
    eml4_n = sum(1 for r in rows if r["gene"] == "EML4")
    alk_n = sum(1 for r in rows if r["gene"] == "ALK")
    print(f"[{build}]   Written to {output_path}  (EML4: {eml4_n}, ALK: {alk_n})")


def main() -> None:
    """Parse arguments and build the requested genome build references."""
    parser = argparse.ArgumentParser(description="Build EML4/ALK exon reference CSVs.")
    group = parser.add_mutually_exclusive_group()
    group.add_argument("--hg38", action="store_true", help="Build hg38 reference only.")
    group.add_argument("--hg19", action="store_true", help="Build hg19 reference only.")
    args = parser.parse_args()

    if args.hg38:
        builds = ["hg38"]
    elif args.hg19:
        builds = ["hg19"]
    else:
        builds = ["hg38", "hg19"]

    for b in builds:
        build_one(b)

    print("Done.")


if __name__ == "__main__":
    main()

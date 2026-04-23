# Fusion Annotator

A web app for classifying EML4-ALK gene fusions into known variant types given chromosomal breakpoint coordinates. Supports both hg19 (GRCh37) and hg38 (GRCh38).

**Try it:** [fusion-annotator.streamlit.app](https://fusion-annotator.streamlit.app)

---

## What it does

Given a fusion name and two breakpoint coordinates, the app maps each breakpoint to an exon in the canonical EML4 or ALK transcript and applies standard variant classification rules. It flags intron-spanning junctions and unexpected exon combinations.

### Classification rules

| EML4 exon | ALK exon | Label |
|---|---|---|
| 13 | 20 | `V1` |
| 20 | 20 | `V2` |
| 6 | 20 | `V3a/b` |
| any | ≠ 20 | `EML4-ALK_NonCanonicalALK` |
| ≠ 13, 20, 6 | 20 | `EML4-ALK_OtherVariant` |

If either breakpoint falls in an intron, `_intron_junction` is appended (e.g. `V3a/b_intron_junction`). Non-EML4-ALK fusions are labeled `Not_EML4-ALK`.

### Reference transcripts

| Gene | Ensembl transcript | Exons |
|---|---|---|
| EML4 | ENST00000318522 | 23 |
| ALK | ENST00000389048 | 29 |

---

## How to use it

Select the genome build (hg19 or hg38) matching your breakpoint coordinates before annotating. The app supports three input modes:

### Single entry

Enter one fusion manually: sample ID (optional), fusion name, Breakpoint A, and Breakpoint B.

### Paste rows

Paste tab- or comma-separated data including a header row. Columns are read by position:

| Position | Field |
|---|---|
| 1st | Sample ID |
| 2nd | Fusion name (e.g. `EML4-ALK`) |
| 3rd | Breakpoint A |
| 4th | Breakpoint B |

Example:

```
Sample	Fusion	BreakpointA	BreakpointB
PT001	EML4-ALK	2_29446394	2_42522656
```

### Upload CSV

Upload a CSV file with the same four-column structure. A 5-row preview is shown before annotating.

Both `2_29446394` and `chr2_29446394` breakpoint formats are accepted.

### Output

Results include the annotated variant type for each row and a summary frequency table. Annotated results can be downloaded as a CSV.

---

## Notes

- Only EML4-ALK fusions are classified; all other fusion names receive `Not_EML4-ALK`.
- Breakpoints that do not fall within any annotated exon or intron of the canonical transcript are reported as unresolved.
- Reference exon coordinates are derived from Ensembl BioMart (GRCh37/GRCh38) for the canonical transcripts listed above.

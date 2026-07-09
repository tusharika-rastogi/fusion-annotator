# Fusion Annotator

A web app for classifying EML4-ALK gene fusions into known variant types from splice-junction breakpoint coordinates reported by RNA-level fusion callers (e.g. STAR-Fusion, Arriba). Supports both hg19 (GRCh37) and hg38 (GRCh38).

DNA-level coordinates from structural variant callers (e.g. DELLY, Manta) may land inside an intron rather than at an exon boundary. The app will still attempt classification by resolving the intronic position to the nearest splice-junction exon, but this use case has not been systematically validated — results will carry the `_intron_junction` suffix (see below).

**Try it:** [fusion-annotator.streamlit.app](https://fusion-annotator.streamlit.app)

---

## What it does

Given a fusion name and two breakpoint coordinates, the app resolves each breakpoint to the nearest splice-junction exon in the canonical EML4 or ALK transcript and applies standard variant classification rules. It flags intron-spanning junctions and unexpected exon combinations. (One exception: the V3a/b check additionally consults a second EML4 transcript for one specific micro-exon case — see below.)

Breakpoints are auto-assigned to EML4 or ALK by proximity to each gene's reference exons, so **it does not matter which breakpoint you enter as "Breakpoint A" vs "Breakpoint B"** — the app tries both orientations and picks whichever one places each coordinate closer to its gene's exons.

### Classification rules

Exon numbers below represent the splice junction, not the raw genomic breakpoint location. For EML4 (plus strand), the listed exon is the last EML4 exon retained in the fusion transcript. For ALK (minus strand), the listed exon is the first ALK exon contributed to the fusion — ALK's exon 1 is at the *high*-coordinate end of chr2, so "exon 20" is genomically upstream of exon 1, not downstream.

| EML4 exon | ALK exon | Label | Notes |
|---|---|---|---|
| 13 | 20 | `V1` | ~33%; V6 (E13;ins69A20) also maps here |
| 20 | 20 | `V2` | ~9% |
| 6 | 20 | `V3a/b` | ~29%; RNA or long-read to distinguish a/b |
| 15 | 20 | `V4` | ~2% |
| 14 | 20 | `V4'` | <1%; ALK bp ~49bp into intron 19 (see below) |
| 14 | 20 | `V7` | <1%; ALK bp ~12bp into intron 19 (see below) |
| 2 | 20 | `V5a/b` | ~2% |
| 18 | 20 | `V5'` | ~2% |
| 17 | 20 | `V8a/b` | <1% |
| any | ≠ 20 | `EML4-ALK_NonCanonicalALK` | |
| other | 20 | `EML4-ALK_OtherVariant` | Novel or uncharacterized EML4 exon |

If either breakpoint does not fall within an annotated exon of the reference transcript, `_intron_junction` is appended (e.g. `V3a/b_intron_junction`). This does **not** mean the RNA read is unspliced or low-confidence — it means the fusion transcript's actual splice junction doesn't match the exon boundaries of the single reference transcript this tool checks against. This is frequently real, reproducible biology: some EML4-ALK isoforms splice at a site inside what the reference transcript calls an intron (e.g. V3b retains a 33-bp EML4 micro-exon not present in the canonical transcript; see Choi et al. 2008, PMID 18593892). **Do not treat `_intron_junction` as a quality filter** — dropping these rows from downstream analysis will disproportionately remove specific variant sub-types rather than removing noise. Non-EML4-ALK fusions are labeled `Not_EML4-ALK`.

#### V3a/b and the micro-exon junction

EML4 exon 6 fused to ALK exon 20 is reported as `V3a/b` regardless of whether the fusion transcript retains a 33-bp micro-exon derived from EML4 intron 6 (i.e. V3a vs V3b are not distinguished as separate labels). However, the app recognizes both known junction positions as clean splice junctions rather than intronic breakpoints, so neither gets the `_intron_junction` suffix:

| EML4 breakpoint (hg19 / hg38) | Result |
|---|---|
| `42491871` / `42264731` (exon 6 end) | `V3a/b` |
| `42492091` / `42264951` (micro-exon end) | `V3a/b` |
| any other intronic position | `V3a/b_intron_junction` (ambiguous) |

The micro-exon junction (from transcript `NM_001410776.1`) is literature-confirmed: Choi et al. 2008 (PMID 18593892) describes the original 33-bp intron-6 insertion; Wang et al. 2022 (PMID 36423218) reports it as "E6ins33;A20"; Hunt et al. 2023 (PMID 37255276) confirms the exact breakpoint coordinate.

#### V4' vs V7 sub-classification

Both V4' (E14;ins11del49A20) and V7 (E14;del12A20) involve EML4 exon 14 fused to ALK exon 20, but differ in how far the ALK breakpoint sits into intron 19. The app uses the exact ALK breakpoint coordinate to sub-classify:

| Distance from ALK exon 20 boundary | Label |
|---|---|
| ≤ 20 bp | `V7` |
| ≥ 40 bp | `V4'` |
| 21–39 bp | `V4'/V7` (ambiguous) |

A warning is emitted when the distance falls in the ambiguous zone. Thresholds are based on published breakpoint offsets (V7 ~12 bp, V4' ~49 bp; [PMC4761370](https://pmc.ncbi.nlm.nih.gov/articles/PMC4761370/)).

### Reference transcripts

| Gene | RefSeq transcript | Exons | Strand | Locus (hg19) |
|---|---|---|---|---|
| EML4 | NM_019063.5 | 23 | + | chr2:42,396,490–42,559,688 |
| ALK | NM_004304.5 | 29 | − | chr2:29,415,640–30,144,432 |

For hg38 coordinates, see `reference/eml4_alk_exons_hg38.tsv`.

The micro-exon junction check for `V3a/b` (see above) additionally uses the micro-exon from transcript `NM_001410776.1` (an EML4 isoform retaining a 33bp exon between exons 6 and 7), but exon numbering elsewhere is unaffected by it.

---

## How to use it

Select the genome build (hg19 or hg38) matching your breakpoint coordinates before annotating. **This choice is not validated against your data** — a build mismatch does not raise an error or warning, it silently produces incorrect exon numbers and variant labels, since the coordinates will simply resolve against the wrong reference's exon boundaries. If your results look biologically implausible (e.g. `NonCanonicalALK` for what should be a common variant), check the build selection first.

Both `2_29446394` and `chr2_29446394` breakpoint formats are accepted in all three input modes below.

The app supports three input modes:

### Single entry

Enter one fusion manually: sample ID (optional), fusion name, Breakpoint A, and Breakpoint B.

### Paste rows

Paste tab- or comma-separated data including a header row. Columns are read by position, **not** by header name — put a header row in even if its exact wording doesn't matter, since the first row is always discarded as a header. If you omit it, your first real fusion row will be silently dropped instead of annotated.

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

A 5-row preview is shown before annotating.

Column order for the two breakpoint fields does not matter — see auto-assignment note above.

### Upload CSV

Upload a CSV file with the same four-column structure (header row required, same silent-drop caveat as above). A 5-row preview is shown before annotating.

### Output

Results are shown on-screen as a table, plus a variant-type frequency summary, and can be downloaded as a CSV. The annotated table (and the downloaded CSV) adds three columns to your original data:

| Column | Meaning |
|---|---|
| `EML4-ALK_VariantType` | The variant label (e.g. `V3a/b`, `V2`, `Not_EML4-ALK`), with `_intron_junction` appended where applicable |
| `Gene_A_Exon` | The **EML4** exon number, regardless of which input column (Breakpoint A or B) that coordinate came from |
| `Gene_B_Exon` | The **ALK** exon number, regardless of which input column that coordinate came from |

**Warnings are not included in the downloaded CSV.** Any warnings (ambiguous V4'/V7 or V3a/b calls, unclassifiable breakpoints, parsing issues) are shown only as an on-screen banner above the results table at annotation time — read them before downloading, since they won't be visible again afterward.

---

## Notes

- Only EML4-ALK fusions are classified; all other fusion names receive `Not_EML4-ALK`.
- Breakpoints on a chromosome not covered by the reference (i.e. not chromosome 2 for EML4 and ALK) are reported as `Unclassified_EML4-ALK`. Any other breakpoint always resolves to some exon (falling back to the nearest exon boundary if intronic) — there is no "unresolved" middle case within chromosome 2.
- Reference exon coordinates are derived from the UCSC RefSeq (`refGene`) table (GRCh37/GRCh38) for the canonical transcripts listed above.

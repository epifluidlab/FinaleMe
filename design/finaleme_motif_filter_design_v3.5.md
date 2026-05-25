# FinaleMe extension: enzyme-targeted cfDNA fragment filter

**Status:** Design draft v3.5
**Changes from v3.4:**
- Added fifth output file `DFFB_dominant` containing fragments where F-III is the argmax F-profile. Two motif set flavors via `--dffb-set`: `canonical` (4 figure-labeled motifs from Zhou et al. Fig 2D) or `fprofile_argmax` (50 motifs, ~65% F-III capture). The DFFB_dominant output is the apoptotic-flux / cell-death readout — distinct from chromatin-accessibility outputs.
- Clarified that `L3_D1_dominant` excludes F-III, F-IV, F-V, AND F-VI (not just F-IV as previously implied). The argmax-based definition was already correct; only the explanatory text was misleading.
- QC PDF page 3 expanded to 3×2 grid for five output motif distributions. Page 4 (F-profile decomposition) now has six bars per F-profile group (input + five outputs).

**Changes from v3.3:**
- QC PDF restructured to include before/after comparisons for the two most diagnostic plots. New layout: page 2 shows the unfiltered 256-bin motif distribution; page 3 shows a 2×2 grid of the four output motif distributions for direct visual comparison; page 4 shows F-profile compositional decomposition with five grouped bars per F-profile (input + four outputs). This validates that each filter produced the expected motif signature rather than just describing the input.

**Changes from v3.2:**
- Added per-sample QC PDF output (§6.4) generated at the end of each run. Six pages: output retention summary, full 256-bin motif distribution with F-profile zones annotated, NNLS-based F-profile compositional decomposition, per-strand motif distribution comparison, QC funnel of filter stages, configuration metadata. Uses matplotlib + seaborn, with scipy.optimize.nnls for the F-profile fit. Disable with `--no-qc-pdf`.

**Changes from v3.1:**
- Reframed as a **standalone script** that imports from FinaleToolkit's API (`finaletoolkit.utils.frag_generator`, `finaletoolkit.utils.low_quality_read_pairs`). Not a modification to FinaleToolkit itself. Distribution as a single Python file alongside FinaleMe, with `finaletoolkit` as a runtime dependency.
- No changes to integration of FinaleToolkit's internals; only to how this tool is packaged and invoked.

**Changes from v3:**
- Strict pair-quality requirements added (§3.4, §3.5). The tool requires both Read 1 and Read 2 to pass MAPQ ≥ threshold. Default mode (`--mate-check require_mq_tag`) requires the BAM to have the `MQ` tag on Read 1 (populated by `samtools fixmate -m`) and aborts with a clear instruction if absent. Alternative mode `fetch_mate` performs per-fragment mate lookup at I/O cost. Multi-mapping (`XA`/`SA` tag) and 5′ soft-clip checks also added.

**Changes from v2 (v3 baseline):**
1. **Strand convention rewritten.** One motif per fragment, dispatched by `read_on_plus`: + strand → forward motif from `[start, start+k]`; − strand → reverse-complemented motif from `[stop-k, stop]`. This is the Read-1-derived 5′ end and is distinct from FinaleToolkit's `both_strands=True` (population statistics) and `both_strands=False` (which drops − strand fragments).
2. **Depletion file rescoped.** No longer just "DFFB-depleted" — now keeps only fragments whose motif is dominated by F-profile I (DNase1L3) or F-profile II (DNase1) among all six F-profiles. Excludes contributions attributed to F-III (DFFB), F-IV (CpG-rich, methylation-correlated), F-V (G-end), F-VI (oxidative). Renamed `L3_D1_dominant`.
3. **Match-policy section removed.** With one motif per fragment, no per-end disagreement to reconcile.

**Purpose:** Add a preprocessing utility to FinaleMe / FinaleToolkit that extracts cfDNA fragments preferentially originating from specific nucleases (DNase1L3, DNase1) based on the Read-1-derived 5′ end 4-mer motif. Output is four filtered fragment files for downstream analysis as a proof-of-concept that enzyme-targeted filtering improves chromatin/methylation inference.

---

## 1. Scope and design principles

**In scope.** A standalone Python script (single file) that:

- Imports `frag_generator`, `low_quality_read_pairs`, and `MIN_QUALITY` from `finaletoolkit.utils.utils` as a runtime dependency.
- Reads cfDNA fragments from BAM/CRAM or tabix-indexed fragment file (`.frag.gz`).
- Extracts a single 5′ end motif per fragment, dispatched by Read 1's mapping strand (§3).
- Classifies each fragment based on its motif against configurable enzyme motif sets.
- Writes four output files in the input format (BAM in → BAM out; frag.gz in → frag.gz out).

**Out of scope (v1).** Modification to FinaleToolkit itself; probabilistic soft assignment; methylation-aware classification; location-prior refinement; empirical re-derivation of motif sets.

**Design principles.**

1. Standalone, single-file script. Drops into a FinaleMe analysis directory without modifying any installed package.
2. FinaleToolkit treated as a library: import its API, call it, don't reimplement.
3. One motif per fragment for classification; clean per-fragment class assignment.
4. Configurable motif sets via TSV file; defaults from Zhou et al. 2023 F-profile attributions.
5. Single-pass over input → four output files.
6. Per-fragment classification preserves paired-end integrity (both Read 1 and Read 2 written to output when the fragment passes).

---

## 2. Inputs

### 2.1 Required arguments

- `--input PATH`: Input file, BAM/CRAM (with index) or `.frag.gz` (with `.tbi`). Format auto-detected by extension.
- `--reference PATH`: Reference genome in `.2bit` format.
- `--output-prefix PATH`: Path prefix for the four output files.

### 2.2 Optional arguments matching FinaleToolkit conventions

- `--motif-config PATH`: TSV file specifying motif sets (default: built-in Zhou et al. F-profile sets).
- `--motif-length INT`: 4 (default) or other k. Must match motif config.
- `--quality-threshold INT`: Minimum MAPQ. Default 30 (`MIN_QUALITY`).
- `--mate-check {require_mq_tag, fetch_mate, skip}`: How to enforce mate MAPQ. Default `require_mq_tag`. See §3.4.
- `--no-check-multimap`: Disable XA/SA tag check (default: enabled).
- `--max-5prime-softclip INT`: Maximum allowed soft-clip at Read 1's 5′ end. Default 5.
- `--fraction-low INT`: Minimum fragment length. Default 50.
- `--fraction-high INT`: Maximum fragment length. Default `None`.
- `--l3-set {ccnn,fprofile_argmax}`: Default `ccnn`. See §4.
- `--d1-set {tgnn,fprofile_argmax}`: Default `tgnn`. See §4.
- `--dffb-set {canonical,fprofile_argmax}`: Default `canonical` (4 figure-labeled motifs). `fprofile_argmax` uses 50-motif F-III-dominant set for meaningful retention. See §4.4.
- `--no-qc-pdf`: Disable QC PDF generation. Default: PDF is generated. See §6.4.
- `--threads INT`: Default 1.
- `--log PATH`: Per-fragment classification log (gzipped TSV). Optional.

Note: no `--both-strands` flag; the strand convention is fixed (one motif per fragment, dispatched by `read_on_plus`).

---

## 3. The 5′ end and strand convention

This implementation uses the **Read-1-derived 5′ end**, dispatched per fragment by the strand of Read 1's mapping. This is implemented on top of FinaleToolkit's `frag_generator`, which already filters out Read 2 (line 313 of `utils.py`: `if ... or read.is_read2: pass`) and emits a `read_on_plus` boolean per fragment.

### 3.1 Per-fragment dispatch

Given the yield from `frag_generator(...)`:
```
(contig, frag_start, frag_stop, mapq, read_on_plus)
```

The 5′ end motif of the fragment is determined by `read_on_plus`:

- **If `read_on_plus` is True** (Read 1 mapped to + strand):
  - Read 1's 5′ end is at the leftmost coordinate `frag_start`.
  - Motif = `refseq.sequence(contig, frag_start, frag_start + k)`, read 5′ → 3′ on + strand.

- **If `read_on_plus` is False** (Read 1 mapped to − strand):
  - Read 1's 5′ end is at the rightmost coordinate `frag_stop` (inclusive of the last base).
  - Genomic 4-mer on the + strand at `[frag_stop - k, frag_stop)`.
  - Motif = reverse complement of that 4-mer, giving the 5′ → 3′ read on the − strand.

This is exactly *one* preserved enzymatic 5′ cut site per fragment.

### 3.2 Relation to FinaleToolkit's existing options

FinaleToolkit's `region_end_motifs(both_strands=True)` (default) extracts **two** motifs per fragment — one from each end — and counts both into population-level statistics. This is appropriate for sample-level motif distributions.

FinaleToolkit's `both_strands=False` extracts only the forward motif and only from fragments where the strand flag is truthy, dropping − strand fragments entirely.

Neither matches what we need for per-fragment classification. The new tool's convention is a third pattern: extract one motif per fragment, with the end choice dispatched by `read_on_plus`. This:
- preserves one motif per fragment (clean class assignment),
- uses every fragment (no − strand drop),
- corresponds to the Read-1-derived 5′ end consistently across BAM and frag.gz inputs.

### 3.3 Motif lookup edge cases

- 4-mer overlaps reference `N` → fragment skipped, increment counter.
- Fragment near chromosome boundary → `py2bit.sequence()` raises `RuntimeError`; caught and skipped.
- Fragment length < k → caller-side check; warn and clamp `fraction_low` to k.

### 3.4 Strict pair-quality requirements (BAM input only)

For per-fragment classification both mates must be reliably positioned. FinaleToolkit's `low_quality_read_pairs` enforces most of what's needed (Read 1 MAPQ ≥ threshold, properly paired, mate not unmapped, not duplicate, not secondary/supplementary, not "both reverse"), but the mate's MAPQ check is conditional on the `MQ` SAM tag being present:

```python
# from utils.py low_quality_read_pairs():
try:
    if read.has_tag("MQ") and read.get_tag("MQ") < min_mapq:
        return True
except Exception:
    pass
return False
```

If the BAM was not processed with `samtools fixmate -m`, the `MQ` tag is absent and the mate MAPQ check is silently skipped — leaving fragments where Read 2 may be low-quality or multi-mapped while Read 1 passes. For per-fragment classification this corrupts `frag_stop`-derived motif extraction on − strand fragments (where `frag_stop` comes from Read 2's mapping).

**The tool requires the mate MAPQ check to be enforced.** Three modes via `--mate-check`:

- **`require_mq_tag` (default).** Check Read 1's `MQ` tag presence on first alignment seen. If absent, abort with: *"Input BAM lacks MQ tag on Read 1 alignments. Run `samtools fixmate -m input.bam fixed.bam` to populate MQ tags, or pass `--mate-check fetch_mate` for per-fragment mate lookup (slower)."*  Once present, FinaleToolkit's existing filter handles it correctly.

- **`fetch_mate`.** Use `sam_file.mate(read1)` per fragment to look up the mate and check its MAPQ directly. Doubles BAM read I/O but works on un-preprocessed BAMs. Acceptable for sample-scale; expensive for deep pools.

- **`skip`.** Skip the mate check entirely. Provided only for diagnostic comparison against legacy results; *not recommended for production.*

### 3.5 Additional pair-quality filters

Beyond mate MAPQ, three further per-fragment filters are applied by default:

1. **Multi-mapping check.** Skip reads with `XA` tag (bwa-mem alternative alignments) or `SA` tag (chimeric alignments). MAPQ ≥ 30 implies uniqueness in practice but this is belt-and-suspenders. Disable with `--no-check-multimap`.

2. **Same-chromosome mates.** Require `read.reference_name == read.next_reference_name`. Translocation / chimeric pairs would yield nonsensical fragment coordinates. Already implicit in `is_proper_pair` for most aligners but enforced explicitly here.

3. **5′ end soft-clip limit.** `--max-5prime-softclip INT` (default 5). Excessive soft-clipping at Read 1's 5′ end means the aligned `reference_start` (or `reference_end` for − strand) may not reflect the true enzymatic cut. Most cfDNA reads have 0–1 bp 5′ soft-clipping; > 5 bp is suspicious. Implementation: parse the leading (or trailing for − strand) CIGAR operation, skip if `S` op length exceeds the threshold.

For frag.gz input these checks are not applicable — the upstream BAM-to-frag.gz conversion is assumed to have applied them. Document this assumption clearly in the user guide so that frag.gz inputs are generated from properly QC'd BAMs.

---

## 4. Motif sets

### 4.1 Source quotes from Zhou et al. 2023 PNAS

- *F-profile I* (attributed to **DNASE1L3**): "F-profile I displayed a predominance of C-end motifs (55%) and was characterized by the 'CC' motifs."
- *F-profile II* (attributed to **DNASE1**): "F-profile II exhibited a major preference for T-end motifs (51%), with a significant enrichment observed for 'TG' motifs."
- *F-profile III* (attributed to **DFFB**): "F-profile III comprised a substantial proportion of A-end motifs (40%) and was characterized by the preference for C and T nucleotides at the third and fourth positions in the 4-mer motifs, respectively, in the 5′ to 3′ direction."

F-profiles IV, V, VI are unattributed in the paper; F-VI was associated with oxidative stress and discriminated HCC at AUC 0.97.

### 4.2 L3 and D1 motif sets — two options each

**Option A: paper-faithful families (CCNN / TGNN), default.** Matches Zhou et al.'s text directly.

**Option B: F-profile argmax sets.** Empirically derived from the F-profile frequency matrix shipped with FinaleToolkit: a motif belongs to L3 if F-I is the largest of the six F-profile contributions for that motif; belongs to D1 if F-II is largest. This captures motifs that are genuinely L3/D1-dominant even when they fall outside CCNN/TGNN, and excludes CCNN motifs that are actually F-IV-dominant (the five CG-containing CCNN motifs CCCG, CCGA, CCGC, CCGG, CCGT have F-IV > F-I).

```
# Default L3 set (--l3-set ccnn): CCNN family, 16 motifs
CCAA CCAC CCAG CCAT CCCA CCCC CCCG CCCT CCGA CCGC CCGG CCGT CCTA CCTC CCTG CCTT

# Alternative L3 set (--l3-set fprofile_argmax): 25 motifs where F-I is dominant
ACAA CAAA CAAC CAAG CAAT CAGA CAGT CATA CATG CCAA CCAC CCAG CCAT CCCA CCCC CCCT
CCTA CCTC CCTG CCTT CTAA CTAT CTGA CTTA TCTT

# Default D1 set (--d1-set tgnn): TGNN family, 16 motifs
TGAA TGAC TGAG TGAT TGCA TGCC TGCG TGCT TGGA TGGC TGGG TGGT TGTA TGTC TGTG TGTT

# Alternative D1 set (--d1-set fprofile_argmax): 53 motifs where F-II is dominant
CGAA CGAT CGCA CGCT CGGT CGTA CGTC CGTG CGTT TAAA TAAC TAAG TAAT TACA TACC TACG
TACT TAGA TAGC TAGG TAGT TATA TATC TATG TATT TCAA TCAC TCAG TCAT TCCA TCCC TCCG
TCCT TCTA TCTC TCTG TGAA TGAC TGAG TGAT TGCA TGCC TGCG TGCT TGGA TGGC TGGG TGGT
TGTA TGTC TGTG TGTT TTCA
```

### 4.3 L3_D1_dominant set (fourth output)

The fourth output retains fragments whose motif is **dominated by F-I or F-II** among all six F-profiles — explicitly excluding F-III (DFFB), F-IV (CpG-rich), F-V (G-end), and F-VI (oxidative). A motif belongs to `L3_D1_dominant` if `argmax(F-I, F-II, F-III, F-IV, F-V, F-VI) ∈ {0, 1}`.

This is the union of the two `fprofile_argmax` sets above — **78 motifs total** (25 F-I-dominant + 53 F-II-dominant).

Per-F-profile coverage of this 78-motif retention set:

| F-profile | Attributed enzyme | Captured by L3_D1_dominant | Excluded |
|---|---|---|---|
| F-I | DNase1L3 | 61.1% | 38.9% |
| F-II | DNase1 | 66.0% | 34.0% |
| F-III | DFFB | 15.2% | **84.8%** |
| F-IV | (CpG-rich) | 28.6% | **71.4%** |
| F-V | (G-end) | 20.4% | **79.6%** |
| F-VI | (oxidative) | 30.6% | **69.4%** |

The trade-off: retains 60-66% of L3 and D1 signal while excluding 70-85% of contributions from each of F-III, F-IV, F-V, and F-VI individually. This is the cleanest available criterion from the F-profile matrix for "exclude non-chromatin-informative fragmentation across the board."

### 4.4 DFFB_dominant set (fifth output)

The fifth output is the F-III analog: fragments whose motif is dominated by F-III (DFFB) among all six F-profiles. A motif belongs to `DFFB_dominant` if `argmax(F-I, F-II, F-III, F-IV, F-V, F-VI) = 2`.

Two flavors of the DFFB set are provided, mirroring the L3 / D1 pattern:

**Option A (`--dffb-set canonical`, default for paper-faithfulness):** The four figure-labeled motifs from Zhou et al. 2023 Figure 2 panel D:
```
ACCC ACCT ACTC ACTT
```

**Option B (`--dffb-set fprofile_argmax`, recommended for meaningful retention):** 50 motifs where F-III is the argmax F-profile:
```
AAAT AACA AACC AACG AACT AATA AATC AATG AATT ACAT ACCA ACCC ACCG ACCT ACTA
ACTC ACTG ACTT AGAT AGCA AGCC AGCG AGCT AGTA AGTC AGTT ATAT ATCC ATCT ATTT
CACA CACC CACG CACT CATC CATT CTCT GACA GACC GACG GACT GATC GATT GCCC GCCT
GCTC GCTT GGCT GGTT GTCT
```

Per-F-profile coverage of the argmax DFFB set (50 motifs):

| F-profile | Captured by DFFB_dominant |
|---|---|
| F-I (L3) | 11.1% |
| F-II (D1) | 9.3% |
| F-III (DFFB) | **65.3%** |
| F-IV | 8.8% |
| F-V | 21.0% |
| F-VI | 26.8% |

The argmax set captures ~65% of DFFB signal with modest leakage from F-V (G-end) and F-VI (oxidative). Expected retention in healthy plasma is ~18% — small enough that the output is enzyme-specific, large enough to be statistically useful.

The DFFB_dominant output is biologically the **apoptotic-flux / cell-death readout** referenced earlier in the conversation. Aggregated per region in clinical samples, it informs which cells are dying (via tissue-of-origin methylation of the retained fragments) — distinct from the chromatin-accessibility readouts of L3 and D1.

### 4.5 Expected per-end retention probabilities

From the F-profile matrix (P(motif | F-profile)):

| Motif set | # motifs | P(end \| F-I) | P(end \| F-II) | P(end \| F-III) | P(end \| F-IV) | P(end \| F-V) | P(end \| F-VI) |
|---|---|---|---|---|---|---|---|
| CCNN | 16 | 0.259 | 0.016 | 0.084 | 0.084 | 0.012 | 0.030 |
| TGNN | 16 | 0.082 | 0.221 | 0.018 | 0.045 | 0.014 | 0.034 |
| F-I-argmax | 25 | 0.343 | 0.057 | 0.045 | 0.014 | 0.014 | 0.048 |
| F-II-argmax | 53 | 0.220 | 0.451 | 0.054 | 0.197 | 0.190 | 0.258 |
| L3_D1_dominant | 78 | 0.563 | 0.508 | 0.099 | 0.211 | 0.204 | 0.306 |
| DFFB-4 (canonical) | 4 | 0.011 | 0.003 | 0.096 | 0.005 | 0.012 | 0.035 |
| DFFB-argmax | 50 | 0.111 | 0.093 | 0.653 | 0.088 | 0.210 | 0.268 |

In healthy plasma (~43% F-I, ~30% F-II, ~12% F-III, ~5% each in F-IV/V/VI), expected per-sample retention rates approximately:

| Output | Default config | F-argmax config |
|---|---|---|
| DNase1L3_enriched | ~13-17% (CCNN) | ~17-22% |
| DNase1_enriched | ~9-12% (TGNN) | ~22-28% |
| DNase1L3_or_DNase1_enriched | ~20-28% | ~38-45% |
| L3_D1_dominant | ~40-50% | same |
| DFFB_dominant | ~2-3% (canonical) | ~18-22% (argmax) |

---

## 5. Algorithm

### 5.1 High-level flow

```
1. Load motif config → build {motif → set of enzyme classes} lookup.
2. Open .2bit reference and four output writers.
3. For each fragment yielded by frag_generator(input, contig, mapq, ...):
     a. Extract one motif via strand dispatch (§3.1).
     b. If motif is None (N or boundary): skip, increment counter.
     c. Determine class memberships from the lookup.
     d. Write to appropriate output files.
4. Index output files; emit summary statistics.
```

### 5.2 Motif extraction

```python
def get_motif(contig, frag_start, frag_stop, read_on_plus, refseq, k=4):
    """Returns the Read-1-derived 5' end motif, or None if invalid."""
    try:
        if read_on_plus:
            kmer = refseq.sequence(contig, int(frag_start), int(frag_start + k))
        else:
            genomic_kmer = refseq.sequence(contig, int(frag_stop - k), int(frag_stop))
            kmer = reverse_complement(genomic_kmer)
    except RuntimeError:
        return None
    if len(kmer) != k or 'N' in kmer.upper():
        return None
    return kmer.upper()
```

### 5.3 Per-fragment classification

```python
def classify(motif, motif_to_classes):
    classes = motif_to_classes.get(motif, set())
    outputs = []
    if 'DNase1L3' in classes:
        outputs.append('DNase1L3_enriched')
    if 'DNase1' in classes:
        outputs.append('DNase1_enriched')
    if 'DNase1L3' in classes or 'DNase1' in classes:
        outputs.append('DNase1L3_or_DNase1_enriched')
    if 'L3_D1_dominant' in classes:
        outputs.append('L3_D1_dominant')
    if 'DFFB_dominant' in classes:
        outputs.append('DFFB_dominant')
    return outputs
```

The motif config TSV assigns motifs to multiple classes as needed. The five classes (`DNase1L3`, `DNase1`, `L3_D1_dominant`, `DFFB_dominant`, and any custom additions) are independent — a motif can belong to multiple classes and contribute to multiple output files. The `DNase1L3_or_DNase1_enriched` and `DFFB_dominant` outputs are mutually exclusive by construction (a motif's argmax F-profile is unique), but `DNase1L3_enriched` and `L3_D1_dominant` overlap substantially (most CCNN motifs are F-I-argmax).

### 5.4 Output writing

- BAM input → four BAM outputs. Use input header. Write both Read 1 and Read 2 when the fragment passes, preserving PE integrity. Index after writing with `pysam.index`.
- frag.gz input → four bgzipped TSV outputs matching the BED3+2 format (`chrom`, `start`, `stop`, `mapq`, `strand`). Tabix-index after writing.

Optionally: write the assigned motif as a BAM tag (`X5:Z:CCAA`) on output reads for downstream tools to use without re-fetching the reference. Off by default; enable with `--emit-motif-tag`.

---

## 6. Outputs

### 6.1 Five primary output files

For `--output-prefix /path/to/sample`:

- `/path/to/sample.DNase1L3_enriched.{bam,frag.gz}`
- `/path/to/sample.DNase1_enriched.{bam,frag.gz}`
- `/path/to/sample.DNase1L3_or_DNase1_enriched.{bam,frag.gz}`
- `/path/to/sample.L3_D1_dominant.{bam,frag.gz}`
- `/path/to/sample.DFFB_dominant.{bam,frag.gz}`

Each fully indexed (`.bai` or `.tbi`).

### 6.2 Summary statistics JSON

`/path/to/sample.motif_filter_stats.json`:

```json
{
  "input": "...",
  "reference": "...",
  "motif_length": 4,
  "l3_set": "ccnn",
  "d1_set": "tgnn",
  "total_fragments_processed": 0,
  "passed_qc": 0,
  "skipped_motif_invalid": 0,
  "skipped_mate_mapq_fail": 0,
  "skipped_multimap": 0,
  "skipped_excessive_softclip": 0,
  "mate_check_mode": "require_mq_tag",
  "motif_counts_plus_strand": {"AAAA": 0, "...": 0},
  "motif_counts_minus_strand": {"AAAA": 0, "...": 0},
  "output_counts": {
    "DNase1L3_enriched": 0,
    "DNase1_enriched": 0,
    "DNase1L3_or_DNase1_enriched": 0,
    "L3_D1_dominant": 0,
    "DFFB_dominant": 0
  },
  "retention_rates": {"...": 0.0},
  "config_hash": "sha256:...",
  "tool_version": "0.1.0",
  "command_line": "..."
}
```

### 6.3 Optional per-fragment log

`/path/to/sample.fragment_log.tsv.gz`: one row per processed fragment with `chrom`, `start`, `stop`, `read_on_plus`, `motif`, `assigned_classes`, `outputs`. Off by default; useful for debugging.

### 6.4 QC summary PDF

`/path/to/sample.motif_filter_qc.pdf` — a one-file visual QC report produced at the end of the run. Includes before/after comparisons for the two most informative diagnostic plots (motif distribution and F-profile decomposition) so the filtering can be validated visually.

**Page 1: Output retention summary.** Bar chart of fragments routed to each of the five outputs, alongside the expected retention range from the F-profile matrix (§4.5). Color bars green if within expected, orange if outside. Includes total input, fragments passing pair-quality, fragments with valid motif, and final per-output counts. Visual sanity check that the filter behaved as expected.

**Page 2: Unfiltered 256-bin motif distribution.** The full 4-mer end-motif histogram for the input sample, x-axis ordered as in Zhou et al. Figure 2 (A-end, C-end, G-end, T-end blocks). Background-shaded regions highlight the L3 (CCNN) zone in red, D1 (TGNN) zone in orange, and DFFB (figure-labeled 4-motif) positions in blue. The labeled top motifs from each F-profile are annotated. This is the baseline for visual comparison to filtered outputs.

**Page 3: Per-output motif distributions (3×2 grid).** Five motif histograms on one page, same x-axis layout as page 2:

- *Top-left: DNase1L3_enriched.* Should show C-end dominance approaching 100% (red zone highly enriched, other zones near-empty). Validates the L3 filter worked.
- *Top-right: DNase1_enriched.* Should show T-end dominance approaching 100% (orange zone highly enriched). Validates the D1 filter worked.
- *Middle-left: DNase1L3_or_DNase1_enriched.* Should show C-end + T-end dominance together (both zones enriched, A-end and G-end depleted).
- *Middle-right: L3_D1_dominant.* Should show the broader F-I-argmax + F-II-argmax signature — similar to the union but with the F-IV/V/VI-leaking CCNN/TGNN motifs removed. The difference from the union plot is subtle but diagnostic.
- *Bottom-left: DFFB_dominant.* Should show A-end dominance with position-3 C and position-4 T preference (matches Zhou et al. Figure 2 panel D). Validates the DFFB filter worked.
- *Bottom-right (sixth panel slot):* Used for a legend/annotation summary, or for the original input distribution (re-shown here for direct adjacency comparison).

Each subplot uses the same y-axis scale relative to its own total (frequency %, not count), so distribution *shapes* are directly comparable across panels.

**Page 4: F-profile compositional decomposition — input vs filtered outputs.** Six grouped bar plots (one per F-profile, I-VI), each with six bars side-by-side: unfiltered input, DNase1L3_enriched, DNase1_enriched, DNase1L3_or_DNase1_enriched, L3_D1_dominant, DFFB_dominant. Computed by non-negative least squares fit:

    motif_distribution_subset = w · F_profile_matrix    subject to w ≥ 0, ∑w = 1

where `F_profile_matrix` is the 6×256 matrix shipped with FinaleToolkit (`end_motif_f_profiles.tsv`).

Expected pattern:
- *Input:* ~40-45% F-I, ~25-35% F-II, ~10-15% F-III, balance in F-IV/V/VI (healthy plasma; cancer will shift this).
- *DNase1L3_enriched:* F-I should approach 80-95%, F-II low, F-III through F-VI low. Confirms L3 enrichment.
- *DNase1_enriched:* F-II should approach 80-95%, F-I low.
- *DNase1L3_or_DNase1_enriched:* F-I + F-II together should approach 90+%.
- *L3_D1_dominant:* Sharper exclusion of F-III through F-VI than the basic enrichment outputs. F-III should be 5-10%, F-IV/V/VI each ≤ 5%.
- *DFFB_dominant:* F-III should approach 60-80% with F-V and F-VI as secondary contributors (10-20% each). F-I and F-II should be low. Confirms DFFB enrichment.

This is the most informative single page for confirming all five filters worked correctly. If any expected enrichment doesn't materialize, the filter or motif config has a bug.

Bootstrap 95% CIs from resampling per-bin counts (1000 iterations). Note in footer: CIs reflect sampling noise only; F-profile model misspecification is not captured.

**Page 5: Per-strand motif distribution (unfiltered).** Two histograms side-by-side: forward-strand-fragment motifs and reverse-strand-fragment motifs (the latter after reverse-complementing). Should be approximately equal — each strand contributes ~50% of fragments with the same enzyme mixing. Large asymmetry indicates strand-specific library prep artifacts. A correlation coefficient or KL divergence between the two distributions is reported as a single number. This is on the unfiltered input only; per-output strand check is excessive.

**Page 6: QC funnel.** Bar chart of fragment counts at each filter stage:
- Total fragments scanned
- Passed pair-quality (MAPQ, MQ tag, flags)
- Passed multi-mapping filter (XA/SA tag check)
- Passed 5′ softclip filter
- Valid motif (not N, not at chromosome boundary)
- Total assigned to ≥ 1 output

Each bar annotated with drop count and percentage from the previous stage.

**Page 7: Configuration and metadata.** Text page with the command-line invocation, motif config used (with hashes), software versions, F-profile reference file path and hash, and a brief notes block. Aids reproducibility.

**Implementation notes.** Motif distributions per output file are accumulated during the single pass over the input (each fragment that goes to output X also increments a 256-bin counter for X). NNLS deconvolution adds ~5 seconds total across the five distributions (input + four outputs) with 1000 bootstraps each. Uses matplotlib + matplotlib.backends.backend_pdf.PdfPages, with seaborn for styling; scipy.optimize.nnls for the F-profile fit.

Disabled via `--no-qc-pdf`. Default: enabled.

---

## 7. Edge cases

| Case | Handling |
|---|---|
| Motif invalid (N or boundary) | Skip fragment; increment counter. |
| Fragment outside length range | Filtered by `frag_generator` upstream. |
| BAM lacks `MQ` tag and `--mate-check require_mq_tag` set | Hard error with instruction to run `samtools fixmate -m`. |
| Read 1 passes MAPQ but Read 2 below threshold | Skip fragment; increment `mate_mapq_fail` counter. |
| Read with `XA` or `SA` tag | Skip fragment if `--no-check-multimap` not set. |
| Read 1 with 5′ soft-clip > `--max-5prime-softclip` | Skip fragment; increment `excessive_softclip` counter. |
| BAM not sorted/indexed | Hard error from FinaleToolkit's reader; surface message. |
| Reference chromosome not in 2bit | Hard error, abort. |
| Motif config has invalid characters | Hard error during config load. |
| Empty enzyme class in config | Allowed (e.g., custom configs for sensitivity testing). |
| Single-end reads | Not supported; FinaleToolkit's `frag_generator` requires paired-end (filters on `is_read2`). |

---

## 8. Integration with FinaleToolkit / FinaleMe

This is a **standalone script** (single Python file, e.g., `motif_filter.py`) that imports from FinaleToolkit but is not part of it. Distribution as a script alongside FinaleMe's own utilities or in a sibling repository.

Imports used from FinaleToolkit:

```python
from finaletoolkit.utils.utils import frag_generator, low_quality_read_pairs, MIN_QUALITY
```

`frag_generator` handles BAM, CRAM, and frag.gz inputs uniformly and yields `(contig, start, stop, mapq, read_on_plus)` tuples — the script only adds motif extraction, classification, and output writing on top.

If `--mate-check fetch_mate` is used, the script needs its own iteration loop over `pysam.AlignmentFile` rather than `frag_generator`, because `frag_generator` doesn't expose the Read 1 object needed to call `.mate()`. Two code paths internally:

- Default path (`require_mq_tag`): uses `frag_generator` with FinaleToolkit's existing `low_quality_read_pairs` (which checks `MQ` tag when present). Cleanest.
- Fetch-mate path: custom iteration over `pysam.AlignmentFile`, calls `low_quality_read_pairs` on Read 1, then `sam_file.mate(read1)` and `low_quality_read_pairs` on the mate, then constructs the fragment tuple manually.

Expected workflow:

```bash
# Step 0 (one-time per BAM): ensure MQ tag is populated
samtools fixmate -m raw.bam fixed.bam
samtools sort -o sorted.bam fixed.bam
samtools index sorted.bam

# Step 1: motif-based filtering using the standalone script
python motif_filter.py \
    --input sorted.bam \
    --reference hg38.2bit \
    --output-prefix sample.filtered

# Step 2: run downstream analyses on each filtered output
for filt in DNase1L3_enriched DNase1_enriched DNase1L3_or_DNase1_enriched L3_D1_dominant; do
    finaleme run --input sample.filtered.${filt}.bam ...
done
```

If the proof-of-concept results justify it, the script can later be upstreamed into FinaleToolkit as a proper subcommand (`finaletoolkit motif-filter`). For now it lives as a script.

---

## 9. Dependencies

- Python 3.10+
- `finaletoolkit >= 0.11.1` (runtime dependency, imported as a library)
- `pysam`, `py2bit`, `numpy`, `scipy` (transitive via finaletoolkit; scipy used for NNLS in QC PDF)
- `matplotlib`, `seaborn` (for QC PDF generation; can be made optional via `--no-qc-pdf`)
- Standard library only beyond that.

Install via:
```
pip install finaletoolkit matplotlib seaborn
# then drop motif_filter.py into your project
```

The QC PDF dependencies (`matplotlib`, `seaborn`) are optional — if `--no-qc-pdf` is passed the script runs without importing them. Documented in the script's docstring.

No separate package distribution needed for v1; the script is self-contained. If later upstreamed into FinaleToolkit, packaging becomes trivial.

---

## 10. Testing

### 10.1 Unit tests

- `get_motif` returns correct motifs for + strand and − strand fragments. Verify reverse complement is computed correctly.
- Boundary cases: chromosome ends, N-containing motifs.
- Motif config loading: well-formed, malformed.
- Classification logic: motifs in single class, multiple classes, no class.

### 10.2 Reference reproducibility tests

- For a test BAM, compare per-motif fragment counts under this tool vs FinaleToolkit's `end_motifs` with `both_strands=False, negative_strand=False` (which uses forward end only on + strand fragments). These should differ in the expected way — our tool additionally captures − strand fragments via the reverse end.

### 10.3 Integration tests

- Small synthetic BAM (~20 fragments with hand-picked + and − strand reads): verify exactly which fragments end up in which outputs.

### 10.4 End-to-end smoke test

- Public small cfDNA dataset (BH01 chromosome subset from FinaleDB): retention rates should fall in ranges per §4.4. L3_D1_dominant should be larger than the union of CCNN+TGNN (since it's broader by definition).

---

## 11. Performance

- Same I/O profile as FinaleToolkit's `end_motifs`: streaming, modest memory.
- Per-fragment work is one motif lookup + classification + conditional file writes.
- Multi-threading via per-chromosome workers (mirror FinaleToolkit's pattern).

---

## 12. Open questions / decisions needed

1. **Default L3/D1 sets:** paper-faithful CCNN/TGNN vs F-profile argmax. CCNN/TGNN matches the paper text and is more interpretable; argmax is more empirically accurate (excludes CCNN motifs that are actually F-IV-dominant). Current default is CCNN/TGNN with argmax as a flag. Run sensitivity analysis during PoC.

2. **L3_D1_dominant always uses F-argmax union:** this is hardcoded since the depletion concept inherently requires the F-profile-derived definition. Acceptable?

3. **F-profile IV/V/VI as user-selectable enzymes:** F-VI was associated with oxidative stress and HCC discrimination at AUC 0.97 per Zhou et al. Could be valuable for cancer-detection use cases. Suggest as v1.1 feature.

4. **BAM motif tag (`X5:Z:CCAA`):** optional flag in v1, would carry the assigned motif on each output read for v2 soft-assignment forward compatibility. Recommend enabling by default.

---

## 13. Implementation timeline

Standalone script scope (no FinaleToolkit modification):

- Core script: motif extraction, classification, BAM/frag.gz output writers, using `finaletoolkit.utils.frag_generator` for the default path: ~2 days.
- `fetch_mate` mode (custom pysam iteration): ~0.5 day.
- 5′ soft-clip and XA/SA filters: ~0.5 day.
- QC PDF (6 pages, NNLS deconvolution, all plotting): ~1 day.
- Tests against synthetic and small public data: ~1 day.
- Documentation, example commands, sample config TSV: ~0.5 day.

Total: ~5-6 days for v1.

---

## 14. Success criteria for the proof-of-concept (downstream analyses)

After running the tool and feeding the four filtered outputs into FinaleMe:

1. **Fragment retention:** Approximately 13-17% (CCNN) or 17-22% (F-I-argmax) for DNase1L3_enriched in healthy plasma; lower in cancer if L3 is suppressed. L3_D1_dominant should retain ~40-50%.

2. **Nucleosome periodicity:** ~167 bp autocorrelation peak should be sharper in DNase1L3_enriched than in unfiltered. Quantify by peak height / width ratio.

3. **TSS depletion sharpness:** NDR depth at active TSSs should be deeper in DNase1L3_enriched. Measure by ratio of fragment density in [−100, +100] vs [−2000, −1000] ∪ [+1000, +2000] bp from TSS.

4. **FinaleMe primary accuracy:** Whatever the current FinaleMe accuracy metric is, compare filtered vs unfiltered runs. The L3_D1_dominant output should also be tested, since it captures more of the chromatin-informative signal at the cost of less enzyme-specificity than the pure CCNN set.

5. **Cancer vs healthy contrast:** Improvement larger in cancer samples than healthy → evidence the mechanism is enzyme-targeting, not generic denoising.

6. **Random-motif control:** Run with a randomly selected 16-motif set (matched retention rate) and verify FinaleMe improvement is substantially smaller than on DNase1L3_enriched. Rules out generic subsampling as the cause.

7. **Comparison across L3 set choices:** CCNN vs F-I-argmax — does the latter (which excludes the 5 CG-containing CCNN motifs that are F-IV-dominant) give cleaner FinaleMe accuracy? Informative for refining the default.

If 4–7 hold, the proof-of-concept supports the enzyme-targeting hypothesis.

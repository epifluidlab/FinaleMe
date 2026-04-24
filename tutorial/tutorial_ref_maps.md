# Building Custom Reference Methylation Maps for CGI+Shore Tissue-of-Origin Analysis

This tutorial explains how to build your own reference methylation atlas restricted to
CpG Island (CGI) and CGI shore regions, and how to use it for tissues-of-origin
deconvolution with FinaleMe-predicted methylation data using `BetaValueDeconvolution`.

> **You usually don't need this.** Pre-built CGI+shore atlases for both hg19 and hg38 are
> downloaded automatically by `./scripts/setup_references.sh` into `data/`:
>
> - `data/Atlas.CGI_shore.U250.l3.hg19.tsv`
> - `data/Atlas.CGI_shore.U250.l3.hg38.tsv`
> - `data/Atlas.pluse_microglia_astrocyte.CGI_shore.U250.l3.hg38.tsv`
>
> Use this tutorial only when you need custom cell types (e.g. tumor subtypes, disease-
> specific cell populations) not covered by the standard 38-cell-type panel, or when you
> want to tune marker selection for a specific cohort. The default Steps 1-6 workflow in
> [tutorial.md](tutorial.md) does not require `wgbstools` or `UXM_deconv`.

## Table of Contents

1. [Background](#1-background)
2. [Prerequisites](#2-prerequisites)
3. [Prepare Reference WGBS Data](#3-prepare-reference-wgbs-data)
4. [Run the CGI+Shore Marker Pipeline](#4-run-the-cgi-shore-marker-pipeline)
5. [Understanding the Pipeline Stages](#5-understanding-the-pipeline-stages)
6. [Use the Atlas for Deconvolution](#6-use-the-atlas-for-deconvolution)
7. [Customization and Advanced Usage](#7-customization-and-advanced-usage)
8. [Output File Reference](#8-output-file-reference)
9. [Training a FinaleMe calibration for `finaleme-too`](#9-training-a-finaleme-calibration-for-finaleme-too)
10. [Troubleshooting](#10-troubleshooting)

---

## 1. Background

### Why CGI+Shore restricted markers?

FinaleMe predicts per-CpG methylation from cell-free DNA (cfDNA) whole-genome sequencing
data. However, FinaleMe predictions are most accurate in **CpG Island (CGI) and CGI shore
regions** due to the distinct methylation patterns and higher CpG density in these regions.

The standard UXM deconvolution (Loyfer et al., Nature 2023) uses genome-wide methylation
markers to estimate tissue-of-origin composition. These markers were selected across the
entire genome, so only a small fraction (~5-10%) fall within CGI+shore regions. When
applied to FinaleMe predictions that are only reliable in CGI+shore regions, most markers
carry noisy or unreliable signal.

To address this, we provide a pipeline (`generate_cgi_shore_markers.py`) that:

1. Defines CGI+shore regions (CGI extended by +/-2kb)
2. Identifies tissue-specific differentially methylated regions within CGI+shore only
3. Builds a CGI+shore atlas using only these regions
4. Produces marker/atlas files directly usable with `BetaValueDeconvolution -markerRegions`

### What is a reference methylation atlas?

A reference atlas is a matrix where:
- **Rows** = marker regions (genomic intervals with tissue-specific methylation)
- **Columns** = cell types / tissues
- **Values** = fraction of unmethylated (U) or methylated (M) reads at each marker in each reference cell type

During deconvolution, the observed methylation pattern in your cfDNA sample is decomposed
as a linear mixture of the reference cell type profiles using non-negative least squares (NNLS).

---

## 2. Prerequisites

### 2.1 Software dependencies

Install the following tools before running the pipeline:

```bash
# 1. wgbs_tools (required for marker selection and atlas building)
git clone https://github.com/nloyfer/wgbs_tools.git
cd wgbs_tools
python setup.py
cd ..

# 2. UXM_deconv (required for atlas building in this custom-map workflow)
git clone https://github.com/nloyfer/UXM_deconv.git
cd UXM_deconv
pip install -r requirements.txt
cd ..

# 3. bedtools (required for region intersection)
# macOS:
brew install bedtools
# Linux:
conda install -c bioconda bedtools

# 4. Python packages
pip install pandas==1.5.3 numpy==1.21.6 matplotlib==3.4.3 scipy==1.9.0
```

### 2.2 Initialize the genome reference

If you haven't already, initialize the wgbs_tools genome reference:

```bash
wgbstools init_genome hg19
# or
wgbstools init_genome hg38
```

This creates the CpG index file required for pat/beta operations.

### 2.3 Directory structure

We recommend the following layout:

```
project/
├── reference_wgbs/           # Reference WGBS data
│   ├── betas/                # .beta files (one per sample)
│   ├── pats/                 # .pat.gz files (one per sample)
│   └── groups.csv            # Sample-to-cell-type mapping
├── results/                  # Pipeline output directory
└── scripts/                  # Contains generate_cgi_shore_markers.py
```

---

## 3. Prepare Reference WGBS Data

The pipeline requires three types of input from reference whole-genome bisulfite sequencing
(WGBS) experiments:

### 3.1 Reference beta files

Beta files store per-CpG methylation levels in a compact binary format. Each file represents
one WGBS sample.

**Format**: Binary file with `NR_SITES` rows x 2 columns of uint8 values:
- Column 0: methylated read count (0-255)
- Column 1: total read count (0-255)

**Where to obtain reference beta files**:

- **Option A: Loyfer et al. atlas** (recommended for standard cell types):
  Download from GEO accession [GSE186458](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE186458).
  Convert BAM files to beta using:
  ```bash
  wgbstools bam2pat sample.bam -o sample
  # This creates sample.pat.gz and sample.beta (or directly use their beta files in GEO)
  ```

- **Option B: Your own WGBS data**:
  If you have WGBS BAM files for your cell types of interest:
  ```bash
  # Convert each BAM to pat + beta format
  wgbstools bam2pat your_sample.bam -o your_sample --genome hg19
  ```

- **Option C: Public WGBS databases**:
  - ENCODE: https://www.encodeproject.org/ (filter for WGBS assay)
  - Roadmap Epigenomics: https://egg2.wustl.edu/roadmap/web_portal/
  - IHEC: https://epigenomesportal.ca/ihec/

### 3.2 Reference pat files

Pat files store per-read methylation patterns, which UXM uses to classify reads as
Unmethylated (U), miXed (X), or Methylated (M).

**Format**: Tab-separated, bgzip-compressed (`.pat.gz`):
```
chr    cpg_index    pattern    count
chr1   1            CTT        5
chr1   4            CCC        3
```

Pat files are generated alongside beta files by `wgbstools bam2pat`.

### 3.3 Groups file

The groups file maps each sample name to its cell type / tissue group. This is a CSV file
with at least two columns:

```csv
name,group
Adipocytes_1,Adipocytes
Adipocytes_2,Adipocytes
Blood_B_1,Blood-B
Blood_B_2,Blood-B
Blood_T_1,Blood-T
Liver_Hep_1,Liver-Hep
Liver_Hep_2,Liver-Hep
Neuron_1,Neuron
Neuron_2,Neuron
```

**Rules**:
- The `name` column must match the basename of the beta/pat files (without the `.beta`
  or `.pat.gz` extension). For example, if the beta file is `/data/betas/Liver_Hep_1.beta`,
  the name should be `Liver_Hep_1`.
- The `group` column defines the cell type label. All samples with the same `group` label
  are treated as biological replicates.
- You can optionally add an `include` column (boolean) to exclude certain samples.
- Having 2+ replicates per group enables statistical testing (t-test) during marker selection.

### 3.4 Example: preparing a small reference panel

Here is a minimal example with 3 cell types and 2 replicates each:

```bash
# Convert BAM files to wgbstools format
mkdir -p reference_wgbs/betas reference_wgbs/pats

for bam in Liver_Hep_1.bam Liver_Hep_2.bam \
           Neuron_1.bam Neuron_2.bam \
           Blood_T_1.bam Blood_T_2.bam; do
    sample=$(basename $bam .bam)
    wgbstools bam2pat $bam -o reference_wgbs/pats/$sample --genome hg19
    # This creates both .pat.gz and .beta
    mv reference_wgbs/pats/${sample}.beta reference_wgbs/betas/
done

# Create groups file
cat > reference_wgbs/groups.csv << 'EOF'
name,group
Liver_Hep_1,Liver-Hep
Liver_Hep_2,Liver-Hep
Neuron_1,Neuron
Neuron_2,Neuron
Blood_T_1,Blood-T
Blood_T_2,Blood-T
EOF
```

---

## 4. Run the CGI+Shore Marker Pipeline

### 4.1 Quick start (minimal command)

```bash
python scripts/generate_cgi_shore_markers.py \
  --genome hg19 \
  --betas reference_wgbs/betas/*.beta \
  --pats reference_wgbs/pats/*.pat.gz \
  --groups reference_wgbs/groups.csv \
  --num-markers 250 \
  --delta-means 0.4 \
  --unmeth-mean-thresh 0.1 \
  --meth-mean-thresh 0.5 \
  --min-cpg 1 \
  --max-cpg 1000 \
  --rlen 3 \
  --threads 10 \
  --out-dir results/cgi_shore_atlas/ \
  --wgbstools-path /path/to/wgbs_tools \
  --uxm-path /path/to/UXM_deconv \
  -v
```

This runs all five stages automatically and produces a ready-to-use atlas at
`results/cgi_shore_atlas/Atlas.CGI_shore.U250.l3.hg19.tsv`.

### 4.2 Full command with all options

```bash
python scripts/generate_cgi_shore_markers.py \
  --genome hg19 \
  --betas reference_wgbs/betas/*.beta \
  --pats reference_wgbs/pats/*.pat.gz \
  --groups reference_wgbs/groups.csv \
  --blocks /path/to/GSE186458_blocks.s205.bed.gz \
  --cgi-bed data/cpgIslandExt.hg19.bed \
  --shore-size 2000 \
  --chrom-sizes data/hg19.chrom.sizes \
  --num-markers 250 \
  --delta-means 0.4 \
  --unmeth-mean-thresh 0.1 \
  --meth-mean-thresh 0.5 \
  --min-cpg 1 \
  --max-cpg 1000 \
  --min-bp 50 \
  --max-bp 5000 \
  --rlen 3 \
  --threads 10 \
  --out-dir results/cgi_shore_atlas/ \
  --wgbstools-path /path/to/wgbs_tools \
  --uxm-path /path/to/UXM_deconv \
  --force \
  -v
```

### 4.3 Command-line options reference

| Option | Tutorial default | Description |
|--------|---------|-------------|
| `--genome` | (required) | Genome build: `hg19` or `hg38` |
| `--betas` | (required) | Reference WGBS beta files (glob pattern) |
| `--pats` | (required) | Reference WGBS pat.gz files (glob pattern) |
| `--groups` | (required) | Groups CSV file (sample-to-cell-type mapping) |
| `--out-dir` | (required) | Output directory for all results |
| `--wgbstools-path` | (required) | Path to wgbs_tools installation |
| `--uxm-path` | (required) | Path to UXM_deconv installation |
| `--cgi-bed` | auto-download | CpG Island BED file from UCSC |
| `--shore-size` | 2000 | Shore extension in bp around each CGI |
| `--chrom-sizes` | none | Chromosome sizes file for boundary capping |
| `--blocks` | optional | Pre-existing blocks file (recommended for reproducibility) |
| `--num-markers` | 250 | Target markers per cell type |
| `--delta-means` | 0.4 | Minimum methylation difference for markers |
| `--unmeth-mean-thresh` | 0.1 | Upper methylation bound for unmethylated markers |
| `--meth-mean-thresh` | 0.5 | Lower methylation bound for methylated markers |
| `--min-cpg` | 1 | Minimum CpGs per candidate block |
| `--max-cpg` | 1000 | Maximum CpGs per candidate block |
| `--min-bp` | 50 | Minimum block length in bp |
| `--max-bp` | 5000 | Maximum block length in bp |
| `--rlen` | 3 | Minimum CpGs per read for U/X/M classification |
| `--threads` | 10 | Number of parallel threads |
| `--force` | false | Overwrite existing output files |
| `--verbose` | false | Print detailed progress information |

---

## 5. Understanding the Pipeline Stages

### Stage 1: Generate CGI+Shore BED Regions

The pipeline first defines the genomic regions where FinaleMe predictions are reliable.

- Downloads UCSC CpG Island annotation (or uses a user-provided file)
- Extends each CGI by `--shore-size` bp in both directions (default: +/-2kb)
- Merges overlapping intervals to avoid double-counting

For hg19, this typically produces ~22,000 merged regions covering ~180 Mb (~6% of the genome).

**Output**: `cgi_plus_shore.{genome}.bed`

### Stage 2: Generate Candidate Blocks

The pipeline needs candidate genomic intervals ("blocks") to test for tissue-specific
methylation. There are two approaches:

- **Default (segmentation)**: Runs `wgbstools segment` on your reference beta files to
  identify regions with coherent methylation patterns, then intersects these with CGI+shore
  regions. This produces biologically meaningful boundaries.

- **Pre-existing blocks**: If you already have a blocks file (e.g., from a previous
  segmentation, or just directly use the hg19/hg38 block files from GSE186458), pass it via `--blocks` and the pipeline will just intersect it with
  CGI+shore.

Blocks are filtered to retain only those with 3-50 CpGs and 50-5000 bp length.

**Output**: `cgi_shore_blocks.{genome}.bed`

### Stage 3: Find Tissue-Specific Markers

This is the core marker selection stage. For each cell type, the pipeline identifies
blocks where that cell type has distinctly different methylation compared to all others.

**Marker types**:
- **U markers**: The target cell type is *unmethylated* while background is methylated
- **M markers**: The target cell type is *methylated* while background is unmethylated

**Selection criteria**:
- Mean methylation difference (delta_means) >= threshold
- t-test p-value <= 0.05
- Sufficient coverage across replicates

**Adaptive threshold relaxation**: If a cell type has fewer than 5 markers at the initial
threshold (default 0.35), the pipeline automatically relaxes to 0.25, then 0.20, then 0.15,
ensuring every cell type gets at least some markers.

**Output**: `markers/Markers.{CellType}.bed` and merged `Markers.CGI_shore.U{N}.{genome}.tsv`

### Stage 4: Build UXM Atlas

The selected markers are converted into a UXM-compatible atlas by computing U/X/M read
fractions for each marker in each reference cell type.

For each marker region and each reference pat file:
1. `wgbstools homog` counts reads that are fully Unmethylated (U), miXed (X), or fully
   Methylated (M)
2. Counts are aggregated across replicates within each cell type
3. Normalized to fractions that sum to 1.0

**Output**: `Atlas.CGI_shore.U{N}.l{rlen}.{genome}.tsv`

### Stage 5: Validation Report

The pipeline generates a summary report including:
- CGI+shore region statistics (count, total coverage)
- Number of candidate blocks
- Markers per cell type (with shortfall warnings)
- Marker quality statistics (mean delta, p-values)
- Overlap analysis with the original genome-wide UXM U250 markers

**Output**: `report.txt`

---

## 6. Use the Atlas for Deconvolution

Once you have the atlas, you can deconvolve any FinaleMe-predicted cfDNA sample.

### 6.1 Generate FinaleMe decode output

First, run FinaleMe decode to produce `*.prediction.bed.gz`:

```bash
JAR="target/FinaleMe-0.63-jar-with-dependencies.jar"

java -Xmx20G -cp "$JAR" \
  edu.northwestern.epifluidlab.finaleme.hmm.FinaleMe \
  results/sample.FinaleMe.model \
  results/sample.cpg_features.hg19.bed.gz \
  results/sample.decode.prediction.bed.gz \
  -decodeModeOnly
```

### 6.2 Run `BetaValueDeconvolution` with the custom atlas

The atlas TSV produced by `generate_cgi_shore_markers.py` has embedded per-cell-type
values, so you can pass it directly with `-refPanel` (no separate `-refBetas`/`-refGroups`
needed). With v0.63 defaults, this also runs stratified bootstrap + permutation:

```bash
JAR="target/FinaleMe-0.63-jar-with-dependencies.jar"

java -Xmx20G -cp "$JAR" \
  edu.northwestern.epifluidlab.finaleme.utils.BetaValueDeconvolution \
  -refPanel results/cgi_shore_atlas/Atlas.CGI_shore.U250.l3.hg19.tsv \
  -cpgIndex data/CpG_index.hg19.bed.gz \
  -output results/sample.deconv.beta.tsv \
  results/sample.decode.prediction.bed.gz
```

If you want the legacy two-file mode (separate marker BED + reference betas):

```bash
java -Xmx20G -cp "$JAR" \
  edu.northwestern.epifluidlab.finaleme.utils.BetaValueDeconvolution \
  -markerRegions results/cgi_shore_atlas/Atlas.CGI_shore.U250.l3.hg19.tsv \
  -refBetas results/cgi_shore_atlas/reference_wgbs/betas/beta_list.txt \
  -refGroups results/cgi_shore_atlas/groups_fixed.csv \
  -cpgIndex data/CpG_index.hg19.bed.gz \
  -output results/sample.deconv.beta.tsv \
  results/sample.decode.prediction.bed.gz
```

Make `beta_list.txt` for that mode:

```bash
ls reference_wgbs/betas/*.beta > results/cgi_shore_atlas/reference_wgbs/betas/beta_list.txt
```

For the full set of `BetaValueDeconvolution` options (bootstrap, permutation, FDR), see
[tutorial.md §8](tutorial.md).

### 6.3 Interpret results

With v0.63 defaults the output is long format, one row per `(sample, cell_type)`:

```
sample   cell_type   proportion   CI_lower   CI_upper   p_value   q_value   significant   p_source      n_replicates
sample1  Blood-T     0.156        0.142      0.171      5.3e-05   5.3e-05   YES           permutation   10000
sample1  Liver-Hep   0.089        0.071      0.104      1.2e-03   1.5e-03   YES           permutation   10000
sample1  Blood-B     0.023        0.000      0.041      0.213     0.213     NO            permutation   10000
...
```

`proportion` values represent the estimated fraction of cfDNA originating from each cell
type and sum approximately to 1.0 across cell types within a sample (minor deviations are
normal due to NNLS renormalization).

For cohort-level differential analysis (e.g. Disease vs Control with covariates), feed
this output to `scripts/too_diff_analysis.py` — see [tutorial.md §9](tutorial.md).


---

## 7. Customization and Advanced Usage

### 7.1 Choosing the number of markers

- **250 markers per cell type** (`--num-markers 250`): Recommended preset in this tutorial for `BetaValueDeconvolution`.
- **25 markers per cell type** (`--num-markers 25`): Use when references are sparse or when you want stricter/high-confidence marker subsets.

### 7.2 Adjusting the shore size

The default shore size is 2kb, a standard definition in the epigenomics literature. You
can adjust it:

```bash
# Narrow: CGI + 1kb shore
--shore-size 1000

# Wide: CGI + 4kb shore
--shore-size 4000
```

Wider shores include more candidate regions but may extend into areas where FinaleMe
predictions are less accurate.

### 7.3 Using custom cell types

To build an atlas for a specific subset of cell types, create a groups file with only
the cell types you need:

```csv
name,group
Blood_T_1,Blood-T
Blood_T_2,Blood-T
Liver_Hep_1,Liver-Hep
Liver_Hep_2,Liver-Hep
Tumor_1,Tumor
Tumor_2,Tumor
```

This is useful when you have custom WGBS data for cell types not in the Loyfer et al.
reference (e.g., tumor subtypes, disease-specific cell populations).

### 7.4 Relaxing marker selection thresholds

If you're working with cell types that have subtle methylation differences in CGI+shore
regions, you can start with a lower threshold:

```bash
--delta-means 0.25
```

The pipeline's adaptive relaxation will still try progressively lower thresholds for cell
types with insufficient markers.

### 7.5 Using pre-existing blocks

If you have previously segmented the genome or want to use a specific block definition:

```bash
python scripts/generate_cgi_shore_markers.py \
  --blocks /path/to/my_blocks.bed \
  --genome hg19 \
  ... # other options
```

The blocks file must have at least 5 columns: `chr`, `start`, `end`, `startCpG`, `endCpG`.

### 7.6 Building for both hg19 and hg38

Run the pipeline separately for each genome build:

```bash
# hg19
python scripts/generate_cgi_shore_markers.py \
  --genome hg19 \
  --betas ref_hg19/betas/*.beta \
  --pats ref_hg19/pats/*.pat.gz \
  --groups ref_hg19/groups.csv \
  --out-dir results/atlas_hg19/ \
  ...

# hg38
python scripts/generate_cgi_shore_markers.py \
  --genome hg38 \
  --betas ref_hg38/betas/*.beta \
  --pats ref_hg38/pats/*.pat.gz \
  --groups ref_hg38/groups.csv \
  --out-dir results/atlas_hg38/ \
  ...
```

---

## 8. Output File Reference

After running the pipeline, the output directory contains:

```
results/cgi_shore_atlas/
├── cgi_plus_shore.hg19.bed                    # Merged CGI+shore regions
├── cgi_shore_blocks.hg19.bed                  # Candidate blocks within CGI+shore
├── markers/                                   # Per-cell-type marker files
│   ├── Markers.Blood-T.bed
│   ├── Markers.Liver-Hep.bed
│   ├── Markers.Neuron.bed
│   ├── ...
│   └── pass_0/                                # Intermediate results per threshold pass
├── Markers.CGI_shore.U250.hg19.tsv            # All markers merged (input to atlas builder)
├── Atlas.CGI_shore.U250.l3.hg19.tsv           # Final atlas (input to BetaValueDeconvolution -markerRegions)
└── report.txt                                 # Validation summary
```

### Key file formats

**Markers file** (`Markers.CGI_shore.U250.hg19.tsv`):

| Column | Description |
|--------|-------------|
| #chr | Chromosome |
| start | Region start (0-based) |
| end | Region end |
| startCpG | First CpG index in region |
| endCpG | Last CpG index in region |
| target | Cell type this marker is specific to |
| region | Region string (e.g., chr1:12345-67890) |
| lenCpG | Number of CpGs in region |
| bp | Region length in bp |
| tg_mean | Mean methylation in target cell type |
| bg_mean | Mean methylation in background |
| delta_means | Difference between bg_mean and tg_mean |
| delta_quants | Quantile-based difference |
| delta_maxmin | Max-min based difference |
| ttest | t-test p-value |
| direction | U (unmethylated in target) or M (methylated in target) |

**Atlas file** (`Atlas.CGI_shore.U250.l3.hg19.tsv`):

| Column | Description |
|--------|-------------|
| chr | Chromosome |
| start | Region start |
| end | Region end |
| startCpG | First CpG index |
| endCpG | Last CpG index |
| target | Cell type this marker is specific to |
| name | Region string |
| direction | U or M |
| CellType1 | U/M fraction for cell type 1 |
| CellType2 | U/M fraction for cell type 2 |
| ... | (one column per cell type) |

---

## 9. Training a FinaleMe calibration for `finaleme-too`

The Python `finaleme-too` package (a separate tool from `BetaValueDeconvolution`) can use a
per-CpG-density calibration model to map FinaleMe predictions onto the same scale as
WGBS truth before running the beta-binomial deconvolution. Training this calibration
needs three inputs:

1. **WGBS counts** — Bis-SNP `*.cpg.6plus2.bed` output, one per sample
2. **Matched FinaleMe counts** — the `*.prediction.bed.gz` output from FinaleMe decode, one per sample
3. **Per-row CpG density** — a TSV with `chrom start end cpg_density`

This section shows how to produce #3 (the `region_annotation.tsv`).

### 9.1 What is `region_annotation.tsv`?

For every row that appears in the matched WGBS/FinaleMe tables, the calibration trainer
needs a `cpg_density` value so it can bin rows into CpG-density classes (one regression
slope/intercept per bin — see math doc §6.1). The density is simply the number of CpGs in
a local window around each row divided by the window size (default 1000 bp).

| Column | Description |
|--------|-------------|
| `chrom` | Chromosome (no `chr` prefix is needed; both conventions are accepted) |
| `start` | Start position (0-based) |
| `end` | End position (half-open) |
| `cpg_density` | (# CpGs in window centered on the row) / window_size |

### 9.2 Three ways to produce it

#### Option A — Let `train-calibration` auto-generate it (simplest)

Pass `--cpg-index` to `finaleme-too train-calibration` and skip `--region-annotation`
entirely. The trainer will compute density on the fly for the unique `(chrom, start, end)`
rows that actually appear in the merged training data, avoiding a genome-wide file.

```bash
finaleme-too train-calibration \
  --matched-wgbs wgbs_samples.tsv \
  --matched-finaleme finaleme_samples.tsv \
  --cpg-index data/CpG_index.hg19.bed.gz \
  --region-annotation-window 1000 \
  --n-bins-candidates 4,6,8 \
  --output calibration_params.json \
  --report calibration_report.json
```

The `--cpg-index` argument accepts the same CpG index BED that `BetaValueDeconvolution`
uses with `-cpgIndex` (downloaded by `scripts/setup_references.sh` into `data/`).

#### Option B — Build it once with `make-region-annotation`

If you want a reusable file (e.g. across many training runs, or to inspect), use the
dedicated subcommand:

```bash
finaleme-too make-region-annotation \
  --regions data/CpG_index.hg19.bed.gz \
  --cpg-index data/CpG_index.hg19.bed.gz \
  --window 1000 \
  --output region_annotation.hg19.tsv
```

Here we pass the CpG index as BOTH the regions and the reference CpG catalog, which
produces a genome-wide per-CpG density file. The output TSV has columns
`chrom  start  end  cpg_density` and is directly consumable by
`finaleme-too train-calibration --region-annotation region_annotation.hg19.tsv`.

If you only care about a specific set of markers (not the full CpG index), pass your
marker BED file to `--regions` instead:

```bash
finaleme-too make-region-annotation \
  --regions my_markers.bed \
  --cpg-index data/CpG_index.hg19.bed.gz \
  --output my_markers.region_annotation.tsv
```

#### Option C — Build it yourself

The file format is trivial. Any script that counts CpGs in a local window around each
row and writes `chrom start end cpg_density` as TSV will work. For example, using
`bedtools` + the CpG index:

```bash
# Compute a 1kb window around each CpG, intersect with the CpG index, count overlaps,
# and divide by 1000. Rough sketch:
awk 'BEGIN{OFS="\t"} {print $1, ($2-500 < 0 ? 0 : $2-500), $3+500}' CpG_index.hg19.bed.gz \
  | bedtools intersect -a - -b CpG_index.hg19.bed.gz -c \
  | awk 'BEGIN{OFS="\t"; print "chrom","start","end","cpg_density"}
         {print $1, $2+500, $2+501, $4/1000}' \
  > region_annotation.hg19.tsv
```

The Python helper is cleaner and handles chrom-prefix normalization automatically.

### 9.3 Which window size should I use?

Default is **1000 bp**. Larger windows smooth out local CpG variation (more stable
density estimates but less resolution); smaller windows pick up tight clusters (more
resolution but noisier). The `--region-annotation-window` flag on both
`make-region-annotation` and `train-calibration` lets you tune this per run.

### 9.4 How many bins should I choose?

`train-calibration` already tunes this for you via leave-one-sample-out cross validation
over `--n-bins-candidates` (default `4,6,8,10,12,16`) and picks the B that minimizes
CV-RMSE. The chosen B and the full candidate table are written to the `--report` JSON.

---

## 10. Troubleshooting

### "wgbstools not found"

Ensure wgbs_tools is installed and the path is correct:
```bash
# Check if wgbstools is accessible
/path/to/wgbs_tools/wgbstools --help

# Pass the correct path
--wgbstools-path /path/to/wgbs_tools
```

### "bedtools not found"

Install bedtools:
```bash
# macOS
brew install bedtools
# Linux
conda install -c bioconda bedtools
```

### Too few markers for some cell types

This is expected when restricting to CGI+shore. The pipeline automatically relaxes
thresholds for underrepresented cell types. Check `report.txt` for details.

If critical cell types have very few markers (<5), consider:
- Increasing `--shore-size` to include more candidate regions
- Lowering `--delta-means` to 0.2 or 0.15
- Adding more reference samples for that cell type

### Memory errors during segmentation

If `wgbstools segment` runs out of memory with many beta files, the pipeline
automatically uses only the first 10 files. You can also provide pre-computed blocks:
```bash
# Segment separately with resource controls
wgbstools segment --betas subset_of_betas/*.beta -o my_blocks.bed -@ 4

# Then use pre-existing blocks
--blocks my_blocks.bed
```

### Slow atlas building

Atlas building calls `wgbstools homog` for each marker x sample combination. Speed up with:
```bash
--threads 10   # Use more threads
```

Results are cached in `UXM_deconv/tmp_dir/`, so re-runs are much faster.

### Atlas has NaN values

NaN values appear when a reference sample has zero coverage at a marker region. This
is handled gracefully by `BetaValueDeconvolution` but may reduce accuracy. Ensure your reference
WGBS data has sufficient genome-wide coverage (>10x recommended).

---

## References

- Loyfer N, et al. "A DNA methylation atlas of normal human cell types." *Nature* 613, 355-364 (2023). https://doi.org/10.1038/s41586-022-05580-6
- wgbs_tools: https://github.com/nloyfer/wgbs_tools
- UXM_deconv: https://github.com/nloyfer/UXM_deconv

# Stage 1: Improved FinaleMe + Fragment-Level TOO Deconvolution for Low-Coverage Samples

## Design Document

**Version:** 2.0
**Date:** April 11, 2026

---

## 1. Overview

### 1.1 Problem Statement

FinaleMe predicts CpG methylation from cfDNA WGS fragmentation patterns using a non-homogeneous HMM. The current model was trained on high-coverage healthy samples and applied to low-coverage samples via a pre-trained model for Viterbi decoding only. Three problems arise at low coverage (~1X), particularly for disease samples:

1. **Coverage feature degeneracy:** At 1X, per-CpG coverage is essentially binary (0 or 1 fragment), making the coverage Z-score uninformative and noisy.
2. **Emission distribution mismatch:** Disease samples (cancer, MS) have altered fragmentation patterns relative to healthy — different nucleosome positioning, different nuclease activity (DNASE1L3, DFFB) — so the healthy-trained emission GMM parameters are misspecified.
3. **Downstream TOO estimation:** The current pipeline aggregates per-CpG binary predictions into window-level continuous methylation, then applies QP deconvolution. At 1X, most windows have very few observations, producing extremely noisy methylation estimates and unreliable TOO fractions.

### 1.2 Goals

- Improve per-CpG methylation prediction by adding 5′ fragment end motif as a new feature.
- Build a coverage-free model that transfers robustly to low-coverage data.
- Implement sample-specific emission adaptation so the model works for non-healthy samples.
- Replace window-level QP deconvolution with fragment-level EM deconvolution that preserves within-fragment co-methylation patterns, accounts for missing/unknown cell types, and provides statistical confidence measures (bootstrap CI, p-values, q-values).

### 1.3 Scope

This document covers Stage 1 only. Stage 2 (direct fragment → TOO supervised model) is deferred to a separate document.

---

## 2. Architecture Overview

```
┌─────────────────────────────────────────────────────────────────────────┐
│                        TRAINING PHASE (High Coverage, Healthy)         │
│                                                                        │
│  BAM/frag.gz + hg19.2bit ──► CpgFeatureMatrixBuilder                  │
│                               (length, dist-to-center,                │
│                                5' end motif from .2bit)               │
│                                    │                                   │
│                                    ▼                                   │
│                            Feature Matrix                              │
│                            (N fragments × K CpGs × 3 features)       │
│                                    │                                   │
│                                    ▼                                   │
│                            FinaleMe.java (Baum-Welch)                 │
│                            - GMM initialization                       │
│                            - Full Baum-Welch (emissions + transitions)│
│                            - Converge by KL divergence                │
│                                    │                                   │
│                                    ▼                                   │
│                            Reference Model File                       │
│                            - Emission GMM params (3 feat × 2 states) │
│                            - Transition matrices (by inter-CpG dist) │
│                            - Initiation probabilities                 │
│                            - Motif score lookup table                 │
└────────────────────────────────────┬────────────────────────────────────┘
                                     │
┌────────────────────────────────────▼────────────────────────────────────┐
│                   ADAPTATION PHASE (Per Sample, Any Coverage)          │
│                                                                        │
│  BAM/frag.gz + hg19.2bit ──► CpgFeatureMatrixBuilder ──► Features     │
│                                    │                                   │
│                                    ▼                                   │
│                            FinaleMe.java (-adaptEmissionOnly)         │
│                            - Load reference model                     │
│                            - Freeze transitions + initiation          │
│                            - Constrained Baum-Welch on emissions only │
│                            - Regularize toward reference (λ-shrinkage)│
│                            - 3-5 iterations max                       │
│                                    │                                   │
│                                    ▼                                   │
│                            Adapted Model + pat.gz output              │
│                            (per-fragment per-CpG methylation states)  │
└────────────────────────────────────┬────────────────────────────────────┘
                                     │
┌────────────────────────────────────▼────────────────────────────────────┐
│              REFERENCE PANEL GENERATION (once per panel)               │
│                                                                        │
│  Reference WGBS pat.gz files ──► generate_reference_panel.py          │
│  (per cell type)                 - Use marker regions from            │
│                                    generate_cgi_shore_markers.py      │
│                                  - Build per-CpG reference matrix     │
│                                    │                                   │
│                                    ▼                                   │
│                            reference_panel.tsv                        │
│                            (CpG × cell type methylation matrix)       │
└────────────────────────────────────┬────────────────────────────────────┘
                                     │
┌────────────────────────────────────▼────────────────────────────────────┐
│                   DECONVOLUTION PHASE (Per Sample)                     │
│                                                                        │
│  pat.gz (decoded) + reference_panel.tsv                               │
│         │                                                              │
│         ▼                                                              │
│  fragment_too_deconvolution.py                                        │
│  - Filter informative fragments                                       │
│  - Per-fragment likelihood vs. reference panel                        │
│  - EM for mixture proportions π (including π_unknown)                 │
│  - Bootstrap CI + LRT p-values + BH q-values                         │
│         │                                                              │
│         ▼                                                              │
│  TOO Estimates per Sample                                             │
│  - π_t ± 95% CI for each cell type + unknown fraction                │
│  - p-value and q-value per cell type                                  │
└─────────────────────────────────────────────────────────────────────────┘
```

---

## 3. Component Specifications

### 3.1 Feature Engineering: Add 5′ End Motif

**File to modify:** `edu.northwestern.epifluidlab.finaleme.utils.CpgFeatureMatrixBuilder`

#### 3.1.1 Motif Extraction

- Extract the 4-mer sequence at the **5′ end only** of each fragment.
- The sequence is obtained from the `.2bit` reference genome file (already the 1st positional argument to `CpgFeatureMatrixBuilder`), NOT from the BAM read sequence. This ensures:
  - Compatibility with `frag.gz` input files (which lack sequence data).
  - Consistent sequence regardless of read quality or soft-clipping.
- Use the fragment's mapped 5′ genomic coordinate to look up the 4-mer from the `.2bit` file. For forward-strand fragments: `[start, start+4)`. For reverse-strand fragments: reverse complement of `[end-4, end)`.
- There are 256 possible 4-mers. Fragments with Ns in the 4-mer are assigned a default motif score (see below).

#### 3.1.2 Motif Score Computation

Convert the categorical 4-mer into a continuous scalar for the GMM emission model:

**Training phase (high-coverage WGBS available):**
- For each 4-mer *s*, compute the empirical methylation rate across all CpGs in all fragments with that 5′ motif:
  ```
  motif_score(s) = mean methylation of CpGs in fragments with 5' motif s
  ```
- This produces a lookup table of 256 motif → score mappings.
- Save this lookup table via `-saveMotifLookup`.

**Training phase (WGS only, no WGBS ground truth):**
- Use the GMM initialization (Step 1 of current model) to assign provisional methylated/unmethylated labels.
- Compute motif scores from these provisional labels.
- The scores refine during Baum-Welch iterations.

**Decoding phase (via `-loadMotifLookup`):**
- Load the motif score lookup table from a previously saved TSV file.
- Each fragment's motif score is looked up by its 5′ 4-mer and Z-score normalized along with the other features.

#### 3.1.3 Feature Vector

The per-CpG feature vector becomes 3-dimensional (previously 4D with coverage, now 3D without):

| Feature | Scope | Description |
|---------|-------|-------------|
| Fragment length | Per-fragment | Total length of the cfDNA fragment |
| Distance to center | Per-CpG | Signed distance from CpG to fragment midpoint |
| 5′ end motif score | Per-fragment | Continuous score derived from 5′ 4-mer lookup |

All three features are Z-score normalized within each BAM/frag.gz file using the per-file mean and standard deviation, consistent with the existing normalization approach.

#### 3.1.4 .2bit File Access

The `.2bit` reference genome is already the 1st positional argument to `CpgFeatureMatrixBuilder`. No new parameter is needed. The existing `org.biojava.nbio.genome.parsers.twobit.TwoBitParser` is already used for sequence lookups (via `binRefParser.loadFragment(offset, length)` in the per-thread processing loop). The motif extraction reuses this same parser instance.

Implementation notes:
- The per-thread `binRefParser` is already set to the current chromosome via `setCurrentSequence(chr)`.
- `loadFragment(offset, length)` uses 0-based offset from the start of the current sequence.
- Fragment coordinates from BAM: `fragStart = Math.min(r.getAlignmentStart(), r.getAlignmentStart() + r.getInferredInsertSize())`.
- For tabix fragment input: start/end come directly from BED columns.

#### 3.1.5 Input/Output Changes to CpgFeatureMatrixBuilder

**New flags:**
- `-useEndMotif`: flag to enable 5′ end motif extraction (default: false for backward compatibility).
- `-noCoverage`: flag to drop coverage from the feature matrix.
- `-saveMotifLookup <path>`: save the 5′ end motif score lookup table (256 4-mer → score mappings) to a TSV file during training.
- `-loadMotifLookup <path>`: load a previously saved motif score lookup table for decoding mode.

**Modified outputs:**
- When `-noCoverage` is set, the coverage column is omitted from the feature matrix.
- When `-useEndMotif` is set, the motif score column is added.
- The column order becomes: fragment_length, distance_to_center, motif_score.

---

### 3.2 Coverage-Free Reference Model Training

**File to modify:** `FinaleMe.java`
(`edu.northwestern.epifluidlab.finaleme.hmm.FinaleMe`)

**Existing USAGE:** `FinaleMe [opts] model input_matrix.txt[.gz] prediction.txt.gz`

**Positional arguments (3, unchanged):**
1. `model` — HMM model file path. Used for both saving (training) and loading (decode). Serialized via Java `ObjectOutputStream` as `BayesianNhmmV5<ObservationVector>`.
2. `input_matrix.txt[.gz]` — feature matrix from `CpgFeatureMatrixBuilder`.
3. `prediction.txt.gz` — output decoded predictions.

**Existing relevant flags:**
- `-lowCoverage` — **already drops coverage** and uses only 2 features (fragLen, DistToCenter). This is the existing mechanism for what we called `-noCoverage`.
- `-decodeModeOnly` — decode without training, loading model from 1st positional arg.
- `-patOutput` — output UXM-compatible `.pat.gz` and `.beta` files. **Already exists.**
- `-cpgIndexFile` — CpG index file required with `-patOutput`. **Already exists.**
- `-bwOutput` — output bigWig files. **Already exists.**
- `-features` — number of features per observation (default: 3).
- `-miniDataPoints` — minimum CpGs per fragment (default: 1; set to 7 for training).
- `-gmm` — use GMM initialization.

#### 3.2.1 New Parameters

| Parameter | Type | Description |
|-----------|------|-------------|
| `-useEndMotif` | Flag | Include 5′ end motif score in emission features |
| `-adaptEmissionOnly` | Flag | Constrained Baum-Welch: freeze transitions, adapt emissions only |
| `-adaptLambda` | Float | Shrinkage regularization toward reference model (default: 0.5) |
| `-adaptMaxIter` | Int | Max iterations for adaptation (default: 5) |

When `-lowCoverage` and `-useEndMotif` are both set, the feature vector becomes: (fragLen, DistToCenter, motif_score) — 3 features. Set `-features 3` accordingly. Without `-useEndMotif`, `-lowCoverage` uses 2 features (fragLen, DistToCenter) as it does currently.

The model file (1st positional arg) handles all HMM parameter persistence (emissions, transitions, initiation) via Java serialization. The motif lookup table is handled separately by `CpgFeatureMatrixBuilder` via `-saveMotifLookup`/`-loadMotifLookup`.

#### 3.2.2 Training Procedure

1. **Input:** High-coverage healthy cfDNA WGS BAM/frag.gz + `.2bit` reference genome.
2. **Feature extraction:** `CpgFeatureMatrixBuilder` with `-useEndMotif -noCoverage -wgsMode -stringentPaired -saveMotifLookup motif_lookup.tsv`.
3. **GMM initialization:** Classify all CpG points into 2 (or 3 with motif) clusters using GMM (flag `-gmm`). Methylated state identified by larger mean distance-to-center.
4. **Baum-Welch:** Full estimation of emissions + transitions + initiation. Fragments with ≥7 CpGs (`-miniDataPoints 7`). Max 50 iterations or KL convergence < 1e-4.
5. **Output:** Model file (1st positional arg) serialized as `BayesianNhmmV5<ObservationVector>` containing all HMM parameters.

#### 3.2.3 Validation Experiments

| Experiment | Data | Metric | Purpose |
|------------|------|--------|---------|
| A | High-cov healthy WGS vs. matched WGBS | auROC (per-CpG binary), Pearson r (1kb windows in CGI/shore) | Measure impact of adding motif / dropping coverage at high coverage |
| B | High-cov healthy WGS downsampled to 1X vs. matched WGBS | Same as A | Confirm 3-feature model outperforms 4-feature at low coverage |
| C | High-cov cancer WGS vs. matched WGBS | Same as A | Baseline for cancer before adaptation (Step 3.3) |

---

### 3.3 Sample-Specific Emission Adaptation

**File to modify:** `FinaleMe.java`

#### 3.3.1 New Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `-adaptEmissionOnly` | Flag | false | Enable constrained Baum-Welch (freeze transitions) |
| `-adaptLambda` | Float | 0.5 | Shrinkage regularization toward reference model (0 = no regularization, 1 = no adaptation) |
| `-adaptMaxIter` | Int | 5 | Maximum Baum-Welch iterations during adaptation |
| `-adaptMinFragments` | Int | 1000 | Minimum fragments with ≥7 CpGs required; if not met, skip adaptation and use reference model directly |

#### 3.3.2 Constrained Baum-Welch Algorithm

```
Input:
  - Reference model M_ref (from Step 3.2, loaded from 1st positional arg model file)
  - Target sample BAM/frag.gz
  - Regularization parameter λ
  - Max iterations T

Algorithm:
  1. Initialize M_current = M_ref
  2. Extract features from target BAM using CpgFeatureMatrixBuilder
     - Z-score normalize using TARGET sample's mean/std (not reference)
  3. For iter = 1 to T:
     a. E-step: compute forward-backward probabilities using M_current
        - Transition matrices = M_ref.transitions (FROZEN)
        - Initiation probs = M_ref.initiation (FROZEN)
        - Emission probs = M_current.emissions
     b. M-step (emissions only):
        - Compute MLE emission parameters θ_MLE from expected sufficient statistics
        - Regularize: θ_new = (1 - λ) · θ_MLE + λ · θ_ref
          For GMM parameters specifically:
            μ_new = (1-λ)·μ_MLE + λ·μ_ref
            Σ_new = (1-λ)·Σ_MLE + λ·Σ_ref
            π_new = (1-λ)·π_MLE + λ·π_ref
     c. Compute KL divergence between M_current and M_new
        - If KL < 1e-4 or ΔKL < 1%: converge, break
     d. M_current.emissions = θ_new
  4. Viterbi decode all fragments in target sample using M_current
  5. Output pat.gz file with per-fragment per-CpG methylation states

Output:
  - pat.gz file: per-CpG binary methylation state for each fragment
  - Adapted model parameters (optional, for inspection/diagnostics)
```

#### 3.3.3 Z-Score Normalization: Critical Detail

The Z-score normalization MUST use the target sample's own feature statistics, not the reference model's. This is because:
- Fragment length distributions can shift between healthy and disease.
- Motif score distributions reflect the target sample's biology.
- Using reference normalization would introduce systematic bias.

However, the GMM parameters from the reference model were estimated under the reference's normalization. Therefore, at adaptation initialization, we must re-center the reference GMM means to account for the normalization shift:

```
μ_ref_adjusted = (μ_ref * σ_ref + mean_ref - mean_target) / σ_target
```

Or more practically: re-normalize the reference GMM parameters into the target sample's Z-score space before starting adaptation.

#### 3.3.4 Handling Very Low Fragment Counts

At 1X WGS, a sample has ~20M total fragments. Fragments with ≥7 CpGs (required for training) will be concentrated in CpG-dense regions (CGIs, shores). Expected count: roughly 50K–200K at 1X, which is sufficient for GMM adaptation with regularization.

If a sample has fewer than `-adaptMinFragments` qualifying fragments (e.g., due to very low coverage or poor quality), skip adaptation and use the reference model directly. Log a warning.

#### 3.3.5 pat.gz Output

FinaleMe already supports `pat.gz` output format, which contains per-fragment per-CpG binary methylation states. This is the native output format for the Viterbi decoding step and serves as the input to the fragment-level TOO deconvolution module (Section 3.5).

No new output flag is needed — the existing pat.gz output pipeline is reused.

#### 3.3.6 Validation Experiments

| Experiment | Data | Comparison | Metric |
|------------|------|------------|--------|
| D | High-cov cancer WGS (full) | Reference-only vs. adapted model vs. full Baum-Welch | auROC, Pearson r vs. WGBS |
| E | High-cov cancer WGS downsampled to 1X | Same as D | Same — expect adapted > reference-only, adapted ≈ full BW |
| F | 1X cancer WGS (real ULP-WGS samples) | Reference-only vs. adapted | TOO correlation with ichorCNA TF |
| G | Vary λ from 0.0 to 1.0 in 0.1 steps | Adapted model at 1X | Pearson r vs. WGBS — find optimal λ |
| H | Vary adaptMaxIter from 1 to 10 | Adapted model at 1X | Pearson r vs. WGBS — confirm 3-5 is sufficient |

---

### 3.4 Reference Panel Generation

**New script:** `generate_reference_panel.py` (extends or complements existing `generate_cgi_shore_markers.py`)

#### 3.4.1 Purpose

Build a per-CpG reference methylation matrix from cell-type-specific WGBS data, restricted to informative marker regions (CGI and CGI shore), for use in fragment-level TOO deconvolution.

#### 3.4.2 Inputs

1. **Cell-type-specific WGBS pat.gz files:** A collection of pat.gz files, one per reference cell type (e.g., Neutrophil, B cell, T cell, Macrophage, Erythroblast, Endothelial, Liver hepatocyte, Oligodendrocyte, Neuron, etc.). Each pat.gz contains per-fragment per-CpG binary methylation states from bulk or sorted-cell WGBS.

2. **Marker region definition:** Cell-type-specific marker regions in CGI and CGI shore, generated by `generate_cgi_shore_markers.py` or a similar script. These are genomic intervals where at least one cell type has a methylation level that is sufficiently distinct from other cell types.

3. **CpG position file:** The standard CG_motif bedgraph file used throughout FinaleMe.

#### 3.4.3 Processing Pipeline

```
Step 1: For each cell type's pat.gz file:
  - Aggregate per-fragment binary methylation states at each CpG site
  - Compute continuous methylation level: β_t(c) = methylated_count / total_count
  - Apply minimum coverage filter: require ≥ N reads at each CpG
    (e.g., N=5 for high-cov, N=3 for low-cov references)

Step 2: Merge into CpG × cell type matrix:
  - Rows: CpG positions present in all (or most) reference cell types
  - Columns: cell types
  - Values: methylation level β_t(c), clipped to [0.01, 0.99]

Step 3: Filter to marker regions:
  - Restrict to CpGs within CGI and CGI shore regions
  - Use generate_cgi_shore_markers.py to identify cell-type-specific DMRs
    OR apply a variance filter: retain CpGs where the variance across
    cell types exceeds a threshold (e.g., top 1% or top 25% most variable,
    depending on coverage regime — consistent with existing TOO pipeline)

Step 4: Output reference panel file
```

#### 3.4.4 Output Format

```
reference_panel.tsv:

chr    position    Neutrophil    B_cell    T_cell    Macrophage    Erythroblast    ...
chr1   10023       0.92          0.05      0.08      0.88          0.91            ...
chr1   10041       0.95          0.08      0.11      0.91          0.93            ...
...
```

#### 3.4.5 Input File Manifest Format

Both `--pat-files` and `--bigwig-files` accept a **manifest file** (TSV) listing input files with their associated metadata. This allows multiple files from the same cell type (biological replicates, different donors) to be merged automatically.

**pat.gz manifest format** (`--pat-files manifest.tsv`):

```
# name          group           file_path
Neutrophil_1    Neutrophil      /data/ref/neutrophil_donor1.pat.gz
Neutrophil_2    Neutrophil      /data/ref/neutrophil_donor2.pat.gz
Bcell_1         B_cell          /data/ref/bcell_donor1.pat.gz
Tcell_1         T_cell          /data/ref/tcell_donor1.pat.gz
Tcell_2         T_cell          /data/ref/tcell_donor2.pat.gz
Oligo_1         Oligodendrocyte /data/ref/oligo_donor1.pat.gz
Neuron_1        Neuron          /data/ref/neuron_donor1.pat.gz
...
```

- **name**: unique identifier for each file (used in logs/diagnostics).
- **group**: cell type label. Files sharing the same group are merged (CpG-level methylated/total counts are pooled) before computing per-CpG methylation levels.
- **file_path**: path to pat.gz file.

**bigWig manifest format** (`--bigwig-files manifest.tsv`):

```
# name          group           methy_bw                            cov_bw
Neutrophil_1    Neutrophil      /data/ref/neutrophil_1.methy.bw     /data/ref/neutrophil_1.cov.bw
Neutrophil_2    Neutrophil      /data/ref/neutrophil_2.methy.bw     /data/ref/neutrophil_2.cov.bw
Bcell_1         B_cell          /data/ref/bcell_1.methy.bw          /data/ref/bcell_1.cov.bw
...
```

- **methy_bw**: bigWig file with methylation count (or methylation level) per CpG.
- **cov_bw**: bigWig file with coverage count per CpG.
- Files sharing the same group are merged by summing methylation counts and coverage counts before computing β values.

#### 3.4.6 Command-Line Interface

```bash
# Using pat.gz inputs (preferred for new panels)
python generate_reference_panel.py \
  --pat-files pat_manifest.tsv \
  --marker-regions cgi_shore_markers.bed \
  --cpg-positions CG_motif.hg19.pos_only.bedgraph \
  --min-coverage 5 \
  --variance-percentile 0.01 \
  --output reference_panel.tsv

# Using bigWig inputs (compatible with existing pipeline)
python generate_reference_panel.py \
  --bigwig-files bigwig_manifest.tsv \
  --marker-regions cgi_shore_markers.bed \
  --cpg-positions CG_motif.hg19.pos_only.bedgraph \
  --min-coverage 5 \
  --variance-percentile 0.01 \
  --output reference_panel.tsv
```

#### 3.4.7 Compatibility with Existing Panels

The script should produce panels compatible with the existing reference methylomes used in FinaleMe's TOO analysis (Moss 2018 panel, New Oligo, Old Oligo). Users can generate custom panels for new applications (e.g., MS-specific panels with neuronal subtypes) using the same pipeline.

---

### 3.5 Fragment-Level TOO Deconvolution

**New module:** `fragment_too_deconvolution.py` (Python)

#### 3.5.1 Input Format

The input is the pat.gz file produced by FinaleMe's Viterbi decoding step (Section 3.3.5). The pat.gz format contains per-fragment per-CpG binary methylation states. The deconvolution module parses this format natively.

#### 3.5.2 Reference Panel

The reference panel from Section 3.4 (reference_panel.tsv). Only CpGs present in BOTH the reference panel and the decoded fragments are used.

#### 3.5.3 Fragment Filtering

Not all fragments are informative for deconvolution. Apply the following filters:

1. **Minimum CpG count:** Fragment must cover ≥3 CpGs present in the reference panel.
2. **Reference panel coverage:** All CpGs in the fragment must have non-NA values in the reference panel.
3. **Informativeness filter:** At least 2 reference cell types must have meaningfully different methylation at the fragment's CpGs. Compute per-fragment informativeness score:
   ```
   I(f) = max_{t1, t2} Σ_i |β_t1(c_i) - β_t2(c_i)| / k
   ```
   where the max is over all cell type pairs and k is the number of CpGs. Discard fragments with I(f) < 0.2 (i.e., all cell types look similar at those CpGs).

#### 3.5.4 Per-Fragment Likelihood Computation

For each fragment f with predicted methylation vector **m** = (m₁, ..., mₖ) and each reference cell type t:

```
log P(f | t) = Σ_i [ m_i · log(β_t(c_i)) + (1 - m_i) · log(1 - β_t(c_i)) ]
```

Clip β_t(c_i) to [0.01, 0.99] to avoid log(0).

Store as an N × T matrix of log-likelihoods, where N = number of informative fragments, T = number of cell types.

#### 3.5.5 Unknown/Missing Cell Type Estimation

The reference panel will never contain all cell types contributing to cfDNA. Fragments originating from cell types not in the panel will be forced into the closest available cell type, biasing proportion estimates. To address this, we add an explicit **unknown cell type** component to the mixture model.

**Uniform background model:**
- Define P(f | unknown) as a uniform likelihood over all possible methylation patterns of length k:
  ```
  log P(f | unknown) = k · log(0.5)
  ```
  This assumes that an unknown cell type produces methylation patterns indistinguishable from random coin flips at the marker CpGs.

**EM with unknown component:**
- The mixture model becomes: P(f) = Σ_t π_t · P(f | t) + π_unknown · P(f | unknown)
- The EM algorithm estimates π_unknown alongside all other proportions, subject to the constraint Σ_t π_t + π_unknown = 1.
- π_unknown captures the fraction of cfDNA that cannot be explained by any reference cell type.

**Interpretation:**
- A high π_unknown (e.g., >0.2) indicates either missing cell types in the reference panel or poor methylation prediction quality.
- π_unknown is reported with bootstrap CI, p-value, and q-value alongside known cell types.
- For downstream analysis, known cell type proportions can be reported both as raw (summing to 1 - π_unknown) and renormalized (summing to 1 after removing unknown).

#### 3.5.6 EM Algorithm for Mixture Proportions

```
Input:
  - Log-likelihood matrix L (N × (T+1)): L[j, t] = log P(f_j | t)
    where t = 1..T are known cell types and t = T+1 is the unknown component
  - Number of cell types T + 1 (including unknown)
  - Convergence threshold ε = 1e-6
  - Maximum iterations = 1000

Algorithm:
  1. Initialize π uniformly: π_t = 1/(T+1) for all t including unknown
  2. Repeat:
     a. E-step: compute responsibilities
        log_numerator[j, t] = log(π_t) + L[j, t]
        γ[j, t] = exp(log_numerator[j, t] - logsumexp(log_numerator[j, :]))
     b. M-step: update proportions
        π_t = Σ_j γ[j, t] / N
     c. Compute log-likelihood:
        LL = Σ_j logsumexp(log_numerator[j, :])
     d. If |ΔLL| < ε: converge, break
  3. Return π, γ, LL

Output:
  - π: estimated cell type proportions ((T+1)-vector, including π_unknown)
  - γ: per-fragment responsibilities (N × (T+1) matrix)
  - LL: final log-likelihood
```

Use logsumexp throughout to maintain numerical stability.

#### 3.5.7 Confidence Estimation: Bootstrap CI

```
Input:
  - Log-likelihood matrix L (N × (T+1))
  - Number of bootstrap replicates B = 1000

Algorithm:
  For b = 1 to B:
    1. Sample N fragment indices with replacement
    2. Run EM on the resampled log-likelihood matrix
    3. Store π_b

Output:
  - For each cell type t (including unknown):
    - π_t: point estimate (from full data EM)
    - 95% CI: [percentile_2.5(π_b_t), percentile_97.5(π_b_t)]
    - Bootstrap standard error: std(π_b_t)
```

#### 3.5.8 Statistical Testing: P-values and Q-values

**Goal:** For each cell type (including unknown), test whether its estimated proportion is significantly greater than zero (i.e., is this tissue type actually contributing to cfDNA?).

**Null hypothesis:** H₀: π_t = 0 (cell type t does not contribute to cfDNA).

**Approach — Likelihood ratio test (default):**

```
For each cell type t:
  1. Full model: EM with all T+1 components → LL_full
  2. Reduced model: EM with cell type t removed (T components) → LL_reduced
  3. LRT statistic: Λ = 2 · (LL_full - LL_reduced)
  4. P-value: from χ²(df=1) distribution (asymptotic)
```

**Alternative approach — Permutation-based p-value (validation / rare cell types):**

```
For each cell type t:
  For p = 1 to P (P=1000):
    1. Permute the rows of column t in L independently
       (breaks the association between fragment methylation patterns
        and cell type t's reference, while preserving other cell types)
    2. Re-run EM on the permuted matrix
    3. Store π_t_null[p]
  p_value[t] = (number of π_t_null[p] ≥ π_t_observed + 1) / (P + 1)
```

The LRT approach is the default for speed. For rare cell types where π_t is near the boundary (0), the χ² approximation can be poor. Use permutation-based p-values as a fallback or validation for cell types with LRT p-values near the significance threshold.

**FDR correction:**
- Apply Benjamini-Hochberg procedure across all T+1 components (including unknown) within each sample.
- Report q-values alongside p-values.

**Output per sample:**

```
cell_type    proportion    proportion_renorm    CI_lower    CI_upper    p_value    q_value    significant
Neutrophil   0.58          0.64                 0.54        0.62        <0.001     <0.001     yes
B_cell       0.07          0.08                 0.04        0.10        <0.001     <0.001     yes
T_cell       0.14          0.15                 0.11        0.17        <0.001     <0.001     yes
Oligodendro  0.03          0.03                 0.01        0.05        0.012      0.024      yes
Neuron       0.005         0.006                0.00        0.02        0.340      0.425      no
Unknown      0.09          —                    0.05        0.13        0.002      0.005      yes
...
```

Where `proportion` sums to 1.0 across all rows (including unknown), and `proportion_renorm` is rescaled to sum to 1.0 across known cell types only (unknown excluded).

#### 3.5.9 Implementation Details

**Dependencies:** numpy, scipy, pandas, statsmodels (for FDR), joblib (for parallel bootstrap/permutation), pysam or tabix (for pat.gz reading).

**Performance considerations:**
- At 1X, after filtering, expect ~100K–1M informative fragments. The log-likelihood matrix (N × (T+1)) fits in memory for T ≤ 50 cell types.
- EM converges fast (typically <100 iterations). Each iteration is O(N × T).
- Bootstrap is embarrassingly parallel. With B=1000 and joblib, runtime is ~5-10 minutes per sample on a single node.
- LRT requires T+1 additional EM runs (one per cell type removed). Fast and deterministic.
- Permutation-based p-values (if used) add P × T EM runs. With P=1000 and T=10, this is 10,000 EM runs — still feasible (~30 min) but LRT is preferred for speed.

**Command-line interface:**

```bash
python fragment_too_deconvolution.py \
  --pat-gz decoded_sample.pat.gz \
  --reference-panel reference_panel.tsv \
  --min-cpgs 3 \
  --informativeness-threshold 0.2 \
  --estimate-unknown \
  --bootstrap-replicates 1000 \
  --test-method lrt \
  --permutation-replicates 1000 \
  --fdr-threshold 0.05 \
  --output-prefix sample_001_too \
  --n-jobs 8
```

**Output files:**
- `{prefix}_proportions.tsv`: cell type proportions (raw + renormalized) with CI, p-values, q-values. Includes the unknown component.
- `{prefix}_fragment_responsibilities.tsv.gz`: per-fragment γ matrix (for downstream analysis or Stage 2 training data).
- `{prefix}_diagnostics.json`: EM convergence info, number of fragments used, number filtered, informativeness distribution, log-likelihood, π_unknown estimate.

---

## 4. Data Flow Summary

```
                    REFERENCE MODEL TRAINING (once, high-coverage healthy)
                    ======================================================

   healthy.bam/frag.gz + hg19.2bit
         │
         ▼
   CpgFeatureMatrixBuilder
   hg19.2bit cpg_list.bed all_cpg.bed healthy.bam cpg_detail.txt.gz
   -useEndMotif -noCoverage -wgsMode -stringentPaired
   -saveMotifLookup motif_lookup.tsv
         │
         ├──► cpg_detail.txt.gz (length, dist-to-center, motif_score)
         └──► motif_lookup.tsv (256 × 1)
               │
               ▼
         FinaleMe.java
         reference_model.ser cpg_detail.txt.gz prediction.txt.gz
         -miniDataPoints 7 -gmm -lowCoverage -useEndMotif
         -patOutput -cpgIndexFile CG_motif.bedgraph
               │
               ▼
         reference_model.ser
         (GMM params, transitions, initiation, motif lookup, norm stats)


                    REFERENCE PANEL GENERATION (once per cell type panel)
                    =====================================================

   Cell-type WGBS pat.gz files + CGI/shore marker regions
         │
         ▼
   generate_reference_panel.py
   --pat-files *.pat.gz
   --marker-regions cgi_shore_markers.bed
   (marker regions from generate_cgi_shore_markers.py)
         │
         ▼
   reference_panel.tsv
   (CpG × cell type methylation matrix, filtered to markers)


                    PER-SAMPLE INFERENCE (each target sample)
                    =========================================

   sample.bam/frag.gz + hg19.2bit
         │
         ▼
   CpgFeatureMatrixBuilder
   hg19.2bit cpg_list.bed all_cpg.bed sample.bam cpg_detail.txt.gz
   -useEndMotif -noCoverage -wgsMode -stringentPaired
   -loadMotifLookup motif_lookup.tsv
         │
         ▼
   FinaleMe.java
   reference_model.ser cpg_detail.txt.gz prediction.txt.gz
   -decodeModeOnly -adaptEmissionOnly -adaptLambda 0.5 -adaptMaxIter 5
   -lowCoverage -useEndMotif
   -patOutput -cpgIndexFile CG_motif.bedgraph
         │
         ├──► adapted_model.ser (optional, diagnostics)
         └──► sample.pat.gz (per-fragment per-CpG methylation)
               │
               ▼
         fragment_too_deconvolution.py
         --pat-gz sample.pat.gz
         --reference-panel reference_panel.tsv
         --estimate-unknown
         --bootstrap-replicates 1000 --test-method lrt
               │
               ├──► sample_too_proportions.tsv
               ├──► sample_fragment_responsibilities.tsv.gz
               └──► sample_diagnostics.json
```

---

## 5. Java Code Modification Details

### 5.1 CpgFeatureMatrixBuilder.java

**Package:** `edu.northwestern.epifluidlab.finaleme.utils`

**USAGE string:** `CpgFeatureMatrixBuilder [opts] hg19.2bit cpg_list.bed all_cpg.bed wgs.bam|fragments.tsv.gz cpg_detail.txt.gz`

**Positional arguments (5, unchanged):**
1. `hg19.2bit` — reference genome (refFile), parsed by `org.biojava.nbio.genome.parsers.twobit.TwoBitParser`
2. `cpg_list.bed` — CpG coordinate bedgraph (cpgListFile)
3. `all_cpg.bed` — all CpG positions for distance calculations (allCpgFile)
4. `wgs.bam|fragments.tsv.gz` — input BAM or tabix-indexed fragment file (wgsBamFile); frag.gz already supported via `-fragmentInputTabix` auto-detection
5. `cpg_detail.txt.gz` — output gzipped feature matrix (detailFile)

**Existing output columns:** `chr, start, end, readName, FragLen, Frag_strand, methy_stat, Norm_Frag_cov, baseQ, Offset_frag, Dist_frag_end` + optional tracks (overlapRegions, distantRegions, valueBeds, valueWigs, kmers).

**Existing .2bit usage:** The `TwoBitParser` (biojava) is already instantiated per-thread as `binRefParser` and used for sequence lookups via `binRefParser.loadFragment(offset, length)`. The parser is set to the current chromosome via `binRefParser.setCurrentSequence(chr)`. This infrastructure is reused for motif extraction — no new dependency needed.

**Existing kmer infrastructure:** The class already has `-kmerLen`, `-kmerString`, `-useFragBaseKmer`, and `CcInferenceUtils.kmerFreqSearch()` for computing k-mer frequencies in reference or fragment sequences. The end motif feature is conceptually different (a single 4-mer at the 5' end mapped to a score) but can follow similar patterns.

**No existing `-saveModel`/`-loadModel`:** These flags do not exist in `CpgFeatureMatrixBuilder`. Model save/load is handled by `FinaleMe.java`. The motif lookup table should be persisted via dedicated flags (`-saveMotifLookup`/`-loadMotifLookup`) since it is computed during feature extraction.

**New @Option fields:**
```java
@Option(name="-useEndMotif", usage="add 5' end 4-mer motif score as a feature. Default: false")
public boolean useEndMotif = false;

@Option(name="-noCoverage", usage="drop normalized coverage from features. Default: false")
public boolean noCoverage = false;

@Option(name="-saveMotifLookup", usage="save motif score lookup table (256 entries) to file. Default: null")
public String saveMotifLookup = null;

@Option(name="-loadMotifLookup", usage="load motif score lookup table from file. Default: null")
public String loadMotifLookup = null;
```

**New class member:**
```java
private Map<String, Double> motifScoreMap;  // 256 4-mer -> methylation rate, built or loaded
```

**New method: `extractFivePrimeMotif()`**

Uses the existing per-thread `binRefParser` (TwoBitParser) already set to the current chromosome:

```java
/**
 * Extract 4-mer at the 5' end of the fragment from the .2bit reference.
 *
 * Uses the existing binRefParser (TwoBitParser) already set to the
 * current chromosome in the per-thread processing loop.
 *
 * For + strand fragments: 5' end is at fragStart, extract [fragStart, fragStart+4)
 * For - strand fragments: 5' end is at fragEnd, extract reverse complement of [fragEnd-4, fragEnd)
 *
 * Fragment coordinates are derived from SAMRecord:
 *   fragStart = Math.min(r.getAlignmentStart(), r.getAlignmentStart() + r.getInferredInsertSize())
 *   fragEnd = Math.max(r.getAlignmentStart(), r.getAlignmentStart() + r.getInferredInsertSize())
 * For tabix fragment input: directly from the BED start/end columns.
 *
 * @param binRefParser  TwoBitParser set to current chromosome
 * @param fragStart     0-based fragment start (genomic)
 * @param fragEnd       0-based fragment end (genomic)
 * @param negStrand     whether fragment is on negative strand
 * @return 4-mer string (uppercase), or "NNNN" if lookup fails or contains N
 */
private String extractFivePrimeMotif(TwoBitParser binRefParser, int fragStart, int fragEnd, boolean negStrand) {
    try {
        String seq;
        if (!negStrand) {
            // + strand: 5' is at fragStart
            seq = binRefParser.loadFragment(fragStart - 1, 4);  // loadFragment is 0-based offset from seq start
        } else {
            // - strand: 5' is at fragEnd, take reverse complement
            seq = binRefParser.loadFragment(fragEnd - 4, 4);
            byte[] seqBytes = seq.getBytes();
            SequenceUtil.reverseComplement(seqBytes);
            seq = new String(seqBytes);
        }
        seq = seq.toUpperCase();
        return seq.contains("N") ? "NNNN" : seq;
    } catch (Exception e) {
        return "NNNN";
    }
}
```

**Integration point in the per-CpG-per-fragment output loop:**

In the existing code, after computing `fragLen`, `cpgOffset`, `distToFragEnd`, and `normalizedFragCov`, and before the `binWriter.write(...)` call:

```java
// Extract 5' end motif score (new)
double motifScore = 0.0;
if (useEndMotif) {
    String motif = extractFivePrimeMotif(binRefParser, fragStart, fragEnd, negStrand);
    motifScore = motifScoreMap.getOrDefault(motif, 0.5);  // default 0.5 for unknown
}
```

**Modified output line (when `-noCoverage` and `-useEndMotif` are set):**

The `Norm_Frag_cov` column is omitted, and a `motif_score` column is appended:
```java
if (noCoverage) {
    // omit normalizedFragCov from output
} else {
    binWriter.write("\t" + String.format("%.6f", normalizedFragCov));
}
if (useEndMotif) {
    binWriter.write("\t" + String.format("%.6f", motifScore));
}
```

The header line is modified correspondingly.

**Motif score table construction (training mode, `-saveMotifLookup`):**

Two-pass approach within the existing parallel processing:
1. **First pass:** For each fragment-CpG pair, record the 5' 4-mer and the `methy_stat` (m/u). Accumulate counts per 4-mer in a thread-safe `ConcurrentHashMap<String, long[]>` where `long[0]` = methylated count, `long[1]` = total count.
2. **After all bins complete:** Compute `motifScoreMap.put(kmer, (double)methyCount / totalCount)` for each of the 256 4-mers. Write to TSV file.

Alternatively, if a single-pass is preferred: use the GMM provisional labels (from the existing wgsMode logic where `methy_stat` is set to '.' for WGS) — but this requires coordination with FinaleMe.java's GMM initialization. A simpler approach for WGS mode: run CpgFeatureMatrixBuilder without motif first, then after FinaleMe GMM initialization provides provisional labels, compute the motif lookup table and re-run with `-loadMotifLookup`.

**Modified feature output:**
- When `-noCoverage` is set, omit the coverage column from the feature matrix.
- When `-useEndMotif` is set, add motif score column.
- The column order becomes: fragment_length, distance_to_center, motif_score.

**Motif score initialization (training mode, `-saveMotifLookup`):**
- First pass: count methylated and total CpGs per 4-mer (using GMM provisional labels or WGBS ground truth if available).
- Compute motif_score(s) = methylated_count(s) / total_count(s).
- Save motif lookup table and feature normalization stats into model file.

**Motif score application (decode mode, `-loadMotifLookup`):**
- Load motif lookup table from the model file.
- Look up each fragment's 5′ 4-mer and assign the corresponding score.

### 5.2 FinaleMe.java

**Package:** `edu.northwestern.epifluidlab.finaleme.hmm.FinaleMe`

**Existing USAGE:** `FinaleMe [opts] model input_matrix.txt[.gz] prediction.txt.gz`

**Model persistence:** The model (1st positional arg) is a serialized `BayesianNhmmV5<ObservationVector>` object via Java `ObjectOutputStream`. The class includes legacy package remapping (`LegacyPackageObjectInputStream`) to handle models saved under the old `org.cchmc.epifluidlab` package.

**Existing `-lowCoverage` flag:** When set, the feature vector drops coverage and uses only 2D: `{fragLen, DistToCenter}`. This is the mechanism we use for coverage-free models. With the addition of end motif, `-lowCoverage -useEndMotif` produces 3D: `{fragLen, DistToCenter, motifScore}`.

**New @Option fields:**
```java
@Option(name="-useEndMotif", usage="include 5' end motif score as a third feature in lowCoverage mode. Default: false")
public boolean useEndMotif = false;

@Option(name="-adaptEmissionOnly", usage="constrained Baum-Welch: freeze transitions, adapt emissions only. Default: false")
public boolean adaptEmissionOnly = false;

@Option(name="-adaptLambda", usage="shrinkage regularization toward reference model (0=no reg, 1=no adapt). Default: 0.5")
public double adaptLambda = 0.5;

@Option(name="-adaptMaxIter", usage="max Baum-Welch iterations during adaptation. Default: 5")
public int adaptMaxIter = 5;

@Option(name="-adaptMinFragments", usage="minimum fragments with >=miniDataPoints CpGs to attempt adaptation. Default: 1000")
public int adaptMinFragments = 1000;
```

**Key code change in `processMatrixFile()` — add motif feature to observation vector:**

In the existing code, when `-lowCoverage` is set, the value array is:
```java
value = new double[]{
    (fragLen-stats[0].getMean())/stats[0].getStandardDeviation(),
    (DistToCenter-stats[2].getMean())/stats[2].getStandardDeviation(),
};
```

With `-useEndMotif`, this becomes:
```java
if (lowCoverage && useEndMotif) {
    value = new double[]{
        (fragLen-stats[0].getMean())/stats[0].getStandardDeviation(),
        (DistToCenter-stats[2].getMean())/stats[2].getStandardDeviation(),
        (motifScore-stats[3].getMean())/stats[3].getStandardDeviation(),
    };
}
```

This requires `CpgFeatureMatrixBuilder` to output the motif score column in the feature matrix (in place of or after the coverage column), and `parseLine()` to read it.

**Modified Baum-Welch loop (when `-adaptEmissionOnly` is set):**

```java
// Pseudocode for constrained Baum-Welch
HmmModel currentModel = deepCopy(referenceModel);

for (int iter = 0; iter < adaptMaxIter; iter++) {
    // E-step: forward-backward with current emissions + FROZEN transitions
    SufficientStats stats = forwardBackward(fragments, currentModel);

    // M-step: update emissions only
    GmmParams mleEmissions = estimateEmissions(stats);

    // Regularize toward reference
    GmmParams newEmissions = interpolate(mleEmissions, referenceModel.emissions, adaptLambda);
    // interpolate: (1 - λ) * MLE + λ * reference, applied to μ, Σ, π separately

    // Check convergence
    double kl = klDivergence(currentModel, newEmissions);
    currentModel.emissions = newEmissions;
    // transitions and initiation remain referenceModel's values

    if (kl < 1e-4) break;
}

// Viterbi decode with adapted model, output pat.gz
viterbiDecode(allFragments, currentModel);
```

**Model persistence (1st positional arg):**

The model is already serialized via Java `ObjectOutputStream` as `BayesianNhmmV5<ObservationVector>`. No change to the serialization mechanism is needed. For adaptation, the reference model is loaded from the 1st positional arg, adapted in memory, and the adapted model can be saved to a different path for diagnostics.

---

## 6. Validation Plan

### 6.1 Datasets

| Dataset | Coverage | Condition | Matched WGBS | Purpose |
|---------|----------|-----------|--------------|---------|
| HD_45 | ~30X | Healthy | Yes (~15X) | Reference model training + validation |
| Cancer_01 | ~16-39X | Prostate cancer | Yes (~10-15X) | Cancer adaptation validation |
| ULP-WGS cohort | ~0.1X | 12 healthy, 22 breast, 43 prostate | Matched ULP-WGBS | Low-coverage TOO validation |
| Downsampled HD_45 | 1X, 2X, 5X | Healthy | Same as HD_45 | Coverage titration |
| Downsampled Cancer_01 | 1X, 2X, 5X | Cancer | Same as Cancer_01 | Coverage titration + adaptation |

### 6.2 Experiments

#### Experiment 1: End Motif Feature Impact (Section 3.1)

- Train on HD_45 with 3 feature sets: (A) original 3 features (length, dist, coverage), (B) new 3 features (length, dist, motif), (C) 4 features (length, dist, coverage, motif).
- Evaluate per-CpG auROC on held-out fragments (by chromosome).
- Report improvement by CpG density bin (≥5, ≥3, ≥1 CpGs per fragment).
- Expected: motif adds ~1-3% auROC improvement, most in CpG-sparse regions.

#### Experiment 2: Coverage Dropout at Low Coverage (Section 3.2)

- Downsample HD_45 to 1X, 2X, 5X.
- Apply models A, B, C from Experiment 1.
- Compare Pearson r of predicted vs. WGBS methylation at 1kb CGI/shore windows.
- Expected: model B (no coverage) outperforms A (with coverage) at 1X; C performs best at high coverage.

#### Experiment 3: Emission Adaptation for Cancer (Section 3.3)

- On Cancer_01 at full coverage and at 1X downsampled:
  - (i) Reference model (healthy, no adaptation)
  - (ii) Adapted model (λ = 0.5, 5 iterations)
  - (iii) Full Baum-Welch on cancer sample alone
- Compare methylation prediction accuracy vs. matched WGBS.
- λ sweep: 0.0, 0.1, 0.2, ..., 1.0 at 1X to find optimal regularization.

#### Experiment 4: Fragment-Level TOO vs. Window-Level QP (Section 3.5)

- On all samples with matched WGBS:
  - (a) Window-level QP from WGBS ground truth (current gold standard)
  - (b) Window-level QP from FinaleMe-predicted methylation (current approach)
  - (c) Fragment-level EM from FinaleMe-predicted methylation (new approach)
- Compare TOO proportions: Pearson r between (a) and (b), between (a) and (c).
- At 1X: (c) should outperform (b) because it preserves co-methylation structure.

#### Experiment 5: Unknown Cell Type Estimation (Section 3.5.5)

- On high-coverage healthy samples where TOO is well-characterized: verify that π_unknown is small and consistent with expected missing cell types.
- On cancer samples: verify that adding tumor-relevant cell types to the panel reduces π_unknown.
- Artificially remove a known major cell type (e.g., Neutrophil) from the reference panel and confirm that π_unknown increases proportionally.

#### Experiment 6: Statistical Reliability (Section 3.5.8)

- On ULP-WGS cancer samples: report p-values and q-values for tumor-related cell types.
- Correlate significance (q < 0.05) with ichorCNA tumor fraction.
- Expected: samples with high TF show significant tumor cell type contributions; healthy samples show non-significant tumor contributions.
- Bootstrap CI width as a function of coverage (from downsampling experiments).
- Compare LRT p-values vs. permutation p-values for rare cell types to validate asymptotic approximation.

---

## 7. Deliverables

| Deliverable | Format | Description |
|-------------|--------|-------------|
| Modified CpgFeatureMatrixBuilder.java | Java | 5′ end motif extraction from .2bit (1st positional arg), coverage dropout, model save/load |
| Modified FinaleMe.java | Java | Constrained BW adaptation, model save/load integration |
| generate_reference_panel.py | Python | Build per-CpG reference matrix from pat.gz + marker regions (from generate_cgi_shore_markers.py) |
| fragment_too_deconvolution.py | Python | Fragment-level EM with unknown component, bootstrap CI, LRT/permutation p-values, FDR q-values |
| Reference model | Binary | Trained on HD_45, 3 features (length, dist, motif), ready for distribution |
| Reference panel(s) | TSV | Generated from Moss 2018 / New Oligo / Old Oligo pat.gz + marker regions |
| Validation figures | PDF/PNG | Experiments 1-6 results |
| Updated README.md | Markdown | Usage instructions for new parameters and workflow |

---

## 8. Risk Assessment

| Risk | Impact | Mitigation |
|------|--------|------------|
| End motif adds minimal predictive value at the HMM level | Low impact on Stage 1 overall; TOO deconvolution is the main contribution | Proceed with motif feature but make it optional; the coverage-free + adaptation improvements are independent |
| Adaptation destabilizes at very low coverage (<0.5X) | Some samples may get worse predictions than reference model | adaptMinFragments threshold + fallback to reference model; log diagnostics |
| Fragment-level EM does not outperform window-level QP | Undermines the framing but is still a valid negative result | Report both methods; the statistical testing (p/q values) is still a contribution even if proportions are similar |
| π_unknown absorbs signal from real cell types due to noisy methylation predictions | Overestimates unknown fraction, underestimates known cell types | Validate on high-coverage data where predictions are accurate; consider making unknown component optional (--estimate-unknown flag) |
| LRT asymptotic approximation is poor for rare cell types | P-values may be miscalibrated for π_t near 0 | Validate LRT p-values against permutation p-values on a subset of samples; use permutation as fallback |
| Reference panel from pat.gz has different CpG coverage than reference from bigWig | Inconsistent marker sets across input formats | Support both pat.gz and bigWig inputs in generate_reference_panel.py; validate cross-format concordance |

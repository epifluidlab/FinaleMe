# FinaleMe TOO (Tissues-of-Origin) Module — Architecture Design v2

## Package: `finaleme.too`

**Version**: 0.2-DRAFT  
**Authors**: Yaping Liu Lab, Northwestern University Feinberg School of Medicine  
**Date**: April 2026  
**Status**: Design Document v2  
**Language**: Python (with optional Cython/numba acceleration for fragment-level likelihood)

---

## 1. Overview

The TOO module performs cell-free DNA (cfDNA) tissue-of-origin deconvolution from methylation data. It supersedes the binary-threshold approach in `BetaValueDeconvolution` with a principled statistical framework featuring:

- Beta-binomial observation model replacing binary 0.1 thresholding
- Coverage-aware weighting throughout the pipeline
- Dual-mode support: WGBS/EM-seq (direct) and FinaleMe (indirect) measurement
- Three-tier coverage strategy (high / low / ultra-low)
- Region-specific calibration for FinaleMe predictions with quantitative goodness metrics
- Bootstrapped and Bayesian uncertainty quantification with per-cell-type reliability p-values
- Missing/unknown cell-type estimation (always on)
- Covariate correction (technical + biological)
- Cross-sample statistical testing: compositional regression (ILR) or Bayesian posterior probability, with FDR control
- Multi-group and pairwise contrast comparisons

---

## 2. Pipeline Overview

```
┌─────────────────────────────────────────────────────────────────┐
│                        INPUT LAYER                              │
│  Sample sheet + Methylation data + Reference panel              │
│  Mode: {WGBS, FINALEME}   Coverage tier: {HIGH, LOW, ULTRALOW} │
└────────────────────────┬────────────────────────────────────────┘
                         │
                         ▼
┌─────────────────────────────────────────────────────────────────┐
│               STAGE 1: PREPROCESSING                            │
│  1a. Input format detection & parsing                           │
│  1b. Coverage assessment → tier assignment                      │
│  1c. Mode-specific calibration (FinaleMe only) + goodness eval  │
│  1d. Technical batch correction (pre-deconvolution)             │
│  1e. Marker filtering & selection                               │
│  1f. Effective coverage computation                             │
└────────────────────────┬────────────────────────────────────────┘
                         │
                         ▼
┌─────────────────────────────────────────────────────────────────┐
│            STAGE 2: OBSERVATION MODEL                           │
│  Beta-binomial likelihood per marker region                     │
│  Inputs: methylated count k_i, total count n_i, dispersion φ_i │
│  Output: posterior distribution on true methylation per marker  │
└────────────────────────┬────────────────────────────────────────┘
                         │
                         ▼
┌─────────────────────────────────────────────────────────────────┐
│             STAGE 3: DECONVOLUTION                              │
│  Constrained optimization (MLE) or Bayesian (Dirichlet-MCMC)   │
│  Coverage-weighted objective function                           │
│  Cell-type-balanced marker weighting                            │
│  Unknown/missing cell-type component (always on)                │
│  Per-cell-type reliability p-values                             │
└────────────────────────┬────────────────────────────────────────┘
                         │
                         ▼
┌─────────────────────────────────────────────────────────────────┐
│        STAGE 4: UNCERTAINTY QUANTIFICATION                      │
│  Bootstrap CIs (fast) or Bayesian posterior (full)              │
│  Per-cell-type proportion credible intervals                    │
│  Reliability p-value per cell type per sample                   │
└────────────────────────┬────────────────────────────────────────┘
                         │
                         ▼
┌─────────────────────────────────────────────────────────────────┐
│        STAGE 5: POST-DECONVOLUTION ANALYSIS                     │
│  5a. Biological covariate adjustment                            │
│  5b. Cross-sample statistical testing (ILR regression or        │
│      Bayesian posterior probability) with FDR                   │
│  5c. Multi-group and pairwise contrast comparisons              │
│  5d. Missing cell-type residual analysis + NMF discovery        │
└────────────────────────┬────────────────────────────────────────┘
                         │
                         ▼
┌─────────────────────────────────────────────────────────────────┐
│                     OUTPUT LAYER                                │
│  Per-sample proportions + CIs + reliability p-values            │
│  Cohort summary with p-values                                   │
│  Multi-group + pairwise comparison results                      │
│  Unknown fraction estimates                                     │
│  Calibration quality report (FinaleMe mode)                     │
│  QC flags                                                       │
└─────────────────────────────────────────────────────────────────┘
```

---

## 3. Input Specification

### 3.1 Sample Sheet (TSV)

Required columns:
| Column | Type | Description |
|--------|------|-------------|
| `sample_id` | string | Unique sample identifier |
| `methylation_file` | path | Path to methylation data (see §3.2 for formats) |
| `mode` | enum | `WGBS` or `FINALEME` |
| `input_format` | enum | See §3.2. Auto-detected if omitted |
| `group` | string | Biological group label (for cohort imputation and comparisons) |

Optional columns — technical (pre-deconvolution batch correction):
| Column | Type | Stage |
|--------|------|-------|
| `extraction_batch` | string | Batch correction |
| `library_date` | date | Batch correction |
| `sequencing_date` | date | Batch correction |
| `institute` | string | Batch correction |
| `plasma_draw_date` | date | Batch correction |

Optional columns — biological/clinical (post-deconvolution covariate adjustment):
| Column | Type | Stage |
|--------|------|-------|
| `age` | numeric | Covariate adjustment |
| `sex` | enum | Covariate adjustment |
| `race` | string | Covariate adjustment |
| `ethnicity` | string | Covariate adjustment |
| `treatment` | string | User-configurable |
| `treatment_efficacy` | string | User-configurable |
| `mutation_status` | string | User-configurable |
| `disease_duration` | numeric | Covariate adjustment |
| `follow_up_time` | numeric | Covariate adjustment |
| `smoking_status` | enum | Covariate adjustment |
| `alcohol_usage` | enum | Covariate adjustment |
| `bmi` | numeric | Covariate adjustment |

Optional columns — QC/technical metadata:
| Column | Type | Stage |
|--------|------|-------|
| `cfdna_concentration_ng_ml` | numeric | Dispersion model, QC |
| `draw_time` | time | Covariate adjustment |
| `hemolysis_flag` | boolean | QC flag |
| `fasting_status` | boolean | Covariate adjustment |

### 3.2 Methylation Data Formats

**FinaleMe mode (`mode=FINALEME`):**

| Format ID | Description | Columns |
|-----------|-------------|---------|
| `finaleme_bed` | FinaleMe prediction.bed.gz output (default) | `chrom, start, end, methylated_count, total_count, predicted_beta` |

**WGBS mode (`mode=WGBS`):**

| Format ID | Description | Columns |
|-----------|-------------|---------|
| `bissnp_6plus2` | Bis-SNP 6+2 BED format (default for WGBS) | Standard 6-column BED + `methylation_pct` (col 7) + `total_count` (col 8). Methylated count derived: `round(methylation_pct/100 * total_count)` |
| `wgbstools_beta` | wgbstools .beta file | Binary format, parsed via wgbstools API |
| `custom_bed` | User-defined column mapping | Any BED-like file. User specifies `--meth-col` and `--total-col` (1-indexed column numbers for methylated count and total count) |

**Auto-detection logic:**
```
if file ends with .beta → wgbstools_beta
elif mode == FINALEME → finaleme_bed
elif file has 8+ columns and col7 looks like percentage (0-100 range) → bissnp_6plus2
else → require explicit --input-format or --meth-col/--total-col
```

### 3.3 Reference Panel Formats

Two supported input modes:

**Mode A: Matrix file (simple)**
```
chrom  start  end  celltype_1  celltype_2  ...  celltype_K
chr1   1000   1500  0.85        0.12        ...  0.45
```

**Mode B: Beta file list + group file (compatible with BetaValueDeconvolution)**

```bash
--ref-betas beta_list.txt \        # one .beta file path per line
--ref-groups groups_pat_ref.csv \   # CSV: marker_name, cell_type_label
--cpg-index CpG.bed.gz             # CpG coordinate index
```

This matches the existing `-refBetas`, `-refGroups`, `-cpgIndex` interface in `BetaValueDeconvolution`. The TOO module parses these into the internal matrix representation.

Optional companion: reference coverage matrix (same dimensions as reference matrix) for reference uncertainty weighting.

### 3.4 Marker Regions

Defines which genomic regions are used for deconvolution. Three supported formats:

**Format A: BED file (simple)**
```
chrom  start  end  [marker_name]
chr1   1000   1500  marker_001
chr1   5000   5800  marker_002
```

Standard 3- or 4-column BED. Methylation data and reference panel values are aggregated within these regions.

**Format B: UXM atlas file (compatible with UXM_deconv)**

The atlas format used by UXM_deconv / wgbstools, where each row defines a marker region with its CpG coordinates and cell-type-specific U/M (unmethylated/methylated) patterns. **Note**: the U/M ratios in the atlas are not methylation levels — they are binary-like markers used by UXM's unique-marker algorithm. For TOO deconvolution, only the **marker region coordinates** are extracted from the atlas file; a separate reference panel with actual methylation levels must still be provided via `--reference-panel`, `--ref-betas`, or another source.

```bash
--marker-regions markers.atlas.gz \     # UXM atlas: coordinates only
--reference-panel reference_betas.tsv   # separate methylation reference required
```

**Format C: BetaValueDeconvolution-compatible (via --marker-regions flag)**

Compatible with `-markerRegions` in `BetaValueDeconvolution`. BED format with marker region coordinates, used together with `--ref-betas`/`--ref-groups`/`--cpg-index`.

```bash
--marker-regions marker_regions.bed \
--ref-betas beta_list.txt \
--ref-groups groups_pat_ref.csv \
--cpg-index CpG.bed.gz
```

**Resolution logic:**
```
if --marker-regions is .atlas.gz → parse as UXM atlas; extract marker coordinates only
                                    (separate reference panel still required)
elif --marker-regions is .bed → use as region coordinates; require separate reference panel
elif --reference-panel provided (matrix with chrom/start/end) → use reference coordinates as markers
else → error: marker regions must be specified
```

CLI flags:
```bash
--marker-regions <path>         # BED or UXM atlas file
--marker-format auto            # auto, bed, uxm_atlas (auto-detect by extension)
```

### 3.5 Region Annotation (shipped with TOO)

Pre-computed per-marker-region annotations used for FinaleMe calibration (CpG density binning):
```
chrom  start  end  cpg_count  cpg_density  region_class  mean_cpg_spacing
```

Where `region_class` ∈ {CGI, shore, shelf, open_sea}.

If not provided, the TOO module computes these annotations automatically from the marker regions and a CpG index file (if available) or from the reference genome.

---

## 4. Coverage Tier System

### 4.1 Tier Definitions

| Tier | Global Coverage | Behavior |
|------|----------------|----------|
| `HIGH` | >10× | All markers usable; per-sample dispersion; no imputation |
| `LOW` | 0.5–10× | Reduced marker set; optional cohort imputation; smoothed coverage |
| `ULTRALOW` | <0.5× | Aggressive marker pruning; cohort imputation required; fragment-level likelihood |

### 4.2 Tier Assignment

Per-sample based on genome-wide mean coverage, modulated per-marker by **effective coverage**:

```
effective_coverage(marker_i, sample_s) = observed_fragments(i, s) / expected_fragments(i, s)
```

A marker with effective coverage below a tier-specific threshold is treated as belonging to the next lower tier for that sample.

### 4.3 Per-Tier Marker Filtering

| Tier | Min reads per marker | Marker selection strategy |
|------|---------------------|--------------------------|
| `HIGH` | ≥3 | All markers in reference panel |
| `LOW` | ≥2 | Top-K most informative per cell type; prefer larger regions |
| `ULTRALOW` | ≥1 (or fragment-level) | Only largest/most discriminative markers; fragment-level fallback |

### 4.4 Per-Tier Imputation Strategy

| Tier | Imputation | Source |
|------|-----------|--------|
| `HIGH` | None | — |
| `LOW` | Marker-level, for markers below local coverage threshold | Same-group cohort (from sample sheet) |
| `ULTRALOW` | Marker-level, required for most markers | Same-group cohort; external reference if no cohort |

**Critical constraint**: never impute across comparison groups.

---

## 5. Stage 1: Preprocessing

### 5.1 Input Parsing

```python
class MethylationLoader:
    """Unified loader supporting all input formats."""

    @staticmethod
    def load(filepath: str, mode: MeasurementMode,
             input_format: Optional[str] = None,
             meth_col: Optional[int] = None,
             total_col: Optional[int] = None) -> MarkerObservations:
        """
        Returns MarkerObservations with fields:
          - chrom, start, end (genomic coordinates)
          - k (methylated count per marker)
          - n (total count per marker)
          - predicted_beta (FinaleMe mode only, raw prediction before calibration)
        """
        if input_format is None:
            input_format = auto_detect_format(filepath, mode)

        if input_format == 'finaleme_bed':
            return _parse_finaleme_bed(filepath)
        elif input_format == 'bissnp_6plus2':
            return _parse_bissnp(filepath)
        elif input_format == 'wgbstools_beta':
            return _parse_wgbstools(filepath)
        elif input_format == 'custom_bed':
            return _parse_custom(filepath, meth_col, total_col)


class ReferencePanelLoader:
    """Supports matrix file or beta-list + groups."""

    @staticmethod
    def load_matrix(filepath: str) -> ReferencePanel: ...

    @staticmethod
    def load_beta_list(beta_list: str, groups_file: str,
                       cpg_index: str) -> ReferencePanel:
        """Compatible with BetaValueDeconvolution -refBetas/-refGroups/-cpgIndex."""
        ...


class MarkerRegionsLoader:
    """Load marker regions from BED or UXM atlas format."""

    @staticmethod
    def load(filepath: str, marker_format: str = 'auto') -> MarkerRegions:
        """
        Returns MarkerRegions with fields:
          - chrom, start, end (genomic coordinates per marker)
          - marker_name (optional)
        Note: UXM atlas files contain U/M ratio patterns, NOT methylation
        levels. Only the marker coordinates are extracted. A separate
        reference panel with actual methylation beta values is always
        required for TOO deconvolution.
        """
        if marker_format == 'auto':
            marker_format = _detect_marker_format(filepath)

        if marker_format == 'uxm_atlas':
            return _parse_uxm_atlas_coordinates(filepath)  # coordinates only
        elif marker_format == 'bed':
            return _parse_marker_bed(filepath)

    @staticmethod
    def _detect_marker_format(filepath: str) -> str:
        if filepath.endswith('.atlas.gz') or filepath.endswith('.atlas'):
            return 'uxm_atlas'
        return 'bed'
```

### 5.2 Coverage Assessment

```python
def assign_coverage_tier(sample: MarkerObservations) -> CoverageTier:
    mean_cov = sample.total_fragments * FRAGMENT_LENGTH / GENOME_SIZE
    if mean_cov > 10:
        return CoverageTier.HIGH
    elif mean_cov >= 0.5:
        return CoverageTier.LOW
    else:
        return CoverageTier.ULTRALOW
```

### 5.3 FinaleMe-Specific Calibration + Goodness Evaluation

#### 5.3.1 Calibration Model

Marker regions binned by CpG density into B bins (default B=8, tunable). For each bin b:

```
logit(μ_calibrated) = a_b · logit(μ_FinaleMe) + c_b
φ_b = exp(d_b)
```

Parameters `{a_b, c_b, d_b}` pre-trained from matched WGBS/FinaleMe cfDNA samples.

#### 5.3.2 Calibration Goodness Metrics

The calibration quality is evaluated at two levels — during training (to select optimal B and validate the model) and at inference time (to flag samples where calibration may be unreliable).

**Training-time metrics (reported in calibration config):**

| Metric | Definition | Good | Concerning |
|--------|-----------|------|------------|
| Per-bin R² | Variance explained by calibration regression within each bin | >0.7 | <0.5 |
| Per-bin MAE | Mean absolute error of calibrated vs. true (WGBS) methylation | <0.10 | >0.20 |
| Per-bin calibration slope | Slope a_b of logit-logit regression | 0.8–1.2 | <0.5 or >1.5 |
| Cross-validated RMSE | Leave-one-sample-out RMSE across all bins | <0.12 | >0.20 |
| Dispersion consistency | Coefficient of variation of φ_b across CV folds | <0.3 | >0.5 |
| Hosmer-Lemeshow test p | Goodness-of-fit test per bin (group observed vs. expected) | >0.05 | <0.01 |

**Inference-time metrics (per sample, reported in QC output):**

| Metric | Definition | Purpose |
|--------|-----------|---------|
| Prediction range coverage | Fraction of markers whose FinaleMe prediction falls within the training range of its bin | Detects extrapolation |
| Residual distribution test | KS test comparing this sample's calibration residuals to the training residual distribution | Detects distributional shift |
| Bin coverage balance | Fraction of bins with ≥10 markers available | Detects uneven marker representation |

**Report format (`calibration_report.json`):**
```json
{
  "calibration_version": "1.0",
  "n_training_samples": 50,
  "n_bins": 8,
  "overall_metrics": {
    "cross_validated_rmse": 0.098,
    "mean_r_squared": 0.82,
    "hosmer_lemeshow_p_min": 0.12
  },
  "per_bin_metrics": [
    {
      "bin_id": 0,
      "cpg_density_range": [0.0, 0.01],
      "n_markers": 5200,
      "r_squared": 0.58,
      "mae": 0.15,
      "slope": 0.72,
      "log_dispersion": 1.2,
      "hosmer_lemeshow_p": 0.08
    },
    ...
  ],
  "per_sample_qc": [
    {
      "sample_id": "S001",
      "prediction_range_coverage": 0.95,
      "residual_ks_p": 0.42,
      "bin_coverage_balance": 1.0,
      "calibration_flag": "PASS"
    },
    ...
  ]
}
```

#### 5.3.3 Tuning the Number of Bins B

The optimal B is selected during calibration training via cross-validation:

```python
def tune_n_bins(matched_wgbs, matched_finaleme, B_candidates=[4, 6, 8, 10, 12, 16]):
    """Select B that minimizes cross-validated RMSE."""
    results = []
    for B in B_candidates:
        cv_rmse = cross_validate_calibration(matched_wgbs, matched_finaleme, n_bins=B)
        aic = compute_aic(matched_wgbs, matched_finaleme, n_bins=B)
        results.append({'B': B, 'cv_rmse': cv_rmse, 'aic': aic,
                        'n_params': B * 3})  # 3 params per bin: a, c, d
    # Select B minimizing CV-RMSE; use AIC to penalize overfitting
    # Report all candidates so user can make informed choice
    return results
```

Considerations: too few bins → underfitting (open-sea and CGI lumped together); too many bins → overfitting (insufficient markers per bin for stable estimation). Default B=8 is a starting point; the training pipeline reports all candidates.

### 5.4 Technical Batch Correction

Applied **before deconvolution** at the marker methylation level. ComBat-style location/scale adjustment. Only applied if ≥2 levels per batch variable and ≥5 samples per level.

### 5.5 Marker Selection

Cell-type-balanced selection: top-N most discriminative markers per cell type. Optional `--strict-regions CGI+shore` flag restricts to CGI and shore markers only.

### 5.6 Effective Coverage Computation

Per-marker, per-sample. High tier uses raw counts; low tier uses smoothed; ultra-low defers to fragment-level.

---

## 6. Stage 2: Observation Model

### 6.1 Beta-Binomial Likelihood

For each marker region i:
```
k_i ~ BetaBinomial(n_i, α_i = μ_i·φ_i, β_i = (1-μ_i)·φ_i)
```

Where μ_i = Σ_j w_j · r_ij + w_unknown · 0.5.

### 6.2 Mode-Dependent Dispersion

| Mode | Tier | φ source | Typical range |
|------|------|----------|---------------|
| WGBS | HIGH | Per-sample from data | 50–200 |
| WGBS | LOW | Per-cohort | 20–100 |
| WGBS | ULTRALOW | Fixed conservative | ~10 |
| FinaleMe | HIGH | Calibration config, per density bin | 5–30 |
| FinaleMe | LOW | Calibration config, reduced | 3–15 |
| FinaleMe | ULTRALOW | Fixed low | 2–5 |

---

## 7. Stage 3: Deconvolution

### 7.1 MLE Deconvolution

```
maximize   Σ_i ω_i · ℓ_i(w)
subject to w_j ≥ 0, Σ w_j + w_unknown = 1
```

Unknown component (flat profile at 0.5) is **always included** — its estimated proportion w_unknown represents the fraction of signal not explained by the reference panel.

### 7.2 Bayesian Deconvolution (optional `--bayesian`)

Dirichlet prior + beta-binomial likelihood, sampled via MCMC.

### 7.3 Per-Cell-Type Reliability P-Values

Two complementary p-values are computed for each cell type j in each sample:

**P_goodness (marker-level goodness-of-fit):**

For each cell type j, select its top-discriminative markers. Compute a chi-squared-like statistic comparing observed methylation at these markers to the prediction under the estimated mixture:

```python
def compute_p_goodness(w, reference, observations, dispersions, cell_type_j):
    """How well do this cell type's discriminative markers support the estimate?"""
    marker_indices = get_discriminative_markers(reference, cell_type_j, top_n=50)
    chi2_stat = 0.0
    for i in marker_indices:
        k_i, n_i = observations[i]
        mu_i = reference[i, :] @ w
        expected_k = mu_i * n_i
        variance = n_i * mu_i * (1 - mu_i) * (n_i + dispersions[i]) / (1 + dispersions[i])
        chi2_stat += (k_i - expected_k)**2 / max(variance, 1e-10)
    df = len(marker_indices) - 1
    p_goodness = 1 - chi2.cdf(chi2_stat, df)
    return p_goodness  # high p → good fit → trustworthy estimate
```

Interpretation: p_goodness > 0.05 → the estimated proportion is consistent with the observed data at this cell type's markers. p_goodness < 0.01 → the model fit is poor for this cell type, suggesting the estimate may be unreliable.

**P_detection (proportion stability):**

From the bootstrap distribution of w_j:

```python
def compute_p_detection(bootstrap_proportions_j, noise_floor=0.001):
    """Is the cell type reliably detected above noise?"""
    n_below_noise = np.sum(bootstrap_proportions_j < noise_floor)
    p_detection = 1 - (n_below_noise / len(bootstrap_proportions_j))
    return p_detection  # high p → reliably detected
```

Interpretation: p_detection > 0.95 → the cell type is reliably detected in this sample. p_detection < 0.10 → the estimated proportion is unstable and may be noise.

**Combined reliability assessment:**

| p_goodness | p_detection | Interpretation |
|-----------|-------------|----------------|
| >0.05 | >0.95 | High confidence: well-supported, stable estimate |
| >0.05 | <0.50 | Low confidence: model fits but proportion near noise floor |
| <0.01 | >0.95 | Caution: stably detected but model fit is poor (possible reference mismatch) |
| <0.01 | <0.50 | Unreliable: poor fit and unstable |

---

## 8. Stage 4: Uncertainty Quantification

### 8.1 Bootstrap CIs

Resample markers with replacement, re-run deconvolution N times (default 1000), compute percentile CIs. Reliability p-values computed from bootstrap distribution (see §7.3).

### 8.2 Bayesian Posterior

Dirichlet-BetaBinomial model via MCMC. CIs from posterior quantiles. Reliability from posterior mass above noise floor.

### 8.3 Covariate-Adjusted CIs

Covariate regression included within each bootstrap iteration to propagate adjustment uncertainty.

---

## 9. Stage 5: Post-Deconvolution Analysis

### 9.1 Biological Covariate Adjustment

Applied **after deconvolution** using ILR transform. Covariates `treatment`, `treatment_efficacy`, `mutation_status` are user-configurable (include/exclude). All other covariates in the sample sheet are included by default if present.

### 9.2 Cross-Sample Statistical Testing

**Default: Compositional regression in ILR space** (most principled for proportion data):

```python
def compositional_regression_test(proportions, group_labels, covariates,
                                   cell_type_names, contrasts):
    """
    Fit multivariate regression in ILR space, test group coefficients.

    proportions: (S x K) matrix of estimated proportions
    group_labels: (S,) group assignments
    covariates: (S x C) covariate matrix
    contrasts: list of (group_A, group_B) pairs to test
    """
    # Transform to ILR space
    ilr_coords = ilr_transform(proportions)  # (S x K-1)

    results = []
    for cell_type_j in range(K):
        # Univariate regression per cell type (in ILR balances involving j)
        # Design matrix: group indicators + covariates
        X = build_design_matrix(group_labels, covariates, contrasts)
        y = ilr_coords[:, j]

        # Fit OLS with bootstrap standard errors
        # (bootstrap SEs account for deconvolution uncertainty)
        model = OLS(y, X).fit()

        for contrast_name, contrast_vec in contrasts.items():
            effect = contrast_vec @ model.params
            se = np.sqrt(contrast_vec @ model.cov_params() @ contrast_vec)
            z = effect / se
            p = 2 * (1 - norm.cdf(abs(z)))
            results.append({
                'cell_type': cell_type_names[j],
                'contrast': contrast_name,
                'effect_size': effect,
                'se': se,
                'z_statistic': z,
                'p_value': p
            })

    # FDR correction across all cell_type × contrast combinations
    all_pvals = [r['p_value'] for r in results]
    _, adjusted = fdrcorrection(all_pvals, alpha=fdr_alpha)
    for r, adj_p in zip(results, adjusted):
        r['adjusted_p_value'] = adj_p

    return results
```

**Bayesian mode: Posterior probability of difference**

```python
def bayesian_group_comparison(posterior_samples_by_group, contrasts):
    """
    posterior_samples_by_group: {group: {sample_id: (n_draws x K) array}}
    """
    results = []
    for cell_type_j in range(K):
        for group_A, group_B in contrasts:
            # Pool posterior draws within each group
            draws_A = aggregate_group_posterior(posterior_samples_by_group[group_A], j)
            draws_B = aggregate_group_posterior(posterior_samples_by_group[group_B], j)
            # P(group_A > group_B) from paired draws
            prob = np.mean(draws_A > draws_B)
            # Convert to two-sided pseudo p-value for FDR compatibility
            p_pseudo = 2 * min(prob, 1 - prob)
            results.append({
                'cell_type': cell_type_names[j],
                'contrast': f'{group_A}_vs_{group_B}',
                'prob_A_gt_B': prob,
                'pseudo_p_value': p_pseudo,
                'mean_A': np.mean(draws_A),
                'mean_B': np.mean(draws_B)
            })
    # FDR on pseudo p-values
    ...
    return results
```

**Why not Wilcoxon (legacy option only)?**

Wilcoxon rank-sum ignores the per-sample uncertainty estimates and treats each proportion as a fixed number. Since the TOO module produces full bootstrap or posterior distributions, this discards valuable information. Compositional regression in ILR space is more powerful because it accounts for the simplex constraint, incorporates covariates directly, and can be combined with bootstrap SEs that propagate deconvolution uncertainty. Wilcoxon is retained as `--test-method wilcoxon` for backward compatibility and for quick exploratory analysis.

### 9.3 Multi-Group and Pairwise Comparisons

The `--group-comparison` flag supports multiple syntaxes:

```bash
# All pairwise comparisons across all groups
--group-comparison all

# Specific contrasts
--group-comparison RRMS:SPMS,RRMS:PPMS,SPMS:PPMS

# One-vs-rest
--group-comparison SPMS:rest

# Multi-group omnibus test + all pairwise post-hoc
--group-comparison omnibus+pairwise
```

**Omnibus test**: For multi-group comparison (>2 groups), first run a global test (MANOVA in ILR space or Kruskal-Wallis per cell type) to test whether any group differs. If significant, proceed to pairwise contrasts with FDR correction across all contrasts × cell types.

**Output structure**: The group_comparison output contains separate rows for each contrast × cell type combination, enabling both per-contrast and per-cell-type views.

### 9.4 Missing Cell-Type Analysis

Unknown component proportion w_unknown is always estimated (§7.1). Additional residual analysis and NMF discovery for cohort-level identification of unknown components.

---

## 10. Output Specification

### 10.1 Per-Sample Output (`{sample_id}.too.tsv`)

```
cell_type    proportion    ci_lower    ci_upper    p_goodness    p_detection    reliability    n_markers    mean_dispersion
Neutrophil   0.452         0.410       0.491       0.72          0.99           HIGH           1200         45.2
Monocyte     0.082         0.065       0.101       0.45          0.98           HIGH           980          38.7
Oligodendro  0.008         0.002       0.018       0.31          0.62           LOW            340          12.5
Neuron       0.002         0.000       0.008       0.15          0.28           UNRELIABLE     280          10.1
...
Unknown      0.035         0.012       0.061       —             0.88           MODERATE       —            —
```

Reliability categories derived from p_goodness × p_detection (see §7.3 table):
- `HIGH`: p_goodness > 0.05 AND p_detection > 0.95
- `MODERATE`: p_goodness > 0.05 AND p_detection ∈ [0.50, 0.95]
- `LOW`: p_goodness > 0.05 AND p_detection < 0.50, OR p_goodness < 0.05 AND p_detection > 0.50
- `UNRELIABLE`: p_goodness < 0.01 AND p_detection < 0.50

### 10.2 Cohort Summary (`cohort_proportions.tsv`)

```
sample_id  group  coverage_tier  Neutrophil  Neutrophil_ci_lo  Neutrophil_ci_hi  Neutrophil_p_goodness  Neutrophil_p_detection  Neutrophil_reliability  ...  Unknown  Unknown_ci_lo  Unknown_ci_hi  qc_flags
S001       RRMS   HIGH           0.452       0.410             0.491             0.72                   0.99                    HIGH                    ...  0.035    0.012          0.061          NONE
S002       SPMS   LOW            0.389       0.340             0.441             0.55                   0.97                    HIGH                    ...  0.058    0.025          0.098          NONE
```

### 10.3 Group Comparison Output (`group_comparison.tsv`)

Supports both omnibus and pairwise contrasts:

```
test_type   contrast        cell_type     mean_A  mean_B  effect_size  se      statistic   p_value    adjusted_p   significant
omnibus     all_groups      Neutrophil    —       —       —            —       12.45       0.0020     0.0140       TRUE
omnibus     all_groups      Oligodendro   —       —       —            —       18.72       0.0001     0.0014       TRUE
pairwise    RRMS_vs_SPMS    Neutrophil    0.452   0.389   0.063        0.018   3.50        0.0005     0.0035       TRUE
pairwise    RRMS_vs_SPMS    Oligodendro   0.012   0.028   -0.016       0.005   -3.20       0.0014     0.0070       TRUE
pairwise    RRMS_vs_PPMS    Oligodendro   0.012   0.032   -0.020       0.006   -3.33       0.0009     0.0054       TRUE
pairwise    SPMS_vs_PPMS    Oligodendro   0.028   0.032   -0.004       0.007   -0.57       0.5680     0.6820       FALSE
bayesian    RRMS_vs_SPMS    Neutrophil    0.452   0.389   —            —       —           0.0010     0.0070       TRUE
```

Columns adapt to test method: frequentist tests include effect_size/se/statistic; Bayesian includes prob_A_gt_B.

### 10.4 Calibration Quality Report (`calibration_report.json`)

See §5.3.2 for full schema.

### 10.5 Residual Analysis Output (`residual_analysis.tsv`)

```
sample_id  unexplained_fraction  mean_residual  residual_sd  qc_flag
```

### 10.6 QC Summary (`qc_summary.tsv`)

```
sample_id  coverage_tier  mean_coverage  n_markers_used  pct_imputed  wbc_fraction  unknown_fraction  calibration_flag  hemolysis  overall_qc
S001       HIGH           25.3           2400             0.0          0.78          0.035             PASS              FALSE      PASS
S002       LOW            3.1            1800             12.5         0.72          0.058             PASS              FALSE      PASS
S003       ULTRALOW       0.3            450              68.2         0.81          0.142             WARN              FALSE      WARN
```

---

## 11. Configuration

### 11.1 Command-Line Interface

```bash
finaleme-too run \
  --sample-sheet samples.tsv \
  --output-dir results/ \
  # --- Reference panel (choose one) ---
  --reference-panel reference.tsv \
  # --- OR BetaValueDeconvolution-compatible format ---
  --ref-betas beta_list.txt \
  --ref-groups groups_pat_ref.csv \
  --cpg-index CpG.bed.gz \
  # --- OR UXM atlas (marker coordinates only; reference panel still needed) ---
  --marker-regions markers.atlas.gz \
  # --- Marker regions (required unless using atlas or matrix with coordinates) ---
  --marker-regions marker_regions.bed \
  --marker-format auto \
  # --- Mode & format ---
  --mode FINALEME \
  --input-format finaleme_bed \
  --calibration calibration_params.json \
  # --- OR for WGBS ---
  --mode WGBS \
  --input-format bissnp_6plus2 \
  # --- OR custom column mapping ---
  --input-format custom_bed --meth-col 5 --total-col 6 \
  # --- Marker selection ---
  --region-annotation regions.tsv \
  --strict-regions CGI+shore \
  --n-markers-per-type 500 \
  # --- Model ---
  --n-bootstrap 1000 \
  --bayesian \
  --bayesian-n-samples 5000 \
  # --- Covariates ---
  --batch-correct extraction_batch,library_date \
  --adjust-covariates age,sex,race,disease_duration,smoking_status \
  --configurable-covariates treatment,treatment_efficacy \
  # --- Testing ---
  --test-method ilr_regression \
  --group-comparison omnibus+pairwise \
  --fdr-method BH \
  --fdr-alpha 0.05 \
  # --- Performance ---
  --threads 8 \
  --min-coverage 3 \
  --coverage-cap 50

# Calibration training subcommand
finaleme-too train-calibration \
  --matched-wgbs wgbs_samples.tsv \
  --matched-finaleme finaleme_samples.tsv \
  --marker-regions marker_regions.bed \
  --region-annotation regions.tsv \
  --n-bins-candidates 4,6,8,10,12,16 \
  --output calibration_params.json \
  --report calibration_report.json \
  --threads 8
```

### 11.2 Configuration File (`too_config.yaml`)

```yaml
model:
  observation: beta_binomial
  deconvolution: mle               # or "bayesian"
  unknown_component: true           # always on; flag kept for documentation
  fragment_level: auto              # auto enables for ULTRALOW tier

coverage:
  tier_high: 10
  tier_low: 0.5
  min_reads: 3
  coverage_cap: 50

calibration:
  mode: FINALEME
  calibration_file: calibration_params.json
  n_density_bins: 8                 # tunable; training selects optimal

markers:
  marker_regions: null              # BED or UXM atlas file (required unless reference has coords)
  marker_format: auto               # auto, bed, uxm_atlas
  region_annotation: null           # optional; auto-computed from marker regions if absent
  strict_regions: null              # or "CGI+shore"
  n_per_type: 500
  specificity_method: "entropy"

uncertainty:
  method: bootstrap                 # or "bayesian" or "both"
  n_bootstrap: 1000
  ci_level: 0.95
  noise_floor: 0.001               # for p_detection
  bayesian_prior_alpha: 1.0
  bayesian_n_samples: 5000

batch_correction:
  technical_covariates: []
  min_levels: 2
  min_samples_per_level: 5

covariate_adjustment:
  biological_covariates: []
  user_configurable: [treatment, treatment_efficacy, mutation_status]
  transform: ILR

testing:
  method: ilr_regression            # or "bayesian_posterior", "wilcoxon" (legacy)
  group_comparison: omnibus+pairwise
  fdr_method: BH
  fdr_alpha: 0.05

input:
  format: auto                      # or finaleme_bed, bissnp_6plus2, wgbstools_beta, custom_bed
  meth_col: null                    # for custom_bed
  total_col: null

reference:
  format: matrix                    # or "beta_list"
  reference_panel: null
  ref_betas: null
  ref_groups: null
  cpg_index: null
  coverage_matrix: null             # optional reference coverage

qc:
  max_wbc_fraction: 0.95
  max_unknown_fraction: 0.30
  max_residual_variance: 0.40
```

---

## 12. Package Structure

```
finaleme-too/
├── pyproject.toml
├── README.md
├── finaleme_too/
│   ├── __init__.py
│   ├── cli.py                         # Click-based CLI entry point
│   ├── pipeline.py                    # TOOPipeline orchestrator
│   ├── config.py                      # TOOConfig dataclass + YAML loading
│   │
│   ├── io/
│   │   ├── __init__.py
│   │   ├── sample_sheet.py            # SampleSheet parser/validator
│   │   ├── methylation_loader.py      # Multi-format methylation data loader
│   │   ├── reference_panel.py         # ReferencePanel: matrix + beta-list modes
│   │   ├── marker_regions.py          # MarkerRegions: BED + UXM atlas loader
│   │   └── output_writer.py           # All output TSV/JSON writers
│   │
│   ├── preprocessing/
│   │   ├── __init__.py
│   │   ├── coverage.py                # CoverageTierAssigner, EffectiveCoverage
│   │   ├── calibration.py             # FinaleMe calibration: apply + train
│   │   ├── calibration_eval.py        # Calibration goodness metrics + bin tuning
│   │   ├── batch_correction.py        # ComBat-style technical correction
│   │   ├── marker_selection.py        # BalancedMarkerSelector
│   │   └── imputation.py             # CohortImputer
│   │
│   ├── core/
│   │   ├── __init__.py
│   │   ├── observation_model.py       # BetaBinomialModel
│   │   ├── deconvolution.py           # MLEDeconvolver, BayesianDeconvolver
│   │   ├── reliability.py             # p_goodness, p_detection computation
│   │   ├── uncertainty.py             # BootstrapCI, BayesianPosterior
│   │   └── fragment_likelihood.py     # FragmentLevelDeconvolver (ULTRALOW)
│   │
│   ├── postprocessing/
│   │   ├── __init__.py
│   │   ├── covariate_adjustment.py    # ILR-based biological adjustment
│   │   ├── statistical_testing.py     # ILR regression, Bayesian comparison
│   │   ├── group_comparison.py        # Omnibus, pairwise, one-vs-rest contrasts
│   │   ├── residual_analysis.py       # Missing cell-type estimation, NMF
│   │   └── qc.py                      # QC flag generation
│   │
│   ├── data/
│   │   ├── default_calibration.json   # Shipped default calibration params
│   │   └── region_annotations/        # Pre-computed CGI/shore/shelf annotations
│   │
│   └── utils/
│       ├── __init__.py
│       ├── transforms.py              # ILR/ALR/CLR transforms
│       ├── beta_binomial.py           # Beta-binomial distribution utilities
│       └── parallel.py               # joblib-based sample parallelism
│
├── tests/
│   ├── test_observation_model.py
│   ├── test_deconvolution.py
│   ├── test_calibration.py
│   ├── test_statistical_testing.py
│   ├── test_io_formats.py
│   └── test_pipeline_integration.py
│
└── scripts/
    ├── train_calibration.py           # Standalone calibration training script
    └── convert_reference.py           # Convert between reference panel formats
```

---

## 13. Mode Comparison: WGBS vs. FinaleMe

| Aspect | WGBS/EM-seq Mode | FinaleMe Mode |
|--------|------------------|---------------|
| Input formats | bissnp_6plus2, wgbstools_beta, custom_bed | finaleme_bed |
| Calibration | None (direct measurement) | Region-specific calibration required |
| Dispersion φ | High (50–200) | Lower (2–30) |
| Noise model | Sampling + bisulfite conversion | Sampling + prediction model error |
| GC correction | Apply before deconvolution | Not applied (in FinaleMe features) |
| Fragment length correction | Apply as covariate | Not applied (in FinaleMe features) |
| CI width | Narrower | Wider |
| Calibration QC | N/A | Per-sample calibration flag |

---

## 14. Implementation Notes

### 14.1 Why Python

- Statistical ecosystem: scipy, statsmodels, emcee, sklearn provide all required machinery
- Community expectations: bioinformatics tools of this type are predominantly Python
- Development speed: substantially faster iteration than Java for statistical modeling
- Parallelism: `joblib` for sample-level parallelism (embarrassingly parallel across samples)
- Performance-critical paths: fragment-level likelihood (ULTRALOW) can use numba JIT or Cython
- Integration: consumes FinaleMe Java output files; no need to share a runtime

### 14.2 Multi-Threading Strategy (Per-Stage)

All stages support `--threads N`. The parallelization granularity varies by stage:

| Stage | Parallel Over | Granularity | Mechanism | Notes |
|-------|--------------|-------------|-----------|-------|
| **1a. Input parsing** | Samples | 1 thread per sample | joblib | I/O-bound; parallel file reading |
| **1b. Coverage assessment** | Samples | 1 thread per sample | joblib | Fast, lightweight |
| **1c. Calibration (FinaleMe)** | Samples | 1 thread per sample | joblib | Apply pre-trained model; independent per sample |
| **1c. Calibration training** | CV folds × bin candidates | 1 thread per (fold, B) pair | joblib | One-time; embarrassingly parallel |
| **1d. Batch correction** | Markers | 1 thread per marker chunk | joblib | Per-marker regression; chunk into blocks of ~500 |
| **1e. Marker selection** | Cell types | 1 thread per cell type | joblib | Specificity computation per cell type |
| **1f. Effective coverage** | Samples | 1 thread per sample | joblib | Independent per sample |
| **2. Observation model** | Samples | 1 thread per sample | numpy vectorized + joblib | Dispersion estimation parallelized across samples |
| **3. Deconvolution (MLE)** | Samples | 1 thread per sample | joblib | Each sample's optimization is independent |
| **3. Deconvolution (Bayesian)** | Samples | 1 thread per sample | joblib wrapping emcee | emcee has internal parallelism; 1 sample per thread |
| **4. Bootstrap CIs** | Bootstrap iterations within sample, AND across samples | Nested parallelism | joblib (outer: samples) + numpy (inner: vectorized bootstrap) | Dominant cost; outer parallelism preferred when S > threads |
| **4. Reliability p-values** | Samples | 1 thread per sample | joblib | Computed alongside bootstrap |
| **5a. Covariate adjustment** | Cell types | 1 thread per cell type | joblib | ILR regression per cell type |
| **5b. Statistical testing** | Cell types × contrasts | 1 thread per (cell_type, contrast) | joblib | Independent tests |
| **5c. Residual analysis** | Samples | 1 thread per sample | joblib | Independent per sample |
| **5d. NMF discovery** | — | Internal sklearn threading | sklearn n_jobs | NMF on full residual matrix |

**Parallelism architecture:**

```python
from joblib import Parallel, delayed

class TOOPipeline:
    def __init__(self, config: TOOConfig):
        self.n_threads = config.threads

    def run(self, sample_sheet, reference):
        # Stage 1: Preprocessing (parallel over samples)
        samples = Parallel(n_jobs=self.n_threads)(
            delayed(self._preprocess_sample)(s, reference)
            for s in sample_sheet.samples
        )

        # Stage 1d: Batch correction (parallel over marker chunks)
        if self.config.batch_correction.technical_covariates:
            samples = self._batch_correct_parallel(samples)

        # Stage 2+3+4: Deconvolution + bootstrap (parallel over samples)
        results = Parallel(n_jobs=self.n_threads)(
            delayed(self._deconvolve_with_uncertainty)(s, reference)
            for s in samples
        )

        # Stage 5: Post-deconvolution (parallel over cell types × contrasts)
        results = self._postprocess_parallel(results, sample_sheet)
        return results

    def _deconvolve_with_uncertainty(self, sample, reference):
        """Per-sample: deconvolution + bootstrap + reliability p-values.
        This is the most expensive step; each call is single-threaded internally
        but multiple samples run in parallel via joblib."""
        obs_model = self._build_observation_model(sample)
        w_hat = self.deconvolver.solve(obs_model, reference, sample.observations)

        # Bootstrap (vectorized numpy operations within each iteration)
        bootstrap_results = []
        for b in range(self.config.uncertainty.n_bootstrap):
            idx = np.random.choice(M, size=M, replace=True)
            w_b = self.deconvolver.solve(obs_model, reference, sample.observations, idx)
            bootstrap_results.append(w_b)

        # Reliability p-values
        p_goodness = compute_p_goodness(w_hat, reference, sample, obs_model)
        p_detection = compute_p_detection(np.array(bootstrap_results))

        return DeconvolutionResult(...)

    def _postprocess_parallel(self, results, sample_sheet):
        """Parallel statistical testing over cell types and contrasts."""
        contrasts = self._build_contrasts(sample_sheet)
        test_results = Parallel(n_jobs=self.n_threads)(
            delayed(self._test_one_contrast)(results, cell_type, contrast)
            for cell_type in range(K)
            for contrast in contrasts
        )
        # FDR correction (sequential, fast)
        test_results = apply_fdr(test_results)
        return test_results
```

**Thread allocation guidance:**

| Scenario | Recommended --threads | Rationale |
|----------|----------------------|-----------|
| Many samples (>50), moderate markers | S (up to CPU count) | Parallelize over samples |
| Few samples (<10), large bootstrap | min(S, CPU count) | Can't exceed sample count |
| Bayesian MCMC mode | S (each sample uses emcee internally) | emcee handles walkers internally |
| Calibration training | n_folds × len(B_candidates) | Embarrassingly parallel |
| Ultra-low coverage, fragment-level | S (numba JIT handles inner loop) | Fragment likelihood is single-threaded but JIT-fast |

### 14.3 Computational Complexity

| Stage | Per-sample cost | Parallelizable? | Notes |
|-------|----------------|-----------------|-------|
| Input parsing | O(M) per file | Yes (across samples) | I/O-bound |
| Calibration (apply) | O(M) | Yes (across samples) | Fast lookup per marker |
| Batch correction | O(M·S) total | Yes (across marker chunks) | Regression per marker |
| Marker selection | O(M·K) | Yes (across cell types) | Specificity scoring |
| Observation model | O(M) per sample | Yes (across samples) | Vectorized numpy |
| Deconvolution (MLE) | O(M·K·I) per sample | Yes (across samples) | I ≈ 100 iterations |
| Bootstrap (1000 iter) | 1000 × O(M·K·I) | Yes (across samples) | **Dominant cost** |
| Bayesian MCMC | O(M·K·N_samples) | Yes (across samples) | Slower; 5000 samples |
| Reliability p-values | O(M·K) per sample | Yes (across samples) | Computed with bootstrap |
| ILR regression | O(S·K·C) total | Yes (across cell types × contrasts) | Fast |
| NMF discovery | O(S·M·n_components·I) | Internal sklearn | One-time cohort analysis |

Estimated wall time for 200 samples, 2000 markers, 20 cell types, 1000 bootstrap, 8 threads: ~5–15 minutes.

---

## 15. Implementation Priority

| Priority | Component | Rationale |
|----------|-----------|-----------|
| P0 | Multi-format I/O (§3.2, §3.3, §12) | Foundation; must support all input formats |
| P0 | Beta-binomial observation model | Core improvement |
| P0 | MLE deconvolution + unknown component | Core engine |
| P0 | WGBS and FinaleMe dual-mode | Required for both measurement types |
| P0 | Per-cell-type reliability p-values | Core quality metric |
| P1 | Bootstrap CIs | Essential for downstream statistics |
| P1 | Coverage tier system + marker filtering | Required for low-coverage support |
| P1 | FinaleMe calibration model + goodness eval | Critical for FinaleMe accuracy |
| P1 | ILR compositional regression testing | Default statistical test |
| P1 | Multi-group + pairwise comparisons | Required for cohort analysis |
| P2 | Technical batch correction | Important for multi-site cohorts |
| P2 | Biological covariate adjustment | Important for confounded cohorts |
| P2 | Cohort imputation (LOW/ULTRALOW) | Required for low-coverage |
| P2 | Calibration bin tuning pipeline | One-time training utility |
| P3 | Bayesian posterior (MCMC) | Full uncertainty; computationally expensive |
| P3 | Fragment-level likelihood (ULTRALOW) | Niche use case |
| P3 | NMF residual discovery | Exploratory |
| P3 | Bayesian group comparison | Alternative to ILR regression |

---

## 16. Dependencies

```
# Core
numpy >= 1.24
scipy >= 1.10
pandas >= 2.0
scikit-learn >= 1.2
statsmodels >= 0.14
pyyaml >= 6.0
click >= 8.0          # CLI framework
joblib >= 1.3         # parallelism
tqdm >= 4.65          # progress bars

# Optional
emcee >= 3.1          # MCMC (Bayesian mode)
numba >= 0.58         # JIT for fragment-level likelihood
arviz >= 0.15         # posterior diagnostics
```

# TOO Module — Mathematical Formulation v3

## 1. Notation

| Symbol | Description |
|--------|-------------|
| M | Number of marker regions |
| K | Number of cell types in reference panel |
| S | Number of samples |
| G | Number of groups |
| B | Number of context bins (FinaleMe mode) |
| w = (w_1, ..., w_K, w_0) | Cell-type proportions + unknown, Σw = 1 |
| w_0 | Unknown/missing cell-type fraction |
| r_{ij} | Reference methylation for marker i, cell type j |
| R | Reference matrix (M × K) |
| k_i | Observed methylated read count at marker i (WGBS mode) |
| n_i | Total read count at marker i (WGBS mode) |
| s_i | Called state at marker i: U, M, or Ambiguous (FinaleMe mode) |
| φ_i | Dispersion parameter at marker i (WGBS mode) |
| μ_i | Expected methylation at marker i under mixture (WGBS mode) |
| τ_low_b, τ_high_b | Binarization thresholds for context bin b (FinaleMe mode) |
| ε_U_b | Error rate P(true=M \| called=U) in bin b (FinaleMe mode) |
| ε_M_b | Error rate P(true=U \| called=M) in bin b (FinaleMe mode) |
| c_i | Genomic context features for marker i |
| ω_i | Marker weight in objective function |

---

## 2. Observation Model — WGBS Mode (Beta-Binomial)

### 2.1 Expected Methylation Under Mixture

With unknown component (always included):
```
μ_i = Σ_{j=1}^{K} w_j · r_{ij} + w_0 · 0.5

where w_0 = 1 - Σ_{j=1}^{K} w_j
```

### 2.2 Beta-Binomial Likelihood

```
P(k_i | n_i, μ_i, φ_i) = C(n_i, k_i) · B(k_i + α_i, n_i - k_i + β_i) / B(α_i, β_i)

where:
  α_i = μ_i · φ_i
  β_i = (1 - μ_i) · φ_i
```

Log-likelihood per marker:
```
ℓ_i(w) = log Γ(n_i + 1) - log Γ(k_i + 1) - log Γ(n_i - k_i + 1)
        + log Γ(k_i + α_i) + log Γ(n_i - k_i + β_i) - log Γ(n_i + α_i + β_i)
        - log Γ(α_i) - log Γ(β_i) + log Γ(α_i + β_i)
```

### 2.3 Dispersion

Variance of the beta-binomial:
```
Var(k_i) = n_i · μ_i · (1 - μ_i) · (n_i + φ_i) / (1 + φ_i)
```

- φ → ∞: approaches binomial variance
- φ → 0: maximum overdispersion

### 2.4 Reference Uncertainty (WGBS Mode)

With reference coverage n_ref_{ij} available:
```
φ_{effective,i,j} = φ_i · min(1, n_ref_{ij} / n_ref_cap)
```

---

## 2B. Observation Model — FinaleMe Mode (Binarization with Error Rates)

### 2B.1 State Assignment

For marker i with FinaleMe predicted beta β̂_i in context bin b:
```
s_i = U        if β̂_i < τ_low_b
s_i = M        if β̂_i > τ_high_b
s_i = Ambiguous if τ_low_b ≤ β̂_i ≤ τ_high_b   (excluded from likelihood)
```

### 2B.2 Reference Binarization

For reference entry r_{ij}:
```
r_binary_{ij} = 0 (U)   if r_{ij} < 0.2
r_binary_{ij} = 1 (M)   if r_{ij} > 0.8
r_binary_{ij} = r_{ij}  if 0.2 ≤ r_{ij} ≤ 0.8  (soft assignment)
```

### 2B.3 Expected State Probabilities Under Mixture

```
P(true_i = U | w) = Σ_{j=1}^{K} w_j · (1 - r_binary_{ij}) + w_0 · 0.5
P(true_i = M | w) = Σ_{j=1}^{K} w_j · r_binary_{ij} + w_0 · 0.5
```

### 2B.4 Observation Likelihood with Error Rates

For marker i in context bin b with called state s_i:
```
P(s_i = U | w, ε) = (1 - ε_U_b) · P(true_i = U | w) + ε_U_b · P(true_i = M | w)
P(s_i = M | w, ε) = (1 - ε_M_b) · P(true_i = M | w) + ε_M_b · P(true_i = U | w)
```

Log-likelihood per marker:
```
ℓ_i(w) = log P(s_i | w, ε)
```

### 2B.5 Design Rationale

The binarization model generalizes the proven binary threshold approach in BetaValueDeconvolution:

| Original approach | Generalized approach |
|-------------------|---------------------|
| Fixed thresholds: 0.1, 0.5 | Learned thresholds τ_low_b, τ_high_b per context bin |
| CGI+shore only | All regions where FinaleMe is discriminative (data-driven) |
| No error modeling | Explicit error rates ε_U, ε_M per context bin |
| Binary: U or M | Three states: U, M, Ambiguous (excluded) |
| No uncertainty from binarization | Error rates propagate into likelihood |

This works because FinaleMe provides reliable *binary classifications* at well-separated regions, not reliable continuous values. Continuous calibration cannot improve discrimination — it can only adjust bias. The binarization model matches what FinaleMe can actually deliver.

---

## 3. MLE Deconvolution

### 3.1 Weighted Objective

```
maximize   L(w) = Σ_{i=1}^{M} ω_i · ℓ_i(w)

subject to:
  w_j ≥ 0          for j = 0, 1, ..., K
  Σ_{j=0}^{K} w_j = 1
```

Where ℓ_i(w) is mode-dependent:
- WGBS: beta-binomial log-likelihood (§2.2)
- FinaleMe: binomial log-likelihood with error rates (§2B.4)

### 3.2 Marker Weights

```
ω_i = coverage_weight_i · balance_weight_i
```

WGBS: `coverage_weight_i = min(n_i, n_cap) / n_cap`
FinaleMe: `coverage_weight_i = min(n_fragments_i, n_cap) / n_cap`

Balance weight normalizes across cell-type marker groups to prevent dominant cell types from overwhelming the objective.

### 3.3 Gradient

**WGBS mode:**
```
∂L/∂w_j = Σ_i ω_i · φ_i · r_{ij} · [ψ(k_i + α_i) - ψ(n_i - k_i + β_i) - ψ(α_i) + ψ(β_i)]

where r_{i,0} = 0.5 for the unknown component
```

ψ(·) is the digamma function.

**FinaleMe mode:**
```
∂L/∂w_j = Σ_i ω_i · (1 / P(s_i | w, ε)) · ∂P(s_i | w, ε)/∂w_j
```

For called state U:
```
∂P(s_i = U | w, ε)/∂w_j = (1 - ε_U_b) · (1 - r_binary_{ij}) - ε_U_b · r_binary_{ij}
                          (substituting -r_binary for ∂P(true_U)/∂w_j component)
```

Symmetrically for called state M. Gradient is straightforward because the likelihood is a simple linear function of w inside a log.

---

## 4. Bayesian Deconvolution

### 4.1 Model

**WGBS mode:**
```
w ~ Dirichlet(α_0)
k_i | w ~ BetaBinomial(n_i, μ_i(w) · φ_i, (1 - μ_i(w)) · φ_i)

P(w | data) ∝ Dir(w | α_0) · Π_i BetaBinom(k_i | n_i, μ_i(w), φ_i)
```

**FinaleMe mode:**
```
w ~ Dirichlet(α_0)
s_i | w ~ Bernoulli with error rates (§2B.4)

P(w | data) ∝ Dir(w | α_0) · Π_{i: s_i ∈ {U,M}} P(s_i | w, ε)
```

Sampled via MCMC (emcee or NUTS).

### 4.2 Posterior Summaries

Point estimate: posterior mean.
Credible interval: [quantile(α/2), quantile(1 - α/2)].

---

## 5. Per-Cell-Type Reliability P-Values

### 5.1 P_goodness: Marker-Level Goodness-of-Fit

For cell type j, select its D_j most discriminative markers.

**WGBS mode (chi-squared on counts):**
```
χ²_j = Σ_{i ∈ D_j} (k_i - n_i · μ̂_i)² / V̂_i

where:
  μ̂_i = Σ_j ŵ_j · r_{ij} + ŵ_0 · 0.5
  V̂_i = n_i · μ̂_i · (1 - μ̂_i) · (n_i + φ_i) / (1 + φ_i)

p_goodness_j = 1 - F_{χ²}(χ²_j ; df = |D_j| - 1)
```

**FinaleMe mode (binomial concordance test):**

Count how many U/M state calls at this cell type's discriminative markers agree with the expected state under the estimated mixture:

```
For each i ∈ D_j where s_i ∈ {U, M}:
  expected_state_i = U  if P(true_i = U | ŵ) > 0.5,  else M
  concordant_i = I(s_i == expected_state_i)

n_concordant = Σ concordant_i
n_tested = |{i ∈ D_j : s_i ∈ {U, M}}|
p_expected = 1 - mean(ε_U_b, ε_M_b)   (expected accuracy given error rates)

p_goodness_j = BinomialTest(n_concordant, n_tested, p_expected, alternative='less')
```

Interpretation:
- p_goodness > 0.05: estimated proportion is consistent with observed marker data
- p_goodness < 0.01: poor fit, estimate may be unreliable

### 5.2 P_detection: Proportion Stability

From B bootstrap replicates of ŵ_j:

```
p_detection_j = 1 - (1/B) · Σ_{b=1}^{B} I(ŵ_j^{(b)} < ε)

where ε = noise floor (default 0.001)
```

Interpretation:
- p_detection > 0.95: cell type reliably detected above noise
- p_detection < 0.10: proportion is unstable, likely noise

### 5.3 Combined Reliability

| p_goodness | p_detection | Category |
|-----------|-------------|----------|
| > 0.05 | > 0.95 | HIGH |
| > 0.05 | [0.50, 0.95] | MODERATE |
| > 0.05 AND p_det < 0.50, OR p_good < 0.05 AND p_det > 0.50 | | LOW |
| < 0.01 | < 0.50 | UNRELIABLE |

---

## 6. FinaleMe Context-Dependent Binarization Model

### 6.1 Context Bin Assignment

Markers are assigned to B bins based on genomic context features (CpG density, region class, CpG spacing). Default B=8.

### 6.2 Threshold Optimization

For each context bin b, using matched WGBS/WGS cfDNA data:

Define ground truth binarization from WGBS:
```
true_state_i = U   if β_WGBS_i < 0.2
true_state_i = M   if β_WGBS_i > 0.8
true_state_i = intermediate   otherwise
```

Find optimal thresholds by maximizing classification accuracy on non-ambiguous calls:
```
(τ_low_b*, τ_high_b*) = argmax_{τ_l, τ_h}  accuracy(τ_l, τ_h, b)

accuracy(τ_l, τ_h, b) = [Σ_i I(β̂_i < τ_l ∧ true_i = U) + Σ_i I(β̂_i > τ_h ∧ true_i = M)]
                         / [Σ_i I(β̂_i < τ_l) + Σ_i I(β̂_i > τ_h)]
```

Search grid: τ_l ∈ [0.05, 0.30], τ_h ∈ [0.35, 0.70].

### 6.3 Error Rate Estimation

At optimal thresholds:
```
ε_U_b = P(true = M | called U) = Σ_i I(β̂_i < τ_low_b ∧ true_i = M) / Σ_i I(β̂_i < τ_low_b)

ε_M_b = P(true = U | called M) = Σ_i I(β̂_i > τ_high_b ∧ true_i = U) / Σ_i I(β̂_i > τ_high_b)
```

### 6.4 Usability Criterion

A context bin b is usable for deconvolution if:
```
usable_b = (ε_U_b < ε_max) ∧ (ε_M_b < ε_max) ∧ (n_called_U_b ≥ 20) ∧ (n_called_M_b ≥ 20)
```

Where ε_max is the maximum tolerable error rate (default 0.15). Markers in non-usable bins are excluded from the likelihood — this is the data-driven generalization of restricting to CGI+shore.

### 6.5 Cross-Validation

With few matched samples (3–5), CV is performed across regions (not samples):

Chromosome-blocked K-fold CV (default K=10):
```
For each fold f:
  held_out_chroms_f ⊂ {chr1, ..., chr22, chrX}
  train on markers from remaining chromosomes
  evaluate accuracy, ε_U, ε_M on held-out chromosomes

CV_accuracy = mean(accuracy_f across folds)
```

### 6.6 Optimal Bin Count Selection

For each candidate B:
```
score(B) = CV_accuracy(B) · n_usable_markers(B) / n_total_markers
```

Balances classification accuracy with marker coverage. Too few bins (B=4) may lump CGI with shelf; too many (B=16) may have insufficient markers per bin for stable error estimation.

### 6.7 Goodness Metrics

**Training-time:**

| Metric | Definition |
|--------|-----------|
| Per-bin accuracy | Classification accuracy at optimal thresholds |
| Per-bin ε_U, ε_M | Error rates |
| Per-bin fraction called | Fraction of markers assigned U or M (not ambiguous) |
| CV accuracy | Chromosome-blocked cross-validated accuracy |
| Total usable markers | Sum across all usable bins |

**Inference-time (per sample):**

| Metric | Definition |
|--------|-----------|
| Fraction called | Fraction of markers receiving U or M call |
| State distribution shift | KL divergence of U:M ratio vs. training distribution |
| Bin balance | Fraction of usable bins with ≥10 called markers |

---

## 7. Compositional Regression (ILR) for Group Comparison

### 7.1 ILR Transform

The isometric log-ratio maps (K+1)-part composition to K-dimensional unconstrained space:

```
y = Ψᵀ · log(w)
```

Where Ψ is the (K+1) × K contrast matrix (Helmert sub-composition basis).

### 7.2 Per-Cell-Type Regression

For cell type j (using the ILR balance that most directly involves cell type j):

```
y_{j,s} = β_0 + Σ_g β_g · I(group_s = g) + Σ_c β_c · covariate_{c,s} + ε_s
```

### 7.3 Testing Group Effects

**Pairwise contrast (group A vs. group B):**

```
H_0: β_A = β_B

Test statistic: z = (β̂_A - β̂_B) / SE(β̂_A - β̂_B)
p = 2 · (1 - Φ(|z|))
```

Standard errors can be:
- OLS standard errors (fast, assumes normality)
- Bootstrap standard errors (robust, propagates deconvolution uncertainty):
  ```
  SE_boot = sd over bootstrap iterations of (β̂_A^{(b)} - β̂_B^{(b)})
  ```

**Omnibus test (all groups differ?):**

MANOVA on ILR coordinates:
```
H_0: μ_{ILR,g1} = μ_{ILR,g2} = ... = μ_{ILR,gG}

Wilks' Λ, Pillai's trace, or Hotelling-Lawley trace
```

Or per-cell-type Kruskal-Wallis as a simpler alternative.

### 7.4 FDR Correction

Benjamini-Hochberg across all (cell_type × contrast) pairs:
```
Total tests = K · C    where C = number of contrasts
Adjusted p_i = min(1, p_i · (K·C) / rank_i)
```

### 7.5 Why ILR Over Wilcoxon

| Criterion | ILR Regression | Wilcoxon | Direct Likelihood Ratio Test |
|-----------|---------------|----------|------------------------------|
| Respects simplex constraint | Yes | No | No (marker-level) |
| Uses per-sample CI/uncertainty | Yes (bootstrap SE) | No | Partially |
| Handles covariates | Yes (in model) | No (separate adjustment) | Yes |
| Multi-group omnibus | Yes (MANOVA) | No (pairwise only) | No |
| Tests proportions directly | Yes | Yes | No (tests markers) |
| Compositional coherence | Yes | No | N/A |

Direct likelihood ratio testing operates at the marker level and answers "do methylation patterns differ?" rather than "do estimated proportions differ?" — a related but different question. ILR regression on estimated proportions directly answers the proportion question while respecting the compositional nature of the data.

---

## 8. Bayesian Group Comparison

### 8.1 Posterior Probability of Difference

For cell type j, groups A and B:

Draw paired posterior samples of group-level mean proportions:
```
w̄_j^A,(t) = (1/|A|) Σ_{s ∈ A} w_j,s^{(t)}
w̄_j^B,(t) = (1/|B|) Σ_{s ∈ B} w_j,s^{(t)}
```

```
P(A > B for cell type j) = (1/T) Σ_t I(w̄_j^{A,(t)} > w̄_j^{B,(t)})
```

Pseudo p-value for FDR compatibility:
```
p_pseudo = 2 · min(P(A > B), 1 - P(A > B))
```

---

## 9. Imputation

### 9.1 Same-Group Cohort Imputation

For marker i in sample s with n_{i,s} < threshold:
```
μ̂_{i,s} = Σ_{s' ∈ G(s)} v_{s'} · μ_{i,s'} / Σ_{s'} v_{s'}

v_{s'} = n_{i,s'} · I(n_{i,s'} ≥ threshold)
```

Imputed marker handling:

WGBS mode — dispersion inflated to reflect imputation uncertainty:
```
φ_{imputed,i} = φ_i · |G|_eff / (|G|_eff + 1)
```

FinaleMe mode — imputed markers use the cohort consensus state call. If the consensus is ambiguous (no clear U/M majority across same-group samples), the marker is excluded.

### 9.2 Constraints

- Never impute across comparison groups
- Require ≥3 same-group samples with coverage at marker i

---

## 10. Technical Batch Correction

**WGBS mode** — ComBat-style per marker on continuous methylation:
```
Y*_{i,s} = (Y_{i,s} - γ̂_{b(s),i}) / δ̂_{b(s),i} · δ̂_pool + ᾱ_i
```

**FinaleMe mode** — batch correction is applied to predicted beta values *before* binarization. This ensures that batch-shifted predictions don't cause systematic miscalls. The same ComBat-style correction is used on the continuous FinaleMe predictions, then binarization thresholds are applied to the corrected values.

Applied before deconvolution in both modes.

---

## 11. Fragment-Level Likelihood (Ultra-Low Coverage)

For fragment f covering CpGs {c_1, ..., c_L}:
```
P(f | j) = Π_{l=1}^L r_{c_l,j}^{m_l} · (1 - r_{c_l,j})^{1-m_l}

P(f | w) = Σ_j w_j · P(f | j)

L(w) = Σ_f log P(f | w)
```

EM algorithm:
```
E: γ_{f,j} = w_j · P(f | j) / Σ_{j'} w_{j'} · P(f | j')
M: w_j = (1/F) Σ_f γ_{f,j}
```

---

## 12. Summary of P-Values in the System

| P-value | Level | Computed From | Interpretation |
|---------|-------|---------------|----------------|
| p_goodness | Per sample, per cell type | WGBS: χ² on counts; FinaleMe: binomial concordance test | Model fit quality for this cell type |
| p_detection | Per sample, per cell type | Bootstrap proportion stability | Is this cell type detectable above noise? |
| p_contrast | Per cell type, per group contrast | ILR regression Wald test or Bayesian posterior | Does this cell type differ between groups? |
| p_omnibus | Per cell type, all groups | MANOVA or Kruskal-Wallis | Does any group differ for this cell type? |
| adjusted_p | Per cell type × contrast | BH FDR correction | Multiple-testing-corrected significance |

### Binarization Quality Metrics (FinaleMe mode, not p-values but reported alongside)

| Metric | Level | Definition |
|--------|-------|-----------|
| Per-bin accuracy | Per context bin (training) | Classification accuracy at learned thresholds |
| ε_U, ε_M | Per context bin (training) | Misclassification rates |
| CV accuracy | Overall (training) | Chromosome-blocked cross-validated accuracy |
| Fraction called | Per sample (inference) | Fraction of markers receiving U/M call |
| State distribution shift | Per sample (inference) | KL divergence vs. training U:M ratio |

# TOO Module — Mathematical Formulation v2

## 1. Notation

| Symbol | Description |
|--------|-------------|
| M | Number of marker regions |
| K | Number of cell types in reference panel |
| S | Number of samples |
| G | Number of groups |
| w = (w_1, ..., w_K, w_0) | Cell-type proportions + unknown, Σw = 1 |
| w_0 | Unknown/missing cell-type fraction |
| r_{ij} | Reference methylation for marker i, cell type j |
| R | Reference matrix (M × K) |
| k_i | Observed methylated read count at marker i |
| n_i | Total read count at marker i |
| φ_i | Dispersion parameter at marker i |
| μ_i | Expected methylation at marker i under mixture |
| c_i | Genomic context features for marker i |
| ω_i | Marker weight in objective function |

---

## 2. Observation Model

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

### 2.4 Context-Dependent Dispersion (FinaleMe Mode)

For marker i in CpG density bin b:
```
φ_i = exp(d_b) · g(n_i)

g(n_i) = min(n_i, n_cap) / n_cap
```

### 2.5 Reference Uncertainty

With reference coverage n_ref_{ij} available:
```
φ_{effective,i,j} = φ_i · min(1, n_ref_{ij} / n_ref_cap)
```

---

## 3. MLE Deconvolution

### 3.1 Weighted Objective

```
maximize   L(w) = Σ_{i=1}^{M} ω_i · ℓ_i(w)

subject to:
  w_j ≥ 0          for j = 0, 1, ..., K
  Σ_{j=0}^{K} w_j = 1
```

### 3.2 Marker Weights

```
ω_i = min(n_i, n_cap) / n_cap  ·  balance_weight_i
```

Balance weight normalizes across cell-type marker groups to prevent dominant cell types from overwhelming the objective.

### 3.3 Gradient

```
∂L/∂w_j = Σ_i ω_i · φ_i · r_{ij} · [ψ(k_i + α_i) - ψ(n_i - k_i + β_i) - ψ(α_i) + ψ(β_i)]

where r_{i,0} = 0.5 for the unknown component
```

ψ(·) is the digamma function.

---

## 4. Bayesian Deconvolution

### 4.1 Model

```
w ~ Dirichlet(α_0)
k_i | w ~ BetaBinomial(n_i, μ_i(w) · φ_i, (1 - μ_i(w)) · φ_i)

P(w | data) ∝ Dir(w | α_0) · Π_i BetaBinom(k_i | n_i, μ_i(w), φ_i)
```

Sampled via MCMC (emcee or NUTS).

### 4.2 Posterior Summaries

Point estimate: posterior mean.
Credible interval: [quantile(α/2), quantile(1 - α/2)].

---

## 5. Per-Cell-Type Reliability P-Values

### 5.1 P_goodness: Marker-Level Goodness-of-Fit

For cell type j, select its D_j most discriminative markers. Compute:

```
χ²_j = Σ_{i ∈ D_j} (k_i - n_i · μ̂_i)² / V̂_i

where:
  μ̂_i = Σ_j ŵ_j · r_{ij} + ŵ_0 · 0.5
  V̂_i = n_i · μ̂_i · (1 - μ̂_i) · (n_i + φ_i) / (1 + φ_i)

p_goodness_j = 1 - F_{χ²}(χ²_j ; df = |D_j| - 1)
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

## 6. FinaleMe Calibration Model

### 6.1 Calibration Regression

For density bin b, using matched WGBS/FinaleMe samples:
```
logit(β_WGBS) = a_b · logit(β_FinaleMe) + c_b + ε_b
```

### 6.2 Dispersion Estimation

Per-bin dispersion φ_b estimated by maximum likelihood of the beta-binomial on calibration residuals.

### 6.3 Goodness-of-Fit Metrics

**Per-bin R²:**
```
R²_b = 1 - SS_res_b / SS_tot_b
```

**Per-bin MAE:**
```
MAE_b = (1/N_b) Σ |β_calibrated - β_WGBS|
```

**Hosmer-Lemeshow test:**
Within each bin, further divide markers into G decile groups by predicted methylation. For each group g:
```
O_g = observed mean methylation
E_g = predicted mean methylation (calibrated)
HL = Σ_g n_g · (O_g - E_g)² / (E_g · (1 - E_g))

p_HL = 1 - F_{χ²}(HL; df = G - 2)
```

**Cross-validated RMSE:**
Leave-one-sample-out:
```
RMSE_CV = sqrt((1/S) Σ_s (1/M_s) Σ_i (β_calibrated_{i,s}^{(-s)} - β_WGBS_{i,s})²)
```

### 6.4 Optimal Bin Selection

For each candidate B:
```
AIC(B) = -2 · L(B) + 2 · (3B)        # 3 params per bin
BIC(B) = -2 · L(B) + log(N) · (3B)
```

Select B minimizing CV-RMSE; use AIC/BIC as tiebreaker. Report all candidates.

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

| Criterion | ILR Regression | Wilcoxon | Beta-Binomial Regression |
|-----------|---------------|----------|-------------------------|
| Respects simplex constraint | Yes | No | No (marker-level) |
| Uses per-sample CI/uncertainty | Yes (bootstrap SE) | No | Partially |
| Handles covariates | Yes (in model) | No (separate adjustment) | Yes |
| Multi-group omnibus | Yes (MANOVA) | No (pairwise only) | No |
| Tests proportions directly | Yes | Yes | No (tests markers) |
| Compositional coherence | Yes | No | N/A |

Beta-binomial regression operates at the marker level and answers "do methylation patterns differ?" rather than "do estimated proportions differ?" — a related but different question. ILR regression on estimated proportions directly answers the proportion question while respecting the compositional nature of the data.

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

Imputed dispersion (inflated):
```
φ_{imputed,i} = φ_i · |G|_eff / (|G|_eff + 1)
```

### 9.2 Constraints

- Never impute across comparison groups
- Require ≥3 same-group samples with coverage at marker i

---

## 10. Technical Batch Correction

ComBat-style per marker:
```
Y*_{i,s} = (Y_{i,s} - γ̂_{b(s),i}) / δ̂_{b(s),i} · δ̂_pool + ᾱ_i
```

Applied before deconvolution.

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
| p_goodness | Per sample, per cell type | χ² test on discriminative markers | Model fit quality for this cell type |
| p_detection | Per sample, per cell type | Bootstrap proportion stability | Is this cell type detectable above noise? |
| p_contrast | Per cell type, per group contrast | ILR regression Wald test or Bayesian posterior | Does this cell type differ between groups? |
| p_omnibus | Per cell type, all groups | MANOVA or Kruskal-Wallis | Does any group differ for this cell type? |
| adjusted_p | Per cell type × contrast | BH FDR correction | Multiple-testing-corrected significance |
| p_HL (calibration) | Per density bin | Hosmer-Lemeshow test | Calibration quality per bin |
| p_KS (calibration QC) | Per sample | KS test on calibration residuals | Sample-level calibration quality |

# FinaleMe TOO - Hierarchical Binarization Model (Unified MLE + Bayesian)

## 1. Purpose

This document specifies a single hierarchical observation model for FinaleMe binarization mode that is shared by:

- point estimation in `MLEDeconvolver`
- posterior inference in `BayesianDeconvolver`

The goal is to replace composite mixing of separate terms with one coherent per-marker joint likelihood that uses both:

- count evidence (`k_i`, `n_i`)
- state-call evidence (`called_state_i`)

while preserving practical robustness, speed, and backward compatibility.

## 2. Problem Statement

Current FinaleMe binarization mode uses both state and count information, but in a weighted hybrid objective. This is practical and works well, yet it is not a fully generative joint model.

Main limitation:

- `called_state` is derived from methylation signal related to the same underlying fragments as `k/n`.
- Naively multiplying independent likelihood terms can double-count evidence.

Design objective:

- introduce one explicit latent-variable model for marker-level methylation and state-calling,
- use the same log-likelihood implementation in both MLE and Bayesian paths,
- keep a dependency-control parameter for robustness.

## 3. Scope

In scope:

- FinaleMe mode only (`MeasurementMode.FINALEME`)
- marker-aggregated deconvolution (`finaleme-too run`)
- MLE (`SLSQP`) and Bayesian (`emcee`) deconvolution
- training and serialization extensions for binarization parameters

Out of scope (first implementation):

- WGBS mode changes
- group-level testing redesign
- fragment-level model redesign

## 4. Notation

For marker `i = 1..M` and cell types `j = 1..K`, with unknown component `j = K+1`:

- `w_j >= 0`, `sum_j w_j = 1` : mixture proportions
- `r_ij` : reference methylation for marker `i`, cell type `j`
- `mu_i(w) = sum_j r_ij * w_j` : expected marker methylation under mixture
- `k_i`, `n_i` : methylated count and total count
- `phi_i` : dispersion (beta-binomial precision)
- `b_i` : context bin index
- `c_i` : observed state call in `{U, M, A, X}`
  - `U`: unmethylated call
  - `M`: methylated call
  - `A`: ambiguous
  - `X`: excluded/unusable

## 5. Hierarchical Generative Model

### 5.1 Latent methylation probability

For each marker:

`p_i | w ~ Beta(alpha_i, beta_i)`

where:

- `alpha_i = mu_i(w) * phi_i`
- `beta_i  = (1 - mu_i(w)) * phi_i`

### 5.2 Counts channel

`k_i | p_i, n_i ~ Binomial(n_i, p_i)`

Marginally this is beta-binomial.

### 5.3 State-call channel

Observed call is generated from latent `p_i` and bin-specific call model:

`c_i | p_i, b_i ~ Categorical( q_b_i(p_i) )`

where `q_b(p)` returns probabilities for `(U, M, A, X)`.

### 5.4 Joint marker likelihood

Per marker:

`L_i(w) = integral_0^1 Binomial(k_i | n_i, p) * Beta(p | alpha_i, beta_i) * q_{b_i,c_i}(p) dp`

Total weighted log-likelihood:

`log L(w) = sum_i omega_i * log(L_i(w))`

This is the one shared likelihood for both MLE and Bayesian.

## 6. Call Model Parameterization (q_b)

Two staged parameterizations are specified.

### 6.1 V1 (recommended initial rollout): piecewise-constant calibrated call model

For each bin `b`, train:

- `tau_low_b`, `tau_high_b`
- confusion/error probabilities for low/mid/high p zones

Define zones by latent `p`:

- low: `p < tau_low_b`
- mid: `tau_low_b <= p <= tau_high_b`
- high: `p > tau_high_b`

Then `q_b(p)` is piecewise constant by zone using trained confusion matrix rows.

Pros:

- easy to train from matched data
- stable and interpretable
- compatible with existing threshold/error-rate framework

### 6.2 V2 (future): smooth monotonic call model

Replace piecewise constants with monotonic logistic/spline probabilities:

- `q_b,U(p)` decreases with `p`
- `q_b,M(p)` increases with `p`
- `q_b,A(p)` peaked near transition region
- `q_b,X(p)` from unusable/noise regime

Pros:

- differentiable and more realistic
- better calibration in transition regions

## 7. Dependence Control (double-count mitigation)

Because `c_i` and `k_i` arise from related evidence, introduce:

- `call_weight` in `[0, 1]`

Adjusted marker likelihood contribution:

`log L_i_adj(w) = log integral [Binom * Beta * q^{call_weight}] dp`

Practical interpretation:

- `call_weight = 1.0` : full joint call-channel contribution
- `call_weight < 1.0` : temper call channel to reduce overconfidence
- tune by chromosome-blocked CV on matched data

This keeps a single hierarchical model while controlling dependence bias.

## 8. Inference

### 8.1 Shared log-likelihood engine

Create one vectorized implementation:

- input: `w`, marker arrays (`k,n,phi,bin,call,weights`), reference matrix
- output: `log L(w)`
- optional output: gradient for MLE

Both MLE and Bayesian must call this same engine.

### 8.2 MLE path

Objective:

`maximize log L(w)`
subject to simplex constraints.

Optimizer:

- continue SLSQP
- initialize from current MLE solution or uniform

Gradient strategy:

- Phase 1: robust finite-difference gradient around shared likelihood engine
- Phase 2: analytic/semi-analytic gradient for speed

### 8.3 Bayesian path

Posterior:

`P(w | data) proportional to Dirichlet(w | alpha0) * L(w)`

- keep current softmax parameterization in unconstrained space
- use same shared log-likelihood engine
- posterior mean can be point estimate when `model.deconvolution == bayesian`

## 9. Numerical Integration Strategy

The marker integral has no simple closed form for general `q_b(p)`.

Recommended implementation:

- Gauss-Jacobi quadrature on `[0,1]` (or fixed transformed Gauss-Legendre)
- default nodes per marker: 24 (configurable)
- vectorized over markers

Performance tricks:

- cache quadrature nodes/weights globally
- precompute bin-specific `q_b(p_nodes)` lookup tables
- chunk markers to control memory

Stability:

- clip probabilities with `eps = 1e-15`
- log-sum-exp style accumulation when needed

## 10. Data Model and Serialization Changes

Extend binarization parameter schema (new version):

- `model_version: "hierarchical_v1"`
- existing fields retained (`tau_low`, `tau_high`, `usable`, etc.)
- new fields:
  - `call_model_type` (`piecewise_v1`)
  - `call_model_params` (per-bin zone confusion matrices)
  - `call_weight_default`

Backward compatibility:

- old params load in compatibility mode (`legacy_v3`)
- hierarchical path requires new fields or auto-upgrade defaults

### 10.1 JSON Schema for `call_model_params` (exact draft)

The block below is a machine-readable draft for the new hierarchical fields.
It is designed to sit at the top level of the existing binarization JSON
(`tau_low`, `tau_high`, `eps_U`, `eps_M`, `usable`, etc. remain unchanged).
Canonical file path in this repo:
`design/schema/binarization_hierarchical_v1.schema.json`.

```json
{
  "$schema": "https://json-schema.org/draft/2020-12/schema",
  "$id": "https://finaleme-too/schema/binarization_hierarchical_v1.json",
  "title": "FinaleMe TOO Binarization Parameters (hierarchical_v1)",
  "type": "object",
  "required": [
    "model_version",
    "n_bins",
    "call_model_type",
    "call_weight_default",
    "call_model_params"
  ],
  "properties": {
    "model_version": {
      "type": "string",
      "const": "hierarchical_v1"
    },
    "n_bins": {
      "type": "integer",
      "minimum": 1
    },
    "call_model_type": {
      "type": "string",
      "const": "piecewise_v1"
    },
    "call_weight_default": {
      "type": "number",
      "minimum": 0.0,
      "maximum": 1.0
    },
    "call_model_params": {
      "type": "object",
      "required": ["states", "zones", "per_bin"],
      "properties": {
        "states": {
          "type": "array",
          "const": ["U", "M", "A", "X"]
        },
        "zones": {
          "type": "array",
          "const": ["low", "mid", "high"]
        },
        "per_bin": {
          "type": "array",
          "minItems": 1,
          "items": {
            "type": "object",
            "required": ["bin_index", "zone_confusion"],
            "properties": {
              "bin_index": {
                "type": "integer",
                "minimum": 0
              },
              "zone_confusion": {
                "type": "object",
                "required": ["low", "mid", "high"],
                "properties": {
                  "low": {
                    "$ref": "#/$defs/state_prob_row"
                  },
                  "mid": {
                    "$ref": "#/$defs/state_prob_row"
                  },
                  "high": {
                    "$ref": "#/$defs/state_prob_row"
                  }
                },
                "additionalProperties": false
              },
              "n_train_low": {
                "type": "integer",
                "minimum": 0
              },
              "n_train_mid": {
                "type": "integer",
                "minimum": 0
              },
              "n_train_high": {
                "type": "integer",
                "minimum": 0
              },
              "fitted": {
                "type": "boolean"
              }
            },
            "additionalProperties": false
          }
        }
      },
      "additionalProperties": false
    }
  },
  "$defs": {
    "state_prob_row": {
      "type": "object",
      "required": ["U", "M", "A", "X"],
      "properties": {
        "U": {
          "type": "number",
          "minimum": 0.0,
          "maximum": 1.0
        },
        "M": {
          "type": "number",
          "minimum": 0.0,
          "maximum": 1.0
        },
        "A": {
          "type": "number",
          "minimum": 0.0,
          "maximum": 1.0
        },
        "X": {
          "type": "number",
          "minimum": 0.0,
          "maximum": 1.0
        }
      },
      "additionalProperties": false
    }
  },
  "additionalProperties": true
}
```

Additional validation rules to enforce in loader code (not expressible
cleanly in plain JSON schema):

- `len(call_model_params.per_bin) == n_bins`
- `bin_index` values are unique and exactly cover `[0, n_bins-1]`
- each probability row sums to `1.0 +/- 1e-6`
- if a bin has insufficient training support, `fitted=false` and fallback
  policy is applied during inference

### 10.2 Example payload (`hierarchical_v1`)

```json
{
  "binarization_version": "2.0",
  "model_version": "hierarchical_v1",
  "n_bins": 2,
  "n_region_classes": 1,
  "density_sub_bins_per_class": 2,
  "region_class_order": ["open_sea"],
  "density_edges": [
    [0.0, 0.02, 1.0],
    [0.0, 0.02, 1.0]
  ],
  "tau_low": [0.10, 0.12],
  "tau_high": [0.90, 0.88],
  "eps_U": [0.03, 0.05],
  "eps_M": [0.02, 0.04],
  "usable": [true, true],
  "call_model_type": "piecewise_v1",
  "call_weight_default": 0.70,
  "call_model_params": {
    "states": ["U", "M", "A", "X"],
    "zones": ["low", "mid", "high"],
    "per_bin": [
      {
        "bin_index": 0,
        "zone_confusion": {
          "low": {
            "U": 0.9500,
            "M": 0.0100,
            "A": 0.0350,
            "X": 0.0050
          },
          "mid": {
            "U": 0.2600,
            "M": 0.2400,
            "A": 0.4700,
            "X": 0.0300
          },
          "high": {
            "U": 0.0150,
            "M": 0.9550,
            "A": 0.0250,
            "X": 0.0050
          }
        },
        "n_train_low": 1240,
        "n_train_mid": 380,
        "n_train_high": 1195,
        "fitted": true
      },
      {
        "bin_index": 1,
        "zone_confusion": {
          "low": {
            "U": 0.9100,
            "M": 0.0200,
            "A": 0.0600,
            "X": 0.0100
          },
          "mid": {
            "U": 0.3000,
            "M": 0.2800,
            "A": 0.3900,
            "X": 0.0300
          },
          "high": {
            "U": 0.0250,
            "M": 0.9150,
            "A": 0.0500,
            "X": 0.0100
          }
        },
        "n_train_low": 840,
        "n_train_mid": 260,
        "n_train_high": 790,
        "fitted": true
      }
    ]
  }
}
```

For production (`n_bins=8` or larger), keep the same structure and provide
one `per_bin` entry per bin.

## 11. Config and CLI Additions

Config (`TOOConfig.model` and `TOOConfig.binarization`):

- `finaleme_observation_model: legacy_v3 | hierarchical_v1`
- `call_weight: float` (default from params)
- `quadrature_points: int` (default 24)

CLI additions:

- `--finaleme-observation-model hierarchical_v1`
- `--call-weight 0.7` (optional override)
- `--quadrature-points 24`

## 12. Pipeline Integration Plan

### 12.1 Core code touch points

- `finaleme_too/core/observation_model_binarization.py`
  - include call-model-ready fields in observation object
- `finaleme_too/core/deconvolution.py`
  - add shared hierarchical likelihood engine
  - MLE and Bayesian both call the same function
- `finaleme_too/preprocessing/binarization.py`
  - train/load/save new call-model parameters
- `finaleme_too/config.py`, `finaleme_too/cli.py`
  - new config/CLI options

### 12.2 Reliability and QC

Keep existing reliability outputs for compatibility, and add:

- posterior predictive calibration diagnostics (optional report)
- count-channel vs call-channel contribution summary

## 13. Validation Plan

### 13.1 Unit tests

- likelihood sanity checks on synthetic markers
- MLE recovers known mixtures under simulated data
- Bayesian posterior mean close to MLE under weak prior
- backward compatibility with legacy params

### 13.2 Calibration tests

On matched WGBS/FinaleMe data:

- calibration curve of predicted vs observed call frequencies
- CI coverage of known synthetic truth
- compare unknown fraction stability vs legacy

### 13.3 Acceptance criteria

For hierarchical_v1 to become default:

- no regression in core deconvolution tests
- improved or equal RMSE on matched benchmark cohorts
- improved uncertainty calibration (credible interval coverage)
- runtime increase <= 1.8x at default quadrature points

## 14. Rollout Strategy

Phase 0 (design + schema):

- finalize equations, schema, and config keys

Phase 1 (feature flag implementation):

- implement `hierarchical_v1` behind explicit flag
- keep `legacy_v3` default

Phase 2 (benchmark + tuning):

- tune `call_weight` and quadrature defaults by CV
- add tutorial notes and interpretation guidance

Phase 3 (default switch):

- switch default to `hierarchical_v1` only after acceptance criteria pass
- keep legacy mode for reproducibility

## 15. Risks and Mitigations

Risk: overconfidence from dependence between call and count channels.

- Mitigation: `call_weight` tempering + CV tuning + posterior predictive checks.

Risk: runtime overhead from per-marker integration.

- Mitigation: vectorized quadrature, cached lookups, marker chunking.

Risk: identifiability at ultra-low coverage.

- Mitigation: stronger prior, marker QC, preserve unknown component, coverage-tier-specific safeguards.

## 16. Open Decisions

1. Should `call_weight` be global, tier-specific, or bin-specific?
2. Should ambiguous/excluded calls contribute weakly (instead of being dropped)?
3. Should V1 use finite-difference gradient only, or implement analytic gradient immediately?
4. Should hierarchical training be integrated into `train-binarization` or added as a separate subcommand?

## 17. Summary

This spec defines a single hierarchical FinaleMe observation model where both count and state-call evidence are modeled jointly through a latent methylation probability per marker. The same likelihood is used in MLE and Bayesian deconvolution, with explicit dependence control (`call_weight`) and a practical rollout path that preserves backward compatibility.

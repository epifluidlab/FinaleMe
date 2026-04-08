"""finaleme_too — Tissue-of-origin deconvolution for cfDNA methylation data.

Dual observation model:

* **WGBS / EM-seq mode** — beta-binomial likelihood on per-marker read counts
  with per-sample dispersion MLE, coverage-weighted objective, three coverage
  tiers, bootstrap and Bayesian uncertainty, per-cell-type reliability
  p-values, and downstream compositional regression in ILR space.
* **FinaleMe mode** — context-dependent binarization with learned per-bin
  thresholds ``(τ_low, τ_high)`` and error rates ``(ε_U, ε_M)`` classifying
  predicted methylation into U / M / Ambiguous / Excluded states per
  region-class × density-sub-bin, fed into a binomial-with-error-rates SLSQP
  solver.

See design/TOO_ARCHITECTURE_v3.md and design/TOO_MATH_FORMULATION_v3.md for
the specification.
"""

from finaleme_too._version import __version__

__all__ = ["__version__"]

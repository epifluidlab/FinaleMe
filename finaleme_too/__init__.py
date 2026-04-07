"""finaleme_too — Tissue-of-origin deconvolution for cfDNA methylation data.

Implements a beta-binomial observation model with three coverage tiers, dual
WGBS/FinaleMe modes, FinaleMe calibration, bootstrap and Bayesian uncertainty,
per-cell-type reliability p-values, and downstream compositional regression.

See design/TOO_ARCHITECTURE_v2.md and design/TOO_MATH_FORMULATION_v2.md for the
specification.
"""

from finaleme_too._version import __version__

__all__ = ["__version__"]

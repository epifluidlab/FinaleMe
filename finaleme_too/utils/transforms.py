"""Compositional data transforms (CLR / ALR / ILR) and their inverses.

Math doc §7.1: ILR uses a Helmert sub-composition basis. We use the Egozcue
orthonormal basis as implemented by the standard formula.
"""

from __future__ import annotations

import numpy as np


def _safe_log(x: np.ndarray, eps: float = 1e-12) -> np.ndarray:
    return np.log(np.maximum(x, eps))


def clr_transform(x: np.ndarray) -> np.ndarray:
    """Centered log-ratio transform.

    Input shape: (..., D). Output shape: (..., D). Sum of CLR coords is 0.
    """
    log_x = _safe_log(x)
    return log_x - np.mean(log_x, axis=-1, keepdims=True)


def clr_inverse(y: np.ndarray) -> np.ndarray:
    exp_y = np.exp(y)
    return exp_y / np.sum(exp_y, axis=-1, keepdims=True)


def helmert_basis(d: int) -> np.ndarray:
    """Egozcue/Helmert orthonormal basis for ILR with d parts → (d-1) coords.

    Returns a (d, d-1) matrix V such that V.T @ clr(x) gives the ILR coords.
    """
    if d < 2:
        raise ValueError("ILR requires at least 2 parts")
    V = np.zeros((d, d - 1), dtype=np.float64)
    for i in range(d - 1):
        n_i = i + 1  # parts in the numerator group
        # column i = scale * [1/n_i, ..., 1/n_i,    -1, 0, ..., 0]
        # the first n_i entries are 1/n_i, the (n_i+1)-th entry is -1, rest 0
        scale = np.sqrt(n_i / (n_i + 1.0))
        V[:n_i, i] = 1.0 / n_i
        V[n_i, i] = -1.0
        V[:, i] *= scale
    return V


def ilr_transform(x: np.ndarray) -> np.ndarray:
    """Isometric log-ratio transform (Egozcue 2003).

    Input shape (..., D). Output shape (..., D-1).
    """
    arr = np.asarray(x, dtype=np.float64)
    d = arr.shape[-1]
    V = helmert_basis(d)
    return clr_transform(arr) @ V


def ilr_inverse(y: np.ndarray, n_parts: int | None = None) -> np.ndarray:
    """Inverse ILR transform: (D-1)-dim coords → D-part composition."""
    arr = np.asarray(y, dtype=np.float64)
    if n_parts is None:
        n_parts = arr.shape[-1] + 1
    V = helmert_basis(n_parts)
    return clr_inverse(arr @ V.T)


def alr_transform(x: np.ndarray, ref_index: int = -1) -> np.ndarray:
    """Additive log-ratio with the chosen reference part."""
    arr = np.asarray(x, dtype=np.float64)
    log_x = _safe_log(arr)
    ref = log_x[..., ref_index : ref_index + 1] if ref_index != -1 else log_x[..., -1:]
    return np.delete(log_x - ref, ref_index, axis=-1)


__all__ = [
    "alr_transform",
    "clr_inverse",
    "clr_transform",
    "helmert_basis",
    "ilr_inverse",
    "ilr_transform",
]

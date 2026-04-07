"""Thin wrapper around joblib for sample-level parallelism."""

from __future__ import annotations

from typing import Callable, Iterable, TypeVar

from joblib import Parallel, delayed

T = TypeVar("T")
R = TypeVar("R")


def parallel_map(
    func: Callable[[T], R],
    items: Iterable[T],
    n_jobs: int = 1,
    backend: str = "loky",
    verbose: int = 0,
) -> list[R]:
    """Apply ``func`` to every element of ``items`` in parallel.

    Falls back to a serial loop if ``n_jobs == 1`` to avoid joblib overhead
    and to make debugging easier.
    """
    items_list = list(items)
    if n_jobs <= 1 or len(items_list) <= 1:
        return [func(item) for item in items_list]
    return Parallel(n_jobs=n_jobs, backend=backend, verbose=verbose)(
        delayed(func)(item) for item in items_list
    )


__all__ = ["parallel_map"]

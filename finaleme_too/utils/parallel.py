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

    Falls back to a serial loop if ``n_jobs == 1`` or there is at most one
    item (avoids joblib overhead and makes debugging easier).

    Backend choice:
      * ``loky`` (default) — process pool. Best when each task is large
        and the arguments are small (e.g. per-sample deconvolution).
      * ``threading`` — thread pool. Best when the inner work is pure
        numpy / scipy that releases the GIL and the shared arguments are
        large (e.g. CV folds over multi-megabyte arrays). No pickling
        overhead, so big arrays are shared by reference.
    """
    items_list = list(items)
    if n_jobs <= 1 or len(items_list) <= 1:
        return [func(item) for item in items_list]
    return Parallel(n_jobs=n_jobs, backend=backend, verbose=verbose)(
        delayed(func)(item) for item in items_list
    )


__all__ = ["parallel_map"]

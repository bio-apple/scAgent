"""Parallel map for marker scoring / per-sample jobs. joblib if present, else serial."""

from __future__ import annotations

from collections.abc import Callable, Iterable
from typing import TypeVar

from scagent.config import performance_params
from scagent.logutil import get_logger

log = get_logger("parallel")
T = TypeVar("T")
R = TypeVar("R")


def n_jobs(override: int | None = None) -> int:
    if override is not None:
        return int(override)
    return int(performance_params()["n_jobs"])


def apply_scanpy_n_jobs(jobs: int | None = None) -> int:
    j = n_jobs(jobs)
    try:
        import scanpy as sc

        sc.settings.n_jobs = j
    except Exception:
        pass
    return j


def map_parallel(fn: Callable[[T], R], items: Iterable[T], *, jobs: int | None = None) -> list[R]:
    seq = list(items)
    j = n_jobs(jobs)
    if j == 1 or len(seq) <= 1:
        return [fn(x) for x in seq]
    try:
        from joblib import Parallel, delayed

        nj = j if j > 0 else -1
        log.info("joblib Parallel n_jobs=%s n_items=%s", nj, len(seq))
        return list(Parallel(n_jobs=nj, prefer="threads")(delayed(fn)(x) for x in seq))
    except ImportError:
        log.debug("joblib missing; serial map")
        return [fn(x) for x in seq]

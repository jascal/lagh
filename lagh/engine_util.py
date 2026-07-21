"""STLSQ support discovery -- the scaling doctrine's workhorse."""

from __future__ import annotations

import numpy as np

from .base import lstsq

THRESHOLDS = np.logspace(-8, 0, 20)


def stlsq_supports(M: np.ndarray, y: np.ndarray,
                   thresholds=THRESHOLDS) -> set[tuple[int, ...]]:
    out: set[tuple[int, ...]] = set()
    n_terms = M.shape[1]
    for thr in thresholds:
        active = np.ones(n_terms, dtype=bool)
        for _ in range(12):
            c = lstsq(M[:, active], y)
            if c is None:
                break
            full = np.zeros(n_terms)
            full[active] = c
            keep = np.abs(full) >= thr
            if keep.sum() == 0 or np.array_equal(keep, active):
                active = keep if keep.sum() else active
                break
            active = keep
        idx = tuple(np.flatnonzero(active).tolist())
        if idx:
            out.add(idx)
    return out

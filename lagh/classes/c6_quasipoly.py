"""C6: the quasi-polynomial tier, promoted into the curriculum.

Reached by escalation after C1-C5 (float) fail on an INTEGER-LATTICE target -- one
whose inputs and outputs are all integers. Float tiers structurally cannot certify
exact-integer data (a float fit never hits machine-precision integer equality), so
this exact-arithmetic tier is the honest terminus for integer laws.

Currently 1-D in the dilation parameter (Ehrhart L_P(t)); the recovery itself is in
lagh/quasipoly.py. This module is the curriculum adapter: detection + invocation.
"""

from __future__ import annotations

import numpy as np

from ..quasipoly import recover

TIER = 6


def is_integer_lattice(X: np.ndarray, y: np.ndarray) -> bool:
    """1-D input, all inputs and outputs integer-valued to machine precision."""
    X = np.asarray(X, float)
    y = np.asarray(y, float).ravel()
    if X.ndim != 2 or X.shape[1] != 1:
        return False
    return bool(np.allclose(X, np.round(X)) and np.allclose(y, np.round(y)))


def recover_integer(ts_all, Ls_all, *, period_max: int = 12, degree_max: int = 4):
    """Pool integer (t, L) pairs and hand to the exact quasi-polynomial recovery.
    The recovery does its own per-class self-split, so pre-splitting is neither
    needed nor wanted."""
    order = np.argsort(ts_all)
    ts = [int(round(ts_all[i])) for i in order]
    Ls = [int(round(Ls_all[i])) for i in order]
    # de-duplicate t (the pooled splits may repeat)
    seen, ts_u, Ls_u = set(), [], []
    for t, L in zip(ts, Ls):
        if t not in seen:
            seen.add(t)
            ts_u.append(t)
            Ls_u.append(L)
    return recover(ts_u, Ls_u, period_max=period_max, degree_max=degree_max)

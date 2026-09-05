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

from ..certify import MACHINE_REL
from ..quasipoly import recover

TIER = 6


def lattice_deviation(X: np.ndarray, y: np.ndarray) -> tuple[float, float]:
    """(max |x - round(x)|, max |y - round(y)|): the quantization the integer
    tier would apply, reported so a certificate can state it."""
    X = np.asarray(X, float)
    y = np.asarray(y, float).ravel()
    dx = float(np.max(np.abs(X - np.round(X)))) if X.size else 0.0
    dy = float(np.max(np.abs(y - np.round(y)))) if y.size else 0.0
    return dx, dy


def _on_lattice(v: np.ndarray) -> bool:
    v = np.asarray(v, float).ravel()
    if not np.all(np.isfinite(v)):
        return False
    return bool(np.all(np.abs(v - np.round(v)) <= MACHINE_REL * np.maximum(1.0, np.abs(v))))


def is_integer_lattice(X: np.ndarray, y: np.ndarray) -> bool:
    """1-D input, all inputs and outputs integer-valued to MACHINE precision:
    |v - round(v)| <= MACHINE_REL * max(1, |v|), the float64-representation
    term of the engine's own error model -- an explicit, scale-aware
    eligibility, never a default tolerance.

    `np.allclose` (rtol 1e-5, atol 1e-8) was the measured hole (jascal/lagh#3):
    at |y| ~ 1e5 its tolerance is ~1, so observations visibly off the lattice
    (by 0.1) were declared integer, ROUNDED by `recover_integer`, and the
    rounded data certified at alpha ~ 1e-211 while the law missed the real
    observations by 0.098 against their declared band of 2e-8. Eligibility is
    now machine-precision only, and the engine re-checks any C6 law against
    the ORIGINAL observations at the declared band before certifying."""
    X = np.asarray(X, float)
    y = np.asarray(y, float).ravel()
    if X.ndim != 2 or X.shape[1] != 1:
        return False
    return _on_lattice(X[:, 0]) and _on_lattice(y)


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

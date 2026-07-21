"""C3: power laws by log-log fit, exponents snapped to small rationals."""

from __future__ import annotations

from fractions import Fraction

import numpy as np
import sympy as sp

from ..base import Candidate, lstsq

TIER = 3


def candidates(ctx) -> list[Candidate]:
    X, y = ctx.X_fit, np.asarray(ctx.y_fit, float).ravel()
    if np.any(y <= 0) or np.any(X <= 0):
        return []
    L = np.column_stack([np.ones(len(X)), np.log(X)])
    c = lstsq(L, np.log(y))
    if c is None:
        return []
    out = []
    for cap in (1, 2, 4):
        exps = [Fraction(float(a)).limit_denominator(cap) for a in c[1:]]
        with np.errstate(all="ignore"):
            base = np.prod([X[:, i] ** float(e) for i, e in enumerate(exps)], axis=0)
        if not np.all(np.isfinite(base)) or np.sum(base**2) == 0:
            continue
        k = float(np.dot(base, y) / np.dot(base, base))
        expr = sp.Float(k)
        for i, e in enumerate(exps):
            expr = expr * ctx.syms[i] ** sp.Rational(e.numerator, e.denominator)
        out.append(Candidate(expr=expr, complexity=int(sp.count_ops(expr)),
                             channel="c3-powerlaw"))
    return out

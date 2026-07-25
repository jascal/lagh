"""C3: power laws by log-log fit, exponents snapped to small rationals."""

from __future__ import annotations

from fractions import Fraction

import numpy as np
import sympy as sp

from ..base import Candidate, lstsq

TIER = 3


def candidates(ctx) -> list[Candidate]:
    X, y = ctx.X_fit, np.asarray(ctx.y_fit, float).ravel()
    # CAP-N (LLMSRBENCH_DEV.md): an all-negative target is a monomial with a
    # constant sign -- fit |y|, restore the sign. 8/19 mined benchmark cells
    # missed ONLY for this.
    sign = 1
    if np.all(y < 0):
        sign, y = -1, -y
    if np.any(y <= 0) or np.any(X <= 0):
        return []
    L = np.column_stack([np.ones(len(X)), np.log(X)])
    c = lstsq(L, np.log(y))
    if c is None:
        return []
    out = []
    # CAP-A: denominator caps for exponent snapping. Extended 4 -> {3,5,10} so
    # denom-5/10/3 rationals snap exactly (x^3.4=17/5, x^-0.3=-3/10, x^-10/3);
    # capped at 4 they mis-snapped and no power law certified (hooke hard cells).
    # Each cap adds ONE checked candidate (no combinatorial inflation), and the
    # exponent is data-driven from the log-log slope then verified -- a wrong snap
    # fails certification, so the checker bounds the added exponents.
    for cap in (1, 2, 3, 4, 5, 10):
        exps = [Fraction(float(a)).limit_denominator(cap) for a in c[1:]]
        with np.errstate(all="ignore"):
            base = np.prod([X[:, i] ** float(e) for i, e in enumerate(exps)], axis=0)
        if not np.all(np.isfinite(base)) or np.sum(base**2) == 0:
            continue
        k = float(np.dot(base, y) / np.dot(base, base))
        expr = sp.Float(sign * k)
        for i, e in enumerate(exps):
            expr = expr * ctx.syms[i] ** sp.Rational(e.numerator, e.denominator)
        out.append(Candidate(expr=expr, complexity=int(sp.count_ops(expr)),
                             channel="c3-powerlaw"))
    return out

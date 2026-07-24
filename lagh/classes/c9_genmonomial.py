"""CAP-F: generalized monomials -- transcendental(input) x power products.

One restricted structure:  y = C * prod_j f_j(x_j)^{p_j},  f_j in {x, e^x, ln x}.

Taking logs it is a LINEAR fit per f-assignment:
    log y = c0 + sum_j p_j * g_j(x_j),   g_j in {log x_j, x_j, log(ln x_j)}.

Enumerating assignments is bounded (3^dim - 1 fits, dim <= 3 by registered bound;
the all-identity assignment is plain C3 and skipped). Exponents snap to small
rationals, the coefficient goes through the escalating snap, and every produced law
is verified by the unchanged checker. Reaches e^gamma- and ln(gamma)-factor laws
(sound-speed hard cells) that no pure power law can express.
"""

from __future__ import annotations

from fractions import Fraction
from itertools import product

import numpy as np
import sympy as sp

from ..base import Candidate, lstsq, snap

TIER = 4
_EXP_CAP = 10
_FIT_TOL = 1e-8


def candidates(ctx) -> list[Candidate]:
    X = np.asarray(ctx.X_fit, float)
    y = np.asarray(ctx.y_fit, float).ravel()
    dim = X.shape[1]
    if dim > 3 or len(X) < 8 or np.any(X <= 0) or np.any(y <= 0):
        return []
    ly = np.log(y)
    # per input: the usable g-columns (None = f inadmissible on this data)
    G = {"id": [np.log(X[:, j]) for j in range(dim)],
         "exp": [X[:, j] if np.all(X[:, j] < 100) else None for j in range(dim)],
         "log": [np.log(np.log(X[:, j])) if np.all(X[:, j] > 1.0 + 1e-9) else None
                 for j in range(dim)]}
    out: list[Candidate] = []
    for assign in product(("id", "exp", "log"), repeat=dim):
        if all(a == "id" for a in assign):
            continue                      # plain C3's job
        cols = [G[a][j] for j, a in enumerate(assign)]
        if any(c is None for c in cols):
            continue
        L = np.column_stack([np.ones(len(X))] + cols)
        c = lstsq(L, ly)
        if c is None:
            continue
        resid = float(np.sqrt(np.mean((L @ c - ly) ** 2)))
        if resid > _FIT_TOL:
            continue
        exps = [Fraction(float(a)).limit_denominator(_EXP_CAP) for a in c[1:]]
        C = snap(float(np.exp(c[0])))
        expr: sp.Expr = sp.Rational(C.numerator, C.denominator)
        for j, (a, e) in enumerate(zip(assign, exps)):
            if e == 0:
                continue
            base = {"id": ctx.syms[j], "exp": sp.exp(ctx.syms[j]),
                    "log": sp.log(ctx.syms[j])}[a]
            expr *= base ** sp.Rational(e.numerator, e.denominator)
        out.append(Candidate(expr=expr, complexity=int(sp.count_ops(expr)) + 2,
                             channel=f"c9-genmono-{'-'.join(assign)}"))
    return out

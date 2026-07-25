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
_FIT_TOL = 1e-3   # loose: the checker decides (1e-8 was dead at sigma_rep=1e-4)
_CONJ_TOL = 5e-2      # loose: a conjecture is a labeled guess, not a certificate


def candidates(ctx) -> list[Candidate]:
    X = np.asarray(ctx.X_fit, float)
    y = np.asarray(ctx.y_fit, float).ravel()
    dim = X.shape[1]
    sign = 1                          # CAP-N: constant-sign closure, as in C3
    if np.all(y < 0):
        sign, y = -1, -np.asarray(y, float)
    if dim > 3 or len(X) < 8 or np.any(X <= 0) or np.any(y <= 0):
        return []
    ly = np.log(y)
    # per input: the usable g-columns (None = f inadmissible on this data)
    G = {"id": [np.log(X[:, j]) for j in range(dim)],
         "exp": [X[:, j] if np.all(X[:, j] < 100) else None for j in range(dim)],
         "log": [np.log(np.log(X[:, j])) if np.all(X[:, j] > 1.0 + 1e-9) else None
                 for j in range(dim)],
         # CAP-T companion: exp(p/x_j) factors (Arrhenius) -- the g-column is 1/x
         "invexp": [1.0 / X[:, j] for j in range(dim)]}
    out: list[Candidate] = []
    for assign in product(("id", "exp", "log", "invexp"), repeat=dim):
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
        expr: sp.Expr = sign * sp.Rational(C.numerator, C.denominator)
        for j, (a, e) in enumerate(zip(assign, exps)):
            if e == 0:
                continue
            base = {"id": ctx.syms[j], "exp": sp.exp(ctx.syms[j]),
                    "log": sp.log(ctx.syms[j]),
                    "invexp": sp.exp(1 / ctx.syms[j])}[a]
            expr *= base ** sp.Rational(e.numerator, e.denominator)
        out.append(Candidate(expr=expr, complexity=int(sp.count_ops(expr)) + 2,
                             channel=f"c9-genmono-{'-'.join(assign)}"))
    return out


def conjecture(X, y):
    """CAP-T (LLMSRBENCH_DEV.md): the best UNSNAPPED generalized-monomial fit --
    float exponents kept -- for the continuous-parameter law class (Arrhenius
    x^p*exp(-E/x)). Returns (expr_str, log_resid) or None. NEVER certified; the
    caller labels it a conjecture. No per-target choices: every assignment is
    tried, lowest log-residual wins, one loose global gate."""
    X = np.atleast_2d(np.asarray(X, float))
    y = np.asarray(y, float).ravel()
    m = np.isfinite(y) & np.all(np.isfinite(X), axis=1)
    X, y = X[m], y[m]
    dim = X.shape[1]
    if dim > 3 or len(X) < 10 or np.any(X <= 0) or np.any(y <= 0):
        return None
    ly = np.log(y)
    G = {"id": [np.log(X[:, j]) for j in range(dim)],
         "exp": [X[:, j] if np.all(X[:, j] < 100) else None for j in range(dim)],
         "log": [np.log(np.log(X[:, j])) if np.all(X[:, j] > 1.0 + 1e-9) else None
                 for j in range(dim)],
         "invexp": [1.0 / X[:, j] for j in range(dim)]}
    syms = [sp.Symbol(f"x_{i}") for i in range(dim)]
    best = None
    for assign in product(("id", "exp", "log", "invexp"), repeat=dim):
        cols = [G[a][j] for j, a in enumerate(assign)]
        if any(c is None for c in cols):
            continue
        L = np.column_stack([np.ones(len(X))] + cols)
        c = lstsq(L, ly)
        if c is None:
            continue
        resid = float(np.sqrt(np.mean((L @ c - ly) ** 2)))
        if best is None or resid < best[0]:
            best = (resid, assign, c)
    if best is None or best[0] > _CONJ_TOL:
        return None
    resid, assign, c = best
    expr: sp.Expr = sp.Float(float(np.exp(c[0])))
    for j, (a, p) in enumerate(zip(assign, c[1:])):
        if abs(float(p)) < 1e-9:
            continue
        base = {"id": syms[j], "exp": sp.exp(syms[j]), "log": sp.log(syms[j]),
                "invexp": sp.exp(1 / syms[j])}[a]
        expr *= base ** sp.Float(float(p))
    return str(expr), resid

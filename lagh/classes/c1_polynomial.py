"""C1: complete polynomial monomials to degree 3, plus guarded unary transcendentals.

Degree-3 completeness (x_i**2*x_j, x_i*x_j*x_k) is a lesson: the predecessor could
express x**3 and x*y but not x**2*y -- an arbitrary hole that cost a recovery.
"""

from __future__ import annotations

from itertools import combinations

import numpy as np

from ..base import ALWAYS, Term, finite_guard

TIER = 1


def terms(dim: int, X_fit=None, y_fit=None, X_cert=None) -> list[Term]:
    out = [Term("1", lambda X: np.ones(len(X)), ALWAYS, 1)]
    for j in range(dim):
        for p in (1, 2, 3):
            nm = f"x_{j}" if p == 1 else f"x_{j}**{p}"
            out.append(Term(nm, (lambda j, p: lambda X: X[:, j] ** p)(j, p),
                            ALWAYS, 1 + (p > 1)))
    for j, k in combinations(range(dim), 2):
        out.append(Term(f"x_{j}*x_{k}",
                        (lambda j, k: lambda X: X[:, j] * X[:, k])(j, k), ALWAYS, 3))
    for j in range(dim):
        for k in range(dim):
            if j != k:
                out.append(Term(f"x_{j}**2*x_{k}",
                                (lambda j, k: lambda X: X[:, j]**2 * X[:, k])(j, k),
                                ALWAYS, 4))
    for i, j, k in combinations(range(dim), 3):
        out.append(Term(f"x_{i}*x_{j}*x_{k}",
                        (lambda i, j, k: lambda X:
                         X[:, i] * X[:, j] * X[:, k])(i, j, k), ALWAYS, 4))
    for j in range(dim):
        for name, f in (("sin", np.sin), ("cos", np.cos), ("exp", np.exp),
                        ("log", np.log), ("sqrt", np.sqrt)):
            fn = (lambda j, f: lambda X: f(X[:, j]))(j, f)
            out.append(Term(f"{name}(x_{j})", fn, finite_guard(fn), 2))
    return out

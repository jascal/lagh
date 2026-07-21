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
    # FRACTIONAL-power features (Half-1: rational exponents), guarded x>0. Bounded to
    # half/third integers that cover physics (^1.5, ^2.5, ^0.5, ^-1.5, cube roots).
    # These stay in the EXACT-rational class -- the irrational boundary is not crossed.
    from fractions import Fraction as _F
    FRAC_EXP = [_F(1,2), _F(3,2), _F(5,2), _F(-1,2), _F(-3,2),
                _F(1,3), _F(2,3), _F(4,3)]
    def _pos(j):
        return (lambda X, j=j: bool(np.all(X[:, j] > 1e-9)))
    for j in range(dim):
        for e in FRAC_EXP:
            ef = float(e)
            fn = (lambda X, j=j, ef=ef: X[:, j] ** ef)
            nm = f"x_{j}**({e.numerator}/{e.denominator})"
            out.append(Term(nm, fn, _pos(j), 2 + abs(e.numerator)))
    # input * fractional-power-of-another-input (decay: lambda * t^1.5) -- half only,
    # to bound growth
    for j in range(dim):
        for k in range(dim):
            if j == k:
                continue
            for e in (_F(1,2), _F(3,2), _F(1,3), _F(2,3)):
                ef = float(e)
                fn = (lambda X, j=j, k=k, ef=ef: X[:, j] * X[:, k] ** ef)
                nm = f"x_{j}*x_{k}**({e.numerator}/{e.denominator})"
                out.append(Term(nm, fn, _pos(k), 3))
    # CAP-B: trig-MONOMIAL products sin^a(x_k) cos^b(x_k), standalone and times another
    # input x_j. A general structural class: any law that is a (product of an input and a)
    # trigonometric monomial of another input -- (alpha sin + beta cos)^2 expands to the
    # (2,0)/(0,2)/(1,1) block; tan^2 = (2,-2); cot^2 = (-2,2); sin^2/cos^3 = (2,-3). The
    # ratio blocks are guarded so the denominator stays bounded away from zero on the box.
    # REGISTERED BOUND (NEWTONBENCH_GAP_PLAN.md): emitted only for dim<=2. These features
    # are admissible on any input (cos(mass) is defined), so at dim>=3 they roughly double
    # the term count and blow up the C2 rational search; every trig-shaped law in the gap
    # set is dim-2 (input x trig-of-angle). A dim>=3 trig law would be silently missed --
    # this is a stated cap, not a claim of coverage.
    TRIG_AB = [(2, 0), (0, 2), (1, 1), (2, -2), (-2, 2), (2, -3)] if dim <= 2 else []
    def _cos_ok(k):
        return (lambda X, k=k: bool(np.all(np.abs(np.cos(X[:, k])) > 1e-6)))
    def _sin_ok(k):
        return (lambda X, k=k: bool(np.all(np.abs(np.sin(X[:, k])) > 1e-6)))
    def _trig(a, b):
        def f(X, k):
            v = np.ones(len(X))
            if a: v = v * np.sin(X[:, k]) ** a
            if b: v = v * np.cos(X[:, k]) ** b
            return v
        return f
    for k in range(dim):
        for a, b in TRIG_AB:
            base = _trig(a, b)
            # guard: negative sin/cos power needs that function bounded away from 0
            if b < 0 and a < 0:
                guard = lambda X, ck=_cos_ok(k), sk=_sin_ok(k): ck(X) and sk(X)
            elif b < 0:
                guard = _cos_ok(k)
            elif a < 0:
                guard = _sin_ok(k)
            else:
                guard = ALWAYS
            nm = f"sin(x_{k})**{a}*cos(x_{k})**{b}"
            out.append(Term(nm, (lambda base=base, k=k: lambda X: base(X, k))(),
                            guard, 3 + abs(a) + abs(b)))
            for j in range(dim):
                if j == k:
                    continue
                nmj = f"x_{j}*sin(x_{k})**{a}*cos(x_{k})**{b}"
                out.append(Term(nmj,
                                (lambda base=base, j=j, k=k:
                                 lambda X: X[:, j] * base(X, k))(),
                                guard, 4 + abs(a) + abs(b)))
    return out

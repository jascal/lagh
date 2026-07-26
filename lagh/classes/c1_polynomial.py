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
    # CAP-C: higher cross-degree monomials x_j**a * x_k**b (a+b in {4,5}) and pure
    # 4th/5th powers. REGISTERED BOUND (NEWTONBENCH_GAP_PLAN.md): dim<=3 only -- the
    # whole registered gap set (gravity/coulomb charge-sum expansions) is dim-3, and
    # at higher arity these inflate the C2 numerator pools quadratically.
    if dim <= 3:
        for j in range(dim):
            for p in (4, 5):
                out.append(Term(f"x_{j}**{p}",
                                (lambda j, p: lambda X: X[:, j] ** p)(j, p),
                                ALWAYS, 1 + p))
        for j in range(dim):
            for k in range(dim):
                if j == k:
                    continue
                for a, b in ((3, 1), (4, 1), (3, 2), (2, 2)):
                    if (a, b) == (2, 2) and j > k:
                        continue                       # symmetric pair once
                    out.append(Term(f"x_{j}**{a}*x_{k}**{b}",
                                    (lambda j, k, a, b:
                                     lambda X: X[:, j]**a * X[:, k]**b)(j, k, a, b),
                                    ALWAYS, 2 + a + b))
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
    # CAP-A part 2: denom-5/10/3 fractional powers as ADDITIVE summands (hooke's
    # multi-term sums x^{17/5}+x^{1/2}+x^{-10/3}). Deferred until the
    # exact-coefficient gate existed (CAP-E lesson); now bounded to dim<=2 -- the
    # registered gap cells are dim-1 -- so the global feature inflation risk the
    # deferral named stays contained.
    if dim <= 2:
        for j in range(dim):
            for e in (_F(17, 5), _F(-3, 10), _F(-10, 3)):
                ef = float(e)
                fn = (lambda X, j=j, ef=ef: X[:, j] ** ef)
                out.append(Term(f"x_{j}**({e.numerator}/{e.denominator})", fn,
                                _pos(j), 2 + abs(e.numerator)))
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
    # CAP-P (LLMSRBENCH_DEV.md): damped / saturating library features. Measured
    # need: oscillator-family laws are sparse LINEAR combos over exactly these --
    # exp(-|x|) damping, log(|x|+1) soft saturation, |x|^{1/3} sublinear response
    # -- plus input x feature products. Bounded dim<=3; all guarded finite.
    if dim <= 3:
        _damped = [("exp(-Abs(x_{j}))", lambda X, j: np.exp(-np.abs(X[:, j]))),
                   ("log(Abs(x_{j}) + 1)",
                    lambda X, j: np.log(np.abs(X[:, j]) + 1.0)),
                   ("Abs(x_{j})**(1/3)",
                    lambda X, j: np.abs(X[:, j]) ** (1.0 / 3.0))]
        for j in range(dim):
            for nm_t, f in _damped:
                fn = (lambda X, j=j, f=f: f(X, j))
                out.append(Term(nm_t.format(j=j), fn, ALWAYS, 3))
                # SELF-products included (i == j): x*exp(-|x|) is the canonical
                # damped-oscillator term; skipping it was the measured PO1 gap
                for i in range(dim):
                    fni = (lambda X, i=i, j=j, f=f: X[:, i] * f(X, j))
                    out.append(Term(f"x_{i}*" + nm_t.format(j=j), fni, ALWAYS, 4))
    # CAP-R (LLMSRBENCH_DEV.md): plain ratio monomials x_i/x_j. The affine-
    # denominator rationals (x_0/(x_1*(x_2+1))) need P-terms like x_0/x_1 in the
    # C2 implicit pass; without them four verified benchmark laws were
    # unreachable. Bounded: dim<=5, n(n-1) terms, positive-input guard.
    if dim <= 5:
        for j in range(dim):
            for k in range(dim):
                if j != k:
                    fn = (lambda X, j=j, k=k: X[:, j] / X[:, k])
                    out.append(Term(f"x_{j}/x_{k}", fn, _pos(k), 3))
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
    # Bound lifted 2 -> 3 (LLMSRBENCH_DEV.md): at dim 3 the features are emitted
    # only for ANGLE-PLAUSIBLE columns (range within (0, 2*pi)), so mass/charge
    # columns don't inflate the search -- the sec-form laws x_0/(x_1*cos(x_2))
    # were unreachable under the old blanket dim<=2 bound.
    TRIG_AB = [(2, 0), (0, 2), (1, 1), (2, -2), (-2, 2), (2, -3), (0, -1),
               (-1, 0)] if dim <= 3 else []
    def _angleish(k):
        if X_fit is None or dim <= 2:
            return True
        col = np.asarray(X_fit, float)[:, k]
        return bool(np.all((col > 0) & (col < 2 * np.pi)))
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
        if not _angleish(k):
            continue
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

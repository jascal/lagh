"""C7: the Lévy-exponent grammar (specialized domain tier, docs/TESTBED_LEVY.md).

The Lévy-Khintchine exponent of a SYMMETRIC infinitely-divisible law is one of a small,
specific set of forms (plus independent sums). The general fractional-power library
(C1-C6) is too rich for the narrow achievable CF range and admits impostors
(x^(4/3) approximating a Gaussian quadratic), so this tier RESTRICTS the candidate set
to the Lévy forms -- exactly as C6 restricts to quasi-polynomials on integer lattices.

L(u) = log|phi_hat(u)|. Symmetric building blocks:

    Gaussian            -a * u^2              a > 0
    compound Poisson     b * (cos(w u) - 1)   b > 0, symmetric jumps at +-1/w
    variance-gamma/Gamma -c * log(1 + u^2/d)  c,d > 0
    symmetric stable    -e * |u|^alpha        alpha in (0,2], RATIONAL only

A Lévy exponent is an independent SUM of these (superposition of independent parts).
Recovery fits coefficients by linear least squares over the chosen feature subset,
inner nonlinear params (w in cos, d in log, alpha in |u|^p) by a bounded grid snapped
to small rationals -- the same solve-outside/grid-inside discipline as C4, no free
irrational parameter admitted (the irrational-alpha boundary is preserved).
"""

from __future__ import annotations

from fractions import Fraction
from itertools import combinations

import numpy as np
import sympy as sp

from ..base import Candidate, lstsq, snap_all, snap_small

TIER = 7

# nonlinear inner-parameter grids (snapped to small rationals -> exact, no irrationals)
W_GRID = [Fraction(k, 2) for k in range(1, 13)]                    # cos(w u)
D_GRID = [Fraction(k, 2) for k in range(1, 21)]                    # log(1 + u^2/d)
ALPHA_GRID = [Fraction(p, q) for q in (1, 2, 3) for p in range(1, 2 * q + 1)
              if 0 < Fraction(p, q) <= 2]                          # |u|^alpha, rational


def is_levy_domain(X: np.ndarray, y: np.ndarray) -> bool:
    """1-D nonpositive-valued input (log|phi| <= 0 for a valid CF), the CF signature."""
    X = np.asarray(X, float)
    y = np.asarray(y, float).ravel()
    return bool(X.ndim == 2 and X.shape[1] == 1 and np.all(X[:, 0] > 0)
                and np.median(y) <= 1e-9)


def _features(u: np.ndarray):
    """(name, column, sympy-builder) for every Lévy building block over its grid."""
    x = u[:, 0]
    feats = []
    # Gaussian: -u^2  (sign folded into the fitted coefficient, constrained >=0 later)
    feats.append(("gauss", -x**2, lambda c: -c * sp.Symbol("x_0")**2))
    # compound Poisson: cos(w u) - 1
    for w in W_GRID:
        wf = float(w)
        feats.append((f"cpois_{w}", np.cos(wf * x) - 1.0,
                      (lambda c, w=w: c * (sp.cos(sp.Rational(w.numerator, w.denominator)
                                                  * sp.Symbol("x_0")) - 1))))
    # variance-gamma / Gamma: -log(1 + u^2/d)
    for d in D_GRID:
        df = float(d)
        feats.append((f"vg_{d}", -np.log(1.0 + x**2 / df),
                      (lambda c, d=d: -c * sp.log(1 + sp.Symbol("x_0")**2
                                                  / sp.Rational(d.numerator, d.denominator)))))
    # symmetric stable: -|u|^alpha
    for a in ALPHA_GRID:
        af = float(a)
        feats.append((f"stable_{a}", -np.abs(x)**af,
                      (lambda c, a=a: -c * sp.Abs(sp.Symbol("x_0"))
                       ** sp.Rational(a.numerator, a.denominator))))
    return feats


def candidates(ctx) -> list[Candidate]:
    """Lévy-exponent candidates: single blocks and sums of <=2 blocks, coefficients by
    least squares, prefiltered to near-exact before symbolic construction."""
    u = np.asarray(ctx.X_fit, float)
    y = np.asarray(ctx.y_fit, float).ravel()
    uv = np.asarray(ctx.X_sel, float)
    yv = np.asarray(ctx.y_sel, float).ravel()
    if not is_levy_domain(u, y):
        return []
    yscale = float(np.sqrt(np.mean(yv**2))) + 1e-300
    # statistical prefilter scale: the val residual of a correct fit is the sampling
    # noise ~se, not machine precision. Gate on the se scale (fed via ctx) so
    # statistically-valid Levy fits are not rejected before the real se-based check.
    se_scale = float(getattr(ctx, "se_scale", 0.0))
    gate = max(5.0 * se_scale, 1e-4 * yscale)
    feats = _features(u)
    feats_v = _features(uv)
    cols = np.column_stack([f[1] for f in feats])
    cols_v = np.column_stack([f[1] for f in feats_v])

    out: list[Candidate] = []
    # single blocks and pairs (independent superposition of <=2 Lévy parts)
    idx = list(range(len(feats)))
    supports = [(i,) for i in idx] + list(combinations(idx, 2))
    for sup in supports:
        c = lstsq(cols[:, list(sup)], y)
        if c is None or np.any(c < -1e-9):           # coefficients must be >= 0 (valid Lévy)
            continue
        pred = cols_v[:, list(sup)] @ c
        vr = float(np.sqrt(np.mean((pred - yv) ** 2)))
        if vr > gate:
            continue
        ch = snap_all(c)
        expr = sp.Integer(0)
        for k, i in enumerate(sup):
            expr = expr + feats[i][2](sp.Rational(ch[k].numerator, ch[k].denominator))
        out.append(Candidate(expr=sp.nsimplify(expr, rational=False),
                             complexity=2 * len(sup), channel="c7-levy"))
    return out

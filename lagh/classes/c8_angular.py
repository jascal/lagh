"""CAP-G: the angular / inverse-trig class (NEWTONBENCH_GAP_PLAN.md).

One restricted structure: an ANGLE output whose sine/cosine/tangent is an exact
power-law monomial of the inputs times a trig of one angular input:

    y = deg(T^-1( C * prod_j x_j^{p_j} * G(rad(x_k))^r ))        (degrees convention)
    y =     T^-1( C * prod_j x_j^{p_j} * G(x_k)^r )              (radians convention)

with T, G in {sin, cos, tan} and small rational exponents. Recovery is a log-log
least-squares in the T-transformed target space -- 3 x 3 x dim fits per convention,
no feature-library inflation, every produced law verified by the unchanged checker.

Angle conventions are guarded by RANGE PLAUSIBILITY, never assumed: the degrees
branch fires only when the candidate angular input lies in (0, 360) and the target
in (0, 90); radians only when input in (0, 2*pi) and target in (0, pi/2). A law is
emitted with sympy pi (exact) -- the irrational-constant boundary is not crossed.
"""

from __future__ import annotations

from fractions import Fraction

import numpy as np
import sympy as sp

from ..base import Candidate, lstsq, snap

TIER = 5                 # participates once transform-tier escalation is reached
RAW_TARGET_ONLY = True   # applies its own target transform; never re-applied to a
                         # C5-transformed target

_EXP_CAP = 10            # exponent snap: same denominator family as C3/CAP-A
_FIT_TOL = 1e-8          # log-space residual gate before constructing a candidate


def _deg2rad(v):
    return v * np.pi / 180.0


_CONV = {
    # name: (to_rad(values), plausible_input, plausible_target, sym_angle(x))
    "deg": (_deg2rad,
            lambda x: np.all((x > 0) & (x < 360)),
            lambda y: np.all((y > 0) & (y < 90)),
            lambda s: sp.pi * s / 180),
    "rad": (lambda v: v,
            lambda x: np.all((x > 0) & (x < 2 * np.pi)),
            lambda y: np.all((y > 0) & (y < np.pi / 2)),
            lambda s: s),
}

_TRIG = {"sin": (np.sin, sp.sin, sp.asin), "cos": (np.cos, sp.cos, sp.acos),
         "tan": (np.tan, sp.tan, sp.atan)}


def candidates(ctx) -> list[Candidate]:
    X = np.asarray(ctx.X_fit, float)
    y = np.asarray(ctx.y_fit, float).ravel()
    if len(X) < 8 or np.any(X <= 0):
        return []
    dim = X.shape[1]
    out: list[Candidate] = []
    for conv, (to_rad, x_ok, y_ok, sym_angle) in _CONV.items():
        if not y_ok(y):
            continue
        yr = to_rad(y)
        for T_name, (T_np, _, T_inv_sym) in _TRIG.items():
            ty = T_np(yr)
            if not (np.all(np.isfinite(ty)) and np.all(ty > 0)):
                continue
            lty = np.log(ty)
            for k in range(dim):
                if not x_ok(X[:, k]):
                    continue
                xr = to_rad(X[:, k])
                for G_name, (G_np, G_sym, _) in _TRIG.items():
                    g = G_np(xr)
                    if not (np.all(np.isfinite(g)) and np.all(g > 1e-9)):
                        continue
                    # log T(y) = c0 + sum_{j!=k} p_j log x_j + r log G(x_k)
                    cols = [np.ones(len(X))]
                    cols += [np.log(X[:, j]) for j in range(dim) if j != k]
                    cols.append(np.log(g))
                    L = np.column_stack(cols)
                    c = lstsq(L, lty)
                    if c is None:
                        continue
                    resid = float(np.sqrt(np.mean((L @ c - lty) ** 2)))
                    if resid > _FIT_TOL:
                        continue
                    exps = [Fraction(float(a)).limit_denominator(_EXP_CAP)
                            for a in c[1:]]
                    C = snap(float(np.exp(c[0])))
                    inner = sp.Rational(C.numerator, C.denominator)
                    others = [j for j in range(dim) if j != k]
                    for j, e in zip(others, exps[:-1]):
                        if e != 0:
                            inner *= ctx.syms[j] ** sp.Rational(e.numerator,
                                                                e.denominator)
                    r = exps[-1]
                    if r != 0:
                        inner *= G_sym(sym_angle(ctx.syms[k])) ** \
                            sp.Rational(r.numerator, r.denominator)
                    expr = T_inv_sym(inner)
                    if conv == "deg":
                        expr = 180 / sp.pi * expr
                    out.append(Candidate(
                        expr=expr, complexity=int(sp.count_ops(expr)) + 4,
                        channel=f"c8-angular-{conv}-{T_name}-{G_name}"))
    return out

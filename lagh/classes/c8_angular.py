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
_FIT_TOL = 1e-3   # loose proposal gate: the CHECKER decides; 1e-8 silently killed every
                  # candidate on sigma_rep-quantized data (measured on LLM-SRBench)


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
    if len(X) < 8:
        return []
    dim = X.shape[1]
    out: list[Candidate] = []
    all_pos = bool(np.all(X > 0))
    # CAP-G2 (LLMSRBENCH_DEV.md): identity-target angular monomials
    # y = C * prod x_j^{p_j} * G(x_k)^r -- the sec/csc product family
    # (x_0/(x_1*cos(x_2))) that per-input feature products cannot span. Same
    # log-log machinery with one trig column; CAP-N sign closure applies.
    # Mixed-sign closure: a sign-changing trig FACTOR makes y mixed-sign
    # (-x_0*cos(x_2)/x_1 as theta crosses pi/2). Fit magnitudes -- log|y| against
    # log|G| -- and emit candidates with SIGNED G and both coefficient signs; the
    # checker resolves which sign pattern (if either) actually certifies.
    # Inputs may be MIXED-SIGN too (a pivot-inverted signed quantity): fit
    # log|x_j| magnitudes, then only emit when every mixed-sign column's exponent
    # snaps to an INTEGER -- sign then flows correctly through odd/even powers,
    # and the checker resolves the overall sign pair.
    mnz = (np.abs(y) > 0) & np.all(np.abs(X) > 1e-12, axis=1)
    if mnz.sum() >= 8:
        Xm = X[mnz]
        mixed = [bool(np.any(Xm[:, j] <= 0)) for j in range(dim)]
        lya = np.log(np.abs(y[mnz]))
        for conv, (to_rad, x_ok, _yok, sym_angle) in _CONV.items():
            for k in range(dim):
                if mixed[k] or not x_ok(Xm[:, k]):
                    continue
                xr = to_rad(Xm[:, k])
                for G_name, (G_np, G_sym, _) in _TRIG.items():
                    g = np.abs(G_np(xr))
                    if not (np.all(np.isfinite(g)) and np.all(g > 1e-9)):
                        continue
                    cols = [np.ones(len(Xm))]
                    cols += [np.log(np.abs(Xm[:, j]))
                             for j in range(dim) if j != k]
                    cols.append(np.log(g))
                    c = lstsq(np.column_stack(cols), lya)
                    if c is None:
                        continue
                    resid = float(np.sqrt(np.mean(
                        (np.column_stack(cols) @ c - lya) ** 2)))
                    if resid > _FIT_TOL:
                        continue
                    exps = [Fraction(float(a)).limit_denominator(_EXP_CAP)
                            for a in c[1:]]
                    others_j = [j for j in range(dim) if j != k]
                    if any(mixed[j] and e.denominator != 1
                           for j, e in zip(others_j, exps[:-1])):
                        continue          # fractional power of a signed column
                    C0 = float(np.exp(c[0]))
                    if not np.isfinite(C0):
                        continue
                    C = snap(C0)
                    base = sp.Rational(C.numerator, C.denominator)
                    others = [j for j in range(dim) if j != k]
                    for j, e in zip(others, exps[:-1]):
                        if e != 0:
                            base *= ctx.syms[j] ** sp.Rational(e.numerator,
                                                               e.denominator)
                    r = exps[-1]
                    if r != 0:
                        base *= G_sym(sym_angle(ctx.syms[k])) ** \
                            sp.Rational(r.numerator, r.denominator)
                    for sgn in (1, -1):
                        out.append(Candidate(
                            expr=sgn * base,
                            complexity=int(sp.count_ops(base)) + 3,
                            channel=f"c8-idmono-{conv}-{G_name}"))
    for conv, (to_rad, x_ok, y_ok, sym_angle) in _CONV.items():
        if not all_pos or not y_ok(y):
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
                    C0 = float(np.exp(c[0]))
                    if not np.isfinite(C0):
                        continue
                    C = snap(C0)
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

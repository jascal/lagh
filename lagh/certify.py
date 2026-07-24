"""The fixed honesty core: exhaustive certification, vacuity, coherence, abstention.

Never per-target, never per-class (docs/DISCOVERER.md 3). The error model assembles
the four terms that were each, individually, a measured correction in the predecessor:

    eps_k =   kappa * sigma * prop_k          declared/replication-estimated noise
            + lam_B * se_k                    observation error bar (trajectory mode)
            + MACHINE_REL * |y_k|             float64 representation
            + floor_abs                       instrument/integrator floor

kappa = lam_B = 4, z-score-anchored against REACHABLE truth (protocol-fitted
coefficients -- never an unobtainable exact law). The noise propagation attachment
point (state vs output) is part of the target declaration, never assumed.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from enum import Enum

import numpy as np
import sympy as sp

from .base import eval_expr

KAPPA = 4
LAM_B = 4
MACHINE_REL = 1e3 * np.finfo(float).eps        # ~2.2e-13
TAU = 0.05                                     # material-difference scale
N_PROBE = 300


class Abstain(str, Enum):
    DOMAIN = "domain"          # query outside certified cover
    STRUCTURAL = "structural"  # materially different laws indistinguishable
    NOISE = "noise"            # tolerance swallows the signal (vacuity)
    SURROGATE = "surrogate"    # rests on high-disagreement surrogate region
    NUMERICAL = "numerical"    # law undefined/unstable inside its own domain
    RANGE = "range"            # sampled box carries no signal above the floor
    PARAMETRIC = "parametric"  # exact rational params not pinned within the noise band


@dataclass
class Certificate:
    certified: bool
    nmiss: int
    nuncov: int
    domain_size: int
    state_bounds: list
    law: str
    abstain: str | None = None
    notes: list = field(default_factory=list)

    def one_line(self) -> str:
        if not self.certified:
            return f"ABSTAIN[{self.abstain}] |D|={self.domain_size} nmiss={self.nmiss}"
        return f"CERTIFIED over |D|={self.domain_size}: {self.law}"


def epsilon(y: np.ndarray, *, sigma: float = 0.0, prop: np.ndarray | None = None,
            se: np.ndarray | None = None, floor_abs: float = 1e-12) -> np.ndarray:
    y = np.asarray(y, float).ravel()
    eps = MACHINE_REL * np.abs(y) + floor_abs
    if sigma > 0:
        eps = eps + KAPPA * sigma * (np.abs(y) if prop is None
                                     else np.asarray(prop, float).ravel())
    if se is not None:
        eps = eps + LAM_B * np.asarray(se, float).ravel()
    return eps


def check(expr, syms, X: np.ndarray, y: np.ndarray, eps: np.ndarray) -> dict:
    """Exhaustive: every point checked; wrong-value and no-value kept distinct."""
    y = np.asarray(y, float).ravel()
    pred = eval_expr(expr, syms, X)
    if pred is None:
        return {"certified": False, "nmiss": 0, "nuncov": len(X)}
    uncov = ~np.isfinite(pred)
    miss = (np.abs(pred - y) > eps) & ~uncov
    return {"certified": bool(miss.sum() == 0 and uncov.sum() == 0),
            "nmiss": int(miss.sum()), "nuncov": int(uncov.sum())}


def vacuous(syms, X: np.ndarray, y: np.ndarray, eps: np.ndarray) -> bool:
    """If the ZERO law certifies, eps >= |y| everywhere: certification proves nothing
    at this noise level. Self-detecting by the certificate's own semantics -- no
    threshold exists to tune."""
    return check(sp.Integer(0), syms, X, y, eps)["certified"]


def sample_box(X: np.ndarray, n: int = N_PROBE, seed: int = 0,
               extend: float = 0.0) -> np.ndarray:
    """Uniform samples of the state box, optionally EXTENDED beyond it by `extend`
    fraction of the range on each side. Coherence uses extend>0 so two laws that agree
    only on the thin sampled tube but diverge as functions are flagged as different --
    the thin-domain under-determination fix (MDBench C9 / Levy narrow CF range)."""
    rng = np.random.default_rng(seed)
    lo, hi = X.min(axis=0), X.max(axis=0)
    if extend <= 0:
        return rng.uniform(lo, hi, size=(n, X.shape[1]))
    # extend in LOG space (all lagh sampling is log-uniform/positive), so the extension
    # stays in the valid domain -- an additive extension went negative and broke
    # fractional-power/log/sqrt laws (measured).
    pos = lo > 0
    out = rng.uniform(lo, hi, size=(n, X.shape[1]))
    if np.any(pos):
        llo, lhi = np.log(np.where(pos, lo, 1)), np.log(np.where(pos, hi, 1))
        span = lhi - llo
        ext = rng.uniform(llo - extend * span, lhi + extend * span, size=(n, X.shape[1]))
        out = np.where(pos, np.exp(ext), out)
    return out


def max_divergence(expr_a, expr_b, syms, P: np.ndarray, yscale: float) -> float:
    """Compared on the FINITE OVERLAP of the probe: an inverse-trig law is genuinely
    undefined on part of an extended box (asin beyond its branch) and that is not
    evidence of difference -- but fewer than 8 comparable points is no evidence of
    sameness either, so tiny overlaps stay 'different' (conservative -> abstain)."""
    va, vb = eval_expr(expr_a, syms, P), eval_expr(expr_b, syms, P)
    if va is None or vb is None:
        return float("inf")
    d = np.abs(va - vb)
    d = d[np.isfinite(d)]
    return float(d.max()) / max(yscale, 1e-300) if d.size >= 8 else float("inf")


def pinned(expr, syms, X: np.ndarray, y: np.ndarray, eps: np.ndarray,
           P: np.ndarray, yscale: float, sigma: float, tau: float = TAU) -> bool:
    """Parametric-uncertainty abstention (docs/RNOISE_STUDY.md).

    Coherence catches FUNCTIONAL under-determination (rival materially-different laws).
    Measurement noise instead creates PARAMETRIC mis-determination: one law whose exact
    rational coefficient/exponent was snapped from a noisy fit, with no functional rival
    generated -- so it certifies as exact when a *nearby* rational fits the noisy data
    just as well. This tests the snap for that: for each rational parameter, is there a
    materially-different neighbour-rational, within the declared noise band, that ALSO
    certifies? If so the exact value is not identified at this noise -> not pinned.

    On CLEAN data (sigma<=0) this is a strict no-op (always pinned) -- the deterministic
    behaviour, and the zero-wrong-on-clean record, are preserved by construction. The
    gate only ever REMOVES certifications under declared noise; it can never add one.
    """
    if not np.isfinite(sigma) or sigma <= 0 or expr is None:
        return True
    from fractions import Fraction
    y = np.asarray(y, float).ravel()
    ev_P = eval_expr(expr, syms, P)                   # original law on the probe box
    if ev_P is None:
        return True
    rats = sorted({a for a in expr.atoms(sp.Rational) if a != 0},
                  key=lambda a: (a.q, abs(a.p)))
    for r in rats:
        v = float(r)
        if v == 0:
            continue
        seen = set()
        for kmul in (1.0, 2.0, 4.0):
            for s in (1, -1):
                vp = v * (1 + s * kmul * sigma)
                for md in (max(2, r.q), 2 * r.q, 4 * r.q):
                    rp = Fraction(vp).limit_denominator(int(md))
                    key = (rp.numerator, rp.denominator)
                    if key in seen:
                        continue
                    seen.add(key)
                    if (rp.numerator, rp.denominator) == (r.p, r.q):
                        continue
                    alt = expr.xreplace({r: sp.Rational(rp.numerator, rp.denominator)})
                    if alt == expr:
                        continue
                    av = eval_expr(alt, syms, X)      # neighbour on the cert points
                    if av is None or not np.all(np.isfinite(av)):
                        continue
                    d2 = float(np.dot(av, av))
                    if d2 <= 0:
                        continue
                    # refit ONE overall scale: a different exponent needs its own
                    # coefficient to compete -- without this the neighbour never fits
                    # and every wrong snap looks pinned.
                    alpha = float(np.dot(av, y) / d2)
                    if not np.isfinite(alpha) or alpha == 0:
                        continue
                    ap = eval_expr(alt, syms, P)
                    if ap is None:
                        continue
                    diff = np.abs(alpha * ap - ev_P)
                    diff = diff[np.isfinite(diff)]
                    if not diff.size or float(diff.max()) / max(yscale, 1e-300) <= tau:
                        continue                     # same function at best scale
                    if np.all(np.abs(alpha * av - y) <= eps):
                        return False                 # a different rational also fits
    return True


def float_pinned(expr, syms, X: np.ndarray, y: np.ndarray,
                 eps: np.ndarray, sigma: float = 0.0) -> tuple[bool, sp.Expr]:
    """The exact-coefficient gate (NEWTONBENCH_GAP_PLAN.md CAP-E lesson).

    A winner carrying a raw `Float` coefficient is claiming exactness it cannot show:
    at loose epsilon (extreme output scale -- the reverted BE transform) a coefficient
    3e-5 off the truth certified anyway. The gate demands the data actually PIN each
    float: perturbing it by ±1e-5/±1e-4 relative must BREAK certification. If a
    perturbed value also certifies, the coefficient is not identified at certification
    precision and the exact claim is unfounded -> reject (parametric abstain).

    A pinned float is then snapped to the shortest-decimal exact rational that still
    certifies (6.674e-5 -> 3337/50000000), so certified laws carry exact rationals;
    a pinned float with no certifying decimal snap is kept as-is (identified to
    certification precision, honestly not a short decimal).

    Gates only ever REMOVE or strengthen certifications -- never add one.

    Rational atoms with denominator > 10^6 (the escalating-snap outputs) carry the
    same exposure as raw Floats -- a huge-denominator rational is a float in a
    costume -- so they are gated identically.

    CLEAN-DATA ONLY (the mirror of pinned()'s clean no-op): under declared noise the
    epsilon is sigma-widened, so a +-1e-5 coefficient perturbation can never break
    certification and the rejection rule would kill every float-carrying candidate
    (measured: 60 dB gravity abstained). Under noise, coefficients are inherently
    noise-limited floats (RNOISE_STUDY.md corrected verdict -- the claim is
    STRUCTURAL there) and decimal-snapping would fabricate false exactness, so the
    gate stands down entirely and pinned() owns the parametric question.
    """
    from fractions import Fraction
    if sigma > 0:
        return True, expr
    floats = sorted((a for a in expr.atoms(sp.Float, sp.Rational)
                     if isinstance(a, sp.Float)
                     or (a.q > 10**6 and not isinstance(a, sp.Integer))),
                    key=lambda f: abs(float(f)))
    if not floats:
        return True, expr
    for f in floats:
        v = float(f)
        if v == 0:
            continue
        for rel in (1e-5, 1e-4):
            for s in (1, -1):
                alt = expr.xreplace({f: sp.Float(v * (1 + s * rel))})
                if alt == expr:
                    continue
                if check(alt, syms, X, y, eps)["certified"]:
                    return False, expr
    snapped = expr
    for f in floats:
        v = float(f)
        for digits in range(1, 13):
            fr = Fraction(f"{v:.{digits}e}")
            cand = snapped.xreplace({f: sp.Rational(fr.numerator, fr.denominator)})
            if cand == snapped:
                break
            if check(cand, syms, X, y, eps)["certified"]:
                snapped = cand
                break
    return True, snapped


def coherent(certifying: list, syms, P: np.ndarray, yscale: float,
             tau: float = TAU) -> list:
    """Greedy tau-clustering of certifying laws into functional equivalence classes.
    ALL complexities are rivals -- parsimony may not veto (the equal-complexity rule
    was a measured clean-data artifact)."""
    reps: list[tuple] = []
    for cand in certifying:
        v = eval_expr(cand.expr, syms, P)
        # a certified law may be legitimately undefined on part of the EXTENDED
        # probe (inverse-trig branch bounds); it stays a rival as long as enough of
        # the probe evaluates. Requiring all-finite silently deleted the true law
        # from coherence -> 0 classes -> a false structural abstain (measured on
        # snell asin/acos cells).
        if v is None or np.isfinite(v).sum() < max(8, len(P) // 20):
            continue
        for i, (rexpr, members) in enumerate(reps):
            if max_divergence(cand.expr, rexpr, syms, P, yscale) <= tau:
                members.append(cand)
                break
        else:
            reps.append((cand.expr, [cand]))
    return reps

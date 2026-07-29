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
    alpha_log10: float | None = None   # log10 of the significance bound |H|*q^h
    n_hypotheses: int | None = None    # |H|: candidates actually checked this run
    partial: dict | None = None        # PARTIAL DETERMINATION, first-class: what
                                       # every law consistent with the data at
                                       # this band agrees on, reported ALONGSIDE
                                       # an abstain so the determined part is not
                                       # discarded (see invariant_content)

    def one_line(self) -> str:
        if not self.certified:
            return f"ABSTAIN[{self.abstain}] |D|={self.domain_size} nmiss={self.nmiss}"
        return f"CERTIFIED over |D|={self.domain_size}: {self.law}"


def epsilon(y: np.ndarray, *, sigma: float = 0.0, prop: np.ndarray | None = None,
            se: np.ndarray | None = None, floor_abs: float = 1e-12,
            hard: np.ndarray | None = None) -> np.ndarray:
    """`hard` is a per-point DETERMINISTIC error bound (computed, not declared as
    a scale) and enters with coefficient 1: the coverage factors KAPPA/LAM_B
    exist to turn a stochastic scale into a band, and multiplying a computed
    bound by 4 would loosen epsilon -- the direction that admits impostors (the
    loose-eps lesson). Introduced for weak-form patch bounds (weakform.py)."""
    y = np.asarray(y, float).ravel()
    eps = MACHINE_REL * np.abs(y) + floor_abs
    if sigma > 0:
        eps = eps + KAPPA * sigma * (np.abs(y) if prop is None
                                     else np.asarray(prop, float).ravel())
    if se is not None:
        eps = eps + LAM_B * np.asarray(se, float).ravel()
    if hard is not None:
        eps = eps + np.asarray(hard, float).ravel()
    return eps


def band(eps, expr=None) -> np.ndarray:
    """Resolve an epsilon argument to a per-point array.

    `eps` is either a fixed array — the standard model, where the inputs are
    exact and only y carries error — or a CALLABLE assembled PER CANDIDATE.
    The per-candidate form exists because in the weak form (weakform.py) the
    design-matrix columns are themselves noisy functionals of the same field as
    the target: the residual of y = sum_k c_k X_k carries
    delta_y - sum_k c_k delta_k, whose band depends on the c_k and, because all
    of them are linear functionals of ONE noise realization, is a single norm
    over the combined weight vector rather than a sum of per-column norms.
    A band assembled without the coefficients has to bound them instead, which
    is both looser and an assumption.

    Every gate that re-checks a MODIFIED law (float_pinned's perturbations,
    reduce_to_minimal's drops) therefore gets the band of the law it is actually
    checking, which is the point: the band is a property of the claim.
    """
    e = eps(expr) if callable(eps) else eps
    return np.asarray(e, float).ravel()


def free_dof(expr) -> int:
    """Free numeric parameters of a law: every numeric atom except 0/±1 (the
    conservative-by-design count -- see significance_log10). Non-sympy law
    objects (C6 QuasiPoly) report 0."""
    if expr is None or not hasattr(expr, "atoms"):
        return 0
    return sum(1 for a in expr.atoms(sp.Number) if a not in (0, 1, -1))


def significance_log10(expr, y: np.ndarray, eps: np.ndarray,
                       n_hypotheses: int, value_range: float | None = None
                       ) -> float:
    """log10 of the a-priori false-certification bound alpha <= |H| * prod(q_k)
    (DIRECTION_SIGNIFICANCE.md). Under a uniform-in-range null, a point chance-
    matches with q_k = 2*eps_k / range(y); a candidate with dof free numeric
    parameters is pinned by dof points, so only the h = n - dof held-out points
    count -- conservatively the h points with the LARGEST q_k. This bounds
    CHANCE agreement only: it says the target is not class-null-random, NOT that
    the certified form is unique (structure uniqueness is coherence's job, and
    the approximant-impostor boundary is unaffected)."""
    y = np.asarray(y, float).ravel()
    eps = band(eps, expr)
    # value_range override: a ZERO-VALUED identity's own range is its residual
    # noise, which makes q vacuous; the chance-agreement range is the scale of
    # the CONSTITUENT terms (measured: color identity got alpha > 1)
    R = float(value_range) if value_range else float(np.max(y) - np.min(y))
    if R <= 0 or len(y) == 0:
        return 0.0
    # (C6 QuasiPoly and other non-sympy law objects: dof=0, conservative)
    h = max(0, len(y) - free_dof(expr))
    if h == 0:
        return 0.0
    q = np.minimum(1.0, 2.0 * eps / R)
    q = np.sort(q)[-h:]                       # largest-q points: upper bound
    with np.errstate(divide="ignore"):
        log_q = np.log10(np.maximum(q, 1e-300))
    return float(np.log10(max(n_hypotheses, 1)) + np.sum(log_q))


def check(expr, syms, X: np.ndarray, y: np.ndarray, eps) -> dict:
    """Exhaustive: every point checked; wrong-value and no-value kept distinct.
    `eps` may be a per-candidate band (see `band`)."""
    y = np.asarray(y, float).ravel()
    eps = band(eps, expr)
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
                    if np.all(np.abs(alpha * av - y) <= band(eps, alpha * alt)):
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

    Under declared noise the perturbation is SIGMA-SCALED instead of standing down
    (the stand-down left a hole: with the prefilter also sigma-widened, a
    dyadic-garbage polynomial overfit CERTIFIED on quantized benchmark data --
    the measured PO12 confident-wrong). At +-10/30 sigma: a REAL coefficient
    (determined to ~sigma/sqrt(n)) shifts predictions by ~10 sigma|y| > the 4
    sigma|y| epsilon and breaks certification -> pinned; a junk term contributing
    below epsilon is invisible to its own perturbation -> unpinned -> rejected.
    Conservative direction only: a true-but-marginal term can be rejected
    (abstain), never accepted wrongly. Decimal-snapping stays clean-data-only
    (snapping under noise would fabricate false exactness).
    """
    from fractions import Fraction
    noisy = sigma > 0
    rels = (10.0 * sigma, 30.0 * sigma) if noisy else (1e-5, 1e-4)
    floats = gated_atoms(expr)
    if not floats:
        return True, expr
    for f in floats:
        v = float(f)
        if v == 0:
            continue
        for rel in rels:
            for s in (1, -1):
                alt = expr.xreplace({f: sp.Float(v * (1 + s * rel))})
                if alt == expr:
                    continue
                if check(alt, syms, X, y, eps)["certified"]:
                    return False, expr
    if noisy:
        return True, expr          # pinned; no decimal-snap under noise
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


def minimal(expr, syms, X: np.ndarray, y: np.ndarray, eps: np.ndarray) -> bool:
    """Parsimony enforcement for noise-regime winners: a JUNK term is one you can
    DROP with certification intact (its whole contribution sits inside epsilon);
    every term of a true minimal law is load-bearing. Replaces the sigma-scaled
    perturbation gate, which wrongly rejected real small terms (a 1%-of-y term
    shifts predictions ~0.1 sigma under its own perturbation -- invisible).
    Reject-only; single-term and non-Add winners are trivially minimal."""
    e = sp.expand(expr)
    if not e.is_Add:
        return True
    terms = e.as_ordered_terms()
    if len(terms) < 2:
        return True
    for i in range(len(terms)):
        reduced = sp.Add(*[t for j, t in enumerate(terms) if j != i])
        if check(reduced, syms, X, y, eps)["certified"]:
            return False               # a smaller law also certifies -> not minimal
    return True


def free_atoms(expr) -> list:
    """Every numeric atom a law's parameters live in (0/+-1 excluded, as in
    free_dof). Under declared noise a MODEST-denominator rational is no more
    supported than a float -- the snapper produced it from noisy data -- so the
    interval report covers all of them, not just the ones `gated_atoms` flags."""
    if expr is None or not hasattr(expr, "atoms"):
        return []
    return sorted((a for a in expr.atoms(sp.Number) if a not in (0, 1, -1)),
                  key=lambda f: abs(float(f)))


def gated_atoms(expr) -> list:
    """The numeric atoms an exactness claim rests on: raw Floats, and Rationals
    with denominator > 10^6 (a huge-denominator rational is a float in a
    costume). Shared by `float_pinned` and `parameter_interval` so the gate and
    the interval always talk about the same numbers."""
    if expr is None or not hasattr(expr, "atoms"):
        return []
    return sorted((a for a in expr.atoms(sp.Float, sp.Rational)
                   if isinstance(a, sp.Float)
                   or (a.q > 10**6 and not isinstance(a, sp.Integer))),
                  key=lambda f: abs(float(f)))


def parameter_interval(expr, syms, X: np.ndarray, y: np.ndarray, eps, atom,
                       *, max_rel: float = 0.5, iters: int = 40):
    """The interval of values for `atom` over which the law STILL CERTIFIES,
    every other parameter held — measured by bisection on the certification
    predicate itself, not by a linearization or a covariance.

    This is what a noisy measurement actually determines about a physical
    coefficient. Demanding an exact rational instead is right for a definitional
    identity and wrong for a diffusivity: measured on weak-form heat rows, the
    certifying interval for nu is 0.1 +- 1.5e-8 / 2.1e-6 / 2.6e-4 at field noise
    1e-8 / 1e-6 / 1e-4 -- linear in sigma, centred on the truth, and never an
    exact rational. Returns (lo, hi), or None when the parameter is not bounded
    within `max_rel` of its own value (the data does not determine it at all).
    """
    v0 = float(atom)
    if v0 == 0 or not np.isfinite(v0):
        return None

    def ok(v):
        return check(expr.xreplace({atom: sp.Float(v)}), syms, X, y,
                     eps)["certified"]

    if not ok(v0):
        return None
    out = []
    for sign in (-1.0, 1.0):
        step, far = abs(v0) * 1e-9, None
        while step <= abs(v0) * max_rel:
            if not ok(v0 + sign * step):
                far = v0 + sign * step
                break
            step *= 2.0
        if far is None:
            return None                     # unbounded within max_rel
        near = v0 + sign * step / 2.0
        for _ in range(iters):
            mid = 0.5 * (near + far)
            if ok(mid):
                near = mid
            else:
                far = mid
        out.append(near)
    return (min(out), max(out))


def reduce_to_minimal(expr, syms, X: np.ndarray, y: np.ndarray,
                      eps: np.ndarray):
    """Parsimony REPAIR (supersedes minimal-as-veto): simplify identity disguises
    (sin^2+cos^2 pairs jointly encoding a monomial defeat term-drop testing),
    then iteratively drop any term whose removal keeps full-data certification.
    Returns the reduced expr -- the input itself when already minimal. A true
    law is a fixed point; a vestigial/overfit term is pruned instead of
    poisoning the candidate."""
    e = sp.expand(expr)
    if e.has(sp.sin, sp.cos, sp.tan):
        try:
            e = sp.expand(sp.trigsimp(e))
        except Exception:                                      # noqa: BLE001
            pass
    changed = True
    while changed and e.is_Add:
        changed = False
        for t in e.as_ordered_terms():
            red = sp.expand(e - t)
            if red != 0 and check(red, syms, X, y, eps)["certified"]:
                e = red
                changed = True
                break
    return e


ARBITRATION_MARGIN_LOG10 = 30.0
ARBITRATION_RIVAL_EVIDENCE_MAX = 0.10


def arbitrate_significance(classes: list, y: np.ndarray, eps: np.ndarray,
                           n_hypotheses: int):
    """Significance arbitration at the Müntz boundary (MUNTZ_ARBITRATION.md,
    registered 2026-07-28, AMENDED the same day after validation).

    Rival classes all certify — they agree everywhere the certificate's domain
    claim applies — and their per-representative chance-fit bounds differ by
    q^(dof gap). Two conditions must BOTH hold before a winner is crowned:

    1. MARGIN: exactly one class's alpha beats every rival's by at least
       ARBITRATION_MARGIN_LOG10 orders of magnitude.
    2. RIVAL EVIDENCE: every defeated rival is evidence-starved — it retains
       under ARBITRATION_RIVAL_EVIDENCE_MAX of the certification sample as
       held-out points (h/n, h = n - free_dof), i.e. it is an INTERPOLATION of
       the data, not a law the data constrains. This is the H1a mechanism
       (`_significance_gate`) applied to rivals rather than the winner.

    Condition 2 is what makes condition 1 sound, and it was added because
    condition 1 alone was measured UNSOUND: on `rational-d1` (reach audit) the
    dense channel proposed two fractional-power twins at tier 1 — dof 34 and 45
    of n=80, h/n 0.57 and 0.44 — neither of them the truth (2x+1)/(x+3), which
    lives in a channel the escalation never reached. Their alpha gap was 125
    orders, so the margin alone crowned an approximant that deviates 2e-5 from
    the truth just outside the sampled box: a form-overclaim. In the contests
    arbitration was built for (`mixed-4term-d2`, `sparse6-d2`) the defeated
    rival is a 66/68-term interpolation with h/n = 0.01 and -0.01 — no held-out
    evidence at all — while the winner carries dof 1-2 with h/n ~ 0.98.
    alpha bounds CHANCE agreement, never structure (see significance_log10);
    dismissing a rival the data never constrained is the only structural claim
    it licenses.

    Returns (winning class, note), or None — the structural abstain stands.
    """
    if len(classes) < 2:
        return None
    n = len(np.asarray(y, float).ravel())
    if n == 0:
        return None
    scored = []
    for cl in classes:
        rep = min(cl[1], key=lambda z: z.complexity)
        a = significance_log10(rep.expr, y, eps, n_hypotheses)
        scored.append((a, max(0, n - free_dof(rep.expr)) / n, cl))
    scored.sort(key=lambda t: t[0])
    best, second = scored[0][0], scored[1][0]
    if second - best < ARBITRATION_MARGIN_LOG10:
        return None
    rival_h = [h for _, h, _ in scored[1:]]
    if max(rival_h) >= ARBITRATION_RIVAL_EVIDENCE_MAX:
        return None                      # a rival the data actually constrains
    note = (f"significance arbitration: {len(classes)} rival classes "
            f"certify; winner alpha <= 1e{best:.0f} (held-out fraction "
            f"{scored[0][1]:.2f}) vs nearest rival 1e{second:.0f} (margin "
            f"{second - best:.0f} >= {ARBITRATION_MARGIN_LOG10:.0f} orders); "
            f"every rival is evidence-starved (max held-out fraction "
            f"{max(rival_h):.2f} < {ARBITRATION_RIVAL_EVIDENCE_MAX:g}) -- "
            "the rivals interpolate the certification sample rather than "
            "being constrained by it; the certificate remains a domain claim")
    return scored[0][2], note


def refit_minimal(expr, syms, X: np.ndarray, y: np.ndarray,
                  eps: np.ndarray):
    """Refit-based parsimony collapse (the loose-floor lesson, Gaia C0): a
    multi-term lstsq fit whose OWN basis admits a sub-support already fitting
    within epsilon is carrying unsupported terms -- delicately canceling
    coefficients that term-by-term dropping (reduce_to_minimal) cannot prune.
    Greedy forward selection over the expr's basis functions, REFITTING at each
    step, stopped at the first support whose refit certifies. Returns the
    refit expr (float coefficients -- the winner gate pins them), or the input
    when it is not a plain linear combination or no smaller support certifies."""
    e = sp.expand(expr)
    if not e.is_Add or len(e.args) < 3:
        return expr
    basis = []
    for t in e.as_ordered_terms():
        coeff, g = t.as_independent(*syms)
        if g.free_symbols - set(syms):
            return expr
        if g not in basis:
            basis.append(g)
    if len(basis) < 3:
        return expr
    try:
        fns = [sp.lambdify(syms, g, "numpy") for g in basis]
        cols = []
        for f in fns:
            v = np.broadcast_to(np.asarray(f(*[X[:, j] for j in range(X.shape[1])]),
                                           float), (len(X),)).astype(float)
            cols.append(v)
        A = np.column_stack(cols)
    except Exception:                                          # noqa: BLE001
        return expr
    if not np.all(np.isfinite(A)):
        return expr
    sel: list[int] = []
    for _ in range(len(basis) - 1):          # strictly smaller supports only
        best = None
        for j in range(len(basis)):
            if j in sel:
                continue
            Aj = A[:, sel + [j]]
            w, *_ = np.linalg.lstsq(Aj, y, rcond=None)
            if callable(eps):        # per-candidate band: score against the band
                e_j = band(eps, sum((sp.Float(c) * basis[k]      # of THIS support
                                     for c, k in zip(w, sel + [j])), sp.S.Zero))
            else:
                e_j = eps            # fixed band: no expression to build
            r = float(np.max(np.abs(y - Aj @ w) - e_j))
            if best is None or r < best[0]:
                best = (r, j, w)
        if best is None:
            break
        sel.append(best[1])
        if best[0] <= 0:
            reduced = sum((sp.Float(c) * basis[k]
                           for c, k in zip(best[2], sel)), sp.S.Zero)
            if check(reduced, syms, X, y, eps)["certified"]:
                return reduced
            break
    return expr


def input_constraints(X: np.ndarray, syms, tol_rel: float = 1e-10,
                      max_constraints: int = 2) -> list:
    """Detect MACHINE-EXACT low-degree polynomial constraints the input data
    satisfies (the constrained-input coherence closure, CASE_STUDY_GAIA_P3.md):
    an SVD null direction of the quadratic feature matrix with singular value
    <= tol_rel * sigma_max is a constraint g(x) = 0 holding on the data
    manifold. Machine tolerance ONLY: statistical correlations (e.g. the RRab
    P-phi31 ridge at rho 0.6) must never trigger. Returns sympy exprs with
    rationalized coefficients."""
    n, d = X.shape
    if n < 3 * (1 + d + d * (d + 1) // 2):
        return []
    feats = [np.ones(n)]
    terms = [sp.S.One]
    for i in range(d):
        feats.append(X[:, i])
        terms.append(syms[i])
    for i in range(d):
        for j in range(i, d):
            feats.append(X[:, i] * X[:, j])
            terms.append(syms[i] * syms[j])
    M = np.column_stack(feats)
    scale = np.sqrt(np.mean(M ** 2, axis=0))
    scale[scale == 0] = 1.0
    _, s, vt = np.linalg.svd(M / scale, full_matrices=False)
    out = []
    from fractions import Fraction
    for k in range(len(s) - 1, -1, -1):
        if s[k] > tol_rel * s[0] or len(out) >= max_constraints:
            break
        coef = vt[k] / scale
        coef = coef / np.max(np.abs(coef))
        expr = sp.S.Zero
        for c, t in zip(coef, terms):
            if abs(c) < 1e-9:
                continue
            fr = Fraction(float(c)).limit_denominator(10 ** 6)
            if abs(float(fr) - float(c)) > 1e-9 * max(1.0, abs(float(c))):
                expr = None
                break
            expr = expr + sp.Rational(fr.numerator, fr.denominator) * t
        if expr is not None and expr != 0:
            out.append(sp.expand(expr))
    return out


def reduce_mod_constraints(expr, syms, constraints: list):
    """Canonical form of a POLYNOMIAL expr modulo the constraint ideal
    (sympy.reduced remainder). Non-polynomial exprs are returned unchanged."""
    try:
        if not expr.is_polynomial(*syms):
            return expr
        _, r = sp.reduced(sp.expand(expr), [sp.expand(g) for g in constraints],
                          *syms)
        return sp.expand(r)
    except Exception:                                          # noqa: BLE001
        return expr


def invariant_content(certifying: list, syms, names=None) -> dict:
    """What EVERY law consistent with the data at the declared band agrees on.

    A structural abstain currently discards the part that WAS determined.
    Measured on PDEBench CFD: 662 materially different classes certified, the
    verdict was `ABSTAIN[structural]`, nothing was reported about the
    coefficients -- and the truth check knew the stated law sat at 0.003 of its
    band. Something was strongly determined and the verdict threw it away.

    The claim this returns is deliberately about the VOCABULARY, the DATA and the
    BAND rather than about nature:

        every law in the declared vocabulary consistent with these observations
        at this declared band has coefficient c_k in [lo, hi]

    which is checkable by construction, needs no assumption that the truth is in
    the certifying set, and CANNOT weaken the zero-confident-wrong record: the
    shared content is a weaker statement than any member that already certifies,
    and it is reported ALONGSIDE the abstain, never as a certificate.

    Two things it reports beyond the intervals, and they are what a reader
    actually wants from an under-determined fit: a term that appears in EVERY
    consistent law is REQUIRED, and one that appears in none is EXCLUDED.

    Linear-in-the-columns candidates only -- the case the declared-basis path
    produces. Anything else is skipped rather than guessed at.
    """
    cols = {s: [] for s in syms}
    n = 0
    for cand in certifying:
        try:
            e = sp.expand(sp.sympify(cand.expr))
        except Exception:                                      # noqa: BLE001
            continue
        if e.free_symbols - set(syms):
            continue
        got, rest = {}, e
        ok = True
        for s in syms:
            c = e.coeff(s)
            if c.free_symbols:                     # not linear in the columns
                ok = False
                break
            got[s] = float(c)
            rest = rest - c * s
        if not ok or sp.simplify(rest).free_symbols:
            continue
        n += 1
        for s in syms:
            cols[s].append(got[s])
    if not n:
        return {"n_certifying_read": 0}
    label = {s: (names[i] if names and i < len(names) else str(s))
             for i, s in enumerate(syms)}
    out, required, excluded = {}, [], []
    for s in syms:
        v = np.array(cols[s], float)
        lo, hi = float(v.min()), float(v.max())
        out[label[s]] = {"lo": lo, "hi": hi, "span": hi - lo,
                         "always_present": bool(np.all(v != 0.0)),
                         "never_present": bool(np.all(v == 0.0))}
        if out[label[s]]["always_present"]:
            required.append(label[s])
        if out[label[s]]["never_present"]:
            excluded.append(label[s])
    order = sorted(out, key=lambda k: out[k]["span"])
    return {"n_certifying_read": n, "coefficients": out,
            "required_terms": required, "excluded_terms": excluded,
            "tightest": order[:3], "loosest": order[-3:],
            "claim": ("every law in the declared vocabulary consistent with "
                      "these observations at this declared band has each "
                      "coefficient inside the stated range; terms listed as "
                      "required appear in all of them and terms listed as "
                      "excluded in none. This is a statement about the "
                      "vocabulary, the data and the band -- NOT a certificate")}


def determination(entries: list, *, status: str, note: str = "") -> dict:
    """ONE vocabulary for partial determination, whatever produced it.

    lagh states partial determination in at least five places and, until this
    existed, in five unrelated encodings: interval-parameter certificates put it
    in a note string, state certificates in a per-mode dict, a structural abstain
    in a count, the constrained-input closure in free text, and the two-strata
    campaigns in a convention. A reader could not compare them and the machinery
    could not compose them -- while `Certificate` was `certified: bool` plus an
    abstain enum, a binary with an excuse.

    `entries` is a list of (name, lo, hi) with lo == hi meaning EXACT and either
    bound None meaning unconstrained in that direction. Each becomes:

      exact         -- pinned to a value
      interval      -- bounded, and RESOLVED when the interval excludes zero
                       (the mode is present, not merely bounded)
      unconstrained -- the data does not determine it, and says so

    `status` names what produced the record (`certified`, `structural-abstain`,
    `state`), because a range from a single certified law and a range over a
    certifying SET are different claims and must not be silently merged.
    """
    out, kinds = {}, {"exact": [], "interval": [], "unconstrained": []}
    for name, lo, hi in entries:
        if lo is None or hi is None or not np.isfinite(lo) or not np.isfinite(hi):
            k = "unconstrained"
            rec = {"kind": k, "lo": None, "hi": None, "resolved": False}
        elif lo == hi:
            k = "exact"
            rec = {"kind": k, "lo": lo, "hi": hi, "span": 0.0,
                   "resolved": bool(lo != 0.0)}
        else:
            k = "interval"
            rec = {"kind": k, "lo": lo, "hi": hi, "span": hi - lo,
                   "resolved": bool(not (lo <= 0.0 <= hi))}
        out[str(name)] = rec
        kinds[k].append(str(name))
    return {"status": status, "components": out,
            "exact": kinds["exact"], "interval": kinds["interval"],
            "unconstrained": kinds["unconstrained"],
            "n_resolved": sum(1 for r in out.values() if r["resolved"]),
            "note": note}


def coherent(certifying: list, syms, P: np.ndarray, yscale: float,
             tau: float = TAU, n_evidence: int | None = None) -> list:
    """Greedy tau-clustering of certifying laws into functional equivalence classes.
    ALL complexities are rivals -- parsimony may not veto (the equal-complexity rule
    was a measured clean-data artifact).

    `n_evidence` (the certification sample size) enables an EXACT early exit, not
    an approximation. Arbitration can crown a winner only when every defeated
    rival is evidence-starved (held-out fraction h/n < ARBITRATION_RIVAL_EVIDENCE_MAX
    -- see arbitrate_significance). So the moment TWO classes are found whose
    representatives each retain h/n >= that bar, at most one of them can be the
    winner and the other is a rival that is not evidence-starved: arbitration
    must decline and the verdict is a structural abstain, whatever the remaining
    candidates would have clustered into. Stopping there changes no verdict.

    It matters because this clustering is PAIRWISE -- quadratic in the certifying
    set -- and a loose declared band inflates that set: measured on PDEBench CFD,
    where a 1e-2 declaration produced 662 classes for one equation and 1553 s of
    clustering to reach an abstain the second evidence-bearing class had already
    settled. Cost rising with the messiness of the data is backwards.
    """
    reps: list[tuple] = []
    with_evidence = 0
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
            if n_evidence:
                h = max(0, n_evidence - free_dof(cand.expr)) / n_evidence
                if h >= ARBITRATION_RIVAL_EVIDENCE_MAX:
                    with_evidence += 1
                    if with_evidence >= 2:
                        return reps        # the structural verdict is settled
    return reps

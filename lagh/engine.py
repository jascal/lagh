"""The discovery engine: escalate through the curriculum, certify, decide.

Escalation is MDL-ordered and happens ONLY when the certifying set is empty
(docs/DISCOVERER.md 4). A nonempty-but-incoherent set is a verdict, not a reason
to escalate: the data does not identify a structure at that tier's reach, and
richer tiers only widen the ambiguity.
"""

from __future__ import annotations

from dataclasses import dataclass
from itertools import combinations

import numpy as np
import sympy as sp

from .base import ALWAYS, Candidate, admissible, design_matrix, lstsq
from .base import Term as BaseTerm
from .base import eval_expr, snap_all, to_expr
from .certify import (MACHINE_REL, Abstain, Certificate,
                      arbitrate_significance, check, coherent, epsilon,
                      float_pinned, free_atoms, gated_atoms,
                      input_constraints, invariant_content, minimal,
                      parameter_interval, pinned,
                      reduce_mod_constraints, reduce_to_minimal,
                      refit_minimal, sample_box, significance_log10, vacuous)
from .classes import CURRICULUM, c5_transforms, c6_quasipoly, c7_levy
from .engine_util import stlsq_supports

PREFILTER_REL = 1e-6
MIN_TOTAL_POINTS = 8   # never certify on fewer TOTAL valid points (significance floor)


@dataclass
class Ctx:
    syms: list
    terms: list
    X_fit: np.ndarray
    y_fit: np.ndarray
    X_sel: np.ndarray
    y_sel: np.ndarray
    sigma: float = 0.0    # declared noise: prefilters must not out-tighten epsilon
    deep_supports: bool = True   # CAP-Q sparse-support search (raw target only;
                                 # transform passes set False)
    band_sel: np.ndarray | None = None   # declared per-row band on the SELECT
                                 # split, when the caller has one (weak-form
                                 # rows): the prefilter is stated in relative-y
                                 # units, which such rows do not obey
    linear_basis: bool = False   # the claim is a LINEAR COMBINATION of the
                                 # declared columns (a weak-form PDE always is),
                                 # so proposal is exhaustive over supports and
                                 # the engine builds no products/powers of them
    n_cert: int | None = None    # certification sample size: bounds the largest
                                 # support worth proposing (dof >= n_cert has no
                                 # held-out evidence and is always demoted)


@dataclass
class Result:
    certificate: Certificate
    expr: sp.Expr | None
    tier: int
    n_candidates: int

    @property
    def abstained(self) -> bool:
        return not self.certificate.certified



ALPHA_CERT_MAX_LOG10 = -6.0


def _significance_gate(cert: Certificate) -> Certificate:
    """Significance is part of certification (H1a closure; measured: a 35-term
    interpolation of 12 binned points 'certified' with its own alpha bound at
    1 -- dof >= points, zero held-out evidence). A certificate whose chance-fit
    bound alpha is not small is NO certificate: demote to a NOISE abstain with
    the bound in the note. Threshold registered at 1e-6; every legitimate
    certificate the program has produced sits at 1e-14 or below."""
    if cert.certified and cert.alpha_log10 is not None \
            and cert.alpha_log10 > ALPHA_CERT_MAX_LOG10:
        cert.certified = False
        cert.abstain = Abstain.NOISE.value
        cert.notes.append(
            f"significance gate: alpha_log10 {cert.alpha_log10:.3f} > "
            f"{ALPHA_CERT_MAX_LOG10:g} -- the domain carries too little "
            f"evidence for this form's dof")
    return cert


LINEAR_BASIS_BUDGET = 20000     # exhaustive-support budget (see _basis_supports)


def _basis_supports(T: int, budget: int = LINEAR_BASIS_BUDGET,
                    max_size: int | None = None) -> tuple:
    """Every support up to the largest size the budget affords, EXHAUSTIVELY.

    `max_size` is bounded by the CERTIFICATION SPLIT, and that bound is exact
    rather than a budget heuristic: significance counts only the held-out points
    h = n - dof, so a support with dof >= n_cert has h = 0, alpha_log10 = 0, and
    is demoted by `_significance_gate` no matter how well it fits. Proposing such
    supports cannot produce a certificate -- it can only enlarge the certifying
    set that coherence must then cluster pairwise. Measured on PDEBench CFD: 5
    certification rows against 12 declared columns, 8191 supports enumerated, 662
    materially different classes, 1553 s to reach an abstain.

    For a DECLARED linear basis the candidate space is the library itself, so
    proposal should be complete rather than greedy: measured on the PDE system
    stages, a 5-term true support over 8 declared columns was proposed by
    neither STLSQ, OMP, nor the size-5-from-top-6-singles rule, and the engine
    returned ZERO candidates while the truth sat 4 orders inside its own band.
    Greedy proposal is the right tool for an open-ended library; for a registered
    one it is a hole.

    Returns (supports, max_exhaustive_size) -- the size is reported, never
    silently capped: above it the greedy channels still contribute.
    """
    from math import comb
    top = T if max_size is None else max(1, min(T, int(max_size)))
    sups, count, reached = set(), 0, 0
    for size in range(1, top + 1):
        c = comb(T, size)
        if count + c > budget:
            break
        sups |= set(combinations(range(T), size))
        count += c
        reached = size
    return sups, reached


def _linear_candidates(ctx: Ctx) -> list[Candidate]:
    M_tr = design_matrix(ctx.terms, ctx.X_fit)
    M_va = design_matrix(ctx.terms, ctx.X_sel)
    y_tr = np.asarray(ctx.y_fit, float).ravel()
    y_va = np.asarray(ctx.y_sel, float).ravel()
    yscale = float(np.sqrt(np.mean(y_va**2))) + 1e-300
    sups = stlsq_supports(M_tr, y_tr)
    if ctx.linear_basis:
        # bounded by the certification split, exactly (see _basis_supports)
        cap = None if ctx.n_cert is None else ctx.n_cert - 1
        sups |= _basis_supports(len(ctx.terms), max_size=cap)[0]
    sups |= {c for size in (1, 2) for c in combinations(range(len(ctx.terms)), size)}
    # Targeted size-3..5 supports from the best singles -- the same measured
    # lesson as C2's pairs/triples fix, finally applied to the linear channel:
    # STLSQ cannot isolate sparse multi-term supports in a collinear library, so
    # 3+-term LINEAR laws (the oscillator family) were unreachable by proposal.
    # Bounded: C(10,3)+C(8,4)+C(6,5) = 196 extra lstsq.
    singles = []
    for k in range(len(ctx.terms)):
        c1 = lstsq(M_tr[:, [k]], y_tr)
        if c1 is None:
            continue
        rr = y_tr - M_tr[:, [k]] @ c1
        singles.append((float(rr @ rr), k))
    singles.sort()
    top = [k for _, k in singles[:10]]
    sups |= {tuple(sorted(t)) for t in combinations(top, 3)}
    sups |= {tuple(sorted(t)) for t in combinations(top[:8], 4)}
    sups |= {tuple(sorted(t)) for t in combinations(top[:6], 5)}
    # CAP-T (dim-2 quadratic reach closure, CASE_STUDY_GAIA_P2.md): the COMPLETE
    # low-degree polynomial supports, proposed unconditionally. The Gaia RRab
    # pipeline formula {1, x0, x1, x0*x1, x1**2} was unreachable: size-5
    # supports came only from top-6 singles (fragile under collinearity), OMP's
    # single greedy path missed it, and CAP-Q's exhaustive stops at size 4. A
    # complete quadric/cubic support is one lstsq; refit_minimal / the
    # minimality repair prunes the terms the data does not need.
    name_idx = {t.name: i for i, t in enumerate(ctx.terms)}
    dim_ = ctx.X_fit.shape[1]
    if dim_ <= 4:
        quad = (["1"] + [f"x_{j}" for j in range(dim_)]
                + [f"x_{j}**2" for j in range(dim_)]
                + [f"x_{j}*x_{k}" for j, k in combinations(range(dim_), 2)])
        sup = tuple(sorted(name_idx[n] for n in quad if n in name_idx))
        if len(sup) >= 3:
            sups.add(sup)
    if dim_ <= 2:
        cub = (["1"] + [f"x_{j}" for j in range(dim_)]
               + [f"x_{j}**{p}" for j in range(dim_) for p in (2, 3)]
               + [f"x_{j}*x_{k}" for j, k in combinations(range(dim_), 2)]
               + [f"x_{j}**2*x_{k}" for j in range(dim_)
                  for k in range(dim_) if j != k])
        sup = tuple(sorted(name_idx[n] for n in cub if n in name_idx))
        if len(sup) >= 3:
            sups.add(sup)
    # Orthogonal matching pursuit: greedy forward selection against the residual,
    # every prefix emitted as a support. Top-N single rankings miss collinear
    # members of a sparse sum (measured: a 4-term oscillator support was never
    # proposed); OMP is the standard sparse-recovery tool for exactly this.
    omp: list[int] = []
    resid = y_tr.copy()
    for _ in range(6):
        with np.errstate(all="ignore"):
            norms = np.sqrt(np.einsum("ij,ij->j", M_tr, M_tr))
            scores = np.abs(M_tr.T @ resid) / np.where(norms > 0, norms, np.inf)
        scores[omp] = -1
        k = int(np.argmax(scores))
        if scores[k] <= 0:
            break
        omp.append(k)
        c_o = lstsq(M_tr[:, omp], y_tr)
        if c_o is None:
            break
        resid = y_tr - M_tr[:, omp] @ c_o
        sups.add(tuple(sorted(omp)))
    # CAP-Q (LLMSRBENCH_DEV.md): correlation-pruned exhaustive supports + swap
    # refinement. In a collinear library every greedy/ranked proposal missed the
    # true sparse sum (measured PO33: zero candidates); exhaustive size-3/4 over
    # cluster REPRESENTATIVES is tractable, and swaps recover the case where the
    # representative is the wrong cluster member. Raw-target pass only.
    T = len(ctx.terms)
    if ctx.deep_supports and T >= 30:
        with np.errstate(all="ignore"):
            Mn = M_tr - M_tr.mean(0)
            sd = np.sqrt(np.einsum("ij,ij->j", Mn, Mn))
            Mn = Mn / np.where(sd > 0, sd, np.inf)
            Corr = np.abs(Mn.T @ Mn)
        reps: list[int] = []
        cluster: dict[int, list[int]] = {}
        for _, k in singles:                    # best-single order
            for r in reps:
                if Corr[k, r] > 0.995:
                    cluster.setdefault(r, []).append(k)
                    break
            else:
                reps.append(k)
        reps = reps[:22]
        scored = []
        for size in (3, 4):
            for s in combinations(reps, size):
                c_s = lstsq(M_tr[:, list(s)], y_tr)
                if c_s is None:
                    continue
                rr = y_tr - M_tr[:, list(s)] @ c_s
                scored.append((float(rr @ rr), s))
        scored.sort(key=lambda z: z[0])
        best = [s for _, s in scored[:12]]
        for s in best:
            sups.add(tuple(sorted(s)))
            sl = list(s)
            for i, k in enumerate(sl):
                for alt in cluster.get(k, [])[:6]:
                    sups.add(tuple(sorted(sl[:i] + [alt] + sl[i + 1:])))
    out = []
    for sup in sups:
        cols = list(sup)
        c = lstsq(M_tr[:, cols], y_tr)
        if c is None:
            continue
        pred = M_va[:, cols] @ c
        if not np.all(np.isfinite(pred)):
            continue
        vr = float(np.sqrt(np.mean((pred - y_va) ** 2)))
        # sigma-aware, REPRESENTATION-SCALE ONLY (LLMSRBENCH_DEV.md): a clean-data
        # 1e-6 gate rejected every multi-term law on quantized data (sigma_rep ~
        # 1e-4) before certification. The widening is CAPPED at 3e-4: at
        # statistical noise (>=1e-3) it admitted impostors without adding
        # recoveries (measured RNOISE regression -- the noisy-recoverable set
        # lives in the unprefiltered channels). sigma=0 -> unchanged.
        lim = max(PREFILTER_REL, min(3.0 * ctx.sigma, 3e-4)) * yscale
        if ctx.band_sel is not None:
            # A DECLARED band says what residual the error model itself expects;
            # a relative-y threshold cannot know that (measured: weak-form rows
            # at field-sigma 1e-6 carry ~1e-3 relative row error, so every
            # support -- including the truth -- was prefiltered away and the
            # engine proposed NOTHING). Widening a PROPOSAL filter is sound in
            # the one direction that matters: it only ever admits more
            # candidates for an exact check, never certifies one.
            lim = max(lim, 3.0 * float(np.sqrt(np.mean(
                np.asarray(ctx.band_sel, float) ** 2))))
        if vr > lim:
            continue
        sub = [ctx.terms[i] for i in cols]
        snapped = snap_all(c)
        expr = to_expr(sub, snapped)
        out.append(Candidate(expr=expr, complexity=sum(t.complexity for t in sub),
                             channel="linear", val_residual=vr))
        # ZERO-CROSSING targets: epsilon at the crossings is the ABSOLUTE floor
        # while snap error scales with the TERM (measured: a 1e-13-relative snap
        # missed certification on Lotka-Volterra rates). When snapping loses
        # precision, also emit the raw-float variant -- the clean-data
        # coefficient gate pins and decimal-snaps float winners soundly.
        # RAW-target pass only (ctx.deep_supports False marks the transform
        # pass): under a log-family transform these linear coefficients become
        # EXPONENTS on inversion, and a float variant there certifies x**Float(e)
        # -- breaking the irrational-exponent wedge (caught by the submit test)
        if ctx.deep_supports and \
                any(abs(float(s_) - ci) > 1e-13 * (abs(ci) + 1e-300)
                    for s_, ci in zip(snapped, c)):
            expr_f = sp.Integer(0)
            for t_, ci in zip(sub, c):
                expr_f = expr_f + sp.Float(ci) * t_.sympy()
            out.append(Candidate(expr=expr_f,
                                 complexity=sum(t.complexity for t in sub),
                                 channel="linear", val_residual=vr))
    return out


def _tier_candidates(tier: int, syms, dim, X_fit, y_fit, X_sel, y_sel,
                     X_cert, sigma: float = 0.0,
                     band_sel: np.ndarray | None = None,
                     linear_basis: bool = False) -> list[Candidate]:
    """All candidates available at `tier`, lower tiers included.

    `linear_basis` restricts the space to linear combinations of the DECLARED
    columns: no powers, no products, no transforms. That is not a shortcut, it
    is what the declaration says -- a weak-form PDE law IS a linear combination
    of its library's columns, and letting the general curriculum build
    x_0**3*x_1 out of patch integrals produces laws like u_xx*[1]**(3/2), the
    measured C1b failure. It also shrinks |H|, which tightens alpha honestly.
    """
    if linear_basis:
        base = [BaseTerm("1", lambda X: np.ones(len(X)), ALWAYS, 1)]
        base += [BaseTerm(f"x_{j}", (lambda j: lambda X: X[:, j])(j), ALWAYS, 1)
                 for j in range(dim)]
        ctx = Ctx(syms, admissible(base, X_fit, X_cert), X_fit, y_fit, X_sel,
                  y_sel, sigma, band_sel=band_sel, linear_basis=True,
                  n_cert=len(X_cert))
        return _linear_candidates(ctx)
    active = [(t, m) for t, m in CURRICULUM if t <= tier]
    base_terms = []
    for t, mod in active:
        if hasattr(mod, "terms"):
            base_terms += mod.terms(dim, X_fit, y_fit, X_cert)
    terms = admissible(base_terms, X_fit, X_cert)
    ctx = Ctx(syms, terms, X_fit, y_fit, X_sel, y_sel, sigma, band_sel=band_sel)
    cands = _linear_candidates(ctx)
    for t, mod in active:
        if hasattr(mod, "candidates"):
            cands += mod.candidates(ctx)
    if tier >= 5:
        for tname, ty_fit, inv in c5_transforms.transforms(y_fit):
            try:
                ty_sel = c5_transforms.apply(tname, y_sel)
            except Exception:                                 # noqa: BLE001
                continue
            if not (np.all(np.isfinite(ty_fit)) and np.all(np.isfinite(ty_sel))):
                continue
            # features must be selected against the TRANSFORMED target: exp(u)
            # amplitudes that explain 1/y are invisible against y (measured in the
            # predecessor -- the BE shape fails without this)
            t_terms = []
            for t2, mod2 in active:
                if hasattr(mod2, "terms"):
                    t_terms += mod2.terms(dim, X_fit, ty_fit, X_cert)
            t_terms = admissible(t_terms, X_fit, X_cert)
            tctx = Ctx(syms, t_terms, X_fit, ty_fit, X_sel, ty_sel, sigma,
                       deep_supports=False)
            inner = _linear_candidates(tctx)
            for t2, mod2 in active:
                if hasattr(mod2, "candidates") and \
                        not getattr(mod2, "RAW_TARGET_ONLY", False):
                    inner += mod2.candidates(tctx)
            for c in inner:
                try:
                    expr = inv(c.expr)
                except Exception:                             # noqa: BLE001
                    continue
                if expr.has(sp.zoo, sp.oo, -sp.oo, sp.nan):
                    continue
                cands.append(Candidate(expr=expr, complexity=c.complexity + 1,
                                       channel=f"t-{tname}"))
    return cands


def discover(X_fit, y_fit, X_sel, y_sel, X_cert, y_cert, *,
             sigma: float = 0.0, se_cert=None, floor_abs: float = 1e-12,
             max_tier: int = 7, hard_cert=None, eps_model=None,
             declared_basis: bool = False, band_sel=None,
             linear_basis: bool = False) -> Result:
    """propose -> certify -> vacuity -> coherence -> answer or abstain.

    Splits must be disjoint: fit, select, certify. Certification is exhaustive on
    (X_cert, y_cert) at the assembled epsilon.

    `declared_basis` says the input columns are a REGISTERED term library whose
    members are the claim's vocabulary (the weak-form PDE library), not an
    open-ended basis a dense channel might fit slop with. It only lifts the
    under-noise dense-channel skip below, whose measured justification (PO12/
    PO40: Taylor-slop approximants certify when the true support goes
    unproposed) is an argument about UNREGISTERED bases and does not transfer:
    with a declared library the candidate space IS the library, and the
    corresponding failure mode there is on-shell degeneracy, which the
    multi-solution holdout closes instead. Default False -- no existing caller
    changes behaviour, and switching it on is a declaration the caller makes.

    `linear_basis` goes one step further and says the claim is a LINEAR
    COMBINATION of those columns -- which a weak-form PDE law always is, since
    every candidate term is already a column. The engine then proposes supports
    EXHAUSTIVELY over the library instead of greedily, and builds no powers,
    products or transforms of the columns. Both halves matter and both were
    measured: the greedy channels proposed nothing at all for a 5-term true
    support over 8 declared columns (PDE system stage 2, where the truth sat 4
    orders inside its own band), and the general curriculum's products of patch
    integrals are the C1b `u_xx*[1]**(3/2)` failure. Default False; it is a
    declaration about the claim, not a search-budget knob.
    """
    X_fit = np.asarray(X_fit, float)
    X_cert = np.asarray(X_cert, float)
    y_cert = np.asarray(y_cert, float).ravel()
    # A non-finite observation is "the oracle declined here" (e.g. total internal
    # reflection returns NaN) -- it is NOT evidence and must not poison the fits:
    # one NaN row made every lstsq return None, so whole in-class cells abstained.
    # Dropping it up front is sound: the certification domain is the finite rows'
    # bounds, and check() then never sees a NaN it would vacuously "cover".
    X_sel = np.atleast_2d(np.asarray(X_sel, float))
    y_fit = np.asarray(y_fit, float).ravel()
    y_sel = np.asarray(y_sel, float).ravel()
    mf = np.isfinite(y_fit) & np.all(np.isfinite(X_fit), axis=1)
    ms = np.isfinite(y_sel) & np.all(np.isfinite(X_sel), axis=1)
    mc = np.isfinite(y_cert) & np.all(np.isfinite(X_cert), axis=1)
    X_fit, y_fit, X_sel, y_sel = X_fit[mf], y_fit[mf], X_sel[ms], y_sel[ms]
    X_cert, y_cert = X_cert[mc], y_cert[mc]
    if se_cert is not None:
        se_cert = np.asarray(se_cert, float).ravel()[mc]
    if hard_cert is not None:
        hard_cert = np.asarray(hard_cert, float).ravel()[mc]
    dim = X_fit.shape[1]
    if min(len(y_fit), len(y_sel), len(y_cert)) < 2:
        cert = Certificate(False, 0, 0, int(len(y_cert)), [], "",
                           abstain=Abstain.RANGE.value,
                           notes=["a split is empty after dropping non-finite rows"])
        return Result(cert, None, 0, 0)
    syms = sp.symbols([f"x_{i}" for i in range(dim)])
    if dim == 1:
        syms = [syms] if not isinstance(syms, (list, tuple)) else list(syms)
    syms = list(syms)
    # ERRORS-IN-VARIABLES REGIME: when the INPUTS carry error too (weak-form patch
    # integrals are noisy functionals of the same field as the target), the band
    # depends on the candidate's own coefficients and cannot be assembled once.
    # `eps_model` is then a callable band (certify.band) and every gate that
    # checks a law -- including the perturbed laws float_pinned invents -- gets
    # the band of the law it is actually checking.
    eps = eps_model if eps_model is not None else \
        epsilon(y_cert, sigma=sigma, se=se_cert, floor_abs=floor_abs,
                hard=hard_cert)
    bounds = [(float(X_cert[:, j].min()), float(X_cert[:, j].max()))
              for j in range(dim)]
    # full-data view for the MINIMALITY gate (split-myopia fix, LLMSRBENCH_DEV.md):
    # a real term material only in a thin region must be visible to the
    # droppable-term test, so it runs on fit+sel+cert at the full-data epsilon
    X_all_m = np.vstack([X_fit, X_sel, X_cert])
    y_all_m = np.concatenate([y_fit, y_sel, y_cert])
    # the minimality gate sees fit+sel+cert, for which only the CERT rows carry a
    # computed hard bound; the patch-max stands in there (uniform patch families
    # make it a close stand-in, and the winner must still certify at the exact
    # per-row eps above)
    eps_all = epsilon(y_all_m, sigma=sigma, floor_abs=floor_abs,
                      hard=(None if hard_cert is None
                            else np.full(len(y_all_m), float(np.max(hard_cert)))))
    # LOOSE-FLOOR REGIME (CASE_STUDY_GAIA_C0.md registered issue): when a
    # declared absolute floor dominates the machine term, the per-candidate
    # float_pinned gate INVERTS -- the true law's coefficients are honestly
    # unpinned at loose epsilon (rejected) while delicately-canceling multi-term
    # approximants are hyper-sensitive to perturbation (accepted). Measured on
    # Gaia C0 at floor 2e-4: the 2-term truth gated out, a 14-term whale
    # certified alone. Same lesson as the noise-gate move at PO12: coherence
    # must see the full certifying set; the gate moves to the winner.
    med_y = float(np.median(np.abs(y_cert)))
    floor_dom = sigma <= 0 and floor_abs > max(1e-9, 100 * MACHINE_REL * med_y)

    # minimum-domain guard: certifying on too few TOTAL valid points is not
    # significant (a constant over 4 overflow-artifact points produced the only
    # confident-wrong in the program's history). The empirical floor the
    # significance direction formalizes: |H|*q^h is meaningless when total evidence
    # is tiny. Counts fit+select+cert so a small cert split alone does not trip it.
    n_total = (int(np.isfinite(np.asarray(y_fit, float)).sum())
               + int(np.isfinite(np.asarray(y_sel, float)).sum())
               + int(np.isfinite(y_cert).sum()))
    if n_total < MIN_TOTAL_POINTS:
        cert = Certificate(False, 0, 0, len(X_cert), bounds, "",
                           abstain=Abstain.RANGE.value,
                           notes=[f"only {n_total} total valid points "
                                  f"(< {MIN_TOTAL_POINTS}); too thin to certify"])
        return Result(cert, None, 0, 0)

    # vacuity first: if the zero law certifies, nothing here can be evidence
    if vacuous(syms, X_cert, y_cert, eps):
        cert = Certificate(False, 0, 0, len(X_cert), bounds, "0",
                           abstain=Abstain.NOISE.value,
                           notes=["VACUOUS: eps swallows the signal"])
        return Result(cert, None, 0, 0)

    # coherence probe EXTENDS beyond the cert box (thin-domain under-determination:
    # impostors that agree on the sampled tube but diverge as functions are caught)
    P = sample_box(X_cert, extend=0.5)
    yscale = float(np.sqrt(np.mean(y_cert**2)))

    # C7: the Lévy-exponent grammar. On a CF-domain target (1-D positive input,
    # nonpositive output = log|phi|), the general library under-determines the exponent
    # (fractional-power impostors); this restricted grammar makes the true form unique.
    # Tried BEFORE the general fallthrough but only on its domain.
    if max_tier >= 7 and dim == 1 and not linear_basis \
            and c7_levy.is_levy_domain(X_fit, y_fit):
        class _Ctx:
            pass
        cx = _Ctx(); cx.X_fit, cx.y_fit = X_fit, y_fit
        cx.X_sel, cx.y_sel = np.asarray(X_sel, float), np.asarray(y_sel, float)
        cx.se_scale = float(np.mean(se_cert)) if se_cert is not None else 0.0
        lc = c7_levy.candidates(cx)
        certifying = [c for c in lc
                      if check(c.expr, syms, X_cert, y_cert, eps)["certified"]]
        if certifying:
            classes = coherent(certifying, syms, sample_box(X_cert, extend=0.5),
                               yscale)
            if len(classes) == 1:
                w = min(classes[0][1], key=lambda z: z.complexity)
                cert = Certificate(True, 0, 0, len(X_cert), bounds, str(w.expr),
                                   notes=["Lévy exponent (C7)"],
                                   alpha_log10=significance_log10(
                                       w.expr, y_cert, eps, len(lc)),
                                   n_hypotheses=len(lc))
                cert = _significance_gate(cert)
                return Result(cert, w.expr if cert.certified else None,
                              7, len(lc))
            cert = Certificate(False, 0, 0, len(X_cert), bounds, "",
                               abstain=Abstain.STRUCTURAL.value,
                               notes=[f"{len(classes)} Lévy exponents certify"])
            return Result(cert, None, 7, len(lc))

    # CAP-S: cost-aware cheap pre-pass (LLMSRBENCH_DEV.md, registered). At dim>=3
    # a time-budgeted run can die inside C2's implicit enumeration before the
    # cheap closed-form classes (C3/C9/C8, and C3 under target transforms) ever
    # run -- 12/19 LLM-verified benchmark certificates were plain C3 monomials the
    # loop never reached. The pre-pass proposes only those (~50 log-fits) and
    # applies the FULL verdict machinery (certification, coefficient gate,
    # coherence on the extended probe, pinned); a unique certifying class returns,
    # anything else falls through to the untouched escalation loop.
    if dim >= 3 and max_tier >= 3 and not floor_dom and not linear_basis:
        from .classes import c3_powerlaw, c8_angular, c9_genmonomial
        pctx = Ctx(syms, [], X_fit, y_fit, X_sel, y_sel, sigma)
        pcands = (c3_powerlaw.candidates(pctx) + c9_genmonomial.candidates(pctx)
                  + c8_angular.candidates(pctx))
        for tname, ty_fit, inv in c5_transforms.transforms(y_fit):
            try:
                ty_sel = c5_transforms.apply(tname, y_sel)
            except Exception:                                 # noqa: BLE001
                continue
            if not (np.all(np.isfinite(ty_fit)) and np.all(np.isfinite(ty_sel))):
                continue
            tctx = Ctx(syms, [], X_fit, ty_fit, X_sel, ty_sel, sigma)
            for c in c3_powerlaw.candidates(tctx):
                try:
                    expr = inv(c.expr)
                except Exception:                             # noqa: BLE001
                    continue
                if expr.has(sp.zoo, sp.oo, -sp.oo, sp.nan):
                    continue
                pcands.append(Candidate(expr=expr, complexity=c.complexity + 1,
                                        channel=f"prepass-t-{tname}"))
        pcert = []
        for c in sorted(pcands, key=lambda z: z.complexity):
            if check(c.expr, syms, X_cert, y_cert, eps)["certified"]:
                if sigma <= 0:      # same ordering rule as the main loop: under
                    ok, gated = float_pinned(c.expr, syms, X_cert, y_cert, eps,
                                             sigma)
                    if not ok:      # noise, coherence sees the full set and the
                        continue    # winner is gated below
                    c.expr = gated
                pcert.append(c)
        if pcert:
            classes = coherent(pcert, syms, P, yscale)
            if len(classes) == 1:
                w = min(classes[0][1], key=lambda z: z.complexity)
                if pinned(w.expr, syms, X_cert, y_cert, eps, P, yscale, sigma) \
                        and (sigma <= 0
                             or check(reduce_to_minimal(w.expr, syms, X_all_m,
                                                        y_all_m, eps_all),
                                      syms, X_cert, y_cert, eps)["certified"]):
                    cert = Certificate(True, 0, 0, len(X_cert), bounds,
                                       str(w.expr), notes=["CAP-S cheap pre-pass"],
                                       alpha_log10=significance_log10(
                                           w.expr, y_cert, eps, len(pcands)),
                                       n_hypotheses=len(pcands))
                    cert = _significance_gate(cert)
                    return Result(cert, w.expr if cert.certified else None,
                                  3, len(pcands))
            # ambiguity or unpinned -> the full loop decides (conservative)

    total = 0
    tiers = [1] if linear_basis else [t for t, _ in CURRICULUM if t <= max_tier]
    for tier in tiers:
        cands = _tier_candidates(tier, syms, dim, X_fit, y_fit, X_sel, y_sel,
                                 X_cert, sigma, band_sel=band_sel,
                                 linear_basis=linear_basis)
        total += len(cands)
        certifying = []
        for c in sorted(cands, key=lambda z: z.complexity):
            r = check(c.expr, syms, X_cert, y_cert, eps)
            if r["certified"]:
                # exact-coefficient gate (CAP-E lesson), per CANDIDATE on CLEAN
                # data: dyadic-garbage overfits certify at floor-dominated eps and
                # would poison coherence. Under DECLARED NOISE the gate moves to
                # the WINNER instead (below): gating candidates first removed the
                # true-but-marginal rival that coherence needed to flag admitted
                # impostors, and impostors then certified ALONE (measured: RNOISE
                # 60 dB structural-CW 2 -> 8). Coherence must see the full
                # certifying set under noise; the winner still cannot carry
                # unpinned coefficients.
                if sigma <= 0:
                    if not floor_dom:
                        ok, gated = float_pinned(c.expr, syms, X_cert, y_cert,
                                                 eps, sigma)
                        if not ok:
                            continue
                        c.expr = gated
                    # floor-dominated: keep the candidate ungated so coherence
                    # sees the true rival; the winner is gated below. First
                    # collapse unsupported basis terms (refit parsimony) so a
                    # full-basis fit whose sub-support suffices rejoins the
                    # class of the law it actually is
                    else:
                        c.expr = refit_minimal(c.expr, syms, X_all_m, y_all_m,
                                               eps_all)
                elif (c.channel in ("linear", "c2-pure", "c2-implicit")
                      and not declared_basis):
                    # SIGNIFICANCE BOUNDARY (measured PO12/PO40): at envelope
                    # epsilon on a bounded box, the DENSE channels certify
                    # Taylor-slop approximants of smooth laws whenever the true
                    # support goes unproposed -- an impostor class no per-
                    # candidate gate can close (that needs |H|*q^h accounting,
                    # DIRECTION_SIGNIFICANCE.md). Under declared noise these
                    # channels are EMPIRICAL-only: their fits reach track B via
                    # the fit scout; certification stays with the small-class
                    # closed-form channels.
                    continue
                certifying.append(c)
        if not certifying:
            continue                              # escalate: reach, not ambiguity
        # SCOPED to the declared-linear-basis path on purpose. The early exit
        # is sound for the verdict (see certify.coherent), but arbitration
        # scores a class by its MIN-COMPLEXITY representative while the exit
        # tests the member that opened the class, and complexity is not dof --
        # so on an open-ended library a class could be evidence-bearing when
        # tested and acquire a lower-complexity, higher-dof representative
        # afterwards. That could turn an arbitrated winner into a structural
        # abstain: conservative, never a false certificate, but a REACH
        # regression. The linear-basis path is where the pathology lives (a
        # loose band inflating the certifying set) and where the verdicts were
        # checked against the slow path; everything else keeps the exhaustive
        # clustering and is bit-identical.
        classes = coherent(certifying, syms, P, yscale,
                           n_evidence=len(y_cert) if linear_basis else None)
        constraint_note = None
        constraints = []
        if len(classes) > 1:
            # CONSTRAINED-INPUT COHERENCE (closure of the GAIA_P3 registered
            # issue): when the inputs satisfy machine-exact polynomial
            # constraints, the probe box leaves the data manifold and rivals
            # differing only by multiples of the constraint LOOK materially
            # different while being identical everywhere the certificate's
            # domain claim applies. Detect the constraint; re-run coherence
            # with the DATA as the probe; a single on-manifold class is a
            # verdict, with the winner canonicalized modulo the constraint.
            constraints = input_constraints(X_all_m, syms)
            if constraints:
                mclasses = coherent(certifying, syms, X_all_m, yscale,
                                    n_evidence=(len(y_cert) if linear_basis
                                                else None))
                if len(mclasses) == 1:
                    classes = mclasses
                    constraint_note = (
                        "domain-restricted certificate: inputs satisfy "
                        + "; ".join(f"{sp.sstr(g)} = 0" for g in constraints)
                        + " (machine-exact); rivals coincide on the "
                        "constraint variety and the representative is "
                        "canonicalized modulo it")
        arb_note = interval_note = None
        if len(classes) > 1:
            # MUNTZ_ARBITRATION.md (registered): the program's own alpha
            # currency arbitrates when one class's chance-fit bound beats
            # every rival's by >= 30 orders of magnitude
            arb = arbitrate_significance(classes, y_cert, eps, total)
            if arb is not None:
                classes = [arb[0]]
                arb_note = arb[1]
        if len(classes) == 1:
            winner = min(classes[0][1], key=lambda z: z.complexity)
            if constraint_note:
                red = reduce_mod_constraints(winner.expr, syms, constraints)
                if red != 0 and check(red, syms, X_cert, y_cert,
                                      eps)["certified"]:
                    # the ideal reduction leaves float dust (~1e-14 residue
                    # terms); the minimality repair sweeps them
                    winner.expr = reduce_to_minimal(red, syms, X_all_m,
                                                    y_all_m, eps_all)
            # parametric-uncertainty gate: under declared noise, abstain if the winner's
            # exact rational params are not pinned (a neighbour-rational within the noise
            # band also certifies). No-op on clean data -- preserves the deterministic
            # zero-wrong record; only ever removes a noisy false-exact certificate.
            if not pinned(winner.expr, syms, X_cert, y_cert, eps, P, yscale, sigma):
                cert = Certificate(False, 0, 0, len(X_cert), bounds, str(winner.expr),
                                   abstain=Abstain.PARAMETRIC.value,
                                   notes=["exact rational parameters not pinned within "
                                          f"the noise band (sigma={sigma:g})"])
                return Result(cert, None, tier, total)
            if sigma > 0:
                winner.expr = reduce_to_minimal(winner.expr, syms, X_all_m,
                                                y_all_m, eps_all)
                if declared_basis:
                    # The dense channel is only reachable under noise because a
                    # declared basis says so, and its candidates arrive with the
                    # fitted coefficient SNAPPED (measured: 90360785/903607052
                    # for a true 1/10) -- a huge-denominator rational is a float
                    # in a costume, which is exactly what this gate exists to
                    # catch. Under noise it is sigma-scaled and never snaps.
                    ok, gated = float_pinned(winner.expr, syms, X_cert, y_cert,
                                             eps, sigma)
                    if ok:
                        # Report the interval EVEN WHEN the gate passes: whether
                        # a coefficient counts as "pinned" is decided by
                        # float_pinned's sigma-scaled perturbation (+-10/30 sigma
                        # relative), which is not derived from this data, while
                        # the interval is. Measured: advection at sigma = 1e-3
                        # passed the gate with a value 1.4e-4 off the truth and
                        # said nothing about its own precision.
                        winner.expr = gated
                        told = []
                        for a_ in free_atoms(winner.expr):
                            iv = parameter_interval(winner.expr, syms, X_cert,
                                                    y_cert, eps, a_)
                            if iv is not None:
                                told.append(f"{a_} in [{iv[0]:.10g}, "
                                            f"{iv[1]:.10g}]")
                        if told:
                            interval_note = ("parameter precision at "
                                             f"sigma={sigma:g}: certifies for "
                                             + "; ".join(told))
                    else:
                        # INTERVAL-PARAMETER CERTIFICATE. An exact rational is
                        # the right claim for a definitional identity and the
                        # wrong one for a measured physical coefficient: a
                        # diffusivity has no reason to be rational, so "not
                        # pinned to an exact value" is not the same as "not
                        # determined". Measure what the data DOES determine --
                        # the interval over which the law still certifies -- and
                        # claim that. An unbounded interval still abstains: then
                        # the parameter really is undetermined.
                        ivs, centred = [], winner.expr
                        for a_ in free_atoms(winner.expr):
                            iv = parameter_interval(winner.expr, syms, X_cert,
                                                    y_cert, eps, a_)
                            if iv is None:
                                ivs = None
                                break
                            ivs.append((str(a_), float(iv[0]), float(iv[1])))
                            centred = centred.xreplace(
                                {a_: sp.Float(0.5 * (iv[0] + iv[1]))})
                        if not ivs:
                            cert = Certificate(
                                False, 0, 0, len(X_cert), bounds,
                                str(winner.expr),
                                abstain=Abstain.PARAMETRIC.value,
                                notes=["declared-basis winner gate: parameters "
                                       f"not determined at sigma={sigma:g}"])
                            return Result(cert, None, tier, total)
                        winner.expr = centred
                        interval_note = (
                            "interval-parameter certificate: no exact value is "
                            "claimed; the law certifies for every parameter in "
                            + "; ".join(f"{v} in [{lo:.10g}, {hi:.10g}]"
                                        for v, lo, hi in ivs))
            elif floor_dom:
                winner.expr = reduce_to_minimal(winner.expr, syms, X_all_m,
                                                y_all_m, eps_all)
                ok, gated = float_pinned(winner.expr, syms, X_cert, y_cert,
                                         eps, sigma)
                if not ok:
                    cert = Certificate(False, 0, 0, len(X_cert), bounds,
                                       str(winner.expr),
                                       abstain=Abstain.PARAMETRIC.value,
                                       notes=["loose-floor winner gate: "
                                              "coefficients not pinned at the "
                                              "declared floor"])
                    return Result(cert, None, tier, total)
                winner.expr = gated
            cert = Certificate(True, 0, 0, len(X_cert), bounds, str(winner.expr),
                               alpha_log10=significance_log10(
                                   winner.expr, y_cert, eps, total),
                               n_hypotheses=total)
            if constraint_note:
                cert.notes.append(constraint_note)
            if arb_note:
                cert.notes.append(arb_note)
            if interval_note:
                cert.notes.append(interval_note)
            cert = _significance_gate(cert)
            return Result(cert, winner.expr if cert.certified else None,
                          tier, total)
        cert = Certificate(False, 0, 0, len(X_cert), bounds,
                           str(min(certifying, key=lambda z: z.complexity).expr),
                           abstain=Abstain.STRUCTURAL.value,
                           notes=[f"at least {len(classes)} materially "
                                  f"different classes certify at tier {tier} "
                                  f"(clustering stops once two of them retain "
                                  f"held-out evidence: arbitration cannot crown "
                                  f"a winner past that point)"])
        if linear_basis:
            # the form is under-determined; report what IS determined rather
            # than discarding it with the abstain (gated to the declared-linear
            # path, where a candidate's coefficients are well defined)
            cert.partial = invariant_content(certifying, syms)
            req = cert.partial.get("required_terms") or []
            exc = cert.partial.get("excluded_terms") or []
            if req or exc:
                cert.notes.append(
                    f"partial determination over {cert.partial['n_certifying_read']} "
                    f"consistent laws: required {req or 'none'}, excluded "
                    f"{exc or 'none'}")
        return Result(cert, None, tier, total)

    # C6: escalate to the exact-integer quasi-polynomial tier when the float tiers
    # are exhausted AND the target is an integer lattice. Float tiers structurally
    # cannot certify exact-integer data, so this is the honest terminus, not a
    # fallback -- and it fires only after C1-C5 have genuinely failed (parsimony).
    X_all = np.vstack([X_fit, np.asarray(X_sel, float), X_cert])
    y_all = np.concatenate([np.asarray(y_fit, float).ravel(),
                            np.asarray(y_sel, float).ravel(), y_cert])
    if max_tier >= 6 and dim == 1 and c6_quasipoly.is_integer_lattice(X_all, y_all):
        qr = c6_quasipoly.recover_integer(X_all[:, 0], y_all)
        if qr.certified:
            cert = Certificate(True, 0, 0, qr.domain_size, bounds, str(qr.quasipoly),
                               notes=[qr.note],
                               alpha_log10=significance_log10(
                                   qr.quasipoly, y_all, np.full(len(y_all), 0.5),
                                   max(total, 1)),
                               n_hypotheses=max(total, 1))
            cert = _significance_gate(cert)
            return Result(cert, qr.quasipoly if cert.certified else None,
                          6, total)
        cert = Certificate(False, 0, 0, qr.domain_size, bounds, "",
                           abstain=qr.abstain, notes=[qr.note])
        return Result(cert, None, 6, total)

    cert = Certificate(False, len(X_cert), 0, len(X_cert), bounds, "",
                       abstain=Abstain.STRUCTURAL.value,
                       notes=[f"no law certifies through tier {max_tier}"])
    return Result(cert, None, max_tier, total)

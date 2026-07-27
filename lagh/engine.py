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

from .base import (Candidate, admissible, design_matrix, lstsq, snap_all,
                   to_expr, eval_expr)
from .certify import (Abstain, Certificate, check, coherent, epsilon,
                      float_pinned, minimal, pinned, reduce_to_minimal,
                      sample_box, significance_log10, vacuous)
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


@dataclass
class Result:
    certificate: Certificate
    expr: sp.Expr | None
    tier: int
    n_candidates: int

    @property
    def abstained(self) -> bool:
        return not self.certificate.certified


def _linear_candidates(ctx: Ctx) -> list[Candidate]:
    M_tr = design_matrix(ctx.terms, ctx.X_fit)
    M_va = design_matrix(ctx.terms, ctx.X_sel)
    y_tr = np.asarray(ctx.y_fit, float).ravel()
    y_va = np.asarray(ctx.y_sel, float).ravel()
    yscale = float(np.sqrt(np.mean(y_va**2))) + 1e-300
    sups = stlsq_supports(M_tr, y_tr)
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
        if vr > max(PREFILTER_REL, min(3.0 * ctx.sigma, 3e-4)) * yscale:
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
        if any(abs(float(s_) - ci) > 1e-13 * (abs(ci) + 1e-300)
               for s_, ci in zip(snapped, c)):
            expr_f = sp.Integer(0)
            for t_, ci in zip(sub, c):
                expr_f = expr_f + sp.Float(ci) * t_.sympy()
            out.append(Candidate(expr=expr_f,
                                 complexity=sum(t.complexity for t in sub),
                                 channel="linear", val_residual=vr))
    return out


def _tier_candidates(tier: int, syms, dim, X_fit, y_fit, X_sel, y_sel,
                     X_cert, sigma: float = 0.0) -> list[Candidate]:
    """All candidates available at `tier`, lower tiers included."""
    active = [(t, m) for t, m in CURRICULUM if t <= tier]
    base_terms = []
    for t, mod in active:
        if hasattr(mod, "terms"):
            base_terms += mod.terms(dim, X_fit, y_fit, X_cert)
    terms = admissible(base_terms, X_fit, X_cert)
    ctx = Ctx(syms, terms, X_fit, y_fit, X_sel, y_sel, sigma)
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
             max_tier: int = 7) -> Result:
    """propose -> certify -> vacuity -> coherence -> answer or abstain.

    Splits must be disjoint: fit, select, certify. Certification is exhaustive on
    (X_cert, y_cert) at the assembled epsilon.
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
    eps = epsilon(y_cert, sigma=sigma, se=se_cert, floor_abs=floor_abs)
    bounds = [(float(X_cert[:, j].min()), float(X_cert[:, j].max()))
              for j in range(dim)]
    # full-data view for the MINIMALITY gate (split-myopia fix, LLMSRBENCH_DEV.md):
    # a real term material only in a thin region must be visible to the
    # droppable-term test, so it runs on fit+sel+cert at the full-data epsilon
    X_all_m = np.vstack([X_fit, X_sel, X_cert])
    y_all_m = np.concatenate([y_fit, y_sel, y_cert])
    eps_all = epsilon(y_all_m, sigma=sigma, floor_abs=floor_abs)

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
    if max_tier >= 7 and dim == 1 and c7_levy.is_levy_domain(X_fit, y_fit):
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
                return Result(cert, w.expr, 7, len(lc))
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
    if dim >= 3 and max_tier >= 3:
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
                    return Result(cert, w.expr, 3, len(pcands))
            # ambiguity or unpinned -> the full loop decides (conservative)

    total = 0
    for tier in [t for t, _ in CURRICULUM if t <= max_tier]:
        cands = _tier_candidates(tier, syms, dim, X_fit, y_fit, X_sel, y_sel,
                                 X_cert, sigma)
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
                    ok, gated = float_pinned(c.expr, syms, X_cert, y_cert, eps,
                                             sigma)
                    if not ok:
                        continue
                    c.expr = gated
                elif c.channel in ("linear", "c2-pure", "c2-implicit"):
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
        classes = coherent(certifying, syms, P, yscale)
        if len(classes) == 1:
            winner = min(classes[0][1], key=lambda z: z.complexity)
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
            cert = Certificate(True, 0, 0, len(X_cert), bounds, str(winner.expr),
                               alpha_log10=significance_log10(
                                   winner.expr, y_cert, eps, total),
                               n_hypotheses=total)
            return Result(cert, winner.expr, tier, total)
        cert = Certificate(False, 0, 0, len(X_cert), bounds,
                           str(min(certifying, key=lambda z: z.complexity).expr),
                           abstain=Abstain.STRUCTURAL.value,
                           notes=[f"{len(classes)} materially different classes "
                                  f"certify at tier {tier}"])
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
            return Result(cert, qr.quasipoly, 6, total)
        cert = Certificate(False, 0, 0, qr.domain_size, bounds, "",
                           abstain=qr.abstain, notes=[qr.note])
        return Result(cert, None, 6, total)

    cert = Certificate(False, len(X_cert), 0, len(X_cert), bounds, "",
                       abstain=Abstain.STRUCTURAL.value,
                       notes=[f"no law certifies through tier {max_tier}"])
    return Result(cert, None, max_tier, total)

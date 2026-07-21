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
                      sample_box, vacuous)
from .classes import CURRICULUM, c5_transforms, c6_quasipoly
from .engine_util import stlsq_supports

PREFILTER_REL = 1e-6


@dataclass
class Ctx:
    syms: list
    terms: list
    X_fit: np.ndarray
    y_fit: np.ndarray
    X_sel: np.ndarray
    y_sel: np.ndarray


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
        if vr > PREFILTER_REL * yscale:
            continue
        sub = [ctx.terms[i] for i in cols]
        expr = to_expr(sub, snap_all(c))
        out.append(Candidate(expr=expr, complexity=sum(t.complexity for t in sub),
                             channel="linear", val_residual=vr))
    return out


def _tier_candidates(tier: int, syms, dim, X_fit, y_fit, X_sel, y_sel,
                     X_cert) -> list[Candidate]:
    """All candidates available at `tier`, lower tiers included."""
    active = [(t, m) for t, m in CURRICULUM if t <= tier]
    base_terms = []
    for t, mod in active:
        if hasattr(mod, "terms"):
            base_terms += mod.terms(dim, X_fit, y_fit, X_cert)
    terms = admissible(base_terms, X_fit, X_cert)
    ctx = Ctx(syms, terms, X_fit, y_fit, X_sel, y_sel)
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
            tctx = Ctx(syms, t_terms, X_fit, ty_fit, X_sel, ty_sel)
            inner = _linear_candidates(tctx)
            for t2, mod2 in active:
                if hasattr(mod2, "candidates"):
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
             max_tier: int = 6) -> Result:
    """propose -> certify -> vacuity -> coherence -> answer or abstain.

    Splits must be disjoint: fit, select, certify. Certification is exhaustive on
    (X_cert, y_cert) at the assembled epsilon.
    """
    X_fit = np.asarray(X_fit, float)
    X_cert = np.asarray(X_cert, float)
    y_cert = np.asarray(y_cert, float).ravel()
    dim = X_fit.shape[1]
    syms = sp.symbols([f"x_{i}" for i in range(dim)])
    if dim == 1:
        syms = [syms] if not isinstance(syms, (list, tuple)) else list(syms)
    syms = list(syms)
    eps = epsilon(y_cert, sigma=sigma, se=se_cert, floor_abs=floor_abs)
    bounds = [(float(X_cert[:, j].min()), float(X_cert[:, j].max()))
              for j in range(dim)]

    # vacuity first: if the zero law certifies, nothing here can be evidence
    if vacuous(syms, X_cert, y_cert, eps):
        cert = Certificate(False, 0, 0, len(X_cert), bounds, "0",
                           abstain=Abstain.NOISE.value,
                           notes=["VACUOUS: eps swallows the signal"])
        return Result(cert, None, 0, 0)

    P = sample_box(X_cert)
    yscale = float(np.sqrt(np.mean(y_cert**2)))
    total = 0
    for tier in [t for t, _ in CURRICULUM if t <= max_tier]:
        cands = _tier_candidates(tier, syms, dim, X_fit, y_fit, X_sel, y_sel, X_cert)
        total += len(cands)
        certifying = []
        for c in sorted(cands, key=lambda z: z.complexity):
            r = check(c.expr, syms, X_cert, y_cert, eps)
            if r["certified"]:
                certifying.append(c)
        if not certifying:
            continue                              # escalate: reach, not ambiguity
        classes = coherent(certifying, syms, P, yscale)
        if len(classes) == 1:
            winner = min(classes[0][1], key=lambda z: z.complexity)
            cert = Certificate(True, 0, 0, len(X_cert), bounds, str(winner.expr))
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
                               notes=[qr.note])
            return Result(cert, qr.quasipoly, 6, total)
        cert = Certificate(False, 0, 0, qr.domain_size, bounds, "",
                           abstain=qr.abstain, notes=[qr.note])
        return Result(cert, None, 6, total)

    cert = Certificate(False, len(X_cert), 0, len(X_cert), bounds, "",
                       abstain=Abstain.STRUCTURAL.value,
                       notes=[f"no law certifies through tier {max_tier}"])
    return Result(cert, None, max_tier, total)

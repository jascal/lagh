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
    va, vb = eval_expr(expr_a, syms, P), eval_expr(expr_b, syms, P)
    if va is None or vb is None:
        return float("inf")
    d = np.abs(va - vb)
    d = d[np.isfinite(d)]
    return float(d.max()) / max(yscale, 1e-300) if d.size else float("inf")


def coherent(certifying: list, syms, P: np.ndarray, yscale: float,
             tau: float = TAU) -> list:
    """Greedy tau-clustering of certifying laws into functional equivalence classes.
    ALL complexities are rivals -- parsimony may not veto (the equal-complexity rule
    was a measured clean-data artifact)."""
    reps: list[tuple] = []
    for cand in certifying:
        v = eval_expr(cand.expr, syms, P)
        if v is None or not np.all(np.isfinite(v)):
            continue
        for i, (rexpr, members) in enumerate(reps):
            if max_divergence(cand.expr, rexpr, syms, P, yscale) <= tau:
                members.append(cand)
                break
        else:
            reps.append((cand.expr, [cand]))
    return reps

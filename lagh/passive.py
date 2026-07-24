"""Passive-data mode: discovery on a fixed dataset, no oracle (docs/DIRECTION_PASSIVE.md).

The regime every fixed-dataset SR benchmark scores: the data arrives sampled where the
benchmark chose, and no query can be made. This wrapper supplies the two passive
substitutes for the active loop's machinery -- K deterministic re-splits of the same
data, and a full-data exhaustive gate that keeps the re-splits sound -- around a
byte-identical `discover()`.

Soundness argument for the re-splits: K splits multiply candidate exposure ~K x, but a
law only certifies if it ALSO passes `check()` on EVERY point in the dataset at the
assembled epsilon. Certification is never granted by a re-split, only re-attempted;
the full-data gate can only remove certifications, never add one.

Sigma is declared-only here: no oracle means no replicates. A benchmark that states
its noise level supplies it; clean benchmarks pass 0.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
import sympy as sp

from .certify import Abstain, Certificate, check, epsilon
from .engine import Result, discover


@dataclass
class PassiveResult:
    result: Result
    resplits_tried: int
    full_check_passed: bool | None   # None when nothing ever certified on a split
    abstain_reasons: list = field(default_factory=list)   # one per non-certifying split

    @property
    def certified(self) -> bool:
        return bool(self.result.certificate.certified and self.full_check_passed)


def discover_passive(X, y, *, sigma: float = 0.0, floor_abs: float = 1e-12,
                     max_tier: int = 7, n_resplits: int = 3,
                     seed: int = 0) -> PassiveResult:
    """X (n,d), y (n,): the dataset as handed out. Returns the first split whose
    certification also survives the full-data gate, else the last abstain."""
    X = np.atleast_2d(np.asarray(X, float))
    y = np.asarray(y, float).ravel()
    m = np.isfinite(y) & np.all(np.isfinite(X), axis=1)
    X, y = X[m], y[m]
    if len(X) < 10:
        cert = Certificate(False, 0, 0, int(len(X)), [], "",
                           abstain=Abstain.RANGE.value,
                           notes=[f"only {len(X)} finite points in the dataset; "
                                  "too thin to attempt discovery"])
        return PassiveResult(Result(cert, None, 0, 0), 0, None, [])
    dim = X.shape[1]
    syms = [sp.Symbol(f"x_{i}") for i in range(dim)]
    eps_full = epsilon(y, sigma=sigma, floor_abs=floor_abs)

    last: Result | None = None
    reasons: list = []
    full_ok: bool | None = None
    for k in range(n_resplits):
        idx = np.random.default_rng(seed + k).permutation(len(X))
        a, b = int(0.6 * len(X)), int(0.8 * len(X))
        r = discover(X[idx[:a]], y[idx[:a]], X[idx[a:b]], y[idx[a:b]],
                     X[idx[b:]], y[idx[b:]], sigma=sigma, floor_abs=floor_abs,
                     max_tier=max_tier)
        last = r
        if not r.certificate.certified:
            reasons.append(r.certificate.abstain)
            continue
        if check(r.expr, syms, X, y, eps_full)["certified"]:
            return PassiveResult(r, k + 1, True, reasons)
        # certified on the split but not on all points: a split artifact the gate
        # exists to catch -- demote and keep trying
        full_ok = False
        reasons.append("full-data-gate")
    if last is not None and last.certificate.certified and full_ok is False:
        # never return a certification that failed the gate
        cert = last.certificate
        cert.certified = False
        cert.abstain = cert.abstain or "structural"
        cert.notes.append("passive: certified on a split but failed the "
                          "full-data exhaustive gate")
        last = Result(cert, None, last.tier, last.n_candidates)
    return PassiveResult(last, n_resplits, full_ok, reasons)

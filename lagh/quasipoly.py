"""Exact quasi-polynomial recovery -- the honesty core in pure integer arithmetic.

Reuses the SEMANTICS of lagh's core (exhaustive exact check, parsimony, coherence,
first-class abstention) instantiated with Fraction arithmetic: certification is exact
integer equality, so the epsilon / noise / floor machinery does not exist here.

A quasi-polynomial: L(t) = sum_i c_i(t) t^i, each c_i periodic with some period p.
Represented as one degree-d polynomial per residue class t mod p, exact-Lagrange-fit.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from fractions import Fraction

from .certify import Abstain


@dataclass
class QuasiPoly:
    period: int
    degree: int
    # per-residue Lagrange nodes: class r -> list of (t, L) determining its polynomial
    nodes: dict = field(default_factory=dict)

    def __call__(self, t: int) -> Fraction:
        pts = self.nodes[t % self.period]
        acc = Fraction(0)
        for i, (ti, Li) in enumerate(pts):
            term = Fraction(Li)
            for j, (tj, _) in enumerate(pts):
                if i != j:
                    term *= Fraction(t - tj, ti - tj)
            acc += term
        return acc

    def __str__(self) -> str:
        return f"quasipoly(period={self.period}, degree={self.degree})"


@dataclass
class QPResult:
    certified: bool
    quasipoly: QuasiPoly | None
    domain_size: int
    abstain: str | None = None
    note: str = ""


def _fit_and_check(ts, Ls, p: int, d: int) -> QuasiPoly | None:
    """Per residue class mod p: use the first d+1 points as exact interpolation nodes,
    then certify on ALL REMAINING points of that class (integer equality). Each class
    self-splits, so no global split can starve it -- a class needs >= d+2 points
    (d+1 to fit, >=1 held out). Certifies iff every held-out point in every class
    matches exactly."""
    by_class: dict[int, list] = {r: [] for r in range(p)}
    for t, L in zip(ts, Ls):
        by_class[t % p].append((t, L))
    nodes = {}
    held = 0
    for r in range(p):
        cls = sorted(by_class[r])
        if len(cls) < d + 2:                 # need d+1 to fit AND >=1 to certify
            return None
        nodes[r] = cls[: d + 1]
    qp = QuasiPoly(p, d, nodes)
    for r in range(p):
        for t, L in sorted(by_class[r])[d + 1:]:
            held += 1
            if qp(t) != Fraction(L):
                return None
    return qp if held > 0 else None


def recover(ts, Ls, *, period_max: int = 12, degree_max: int = 4) -> QPResult:
    """Recover the quasi-polynomial or abstain. ts: sorted distinct positive ints.

    The fit/certify split is RANDOM (seeded), never strided: a stride-k holdout
    aligns with period k and starves both of that period's residue classes on one
    side -- measured, it caused period-2/6 laws to over-abstain. A random split has
    no periodic structure to collide with, so every residue class lands on both
    sides in expectation.
    """
    ts = list(ts)
    Ls = list(Ls)
    for p in range(1, period_max + 1):
        # parsimony: smallest period that certifies at ANY degree is the answer
        for d in range(degree_max + 1):
            qp = _fit_and_check(ts, Ls, p, d)
            if qp is not None:
                held = len(ts) - p * (d + 1)
                return QPResult(True, qp, held,
                                note=f"period={p} degree={d}")

    return QPResult(False, None, 0, abstain=Abstain.RANGE.value,
                    note=f"no quasi-polynomial with period<= {period_max}, "
                         f"degree<= {degree_max} certifies within the "
                         f"{len(ts)}-value budget (underdetermined)")

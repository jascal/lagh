"""Ehrhart oracle: exact lattice-point counts of dilated rational polytopes.

L_P(t) = #(t*P intersect Z^d), computed by exact enumeration over the bounding box.
No floats: the inequality test is exact Fraction arithmetic. The polytope is the
scaled standard simplex  P = { x >= 0 : sum_i x_i / a_i <= 1 },  a_i = p_i/q_i rational,
so  t*P = { x >= 0 : sum_i x_i q_i / (t p_i) <= 1 }  and x_i ranges 0..floor(t a_i).
"""

from __future__ import annotations

from fractions import Fraction
from itertools import product

import numpy as np


def random_simplex(dim: int, rng: np.random.Generator,
                   denoms=(2, 3)) -> list[Fraction]:
    """Scaled standard simplex, given by its axis intercepts a_i (rational)."""
    a = []
    for _ in range(dim):
        q = int(rng.choice(denoms))
        p = int(rng.integers(q + 1, 4 * q + 1))       # a_i = p/q > 1, not integer-ish
        a.append(Fraction(p, q))
    return a


def _count_int(weights: list[int], ubs: list[int], B: int) -> int:
    """Exact #{x >= 0 : x_i <= ub_i, sum w_i x_i <= B}. Last axis closed-formed,
    the rest recursed. Pure integers -- no float, no Fraction in the hot loop."""
    if not weights:
        return 1 if B >= 0 else 0
    if len(weights) == 1:
        w, ub = weights[0], ubs[0]
        return 0 if B < 0 else min(ub, B // w) + 1
    w, ub = weights[0], ubs[0]
    if B < 0:
        return 0
    xmax = min(ub, B // w)
    return sum(_count_int(weights[1:], ubs[1:], B - x * w) for x in range(xmax + 1))


def lattice_count(intercepts: list[Fraction], t: int) -> int:
    """Exact #(t*P intersect Z^d) for the scaled simplex with these intercepts.

    Condition sum_i x_i/(t a_i) <= 1 with a_i = p_i/q_i cleared to integers:
    let L = lcm(p_i); w_i = q_i * L/p_i; B = t*L; then sum_i w_i x_i <= B, with
    0 <= x_i <= floor(t*p_i/q_i). Exact and denominator-free."""
    from math import lcm
    ps = [ai.numerator for ai in intercepts]
    qs = [ai.denominator for ai in intercepts]
    L = lcm(*ps) if len(ps) > 1 else ps[0]
    weights = [q * (L // pp) for q, pp in zip(qs, ps)]
    ubs = [(t * pp) // q for pp, q in zip(ps, qs)]
    return _count_int(weights, ubs, t * L)


def make_oracle(intercepts: list[Fraction]):
    def oracle(ts) -> list[int]:
        return [lattice_count(intercepts, int(t)) for t in np.atleast_1d(ts)]
    return oracle

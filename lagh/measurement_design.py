"""Oracle-free design from structural refusal content.

Scores a finite, caller-declared probe set, not a continuous optimum. Candidate
values are predictions, never acquired evidence. No certification is performed.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import sympy as sp

from .base import eval_expr


@dataclass(frozen=True)
class MeasurementChoice:
    point: tuple[float, ...] | None
    predictions: tuple[float, ...]
    disagreement: float
    cost: float | None
    utility: float
    reason: str


def choose_measurement(rivals, probes, *, bounds, costs=1.0) -> MeasurementChoice:
    """Maximize (max prediction - min prediction) / cost on finite probes.

    Bounds and costs must be declared independently of future responses. Ties
    take the first probe in the declared order. Require ALL rivals to evaluate
    finitely: undefined values do not count as observable disagreement. This
    spread heuristic need not separate every pair in a multi-rival set.
    """
    X = np.asarray(probes, dtype=float)
    box = np.asarray(bounds, dtype=float)
    if (X.ndim != 2 or not X.shape[1] or box.shape != (2, X.shape[1])
            or not np.all(np.isfinite(X)) or not np.all(np.isfinite(box))
            or np.any(box[0] > box[1])
            or np.any(X < box[0]) or np.any(X > box[1])):
        raise ValueError("finite probes must lie inside declared (2, d) bounds")
    c = np.asarray(costs, dtype=float)
    if c.ndim == 0:
        c = np.full(len(X), float(c))
    if c.shape != (len(X),) or np.any(~np.isfinite(c)) or np.any(c <= 0):
        raise ValueError("one finite positive declared cost is required per probe")
    rivals = tuple(rivals)
    empty = MeasurementChoice(None, (), 0.0, None, 0.0,
                              "no finite informative query among declared probes")
    if len(rivals) < 2 or not len(X):
        return empty
    syms = [sp.Symbol(f"x_{i}") for i in range(X.shape[1])]
    values = []
    for rival in rivals:
        v = eval_expr(rival, syms, X)
        if v is None or np.asarray(v).shape != (len(X),):
            return empty
        values.append(v)
    V = np.asarray(values, dtype=float)
    finite = np.all(np.isfinite(V), axis=0)
    with np.errstate(all="ignore"):
        spread = np.max(V, axis=0) - np.min(V, axis=0)
        utility = spread / c
    valid = finite & np.isfinite(utility) & (utility > 0)
    if not np.any(valid):
        return empty
    i = int(np.argmax(np.where(valid, utility, -np.inf)))
    return MeasurementChoice(tuple(map(float, X[i])), tuple(map(float, V[:, i])),
                             float(spread[i]), float(c[i]), float(utility[i]),
                             "maximum finite rival spread per cost on declared probes")

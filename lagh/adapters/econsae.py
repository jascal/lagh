"""Oracle adapter for econ-sae (docs/TESTBED_ECONSAE.md).

econ-sae is a published-package dependency, queried as a black box and NEVER modified.
This wraps `Economy.small().step(shock)` into `oracle(shock_matrix) -> aggregate`,
running burn+avg per row to read a steady-state aggregate. The target law bodies (the
emergent aggregates) are unread by construction -- they do not exist as source.

Shock-key mapping (registration amendment 2026-07-21, appended not edited): the
registration named targets conceptually; interface inspection pins them to the real
`_apply_shock` schema:
    E1 fiscal        shock={'transfer_per_hh': g}          -> GDP
    E2 productivity  shock={'productivity_mult': a}        -> price_level
    E3 credit        shock={'interest_rate': s}            -> debt_outstanding
    E4 distribution  shock={'transfer_per_hh','productivity_mult'} -> Gini(HH money)
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

ECONSAE = Path("/home/allans/code/econ-sae")
if str(ECONSAE) not in sys.path:
    sys.path.insert(0, str(ECONSAE))

N_BURN = 200
N_AVG = 50

TARGETS = {
    "E1": {"inputs": ["transfer_per_hh"], "box": ([0.5], [2.0]), "aggregate": "GDP"},
    "E2": {"inputs": ["productivity_mult"], "box": ([0.8], [1.25]),
           "aggregate": "price_level"},
    "E3": {"inputs": ["interest_rate"], "box": ([0.01], [0.1]),
           "aggregate": "debt_outstanding"},
    "E4": {"inputs": ["transfer_per_hh", "productivity_mult"],
           "box": ([0.5, 0.8], [2.0, 1.25]), "aggregate": "gini_money"},
}


def _gini(x: np.ndarray) -> float:
    x = np.sort(np.abs(np.asarray(x, float)))
    n = len(x)
    if n == 0 or x.sum() == 0:
        return 0.0
    return float((2 * np.arange(1, n + 1) - n - 1).dot(x) / (n * x.sum()))


def _steady(inputs, row, aggregate, *, n_burn=N_BURN, n_avg=N_AVG, seed=0) -> float:
    from econsae.simulator.core import Economy
    from econsae.embeddings import COORD_IDX
    e = Economy.small(seed=seed)
    shock = {k: float(v) for k, v in zip(inputs, row)}
    vals = []
    for t in range(n_burn + n_avg):
        state, _, macros = e.step(shock=shock)
        if t < n_burn:
            continue
        if aggregate == "gini_money":
            hh = [a for a in e.households()]
            money = np.array([a.get("money") for a in hh])
            vals.append(_gini(money))
        else:
            vals.append(float(macros[aggregate]))
    return float(np.mean(vals))


def make_oracle(target_id: str):
    """oracle: X (n, d) -> y (n,). Deterministic; per-query seed fixed."""
    spec = TARGETS[target_id]
    inputs, aggregate = spec["inputs"], spec["aggregate"]

    def oracle(X: np.ndarray) -> np.ndarray:
        X = np.atleast_2d(np.asarray(X, float))
        return np.array([_steady(inputs, row, aggregate) for row in X])

    return oracle

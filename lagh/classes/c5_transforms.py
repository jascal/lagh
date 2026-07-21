"""C5: target transforms -- fit t(y) with the whole lower-tier stack, invert.

Reaches sqrt(rational), 1/(e^u - 1), arcsin(linear in sin), and exponential-family
laws while fitting stays a linear solve throughout. Guards are part of each
transform's admissibility; inverses that produce infinities are dropped upstream.
"""

from __future__ import annotations

import numpy as np
import sympy as sp

TIER = 5


def transforms(y: np.ndarray):
    y = np.asarray(y, float).ravel()
    out = []
    if np.all(y > 0):
        out.append(("log", np.log(y), lambda e: sp.exp(e)))
        out.append(("inv", 1.0 / y, lambda e: 1 / e))
    out.append(("square", y * y, lambda e: sp.sqrt(e)))
    if np.all(np.abs(y) <= np.pi / 2):
        out.append(("sin", np.sin(y), lambda e: sp.asin(e)))
    # exp-of-target: reaches log(sum of exps) forms (logaddexp family). exp(y) is a
    # clean C4 sum e^x0 + e^x1; invert with log. Guarded so exp(y) does not overflow.
    if np.all(y <= 30):
        out.append(("exp", np.exp(y), lambda e: sp.log(e)))
    return out


def apply(name: str, y: np.ndarray) -> np.ndarray:
    return {"log": np.log, "inv": lambda v: 1.0 / v, "square": lambda v: v * v,
            "sin": np.sin, "exp": np.exp}[name](np.asarray(y, float))

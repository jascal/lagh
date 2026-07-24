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
    # CAP-D: sqrt-FORWARD (fit sqrt(y), invert by squaring) reaches y = (rational)^2
    # and any y = g(x)^2 where g is lower-tier -- the C2 rational base becomes visible
    # under the root. Guarded y>=0. Complements "square" (which fits y^2 for y=sqrt(.)).
    if np.all(y >= 0):
        out.append(("sqrt", np.sqrt(y), lambda e: e * e))
    if np.all(np.abs(y) <= np.pi / 2):
        out.append(("sin", np.sin(y), lambda e: sp.asin(e)))
    # exp-of-target: reaches log(sum of exps) forms (logaddexp family). exp(y) is a
    # clean C4 sum e^x0 + e^x1; invert with log. Guarded so exp(y) does not overflow.
    if np.all(y <= 30):
        out.append(("exp", np.exp(y), lambda e: sp.log(e)))
    # CAP-E: rational-of-exponential (Bose-Einstein / Fermi-Dirac). np.log1p is the
    # registered stability fix -- log(1/y + 1) cancelled catastrophically at y~1e10
    # and manufactured the program's only confident-wrong; log1p(1/y) does not. The
    # exact-coefficient gate (the other registered precondition) now rejects any
    # loose-epsilon float certificate this transform could produce.
    if np.all(y > 0):
        out.append(("bose", np.log1p(1.0 / y), lambda e: 1 / (sp.exp(e) - 1)))
    if np.all((y > 0) & (y < 1)):
        out.append(("fermi", np.log((1.0 - y) / y), lambda e: 1 / (sp.exp(e) + 1)))
    return out


def apply(name: str, y: np.ndarray) -> np.ndarray:
    return {"log": np.log, "inv": lambda v: 1.0 / v, "square": lambda v: v * v,
            "sqrt": np.sqrt, "sin": np.sin, "exp": np.exp,
            "bose": lambda v: np.log1p(1.0 / v),
            "fermi": lambda v: np.log((1.0 - v) / v)}[name](np.asarray(y, float))

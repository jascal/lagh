"""The lab's hidden ground-truth problems + public cards.

`fn` and `truth` are NEVER sent over MCP -- only `card()` (inputs, domain, hint) and the
sampled `y` values are. Offline scoring compares the model's reported law to `truth`.
"""

from __future__ import annotations

import numpy as np

# id -> (input names, suggested domain [lo, hi], hint, hidden fn, hidden truth, expected)
PROBLEMS: dict[str, dict] = {
    "orbit": {
        "inputs": ["m1", "m2", "r"],
        "domain": [[1.0, 1.0, 0.5], [10.0, 10.0, 5.0]],
        "hint": "attractive force between two masses separated by a distance",
        "fn": lambda X: 6.674e-11 * X[:, 0] * X[:, 1] / X[:, 2] ** 2,
        "truth": "6.674e-11*m1*m2/r**2",
        "expected": "recover (clean power law; a physical constant coefficient)",
    },
    "spring": {
        "inputs": ["k", "x"],
        "domain": [[1.0, 0.1], [10.0, 3.0]],
        "hint": "energy stored when an elastic element is displaced",
        "fn": lambda X: 0.5 * X[:, 0] * X[:, 1] ** 2,
        "truth": "k*x**2/2",
        "expected": "recover (clean; coefficient 1/2)",
    },
    "polarizer": {
        "inputs": ["I0", "theta"],
        "domain": [[1.0, 0.05], [10.0, 1.4]],
        "hint": "intensity transmitted through a polarizer at angle theta (radians)",
        "fn": lambda X: X[:, 0] * np.cos(X[:, 1]) ** 2,
        "truth": "I0*cos(theta)**2",
        "expected": "recover (trig; the half-angle form is exact)",
    },
    "narrow": {
        "inputs": ["x"],
        "domain": [[0.9], [1.1]],  # deliberately too tight -> under-determined
        "hint": "a smooth response; the suggested range may be too narrow to pin it",
        "fn": lambda X: 3.0 * X[:, 0] ** 3 - 2.0 * X[:, 0],
        "truth": "3*x**3 - 2*x",
        "expected": "fit says acquire_more_data; recover after broadening the box",
    },
    "steps": {
        "inputs": ["x"],
        "domain": [[0.5], [6.0]],
        "hint": "a quantity that changes in discrete jumps",
        "fn": lambda X: np.floor(2.0 * X[:, 0]),
        "truth": "floor(2*x)  [non-smooth]",
        "expected": "recover ABSTAINS (non-smooth / not a closed form) -- the honest answer",
    },
    "exotic": {
        "inputs": ["x"],
        "domain": [[1.0], [4.0]],
        "hint": "a power law whose exponent may not be a simple fraction",
        "fn": lambda X: 2.0 * X[:, 0] ** np.e,
        "truth": "2*x**e  [irrational exponent]",
        "expected": "recover ABSTAINS; fit flags a continuum near e -> declare_and_verify",
    },
}

MAX_POINTS = 2000


def card(pid: str) -> dict:
    """Public problem card -- no hidden fn/truth."""
    p = PROBLEMS[pid]
    return {"id": pid, "inputs": p["inputs"], "suggested_domain": p["domain"],
            "hint": p["hint"], "n_inputs": len(p["inputs"])}


def sample(pid: str, X) -> list[float]:
    """Query the hidden oracle at points X (list of rows). Noise-free."""
    if pid not in PROBLEMS:
        raise KeyError(f"unknown problem {pid!r}; known: {sorted(PROBLEMS)}")
    X = np.asarray(X, float)
    if X.ndim == 1:
        X = X[:, None]
    d = len(PROBLEMS[pid]["inputs"])
    if X.shape[1] != d:
        raise ValueError(f"problem {pid!r} expects {d} input(s) per row, got {X.shape[1]}")
    if len(X) > MAX_POINTS:
        raise ValueError(f"too many points ({len(X)} > {MAX_POINTS})")
    y = np.asarray(PROBLEMS[pid]["fn"](X), float).ravel()
    return [float(v) for v in y]


def truth(pid: str) -> str:
    """Hidden true form -- OFFLINE scoring only, never exposed over MCP."""
    return PROBLEMS[pid]["truth"]

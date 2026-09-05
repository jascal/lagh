"""Registered public-surface refusal read; no oracle or paid proposer."""
import json
from pathlib import Path

import numpy as np

from lagh.mcp.core import verify


def cases():
    x = np.linspace(1, 2, 100)
    yield "omitted_quadratic", x, 2*x + .1*x*x
    y = 2*x
    a, b = np.random.default_rng(0).permutation(100)[:2]
    y[a] += .1
    y[b] -= .1*x[a]/x[b]  # preserve training scale; mismatch outside check split
    yield "fitting_rows_only", x, y
    yield "clean_control", x, 2*x


if __name__ == "__main__":
    out = {name: verify(x, y, "x_0", floor_abs=1e-6)
           for name, x, y in cases()}
    path = Path("experiments/results/refusal_measurements.json")
    path.write_text(json.dumps({"tag": "empirical", "cases": out}, indent=2,
                               allow_nan=False) + "\n")
    print(json.dumps({k: {"certified": v["certified"],
                         "abstain": v.get("abstain"),
                         "domain": v.get("measurement", {}).get("domain")}
                      for k, v in out.items()}, indent=2))

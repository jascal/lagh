"""Registered public-surface refusal read; no oracle or paid proposer."""
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import numpy as np

from lagh.mcp.core import verify


def cases():
    x = np.linspace(1, 2, 100)
    yield "omitted_quadratic", x, 2*x + .1*x*x
    y = 2*x
    a, b = 82, 36  # frozen fitting rows under the registered seed-0 split
    y[a] += .1
    y[b] -= .1*x[a]/x[b]  # preserve training scale; mismatch outside check split
    yield "fitting_rows_only", x, y
    yield "clean_control", x, 2*x


if __name__ == "__main__":
    out = {name: verify(x, y, "x_0", floor_abs=1e-6)
           for name, x, y in cases()}
    assert out["fitting_rows_only"]["measurement"]["domain"] == (
        "all finite rows of supplied dataset"), "registered split fixture changed"
    path = Path(__file__).resolve().parent / "results/refusal_measurements.json"
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps({"tag": "empirical", "cases": out}, indent=2,
                               allow_nan=False) + "\n")
    print(json.dumps({k: {"certified": v["certified"],
                         "abstain": v.get("abstain"),
                         "domain": v.get("measurement", {}).get("domain")}
                      for k, v in out.items()}, indent=2))

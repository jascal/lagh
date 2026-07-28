"""The reach-envelope audit (registered): a systematic matrix of
(dimension x form family) cells, each a clean synthetic recovery at machine
floor. Records certified / abstain-with-reason per cell — the measured
proposal-reach envelope, so future gaps are known caps rather than
surprises (the RRab lesson: the binding constraint is PROPOSAL reach).

Every cell: n=400 points, inputs uniform in [0.5, 3] (positive: power/log
safe), exact targets, sigma=0, default floor. A cell FAILS REACH when the
true law certifies nowhere in the escalation (structural 'no law
certifies'); ambiguity abstains are recorded distinctly.
"""
from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).parent.parent.parent))
from lagh.passive import discover_passive  # noqa: E402

OUT = Path("experiments/results/reach_audit.json")


def X(dim, n=400, seed=0):
    """Per-cell seeding: every cell is standalone-reproducible (the shared
    stream made filtered re-runs draw different data than the full run)."""
    return np.random.default_rng(seed).uniform(0.5, 3.0, (n, dim))


CELLS = []


def cell(name, dim, fn):
    CELLS.append((name, dim, fn))


# --- linear family: k-term affine at increasing dim
cell("linear-2term-d1", 1, lambda x: 3 * x[:, 0] - 2)
cell("linear-3term-d2", 2, lambda x: 2 * x[:, 0] - 5 * x[:, 1] + 1)
cell("linear-4term-d3", 3, lambda x: x[:, 0] + 2 * x[:, 1] - 3 * x[:, 2] + 4)
cell("linear-5term-d4", 4,
     lambda x: x[:, 0] - x[:, 1] + 2 * x[:, 2] - 3 * x[:, 3] + 5)
cell("linear-6term-d5", 5,
     lambda x: x[:, 0] + x[:, 1] - x[:, 2] + 2 * x[:, 3] - 2 * x[:, 4] + 3)
cell("linear-7term-d6", 6,
     lambda x: (x[:, 0] - x[:, 1] + x[:, 2] + 2 * x[:, 3] - x[:, 4]
                + 3 * x[:, 5] - 1))

# --- complete quadric at increasing dim (CAP-T covers d<=4)
cell("quadric-full-d2", 2,
     lambda x: (1 + x[:, 0] - x[:, 1] + 2 * x[:, 0]**2 - x[:, 1]**2
                + 3 * x[:, 0] * x[:, 1]))
cell("quadric-full-d3", 3,
     lambda x: (2 - x[:, 0] + x[:, 1] + x[:, 2] + x[:, 0]**2
                + 2 * x[:, 1] * x[:, 2] - x[:, 0] * x[:, 2]))
cell("quadric-full-d4", 4,
     lambda x: (1 + sum(x[:, j] for j in range(4))
                + x[:, 0] * x[:, 1] - x[:, 2] * x[:, 3] + x[:, 3]**2))
cell("quadric-full-d5", 5,
     lambda x: (1 + sum(x[:, j] for j in range(5))
                + x[:, 0] * x[:, 4] - x[:, 1] * x[:, 3] + x[:, 2]**2))

# --- cubic cross terms (CAP-T cubic covers d<=2)
cell("cubic-full-d1", 1, lambda x: x[:, 0]**3 - 2 * x[:, 0]**2 + x[:, 0] - 1)
cell("cubic-cross-d2", 2,
     lambda x: x[:, 0]**3 - x[:, 0]**2 * x[:, 1] + 2 * x[:, 1] - 1)
cell("cubic-cross-d3", 3,
     lambda x: x[:, 0]**2 * x[:, 1] - x[:, 1] * x[:, 2] + x[:, 2]**3)

# --- monomials / power laws
cell("monomial-d2", 2, lambda x: 3.5 * x[:, 0]**2 / x[:, 1])
cell("monomial-d3", 3, lambda x: 2 * x[:, 0] * x[:, 1]**2 / x[:, 2]**3)
cell("monomial-d4", 4,
     lambda x: x[:, 0] * x[:, 1] / (x[:, 2] * x[:, 3]))
cell("monomial-frac-d1", 1, lambda x: 2 * x[:, 0]**1.5)
cell("monomial-frac-d2", 2, lambda x: x[:, 0]**1.5 / x[:, 1]**0.5)

# --- rational forms (affine denominators)
cell("rational-d1", 1, lambda x: (2 * x[:, 0] + 1) / (x[:, 0] + 3))
cell("rational-d2", 2, lambda x: x[:, 0] / (x[:, 1] + 1))
cell("rational-d3", 3, lambda x: (x[:, 0] + x[:, 1]) / (x[:, 2] + 2))

# --- transcendental sums
cell("trig-sum-d1", 1, lambda x: 2 * np.sin(x[:, 0]) + np.cos(x[:, 0]))
cell("trig-prod-d2", 2, lambda x: x[:, 0] * np.sin(x[:, 1]))
cell("trig-cross-d2", 2, lambda x: np.sin(x[:, 0]) * np.cos(x[:, 1]))
cell("exp-decay-d1", 1, lambda x: 3 * np.exp(-x[:, 0]))
cell("exp-decay-d2", 2, lambda x: x[:, 0] * np.exp(-x[:, 1]))
cell("log-sum-d2", 2, lambda x: 2 * np.log(x[:, 0]) + x[:, 1])
cell("mixed-4term-d2", 2,
     lambda x: np.sin(x[:, 0]) + np.log(x[:, 1]) + x[:, 0] * x[:, 1] - 2)

# --- sparse sums at size 3..6 over a wide basis (the oscillator lesson)
cell("sparse3-d1", 1,
     lambda x: 2 * np.sin(x[:, 0]) + x[:, 0]**2 - 3 / x[:, 0])
cell("sparse4-d1", 1,
     lambda x: np.sin(x[:, 0]) + np.cos(x[:, 0]) + x[:, 0] - np.exp(-x[:, 0]))
cell("sparse5-d2", 2,
     lambda x: (np.sin(x[:, 0]) + np.log(x[:, 1]) + x[:, 0]**2
                - x[:, 1] + 2 / x[:, 0]))
cell("sparse6-d2", 2,
     lambda x: (np.sin(x[:, 0]) + np.cos(x[:, 1]) + x[:, 0] * x[:, 1]
                - np.sqrt(x[:, 0]) + np.exp(-x[:, 1]) - 1))

# --- compositions / inner transforms
cell("gaussian-d1", 1, lambda x: np.exp(-(x[:, 0] - 1.5)**2))
cell("inner-linear-trig-d1", 1, lambda x: np.sin(2 * x[:, 0] + 1))
cell("sqrt-sum-d2", 2, lambda x: np.sqrt(x[:, 0]**2 + x[:, 1]**2))
cell("inverse-square-sum-d2", 2, lambda x: 1.0 / (x[:, 0]**2 + x[:, 1]**2))


def main(only=None):
    results = {}
    for name, dim, fn in CELLS:
        if only and only not in name:
            continue
        import zlib
        x = X(dim, seed=zlib.crc32(name.encode()))
        y = fn(x)
        t0 = time.time()
        try:
            r = discover_passive(x, y, sigma=0.0)
            cert = r.certified
            law = str(r.result.expr)[:90] if cert else ""
            ab = r.result.certificate.abstain
            notes = [str(n)[:70] for n in r.result.certificate.notes][:1]
        except Exception as e:                                 # noqa: BLE001
            cert, law, ab, notes = False, "", f"crash: {e}", []
        results[name] = {"dim": dim, "certified": bool(cert), "law": law,
                         "abstain": ab, "note": notes,
                         "seconds": round(time.time() - t0, 1)}
        print(f"{name:26s} {'CERT' if cert else 'MISS':4s} "
              f"{ab or '':12s} {results[name]['seconds']:6.1f}s  {law[:50]}",
              flush=True)
    OUT.write_text(json.dumps(results, indent=1))
    n_cert = sum(1 for v in results.values() if v["certified"])
    print(f"\nREACH: {n_cert}/{len(results)} cells certified -> {OUT}")


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else None)

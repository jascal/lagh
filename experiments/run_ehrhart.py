"""Run the registered Ehrhart targets H1-H3 (docs/TESTBED_EHRHART.md).

Exact integer oracle; exact quasi-polynomial recovery; zero-wrong invariant checked
against an extended-range reference that interpolation cannot game.
"""

from __future__ import annotations

import json
import sys
import time
from fractions import Fraction
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent))

import numpy as np

from lagh.adapters.ehrhart import lattice_count, random_simplex
from lagh.quasipoly import recover

T_MAX = 48
EXT = list(range(T_MAX + 1, T_MAX + 41))       # far beyond fit range
N_PER_DIM = 10


def extended_reference_ok(qp, intercepts) -> bool:
    for t in EXT:
        if qp(t) != Fraction(lattice_count(intercepts, t)):
            return False
    return True


def main() -> int:
    out = Path("experiments/results"); out.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(2026)
    rows = []
    for dim, tag in [(1, "H1"), (2, "H2"), (3, "H3")]:
        for k in range(N_PER_DIM):
            inter = random_simplex(dim, rng)
            ts = list(range(1, T_MAX + 1))
            t0 = time.time()
            Ls = [lattice_count(inter, t) for t in ts]
            r = recover(ts, Ls)
            correct = (extended_reference_ok(r.quasipoly, inter)
                       if r.certified else None)
            rec = {"target": tag, "dim": dim,
                   "intercepts": [str(a) for a in inter],
                   "certified": r.certified, "abstain": r.abstain,
                   "note": r.note, "domain_size": r.domain_size,
                   "ext_ref_correct": correct,
                   "confident_wrong": bool(r.certified and correct is False),
                   "seconds": round(time.time() - t0, 2)}
            rows.append(rec)
            flag = ("CERT " + ("ok" if correct else "**WRONG**")) if r.certified \
                else f"abstain[{r.abstain}]"
            print(f"{tag} {str([str(a) for a in inter]):32s} {flag:16s} {r.note}",
                  flush=True)

    (out / "ehrhart.jsonl").write_text("\n".join(json.dumps(r) for r in rows) + "\n")
    for tag in ("H1", "H2", "H3"):
        g = [r for r in rows if r["target"] == tag]
        c = sum(r["certified"] and r["ext_ref_correct"] for r in g)
        print(f"{tag}: {c}/{len(g)} correct, "
              f"{sum(r['abstain'] is not None for r in g)} abstain")
    cw = sum(r["confident_wrong"] for r in rows)
    print(f"\nconfident-wrong: {cw}  (invariant: must be 0)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

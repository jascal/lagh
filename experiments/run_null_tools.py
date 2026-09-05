"""Null calibration at the TOOL surface (the companion of run_null_calibration.py).

The engine's null sweep (0/200 true-random targets through `discover`) never
exercises the public tools, and the review of 2026-09-04 (jascal/lagh#1-#6)
found five confident-wrong paths in exactly the code it skips: `verify`,
box-search, the C6 tier, state certificates. This sweep puts TRUE-RANDOM
targets through the tools themselves -- `recover` on data, `verify` on a set
of declared forms, `recover` in box-search mode against an oracle that answers
at random, and `certify_state` on a random target -- and counts certificates.
Every count must be zero; a single `proved` falsifies the tool's accounting.

OS-seeded (a true null, never replayable), like the engine sweep. Slow by
nature: a `recover` on random data walks every tier before it abstains.

    .venv/bin/python experiments/run_null_tools.py --trials 20
"""
from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).parent.parent))
from lagh.mcp.core import recover, verify          # noqa: E402
from lagh.statecert import certify_state           # noqa: E402

FORMS = ["x_0", "x_0**2", "1/x_0", "sqrt(x_0)", "log(x_0)", "exp(x_0)", "x_0**E",
         "x_0*x_1", "x_0/x_1", "x_0 + x_1", "0", "1"]


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default="experiments/results/null_tools.jsonl")
    ap.add_argument("--trials", type=int, default=20)
    ap.add_argument("--skip-boxsearch", action="store_true")
    args = ap.parse_args(argv)
    rng = np.random.default_rng()
    rows, certs = [], {"recover": 0, "verify": 0, "boxsearch": 0, "state": 0}
    for t in range(args.trials):
        dim = int(rng.integers(1, 4))
        n = int(rng.integers(8, 61))
        X = np.exp(rng.uniform(np.log(0.5), np.log(10), (n, dim)))
        y = rng.uniform(-1, 1, n) * 10 ** rng.uniform(-3, 3)
        row = {"trial": t, "dim": dim, "n": n}
        t0 = time.time()
        r = recover(X, y)
        row["recover"] = {"certified": r["certified"], "abstain": r.get("abstain"),
                          "law": r.get("law"), "secs": round(time.time() - t0, 1)}
        certs["recover"] += bool(r["certified"])
        vs = {}
        for f in FORMS:
            if dim == 1 and "x_1" in f:
                continue
            v = verify(X, y, f)
            vs[f] = {"certified": v["certified"], "abstain": v.get("abstain")}
            certs["verify"] += bool(v["certified"])
        row["verify"] = vs
        if not args.skip_boxsearch and t % 4 == 0:
            scale = 10 ** rng.uniform(-2, 2)
            oracle = lambda Xq, s=scale: rng.uniform(-1, 1, len(Xq)) * s   # noqa: E731
            t0 = time.time()
            b = recover(oracle=oracle, box=[[0.5] * dim, [10.0] * dim], box_search=True,
                        time_budget_s=60.0)
            row["boxsearch"] = {"certified": b["certified"], "abstain": b.get("abstain"),
                                "heldout_box_ok": b["acquisition"].get("heldout_box_ok"),
                                "secs": round(time.time() - t0, 1)}
            certs["boxsearch"] += bool(b["certified"])
        d = int(rng.integers(1, 5))
        nr = int(rng.integers(d + 8, 60))
        B = rng.normal(size=(nr, d)) * 10 ** rng.uniform(-3, 3, d)
        yr = rng.normal(size=nr)
        eps = 10 ** rng.uniform(-6, -1) * np.ones(nr)
        c = certify_state(B, yr, eps, [f"m{j}" for j in range(d)])
        row["state"] = {"certified": c.certified, "abstain": c.abstain,
                        "alpha_log10": c.alpha_log10}
        certs["state"] += bool(c.certified)
        rows.append(row)
        print(f"[{t + 1:>3}/{args.trials}] dim={dim} n={n:>2} recover={r.get('abstain')} "
              f"({row['recover']['secs']}s) verify certs={sum(v['certified'] for v in vs.values())} "
              f"state={c.abstain or 'CERT'} | running certs: {certs}", flush=True)
    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    with out.open("a") as fh:
        for row in rows:
            fh.write(json.dumps(row) + "\n")
    print("\nfalse certifications:", certs, "over", args.trials, "trials")
    verdict = "PASS" if not any(certs.values()) else "FAIL -- the tool accounting is falsified"
    print(verdict)
    return 0 if verdict == "PASS" else 1


if __name__ == "__main__":
    sys.exit(main())

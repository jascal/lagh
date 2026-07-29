"""Null calibration of the significance bound (DIRECTION_SIGNIFICANCE.md job 1).

Run discovery on TRUE-RANDOM targets (i.i.d. uniform y over a range; inputs
log-uniform as everywhere else). The bound alpha <= |H|*q^h claims a spurious
certification is essentially impossible (alpha ~ 1e-50-ish at machine epsilon);
therefore the observed false-certification count over N trials MUST be zero,
and a single certification falsifies the |H| accounting (the transform/inner
tiers would have been undercounted). Clean sigma=0 regime first.
"""
from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).parent.parent))
from lagh.engine import discover  # noqa: E402

N_TRIALS = 200


def main(argv=None):
    # The re-run under the AMENDED arbitration rule (the evidence-starved rival
    # bar, 2026-07-28) writes beside the original rather than over it: the
    # amendment can only ever REDUCE certifications, so the comparison against
    # the pre-amendment 0/200 is the whole point of running it again.
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default="experiments/results/null_calibration.jsonl")
    ap.add_argument("--trials", type=int, default=N_TRIALS)
    args = ap.parse_args(argv)
    rng = np.random.default_rng()          # OS-seeded: a true null, not replayable
    out = Path(args.out)
    n_trials = args.trials
    rows = []
    certs = 0
    for t in range(n_trials):
        dim = int(rng.integers(1, 4))
        n = int(rng.integers(60, 200))
        lo, hi = 0.5, 10.0
        X = np.exp(rng.uniform(np.log(lo), np.log(hi), (n, dim)))
        scale = 10.0 ** rng.uniform(-3, 3)
        y = rng.uniform(-scale, scale, n)          # i.i.d. -- the null
        a, b = int(0.6 * n), int(0.8 * n)
        t0 = time.time()
        r = discover(X[:a], y[:a], X[a:b], y[a:b], X[b:], y[b:])
        rec = {"trial": t, "dim": dim, "n": n,
               "certified": bool(r.certificate.certified),
               "alpha_log10": r.certificate.alpha_log10,
               "n_hypotheses": r.certificate.n_hypotheses,
               "law": str(r.expr)[:80] if r.expr is not None else None,
               "secs": round(time.time() - t0, 1)}
        rows.append(rec)
        if rec["certified"]:
            certs += 1
            print(f"!! FALSE CERTIFICATION trial {t}: {rec}", flush=True)
        if (t + 1) % 20 == 0:
            print(f"{t+1}/{n_trials} trials, false certs: {certs}", flush=True)
    out.write_text("\n".join(json.dumps(r) for r in rows) + "\n")
    print(f"\nNULL CALIBRATION: {certs}/{n_trials} false certifications "
          f"(bound demands 0; any hit falsifies the |H| accounting)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

"""Measure lagh against NewtonBench vanilla-equation cells (docs, forthcoming reg).

CAPABILITY READ first: easy cells (exposure-disclosed -- wyly N1 read them). This is
not the blind claim; it measures where the full lagh instrument now stands on the real
benchmark laws. Scored by dense-grid reference (oracle only), 1e-6 relative.
"""

from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent))

import numpy as np
import sympy as sp

from lagh.acquisition import run_active
from lagh.adapters.newtonbench import MODULES, available_versions, make_oracle


def dense_ref_ok(expr, oracle, lo, hi, dim, tol=1e-6, n=200):
    if expr is None:
        return None
    rng = np.random.default_rng(4321)
    lo, hi = np.array(lo), np.array(hi)
    X = np.exp(rng.uniform(np.log(np.maximum(lo, 1e-6)), np.log(hi), (n, dim)))
    y = oracle(X)
    syms = [sp.Symbol(f"x_{i}") for i in range(dim)]
    try:
        got = np.broadcast_to(np.asarray(sp.lambdify(syms, expr, "numpy")(*X.T),
                                         float), y.shape)
    except Exception:                                        # noqa: BLE001
        return False
    ok = np.isfinite(got) & np.isfinite(y) & (np.abs(y) > 1e-9)
    if ok.sum() < n // 2:
        return False
    return bool(np.max(np.abs(got[ok] - y[ok]) / np.abs(y[ok])) < tol)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--difficulty", default="easy")
    ap.add_argument("--out", type=Path,
                    default=Path("experiments/results/newtonbench_easy.jsonl"))
    args = ap.parse_args()
    out = args.out; out.parent.mkdir(parents=True, exist_ok=True)
    rows = []
    for module, (inputs, lo, hi) in MODULES.items():
        dim = len(inputs)
        for v in available_versions(module, args.difficulty):
            oracle = make_oracle(module, v)
            t0 = time.time()
            r = run_active(oracle, lo, hi, seed=1)
            expr = r.result.expr
            correct = dense_ref_ok(expr, oracle, lo, hi, dim)
            rec = {"module": module, "version": v, "dim": dim,
                   "certified": r.result.certificate.certified,
                   "abstain": r.result.certificate.abstain,
                   "law": str(expr) if expr is not None else None,
                   "correct": correct,
                   "confident_wrong": bool(r.result.certificate.certified
                                           and correct is False),
                   "tier": r.result.tier, "queries": r.queries_used,
                   "seconds": round(time.time() - t0, 1)}
            rows.append(rec)
            flag = ("CERT " + ("ok" if correct else "**WRONG**")) \
                if r.result.certificate.certified else f"abstain[{r.result.certificate.abstain}]"
            print(f"{module:24s} {v} {flag:16s} tier={r.result.tier} "
                  f"q={r.queries_used} ({rec['seconds']}s)", flush=True)
            if rec["law"]:
                print(f"     {rec['law'][:78]}", flush=True)

    out.write_text("\n".join(json.dumps(r) for r in rows) + "\n")
    rec_ok = sum(r["certified"] and r["correct"] for r in rows)
    cw = sum(r["confident_wrong"] for r in rows)
    n = len(rows)
    # per-module recovery (any version)
    mods = {}
    for r in rows:
        mods.setdefault(r["module"], []).append(r["certified"] and r["correct"])
    modrec = sum(any(v) for v in mods.values())
    print(f"\nrecovered {rec_ok}/{n} tasks | {modrec}/{len(mods)} modules (any version)")
    print(f"confident-wrong {cw} (must be 0)")
    print(f"bars: GPT-5 0.903 / best-published 0.965  ->  we are at {rec_ok/n:.3f}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

"""Re-score every published result that went through the MCP `verify` tool, after
jascal/lagh#4 gave that tool the certification gates discovery had always run
(vacuity, coefficient, parametric, full-data, significance). A changed gate means
a re-scored ledger. No LLM call is made: the verified forms are stored, and the
data is either regenerated from the NewtonBench oracles (machine composite) or
loaded from LLM-SRBench (dev sweep v1), with the same seeds and subsampling the
original runs used.

    .venv/bin/python experiments/reverify_gated_verify.py

writes experiments/results/reverify_gated_verify.json.
"""
from __future__ import annotations

import collections
import json
import os
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).parent.parent))
sys.path.insert(0, "/home/allans/code/llm-srbench")
os.environ.setdefault("LAB_SOURCE", "newtonbench")

OUT = Path("experiments/results/reverify_gated_verify.json")


def machine_composite():
    """Every `proved` law in machine/run_results.json, re-verified on the sample the
    driver's verifying state draws (seed 0 + 2)."""
    from lagh.adapters.newtonbench import MODULES, make_oracle
    from lagh.mcp.core import verify
    from machine.driver import _sample
    rows = json.load(open("machine/run_results.json"))
    out = []
    for r in rows:
        if r["outcome"] != "proved" or not r["law"]:
            continue
        m, d, v = r["pid"].split("/")
        _, lo, hi = MODULES[m]
        X, y = _sample(make_oracle(m, v, d), [list(map(float, lo)), list(map(float, hi))],
                       seed=2)
        vr = verify(X.tolist(), y.tolist(), r["law"])
        out.append({"pid": r["pid"], "cat": r["cat"], "lagh_alone": r["lagh_alone"],
                    "law": r["law"], "now_certified": vr["certified"],
                    "abstain": vr.get("abstain"), "alpha_log10": vr.get("alpha_log10"),
                    "note": vr["note"][:200]})
        print(f"  {r['pid']:<28} {'stands ' if vr['certified'] else 'DEMOTED'} "
              f"{vr.get('abstain') or '':<10} {r['law'][:48]}", flush=True)
    return out


def dev_sweep_v1():
    """The 19 `llm-verified` certificates of dev sweep v1, re-verified on the same
    SUB_CAP subsample at the same declared sigma_rep."""
    import experiments.run_dev_llmsrbench as dev
    from experiments.run_blind_llmsrbench import _train_xy, load_all_problems
    from lagh.mcp.core import verify
    byid = {f"{m}/{p.equation_idx}": p for m, p in load_all_problems()}
    rows = [json.loads(l) for l in open("experiments/results/dev_llmsrbench_v1.jsonl")]
    out = []
    for r in rows:
        if r["channel"] != "llm-verified":
            continue
        p = byid[r["id"]]
        X, y = _train_xy(p)
        if len(X) > dev.SUB_CAP:
            idx = np.sort(np.random.default_rng(0).choice(len(X), dev.SUB_CAP,
                                                          replace=False))
            X, y = X[idx], y[idx]
        m = np.isfinite(y)
        v = verify(X[m].tolist(), y[m].tolist(), r["expr"], sigma=dev.SIGMA_REP)
        out.append({"id": r["id"], "expr": r["expr"], "now_certified": v["certified"],
                    "abstain": v.get("abstain"), "alpha_log10": v.get("alpha_log10"),
                    "law": v.get("law"), "note": v["note"][:200]})
        print(f"  {r['id']:<34} {'stands ' if v['certified'] else 'DEMOTED'} "
              f"alpha<=1e{v.get('alpha_log10', float('nan')):.0f}  {r['expr'][:46]}",
              flush=True)
    return out


def main():
    print("machine composite (machine/run_results.json), proved laws:")
    mc = machine_composite()
    print("\ndev sweep v1, llm-verified certificates:")
    dv = dev_sweep_v1()
    summary = {
        "machine": {"n_proved": len(mc),
                    "stand": sum(o["now_certified"] for o in mc),
                    "demoted": [o["pid"] for o in mc if not o["now_certified"]],
                    "demoted_reasons": dict(collections.Counter(
                        o["abstain"] for o in mc if not o["now_certified"]))},
        "dev_v1": {"n_llm_verified": len(dv),
                   "stand": sum(o["now_certified"] for o in dv),
                   "demoted": [o["id"] for o in dv if not o["now_certified"]],
                   "alpha_log10_range": [min(o["alpha_log10"] for o in dv
                                             if o["alpha_log10"] is not None),
                                         max(o["alpha_log10"] for o in dv
                                             if o["alpha_log10"] is not None)]},
    }
    json.dump({"summary": summary, "machine": mc, "dev_v1": dv},
              open(OUT, "w"), indent=1)
    print("\n" + json.dumps(summary, indent=1))


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Run the LawResearch machine composite over NewtonBench-dev and score it.

No Hermes, no Docker, no subprocess, no timeouts -- the verified state machine bounds the
loop in-process, so every cell terminates in the certify|characterize ladder. The LLM is a
bounded proposer (machine/llm.py) read from machine/.env; unconfigured -> lagh alone.

    .venv/bin/python machine/run_bench.py --subset easy
    .venv/bin/python machine/run_bench.py --subset all
    NEWTONBENCH_DIR=/path .venv/bin/python machine/run_bench.py --cells m4_snell_law/easy/v0
"""
from __future__ import annotations

import argparse
import json
import os
import time
from collections import Counter

import numpy as np
import sympy as sp

os.environ.setdefault("LAB_SOURCE", "newtonbench")
from lagh.adapters.newtonbench import MODULES, make_oracle          # noqa: E402
from machine.driver import run_sync                                 # noqa: E402
from machine.llm import _load_env, propose_forms                     # noqa: E402

_load_env()  # load machine/.env into os.environ BEFORE the have_llm gate reads it


def _all_cells() -> list[str]:
    out = []
    for m in MODULES:
        for d in ("easy", "medium", "hard"):
            for v in ("v0", "v1", "v2"):
                out.append(f"{m}/{d}/{v}")
    return out


def score(pid: str, law: str):
    """Dense-compare a proved law (already in x_0..) to the real oracle. (correct, rel_err)
    or None if it can't be parsed/evaluated (a formatting issue, never a confident-wrong)."""
    m, d, v = pid.split("/")
    inp, lo, hi = MODULES[m]
    dim = len(inp)
    try:
        expr = sp.sympify(law)
    except Exception:                                               # noqa: BLE001
        return None
    orc = make_oracle(m, v, d)
    rng = np.random.default_rng(20260722)
    X = np.exp(rng.uniform(np.log(np.maximum(np.array(lo, float), 1e-9)),
                           np.log(np.array(hi, float)), (200, dim)))
    y = np.asarray(orc(X), float)
    try:
        f = sp.lambdify([sp.Symbol(f"x_{i}") for i in range(dim)], expr, "numpy")
        got = np.broadcast_to(np.asarray(f(*X.T), float), y.shape)
    except Exception:                                               # noqa: BLE001
        return None
    ok = np.isfinite(got) & np.isfinite(y) & (np.abs(y) > 1e-9)
    if ok.sum() < 20:
        return None
    return float(np.median(np.abs(got[ok] - y[ok]) / np.abs(y[ok]))) < 1e-2, None


def baseline_map() -> dict:
    path = os.path.join(os.path.dirname(__file__), "..", "experiments", "results",
                        "newtonbench_all.jsonl")
    out = {}
    if os.path.exists(path):
        for ln in open(path):
            r = json.loads(ln)
            out[f"{r['module']}/{r['difficulty']}/{r['version']}"] = bool(
                r["certified"] and r["correct"])
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--subset", default="easy", help="easy|medium|hard|all")
    ap.add_argument("--cells", default="")
    ap.add_argument("--limit", type=int, default=0)
    args = ap.parse_args()

    ids = _all_cells()
    if args.cells:
        want = set(args.cells.split(","))
        ids = [i for i in ids if i in want]
    elif args.subset != "all":
        ids = [i for i in ids if i.split("/")[1] == args.subset]
    if args.limit:
        ids = ids[: args.limit]

    base = baseline_map()
    have_llm = bool(os.environ.get("LLM_API_KEY") and os.environ.get("LLM_MODEL"))
    print(f"machine composite over {len(ids)} cell(s) [subset={args.subset}] "
          f"-- proposer: {'LLM ' + os.environ.get('LLM_MODEL', '') if have_llm else 'lagh alone'}\n",
          flush=True)
    rows = []
    for k, pid in enumerate(ids, 1):
        m, d, v = pid.split("/")
        inp, lo, hi = MODULES[m]
        oracle = make_oracle(m, v, d)
        box = [list(map(float, lo)), list(map(float, hi))]
        t0 = time.time()
        # always pass propose_form -- it self-degrades to None when unconfigured, so this is
        # a no-op in lagh-alone mode and calls the LLM only when machine/.env is set.
        out = run_sync(pid, oracle, box, propose_fn=propose_forms)
        law = out.get("law", "")
        if out["outcome"] == "proved" and law:
            sc = score(pid, law)
            cat = "PARSE_ERROR" if sc is None else ("CORRECT" if sc[0] else "CONFIDENT_WRONG")
        else:
            cat = "HONEST_ABSTAIN"
        rows.append({"pid": pid, "cat": cat, "outcome": out["outcome"], "law": law,
                     "class": out.get("characterization", ""), "lagh_alone": base.get(pid),
                     "secs": round(time.time() - t0, 1)})
        tail = law if law else f"[{out.get('characterization', '')}]"
        print(f"  [{k:>3}/{len(ids)}] {pid:<28} {cat:<16} {out['outcome']:<9} {tail[:34]}",
              flush=True)

    c = Counter(r["cat"] for r in rows)
    base_n = sum(1 for r in rows if r["lagh_alone"])
    gains = [r["pid"] for r in rows if r["cat"] == "CORRECT" and not r["lagh_alone"]]
    losses = [r["pid"] for r in rows if r["cat"] in ("CONFIDENT_WRONG",) and r["lagh_alone"]]
    print("\n" + "=" * 52 + f"\nSCORECARD ({len(rows)} cells)\n" + "=" * 52)
    for kk in ("CORRECT", "HONEST_ABSTAIN", "CONFIDENT_WRONG", "PARSE_ERROR"):
        print(f"  {kk:<16} {c[kk]}")
    print(f"\n  composite (machine+lagh) certified-correct: {c['CORRECT']}/{len(rows)}")
    print(f"  lagh-alone baseline on same cells:          {base_n}/{len(rows)}")
    print(f"  gains (composite got, baseline missed):     {len(gains)}  {gains}")
    print(f"\n  CONFIDENT-WRONG must be 0 -> {'OK' if c['CONFIDENT_WRONG'] == 0 else 'INVARIANT BREACH'}")
    json.dump(rows, open(os.path.join(os.path.dirname(__file__), "run_results.json"), "w"), indent=1)


if __name__ == "__main__":
    main()

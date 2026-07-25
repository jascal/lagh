"""Grammar-only re-probe of the LLM-verified benchmark cells (no LLM anywhere):
did CAP-S / CAP-R / CAP-B-lift make the grammar reach what only the proposer
could reach before? Scored per the registered predictions in LLMSRBENCH_DEV.md.
"""
from __future__ import annotations

import json
import signal
import sys
import time
from pathlib import Path

import numpy as np
import sympy as sp

sys.path.insert(0, str(Path(__file__).parent.parent))
sys.path.insert(0, "/home/allans/code/llm-srbench")

from experiments.run_dev_llmsrbench import (SIGMA_REP, SUB_CAP, _gt_expr_dev)  # noqa: E402
from experiments.rejudge_dev_llmsrbench import equivalent  # noqa: E402


def _h(s, f):
    raise TimeoutError()


def main():
    from experiments.run_blind_llmsrbench import load_all_problems, _train_xy
    from lagh.passive import discover_passive
    d = json.load(open("experiments/results/dev_llmsrbench_v1_scores_v2.json"))
    targets = [r for r in d["rows"]
               if r["track"] == "certified" and r["channel"] == "llm-verified"]
    problems = {f"{m}/{p.equation_idx}": p for m, p in load_all_problems()}
    signal.signal(signal.SIGALRM, _h)
    got = 0
    for r in sorted(targets, key=lambda x: x["id"]):
        p = problems[r["id"]]
        X, y = _train_xy(p)
        if len(X) > SUB_CAP:
            idx = np.sort(np.random.default_rng(0).choice(len(X), SUB_CAP,
                                                          replace=False))
            X, y = X[idx], y[idx]
        t0 = time.time()
        law = None
        try:
            signal.alarm(350)
            pr = discover_passive(X, y, sigma=SIGMA_REP)
            if pr.certified:
                law = str(pr.result.expr)
        except TimeoutError:
            pass
        finally:
            signal.alarm(0)
        verdict = "-"
        if law:
            ok, how = equivalent(law, _gt_expr_dev(p))
            verdict = f"CERT/{'SA-ok' if ok else 'SA-WRONG'}"
            got += ok
        print(f"{r['id']:45s} {verdict:14s} {time.time()-t0:5.0f}s  "
              f"{(law or '')[:60]}", flush=True)
    print(f"\ngrammar-only SA-correct: {got}/{len(targets)} "
          f"(were 19/19 with the proposer)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

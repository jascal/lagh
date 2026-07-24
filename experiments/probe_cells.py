"""Targeted cell probe for capability scoring (gap-plan ground rule: isolated verify
on the predicted cells + regression cells BEFORE any full re-sweep).
Usage: probe_cells.py module/diff/version [...]"""
from __future__ import annotations
import sys, time
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent.parent))
from lagh.acquisition import run_active
from lagh.adapters.newtonbench import MODULES, make_oracle
from experiments.run_newtonbench_all import dense_ok


def main(cells):
    for cell in cells:
        module, diff, v = cell.split("/")
        inputs, lo, hi = MODULES[module]
        oracle = make_oracle(module, v, diff)
        t0 = time.time()
        r = run_active(oracle, lo, hi, seed=1, time_budget_s=150)
        correct = dense_ok(r.result.expr, oracle, lo, hi, len(inputs))
        mk = ("ok" if (r.result.certificate.certified and correct) else
              "WRONG" if (r.result.certificate.certified and correct is False) else
              f"-[{r.result.certificate.abstain}]")
        law = str(r.result.expr)[:70] if r.result.expr is not None else ""
        print(f"{cell:34s} {mk:14s} {time.time()-t0:6.1f}s  {law}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))

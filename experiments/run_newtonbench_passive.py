"""NewtonBench-dev PASSIVE sweep (docs/DIRECTION_PASSIVE.md): one fixed dataset per
cell, no oracle at discovery time, both sampling laws. Scored by the same dense grid
as the active sweep so the two regimes are directly comparable.
DEV MEASUREMENT (STRATEGY.md) -- no win claim."""
from __future__ import annotations
import json, sys, time
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent.parent))
import numpy as np
from lagh.passive import discover_passive
from lagh.adapters.newtonbench import MODULES, available_versions, make_oracle
from experiments.run_newtonbench_all import dense_ok

N_POINTS = 250


def make_dataset(oracle, lo, hi, dim, sampling, seed):
    rng = np.random.default_rng(seed)
    lo, hi = np.asarray(lo, float), np.asarray(hi, float)
    if sampling == "loguniform":
        X = np.exp(rng.uniform(np.log(lo), np.log(hi), (N_POINTS, dim)))
    else:
        X = rng.uniform(lo, hi, (N_POINTS, dim))
    return X, oracle(X)


def main():
    out = Path("experiments/results/newtonbench_passive.jsonl")
    out.parent.mkdir(parents=True, exist_ok=True)
    rows = []
    for sampling in ("loguniform", "uniform"):
        for mi, (module, (inputs, lo, hi)) in enumerate(MODULES.items()):
            dim = len(inputs)
            for di, diff in enumerate(["easy", "medium", "hard"]):
                for vi, v in enumerate(available_versions(module, diff)):
                    oracle = make_oracle(module, v, diff)
                    seed = 100000 + 1000 * mi + 100 * di + vi   # fixed per cell
                    X, y = make_dataset(oracle, lo, hi, dim, sampling, seed)
                    t0 = time.time()
                    r = discover_passive(X, y)
                    expr = r.result.expr if r.certified else None
                    correct = dense_ok(expr, oracle, lo, hi, dim)
                    rec = {"sampling": sampling, "module": module, "difficulty": diff,
                           "version": v, "certified": r.certified, "correct": correct,
                           "abstain": r.result.certificate.abstain,
                           "resplits": r.resplits_tried,
                           "law": str(expr) if expr is not None else None,
                           "confident_wrong": bool(r.certified and correct is False),
                           "seconds": round(time.time() - t0, 1)}
                    rows.append(rec)
                    mk = "ok" if (rec["certified"] and correct) else \
                        ("WRONG" if rec["confident_wrong"] else "-")
                    print(f"{sampling:10s} {module:22s} {diff:6s} {v} {mk:5s} "
                          f"{rec['seconds']}s", flush=True)
        out.write_text("\n".join(json.dumps(r) for r in rows) + "\n")

    # score against the ACTIVE sweep (predictions P1-P4)
    act_path = Path("experiments/results/newtonbench_all.jsonl")
    active_ok = set()
    if act_path.exists():
        for line in act_path.read_text().splitlines():
            a = json.loads(line)
            if a["certified"] and a["correct"]:
                active_ok.add((a["module"], a["difficulty"], a["version"]))
    for sampling in ("loguniform", "uniform"):
        sub = [r for r in rows if r["sampling"] == sampling]
        ok = {(r["module"], r["difficulty"], r["version"])
              for r in sub if r["certified"] and r["correct"]}
        cw = sum(r["confident_wrong"] for r in sub)
        print(f"\n{sampling}: {len(ok)}/108 correct, confident-wrong {cw} (must be 0)")
        lost = sorted(active_ok - ok)
        new = sorted(ok - active_ok)
        print(f"  lost vs active ({len(lost)}): {['/'.join(t) for t in lost]}")
        print(f"  new vs active ({len(new)}): {['/'.join(t) for t in new]}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

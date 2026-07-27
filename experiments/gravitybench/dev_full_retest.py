"""DEV re-measure of the full-obs variant with the native-cadence epoch fix.
The read numbers (63.59%) STAND; this is registered dev work on the diagnosed
mechanism (BLIND_READ_REPORT_GRAVITYBENCH.md)."""
import io, json, signal, sys
from pathlib import Path
import numpy as np, pandas as pd
sys.path.insert(0, str(Path(__file__).parent.parent.parent))
from experiments.gravitybench.run_blind_read import (detect_units, load_instances,
                                                     solve)

def h(s, f): raise TimeoutError()

def main():
    insts = load_instances()
    signal.signal(signal.SIGALRM, h)
    out = Path("experiments/results/gravitybench_dev_full_fix.jsonl")
    done = set()
    if out.exists():
        done = {json.loads(l)["scenario_id"] for l in out.read_text().splitlines()}
    with out.open("a") as fh:
        for inst in insts:
            if inst["scenario_id"] in done: continue
            df = pd.read_csv(io.StringIO(inst["simulation_csv_content"]))
            rec = {"scenario_id": inst["scenario_id"], "task": inst["scenario_name"]}
            try:
                signal.alarm(300)
                ans, val, n_used, state = solve(df, inst["scenario_name"],
                                                detect_units(inst), "full",
                                                variant_epoch_fix=True)
                rec["answer"] = (bool(ans) if isinstance(ans, (bool, np.bool_))
                                 else (float(ans) if ans is not None else None))
            except Exception as e:
                rec["answer"] = None; rec["error"] = str(e)[:100]
            finally:
                signal.alarm(0)
            truth = inst["true_answer"]; thr = inst.get("full_obs_threshold_percent")
            ok = False
            if rec["answer"] is not None and truth is not None:
                if str(truth) in ("True", "False") or isinstance(truth, bool):
                    tb = truth if isinstance(truth, bool) else (str(truth) == "True")
                    ok = bool(rec["answer"]) == tb
                else:
                    tv = float(truth)
                    err = abs(float(rec["answer"]) - tv) / (abs(tv) + 1e-300) * 100
                    ok = (err <= float(thr)) if thr is not None else (err <= 10.0)
            rec["correct"] = ok
            fh.write(json.dumps(rec) + "\n"); fh.flush()
    rows = [json.loads(l) for l in out.read_text().splitlines()]
    c = sum(r["correct"] for r in rows)
    print(f"DEV full-variant with epoch fix: {c}/{len(rows)} = {100*c/len(rows):.2f}% "
          f"(read stood at 63.59%; SOTA 74%)")

if __name__ == "__main__":
    main()

"""THE GRAVITY-BENCH ONE-SHOT READ (docs/BLIND_READ_REGISTRATION_GRAVITYBENCH.md).

206 validation instances x two variants:
- budgeted: the fixed observation policy over an Observe(times) interface
  backed by the instance's simulation table (cubic interpolation, budget 100,
  10 per request -- their row-wise contract);
- full: the entire simulation table handed to system_id (no planner).

Answers in SI. Scoring: per-instance percent thresholds from the dataset
(budget_obs_threshold_percent / full_obs_threshold_percent) against
true_answer. Incremental journal; crash-fixes allowed and logged;
scoring-affecting changes prohibited (registration section 5).
"""
from __future__ import annotations

import io
import json
import signal
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.interpolate import CubicSpline

sys.path.insert(0, str(Path(__file__).parent.parent.parent))
from experiments.gravitybench import astronomer as ast  # noqa: E402
from experiments.gravitybench.driver import obs_to_si, unit_factors  # noqa: E402
from experiments.gravitybench.twin import Twin, system_id  # noqa: E402

OUT = Path("experiments/results/gravitybench_read.jsonl")
# CRASH-FIX (logged, registration 5): the repo's validation.jsonl carries
# placeholder CSV content; the real data is the HF dataset. All 412 first-pass
# solves crashed on the placeholder (zero scoring information revealed).
# BOUNDED-RUNTIME AMENDMENT (logged): sims are ~1e6 rows; the FULL variant
# subsamples uniformly to 2400 rows (uniform cadence preserved -- their own
# expert baseline used 100 uniform samples).
FULL_SUBSAMPLE = 2400

# task_prompt -> units detection: instances carry expected_units; the sim table
# columns are in the scenario's native units. Column names follow Binary.py.


def _timeout(s, f):
    raise TimeoutError()


def load_instances():
    from datasets import load_dataset
    ds = load_dataset("GravityBench/GravityBench")
    split = list(ds.keys())[0]
    return list(ds[split])


def detect_units(inst):
    """Native units from the filename/variation conventions: the HF set stores
    sims in SI unless the variation says otherwise."""
    v = (inst.get("variation_name") or "") + " " + (inst.get("simulation_csv_filename") or "")
    v = v.lower()
    if "cgs" in v:
        return ("s", "cm", "g")
    if "yraumsun" in v or "astronomical" in v:
        return ("yr", "AU", "Msun")
    return ("m", "s", "kg")


def solve(df, task, units, variant):
    lf, tf = unit_factors(units)
    t_nat = df["time"].to_numpy(float)
    maxtime_nat = float(t_nat.max())
    cols = [c for c in df.columns if c != "time"]
    splines = {c: CubicSpline(t_nat, df[c].to_numpy(float)) for c in cols}

    if variant == "full":
        step = max(1, len(df) // FULL_SUBSAMPLE)
        sub = df.iloc[::step]
        obs = obs_to_si({c: sub[c].to_numpy(float) for c in sub.columns}, lf, tf)
        n_used = len(sub)
    else:
        used = {"n": 0}

        def observe_si(times_si):
            tn = np.clip(np.asarray(times_si, float) / tf, 0, maxtime_nat)
            used["n"] += len(tn)
            out = {"time": tn}
            for c in cols:
                out[c] = splines[c](tn)
            return obs_to_si(out, lf, tf)

        obs, _ = ast.plan_and_observe(observe_si, maxtime_nat * tf, budget=100,
                                      per_request=10)
        n_used = used["n"]
    state = system_id(obs)
    tw = Twin(state, float(obs["time"].max()))
    ans = tw.answer(task)
    return ans, tw.validate(obs), n_used, state


def main():
    insts = load_instances()
    done = set()
    if OUT.exists():
        for line in OUT.read_text().splitlines():
            r = json.loads(line)
            done.add((r["scenario_id"], r["variant"]))
    signal.signal(signal.SIGALRM, _timeout)
    with OUT.open("a") as fh:
        for inst in insts:
            df = pd.read_csv(io.StringIO(inst["simulation_csv_content"]))
            task = inst["scenario_name"]
            units = detect_units(inst)
            for variant in ("budget", "full"):
                key = (inst["scenario_id"], variant)
                if key in done:
                    continue
                t0 = time.time()
                rec = {"scenario_id": inst["scenario_id"], "task": task,
                       "variation": inst.get("variation_name"),
                       "variant": variant, "units": list(units)}
                try:
                    signal.alarm(300)
                    ans, val, n_used, state = solve(df, task, units, variant)
                    rec.update(answer=(bool(ans) if isinstance(ans, (bool, np.bool_))
                                       else (float(ans) if ans is not None else None)),
                               twin_validation=float(val), n_obs=n_used,
                               alpha_fit=state.get("p_raw"))
                except TimeoutError:
                    rec.update(answer=None, error="timeout")
                except Exception as e:                          # noqa: BLE001
                    rec.update(answer=None, error=str(e)[:120])
                finally:
                    signal.alarm(0)
                rec["secs"] = round(time.time() - t0, 1)
                fh.write(json.dumps(rec) + "\n")
                fh.flush()
                print(f"{inst['scenario_id']:38s} {variant:6s} "
                      f"{'ok' if rec.get('answer') is not None else rec.get('error','-'):18s}"
                      f" {rec['secs']}s", flush=True)

    # ---- scoring (mechanical, per-instance thresholds from the dataset) ----
    rows = [json.loads(line) for line in OUT.read_text().splitlines()]
    inst_by_id = {i["scenario_id"]: i for i in insts}
    scored = []
    for r in rows:
        inst = inst_by_id[r["scenario_id"]]
        truth = inst["true_answer"]
        thr_key = ("budget_obs_threshold_percent" if r["variant"] == "budget"
                   else "full_obs_threshold_percent")
        thr = inst.get(thr_key)
        ok = False
        if r.get("answer") is not None and truth is not None:
            if isinstance(truth, bool) or str(truth) in ("True", "False"):
                tb = truth if isinstance(truth, bool) else (str(truth) == "True")
                ok = bool(r["answer"]) == tb
            else:
                tv = float(truth)
                err = abs(float(r["answer"]) - tv) / (abs(tv) + 1e-300) * 100.0
                ok = (err <= float(thr)) if thr is not None else (err <= 10.0)
        scored.append({**r, "correct": ok})
    summary = {}
    for variant in ("budget", "full"):
        sub = [r for r in scored if r["variant"] == variant]
        summary[variant] = {"n": len(sub),
                            "correct": sum(r["correct"] for r in sub),
                            "pct": round(100 * sum(r["correct"] for r in sub)
                                         / max(len(sub), 1), 2)}
    Path("experiments/results/gravitybench_read_scores.json").write_text(
        json.dumps({"summary": summary, "rows": scored}, indent=1))
    print(json.dumps(summary, indent=1))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

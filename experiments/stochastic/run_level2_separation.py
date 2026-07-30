"""Level 2, the separation increment: process noise vs measurement noise, MEASURED.

`docs/DIRECTION_ERROR_PROVENANCE.md`'s central question, and the reason Levels 0-2
were justified independently of the benchmark framing: that question had ground truth
NOWHERE in this repo. This run constructs it -- systems where sigma_obs and b are both
known by construction -- and scores the recovery of both.

The mechanism is one polynomial fit in the stride s (`weakform.qv_three_way`):

    sum_i (u[i+s] - u[i])^2  ~=  c + alpha*s + beta*s^2
                                 |    |         `- smooth (differentiable) part
                                 |    `- PROCESS noise: a martingale increment's
                                 |       variance grows linearly in the lag
                                 `- OBSERVATION noise: iid, so E[(e_{i+s}-e_i)^2]
                                    = 2 sigma^2 at EVERY lag

Three sources, three distinct exponents. So sigma_obs is measured rather than
declared, which retires at its root the declaration that produced this arc's only
confident-wrong.

Run: .venv/bin/python experiments/stochastic/run_level2_separation.py
"""
from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).parent.parent.parent))

from experiments.stochastic.generator import (double_well_paths,  # noqa: E402
                                              ou_paths)
from lagh.weakform import QV_STRIDES_3, qv_three_way  # noqa: E402

OUT = Path("experiments/results/stochastic_level2_separation.json")


def sums(u, strides=QV_STRIDES_3):
    return [float(np.sum((u[s:] - u[:-s]) ** 2)) for s in strides]


def measure(u, T):
    """(sigma_obs, implied b^2, info) for one observed path."""
    sc, info = qv_three_way(sums(u), n_increments=len(u) - 1)
    if sc is None:
        return None, None, info
    return info.get("sigma_obs"), sc / T, info


def main():
    t0 = time.time()
    rows = []
    rng = np.random.default_rng(21)

    # ---- a grid over BOTH noise scales, on two process families. sigma_obs = 0 and
    # b = 0 are included because a null in each direction is what makes a rate
    # meaningful -- the same discipline every lagh campaign carries.
    for family in ("ou", "double_well", "deterministic"):
        for b in (0.0, 0.2, 0.7, 1.4):
            for sob in (0.0, 1e-3, 1e-2, 5e-2):
                T, dt, n_traj = 200.0, 1e-3, 1
                if family == "ou":
                    if b == 0.0:
                        continue                       # OU with b = 0 is a fixed point
                    t, X = ou_paths(theta=1.0, b=b, T=T, dt=dt, n_traj=n_traj,
                                    seed=int(1000 * b) + int(1e4 * sob), x0=0.0)
                elif family == "double_well":
                    if b == 0.0:
                        continue
                    t, X = double_well_paths(theta=1.0, b=b, T=T, dt=dt,
                                             n_traj=n_traj, substeps=4,
                                             seed=int(2000 * b) + int(1e4 * sob))
                else:
                    if b != 0.0:
                        continue                       # the NO-process-noise family
                    n = int(round(T / dt)) + 1
                    t = np.arange(n) * dt
                    X = (np.exp(-t) + 0.3 * np.sin(2 * t))[None, :]
                u = X[0] + (sob * rng.standard_normal(X.shape[1]) if sob else 0.0)
                s_hat, b2_hat, info = measure(u, T)
                b2_true = b ** 2
                rows.append({
                    "family": family, "b_true": b, "sigma_obs_true": sob,
                    "sigma_obs_measured": s_hat,
                    "sigma_obs_rel_err": (None if not sob else
                                          abs(s_hat - sob) / sob),
                    "b2_true": b2_true, "b2_measured": b2_hat,
                    "b2_rel_err": (None if b2_true == 0 else
                                   abs(b2_hat - b2_true) / b2_true),
                    "dominant": info.get("dominant"),
                    "separable": info.get("separable"),
                    "process_over_observation": info.get(
                        "process_over_observation"),
                    "fit_residual_rel": info.get("fit_residual_rel"),
                })

    # ---- scoring: how well is each source recovered, and is the verdict right?
    # the headline claim is scoped to the SEPARABLE regime, because outside it the
    # measurement is not merely imprecise but biased in the dangerous direction
    with_obs = [r for r in rows if r["sigma_obs_true"] > 0 and r["separable"]]
    buried = [r for r in rows if r["sigma_obs_true"] > 0 and not r["separable"]]
    with_proc = [r for r in rows if r["b2_true"] > 0]
    no_obs = [r for r in rows if r["sigma_obs_true"] == 0]
    no_proc = [r for r in rows if r["b2_true"] == 0]
    summary = {
        "n_cases": len(rows),
        "sigma_obs_separable_regime": {
            "n": len(with_obs),
            "median_rel_err": float(np.median([r["sigma_obs_rel_err"]
                                               for r in with_obs])),
            "worst_rel_err": float(np.max([r["sigma_obs_rel_err"]
                                           for r in with_obs])),
            "within_10pct": int(sum(r["sigma_obs_rel_err"] < 0.10
                                    for r in with_obs)),
        },
        "sigma_obs_buried_refused": {
            "n": len(buried),
            "worst_rel_err_had_it_been_used": float(np.max(
                [r["sigma_obs_rel_err"] for r in buried])) if buried else None,
            "note": ("refused by the c >= 0.1*alpha bar; these are the cases where "
                     "using the measurement would have over-declared sigma_obs and "
                     "tightened a band"),
        },
        "process_b2": {
            "n": len(with_proc),
            "median_rel_err": float(np.median([r["b2_rel_err"]
                                               for r in with_proc])),
            "worst_rel_err": float(np.max([r["b2_rel_err"] for r in with_proc])),
            "within_10pct": int(sum(r["b2_rel_err"] < 0.10 for r in with_proc)),
        },
        # THE NULLS, in both directions: a rate is meaningless without them
        "null_no_observation_noise": {
            "n": len(no_obs),
            "worst_spurious_sigma_obs": float(np.max([r["sigma_obs_measured"]
                                                      for r in no_obs])),
        },
        "null_no_process_noise": {
            "n": len(no_proc),
            "worst_spurious_b2": float(np.max([r["b2_measured"]
                                               for r in no_proc])),
        },
    }
    res = {"level": 2, "increment": "process vs measurement separation",
           "seconds": round(time.time() - t0, 1),
           "strides": list(QV_STRIDES_3), "cases": rows, "summary": summary}
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(res, indent=1, default=str))
    print(f"wrote {OUT} ({res['seconds']}s)\n")
    sg = summary["sigma_obs_separable_regime"]
    print("sigma_obs recovery, SEPARABLE regime (c >= 0.1*alpha):")
    print("   n=%d  median rel err %.4f  worst %.4f  within 10%%: %d/%d"
          % (sg["n"], sg["median_rel_err"], sg["worst_rel_err"],
             sg["within_10pct"], sg["n"]))
    bu = summary["sigma_obs_buried_refused"]
    print("   REFUSED as buried: %d cases (worst error had it been used: %.2f)"
          % (bu["n"], bu["worst_rel_err_had_it_been_used"] or 0.0))
    print("process b^2 recovery (cases with process noise):")
    print("   n=%d  median rel err %.4f  worst %.4f  within 10%%: %d/%d"
          % (summary["process_b2"]["n"], summary["process_b2"]["median_rel_err"],
             summary["process_b2"]["worst_rel_err"],
             summary["process_b2"]["within_10pct"], summary["process_b2"]["n"]))
    print("NULLS:")
    print("   no observation noise -> worst spurious sigma_obs = %.3e (%d cases)"
          % (summary["null_no_observation_noise"]["worst_spurious_sigma_obs"],
             summary["null_no_observation_noise"]["n"]))
    print("   no process noise     -> worst spurious b^2       = %.3e (%d cases)"
          % (summary["null_no_process_noise"]["worst_spurious_b2"],
             summary["null_no_process_noise"]["n"]))


if __name__ == "__main__":
    main()

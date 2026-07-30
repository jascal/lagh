"""H2a driver: speaks the benchmark's contract, converts units at the boundary,
runs the fixed observation policy + digital twin, returns the answer.

All physics runs in SI internally; observations arrive in the scenario's native
units (AU/yr, CGS, or SI) and are converted ONCE at the boundary. Answers return
in SI (the benchmark validates in SI after its own conversion). The observation
budget is respected exactly; the policy is the registered fixed one.
"""
from __future__ import annotations

import numpy as np

from experiments.gravitybench import astronomer as ast
from experiments.gravitybench.twin import Twin, system_id

AU_M = 1.495978707e11
YR_S = 3.15576e7          # Julian year
LENGTH = {"m": 1.0, "AU": AU_M, "cm": 1e-2}
TIME = {"s": 1.0, "yr": YR_S}
MASS = {"kg": 1.0, "Msun": 1.98892e30, "g": 1e-3}


def unit_factors(units):
    """units: iterable of unit names in any order (their tuples vary)."""
    lf = tf = None
    for u in units:
        if u in LENGTH:
            lf = LENGTH[u]
        if u in TIME:
            tf = TIME[u]
    if lf is None or tf is None:
        raise ValueError(f"unrecognized unit system {units!r}")
    return lf, tf


def obs_to_si(obs, lf, tf):
    out = {"time": np.asarray(obs["time"], float) * tf}
    for k, v in obs.items():
        if k != "time":
            out[k] = np.asarray(v, float) * lf
    return out


def solve_instance(observe_native, maxtime_native, task, units,
                   budget=100, per_request=10):
    """observe_native(times_native)->columns dict in native units. Returns
    {answer (SI), state, twin_validation, n_obs_used}."""
    lf, tf = unit_factors(units)

    used = {"n": 0}

    def observe_si(times_si):
        times_nat = np.asarray(times_si, float) / tf
        times_nat = np.clip(times_nat, 0.0, maxtime_native)
        out_frames = []
        for i in range(0, len(times_nat), per_request):
            chunk = times_nat[i:i + per_request]
            used["n"] += len(chunk)
            out_frames.append(observe_native(chunk))
        merged = {k: np.concatenate([np.asarray(f[k], float) for f in out_frames])
                  for k in out_frames[0]}
        return obs_to_si(merged, lf, tf)

    maxtime_si = maxtime_native * tf
    obs, P0 = ast.plan_and_observe(observe_si, maxtime_si, budget=budget,
                                   per_request=per_request)
    state = system_id(obs)
    tw = Twin(state, maxtime_si)
    ans, val, refusal = tw.gated_answer(task, obs)
    return {"answer": ans, "state": state, "twin_validation": val,
            "refusal": refusal, "n_obs_used": used["n"], "period_est": P0}

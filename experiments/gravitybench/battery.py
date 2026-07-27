"""H2a robustness battery: the astronomer across regimes, DEV orbits only.

Sweeps eccentricity x mass-ratio x window coverage, plus physics variants
(modified gravity, drag, unbound flyby, proper motion) and unit systems
(SI / AU-yr / CGS through the driver's boundary conversion). Scores the CORE
quantities against analytic truth + the twin's own prediction error. The
output table drives planner/estimator fixes BEFORE the sealed read.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).parent.parent.parent))
from experiments.gravitybench.integrator import G_SI, TwoBody, make_circularish  # noqa: E402
from experiments.gravitybench.driver import AU_M, YR_S, solve_instance  # noqa: E402

OUT = Path("experiments/results/gravitybench_battery.jsonl")


def native_observe(tb, lf, tf):
    def fn(times_native):
        return {k: (np.asarray(v) / (tf if k == "time" else lf))
                for k, v in tb.observe(np.asarray(times_native, float) * tf).items()}
    return fn


def run_case(name, tb, P_ref, maxtime_si, checks, units=("m", "s", "kg"),
             lf=1.0, tf=1.0, budget=100):
    tb.run(maxtime_si)
    obs_fn = native_observe(tb, lf, tf)
    row = {"case": name, "budget": budget}
    try:
        first_task = next(iter(checks))
        res = solve_instance(obs_fn, maxtime_si / tf, first_task, units,
                             budget=budget)
        row["twin_validation"] = round(res["twin_validation"], 4)
        row["n_obs"] = res["n_obs_used"]
        from experiments.gravitybench.twin import Twin
        tw = Twin(res["state"], maxtime_si)
        errs = {}
        for task, truth in checks.items():
            got = tw.answer(task)
            if isinstance(truth, bool):
                errs[task] = 0.0 if got == truth else 1.0
            elif got is None or truth is None:
                errs[task] = None if truth is None else 1.0
            elif task == "eccentricity":
                # small-e: score absolute error against a 0.1 floor (their
                # thresholds are percent-based on healthy-e instances)
                errs[task] = round(abs(got - truth) / max(abs(truth), 0.1), 4)
            else:
                errs[task] = round(abs(got - truth) / (abs(truth) + 1e-300), 4)
        row["errors"] = errs
        row["worst"] = max([e for e in errs.values() if e is not None] or [None])
    except Exception as e:                                     # noqa: BLE001
        row["error"] = str(e)[:120]
        row["worst"] = 1.0
    return row


def main():
    rows = []
    # sweep: eccentricity x mass ratio x window
    for ecc in (0.05, 0.3, 0.6, 0.85):
        for q in (1.0, 3.0, 10.0):
            for wins in (1.5, 4.0, 12.0):
                m2 = 1e30
                m1 = q * m2
                a = 1.5e11
                tb, P = make_circularish(m1=m1, m2=m2, a=a, ecc=ecc)
                checks = {"period": P, "total_mass": m1 + m2,
                          "mass_star1": m1, "eccentricity": ecc,
                          "semi_major_axis": a,
                          "apoastron": a * (1 + ecc), "periastron": a * (1 - ecc),
                          "is_bound": True, "kepler_3rd_law": True}
                rows.append(run_case(f"e{ecc}_q{q}_w{wins}", tb, P, wins * P,
                                     checks))
                print(f"{rows[-1]['case']:22s} worst={rows[-1].get('worst')}",
                      flush=True)
    # physics variants
    tb, P = make_circularish(ecc=0.35, mod_gravity_exponent=-2.4)
    rows.append(run_case("modgrav_a0.4", tb, P, 4 * P,
                         {"modified_gravity_power_law": 0.4,
                          "kepler_3rd_law": False}))
    tb, P = make_circularish(ecc=0.35, mod_gravity_exponent=-1.7)
    rows.append(run_case("modgrav_a-0.3", tb, P, 4 * P,
                         {"modified_gravity_power_law": -0.3,
                          "kepler_3rd_law": False}))
    tb, P = make_circularish(ecc=0.3, drag_tau=8 * 6.1e7)
    rows.append(run_case("drag", tb, P, 4 * P, {"linear_drag": 8 * 6.1e7}))
    # unbound flyby: hyperbolic (v > v_escape at start)
    m1, m2 = 2e30, 1e30
    M = m1 + m2
    r0 = 3e11
    v_esc = np.sqrt(2 * G_SI * M / r0)
    # COM-consistent construction so |v2 - v1| is EXACTLY 1.3 x escape speed
    # (two previous attempts hand-built components and accidentally made bound
    # systems -- the astronomer's triplet energies exposed both test bugs)
    vdir = np.array([0.3, -0.9, 0.3])
    vrel_vec = 1.3 * v_esc * vdir / np.linalg.norm(vdir)
    v_com = np.array([1e3, 0.0, 0.0])
    tbu = TwoBody(m1, m2, [-m2 / M * r0, 0, 0], [m1 / M * r0, 0, 0],
                  v_com - m2 / M * vrel_vec, v_com + m1 / M * vrel_vec)
    rows.append(run_case("unbound", tbu, None, 2e8,
                         {"is_bound": False, "total_mass": M, "mass_ratio": 2.0}))
    # proper motion: drifting COM
    tbd, P = make_circularish(ecc=0.3)
    for v in (tbd.v1, tbd.v2):
        v += np.array([3e3, -1e3, 5e2])
    rows.append(run_case("proper_motion", tbd, P, 4 * P,
                         {"period": P, "mass_ratio": 2.0, "total_mass": 3e30}))
    # unit systems through the driver boundary
    tb, P = make_circularish(ecc=0.4)
    rows.append(run_case("units_AUyr", tb, P, 4 * P,
                         {"period": P, "total_mass": 3e30},
                         units=("yr", "AU", "Msun"), lf=AU_M, tf=YR_S))
    tb, P = make_circularish(ecc=0.4)
    rows.append(run_case("units_CGS", tb, P, 4 * P,
                         {"period": P, "total_mass": 3e30},
                         units=("s", "cm", "g"), lf=1e-2, tf=1.0))
    for r in rows[-6:]:
        print(f"{r['case']:22s} worst={r.get('worst')} "
              f"{r.get('error','')}", flush=True)

    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text("\n".join(json.dumps(r) for r in rows) + "\n")
    bad = [r for r in rows if r.get("worst") is None or r["worst"] > 0.15]
    print(f"\nBATTERY: {len(rows)} cases, {len(bad)} above 15% worst-error")
    for r in bad:
        print("  FAIL", r["case"], r.get("errors", r.get("error")))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

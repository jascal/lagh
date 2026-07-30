"""Level 0 of the stochastic suite: calibrate the coverage factor (step 2 of
docs/DIRECTION_STOCHASTIC.md).

Six sections, in the order the registered predictions need them:

  A  kappa CALIBRATION -- the deliverable. Replicate runs at the TRUE law: does
     the exhaustive martingale band hold on every row at the claimed budget?
  B  the reach boundary -- where a drift is vacuous, partially determined, or
     certified, against the arithmetic thresholds `lagh/ito.py` derives.
  C  S1 -- how the interval width scales in the ROW COUNT and in the WINDOW
     LENGTH. These are different answers and the prediction assumed one.
  D  the five tasks, scored through the FROZEN checker (S4, S5).
  E  S6 -- realized quadratic variation on a system with no process noise.
  F  S7 -- alpha against the disjoint window count.

Run: .venv/bin/python experiments/stochastic/run_level0.py
"""
from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).parent.parent.parent))

from experiments.stochastic.generator import (LIBRARY, level0_tasks,  # noqa: E402
                                              ou_paths, paths_for)
from lagh.certify import admissible_interval, coverage_factor  # noqa: E402
from lagh.ito import LAM_QV, ItoBand, build_rows, certify_drift  # noqa: E402
from lagh.stochcheck import (Consumer, Coverage, Declaration,  # noqa: E402
                             Submission, component, rank, score_task,
                             suite_totals)

OUT = Path("experiments/results/stochastic_level0.json")


def _band_of(rows, delta=0.05, sigma_obs=0.0):
    return ItoBand(rows.qv, rows.qv_se, rows.corr_se, rows.quad, len(rows.y),
                   delta=delta, sigma_obs=sigma_obs, y=rows.y,
                   feat_names=list(rows.names))


def _lp(rows, delta=0.05, coeff_max=10.0):
    band = _band_of(rows, delta)
    e0 = band.martingale() + LAM_QV * rows.corr_se + rows.quad[:, 0]
    q = rows.quad[:, 1:].sum(axis=1)
    return admissible_interval(rows.A, rows.y, lambda c: e0 + c * q,
                               coeff_max=coeff_max)


# ---------------------------------------------------------- A: kappa calibration

def section_a(*, n_rep: int = 40, deltas=(0.5, 0.2, 0.05, 0.01)) -> dict:
    """Does the exhaustive band hold at its claimed budget?

    Per replicate: assemble rows from independent seeds, evaluate the residual of
    the TRUE law on every row, and record whether ANY row falls outside its band.
    The empirical failure rate is then compared with delta -- the direct
    calibration the sequencing asks for, on discovery rows rather than in
    simulation.
    """
    out = {}
    for name, mk, truth in (
        ("ou", lambda s: ou_paths(theta=1.0, b=1.4, T=640.0, dt=1e-3,
                                  n_traj=2, seed=1000 + s, x0=0.0), -1.0),
    ):
        rows_cache = []
        for s in range(n_rep):
            t, X = mk(s)
            rows_cache.append(build_rows(t, X, LIBRARY, half=32000))
        per = {}
        for delta in deltas:
            fails, worst = 0, []
            for rows in rows_cache:
                band = _band_of(rows, delta)
                e = band(None)
                j = LIBRARY.index("x")
                resid = rows.y - truth * rows.A[:, j]
                ratio = float(np.max(np.abs(resid) / e))
                worst.append(ratio)
                fails += int(ratio > 1.0)
            per[str(delta)] = {
                "kappa": coverage_factor(len(rows_cache[0].y), delta),
                "n_replicates": n_rep,
                "failures": fails,
                "empirical_rate": fails / n_rep,
                "declared_delta": delta,
                "holds": bool(fails / n_rep <= delta),
                "max_resid_over_band_median": float(np.median(worst)),
                "max_resid_over_band_max": float(np.max(worst)),
            }
        out[name] = {"n_rows": int(len(rows_cache[0].y)), "by_delta": per}
    return out


# ------------------------------------------------------- B: the reach boundary

def section_b() -> dict:
    """OU against theta*L > 2 kappa^2, and GBM against the diffusion intensity."""
    ou = []
    for half, T in ((2000, 640.0), (8000, 640.0), (32000, 640.0),
                    (64000, 1280.0)):
        t, X = ou_paths(theta=1.0, b=1.4, T=T, dt=1e-3, n_traj=8, seed=0,
                        x0=0.0)
        rows = build_rows(t, X, LIBRARY, half=half)
        if len(rows.y) < 6:
            ou.append({"L": 2 * half * 1e-3, "rows": int(len(rows.y)),
                       "verdict": "too few rows"})
            continue
        band = _band_of(rows)
        r = certify_drift(rows, delta=0.05, seed=0)
        bb, _ = _lp(rows)
        lo, hi = bb[LIBRARY.index("x")] if bb else (None, None)
        L = 2 * half * 1e-3
        ou.append({
            "L": L, "theta_L": 1.0 * L, "two_kappa_sq": 2 * band.kappa ** 2,
            "rows": int(len(rows.y)), "n_disjoint": rows.n_disjoint(),
            "median_signal_to_band": r.get("median_signal_to_band"),
            "verdict": ("certified" if r.get("certified")
                        else f"abstain[{r.get('abstain')}]"),
            "x_bound": [lo, hi],
            "x_width": None if lo is None or hi is None else hi - lo,
            "x_resolved": bool(lo is not None and hi is not None
                               and not (lo <= 0.0 <= hi)),
        })
    gbm = []
    from experiments.stochastic.generator import gbm_paths
    for b in (0.3, 0.2, 0.1, 0.05, 0.02, 0.01):
        t, S = gbm_paths(mu=0.8, b=b, T=6.0, dt=1e-4, n_traj=8, seed=1)
        rows = build_rows(t, S, LIBRARY, half=5000)
        if len(rows.y) < 6:
            gbm.append({"b": b, "rows": int(len(rows.y)),
                        "verdict": "too few rows"})
            continue
        r = certify_drift(rows, delta=0.05, seed=0)
        bb, _ = _lp(rows)
        lo, hi = bb[LIBRARY.index("x")] if bb else (None, None)
        gbm.append({
            "b": b, "rows": int(len(rows.y)),
            "median_signal_to_band": r.get("median_signal_to_band"),
            "verdict": ("certified" if r.get("certified")
                        else f"abstain[{r.get('abstain')}]"),
            "law": r.get("law"),
            "x_bound": [lo, hi],
            "x_width": None if lo is None or hi is None else hi - lo,
            "x_resolved": bool(lo is not None and hi is not None
                               and not (lo <= 0.0 <= hi)),
            "x_covers_truth": bool((lo is None or lo <= 0.8)
                                   and (hi is None or 0.8 <= hi)),
        })
    return {"ou_window_length": ou, "gbm_diffusion_intensity": gbm}


# ------------------------------------------------------------------ C: S1

def section_c() -> dict:
    """S1: the interval half-width against the row count and the window length.

    Registered as scaling like the CLT rate sigma/sqrt(N*T). The two halves are
    measured separately because they are NOT the same question: adding rows costs
    kappa(n, delta) = sqrt(2 ln(2n/delta)) in band width, and that growth is
    exactly the size of the averaging gain it would have bought.
    """
    by_rows = []
    for n_traj in (2, 4, 8, 16, 32):
        t, X = ou_paths(theta=1.0, b=1.4, T=640.0, dt=1e-3, n_traj=n_traj,
                        seed=0, x0=0.0)
        rows = build_rows(t, X, LIBRARY, half=32000)
        bb, _ = _lp(rows)
        lo, hi = bb[LIBRARY.index("x")]
        by_rows.append({"n_traj": n_traj, "rows": int(len(rows.y)),
                        "kappa": _band_of(rows).kappa,
                        "x_width": None if lo is None or hi is None else hi - lo})
    by_window = []
    for half, T in ((8000, 160.0), (16000, 320.0), (32000, 640.0),
                    (64000, 1280.0)):
        t, X = ou_paths(theta=1.0, b=1.4, T=T, dt=1e-3, n_traj=8, seed=0,
                        x0=0.0)
        rows = build_rows(t, X, LIBRARY, half=half)
        bb, _ = _lp(rows)
        lo, hi = bb[LIBRARY.index("x")]
        w = None if lo is None or hi is None else hi - lo
        L = 2 * half * 1e-3
        by_window.append({"L": L, "rows": int(len(rows.y)),
                          "kappa": _band_of(rows).kappa, "x_width": w,
                          "width_times_sqrt_L": None if w is None
                          else w * np.sqrt(L)})
    return {"by_row_count_fixed_window": by_rows,
            "by_window_length_fixed_rows": by_window}


# ----------------------------------------------- D: the tasks, through the checker

def section_d(tasks) -> dict:
    """Run every Level 0 task and score it with the FROZEN checker.

    This is also the first time `stochcheck` is fed by a real producer rather than
    by constructed records -- the same role `run_darcy_domains.py` played for the
    domain qualifier.
    """
    halves = {"L0-ou": 32000, "L0-gbm": 5000, "L0-ode-obs": 32000,
              "L0-null-noise": 32000, "L0-null-coarse": 20}
    rows_by_task, results = {}, []
    for task in tasks:
        t, X, info = paths_for(task)
        half = halves[task.task_id]
        sigma_obs = float(task.sampling.get("sigma_obs", 0.0))
        # sigma_obs MUST reach the assembler: it debiases realized quadratic
        # variation, and declaring it only downstream produced the L0-ode-obs
        # confident-wrong. `certify_drift` now raises on the mismatch.
        rows = build_rows(t, X, LIBRARY, half=half, sigma_obs=sigma_obs)
        rows_by_task[task.task_id] = rows
        r = certify_drift(rows, delta=0.05, sigma_obs=sigma_obs, seed=0)
        # BUILD THE SUBMISSION the checker consumes. Declarations name their
        # consumer; the run declares only what it actually used.
        decls = [Declaration(Consumer.DRIFT_BAND,
                             float(np.median(rows.qv)) if len(rows.qv) else 0.0,
                             provenance="measured",
                             note="realized int phi^2 f'^2 d[X], median row")]
        if sigma_obs > 0:
            decls.append(Declaration(Consumer.OBSERVATION, sigma_obs,
                                     provenance="declared"))
        rec = None
        if r.get("partial", {}).get("components"):
            rec = dict(r["partial"])
            rec["components"] = {component("drift", k): v
                                 for k, v in rec["components"].items()}
            for key in ("exact", "interval", "unconstrained"):
                rec[key] = [component("drift", n) for n in rec.get(key, [])]
        cov = Coverage(kappa=r["kappa"], delta=r["delta"], n_rows=r["n_rows"],
                       n_disjoint=r.get("n_disjoint_cert") or r["n_disjoint"],
                       qv_provenance="measured") if r.get("kappa") else None
        if r.get("certified"):
            sub = Submission(
                task_id=task.task_id, kind="answer", record=rec,
                declarations=tuple(decls), coverage=cov,
                alpha_log10=r.get("alpha_log10"),
                alpha_kind=r.get("alpha_kind", ""),
                law=r.get("law", ""), submission_id=task.task_id)
        else:
            # ABSTENTION CARRYING PARTIAL DETERMINATION -- the shape this arc
            # produces. The reason is scored, and the record speaks for the
            # components it determined; the token covers the rest.
            sub = Submission(
                task_id=task.task_id, kind="abstain",
                abstain=r.get("abstain") or "resolution",
                reason_detail="; ".join(r.get("notes", []))[:400],
                record=rec, coverage=cov if rec else None,
                declarations=tuple(decls), submission_id=task.task_id)
        score = score_task(task, [sub])
        results.append({
            "task_id": task.task_id, "system": task.system,
            "verdict": ("certified" if r.get("certified")
                        else f"abstain[{r.get('abstain')}]"),
            "law": r.get("law"), "rows": r["n_rows"],
            "n_disjoint": r["n_disjoint"], "rejected": r["rejected"],
            "kappa": r.get("kappa"), "sigma_effective": r.get("sigma_effective"),
            "alpha_log10": r.get("alpha_log10"),
            "median_signal_to_band": r.get("median_signal_to_band"),
            "median_martingale_share": r.get("median_martingale_share"),
            "submission_kind": sub.kind,
            "score": {
                "confident_wrong": score.n_confident_wrong,
                "exact": score.n_exact, "covered": score.n_covered,
                "informative": score.n_informative,
                "abstained_correctly": score.n_abstained_correctly,
                "missed": score.n_missed, "resolved": score.n_resolved,
                "refusals": score.refusals,
                "abstention": score.abstention,
                "declaration_audit": score.declaration_audit,
                "components": {k: {"outcome": v["outcome"], "truth": v["truth"],
                                   "lo": v["lo"], "hi": v["hi"],
                                   "resolved": v["resolved"],
                                   "informative": v["informative"]}
                               for k, v in score.components.items()},
            },
            "notes": r.get("notes", [])[:4],
        })
        results[-1]["_score_obj"] = score
    scores = [x.pop("_score_obj") for x in results]
    totals = suite_totals(scores)
    return {"tasks": results, "totals": totals,
            "ranking": rank({"lagh-ito-L0": scores})}, rows_by_task


# ------------------------------------------------------------------ E: S6

def section_e() -> dict:
    """Realized quadratic variation where there is NO process noise.

    Observational error contributes 2 sigma_obs^2 per step to [X], so the implied
    diffusion is 2 sigma_obs^2 / dt and DIVERGES as the sampling rate rises. The
    prediction and its consequence: the band is inflated (conservative -- the safe
    direction), and quadratic variation alone cannot separate process noise from
    measurement noise, which is precisely what Level 2 is built to supply.
    """
    from experiments.stochastic.generator import ode_obs_paths
    rows_out = []
    for dt in (4e-3, 2e-3, 1e-3, 5e-4):
        t, X = ode_obs_paths(theta=1.0, sigma_obs=1e-3, T=640.0, dt=dt,
                             n_traj=4, seed=7)
        half = int(round(32.0 / dt))
        # deliberately NOT debiased here: this section measures the raw estimator's
        # contamination, which is the quantity S6 is about
        rows = build_rows(t, X, LIBRARY, half=half)
        if not len(rows.y):
            rows_out.append({"dt": dt, "rows": 0})
            continue
        # implied b^2 from qv, on the f = x rows (f' = 1 there, so
        # qv = int phi^2 b^2 dt and the phi^2 integral is known)
        m = np.array(rows.fname) == "x"
        if not m.any():
            continue
        # int phi^2 dt over a window of half-width at, for psi = (1-s^2)^8
        at = half * dt
        int_phi2 = 0.4331 * at
        implied = float(np.median(rows.qv[m]) / int_phi2)
        rows_out.append({
            "dt": dt, "rows": int(len(rows.y)),
            "implied_b_squared": implied,
            "predicted_2_sigma_obs_sq_over_dt": 2 * (1e-3) ** 2 / dt,
            "ratio": implied / (2 * (1e-3) ** 2 / dt),
            "true_b_squared": 0.0,
        })
    return {"spurious_diffusion_from_observation_noise": rows_out,
            "note": ("no process noise exists in this system; every one of these "
                     "implied diffusions is spurious, and the value scales as "
                     "1/dt exactly as predicted. The band it produces is too "
                     "WIDE, which is the conservative direction")}


# ------------------------------------------------------------------ F: S7

def section_f() -> dict:
    """S7: alpha against the disjoint window count.

    Measured on GBM rather than OU, and the reason is itself a result: for a
    STATIONARY drift the martingale band exceeds the target's own range, so
    q = min(1, 2 eps / range(y)) saturates at 1 and alpha is VACUOUS -- the
    significance bound says nothing in the regime where the drift is pure
    fluctuation. GBM's drift accumulates, q < 1, and the trade S7 predicts is
    visible there.
    """
    from experiments.stochastic.generator import gbm_paths

    def alpha_of(rows, dof=1):
        band = _band_of(rows)
        e = band(None)
        R = float(np.max(rows.y) - np.min(rows.y))
        q = float(np.median(np.minimum(1.0, 2.0 * e / max(R, 1e-300))))
        nd = rows.n_disjoint()
        h = max(0, nd - dof)
        return {"rows": int(len(rows.y)), "n_disjoint": nd,
                "kappa": band.kappa, "median_q": q, "h_disjoint": h,
                "log10_alpha_at_h": float(h * np.log10(max(q, 1e-300))),
                "vacuous": bool(q >= 1.0)}

    gbm = []
    for n_traj in (2, 4, 8, 16):
        t, S = gbm_paths(mu=0.8, b=0.01, T=6.0, dt=1e-4, n_traj=n_traj, seed=1)
        rows = build_rows(t, S, LIBRARY, half=5000)
        gbm.append({"n_traj": n_traj, **alpha_of(rows)})
    ou = []
    for n_traj in (2, 8):
        t, X = ou_paths(theta=1.0, b=1.4, T=640.0, dt=1e-3, n_traj=n_traj,
                        seed=0, x0=0.0)
        rows = build_rows(t, X, LIBRARY, half=32000)
        ou.append({"n_traj": n_traj, **alpha_of(rows)})
    return {"gbm_alpha_vs_disjoint_count": gbm,
            "ou_alpha_is_vacuous": ou,
            "note": ("kappa grows like sqrt(log n) while the exponent grows "
                     "linearly in the disjoint count, so more windows tighten "
                     "alpha faster than they loosen the band -- S7. On a "
                     "stationary drift the band exceeds the target's range, q "
                     "saturates at 1 and alpha is vacuous: the trade is "
                     "measurable only where the drift accumulates")}


def main():
    t0 = time.time()
    tasks = level0_tasks()
    print("A: kappa calibration ...", flush=True)
    a = section_a()
    print("B: reach boundary ...", flush=True)
    b = section_b()
    print("C: S1 scaling ...", flush=True)
    c = section_c()
    print("D: tasks through the frozen checker ...", flush=True)
    d, _ = section_d(tasks)
    print("E: S6 spurious quadratic variation ...", flush=True)
    e = section_e()
    print("F: S7 alpha vs disjoint count ...", flush=True)
    f = section_f()
    res = {"level": 0, "seconds": round(time.time() - t0, 1),
           "library": list(LIBRARY),
           "A_kappa_calibration": a, "B_reach_boundary": b, "C_S1_scaling": c,
           "D_tasks_scored": d, "E_S6_spurious_qv": e, "F_S7_alpha": f}
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(res, indent=1, default=str))
    print(f"\nwrote {OUT}  ({res['seconds']}s)")
    print("\n--- A: does the band hold at its claimed budget? ---")
    for sysname, blk in a.items():
        for dl, r in blk["by_delta"].items():
            print(f"  {sysname} delta={dl:>5} kappa={r['kappa']:.2f} "
                  f"failures {r['failures']}/{r['n_replicates']} "
                  f"(rate {r['empirical_rate']:.3f}) holds={r['holds']} "
                  f"max|r|/band med {r['max_resid_over_band_median']:.3f}")
    print("\n--- D: the five tasks ---")
    for x in d["tasks"]:
        s = x["score"]
        print(f"  {x['task_id']:16s} {x['verdict']:22s} "
              f"CW={s['confident_wrong']} cov={s['covered']} inf={s['informative']} "
              f"abst-ok={s['abstained_correctly']} missed={s['missed']} "
              f"resolved={s['resolved']}")
    print("  TOTALS:", d["totals"])


if __name__ == "__main__":
    main()

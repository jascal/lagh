"""Level 1 of the stochastic suite, scored through the FROZEN checker.

Four systems, drift AND diffusion per component, each verdict submitted as an
abstention-or-answer carrying its partial determination. What this run is evidence
for, precisely:

  * **S5** -- zero confident-wrong across Level 1, which Level 0 established and
    this extends to nonlinear drift, state-dependent diffusion, and a system whose
    state space is disconnected.
  * the frozen checker exercised on the DIFFUSION path for the first time, i.e.
    against a producer emitting two component families at once.

What it is NOT evidence for: the expectations. Those are CALIBRATED from the Level 0
run and the two Level 1 increments, not registered blind, and `level1_tasks` says so
in its docstring. A blind Level 1 would need systems nobody has probed.

Run: .venv/bin/python experiments/stochastic/run_level1.py
"""
from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import numpy as np
import sympy as sp

sys.path.insert(0, str(Path(__file__).parent.parent.parent))

from experiments.stochastic.generator import (DIFF_LIBRARY, LIBRARY,  # noqa: E402
                                              l1_paths, level1_tasks)
from lagh.certify import admissible_functional  # noqa: E402
from lagh.ito import (LAM_QV, ItoBand, build_qv_rows, build_rows,  # noqa: E402
                      certify_diffusion, certify_drift)
from lagh.stochcheck import (Consumer, Coverage, Declaration,  # noqa: E402
                             Submission, rank, score_task, suite_totals)

OUT = Path("experiments/results/stochastic_level1.json")

# Per task: the drift window, the QV window and the w-family. Chosen from the
# measured identifiability arithmetic (theta*L > 2 kappa^2 for the drift, state
# spread for the diffusion).
#
# `drift_max` is NOT here, and that is the point. The diffusion's band needs a bound
# on |a| because E[(dX)^2] = b^2 dt + a^2 dt^2, and the first version of this run
# hand-declared one per task -- all four of which turned out UNDER-declared against
# the truth, GBM by 96x. A hand-declared band input that nobody verified is the
# dangerous direction, and it is the third time this session that pattern appeared.
# So the bound is now DERIVED from the drift's own admissible envelope
# (`drift_envelope_max`): measured QV -> drift band -> drift envelope -> diffusion
# band, one direction, no circularity. When the drift is undetermined the envelope
# is wide and the diffusion band widens with it, which is the honest coupling.
CONFIG = {
    "L1-dw-additive":       dict(half=32000, qv_half=32000, ws=("1", "x", "x**2")),
    "L1-dw-multiplicative": dict(half=32000, qv_half=32000, ws=("1", "x", "x**2")),
    "L1-cir":               dict(half=8000, qv_half=8000, ws=("1", "x", "x**2")),
    "L1-gbm-mult":          dict(half=5000, qv_half=5000, ws=("1", "x", "x**2")),
}


def drift_envelope(rows, X, *, delta=0.05, n_query=61):
    """A POINTWISE bound on |a(x)| over the visited range, across every drift
    consistent with the data -- what the diffusion's band needs, derived rather than
    declared.

    One LP per query state via `certify.admissible_functional`, returned as a
    callable that interpolates the envelope. Pointwise rather than its maximum
    because the leakage term is `sum_i |w_i| a(X_i)^2 dt^2`, and a process spends its
    time where its own drift is SMALL -- most of all a bistable one, whose wells sit
    at a(x) = 0. Bounding a^2 by its worst value anywhere in the range over-declares
    by the square of a ratio that is an order of magnitude here, and that flat bound
    was measured to destroy the diffusion's precision entirely.

    Returns (callable, max) or (None, None) when the envelope is unbounded -- itself
    the answer: an unbounded drift envelope means the leakage cannot be bounded, and
    the run must say so rather than pick a number.
    """
    if len(rows.y) < 6:
        return None, None
    band = ItoBand(rows.qv, rows.qv_se, rows.corr_se, rows.quad, len(rows.y),
                   delta=delta, y=rows.y, feat_names=list(rows.names))
    e0 = band.martingale() + LAM_QV * rows.corr_se + rows.quad[:, 0]
    q = rows.quad[:, 1:].sum(axis=1)
    grid = np.linspace(float(X.min()), float(X.max()), n_query)
    W = np.array([[float(sp.sympify(g).subs(sp.Symbol("x"), z))
                   for g in LIBRARY] for z in grid])
    got, _ = admissible_functional(rows.A, rows.y, lambda c: e0 + c * q, W,
                                   coeff_max=50.0)
    if got is None or any(lo is None or hi is None for lo, hi in got):
        return None, None
    env = np.array([max(abs(lo), abs(hi)) for lo, hi in got])
    return (lambda z: np.interp(np.asarray(z, float), grid, env),
            float(env.max()))


def _record(part_rows, prefix):
    """Re-key a determination record's components into the frozen component
    vocabulary. The drift path names its columns bare; the QV path already
    prefixes them."""
    rec = dict(part_rows)
    comps = {}
    for k, v in rec.get("components", {}).items():
        comps[k if ":" in k else f"{prefix}:{k}"] = v
    rec["components"] = comps
    for key in ("exact", "interval", "unconstrained"):
        rec[key] = [k if ":" in k else f"{prefix}:{k}" for k in rec.get(key, [])]
    return rec


def _merge(a, b):
    """Two records, one submission: the drift's components and the diffusion's."""
    if a is None:
        return b
    if b is None:
        return a
    out = dict(a)
    out["components"] = {**a.get("components", {}), **b.get("components", {})}
    for key in ("exact", "interval", "unconstrained"):
        out[key] = list(a.get(key, [])) + list(b.get(key, []))
    out["status"] = a.get("status", "conjoined")
    out["note"] = ("drift from the Itô weak form, diffusion from realized "
                   "quadratic variation -- two estimators, one submission, "
                   "because they are the efficient estimator for their own part "
                   "(docs/CASE_STUDY_STOCHASTIC_L1.md)")
    return out


def main():
    t0 = time.time()
    rows_out, scores = [], []
    for task in level1_tasks():
        cfg = CONFIG[task.task_id]
        t, X = l1_paths(task)
        # --- the DRIFT, from the Itô weak form
        dr_rows = build_rows(t, X, LIBRARY, half=cfg["half"])
        dr = certify_drift(dr_rows, delta=0.05, seed=0) if len(dr_rows.y) else {}
        # the bound on |a| the diffusion's band needs, DERIVED from the drift's
        # admissible envelope over the visited range
        env_fn, dmax = drift_envelope(dr_rows, X)
        # --- the DIFFUSION, from realized quadratic variation
        qv_rows = build_qv_rows(t, X, DIFF_LIBRARY, half=cfg["qv_half"],
                                ws=cfg["ws"], drift_envelope=env_fn)
        df = certify_diffusion(qv_rows, delta=0.05, drift_max=dmax or 0.0, seed=0) \
            if len(qv_rows.y) and env_fn is not None else {}
        rec = _merge(_record(dr["partial"], "drift") if dr.get("partial") else None,
                     _record(df["partial"], "diffusion") if df.get("partial")
                     else None)
        decls = [
            Declaration(Consumer.DRIFT_BAND,
                        float(np.median(dr_rows.qv)) if len(dr_rows.qv) else 0.0,
                        provenance="measured",
                        note="realized int phi^2 f'^2 d[X], median row"),
            Declaration(Consumer.DIFFUSION_QV,
                        float(np.median(qv_rows.y)) if len(qv_rows.y) else 0.0,
                        provenance="measured",
                        note="realized int phi w d[X], median row; the DIFFUSION's "
                             "own target rather than a band scale"),
        ]
        cov = Coverage(kappa=dr.get("kappa") or df.get("kappa"),
                       delta=0.05,
                       n_rows=(dr.get("n_rows") or 0) + (df.get("n_rows") or 0),
                       n_disjoint=(dr.get("n_disjoint_cert")
                                   or dr.get("n_disjoint") or 0)
                       + (df.get("n_disjoint") or 0),
                       qv_provenance="measured") if rec else None
        certified = bool(dr.get("certified")) and bool(df.get("certified"))
        if certified:
            sub = Submission(task_id=task.task_id, kind="answer", record=rec,
                             declarations=tuple(decls), coverage=cov,
                             law=f"drift {dr.get('law','')} | "
                                 f"diffusion {df.get('law','')}",
                             submission_id=task.task_id)
        else:
            sub = Submission(
                task_id=task.task_id, kind="abstain",
                abstain=(dr.get("abstain") or df.get("abstain") or "structural"),
                reason_detail=(f"drift: {dr.get('abstain')}; "
                               f"diffusion: {df.get('abstain')}")[:400],
                record=rec, coverage=cov, declarations=tuple(decls),
                submission_id=task.task_id)
        score = score_task(task, [sub])
        scores.append(score)
        rows_out.append({
            "task_id": task.task_id, "system": task.system,
            "drift": {"verdict": ("certified" if dr.get("certified")
                                  else f"abstain[{dr.get('abstain')}]"),
                      "law": dr.get("law"), "rows": dr.get("n_rows"),
                      "kappa": dr.get("kappa"),
                      "signal_to_band": dr.get("median_signal_to_band")},
            "diffusion": {"verdict": ("certified" if df.get("certified")
                                      else f"abstain[{df.get('abstain')}]"),
                          "law": df.get("law"), "rows": df.get("n_rows"),
                          "drift_leakage_share":
                              df.get("median_drift_leakage_share"),
                          "leakage_bound": df.get("leakage_bound")},
            "state_range": [float(X.min()), float(X.max())],
            # std, NOT std/|mean|: the first version reported the ratio, which is
            # meaningless for a symmetric double well whose mean is ~0 (it read 74
            # and 238). The diffusion's identifiability tracks the spread relative
            # to where the library terms differ, so both are reported raw.
            "state_std": float(X.std()),
            "drift_envelope_max": dmax,
            "score": {"confident_wrong": score.n_confident_wrong,
                      "exact": score.n_exact, "covered": score.n_covered,
                      "informative": score.n_informative,
                      "abstained_correctly": score.n_abstained_correctly,
                      "missed": score.n_missed, "resolved": score.n_resolved,
                      "refusals": score.refusals,
                      "exceeded_expectation": score.exceeded_expectation},
            "components": {k: {"outcome": v["outcome"], "truth": v["truth"],
                               "lo": v["lo"], "hi": v["hi"],
                               "resolved": v["resolved"]}
                           for k, v in score.components.items()},
        })
        print(f"  {task.task_id:22s} drift {rows_out[-1]['drift']['verdict']:22s} "
              f"diffusion {rows_out[-1]['diffusion']['verdict']:22s} "
              f"CW={score.n_confident_wrong} res={score.n_resolved}", flush=True)
    totals = suite_totals(scores)
    res = {"level": 1, "seconds": round(time.time() - t0, 1),
           "drift_library": list(LIBRARY), "diffusion_library": list(DIFF_LIBRARY),
           "tasks": rows_out, "totals": totals,
           "ranking": rank({"lagh-ito-L1": scores}),
           "note": ("expectations are CALIBRATED from Level 0 and the Level 1 "
                    "increments, not registered blind -- see level1_tasks")}
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(res, indent=1, default=str))
    print(f"\nwrote {OUT} ({res['seconds']}s)")
    print("TOTALS:", totals)


if __name__ == "__main__":
    main()

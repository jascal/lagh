"""Van der Pol in the scored task set: a 2-D state, one submission per EQUATION.

Separate from `run_level1.py` because the shape of the claim is different, not
because the machinery is. A scalar task has one drift and one diffusion; a
2-D system has one drift equation PER COMPONENT, and the frozen interface's answer
to that is `part[index]:term` plus several submissions per task -- both of which
were in the interface from the freeze and neither of which had a producer until now.

Each equation is submitted on its own, with its own coverage statement, and the
checker scores the union. The diffusion is submitted once, from quadratic variation,
per component.

Run: .venv/bin/python experiments/stochastic/run_level1_vdp.py
"""
from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).parent.parent.parent))

from experiments.stochastic.generator import (VDP_FIELDS,  # noqa: E402
                                              VDP_LIBRARIES, vdp_paths, vdp_task)
from lagh.ito import (build_qv_rows, build_rows_nd, certify_diffusion,  # noqa: E402
                      certify_drift)
from lagh.stochcheck import (Consumer, Coverage, Declaration,  # noqa: E402
                             Submission, score_task, suite_totals)

OUT = Path("experiments/results/stochastic_level1_vdp.json")

# One test function per component: f = x reads the x equation (df/dx = 1, df/dy = 0),
# f = y^2/2 reads the y equation. Both were measured per-equation before this run --
# the noise-free x equation carries 30x more signal than its band, the driven y
# equation is vacuous -- so the pairing is chosen, not guessed.
F_FOR = {"x": "x", "y": "y**2/2"}
HALF = 8000
DIFF_LIBRARY = ("1", "x", "x**2")


def _reindex(record, index: int):
    """Re-key `<field>:<term>` (or a bare term) into the frozen
    `<part>[<index>]:<term>` convention. This is the whole of what the component
    INDEX exists for."""
    out = dict(record)
    comps = {}
    for k, v in record.get("components", {}).items():
        term = k.split(":", 1)[1] if ":" in k else k
        comps[f"drift[{index}]:{term}"] = v
    out["components"] = comps
    for key in ("exact", "interval", "unconstrained"):
        out[key] = [f"drift[{index}]:{k.split(':', 1)[1] if ':' in k else k}"
                    for k in record.get(key, [])]
    return out


def _reindex_diffusion(record, index: int):
    out = dict(record)
    out["components"] = {f"diffusion[{index}]:{k.split(':', 1)[1]}": v
                         for k, v in record.get("components", {}).items()}
    for key in ("exact", "interval", "unconstrained"):
        out[key] = [f"diffusion[{index}]:{k.split(':', 1)[1]}"
                    for k in record.get(key, [])]
    return out


def main():
    t0 = time.time()
    task = vdp_task()
    s = task.sampling
    t, X, Y = vdp_paths(mu=s["mu"], b=s["b"], T=s["T"], dt=s["dt"],
                        n_traj=s["n_traj"], seed=s["seed"],
                        substeps=s["substeps"])
    fields = {"x": X, "y": Y}
    subs, report = [], []

    # ---- one DRIFT submission per component equation
    for i, fld in enumerate(VDP_FIELDS):
        rows = build_rows_nd(t, fields, VDP_LIBRARIES, F_FOR[fld], half=HALF)
        r = certify_drift(rows, delta=0.05, seed=0)
        rec = _reindex(r["partial"], i) if r.get("partial") else None
        cov = Coverage(kappa=r["kappa"], delta=r["delta"], n_rows=r["n_rows"],
                       n_disjoint=r.get("n_disjoint_cert") or r["n_disjoint"],
                       qv_provenance="measured") if r.get("kappa") else None
        kind = "answer" if r.get("certified") else "abstain"
        subs.append(Submission(
            task_id=task.task_id, kind=kind, record=rec, coverage=cov,
            abstain=None if kind == "answer" else (r.get("abstain") or "noise"),
            reason_detail=f"component {fld}: {r.get('abstain')}",
            alpha_log10=r.get("alpha_log10"), alpha_kind=r.get("alpha_kind", ""),
            law=r.get("law", ""), submission_id=f"drift-{fld}"))
        report.append({"part": f"drift[{i}] ({fld})", "f": F_FOR[fld],
                       "verdict": ("certified" if r.get("certified")
                                   else f"abstain[{r.get('abstain')}]"),
                       "law": r.get("law"), "rows": r["n_rows"],
                       "signal_to_band": r.get("median_signal_to_band"),
                       "alpha_log10": r.get("alpha_log10")})

    # ---- the DIFFUSION, per component, from quadratic variation
    for i, fld in enumerate(VDP_FIELDS):
        qv = build_qv_rows(t, fields[fld], DIFF_LIBRARY, half=HALF,
                           ws=("1", "x", "x**2"))
        rd = certify_diffusion(qv, delta=0.05, drift_max=8.0, seed=0)
        rec = _reindex_diffusion(rd["partial"], i) if rd.get("partial") else None
        cov = Coverage(kappa=rd["kappa"], delta=rd["delta"], n_rows=rd["n_rows"],
                       n_disjoint=rd["n_disjoint"], qv_provenance="measured") \
            if rd.get("kappa") else None
        kind = "answer" if rd.get("certified") else "abstain"
        subs.append(Submission(
            task_id=task.task_id, kind=kind, record=rec, coverage=cov,
            abstain=None if kind == "answer" else (rd.get("abstain") or "noise"),
            reason_detail=f"diffusion of {fld}: {rd.get('abstain')}",
            declarations=(Declaration(Consumer.DIFFUSION_QV,
                                      float(np.median(qv.y)) if len(qv.y) else 0.0,
                                      provenance="measured"),),
            law=rd.get("law", ""), submission_id=f"diffusion-{fld}"))
        report.append({"part": f"diffusion[{i}] ({fld})",
                       "verdict": ("certified" if rd.get("certified")
                                   else f"abstain[{rd.get('abstain')}]"),
                       "law": rd.get("law"), "rows": rd["n_rows"],
                       "signal_to_band": rd.get("median_signal_to_band")})

    score = score_task(task, subs)
    res = {"task": task.task_id, "state_dim": 2,
           "seconds": round(time.time() - t0, 1),
           "parts": report,
           "score": {"confident_wrong": score.n_confident_wrong,
                     "exact": score.n_exact, "covered": score.n_covered,
                     "informative": score.n_informative,
                     "abstained_correctly": score.n_abstained_correctly,
                     "missed": score.n_missed, "resolved": score.n_resolved,
                     "reach": list(score.reach), "refusals": score.refusals,
                     "exceeded_expectation": score.exceeded_expectation},
           "components": {k: {"outcome": v["outcome"], "truth": v["truth"],
                              "lo": v["lo"], "hi": v["hi"],
                              "resolved": v["resolved"],
                              "informative": v["informative"]}
                          for k, v in score.components.items()},
           "totals": suite_totals([score])}
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(res, indent=1, default=str))
    for p in report:
        print(f"  {p['part']:22s} {p['verdict']:24s} S/B="
              f"{(p.get('signal_to_band') or 0):7.2f}  {p.get('law') or ''}")
    print(f"\nwrote {OUT} ({res['seconds']}s)")
    print("SCORE:", res["score"] | {"refusals": len(score.refusals)})


if __name__ == "__main__":
    main()

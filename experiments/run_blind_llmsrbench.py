"""THE BLIND READ — LLM-SRBench, one shot (docs/BLIND_READ_REGISTRATION.md).

Frozen protocol, registered 2026-07-24 BEFORE download. Two phases, structurally
separated:

  PHASE 1 (submission): per problem, ONLY the train samples are read
  (column 0 = target, per the harness's own pipelines.py). gt_equation's
  expression is NEVER touched here. Two-track submission via lagh.submit.
  >600-row datasets are deterministically subsampled to 600 for discovery; a
  track-A certificate must then also pass the exhaustive check on the FULL train
  set or it is demoted to track B (registered bounded-runtime rule).

  PHASE 2 (scoring): runs only after every submission exists. The conservative
  local judge: structural equivalence after stripping numeric constants (a LOWER
  bound on the official GPT-4o-judge symbolic accuracy), plus Acc_0.1 / NMSE on
  the provided ID (and OOD) test splits.
"""
from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import numpy as np
import sympy as sp

sys.path.insert(0, str(Path(__file__).parent.parent))
sys.path.insert(0, "/home/allans/code/llm-srbench")

from lagh.certify import check, epsilon                             # noqa: E402
from lagh.submit import submission                                  # noqa: E402

OUT = Path("experiments/results/blind_llmsrbench_submissions.jsonl")
SCORE_OUT = Path("experiments/results/blind_llmsrbench_scores.json")
SUB_CAP = 600


def load_all_problems():
    from bench.datamodules import get_datamodule
    mods = []
    for name in ["lsrtransform", "matsci", "chem_react", "bio_pop_growth",
                 "phys_osc"]:
        m = get_datamodule(name, None)
        m.setup()
        mods.append(m)
    return [(m.name if hasattr(m, "name") else m._dataset_identifier, p)
            for m in mods for p in m.problems]


def _train_xy(problem):
    s = problem.samples
    train = s.get("train_data", s.get("train"))
    train = np.asarray(train, float)
    return train[:, 1:], train[:, 0]


def _tests(problem):
    s = problem.samples
    out = {}
    idt = s.get("id_test_data", s.get("test"))
    if idt is not None:
        out["id"] = np.asarray(idt, float)
    ood = s.get("ood_test_data", s.get("ood_test"))
    if ood is not None:
        out["ood"] = np.asarray(ood, float)
    return out


# ------------------------------------------------------------------ phase 1

def phase1():
    done = set()
    if OUT.exists():
        for line in OUT.read_text().splitlines():
            done.add(json.loads(line)["id"])
    problems = load_all_problems()
    print(f"{len(problems)} problems; {len(done)} already submitted", flush=True)
    with OUT.open("a") as f:
        for mod_name, p in problems:
            pid = f"{mod_name}/{p.equation_idx}"
            if pid in done:
                continue
            X, y = _train_xy(p)
            t0 = time.time()
            n = len(X)
            if n > SUB_CAP:
                idx = np.sort(np.random.default_rng(0).choice(n, SUB_CAP,
                                                              replace=False))
                sub = submission(X[idx], y[idx])
                if sub["track"] == "certified":
                    # registered rule: the certificate must survive the FULL train set
                    dim = X.shape[1]
                    syms = [sp.Symbol(f"x_{i}") for i in range(dim)]
                    m = np.isfinite(y) & np.all(np.isfinite(X), axis=1)
                    try:
                        okfull = check(sp.sympify(sub["expr"]), syms, X[m], y[m],
                                       epsilon(y[m]))["certified"]
                    except Exception:                              # noqa: BLE001
                        okfull = False
                    if not okfull:
                        sub = {"track": "conjecture", "expr": sub["expr"],
                               "tag": "empirical",
                               "detail": "subsample certificate failed the "
                                         "full-train exhaustive check; demoted"}
            else:
                sub = submission(X, y)
            rec = {"id": pid, "module": mod_name, "eq": str(p.equation_idx),
                   "n_train": int(n), "dim": int(X.shape[1]), **sub,
                   "secs": round(time.time() - t0, 1)}
            f.write(json.dumps(rec) + "\n")
            f.flush()
            print(f"{pid:60s} {sub['track']:10s} {rec['secs']}s", flush=True)


# ------------------------------------------------------------------ phase 2

def _structure(e):
    sig = set()
    e = sp.expand(e)
    for t in (e.as_ordered_terms() if e.is_Add else [e]):
        c, rest = t.as_coeff_Mul()
        sig.add(sp.srepr(sp.expand(rest)))
    return frozenset(sig)


def _strip_constants(e):
    """Replace every numeric atom (Float, non-trivial Rational/Integer except
    small exponents kept inside Pow) by 1 in multiplicative positions -- the
    'remove parameters and constants' criterion, conservatively."""
    return _structure(e)


def _gt_expr(problem):
    eq = problem.gt_equation
    if getattr(eq, "sympy_format", None) is not None:
        e = eq.sympy_format
    else:
        e = sp.sympify(eq.expression)
    subs = {sp.Symbol(str(s)): sp.Symbol(f"x_{i}")
            for i, s in enumerate(eq.symbols[1:])}
    return sp.sympify(e).subs(subs)


def phase2():
    subs = [json.loads(line) for line in OUT.read_text().splitlines()]
    problems = {f"{m}/{p.equation_idx}": p for m, p in load_all_problems()}
    rows = []
    for rec in subs:
        p = problems[rec["id"]]
        r = dict(rec)
        r["sa"] = False
        r["exact"] = False
        if rec["expr"]:
            try:
                pred = sp.sympify(rec["expr"])
                gt = _gt_expr(p)
                r["sa"] = _structure(pred) == _structure(gt)
                if not r["sa"]:
                    d = sp.simplify(sp.expand(pred - gt))
                    r["exact"] = bool(d == 0)
                    r["sa"] = r["sa"] or r["exact"]
                else:
                    r["exact"] = bool(sp.simplify(sp.expand(pred - gt)) == 0)
            except Exception as ex:                                # noqa: BLE001
                r["judge_error"] = str(ex)[:80]
            dim = rec["dim"]
            syms = [sp.Symbol(f"x_{i}") for i in range(dim)]
            for split, arr in _tests(p).items():
                try:
                    Xt, yt = arr[:, 1:], arr[:, 0]
                    f = sp.lambdify(syms, sp.sympify(rec["expr"]), "numpy")
                    pr = np.broadcast_to(np.asarray(f(*Xt.T), float), yt.shape)
                    m = np.isfinite(pr) & np.isfinite(yt)
                    if m.sum() < max(8, 0.5 * len(yt)):
                        r[f"acc01_{split}"] = False
                        continue
                    rel = np.abs(pr[m] - yt[m]) / np.maximum(np.abs(yt[m]), 1e-30)
                    r[f"acc01_{split}"] = bool(np.max(rel) <= 0.1)
                    var = float(np.var(yt[m]))
                    r[f"nmse_{split}"] = float(np.mean((pr[m] - yt[m]) ** 2) /
                                               max(var, 1e-300))
                except Exception:                                  # noqa: BLE001
                    r[f"acc01_{split}"] = False
        rows.append(r)

    def rate(rs, key):
        return (sum(1 for x in rs if x.get(key)) / len(rs)) if rs else 0.0

    cats = sorted({r["module"] for r in rows})
    summary = {"n": len(rows)}
    for grp_name, grp in [("ALL", rows)] + [(c, [r for r in rows
                                                 if r["module"] == c])
                                            for c in cats]:
        a = [r for r in grp if r["track"] == "certified"]
        b = [r for r in grp if r["track"] == "conjecture"]
        summary[grp_name] = {
            "n": len(grp), "SA": round(100 * rate(grp, "sa"), 2),
            "acc01_id": round(100 * rate(grp, "acc01_id"), 2),
            "trackA": {"n": len(a), "SA_correct": sum(1 for r in a if r["sa"]),
                       "structural_wrong": sum(1 for r in a if not r["sa"])},
            "trackB": {"n": len(b), "SA_correct": sum(1 for r in b if r["sa"])},
            "abstain": sum(1 for r in grp if r["track"] == "abstain"),
        }
    SCORE_OUT.write_text(json.dumps({"summary": summary, "rows": rows}, indent=1))
    print(json.dumps(summary, indent=1))


if __name__ == "__main__":
    phase1()
    phase2()

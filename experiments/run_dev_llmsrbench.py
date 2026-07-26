"""LLM-SRBench DEV sweep v1 (docs/LLMSRBENCH_DEV.md — dev metrics, never a headline).

Changes vs the spent blind protocol, all registered in the dev doc BEFORE running:
sigma_rep = 1e-4 declared globally (float32 + gradient-amplification envelope);
k-form LLM proposer on abstain (verified forms are certificates — same sound
checker; unverified best proposal replaces the log-log conjecture); judge fixed
for functional notation; 8-way parallel; 600 s/problem cap.
"""
from __future__ import annotations

import json
import os
import signal
import sys
import time
from pathlib import Path

import numpy as np
import sympy as sp

sys.path.insert(0, str(Path(__file__).parent.parent))
sys.path.insert(0, "/home/allans/code/llm-srbench")

SIGMA_REP = 1e-4
SUB_CAP = 600
# v1.1: per-STAGE caps (registered amendment in LLMSRBENCH_DEV.md). v1's single
# 600s cap starved the proposer: sigma-widened discovery is slower, so timeouts
# fired before the LLM stage ever ran and emitted bare abstains -- strictly worse
# than the blind run's fallback. Discovery gets 350s; the proposer+verify stage
# gets its own 200s; the log-log fallback is always reachable.
DISCOVER_CAP_S = 350
PROPOSE_CAP_S = 200
WORKERS = 8
_TAG = os.environ.get("DEV_SWEEP_TAG", "v1")
OUT = Path(f"experiments/results/dev_llmsrbench_{_TAG}.jsonl")
SCORE_OUT = Path(f"experiments/results/dev_llmsrbench_{_TAG}_scores.json")


def _timeout_handler(signum, frame):
    raise TimeoutError("per-problem cap")


def _fallback_conjecture(X, y):
    """Cheap conjectures, best first: CAP-T unsnapped generalized monomial
    (continuous exponents -- the Arrhenius class), else the log-log probe.
    Always reachable, never times out."""
    from lagh.classes.c9_genmonomial import conjecture as c9_conjecture
    try:
        c = c9_conjecture(X, y)
        if c is not None:
            return c[0]
    except Exception:                                          # noqa: BLE001
        pass
    m = np.isfinite(y) & np.all(np.isfinite(X), axis=1)
    X, y = X[m], y[m]
    if len(X) < 10 or not (np.all(X > 0) and np.all(y > 0)):
        return None
    from fractions import Fraction
    L = np.column_stack([np.ones(len(X)), np.log(X)])
    c, *_ = np.linalg.lstsq(L, np.log(y), rcond=None)
    expr = sp.Float(float(np.exp(c[0])))
    for i, a in enumerate(c[1:]):
        fr = Fraction(float(a)).limit_denominator(12)
        expr *= sp.Symbol(f"x_{i}") ** sp.Rational(fr.numerator, fr.denominator)
    return str(expr)


def solve_one(args):
    """Worker: one problem -> submission record. Per-STAGE alarm caps (v1.1) so a
    slow discovery can never starve the proposer or the fallback."""
    pid, X, y = args
    from lagh.certify import check, epsilon
    from lagh.characterize import characterize
    from lagh.mcp.core import verify
    from lagh.passive import discover_passive
    from machine.llm import propose_forms

    signal.signal(signal.SIGALRM, _timeout_handler)
    t0 = time.time()
    rec = {"id": pid, "n_train": int(len(X)), "dim": int(X.shape[1]),
           "timed_out": False}
    try:
        n = len(X)
        Xd, yd = X, y
        if n > SUB_CAP:
            idx = np.sort(np.random.default_rng(0).choice(n, SUB_CAP, replace=False))
            Xd, yd = X[idx], y[idx]

        # ---- stage 1: grammar discovery (DISCOVER_CAP_S) ----
        r = None
        try:
            signal.alarm(DISCOVER_CAP_S)
            r = discover_passive(Xd, yd, sigma=SIGMA_REP)
        except TimeoutError:
            rec["timed_out"] = True
        finally:
            signal.alarm(0)
        if r is not None and r.certified:
            syms = [sp.Symbol(f"x_{i}") for i in range(X.shape[1])]
            m = np.isfinite(y) & np.all(np.isfinite(X), axis=1)
            ok = check(r.result.expr, syms, X[m], y[m],
                       epsilon(y[m], sigma=SIGMA_REP))["certified"] if n > SUB_CAP \
                else True
            if ok:
                rec.update(track="certified", expr=str(r.result.expr),
                           tag="empirical-structural", channel="grammar",
                           detail=f"certified at declared sigma_rep={SIGMA_REP}")
                return rec

        # ---- stage 2: characterize -> k-form LLM proposer -> verify (PROPOSE_CAP_S)
        forms = []
        try:
            signal.alarm(PROPOSE_CAP_S)
            abstain = r.result.certificate.abstain if r is not None else "timeout"
            ch = characterize(Xd, yd, sigma=SIGMA_REP, abstain_reason=abstain)
            k = min(24, len(Xd))
            forms = propose_forms({"class": ch.get("class", ""),
                                   "why": ch.get("why", ""),
                                   "X": Xd[:k].tolist(), "y": yd[:k].tolist()},
                                  k=3) or []
            for form in forms:
                v = verify(Xd[np.isfinite(yd)].tolist(),
                           yd[np.isfinite(yd)].tolist(), form, sigma=SIGMA_REP)
                if v.get("certified"):
                    rec.update(track="certified", expr=v["law"],
                               tag="empirical-structural", channel="llm-verified",
                               detail="LLM proposal verified by the sound checker "
                                      f"at sigma_rep={SIGMA_REP}")
                    return rec
        except TimeoutError:
            rec["timed_out"] = True
        except Exception:                                      # noqa: BLE001
            pass
        finally:
            signal.alarm(0)
        if forms:
            rec.update(track="conjecture", expr=forms[0], tag="empirical",
                       channel="llm-conjecture",
                       detail="unverified LLM proposal (best of k)")
            return rec

        # ---- stage 3: always-reachable fallback ----
        conj = _fallback_conjecture(Xd, yd) or \
            (r.result.certificate.law if r is not None else None) or None
        if conj:
            rec.update(track="conjecture", expr=conj, tag="empirical",
                       channel="fit-probe", detail="log-log fallback conjecture")
        else:
            rec.update(track="abstain", expr=None, tag="open", channel="none",
                       detail=(r.result.certificate.abstain if r is not None
                               else "timeout") or "structural")
        return rec
    except Exception as e:                                     # noqa: BLE001
        rec.update(track="abstain", expr=None, tag="open", channel="error",
                   detail=str(e)[:100])
        return rec
    finally:
        signal.alarm(0)
        rec["secs"] = round(time.time() - t0, 1)


# ---------------------------------------------------------------- judge (fixed)

def _gt_expr_dev(problem):
    from sympy.core.function import AppliedUndef
    eq = problem.gt_equation
    # parse with the problem's OWN symbol names bound to Symbols, so names like
    # `gamma`/`beta` do not resolve to sympy FUNCTIONS (measured crash: gamma-1)
    names = {str(s) for s in eq.symbols}
    import re
    expr_s = eq.expression
    keep = {"sqrt", "exp", "log", "sin", "cos", "tan", "asin", "acos", "atan",
            "sinh", "cosh", "tanh", "pi", "E", "Abs"}
    # benchmark data artifact: some expression strings carry a stray `_f` suffix
    # on numeric literals (`0.325..._f`) -- strip it before tokenizing
    expr_s = re.sub(r"(\d)_f\b", r"\1", expr_s)
    # functional notation: A(t) / P(t) / x(t) mean "the observable named A/P/x"
    # (the dataset column), not a call -- collapse f(arg) -> f for non-math names
    expr_s = re.sub(
        r"\b([A-Za-z_]\w*)\s*\(\s*[A-Za-z_]\w*\s*\)",
        lambda m: m.group(1) if m.group(1) not in keep else m.group(0), expr_s)
    names |= set(re.findall(r"[A-Za-z_][A-Za-z_0-9]*", expr_s))
    local = {n: sp.Symbol(n) for n in names if n not in keep}
    e = sp.sympify(expr_s, locals=local)
    e = e.replace(lambda x: isinstance(x, AppliedUndef),
                  lambda x: sp.Symbol(x.func.__name__))
    subs = {sp.Symbol(str(s)): sp.Symbol(f"x_{i}")
            for i, s in enumerate(eq.symbols[1:])}
    return e.subs(subs)


def _structure(e):
    sig = set()
    e = sp.expand(e)
    for t in (e.as_ordered_terms() if e.is_Add else [e]):
        c, rest = t.as_coeff_Mul()
        sig.add(sp.srepr(sp.expand(rest)))
    return frozenset(sig)


def score(rows, problems):
    from experiments.run_blind_llmsrbench import _tests
    out = []
    for rec in rows:
        p = problems[rec["id"]]
        r = dict(rec)
        r["sa"] = False
        if rec.get("expr"):
            try:
                pred = sp.sympify(rec["expr"])
                gt = _gt_expr_dev(p)
                r["sa"] = _structure(pred) == _structure(gt) or \
                    bool(sp.simplify(sp.expand(pred - gt)) == 0)
            except Exception as ex:                            # noqa: BLE001
                r["judge_error"] = str(ex)[:80]
            syms = [sp.Symbol(f"x_{i}") for i in range(rec["dim"])]
            for split, arr in _tests(p).items():
                try:
                    Xt, yt = arr[:, 1:], arr[:, 0]
                    f = sp.lambdify(syms, sp.sympify(rec["expr"]), "numpy")
                    pr = np.broadcast_to(np.asarray(f(*Xt.T), float), yt.shape)
                    mm = np.isfinite(pr) & np.isfinite(yt)
                    if mm.sum() < max(8, 0.5 * len(yt)):
                        r[f"acc01_{split}"] = False
                        continue
                    rel = np.abs(pr[mm] - yt[mm]) / np.maximum(np.abs(yt[mm]), 1e-30)
                    r[f"acc01_{split}"] = bool(np.max(rel) <= 0.1)
                except Exception:                              # noqa: BLE001
                    r[f"acc01_{split}"] = False
        out.append(r)
    return out


def summarize(rows):
    def rate(rs, key):
        return round(100 * sum(1 for x in rs if x.get(key)) / len(rs), 2) if rs else 0.0
    cats = sorted({r["module"] for r in rows})
    summary = {}
    for name, grp in [("ALL", rows)] + [(c, [r for r in rows if r["module"] == c])
                                        for c in cats]:
        a = [r for r in grp if r["track"] == "certified"]
        summary[name] = {
            "n": len(grp), "SA": rate(grp, "sa"), "acc01_id": rate(grp, "acc01_id"),
            "certified": {"n": len(a), "SA_correct": sum(1 for r in a if r["sa"]),
                          "structural_wrong": sum(1 for r in a if not r["sa"]),
                          "llm_verified": sum(1 for r in a
                                              if r["channel"] == "llm-verified")},
            "conjecture": {"n": sum(1 for r in grp if r["track"] == "conjecture"),
                           "SA_correct": sum(1 for r in grp
                                             if r["track"] == "conjecture" and r["sa"])},
            "abstain": sum(1 for r in grp if r["track"] == "abstain"),
            "timeouts": sum(1 for r in grp if r.get("timed_out")),
        }
    return summary


def main():
    from experiments.run_blind_llmsrbench import load_all_problems, _train_xy
    from machine.llm import _load_env
    _load_env()
    problems = {}
    jobs = []
    for mod_name, p in load_all_problems():
        pid = f"{mod_name}/{p.equation_idx}"
        problems[pid] = p
        X, y = _train_xy(p)
        jobs.append((pid, X, y))
    done = {}
    if OUT.exists():
        for line in OUT.read_text().splitlines():
            r = json.loads(line)
            done[r["id"]] = r
    todo = [j for j in jobs if j[0] not in done]
    print(f"{len(jobs)} problems, {len(done)} done, {len(todo)} to run "
          f"({WORKERS} workers)", flush=True)
    import multiprocessing as mp
    ctx = mp.get_context("fork")
    with OUT.open("a") as fh, ctx.Pool(WORKERS) as pool:
        for rec in pool.imap_unordered(solve_one, todo, chunksize=1):
            rec["module"] = rec["id"].split("/")[0]
            fh.write(json.dumps(rec) + "\n")
            fh.flush()
            print(f"{rec['id']:55s} {rec['track']:10s} {rec.get('channel',''):14s} "
                  f"{rec.get('secs','?')}s", flush=True)
    rows = [json.loads(line) for line in OUT.read_text().splitlines()]
    for r in rows:
        r.setdefault("module", r["id"].split("/")[0])
    scored = score(rows, problems)
    summary = summarize(scored)
    SCORE_OUT.write_text(json.dumps({"summary": summary, "rows": scored}, indent=1))
    print(json.dumps(summary, indent=1))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

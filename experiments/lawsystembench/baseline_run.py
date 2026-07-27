"""LawSystemBench v1 baseline: the lagh system discoverer over all 80 problems
(8-way parallel), scored against ground truth: equations via the judge-v4
structural comparator, invariants via affine constant-equivalence on the data,
shared-constant verdict counts (the kill-criterion measurement).
Predictions P1-P4 in docs/LAWSYSTEMBENCH.md are frozen before this runs.
"""
from __future__ import annotations

import json
import signal
import sys
import time
from pathlib import Path

import numpy as np
import sympy as sp

sys.path.insert(0, str(Path(__file__).parent.parent.parent))
sys.path.insert(0, "/home/allans/code/llm-srbench")   # judge helpers only
from experiments.rejudge_dev_llmsrbench import canon, sig, strip_consts  # noqa: E402
from lagh.systems import discover_system  # noqa: E402

IN = Path("experiments/results/lawsystembench_v1.jsonl")
OUT = Path("experiments/results/lawsystembench_v1_baseline.jsonl")


def _h(s, f):
    raise TimeoutError()


def eq_match(pred, gt):
    try:
        return sig(strip_consts(canon(pred))) == sig(strip_consts(canon(gt)))
    except Exception:                                          # noqa: BLE001
        return False


def inv_match(found_exprs, gt_inv, columns):
    """A gt invariant counts as recovered when some found invariant is
    affinely equivalent to it on the data (same conserved level sets)."""
    gt_e = sp.sympify(gt_inv)
    syms = sorted(gt_e.free_symbols, key=str)
    vals_gt = np.asarray(sp.lambdify(syms, gt_e, "numpy")(
        *[np.asarray(columns[str(s)]) for s in syms]), float)
    for f in found_exprs:
        try:
            fe = sp.sympify(f)
            fs = sorted(fe.free_symbols, key=str)
            vals_f = np.asarray(sp.lambdify(fs, fe, "numpy")(
                *[np.asarray(columns[str(s)]) for s in fs]), float)
            A = np.column_stack([np.ones(len(vals_gt)), vals_f])
            c, *_ = np.linalg.lstsq(A, vals_gt, rcond=None)
            pred = A @ c
            scale = float(np.std(vals_gt)) + float(np.abs(vals_gt).mean()) + 1e-12
            if np.max(np.abs(pred - vals_gt)) < 1e-4 * scale and abs(c[1]) > 1e-8:
                return True
        except Exception:                                      # noqa: BLE001
            continue
    return False


def solve_one(p):
    signal.signal(signal.SIGALRM, _h)
    t0 = time.time()
    data = {k: np.asarray(v, float) for k, v in p["columns"].items()}
    rec = {"id": p["id"], "family": p["family"], "tier": p["tier"]}
    try:
        signal.alarm(600)
        cert = discover_system(data, sigma=p["sigma"])
        inv_map = {v: k for k, v in p["mapping"].items()}       # c_i -> name
        sub = {sp.Symbol(c): sp.Symbol(orig) for c, orig in inv_map.items()}

        def rename(e):
            return str(sp.sympify(e).xreplace(sub))

        eqs = {inv_map[t]: rename(e["expr"]) for t, e in cert.equations.items()}
        invs = [rename(iv["expr"]) for iv in cert.invariants]
        gt = p["gt_equations"]
        eq_ok = {t: (t in eqs and eq_match(eqs[t], gt[t])) for t in gt}
        wrong = [t for t in eqs if t in gt and not eq_ok[t]]
        extra = [t for t in eqs if t not in gt]
        cols_named = {inv_map[c]: np.asarray(v, float)
                      for c, v in p["columns"].items()}
        inv_ok = [inv_match(invs, gi, cols_named) for gi in p["gt_invariants"]]
        rec.update(
            all_eqs_ok=all(eq_ok.values()), eq_ok=sum(eq_ok.values()),
            eq_total=len(gt), eq_wrong=len(wrong), eq_extra=len(extra),
            inv_found=len(invs), inv_recovered=sum(inv_ok),
            inv_total=len(p["gt_invariants"]),
            shared_verdicts=len(cert.shared),
            alpha=cert.alpha_log10_total)
    except TimeoutError:
        rec["error"] = "timeout"
    except Exception as e:                                     # noqa: BLE001
        rec["error"] = str(e)[:120]
    finally:
        signal.alarm(0)
    rec["secs"] = round(time.time() - t0, 1)
    return rec


def main():
    problems = [json.loads(l) for l in IN.read_text().splitlines()]
    done = set()
    if OUT.exists():
        done = {json.loads(l)["id"] for l in OUT.read_text().splitlines()}
    todo = [p for p in problems if p["id"] not in done]
    print(f"{len(problems)} problems, {len(todo)} to run", flush=True)
    import multiprocessing as mp
    ctx = mp.get_context("fork")
    with OUT.open("a") as fh, ctx.Pool(8) as pool:
        for rec in pool.imap_unordered(solve_one, todo, chunksize=1):
            fh.write(json.dumps(rec) + "\n")
            fh.flush()
            print(f"{rec['id']:26s} eqs {rec.get('eq_ok','-')}/{rec.get('eq_total','-')}"
                  f" wrong={rec.get('eq_wrong','-')} inv {rec.get('inv_recovered','-')}"
                  f"/{rec.get('inv_total','-')} shared={rec.get('shared_verdicts','-')}"
                  f" {rec.get('secs','')}s {rec.get('error','')}", flush=True)

    rows = [json.loads(l) for l in OUT.read_text().splitlines()]
    ok = [r for r in rows if "error" not in r]
    print("\n== LAWSYSTEMBENCH v1 BASELINE ==")
    for name in ("clean", "noisy"):
        grp = [r for r in ok if r["tier"] == name]
        if not grp:
            continue
        print(f"{name}: n={len(grp)} all-eqs-correct={sum(r['all_eqs_ok'] for r in grp)}"
              f" eq-wrong-total={sum(r['eq_wrong'] for r in grp)}"
              f" invariants {sum(r['inv_recovered'] for r in grp)}"
              f"/{sum(r['inv_total'] for r in grp)}"
              f" shared-verdicts={sum(r['shared_verdicts'] for r in grp)}")
    print("errors:", sum(1 for r in rows if "error" in r))


if __name__ == "__main__":
    main()

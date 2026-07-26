"""Judge v2 re-score of dev sweep v1.1 (no re-discovery, no LLM).

The conservative v1 judge under-counted: (a) split vs fused radicals differ in
srepr on positive domains; (b) noise-limited float coefficients defeat exact
simplify(pred-gt)==0; (c) inner constants (exp(2.31*x) vs exp(c*x)) were never
stripped. v2 = canonicalize under positive symbols + strip ALL fitted-scale
constants + three equivalence probes. Still a LOWER bound on the official
GPT-4o judge; both judges' numbers are reported side by side.
"""
from __future__ import annotations

import json
import signal
import sys
from pathlib import Path

import numpy as np
import sympy as sp

sys.path.insert(0, str(Path(__file__).parent.parent))
sys.path.insert(0, "/home/allans/code/llm-srbench")

import os
_TAG = os.environ.get("DEV_SWEEP_TAG", "v1")
IN = Path(f"experiments/results/dev_llmsrbench_{_TAG}_scores.json")
OUT = Path(f"experiments/results/dev_llmsrbench_{_TAG}_scores_v2.json")

POS = {sp.Symbol(f"x_{i}"): sp.Symbol(f"x_{i}", positive=True) for i in range(12)}


def canon(e):
    e = sp.sympify(e).xreplace(POS)
    # gt strings carry FLOAT exponents (x**0.333333333333333); rationalize the
    # near-rational ones so they compare equal to a grammar's exact Rational(1,3)
    # and so strip_consts never mistakes an exponent for a coefficient
    def _rx(x):
        fr = sp.nsimplify(x.exp, rational=True, tolerance=1e-6)
        if fr.is_Rational and fr.q <= 24 and \
                abs(float(fr) - float(x.exp)) < 1e-6:
            return sp.Pow(x.base, fr)
        return x
    e = e.replace(lambda x: isinstance(x, sp.Pow) and isinstance(x.exp, sp.Float),
                  _rx)
    return sp.powsimp(sp.together(sp.expand(e)), force=True)


def strip_consts(e):
    """Replace every fitted-scale numeric atom by 1: Floats always; rationals and
    integers with any component > 1000 (noise-limited snaps); small rationals
    (exponents like 3/2) survive."""
    subs = {}
    for a in e.atoms(sp.Number):
        if isinstance(a, sp.Float):
            subs[a] = sp.Integer(1)
        elif isinstance(a, sp.Rational) and (abs(a.p) > 1000 or a.q > 1000):
            subs[a] = sp.Integer(1)
    return e.xreplace(subs)


def sig(e):
    """Per-additive-term signature with the ENTIRE numeric-valued factor stripped
    (floats, rationals, pi, sqrt(6), E, ...) -- the official protocol's 'remove
    parameters and constants'. as_independent splits on free symbols, so 1/(8*pi^2)
    and 0.0126651... strip identically."""
    e = sp.expand(e)
    out = set()
    for t in (e.as_ordered_terms() if e.is_Add else [e]):
        if not t.free_symbols:
            out.add("CONST")
            continue
        _, rest = t.as_independent(*t.free_symbols, as_Add=False)
        out.add(sp.srepr(sp.powsimp(sp.expand(rest), force=True)))
    return frozenset(out)


class _TO(Exception):
    pass


def _h(s, f):
    raise _TO()


def equivalent(pred_s, gt):
    prev = signal.getsignal(signal.SIGALRM)   # restore the caller's handler after
    signal.signal(signal.SIGALRM, _h)
    try:
        signal.alarm(30)
        pred = canon(pred_s)
        gtc = canon(gt)
        if sig(strip_consts(pred)) == sig(strip_consts(gtc)):
            return True, "structure"
        q = sp.simplify(pred / gtc)
        if q.is_Number and q != 0:
            return True, "ratio-constant"
        d = sp.simplify(pred - gtc)
        if d.is_Number:
            return True, "diff-constant"
        return False, ""
    except (_TO, Exception):                                   # noqa: BLE001
        return False, "judge-v2-error"
    finally:
        signal.alarm(0)
        signal.signal(signal.SIGALRM, prev)


def main():
    from experiments.run_blind_llmsrbench import load_all_problems
    from experiments.run_dev_llmsrbench import _gt_expr_dev, summarize
    problems = {f"{m}/{p.equation_idx}": p for m, p in load_all_problems()}
    rows = json.load(open(IN))["rows"]
    changed = 0
    for r in rows:
        r["sa_v1"] = r.get("sa", False)
        if r.get("expr") and not r.get("sa"):
            try:
                gt = _gt_expr_dev(problems[r["id"]])
            except Exception:                                  # noqa: BLE001
                continue
            ok, how = equivalent(r["expr"], gt)
            if ok:
                r["sa"] = True
                r["sa_v2_how"] = how
                changed += 1
    summary = summarize(rows)
    OUT.write_text(json.dumps({"summary": summary, "rows": rows}, indent=1))
    print(f"judge v2 upgraded {changed} rows")
    print(json.dumps(summary, indent=1))


if __name__ == "__main__":
    raise SystemExit(main())

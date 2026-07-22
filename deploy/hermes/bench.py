#!/usr/bin/env python3
"""Full-run benchmark driver + scorer for MiniMax + Hermes + lagh.

Runs INSIDE the container. For each cell: build a per-cell prompt (fresh context),
invoke `hermes chat -q`, extract the model's strict RESULT line, then score the reported
law against the REAL oracle (dense comparison) and against lagh-alone (the direct-sweep
baseline). Prints the composite scorecard.

Usage (via run_bench.sh, or directly in the container):
    python deploy/hermes/bench.py --subset easy          # all 36 easy cells
    python deploy/hermes/bench.py --subset all            # 108
    python deploy/hermes/bench.py --cells m0_gravity/easy/v0,m4_snell_law/easy/v0
    python deploy/hermes/bench.py --subset easy --limit 6 --timeout 240
"""
from __future__ import annotations

import argparse
import json
import os
import re
import subprocess
import sys
import time

import numpy as np
import sympy as sp
from sympy.parsing.sympy_parser import (convert_xor,                  # noqa: E402
                                        implicit_multiplication_application,
                                        parse_expr, standard_transformations)

_TF = standard_transformations + (implicit_multiplication_application, convert_xor)

os.environ.setdefault("LAB_SOURCE", "newtonbench")
from lagh.lab import get_source                                       # noqa: E402
from lagh.adapters.newtonbench import MODULES, make_oracle           # noqa: E402

PROMPT = """You are discovering an exact physical law by querying an oracle with certified tools. Work ONLY from the samples -- the hint may be a distractor and the true law is often a VARIANT of the textbook one.

Problem id: {pid}
Inputs (use this exact column order for every row of X): {inputs}
Suggested domain: lo={lo}, hi={hi}
Hint (physical context, may mislead): {hint}

Tools available:
- lab.sample(problem, X): query the oracle. problem={pid!r}. X is a list of rows; each row has {n} number(s) in the input column order above. Returns y for each row.
- lagh.fit(X, y): scout -- conjectures + a diagnosis (pinned / continuum / acquire_more_data). A guess, NOT a certificate.
- lagh.recover(X, y): discover an exact law -> a certificate (certified true, law, strength) OR an abstention (certified false, abstain reason).
- lab.discover(problem): DELEGATE the whole thing to lagh -- it runs its OWN adaptive sampling loop against the oracle and certifies where a hand-picked sample often can't. Your best move when recover abstains.
- lagh.verify(X, y, form): check a form YOU declare -> certificate (pinned or consistent) or abstain.

Procedure:
1. Pick at least 12 points spread across the suggested domain (log-spaced if it spans decades); call lab.sample.
2. lagh.fit to scout. If it says acquire_more_data or continuum, sample a WIDER range (you may go outside the suggested domain) and try again.
3. lagh.recover on your best (X, y).
4. If recover ABSTAINS, call lab.discover(problem) -- lagh runs its full active-acquisition loop and usually certifies. Report ITS certificate (via="recover"). Do NOT give up before trying discover.
5. Only if discover ALSO abstains may you declare a form and lagh.verify it; if still nothing certifies, abstain honestly.

Honesty: report only what the tools certified; do not invent a law; do not assume the textbook form.

End with EXACTLY ONE final line and nothing after it. Write it as RAW TEXT -- no markdown, no bold/italics -- and copy the law VERBATIM from the tool with every * and ** operator intact (e.g. 3*x_0**2, NOT 3x_0^2 and NOT 3x_0*2):
RESULT: {{"status": "certified"|"consistent"|"abstained", "law": "<the exact sympy expr lagh returned in x_0.. variables, or null>", "via": "recover"|"verify"|"none"}}
"""


def cells(args) -> list[str]:
    ids = [c["id"] for c in get_source("newtonbench").problems()]
    if args.cells:
        want = set(args.cells.split(","))
        ids = [i for i in ids if i in want]
    elif args.subset != "all":
        ids = [i for i in ids if i.split("/")[1] == args.subset]
    if args.limit:
        ids = ids[: args.limit]
    return ids


def prompt_for(pid: str) -> str:
    m, d, v = pid.split("/")
    inp, lo, hi = MODULES[m]
    return PROMPT.format(pid=pid, inputs=list(inp), lo=list(map(float, lo)),
                         hi=list(map(float, hi)),
                         hint=f"{m.split('_', 1)[1].replace('_', ' ')} ({d})",
                         n=len(inp))


def run_cell(pid: str, timeout: int, trace: bool = False) -> tuple[str, dict | None]:
    """Invoke hermes headless for one cell; return (raw_output, parsed RESULT|None).
    trace=True uses `hermes chat -q` (full tool transcript, for debugging) instead of
    `hermes -z` (clean final answer, the default -- no markdown to mangle the law)."""
    cmd = (["hermes", "chat", "-q", prompt_for(pid)] if trace
           else ["hermes", "-z", prompt_for(pid)])
    try:
        p = subprocess.run(cmd, capture_output=True, text=True, timeout=timeout)
        raw = (p.stdout or "") + (p.stderr or "")
    except subprocess.TimeoutExpired:
        return "TIMEOUT", None
    except Exception as e:                                            # noqa: BLE001
        return f"ERROR: {e}", None
    m = list(re.finditer(r"RESULT:\s*(\{.*?\})", raw, re.S))
    if not m:
        return raw, None
    try:
        return raw, json.loads(m[-1].group(1))
    except Exception:                                                # noqa: BLE001
        return raw, None


def _to_x(expr, dim):
    """Map any input-name symbols in the reported law onto x_0..x_{dim-1} by position."""
    subs = {}
    for m in MODULES.values():
        for i, name in enumerate(m[0]):
            if i < dim:
                subs[sp.Symbol(name)] = sp.Symbol(f"x_{i}")
    return expr.xreplace(subs)


def score_law(pid: str, law: str):
    """Dense-compare the reported law to the real oracle. Returns (correct, rel_err),
    or None if the law can't be parsed/evaluated at all -- a formatting/capture issue,
    NOT a wrong law (so it is never miscounted as a confident-wrong)."""
    m, d, v = pid.split("/")
    inp, lo, hi = MODULES[m]
    dim = len(inp)
    try:
        expr = _to_x(parse_expr(str(law), transformations=_TF), dim)  # implicit-mult tolerant
    except Exception:                                                # noqa: BLE001
        return None
    orc = make_oracle(m, v, d)
    rng = np.random.default_rng(20260722)
    X = np.exp(rng.uniform(np.log(np.maximum(np.array(lo, float), 1e-9)),
                           np.log(np.array(hi, float)), (200, dim)))
    y = np.asarray(orc(X), float)
    try:
        f = sp.lambdify([sp.Symbol(f"x_{i}") for i in range(dim)], expr, "numpy")
        got = np.broadcast_to(np.asarray(f(*X.T), float), y.shape)
    except Exception:                                                # noqa: BLE001
        return None
    ok = np.isfinite(got) & np.isfinite(y) & (np.abs(y) > 1e-9)
    if ok.sum() < 20:
        return None
    re_ = float(np.median(np.abs(got[ok] - y[ok]) / np.abs(y[ok])))
    return re_ < 1e-2, re_


def lagh_alone_map() -> dict:
    path = os.path.join(os.path.dirname(__file__), "..", "..",
                        "experiments", "results", "newtonbench_all.jsonl")
    out = {}
    if os.path.exists(path):
        for ln in open(path):
            r = json.loads(ln)
            out[f"{r['module']}/{r['difficulty']}/{r['version']}"] = bool(
                r["certified"] and r["correct"])
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--subset", default="easy", help="easy|medium|hard|all")
    ap.add_argument("--cells", default="")
    ap.add_argument("--limit", type=int, default=0)
    ap.add_argument("--timeout", type=int, default=240)
    ap.add_argument("--trace", action="store_true",
                    help="save the FULL tool transcript (hermes chat -q) instead of just "
                         "the final answer -- use to see whether the model called discover")
    ap.add_argument("--out", default=os.path.join(os.path.dirname(__file__), "bench_results"))
    args = ap.parse_args()

    os.makedirs(args.out, exist_ok=True)
    ids = cells(args)
    base = lagh_alone_map()
    print(f"running {len(ids)} cell(s) [subset={args.subset}] through MiniMax+Hermes+lagh\n",
          flush=True)
    rows = []
    for k, pid in enumerate(ids, 1):
        t0 = time.time()
        raw, res = run_cell(pid, args.timeout, args.trace)
        open(os.path.join(args.out, pid.replace("/", "_") + ".txt"), "w").write(raw)
        status = (res or {}).get("status")
        law = (res or {}).get("law")
        via = (res or {}).get("via")
        if res is None:
            cat = "ERROR"
        elif status in ("certified", "consistent") and law and str(law).lower() != "null":
            sc = score_law(pid, law)
            cat = ("PARSE_ERROR" if sc is None
                   else ("CORRECT" if sc[0] else "CONFIDENT_WRONG"))
        else:  # abstained / no law
            cat = "MISS" if base.get(pid) else "HONEST_ABSTAIN"
        rows.append({"pid": pid, "cat": cat, "status": status, "via": via,
                     "law": law, "lagh_alone": base.get(pid), "secs": round(time.time() - t0, 1)})
        print(f"  [{k:>3}/{len(ids)}] {pid:<28} {cat:<16} via={via or '-':<7} "
              f"{(str(law)[:34] if law else '')}", flush=True)

    json.dump(rows, open(os.path.join(args.out, "results.json"), "w"), indent=1)
    from collections import Counter
    c = Counter(r["cat"] for r in rows)
    comp = c["CORRECT"]
    base_n = sum(1 for r in rows if r["lagh_alone"])
    gains = [r["pid"] for r in rows if r["cat"] == "CORRECT" and not r["lagh_alone"]]
    losses = [r["pid"] for r in rows if r["cat"] in ("MISS", "CONFIDENT_WRONG") and r["lagh_alone"]]
    print("\n" + "=" * 52)
    print(f"SCORECARD ({len(rows)} cells)")
    print("=" * 52)
    for k in ("CORRECT", "HONEST_ABSTAIN", "MISS", "CONFIDENT_WRONG", "PARSE_ERROR", "ERROR"):
        print(f"  {k:<16} {c[k]}")
    print(f"\n  composite (MiniMax+lagh) certified-correct: {comp}/{len(rows)}")
    print(f"  lagh-alone (direct sweep) on same cells:    {base_n}/{len(rows)}")
    print(f"  gains (composite got, lagh-alone missed):   {len(gains)}  {gains}")
    print(f"  losses (lagh-alone got, composite missed):  {len(losses)}  {losses}")
    print(f"\n  CONFIDENT-WRONG must be 0 -> {'OK' if c['CONFIDENT_WRONG'] == 0 else 'INVARIANT BREACH'}")
    return 0


if __name__ == "__main__":
    sys.exit(main())

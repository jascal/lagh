"""Registered escalation-rule study (docs/ESCALATION_REGISTRATION.md).

Scores the engine's escalation rule -- `first` (stop at the first tier with a
non-empty certifying set, the shipped default) against `pool` (judge the
highest tier's certifying set, which includes every lower tier's candidates) --
on three banks, recording CORRECTNESS of every certificate, not only whether
one was issued:

  p2     the 16 registered P2 draws of (3x+2)/(x+4) on [.5,3] (seeds 20-35)
  p1     the three distinct P1 seed10 datasets (rational-a/b/c on [.5,3];
         restricted-a shares rational-a's data and is not repeated)
  reach  the 36 cells of experiments/reach/audit.py, unchanged data

Each certified law is evaluated on 400 FRESH in-box points (the domain claim)
and 400 points of an extended box [.25, 6]^d (the form claim). One case per
process call; results are append-only JSON files that are never overwritten.
"""
from __future__ import annotations

import argparse
import json
import time
import zlib
from pathlib import Path

import numpy as np
import sympy as sp

from experiments.reach import audit
from lagh.base import eval_expr
from lagh.passive import discover_passive

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'experiments/results/escalation'
RULES = ('first', 'pool', 'accumulate')
WRONG_REL = 1e-6        # registered: a certificate is wrong beyond this
EXT_BOX = (.25, 6.)
x = sp.Symbol('x_0')


def _rational(expr):
    f = sp.lambdify(x, expr, 'numpy')
    return (lambda X: f(X[:, 0])), expr


def cases():
    out = {}
    for seed in range(20, 36):
        fn, truth = _rational((3*x+2)/(x+4))
        X = np.random.default_rng(seed).uniform(.5, 3., (400, 1))
        out[f'p2-seed{seed}'] = (X, fn, truth)
    for name, truth in (('rational-a', (2*x+1)/(x+3)),
                        ('rational-b', (3*x+2)/(x+4)),
                        ('rational-c', (x+1)/(x+5))):
        fn, _ = _rational(truth)
        out[f'p1-{name}'] = (np.random.default_rng(10).uniform(.5, 3., (400, 1)),
                             fn, truth)
    for name, dim, fn in audit.CELLS:
        out[f'reach-{name}'] = (audit.X(dim, seed=zlib.crc32(name.encode())),
                                fn, None)
    return out


def _rel_error(expr, fn, P):
    syms = [sp.Symbol(f'x_{j}') for j in range(P.shape[1])]
    pred, target = eval_expr(expr, syms, P), fn(P)
    if pred is None or not np.all(np.isfinite(pred)):
        return None
    return float(np.max(np.abs(pred - target)) / np.max(np.abs(target)))


def score(expr, fn, truth, dim, seed):
    rng = np.random.default_rng(seed + 1)
    inbox = _rel_error(expr, fn, rng.uniform(.5, 3., (400, dim)))
    ext = _rel_error(expr, fn, rng.uniform(*EXT_BOX, (400, dim)))
    exact = None if truth is None else bool(sp.cancel(expr - truth) == 0)
    wrong_domain = inbox is None or inbox > WRONG_REL
    wrong_form = exact is False if exact is not None else (
        ext is None or ext > WRONG_REL)
    return {'inbox_rel_error': inbox, 'ext_rel_error': ext, 'exact': exact,
            'wrong_domain': wrong_domain, 'wrong_form': bool(wrong_form)}


def run(name, rule):
    path = OUT / rule / f'{name}.json'
    if path.exists():
        print(f'{rule} {name}: retained', flush=True)
        return
    X, fn, truth = cases()[name]
    y = fn(X)
    t0 = time.time()
    r = discover_passive(X, y, sigma=0.0, escalation=rule)
    record = {'tag': 'empirical', 'case': name, 'rule': rule,
              'certified': bool(r.certified), 'tier': r.result.tier,
              'abstain': r.result.certificate.abstain,
              'law': str(r.result.expr) if r.certified else None,
              'rivals': [str(z) for z in r.result.rivals],
              'seconds': round(time.time() - t0, 1)}
    if r.certified:
        record.update(score(r.result.expr, fn, truth, X.shape[1],
                            zlib.crc32(name.encode())))
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(record, indent=1, allow_nan=False) + '\n')
    print(f"{rule} {name}: certified={record['certified']} "
          f"wrong_form={record.get('wrong_form')} {record['seconds']}s", flush=True)


def summarize():
    rows = [json.loads(p.read_text()) for p in sorted(OUT.glob('*/*.json'))]
    summary = {}
    for rule in RULES:
        mine = [r for r in rows if r['rule'] == rule]
        for bank in ('p2', 'p1', 'reach'):
            b = [r for r in mine if r['case'].startswith(bank + '-')]
            cert = [r for r in b if r['certified']]
            summary[f'{rule}/{bank}'] = {
                'cases': len(b), 'certified': len(cert),
                'wrong_domain': sorted(r['case'] for r in cert if r['wrong_domain']),
                'wrong_form': sorted(r['case'] for r in cert if r['wrong_form']),
                'median_seconds': float(np.median([r['seconds'] for r in b])) if b else None}
    (OUT / 'summary.json').write_text(json.dumps(summary, indent=1) + '\n')
    print(json.dumps(summary, indent=1))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--list', action='store_true')
    ap.add_argument('--case')
    ap.add_argument('--rule', choices=RULES)
    ap.add_argument('--summarize', action='store_true')
    a = ap.parse_args()
    if a.list:
        print('\n'.join(cases()))
    elif a.summarize:
        summarize()
    else:
        run(a.case, a.rule)


if __name__ == '__main__':
    main()

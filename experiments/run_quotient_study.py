"""Registered A6 study: the joint gate in the exact quotient space
(docs/ESCALATION_REGISTRATION.md, A6).

The constrained-input bank `cs`: inputs on a variety with exact constraints,
where `joint_modulo` (A4: reduce, then `reduce_to_minimal` on all rows) and
`joint_quotient` (A6: exclude the ideal directions exactly inside the joint
test) can differ. Four varieties x six laws, 400 rows each, three gates.

Scoring is on the variety: `inbox` = 400 fresh points from the sampling
region (the domain claim), `ext` = 400 fresh points from a wider region of the
SAME variety (the form claim). Off-variety points are never used: a
domain-restricted certificate says nothing there.
"""
from __future__ import annotations

import argparse
import json
import time
import zlib
from pathlib import Path

import numpy as np
import sympy as sp

from lagh.base import eval_expr
from lagh.passive import discover_passive

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'experiments/results/joint_quotient'
GATES = ('marginal', 'joint_modulo', 'joint_quotient')
WRONG_REL = 1e-6


def _sphere(rng, n, wide):
    v = rng.normal(size=(n, 3))
    if not wide:
        v[:, 2] = np.abs(v[:, 2])               # sampled on the upper hemisphere
    return v / np.linalg.norm(v, axis=1, keepdims=True)


def _circle(rng, n, wide):
    t = rng.uniform(0, 2 * np.pi if wide else np.pi / 2, n)
    return np.column_stack([np.cos(t), np.sin(t)])


def _hyperbola(rng, n, wide):
    x = rng.uniform(*((.25, 6.) if wide else (.5, 3.)), n)
    return np.column_stack([x, 1 / x])


def _plane(rng, n, wide):
    lo, hi = (-.5, 1.5) if wide else (.2, 1.)
    a, b = rng.uniform(lo, hi, n), rng.uniform(lo, hi, n)
    return np.column_stack([a, b, 1 - a - b])


VARIETIES = {'sphere': (_sphere, 3), 'circle': (_circle, 2),
             'hyperbola': (_hyperbola, 2), 'plane': (_plane, 3)}
LAWS = {
    'float-linear': lambda X: X @ np.array([0.4559837762, -0.8676661490,
                                            -0.1980763734][:X.shape[1]]),
    'int-linear': lambda X: 2 * X[:, 0] - X[:, 1] + (X[:, 2] if X.shape[1] > 2 else 0),
    'rational': lambda X: 1 / (2 + X[:, 0]),
    'product': lambda X: X[:, 0] * X[:, 1] + 3,
    'exp-bait': lambda X: np.exp(X[:, 0] / 2),
    'sqrt-bait': lambda X: np.sqrt(3 + X[:, 1]),
}


def cases():
    return [f'cs-{v}-{l}' for v in VARIETIES for l in LAWS]


def _data(name, n, wide, salt=0):
    _, v, *rest = name.split('-', 2)
    law = rest[0]
    gen, _dim = VARIETIES[v]
    rng = np.random.default_rng(zlib.crc32(name.encode()) + salt)
    X = gen(rng, n, wide)
    return X, LAWS[law](X)


def _rel(expr, X, y):
    syms = [sp.Symbol(f'x_{j}') for j in range(X.shape[1])]
    p = eval_expr(expr, syms, X)
    if p is None or not np.all(np.isfinite(p)):
        return None
    return float(np.max(np.abs(p - y)) / np.max(np.abs(y)))


def run(name, gate):
    path = OUT / gate / f'{name}.json'
    if path.exists():
        print(f'{gate} {name}: retained', flush=True)
        return
    X, y = _data(name, 400, wide=False)
    t0 = time.time()
    s0 = time.clock_gettime(time.CLOCK_BOOTTIME) - time.monotonic()
    r = discover_passive(X, y, sigma=0.0, coefficient_gate=gate)
    rec = {'tag': 'empirical', 'case': name, 'gate': gate,
           'certified': bool(r.certified), 'tier': r.result.tier,
           'abstain': r.result.certificate.abstain,
           'law': str(r.result.expr) if r.certified else None,
           'notes': [str(n)[:200] for n in r.result.certificate.notes][:3],
           'seconds': round(time.time() - t0, 1),
           'suspended_seconds': round(time.clock_gettime(time.CLOCK_BOOTTIME)
                                      - time.monotonic() - s0, 1)}
    if r.certified:
        Xi, yi = _data(name, 400, wide=False, salt=1)
        Xe, ye = _data(name, 400, wide=True, salt=2)
        rec['inbox_rel_error'] = _rel(r.result.expr, Xi, yi)
        rec['ext_rel_error'] = _rel(r.result.expr, Xe, ye)
        rec['wrong_domain'] = rec['inbox_rel_error'] is None or rec['inbox_rel_error'] > WRONG_REL
        rec['wrong_form'] = rec['ext_rel_error'] is None or rec['ext_rel_error'] > WRONG_REL
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(rec, indent=1, allow_nan=False) + '\n')
    print(f"{gate} {name}: certified={rec['certified']} "
          f"wrong={rec.get('wrong_form') or rec.get('wrong_domain')} {rec['seconds']}s",
          flush=True)


def summarize():
    rows = {g: {p.stem: json.loads(p.read_text()) for p in (OUT / g).glob('*.json')}
            for g in GATES}
    out = {}
    for g in GATES:
        cert = [r for r in rows[g].values() if r['certified']]
        out[g] = {'cases': len(rows[g]), 'certified': len(cert),
                  'wrong': sorted(r['case'] for r in cert
                                  if r['wrong_form'] or r['wrong_domain']),
                  'slept': sorted(r['case'] for r in rows[g].values()
                                  if r['suspended_seconds'] > 1),
                  'median_seconds': float(np.median([r['seconds'] for r in rows[g].values()]))
                  if rows[g] else None}
    q, m = rows['joint_quotient'], rows['joint_modulo']
    keys = sorted(set(q) & set(m))
    out['quotient_vs_modulo'] = {
        'differences': [{'case': k, 'joint_modulo': (m[k]['certified'], m[k]['law']),
                         'joint_quotient': (q[k]['certified'], q[k]['law'])}
                        for k in keys if (m[k]['certified'], m[k]['law'])
                        != (q[k]['certified'], q[k]['law'])],
        'median_ratio': float(np.median([q[k]['seconds'] / max(m[k]['seconds'], .1)
                                         for k in keys])) if keys else None}
    (OUT / 'summary.json').write_text(json.dumps(out, indent=1) + '\n')
    print(json.dumps(out, indent=1))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--list', action='store_true')
    ap.add_argument('--case')
    ap.add_argument('--gate', choices=GATES)
    ap.add_argument('--summarize', action='store_true')
    a = ap.parse_args()
    if a.list:
        print('\n'.join(cases()))
    elif a.summarize:
        summarize()
    else:
        run(a.case, a.gate)


if __name__ == '__main__':
    main()

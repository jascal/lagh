"""Registered default-switch study, FRESH banks (docs/ESCALATION_REGISTRATION.md, A3).

Compares the shipped default coefficient gate (`marginal`) against `joint`, both
under the default `first` escalation, on data that no earlier study has seen:

  fr    24 draws of (a*x+b)/(x+c) on [.5,3], 400 rows, seeds 100-123; a, b, c
        are integers 1..5 from rng(seed), redrawn deterministically when
        a*c == b (the law would collapse to the constant a)
  frch  the 36 reach cells re-drawn with seed crc32(name)+1

Scoring is `run_escalation_study.score`: fresh in-box points (domain claim) and
the extended box [.25, 6]^d (form claim). Results are append-only.
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
from experiments.run_escalation_study import _rational, score, x
from lagh.passive import discover_passive

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'experiments/results/default_switch'
GATES = ('marginal', 'joint')


def _fresh_rational(seed):
    rng = np.random.default_rng(seed)
    while True:
        a, b, c = (int(v) for v in rng.integers(1, 6, 3))
        if a * c != b:                      # a*c == b collapses to the constant a
            return (a*x + b)/(x + c)


def cases():
    out = {}
    for seed in range(100, 124):
        fn, truth = _rational(_fresh_rational(seed))
        out[f'fr-seed{seed}'] = (np.random.default_rng(seed).uniform(.5, 3., (400, 1)),
                                 fn, truth)
    for name, dim, fn in audit.CELLS:
        out[f'frch-{name}'] = (audit.X(dim, seed=zlib.crc32(name.encode()) + 1),
                               fn, None)
    return out


def run(name, gate):
    path = OUT / gate / f'{name}.json'
    if path.exists():
        print(f'{gate} {name}: retained', flush=True)
        return
    X, fn, truth = cases()[name]
    t0 = time.time()
    r = discover_passive(X, fn(X), sigma=0.0, coefficient_gate=gate)
    rec = {'tag': 'empirical', 'case': name, 'gate': gate,
           'truth': None if truth is None else str(truth),
           'certified': bool(r.certified), 'tier': r.result.tier,
           'abstain': r.result.certificate.abstain,
           'law': str(r.result.expr) if r.certified else None,
           'rivals': [str(z) for z in r.result.rivals],
           'seconds': round(time.time() - t0, 1)}
    if r.certified:
        rec.update(score(r.result.expr, fn, truth, X.shape[1],
                         zlib.crc32(name.encode())))
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(rec, indent=1, allow_nan=False) + '\n')
    print(f"{gate} {name}: certified={rec['certified']} "
          f"wrong={rec.get('wrong_form') or rec.get('wrong_domain')} {rec['seconds']}s",
          flush=True)


def summarize():
    rows = {g: {p.stem: json.loads(p.read_text()) for p in (OUT / g).glob('*.json')}
            for g in GATES}
    summary = {}
    for bank in ('fr', 'frch'):
        for g in GATES:
            b = {k: r for k, r in rows[g].items() if k.startswith(bank + '-')}
            cert = [r for r in b.values() if r['certified']]
            summary[f'{g}/{bank}'] = {
                'cases': len(b), 'certified': len(cert),
                'exact': sum(bool(r.get('exact')) for r in cert),
                'wrong': sorted(r['case'] for r in cert
                                if r['wrong_form'] or r['wrong_domain']),
                'median_seconds': float(np.median([r['seconds'] for r in b.values()]))
                if b else None}
        # D6: cost on the cases where the marginal default did NOT certify wrongly
        m, j = rows['marginal'], rows['joint']
        fair = [k for k in m if k.startswith(bank + '-') and k in j and not (
            m[k]['certified'] and (m[k]['wrong_form'] or m[k]['wrong_domain']))]
        ratios = [j[k]['seconds'] / max(m[k]['seconds'], 0.1) for k in fair]
        summary[f'cost/{bank}'] = {
            'fair_cases': len(fair),
            'median_ratio_joint_over_marginal': float(np.median(ratios)) if ratios else None,
            'lost_certificates': sorted(k for k in fair if m[k]['certified']
                                        and not j[k]['certified'])}
    (OUT / 'summary.json').write_text(json.dumps(summary, indent=1) + '\n')
    print(json.dumps(summary, indent=1))


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

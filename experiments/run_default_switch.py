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
from lagh.engine import discover
from lagh.passive import discover_passive

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'experiments/results/default_switch'
GATES = ('marginal', 'joint', 'joint_modulo')


def _fresh_rational(seed):
    rng = np.random.default_rng(seed)
    while True:
        a, b, c = (int(v) for v in rng.integers(1, 6, 3))
        if a * c != b:                      # a*c == b collapses to the constant a
            return (a*x + b)/(x + c)


def cases():
    out = {}
    # A5 fresh banks (registered before any A5 run)
    for seed in range(300, 324):
        fn, truth = _rational(_fresh_rational(seed))
        out[f'fr3-seed{seed}'] = (np.random.default_rng(seed).uniform(.5, 3., (400, 1)),
                                  fn, truth)
    for name, dim, fn in audit.CELLS:
        out[f'frch3-{name}'] = (audit.X(dim, seed=zlib.crc32(name.encode()) + 3),
                                fn, None)
    # A4 fresh banks (registered before any A4 run): new seeds, never run before
    for seed in range(200, 224):
        fn, truth = _rational(_fresh_rational(seed))
        out[f'fr2-seed{seed}'] = (np.random.default_rng(seed).uniform(.5, 3., (400, 1)),
                                  fn, truth)
    for name, dim, fn in audit.CELLS:
        out[f'frch2-{name}'] = (audit.X(dim, seed=zlib.crc32(name.encode()) + 2),
                                fn, None)
    for seed in range(100, 124):
        fn, truth = _rational(_fresh_rational(seed))
        out[f'fr-seed{seed}'] = (np.random.default_rng(seed).uniform(.5, 3., (400, 1)),
                                 fn, truth)
    for name, dim, fn in audit.CELLS:
        out[f'frch-{name}'] = (audit.X(dim, seed=zlib.crc32(name.encode()) + 1),
                               fn, None)
    return out


def split_verdicts(X, fn, truth, gate, name):
    """A5: the per-split verdicts of the passive run, recomputed OUTSIDE the
    timed call with passive's own split procedure (seed 0 + k, 60/20/20).
    A split that certifies a wrong law is a wrong verdict even when the
    passive full-data gate later rejects it."""
    y = fn(X)
    out = []
    for k in range(3):
        idx = np.random.default_rng(k).permutation(len(X))
        a, b = int(0.6 * len(X)), int(0.8 * len(X))
        r = discover(X[idx[:a]], y[idx[:a]], X[idx[a:b]], y[idx[a:b]],
                     X[idx[b:]], y[idx[b:]], sigma=0.0, coefficient_gate=gate)
        rec = {'split': k, 'certified': bool(r.certificate.certified),
               'tier': r.tier,
               'law': str(r.expr) if r.certificate.certified else None}
        if r.certificate.certified:
            sc = score(r.expr, fn, truth, X.shape[1], zlib.crc32(name.encode()))
            rec['wrong'] = bool(sc['wrong_form'] or sc['wrong_domain'])
        out.append(rec)
    return out


def run(name, gate):
    path = OUT / gate / f'{name}.json'
    if path.exists():
        print(f'{gate} {name}: retained', flush=True)
        return
    X, fn, truth = cases()[name]
    t0 = time.time()
    s0 = time.clock_gettime(time.CLOCK_BOOTTIME) - time.monotonic()
    r = discover_passive(X, fn(X), sigma=0.0, coefficient_gate=gate)
    rec = {'tag': 'empirical', 'case': name, 'gate': gate,
           'truth': None if truth is None else str(truth),
           'certified': bool(r.certified), 'tier': r.result.tier,
           'abstain': r.result.certificate.abstain,
           'law': str(r.result.expr) if r.certified else None,
           'rivals': [str(z) for z in r.result.rivals],
           'seconds': round(time.time() - t0, 1),
           'suspended_seconds': round(time.clock_gettime(time.CLOCK_BOOTTIME)
                                      - time.monotonic() - s0, 1)}
    if r.certified:
        rec.update(score(r.result.expr, fn, truth, X.shape[1],
                         zlib.crc32(name.encode())))
    if name.startswith(('fr3-', 'frch3-')):
        rec['splits'] = split_verdicts(X, fn, truth, gate, name)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(rec, indent=1, allow_nan=False) + '\n')
    print(f"{gate} {name}: certified={rec['certified']} "
          f"wrong={rec.get('wrong_form') or rec.get('wrong_domain')} {rec['seconds']}s",
          flush=True)


def summarize(gates=('marginal', 'joint'), banks=('fr', 'frch')):
    rows = {g: {p.stem: json.loads(p.read_text()) for p in (OUT / g).glob('*.json')}
            for g in gates}
    summary = {}
    for bank in banks:
        for g in gates:
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
        m, j = rows[gates[0]], rows[gates[1]]
        fair = [k for k in m if k.startswith(bank + '-') and k in j and not (
            m[k]['certified'] and (m[k]['wrong_form'] or m[k]['wrong_domain']))]
        ratios = [j[k]['seconds'] / max(m[k]['seconds'], 0.1) for k in fair]
        summary[f'cost/{bank}'] = {
            'fair_cases': len(fair),
            'median_ratio_joint_over_marginal': float(np.median(ratios)) if ratios else None,
            'lost_certificates': sorted(k for k in fair if m[k]['certified']
                                        and not j[k]['certified'])}
        # A4's split: same verdict (gate overhead) vs changed verdict
        same = [k for k in fair if m[k]['certified'] == j[k]['certified']
                and m[k].get('law') == j[k].get('law')]
        changed = [k for k in fair if k not in same]
        summary[f'cost/{bank}'].update({
            'same_verdict_cases': len(same),
            'same_verdict_median_ratio': float(np.median(
                [j[k]['seconds'] / max(m[k]['seconds'], 0.1) for k in same])) if same else None,
            'changed_verdict_cases': len(changed),
            'changed_verdict_median_seconds': float(np.median(
                [j[k]['seconds'] for k in changed])) if changed else None,
            'max_seconds_second_gate': max((j[k]['seconds'] for k in j
                                            if k.startswith(bank + '-')), default=None)})
    name = 'summary.json' if gates == ('marginal', 'joint') else f'summary_{gates[1]}.json'
    (OUT / name).write_text(json.dumps(summary, indent=1) + '\n')
    print(json.dumps(summary, indent=1))


def summarize_a5():
    """A5 cost classification: a case is CHANGED if the final verdicts differ
    or marginal certified a wrong law on any split; otherwise SAME."""
    m = {p.stem: json.loads(p.read_text()) for p in (OUT / 'marginal').glob('fr*3-*.json')}
    j = {p.stem: json.loads(p.read_text()) for p in (OUT / 'joint_modulo').glob('fr*3-*.json')}
    summary = {}
    for bank in ('fr3', 'frch3'):
        keys = sorted(k for k in m if k.startswith(bank + '-') and k in j)
        for g, rows in (('marginal', m), ('joint_modulo', j)):
            cert = [rows[k] for k in keys if rows[k]['certified']]
            summary[f'{g}/{bank}'] = {
                'cases': len(keys), 'certified': len(cert),
                'exact': sum(bool(r.get('exact')) for r in cert),
                'wrong': sorted(r['case'] for r in cert
                                if r['wrong_form'] or r['wrong_domain']),
                'split_wrong_cases': sorted(k for k in keys if any(
                    s_.get('wrong') for s_ in rows[k].get('splits', []))),
                'median_seconds': float(np.median([rows[k]['seconds'] for k in keys]))}
        def changed(k):
            a, b = m[k], j[k]
            return (a['certified'] != b['certified'] or a.get('law') != b.get('law')
                    or any(s_.get('wrong') for s_ in a.get('splits', [])))
        same = [k for k in keys if not changed(k)]
        chg = [k for k in keys if changed(k)]
        summary[f'cost/{bank}'] = {
            'same_cases': len(same),
            'same_median_ratio': float(np.median(
                [j[k]['seconds'] / max(m[k]['seconds'], 0.1) for k in same])) if same else None,
            'changed_cases': len(chg),
            'changed_median_seconds': float(np.median([j[k]['seconds'] for k in chg])) if chg else None,
            'max_seconds_joint_modulo': max(j[k]['seconds'] for k in keys),
            'lost_correct_certificates': sorted(
                k for k in keys if m[k]['certified'] and not (m[k]['wrong_form'] or m[k]['wrong_domain'])
                and not j[k]['certified'])}
    (OUT / 'summary_a5.json').write_text(json.dumps(summary, indent=1) + '\n')
    print(json.dumps(summary, indent=1))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--list', action='store_true')
    ap.add_argument('--case')
    ap.add_argument('--gate', choices=GATES + ('joint_modulo_a5',))
    ap.add_argument('--summarize', action='store_true')
    a = ap.parse_args()
    if a.list:
        print('\n'.join(cases()))
    elif a.summarize:
        if a.gate == 'joint_modulo_a5':
            summarize_a5()
        elif a.gate == 'joint_modulo':
            summarize(('marginal', 'joint_modulo'), ('fr2', 'frch2'))
        else:
            summarize()
    else:
        run(a.case, a.gate)


if __name__ == '__main__':
    main()

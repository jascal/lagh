"""Synthetic counter-inputs; no benchmark or tweezers data are opened.
Run with .venv/bin/python -m experiments.run_falsifiability [--baseline].
"""
import argparse
import json
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import sympy as sp

from lagh import certify as c, instrument as i


def run():
    x, z = sp.symbols('x z')
    X = np.linspace(1, 2, 40)[:, None]
    y = X[:, 0]
    tiny = np.full(40, 1e-6)
    rows = []

    def witness(name, supplied, expected, fn, predicate):
        try:
            result = fn()
            met = bool(predicate(result))
        except Exception as exc:
            result, met = {'exception': type(exc).__name__, 'message': str(exc)}, False
        rows.append(dict(gate=name, input=supplied, refusal_expected=expected,
                         observed=str(result), witness_observed=met, evidence='empirical'))

    w = witness
    w('check/wrong value', 'x=linspace(1,2,40), y=2x, law=x, eps=1e-6', '40 misses',
      lambda: c.check(x, [x], X, 2*y, tiny), lambda r: r['nmiss'] == 40)
    w('check/undefined', 'law=log(-x), x in [1,2]', '40 uncovered',
      lambda: c.check(sp.log(-x), [x], X, y, tiny), lambda r: r['nuncov'] == 40)
    for name, yy, ee in [('NaN target', np.full(40, np.nan), tiny),
                         ('NaN band', y, np.full(40, np.nan)),
                         ('infinite band', y, np.full(40, np.inf)),
                         ('negative band', y, -tiny)]:
        w('check/'+name, name+' on 40 exact x rows', 'refuse',
          lambda yy=yy, ee=ee: c.check(x, [x], X, yy, ee), lambda r: not r['certified'])
    w('check/empty', 'X=(0,1), y=eps=[]', 'refuse',
      lambda: c.check(x, [x], X[:0], y[:0], tiny[:0]), lambda r: not r['certified'])
    w('vacuous', 'y=x, eps=3', 'zero law covers: vacuous',
      lambda: c.vacuous([x], X, y, 3), bool)
    w('significance/empty', 'y=[], eps=[]', 'no evidence: log bound 0',
      lambda: c.significance_log10(x, [], [], 1), lambda r: r == 0)
    w('significance/evidence', 'y=x, eps=10', 'no evidence: log bound >=0',
      lambda: c.significance_log10(x, y, 10, 1), lambda r: r >= 0)
    w('float_pinned', 'law=2.1*x, y=2.1*x, eps=.1', 'unidentified coefficient',
      lambda: c.float_pinned(sp.Float(2.1)*x, [x], X, 2.1*y, .1)[0], lambda r: not r)
    w('pinned', 'x in [1,1.001], law=x**(3/2), sigma=.2, eps=.1; probe [1,2]',
      'rival exponent fits narrow domain but diverges on probe',
      lambda: c.pinned(x**sp.Rational(3,2), [x], 1+(X-1)*.001,
                        (1+(y-1)*.001)**1.5, .1, X, 1, .2), lambda r: not r)
    w('minimal', 'law=x+1/1000, y=x, eps=.01', 'junk term removable',
      lambda: c.minimal(x+sp.Rational(1,1000), [x], X, y, .01), lambda r: not r)
    w('parameter_interval', 'law=2*x, y=2*x, eps=100', 'unbounded at max_rel=.5',
      lambda: c.parameter_interval(2*x, [x], X, 2*y, 100, sp.Integer(2)), lambda r: r is None)
    w('admissible_interval/infeasible', 'A=[[1],[1]], y=[0,1], eps=0', 'infeasible',
      lambda: c.admissible_interval([[1],[1]], [0,1], 0), lambda r: r[0] is None)
    w('admissible_interval/unbounded', 'A=[[0]], y=[0], eps=0', 'no coefficient bound',
      lambda: c.admissible_interval([[0]], [0], 0), lambda r: r[0] == [(None,None)])
    w('admissible_interval/coefficient budget', 'A=[[1]], y=[2], eps=lambda cap:cap, cap=1, iters=1',
      'coefficient declaration unverified',
      lambda: c.admissible_interval([[1]], [2], lambda cap: cap, coeff_max=1, iters=1),
      lambda r: not r[1]['coeff_max_verified'])
    w('admissible_functional', 'A=[[1],[1]], y=[0,1], eps=0, W=[[1]]', 'infeasible',
      lambda: c.admissible_functional([[1],[1]], [0,1], 0, [[1]]), lambda r: r[0] is None)
    w('max_divergence/overlap', 'identical x laws, only 7 probes', 'insufficient overlap: infinity',
      lambda: c.max_divergence(x, x, [x], X[:7], 1), np.isinf)
    cand = lambda e: SimpleNamespace(expr=e, complexity=1)
    classes = [(x,[cand(x)]),(2*x,[cand(2*x)])]
    w('arbitrate_significance/margin', 'rivals x and 2x on 40 rows', 'no winner',
      lambda: c.arbitrate_significance(classes, y, tiny, 2), lambda r: r is None)
    many = sum(sp.Rational(k, k+41)*x**k for k in range(2,16))
    w('arbitrate_significance/rival evidence', 'x versus polynomial with 14 coefficients, 40 rows, eps=1e-6',
      'margin alone cannot defeat evidence-bearing rival',
      lambda: c.arbitrate_significance([(x,[cand(x)]),(many,[cand(many)])], y, tiny, 2), lambda r: r is None)
    w('coherent', 'laws x and 2x, probe [1,2], yscale=1', 'two classes',
      lambda: len(c.coherent([cand(x),cand(2*x)], [x], X, 1, n_evidence=40)), lambda r: r == 2)
    w('input_constraints', 'seed 0 independent uniform 40x2 design', 'no exact constraint',
      lambda: c.input_constraints(np.random.default_rng(0).uniform(1,2,(40,2)), [x,z]), lambda r: r == [])
    w('conjoin_determination/domain', 'qualified domains x>0 and x<0', 'refused',
      lambda: c.conjoin_determination([dict(qualifier=c.domain_qualifier(s)) for s in ('x>0','x<0')]),
      lambda r: r['status'] == 'refused')
    w('conjoin_determination/contradiction', 'a in [1,2] AND [3,4]', 'contradiction',
      lambda: c.conjoin_determination([c.determination([('a',lo,hi)], status='test') for lo,hi in [(1,2),(3,4)]]),
      lambda r: r.get('contradiction') == ['a'])
    w('determination/resolution', 'a in [-1,1]', 'not resolved',
      lambda: c.determination([('a',-1,1)], status='test'), lambda r: not r['n_resolved'])
    w('acf_timescale', 'zeros(1000), dt=.001', 'unreadable',
      lambda: i.acf_timescale(np.zeros(1000), .001), lambda r: r['theta'] is None)
    w('axis_gate', 'zeros(1000), dt=.001, theta_ref=1', 'refuse',
      lambda: i.axis_gate(np.zeros(1000), .001, 1), lambda r: not r['passed'])
    for name, axes, verdict in [
        ('unreadable', {'a':{'ratio':None}}, 'unreadable'),
        ('contaminated', {'a':{'ratio':.1,'var_ok':False}}, 'contaminated'),
        ('single', {'a':{'ratio':.7,'var_ok':True}}, 'unresolved'),
        ('spread', {'a':{'ratio':.5,'var_ok':True},'b':{'ratio':.7,'var_ok':True}}, 'unresolved'),
        ('missing variance', {'a':{'ratio':.7},'b':{'ratio':.7}}, 'unresolved'),
        ('missing axis', {'a':{'ratio':.7,'var_ok':True},'b':{'ratio':.7,'var_ok':True},'c':{'ratio':None}}, 'unresolved')]:
        w('gate_record/'+name, repr(axes), verdict, lambda axes=axes: i.gate_record(axes),
          lambda r, verdict=verdict: r['verdict'] == verdict)
    w('retention/fit bins', 'zeros(20), fs=100, fit_range=(10,11)', 'ValueError insufficient bins',
      lambda: retention_bins(), lambda r: r == 'ValueError')
    for ratio in (.9, 3., 1., 1.02):
        w('faxen_height/'+str(ratio), 'radius=1, ratio='+str(ratio)+(', se=.1' if ratio==1.02 else ''),
          'no finite height', lambda ratio=ratio: i.faxen_height(ratio, 1, .1 if ratio==1.02 else None),
          lambda r: r['height_um'] is None)
    for triple, verdict in [((0,1,1),'unreadable'),((.019,.885,3.797),'unattributed'),
                             ((.5,1,2),'unattributed'),((.1,1,1.04),'unattributed'),
                             ((2,3,4),'unattributed')]:
        w('attribute_deviation/'+str(triple), repr(triple), verdict,
          lambda triple=triple: i.attribute_deviation(*triple), lambda r,verdict=verdict: r['verdict']==verdict)
    return rows


def retention_bins():
    try:
        i.retention(np.zeros(20), 100, fit_range=(10,11), alpha=1, f_diode=10)
    except ValueError:
        return 'ValueError'
    return 'accepted'


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--baseline', action='store_true')
    parser.add_argument('--historical-worker', action='store_true', help=argparse.SUPPRESS)
    args = parser.parse_args()
    if args.baseline:
        import io
        import os
        import shutil
        import subprocess
        import sys
        import tarfile
        import tempfile
        revision = '582227fc59f22d6d1a5d09fbec7a575a8a2224df'
        archive = subprocess.check_output(['git', 'archive', revision, 'lagh'])
        with tempfile.TemporaryDirectory() as tmp:
            with tarfile.open(fileobj=io.BytesIO(archive)) as tar:
                tar.extractall(tmp, filter='data')
            output = Path(tmp)/'experiments/results'
            output.mkdir(parents=True)
            env = dict(os.environ, PYTHONPATH=tmp)
            subprocess.run([sys.executable, str(Path(__file__).resolve()),
                            '--historical-worker'], cwd=tmp, env=env, check=True)
            shutil.copyfile(output/'falsifiability_baseline.json',
                            'experiments/results/falsifiability_baseline.json')
        raise SystemExit(0)
    rows = run()
    suffix = '_baseline' if args.historical_worker else ''
    path = Path('experiments/results/falsifiability'+suffix+'.json')
    path.write_text(json.dumps(dict(evidence='empirical', witnesses=rows), indent=2)+'\n')
    print(f'{sum(r["witness_observed"] for r in rows)}/{len(rows)} expected refusals observed; {path}')
    for r in rows:
        if not r['witness_observed']:
            print(r['gate'], r['observed'])

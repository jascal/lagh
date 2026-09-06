"""P1 registered known-truth pilot. No scored seeds are run by this module yet."""
from __future__ import annotations

import argparse
from dataclasses import asdict
import hashlib
import json
from pathlib import Path

import numpy as np
import sympy as sp

from lagh.base import eval_expr
from lagh.certify import epsilon
from lagh.passive import discover_passive
from lagh.refusal_acquisition import RefusalPolicy, fingerprint, run_refusal_acquisition

ROOT = Path(__file__).resolve().parents[1]
x = sp.Symbol('x_0')
CASES = {
    'rational-a': ((2*x+1)/(x+3), ((.005,), (300.,))),
    'rational-b': ((3*x+2)/(x+4), ((.005,), (300.,))),
    'rational-c': ((x+1)/(x+5), ((.005,), (300.,))),
    'restricted-a': ((2*x+1)/(x+3), ((.5,), (3.,))),
}
ARMS = ('guided', 'fixed-x10', 'fixed-ladder', 'native-ladder')


def code_hash():
    h = hashlib.sha256()
    for directory in (ROOT/'lagh',):
        for file in sorted(directory.rglob('*.py')):
            h.update(str(file.relative_to(ROOT)).encode())
            h.update(file.read_bytes())
    h.update(Path(__file__).read_bytes())
    return h.hexdigest()


def score(expr, truth):
    if expr is None:
        return {'exact_truth': None, 'probe_max_error': None, 'probe_checks': 0}
    exact = sp.cancel(expr-truth) == 0
    probes = np.geomspace(.005, 300., 257)[:, None]
    predicted, target = eval_expr(expr, [x], probes), eval_expr(truth, [x], probes)
    finite = predicted is not None and np.all(np.isfinite(predicted))
    error = float(np.max(np.abs(predicted-target))) if finite else None
    return {'exact_truth': bool(exact), 'probe_max_error': error,
            'probe_checks': len(probes), 'probe_all_finite': bool(finite)}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--case', choices=tuple(CASES))
    parser.add_argument('--arm', choices=ARMS)
    args = parser.parse_args()
    directory = ROOT/'experiments/results/acquisition_p1'
    directory.mkdir(exist_ok=True)
    version = code_hash()
    seed = 10
    for name, (truth, bounds) in CASES.items():
        if args.case and args.case != name:
            continue
        X = np.random.default_rng(seed).uniform(.5, 3., (400, 1))
        y = eval_expr(truth, [x], X)
        print(f'{name}: initial discovery on400 design rows', flush=True)
        initial = discover_passive(X, y, sigma=0., seed=0)
        print(f'{name}: initial certified={initial.certified}, '
              f'rivals={len(initial.result.rivals)}', flush=True)
        initial_truth = score(initial.result.expr, truth) if initial.certified else None
        for arm in ARMS:
            if args.arm and args.arm != arm:
                continue
            path = directory/f'{name}_{arm}.json'
            if path.exists():
                data = json.loads(path.read_text())
                if data['code_hash'] != version:
                    raise RuntimeError(f'{path}: code changed; preserve old pilot and register amendment')
                print(f'{name}/{arm}: completed artifact retained', flush=True)
                continue
            planfile = directory/f'{name}_{arm}.plans.jsonl'
            if planfile.exists():
                raise RuntimeError(f'{planfile}: unfinished run; inspect before restarting')
            policy = RefusalPolicy(bounds)
            counter = 0
            def oracle(points):
                nonlocal counter
                counter += len(points)
                return eval_expr(truth, [x], points)
            with planfile.open('x') as stream:
                def record(plan):
                    stream.write(json.dumps(plan, allow_nan=False) + '\n')
                    stream.flush()
                    print(f'{name}/{arm}: {plan["role"]}, '
                          f'next={plan["count"]}, so far={plan["query_start"]}', flush=True)
                out = run_refusal_acquisition(oracle, X, y, initial_box=[[.5], [3.]],
                                              policy=policy, seed=seed, strategy=arm,
                                              initial_result=initial, on_plan=record)
            if out.queries != len(X) + counter:
                raise RuntimeError('actual oracle cost does not match reported cost')
            truth_score = score(out.result.expr, truth) if out.certified else score(None, truth)
            data = {'tag': 'empirical', 'stage': 'pilot P1', 'case': name,
                    'arm': arm, 'seed': seed, 'code_hash': version,
                    'policy': asdict(policy), 'initial_hash': fingerprint(X, y),
                    'initial_structural': out.initial_structural,
                    'initial_rivals': out.initial_rivals, 'initial_certified': initial.certified,
                    'initial_truth': initial_truth,
                    'certified': out.certified, 'queries': out.queries,
                    'actual_new_queries': counter,
                    'law': str(out.result.expr) if out.certified else None,
                    'abstain': out.result.certificate.abstain,
                    'notes': out.result.certificate.notes,
                    'native_verdict': out.native_verdict, 'final_passed': out.final_passed,
                    'domain': out.domain, 'truth_score': truth_score,
                    'history': out.history,
                    'wrong_certified': out.certified and not truth_score['exact_truth'],
                    'initial_representatives': [str(e) for e in initial.result.rivals],
                    'initial_rival_misses': []}
            # Read-only diagnostic, never fed to selection: split representatives
            # need not survive every initial row (both old twins can be wrong).
            from lagh.certify import check
            for e in initial.result.rivals:
                c = check(e, [x], X, y, epsilon(y))
                data['initial_rival_misses'].append({'nmiss': c['nmiss'], 'nuncov': c['nuncov']})
            path.write_text(json.dumps(data, indent=2, allow_nan=False) + '\n')
            print(f'{name}/{arm}: certified={out.certified}, queries={out.queries}, '
                  f'wrong={data["wrong_certified"]}', flush=True)


if __name__ == '__main__':
    main()

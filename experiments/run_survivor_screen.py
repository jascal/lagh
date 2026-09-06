"""Registered P2 design-only screen; all16 seeds retained, no future samples."""
import argparse
import json
from pathlib import Path

import numpy as np
import sympy as sp

from experiments.run_refusal_comparison import code_hash, score
from lagh.certify import check, epsilon
from lagh.passive import discover_passive
from lagh.refusal_acquisition import fingerprint


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--start', type=int, default=20)
    parser.add_argument('--stop', type=int, default=36)
    args = parser.parse_args()
    if not 20 <= args.start < args.stop <= 36:
        raise ValueError('only registered P2 seeds20 through35 may run')
    directory = Path('experiments/results/acquisition_p2')
    directory.mkdir(exist_ok=True)
    version = code_hash()
    x = sp.Symbol('x_0')
    truth = (3*x+2)/(x+4)
    for seed in range(args.start, args.stop):
        path = directory/f'seed_{seed}.json'
        if path.exists():
            if json.loads(path.read_text())['code_hash'] != version:
                raise RuntimeError('code changed; register an amendment instead of replacing results')
            print(f'seed{seed}: completed artifact retained', flush=True)
            continue
        X = np.random.default_rng(seed).uniform(.5, 3., (400, 1))
        y = (3*X[:, 0]+2)/(X[:, 0]+4)
        print(f'seed{seed}:400 design observations', flush=True)
        out = discover_passive(X, y, sigma=0., seed=0)
        records = []
        for expr in out.result.rivals:
            witness = check(expr, [x], X, y, epsilon(y))
            records.append({'expr': str(expr), **dict(witness)})
        n = sum(r['certified'] for r in records)
        eligible = out.result.certificate.abstain == 'structural' and n >= 2
        data = {'tag': 'empirical', 'stage': 'P2 design only', 'seed': seed,
                'code_hash': version, 'design_queries':400, 'final_queries':0,
                'data_hash': fingerprint(X, y), 'certified': out.certified,
                'law': str(out.result.expr) if out.certified else None,
                'truth_score': score(out.result.expr, truth) if out.certified else None,
                'abstain': out.result.certificate.abstain, 'rivals': records,
                'full_data_survivors': n, 'eligible': bool(eligible)}
        path.write_text(json.dumps(data, indent=2, allow_nan=False)+'\n')
        print(f'seed{seed}: certified={out.certified}, rivals={len(records)}, '
              f'survivors={n}, eligible={eligible}', flush=True)


if __name__ == '__main__':
    main()

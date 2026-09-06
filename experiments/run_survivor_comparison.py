"""P3 and separately registered scored trials on the frozen P2 twin bank."""
from __future__ import annotations

import argparse
from dataclasses import asdict
import hashlib
import json
from pathlib import Path

import numpy as np
import sympy as sp

from experiments.run_refusal_comparison import ARMS, ROOT, code_hash, score
from lagh.passive import discover_passive
from lagh.refusal_acquisition import RefusalPolicy, fingerprint, run_refusal_acquisition

DESIGN_SEEDS = (30, 32, 35)
DOMAINS = {'wide': ((.005,), (300.,)), 'restricted': ((.5,), (3.,))}


def study_hash():
    return hashlib.sha256(code_hash().encode()+Path(__file__).read_bytes()).hexdigest()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--stage', choices=('pilot', 'scored'), default='pilot')
    parser.add_argument('--design-seed', type=int, choices=DESIGN_SEEDS)
    parser.add_argument('--domain', choices=tuple(DOMAINS))
    parser.add_argument('--arm', choices=ARMS)
    parser.add_argument('--followup-seed', type=int)
    args = parser.parse_args()
    version = study_hash()
    if args.stage == 'pilot':
        followup_seeds = (99,)
    else:
        path = ROOT/'experiments/acquisition_scored_protocol.json'
        if not path.exists():
            raise RuntimeError('scored protocol has not been frozen; no scored samples generated')
        protocol = json.loads(path.read_text())
        if protocol['code_hash'] != version:
            raise RuntimeError('scored code differs from the frozen protocol')
        followup_seeds = tuple(protocol['followup_seeds'])
    if args.followup_seed is not None:
        if args.followup_seed not in followup_seeds:
            raise ValueError('follow-up seed is not registered for this stage')
        followup_seeds = (args.followup_seed,)
    directory = ROOT/'experiments/results'/('acquisition_p3' if args.stage == 'pilot'
                                          else 'acquisition_scored')
    directory.mkdir(exist_ok=True)
    x = sp.Symbol('x_0')
    truth = (3*x+2)/(x+4)
    for design_seed in DESIGN_SEEDS:
        if args.design_seed is not None and args.design_seed != design_seed:
            continue
        X = np.random.default_rng(design_seed).uniform(.5, 3., (400, 1))
        y = (3*X[:, 0]+2)/(X[:, 0]+4)
        recorded = json.loads((ROOT/f'experiments/results/acquisition_p2/seed_{design_seed}.json').read_text())
        if fingerprint(X, y) != recorded['data_hash'] or not recorded['eligible']:
            raise RuntimeError('initial twin bank does not match the frozen design artifact')
        print(f'design{design_seed}: reconstructing frozen refusal', flush=True)
        initial = discover_passive(X, y, sigma=0., seed=0)
        if [str(e) for e in initial.result.rivals] != [r['expr'] for r in recorded['rivals']]:
            raise RuntimeError('initial rivals changed; do not silently rerun a different experiment')
        for domain, bounds in DOMAINS.items():
            if args.domain and args.domain != domain:
                continue
            for seed in followup_seeds:
                for arm in ARMS:
                    if args.arm and args.arm != arm:
                        continue
                    stem = f'design{design_seed}_{domain}_seed{seed}_{arm}'
                    path = directory/f'{stem}.json'
                    if path.exists():
                        if json.loads(path.read_text())['code_hash'] != version:
                            raise RuntimeError(f'{stem}: preserve old artifact; code changed')
                        print(f'{stem}: completed artifact retained', flush=True)
                        continue
                    plans = directory/f'{stem}.plans.jsonl'
                    if plans.exists():
                        raise RuntimeError(f'{plans}: unfinished plan ledger, inspect before restarting')
                    policy = RefusalPolicy(bounds)
                    counter = 0
                    def oracle(points):
                        nonlocal counter
                        counter += len(points)
                        return (3*points[:, 0]+2)/(points[:, 0]+4)
                    with plans.open('x') as stream:
                        def record(plan):
                            stream.write(json.dumps(plan, allow_nan=False)+'\n')
                            stream.flush()
                            print(f'{stem}: {plan["role"]}, next{plan["count"]}, '
                                  f'spent{plan["query_start"]}', flush=True)
                        out = run_refusal_acquisition(oracle, X, y, initial_box=[[.5], [3.]],
                                                      policy=policy, seed=seed, strategy=arm,
                                                      initial_result=initial, on_plan=record)
                    if out.queries != len(X)+counter:
                        raise RuntimeError('ledger cost differs from actual oracle observations')
                    truth_score = score(out.result.expr, truth) if out.certified else score(None, truth)
                    data = {'tag': 'empirical', 'stage': args.stage,
                            'design_seed': design_seed, 'domain_name': domain,
                            'followup_seed': seed, 'arm': arm, 'code_hash': version,
                            'policy': asdict(policy), 'initial_hash': fingerprint(X,y),
                            'initial_structural': out.initial_structural,
                            'initial_split_rivals': out.initial_rivals,
                            'initial_full_survivors': recorded['full_data_survivors'],
                            'certified': out.certified, 'queries': out.queries,
                            'actual_new_queries': counter,
                            'law': str(out.result.expr) if out.certified else None,
                            'abstain': out.result.certificate.abstain,
                            'notes': out.result.certificate.notes,
                            'native_verdict': out.native_verdict,
                            'native_law': out.native_law,
                            'native_truth': (score(sp.sympify(out.native_law), truth)
                                             if out.native_law else None),
                            'final_passed': out.final_passed, 'domain': out.domain,
                            'truth_score': truth_score, 'history': out.history,
                            'wrong_certified': out.certified and not truth_score['exact_truth']}
                    path.write_text(json.dumps(data, indent=2, allow_nan=False)+'\n')
                    print(f'{stem}: certified={out.certified}, queries={out.queries}, '
                          f'wrong={data["wrong_certified"]}', flush=True)


if __name__ == '__main__':
    main()

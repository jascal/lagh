"""Registered S1 completeness, cost and soundness scoring; no fitting or tuning."""
from __future__ import annotations

from itertools import product
import json
from pathlib import Path

import numpy as np
import sympy as sp

from experiments.run_refusal_comparison import ROOT, score
from lagh.certify import check, epsilon
from lagh.refusal_acquisition import fingerprint

TRUTH = (3*sp.Symbol('x_0')+2)/(sp.Symbol('x_0')+4)


def key(row):
    return (row['design_seed'],row['domain_name'],row['followup_seed'],row['arm'])


def summarize(rows, protocol):
    expected=set(product(protocol['design_seeds'],protocol['domains'],
                         protocol['followup_seeds'],protocol['arms']))
    keys=[key(row) for row in rows]
    if len(keys)!=len(set(keys)) or set(keys)!=expected:
        raise ValueError('missing, duplicated, or unexpected trial/arm; summary is incomplete')
    if any(r['code_hash']!=protocol['code_hash'] for r in rows):
        raise ValueError('trial code differs from registered code')
    truth={key(r): (score(sp.sympify(r['law']),TRUTH)['exact_truth']
                    if r['certified'] else False) for r in rows}
    wrong=[key(r) for r in rows if r['certified'] and not truth[key(r)]]
    bykey={key(r):r for r in rows}
    arms={}
    for arm in protocol['arms']:
        group=[r for r in rows if r['arm']==arm]
        resolved=[r for r in group if r['certified'] and truth[key(r)]]
        costs=[r['queries'] if r['certified'] and truth[key(r)] else
               protocol['unresolved_penalty'] for r in group]
        arms[arm]={'trials':len(group),'resolved':len(resolved),
                   'unresolved':len(group)-len(resolved),
                   'resolved_queries':sorted(r['queries'] for r in resolved),
                   'actual_queries':sum(r['queries'] for r in group),
                   'actual_new_queries':sum(r['actual_new_queries'] for r in group),
                   'mean_censored_cost':float(np.mean(costs)),
                   'native_raw_certificates':sum(r.get('native_verdict') is True for r in group)}
    pairs={}
    for baseline in ('fixed-x10','fixed-ladder','native-ladder'):
        if baseline not in protocol['arms']:
            continue
        diffs=[]
        for design,domain,seed in product(protocol['design_seeds'],protocol['domains'],
                                          protocol['followup_seeds']):
            g=bykey[(design,domain,seed,'guided')]
            b=bykey[(design,domain,seed,baseline)]
            if truth[key(g)] and truth[key(b)]:
                diffs.append(g['queries']-b['queries'])
        pairs[baseline]={'joint_resolutions':len(diffs),
                         'differences':diffs,
                         'median':float(np.median(diffs)) if diffs else None}
    gain=not wrong and all(
        arms['guided']['resolved']>=arms[b]['resolved'] and
        pairs[b]['median'] is not None and pairs[b]['median']<0
        for b in ('fixed-x10','fixed-ladder'))
    return {'tag':'empirical','complete':True,'total_arm_runs':len(rows),
            'arms':arms,'paired_query_differences':pairs,
            'wrong_terminal_certificates':len(wrong),'wrong_trials':wrong,
            'strict_gain_over_both_fixed_strategies':bool(gain)}


def audit_trial(row, plans_path):
    """Reconstruct declared inputs, finite domain and checks from persisted plans."""
    d=row['design_seed']; seed=row['followup_seed']
    initial=np.random.default_rng(d).uniform(.5,3.,(400,1))
    target=lambda X:(3*X[:,0]+2)/(X[:,0]+4)
    if fingerprint(initial,target(initial))!=row['initial_hash']:
        raise ValueError('initial data hash differs')
    plans=[json.loads(line) for line in Path(plans_path).read_text().splitlines() if line]
    observed=[r for r in row['history'] if 'count' in r]
    if len(plans)!=len(observed):
        raise ValueError('plan/observation count differs')
    spent=400; designs=[initial]; finals=[]
    for plan,record in zip(plans,observed):
        if {k:v for k,v in record.items() if k!='y_hash'}!=plan:
            raise ValueError('pre-query plan differs from recorded observation')
        X=np.asarray(plan['X'],float)
        if (X.shape!=(plan['count'],1) or fingerprint(X)!=plan['X_hash']
                or plan['query_start']!=spent or not np.all(np.isfinite(X))):
            raise ValueError('query points, hash, shape or cost inconsistent')
        bounds=np.asarray(row['policy']['admissible'],float)
        # exp(log(endpoint)) can differ by one ULP from the declared endpoint.
        if np.any(X<np.nextafter(bounds[0],-np.inf)) or np.any(X>np.nextafter(bounds[1],np.inf)):
            raise ValueError('measurement outside admissible domain')
        if fingerprint(target(X))!=record['y_hash']:
            raise ValueError('oracle values inconsistent with the owned substrate')
        spent+=len(X)
        if plan['role']=='final':
            if finals or record is not observed[-1]:
                raise ValueError('final observations reused for subsequent design')
            finals.append(X)
        else:
            designs.append(X)
    if spent!=row['queries'] or spent!=400+row['actual_new_queries']:
        raise ValueError('not every oracle observation is charged')
    design=np.vstack(designs)
    if finals:
        box=np.array([design.min(axis=0),design.max(axis=0)])
        final_seed=np.random.SeedSequence(seed).spawn(3)[1]
        expected=np.exp(np.random.default_rng(final_seed).uniform(
            np.log(box[0]),np.log(box[1]),(row['policy']['final_points'],1)))
        if not np.array_equal(expected,finals[0]):
            raise ValueError('final sample differs from the independently declared RNG stream')
    if row['certified']:
        if not finals or row['final_passed'] is not True or row['domain'] is None:
            raise ValueError('certificate lacks independent final evidence')
        if plans[-1].get('frozen_law') != row['law']:
            raise ValueError('terminal law differs from the law frozen before final queries')
        all_X=np.vstack([design,*finals]); all_y=target(all_X)
        expected_domain={'kind':'finite observed rows, not the continuous bounding box',
                         'bounds':[all_X.min(axis=0).tolist(),all_X.max(axis=0).tolist()],
                         'rows':len(all_X),'design_hash':fingerprint(design,target(design)),
                         'final_hash':fingerprint(finals[0],target(finals[0]))}
        if row['domain']!=expected_domain:
            raise ValueError('certificate inherited or misstated its finite domain')
        law=sp.sympify(row['law'])
        if not check(law,[sp.Symbol('x_0')],all_X,all_y,epsilon(all_y))['certified']:
            raise ValueError('candidate fails reconstructed measured observations')
        bank=json.loads((ROOT/f'experiments/results/acquisition_p2/seed_{d}.json').read_text())
        for rival in bank['rivals']:
            if rival['certified']:
                expr=sp.sympify(rival['expr'])
                if expr!=law and check(expr,[sp.Symbol('x_0')],design,target(design),
                                       epsilon(target(design)))['certified']:
                    raise ValueError('old compatible rival was not eliminated by design measurements')
    return {'plans':len(plans),'final_samples':len(finals),'charged_queries':spent}


def main():
    protocol=json.loads((ROOT/'experiments/acquisition_scored_protocol.json').read_text())
    folder=ROOT/'experiments/results/acquisition_scored'
    paths=sorted(folder.glob('*.json'))
    rows=[json.loads(p.read_text()) for p in paths]
    summary=summarize(rows,protocol)  # completeness before claiming any totals
    audits=[audit_trial(r,p.with_suffix('.plans.jsonl')) for r,p in zip(rows,paths)]
    summary['audited_trials']=len(audits)
    summary['prequery_plan_records']=sum(a['plans'] for a in audits)
    summary['independent_final_samples']=sum(a['final_samples'] for a in audits)
    predictions=protocol['predictions']
    summary['predictions']={
        'resolution_counts':all(summary['arms'][a]['resolved']==n for a,n in
                                predictions['resolutions'].items()),
        'minimum_unresolved':summary['arms']['guided']['unresolved']>=predictions['guided_unresolved_minimum'],
        'paired_medians':all(summary['paired_query_differences'][a]['median']==n for a,n in
                             predictions['paired_median_differences'].items()),
        'zero_wrong':summary['wrong_terminal_certificates']==0,
        'expanded_wide_domains':all(r['domain']['bounds'][1][0]>3 for r in rows
                                    if r['certified'] and r['domain_name']=='wide'),
        'successful_query_counts':all(all(q==n for q in summary['arms'][a]['resolved_queries'])
                                      for a,n in predictions['successful_queries'].items())}
    out=ROOT/'experiments/results/acquisition_scored_summary.json'
    out.write_text(json.dumps(summary,indent=2)+'\n')
    print(json.dumps(summary,indent=2))


if __name__=='__main__':
    main()

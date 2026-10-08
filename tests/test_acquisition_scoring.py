"""Registered scoring falsifiers, using synthetic current-runtime observations."""
from copy import deepcopy
import json

import numpy as np
import pytest

from experiments.score_refusal_acquisition import TRUTH, audit_trial, summarize
from lagh.refusal_acquisition import fingerprint


def fixture():
    protocol={'design_seeds':[30],'domains':['wide'],'followup_seeds':[100],
              'arms':['guided','fixed-x10','fixed-ladder'],'code_hash':'test',
              'unresolved_penalty':1801}
    rows=[{'design_seed':30,'domain_name':'wide','followup_seed':100,'arm':arm,
           'code_hash':'test','certified':True,'law':str(TRUTH),'queries':cost,
           'actual_new_queries':cost-400} for arm,cost in
          [('guided',560),('fixed-x10',560),('fixed-ladder',640)]]
    return rows,protocol


def test_tie_is_not_a_strict_gain():
    rows,p=fixture()
    assert summarize(rows,p)['strict_gain_over_both_fixed_strategies'] is False
    # Positive counterpart prevents an always-false gain test.
    rows[0]['queries']=540
    rows[0]['actual_new_queries']=140
    assert summarize(rows,p)['strict_gain_over_both_fixed_strategies'] is True


def test_wrong_law_is_scored_even_if_the_verdict_says_certified():
    rows,p=fixture()
    rows[0]['law']='x_0 + 100'
    out=summarize(rows,p)
    assert out['wrong_terminal_certificates']==1
    assert out['strict_gain_over_both_fixed_strategies'] is False


@pytest.mark.parametrize('damage',['missing','duplicate'])
def test_incomplete_or_duplicate_arms_invalidate_summary(damage):
    rows,p=fixture()
    rows=rows[:-1] if damage=='missing' else rows+[deepcopy(rows[0])]
    with pytest.raises(ValueError,match='missing, duplicated'):
        summarize(rows,p)


def test_cheap_abstention_does_not_win():
    rows,p=fixture()
    rows[0].update(certified=False,law=None,queries=400,actual_new_queries=0)
    out=summarize(rows,p)
    assert out['arms']['guided']['mean_censored_cost']==1801
    assert out['strict_gain_over_both_fixed_strategies'] is False


def synthetic_trial(tmp_path):
    """Independent plans and observations without archived exp/log byte assumptions.

    Exact stream/hash checks still run in audit_trial. This fixture is never
    written to study artifacts or counted as empirical discovery evidence.
    """
    target=lambda X:(3*X[:,0]+2)/(X[:,0]+4)
    initial=np.random.default_rng(30).uniform(.5,3.,(400,1))
    added=np.geomspace(.05,30.,80)[:,None]
    design=np.vstack([initial,added])
    box=np.array([design.min(axis=0),design.max(axis=0)])
    final_seed=np.random.SeedSequence(99).spawn(3)[1]
    final=np.exp(np.random.default_rng(final_seed).uniform(
        np.log(box[0]),np.log(box[1]),(80,1)))
    plans=[]
    history=[]
    for role,X,start in [('design',added,400),('final',final,480)]:
        plan={'role':role,'X':X.tolist(),'X_hash':fingerprint(X),
              'count':len(X),'query_start':start}
        if role=='final':
            plan['frozen_law']=str(TRUTH)
        plans.append(plan)
        history.append({**deepcopy(plan),'y_hash':fingerprint(target(X))})
    all_X=np.vstack([design,final])
    row={'design_seed':30,'followup_seed':99,'initial_hash':fingerprint(initial,target(initial)),
         'policy':{'admissible':[[.05],[30.]],'final_points':80},
         'history':history,'queries':560,'actual_new_queries':160,
         'certified':True,'final_passed':True,'law':str(TRUTH),
         'domain':{'kind':'finite observed rows, not the continuous bounding box',
                   'bounds':[all_X.min(axis=0).tolist(),all_X.max(axis=0).tolist()],
                   'rows':len(all_X),'design_hash':fingerprint(design,target(design)),
                   'final_hash':fingerprint(final,target(final))}}
    path=tmp_path/'plans.jsonl'
    path.write_text('\n'.join(json.dumps(p) for p in plans)+'\n')
    assert audit_trial(row,path)['charged_queries']==560
    return row,path,plans


def test_old_bounds_fail_domain_audit(tmp_path):
    row,plans,_=synthetic_trial(tmp_path)
    assert audit_trial(row,plans)['charged_queries']==560
    row['domain']['bounds']=[[.5],[3.]]
    with pytest.raises(ValueError,match='inherited or misstated'):
        audit_trial(row,plans)


def test_final_observations_cannot_license_a_different_law(tmp_path):
    # A plan for x+100 cannot be credited as final evidence for the true law,
    # even if both the saved plan and observation metadata agree with each other.
    row,_,plans=synthetic_trial(tmp_path)
    plans[-1]['frozen_law']='x_0 + 100'
    next(r for r in row['history'] if r['role']=='final')['frozen_law']='x_0 + 100'
    fake=tmp_path/'plans.jsonl'
    fake.write_text('\n'.join(json.dumps(p) for p in plans)+'\n')
    with pytest.raises(ValueError,match='law frozen before final'):
        audit_trial(row,fake)


def test_final_stream_check_rejects_even_one_ulp_with_matching_hashes(tmp_path):
    row,path,plans=synthetic_trial(tmp_path)
    final=np.asarray(plans[-1]['X'])
    final[0,0]=np.nextafter(final[0,0],np.inf)
    plans[-1].update(X=final.tolist(),X_hash=fingerprint(final))
    # Both plan and observation agree, but the independently seeded stream does not.
    row['history'][-1]={**deepcopy(plans[-1]),
        'y_hash':fingerprint((3*final[:,0]+2)/(final[:,0]+4))}
    path.write_text('\n'.join(json.dumps(p) for p in plans)+'\n')
    with pytest.raises(ValueError,match='independently declared RNG stream'):
        audit_trial(row,path)

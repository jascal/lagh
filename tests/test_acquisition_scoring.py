"""Registered scoring falsifiers, using synthetic results and frozen P3 data."""
from copy import deepcopy
import json

import pytest

from experiments.score_refusal_acquisition import ROOT, TRUTH, audit_trial, summarize


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


def test_old_bounds_fail_domain_audit():
    path=ROOT/'experiments/results/acquisition_p3/design30_wide_seed99_guided.json'
    row=json.loads(path.read_text())
    plans=path.with_suffix('.plans.jsonl')
    assert audit_trial(row,plans)['charged_queries']==560
    row['domain']['bounds']=[[.5],[3.]]
    with pytest.raises(ValueError,match='inherited or misstated'):
        audit_trial(row,plans)


def test_final_observations_cannot_license_a_different_law(tmp_path):
    # A plan for x+100 cannot be credited as final evidence for the true law,
    # even if both the saved plan and observation metadata agree with each other.
    path=ROOT/'experiments/results/acquisition_p3/design30_wide_seed99_guided.json'
    row=json.loads(path.read_text())
    plans=[json.loads(s) for s in path.with_suffix('.plans.jsonl').read_text().splitlines()]
    plans[-1]['frozen_law']='x_0 + 100'
    next(r for r in row['history'] if r['role']=='final')['frozen_law']='x_0 + 100'
    fake=tmp_path/'plans.jsonl'
    fake.write_text('\n'.join(json.dumps(p) for p in plans)+'\n')
    with pytest.raises(ValueError,match='law frozen before final'):
        audit_trial(row,fake)

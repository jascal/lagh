"""Counter-inputs recorded in falsifiability_baseline.json before these fixes."""
import numpy as np
import sympy as sp
import pytest

from lagh import certify as c
from lagh.instrument import gate_record, faxen_height, covariance_timescale


@pytest.mark.parametrize('y,eps', [(np.nan, .1), (1., np.nan), (1., np.inf), (1., -.1)])
def test_invalid_rows_cannot_certify(y, eps):
    # Witness: each invalid value must defeat an otherwise exact x=1 claim.
    r = c.check(sp.Symbol('x'), [sp.Symbol('x')], np.ones((1,1)), [y], eps)
    assert not r['certified'] and r['nuncov'] == 1


def test_empty_domain_has_no_evidence():
    assert not c.check(sp.S.One, [], np.empty((0,0)), [], [])['certified']
    assert c.significance_log10(sp.S.One, [], [], 1) == 0


def test_evidence_comes_from_one_evaluation(monkeypatch):
    # Witness: second evaluation deliberately returns a different, perfect law.
    calls=[]
    def evaluated(*args):
        calls.append(1)
        return np.array([1.,1.]) if len(calls)==1 else np.array([1.,3.])
    monkeypatch.setattr(c, 'eval_expr', evaluated)
    r=c.check(sp.S.One, [], np.empty((2,0)), [1.,3.], .1, row_indices=[8,2])
    m=r.measurement(domain='held-out fixture')['measurement']
    assert len(calls)==1 and not r['certified']
    assert r.residual.tolist()==[0.,2.] and r.row_indices.tolist()==[8,2]
    assert m['row_indices'][0]==2 and m['residual'][0]==2.
    assert m['n_exceeding']==1


def test_common_mode_needs_all_axes_and_variance():
    missing={'a':{'ratio':.7}, 'b':{'ratio':.7}}
    assert gate_record(missing)['verdict']=='unresolved'
    good={k:dict(v,var_ok=True) for k,v in missing.items()}
    assert gate_record(good)['verdict']=='common-mode'
    assert gate_record(dict(good,c={'ratio':None}))['verdict']=='unresolved'


def test_bulk_height_and_unreadable_timescale():
    assert faxen_height(1.,1.)['height_um'] is None
    assert covariance_timescale(np.zeros(10000), .001)['theta'] is None
    assert covariance_timescale(np.arange(100.), .001)['theta'] is None


def test_evidence_survives_caller_mutation():
    # Witness exhibited in residual_mutation_baseline.json before snapshot fix.
    y=np.array([1.,3.]); eps=np.array([.1,.1]); ids=np.array([8,2])
    r=c.check(sp.S.One,[],np.empty((2,0)),y,eps,row_indices=ids)
    y[:]=1; eps[:]=100; ids[:]=[0,1]
    m=r.measurement()['measurement']
    assert m['n_exceeding']==1 and m['row_indices'][0]==2
    assert m['residual'][0]==2 and m['epsilon'][0]==.1

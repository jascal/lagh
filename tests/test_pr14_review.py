"""Counter-input families registered and exhibited before crediting these fixes."""
from pathlib import Path
from types import SimpleNamespace
import importlib.util

import numpy as np
import pytest
import sympy as sp
from lagh import certify as c, engine
from lagh.instrument import faxen_height, gate_record, covariance_timescale
from lagh.mcp import core

_spec=importlib.util.spec_from_file_location('compare',Path(__file__).resolve().parents[1]/'experiments/compare_artifacts.py')
_compare=importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(_compare)


@pytest.mark.parametrize('X',[np.linspace(1,2,5),np.array(1.),np.ones((5,1,1)),np.ones((5,2)),np.full((5,1),'bad')])
def test_invalid_domain_shapes_refuse(X):
    x=sp.Symbol('x_0')
    r=c.check(x,[x],X,2*np.linspace(1,2,5),1e-6)
    assert not r['certified'] and r['nuncov']==5


@pytest.mark.parametrize('indices',[[0,1],[-1,0,1,2,3],[0,0,1,2,3],[0.,1.,2.,3.,4.]])
def test_indices_fail_at_check_boundary(indices):
    with pytest.raises(ValueError,match='row_indices'):
        c.check(sp.S.One,[],np.empty((5,0)),np.ones(5),.1,row_indices=indices)


@pytest.mark.parametrize('F',[np.nextafter(1.,2.),1.+1e-12,1.+1e-9,1.5])
def test_faxen_near_bulk_requires_precision_and_stable_inversion(F):
    unknown=faxen_height(F,1.)
    assert unknown['height_um'] is None and not unknown['height_determined']
    assert np.isfinite(unknown['conditional_height_um'])
    exact=faxen_height(F,1.,0.)
    assert exact['height_determined'] and exact['height_um']==unknown['conditional_height_um']
    # Near bulk h*(F-1) tends to9/16, an independent asymptotic check of
    # the scaled root (the old 1e-9 lower bracket fails the entire tiny family).
    if F-1<1e-6:
        assert exact['height_um']*(F-1)==pytest.approx(9/16,rel=1e-6)
    bulk=faxen_height(F,1.,F-1+.01)
    assert bulk['height_um'] is None and not bulk['height_determined']


@pytest.mark.parametrize('se',[float('nan'),float('inf'),-.1])
def test_invalid_uncertainty_cannot_mean_exact(se):
    r=faxen_height(1.5,1.,se)
    assert r['height_um'] is None and r['conditional_height_um'] is None
    assert 'uncertainty' in r['note']


def test_missing_variance_never_grants_clean():
    assert gate_record({'a':{'ratio':1.},'b':{'ratio':1.}})['verdict']=='unresolved'
    assert gate_record({'a':{'ratio':1.,'var_ok':True}})['verdict']=='clean'


def test_gls_refusal_retains_schema():
    r=covariance_timescale(np.zeros(100),.001)
    assert r['theta'] is None and r['evidence']=='empirical' and r['n_blocks']==0
    assert r['n_lags']==192 and r['amplitude_normalized'] is None


def test_quantized_boundary_has_one_ulp_counter_input():
    v=1.234567890125
    neighbor=np.nextafter(v,np.inf if format(v,'.12g')=='1.23456789012' else -np.inf)
    a,b=(float(format(z,'.12g')) for z in (v,neighbor))
    assert a!=b and _compare.differences(a,b)
    assert not _compare.differences(a,b,significant_digits=12)
    assert _compare.differences(a,a*1.00001,significant_digits=12)


def test_verify_does_not_blame_form_for_nan_band():
    X=np.linspace(1,2,30)
    r=core.verify(X,2*X,'x_0',floor_abs=float('nan'))
    assert not r['certified'] and 'band' in r['note'] and 'diverges' not in r['note']


def test_recover_surfaces_both_refusal_routes(monkeypatch):
    x=sp.Symbol('x_0');X=np.linspace(1,2,30)[:,None];y=2*X[:,0]
    cert=c.Certificate(False,30,0,30,[(1,2)],str(x),abstain='structural')
    c.attach_check_evidence(cert,c.check(x,[x],X,y,.01),domain='fixture')
    result=engine.Result(cert,None,1,1)
    monkeypatch.setattr(core,'discover_passive',lambda *a,**kw:SimpleNamespace(result=result))
    passive=core.recover(X,y,max_tier=1)
    active=SimpleNamespace(result=result,box_final=[[1],[2]],queries_used=30,
                           ledger=SimpleNamespace(spent=30),ranging_trajectory=[])
    monkeypatch.setattr(core,'run_active',lambda *a,**kw:active)
    monkeypatch.setattr(core,'_characterize_oracle',lambda *a,**kw:None)
    acquired=core.recover(oracle=lambda X:2*X[:,0],box=[[1],[2]])
    for r in (passive,acquired):
        assert r['measurement']['n_exceeding']==30 and r['tag']=='open'
        assert r['measurement']['candidate']=='x_0'


@pytest.mark.parametrize('cause',['parametric','structural','C6'])
def test_engine_retains_existing_evidence_for_other_abstains(monkeypatch,cause):
    x=sp.Symbol('x_0')
    X=np.arange(1.,101.)[:,None] if cause=='C6' else np.linspace(1,2,100)[:,None]
    y=2*X[:,0] if cause=='C6' else X[:,0]
    candidates=[] if cause=='C6' else [SimpleNamespace(expr=x,complexity=1,channel='witness')]
    monkeypatch.setattr(engine,'_tier_candidates',lambda *a,**kw:candidates)
    if cause=='parametric':
        monkeypatch.setattr(engine,'pinned',lambda *a,**kw:False)
    elif cause=='structural':
        monkeypatch.setattr(engine,'coherent',lambda *a,**kw:[(x,candidates),(2*x,candidates)])
        monkeypatch.setattr(engine,'input_constraints',lambda *a,**kw:[])
        monkeypatch.setattr(engine,'arbitrate_significance',lambda *a,**kw:None)
    else:
        monkeypatch.setattr(engine.c6_quasipoly,'is_integer_lattice',lambda *a:True)
        monkeypatch.setattr(engine.c6_quasipoly,'recover_integer',lambda *a:SimpleNamespace(
            certified=True,quasipoly=x,domain_size=100,note='deliberately wrong integer candidate'))
    r=engine.discover(X[:60],y[:60],X[60:80],y[60:80],X[80:],y[80:],max_tier=6 if cause=='C6' else 1)
    assert not r.certificate.certified
    assert r.certificate.measurement['candidate']=='x_0'
    if cause=='C6':
        assert r.certificate.measurement['domain']=='all engine rows'
        assert r.certificate.measurement['n_exceeding']==100
    else:
        assert r.certificate.abstain==cause
        assert r.certificate.measurement['n_exceeding']==0
        assert 'not necessarily' in r.certificate.measurement['role']

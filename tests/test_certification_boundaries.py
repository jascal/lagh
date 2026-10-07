"""Measured BND1 failures: path parity and the meaning of parameter slices."""
import numpy as np
import sympy as sp
import pytest

from experiments.boundary_witnesses import interval_fixture, prepass_result
from lagh.certify import check, epsilon, parameter_interval


def test_prepass_honors_the_selected_joint_coefficient_gate():
    assert prepass_result("marginal").certificate.certified
    assert not prepass_result("joint_quotient").certificate.certified


def test_parameter_interval_endpoints_are_checked_even_below_initial_step():
    x = sp.Symbol("x_0")
    X = np.linspace(.5, 3., 80)[:, None]
    a = sp.Rational(13, 10)
    y = 1.3 * X[:, 0]
    eps = epsilon(y)
    iv = parameter_interval(a*x, [x], X, y, eps, a)
    assert iv is not None
    assert iv[0] <= 1.3 <= iv[1]
    for endpoint in iv:
        assert check(sp.Float(endpoint)*x, [x], X, y, eps)["certified"]


def test_coordinate_slices_are_not_joint_boxes_or_marginal_coverage():
    X, y, eps, expr, (a, b) = interval_fixture()
    syms = list(sp.symbols("x_0:2"))
    ia = parameter_interval(expr, syms, X, y, eps, a)
    ib = parameter_interval(expr, syms, X, y, eps, b)
    for atom, iv in ((a, ia), (b, ib)):
        assert check(expr.xreplace({atom: sp.Float(iv[1])}), syms, X, y, eps)["certified"]
    corner = expr.xreplace({a: sp.Float(ia[1]), b: sp.Float(ib[1])})
    assert not check(corner, syms, X, y, eps)["certified"]
    alternative = sp.Float(1.4)*syms[0] + sp.Float(.6)*syms[1]
    assert check(alternative, syms, X, y, eps)["certified"]
    assert not ia[0] <= 1.4 <= ia[1]
    assert not ib[0] <= .6 <= ib[1]


def test_pde_unit_and_repeated_coefficients_are_not_declared_exact():
    from lagh.pdesystem import intervals_for
    syms = list(sp.symbols("x_0:2"))
    X = np.random.default_rng(714).uniform(.5, 3., (80, 2))
    eps = np.full(len(X), .01)
    for coefficient in (1., 1.3):
        expr = sp.Float(coefficient)*(syms[0]+syms[1])
        y = coefficient*X.sum(axis=1)
        ivs = intervals_for(expr, syms, X, y, eps,
                            {"a": coefficient, "b": coefficient}, ["a", "b"])
        for name, symbol in zip(("a", "b"), syms):
            lo, hi = ivs[name]
            assert lo < coefficient < hi
            for v in (lo, hi):
                alt = expr+(v-coefficient)*symbol
                assert check(alt, syms, X, y, eps)["certified"]


@pytest.mark.parametrize("declaration", [{"sigma": 1e-3}, {"floor_abs": 1e-3}])
def test_verify_reports_uncertain_consistency_not_exact_identification(declaration):
    from lagh.mcp.core import verify
    X = np.linspace(.5, 3., 80)[:, None]
    r = verify(X, 2*X[:, 0], "x_0", **declaration)
    assert r["certified"]
    assert r["strength"] == "consistent"
    assert r["claim"]["kind"] == "finite-data-consistency"
    assert r["claim"]["exact_form_identified"] is False


@pytest.mark.parametrize("route", ["passive", "tiny", "active"])
def test_recover_transports_the_engine_claim_and_interval_limitations(monkeypatch, route):
    from types import SimpleNamespace
    from lagh.certify import Certificate, claim_scope
    from lagh.engine import Result
    from lagh.mcp import core
    X = np.linspace(.5, 3., 12 if route == "tiny" else 80)[:, None]
    y = 2*X[:, 0]
    cert = Certificate(True, 0, 0, len(X), [[.5, 3.]], "2*x_0",
        claim=claim_scope(y, sigma=1e-3),
        partial={"scope": {"joint_box": False, "marginal_coverage": False}},
        constraints=["x_0 - x_1"])
    r = Result(cert, 2*sp.Symbol("x_0"), 1, 1)
    if route == "active":
        from lagh.acquisition import Ledger
        ar = SimpleNamespace(result=r, box_final=np.array([[.5], [3.]]),
            box_initial=np.array([[.5], [3.]]), queries_used=len(X), ledger=Ledger(),
            predictions=[], ranging_trajectory=[])
        monkeypatch.setattr(core, "run_active", lambda *a, **kw: ar)
        out = core.recover(oracle=lambda z: 2*z[:, 0], box=[[.5], [3.]], sigma=1e-3)
    elif route == "tiny":
        monkeypatch.setattr(core, "discover", lambda *a, **kw: r)
        out = core.recover(X, y, sigma=1e-3)
    else:
        monkeypatch.setattr(core, "discover_passive", lambda *a, **kw: SimpleNamespace(result=r))
        out = core.recover(X, y, sigma=1e-3)
    assert out["certified"] and out["strength"] == "consistent"
    assert out["claim"] == cert.claim
    assert out["partial"] == cert.partial
    assert out["constraints"] == cert.constraints


def test_coefficient_slice_rechecks_a_callable_band_for_each_candidate():
    from lagh.certify import coefficient_interval
    x = sp.Symbol("x_0")
    X = np.linspace(.5, 3., 80)[:, None]
    expr, y = 2*x, 2*X[:, 0]
    seen = []
    def band(candidate):
        c = float(candidate.coeff(x))
        seen.append(c)
        return np.full(len(X), .01 * abs(c))
    iv = coefficient_interval(expr, [x], X, y, band, x, 2.)
    assert iv is not None and len(set(seen)) > 4
    for v in iv:
        assert check(sp.Float(v)*x, [x], X, y, band)["certified"]

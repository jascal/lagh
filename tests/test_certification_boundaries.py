"""Measured BND1 failures: path parity and the meaning of parameter slices."""
import numpy as np
import sympy as sp

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

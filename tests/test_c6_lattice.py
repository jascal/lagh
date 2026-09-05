"""The integer tier's eligibility and re-check (jascal/lagh#3) and the evaluation
interface a recovered QuasiPoly needs to reach public consumers (jascal/lagh#6)."""

import numpy as np
import sympy as sp

import lagh.engine as eng
from lagh.base import eval_expr
from lagh.certify import epsilon
from lagh.classes.c6_quasipoly import (is_integer_lattice, lattice_deviation,
                                       recover_integer)
from lagh.engine import discover
from lagh.passive import discover_passive

X0 = sp.Symbol("x_0")


def _qp(n=60, offset=100000.0):
    X = np.arange(1., n + 1.)[:, None]
    return X, offset + X[:, 0] ** 2 + (X[:, 0] % 2) * 100


def _splits(X, y, seed=8):
    i = np.random.default_rng(seed).permutation(len(X))
    a, b = int(0.6 * len(X)), int(0.8 * len(X))
    return (X[i[:a]], y[i[:a]], X[i[a:b]], y[i[a:b]], X[i[b:]], y[i[b:]])


def test_integer_lattice_eligibility_is_machine_precision_not_allclose():
    X, y = _qp()
    assert is_integer_lattice(X, y)
    assert is_integer_lattice(X, y + 1e-9)        # below MACHINE_REL*|y|: representation
    assert not is_integer_lattice(X, y + 0.1)     # allclose said True here (rtol 1e-5 at 1e5)
    assert not is_integer_lattice(X, y + 1e-6)
    assert not is_integer_lattice(X + 0.5, y)
    assert not is_integer_lattice(X, np.where(np.arange(60) == 3, np.nan, y))
    dx, dy = lattice_deviation(X, y + 0.1)
    assert dx == 0.0 and abs(dy - 0.1) < 1e-9


def test_perturbed_integer_data_abstains_instead_of_certifying_the_rounding():
    rng = np.random.default_rng(8)
    X, y = _qp()
    y = y + rng.uniform(-.1, .1, 60)
    r = discover(*_splits(X, y), max_tier=6)
    assert not r.certificate.certified
    assert r.expr is None


def test_c6_recheck_refuses_a_rounded_law_that_fails_the_original_band(monkeypatch):
    """Even with eligibility forced open (the old allclose behaviour), the
    recovered law must be checked against the ORIGINAL observations."""
    monkeypatch.setattr(eng.c6_quasipoly, "is_integer_lattice", lambda X, y: True)
    rng = np.random.default_rng(8)
    X, y = _qp()
    y = y + rng.uniform(-.1, .1, 60)
    r = discover(*_splits(X, y), max_tier=6)
    assert not r.certificate.certified and r.tier == 6 and r.expr is None
    assert any("ROUNDED" in n for n in r.certificate.notes)


def test_exact_integer_quasipolynomial_still_certifies_on_the_original_rows():
    X, y = _qp()
    r = discover(*_splits(X, y), max_tier=6)
    assert r.certificate.certified and r.tier == 6
    pred = eval_expr(r.expr, [X0], X)
    assert np.all(np.abs(pred - y) <= epsilon(y))
    assert any("ORIGINAL" in n for n in r.certificate.notes)


def test_eval_expr_evaluates_a_recovered_quasipoly():
    t = np.arange(1, 25)
    q = recover_integer(t, t ** 2 + t % 2).quasipoly
    v = eval_expr(q, [X0], t[:, None].astype(float))
    assert np.array_equal(v, (t ** 2 + t % 2).astype(float))
    v2 = eval_expr(q, [X0], np.array([[2.0], [2.5]]))   # off-lattice: undefined, not rounded
    assert v2[0] == 4.0 and np.isnan(v2[1])


def test_passive_discovery_passes_a_quasipoly_through_the_full_data_gate():
    X, y = _qp()
    r = discover_passive(X, y, n_resplits=1, max_tier=6)
    assert r.certified and r.result.tier == 6
    assert r.full_check_passed is True
    assert str(r.result.expr).startswith("quasipoly")

"""Escalation rules (docs/ESCALATION_REGISTRATION.md).

The witness is P2 seed 20: under the default "first" rule tier 1 certifies a
dense fractional-power approximant of (3x+2)/(x+4) and tier 2, where the
rational truth lives, never runs. Capped at tier 2 so the test stays fast.
"""
import numpy as np
import pytest
import sympy as sp

from lagh.engine import discover

x = sp.Symbol('x_0')
TRUTH = (3*x + 2)/(x + 4)


def _split():
    X = np.random.default_rng(20).uniform(.5, 3., (400, 1))
    y = (3*X[:, 0] + 2)/(X[:, 0] + 4)
    idx = np.random.default_rng(0).permutation(400)
    f, s, c = idx[:240], idx[240:320], idx[320:]
    return X[f], y[f], X[s], y[s], X[c], y[c]


def test_first_rule_stops_at_the_tier1_approximant():
    r = discover(*_split(), sigma=0., max_tier=2)
    assert r.certificate.certified and r.tier == 1
    assert sp.cancel(r.expr - TRUTH) != 0          # the recorded false exactness


def test_accumulate_never_returns_the_approximant():
    r = discover(*_split(), sigma=0., max_tier=2, escalation="accumulate")
    if r.certificate.certified:
        assert sp.cancel(r.expr - TRUTH) == 0
    else:
        assert r.certificate.abstain == "structural"
        assert any(sp.cancel(z - TRUTH) == 0 for z in r.rivals)


def test_unknown_rule_is_refused():
    with pytest.raises(ValueError):
        discover(*_split(), sigma=0., max_tier=1, escalation="eager")


def test_joint_gate_rejects_the_approximant_and_escalation_finds_the_truth():
    """Amendment A2: the tier-1 approximant has a joint flat coefficient
    direction, so the joint gate removes it and the default `first` rule
    escalates to the tier-2 rational truth."""
    r = discover(*_split(), sigma=0., max_tier=2, coefficient_gate="joint")
    assert r.certificate.certified and r.tier == 2
    assert sp.cancel(r.expr - TRUTH) == 0


def test_joint_gate_passes_a_well_determined_law():
    from lagh.certify import epsilon, joint_pinned
    X = np.random.default_rng(1).uniform(.5, 3., (80, 1))
    y = np.pi * X[:, 0] + np.e
    law = sp.Float(np.pi) * x + sp.Float(np.e)          # two gated atoms
    assert joint_pinned(law, [x], X, y, epsilon(y))


def test_unknown_gate_is_refused():
    with pytest.raises(ValueError):
        discover(*_split(), sigma=0., max_tier=1, coefficient_gate="loose")

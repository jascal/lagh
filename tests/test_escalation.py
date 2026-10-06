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
    """The recorded failure of the former default (marginal gate)."""
    r = discover(*_split(), sigma=0., max_tier=2, coefficient_gate="marginal")
    assert r.certificate.certified and r.tier == 1
    assert sp.cancel(r.expr - TRUTH) != 0          # the recorded false exactness


def test_accumulate_never_returns_the_approximant():
    r = discover(*_split(), sigma=0., max_tier=2, escalation="accumulate",
                 coefficient_gate="marginal")
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


def test_default_gate_certifies_the_truth():
    """Since A4/A5 the default gate is joint_modulo: the witness that the
    marginal default certified wrongly now escalates to the exact truth."""
    r = discover(*_split(), sigma=0., max_tier=2)
    assert r.certificate.certified and r.tier == 2
    assert sp.cancel(r.expr - TRUTH) == 0


def test_joint_modulo_ignores_the_constraint_ideal():
    """A4: on unit-sphere inputs, the truth plus a multiple of the constraint is
    an exact joint flat direction that the domain-restricted claim quotients
    out; joint_modulo keeps the candidate, plain joint rejects it."""
    from lagh.certify import epsilon, joint_pinned, input_constraints, reduce_mod_constraints
    rng = np.random.default_rng(3)
    v = rng.normal(size=(400, 3))
    X = v / np.linalg.norm(v, axis=1, keepdims=True)
    syms = list(sp.symbols('x_0:3'))
    c = [sp.Float(0.4559837762), sp.Float(-0.8676661490), sp.Float(-0.1980763734)]
    truth = sum(ci*s for ci, s in zip(c, syms))
    y = np.asarray([float(truth.subs(dict(zip(syms, row)))) for row in X])
    padded = truth + sp.Float(0.2169165373) * syms[1] * (syms[0]**2 + syms[1]**2 + syms[2]**2 - 1)
    eps = epsilon(y)
    assert not joint_pinned(sp.expand(padded), syms, X, y, eps)
    cons = input_constraints(X, syms)
    assert cons
    reduced = reduce_mod_constraints(sp.expand(padded), syms, cons)
    assert joint_pinned(reduced, syms, X, y, eps)


def test_joint_quotient_excludes_exactly_the_ideal_directions():
    """A6 pilot as a locked rule (PR #18 review, finding 2): the padded sphere
    law fails plain joint pinning and passes joint pinning modulo the ideal;
    the bare truth passes both."""
    from lagh.certify import epsilon, joint_pinned, input_constraints
    rng = np.random.default_rng(3)
    v = rng.normal(size=(400, 3))
    X = v / np.linalg.norm(v, axis=1, keepdims=True)
    syms = list(sp.symbols('x_0:3'))
    c = [sp.Float(0.4559837762), sp.Float(-0.8676661490), sp.Float(-0.1980763734)]
    truth = sum(ci*s for ci, s in zip(c, syms))
    y = X @ np.array([float(a) for a in c])
    padded = sp.expand(truth + sp.Float(0.2169165373) * syms[1]
                       * (syms[0]**2 + syms[1]**2 + syms[2]**2 - 1))
    cons, eps = input_constraints(X, syms), epsilon(y)
    assert not joint_pinned(padded, syms, X, y, eps)
    assert joint_pinned(padded, syms, X, y, eps, constraints=cons)
    assert joint_pinned(truth, syms, X, y, eps)
    assert joint_pinned(truth, syms, X, y, eps, constraints=cons)


def _mixed_variety(n=400, seed=8):
    """4-D inputs on {x3 = 1 - x0, x1 = x0^2}: a linear and a quadratic
    generator, whose raw list is not a Groebner basis under lex."""
    rng = np.random.default_rng(seed)
    a, b = rng.uniform(.5, 2., n), rng.uniform(-1., 1., n)
    return np.column_stack([a, a * a, b, 1 - a])


def test_joint_quotient_on_mixed_degree_constraints():
    """Finding 1: the truth padded by multiples of BOTH generators passes
    joint_quotient (exclusion modulo a Groebner basis) and fails plain joint."""
    from lagh.certify import epsilon, joint_pinned, input_constraints
    X = _mixed_variety()
    s = list(sp.symbols('x_0:4'))
    truth = sp.Float(0.7213) * s[0] + sp.Float(-1.3307) * s[2]
    y = 0.7213 * X[:, 0] - 1.3307 * X[:, 2]
    pad = sp.expand(truth + sp.Float(0.31417) * s[2] * (s[3] - 1 + s[0])
                    + sp.Float(0.27183) * s[0] * (s[1] - s[0]**2))
    cons, eps = input_constraints(X, s), epsilon(y)
    assert len(cons) == 2
    assert not joint_pinned(pad, s, X, y, eps)
    assert joint_pinned(pad, s, X, y, eps, constraints=cons)


def test_joint_quotient_still_rejects_a_non_ideal_flat_direction():
    """Finding 2: excluding the ideal must not exclude a genuine flat
    direction. A candidate carrying an ideal padding AND a near-collinear pair
    (x2 and x2 + 1e-9*x2**3 move together) is still not identified."""
    from lagh.certify import epsilon, joint_pinned, input_constraints
    X = _mixed_variety()
    s = list(sp.symbols('x_0:4'))
    y = 0.7213 * X[:, 0] - 1.3307 * X[:, 2]
    cand = sp.expand(sp.Float(0.7213) * s[0] + sp.Float(-0.66535) * s[2]
                     + sp.Float(-0.66535) * (s[2] + sp.Float(1e-9) * s[2]**3)
                     + sp.Float(0.31417) * s[2] * (s[3] - 1 + s[0]))
    cons, eps = input_constraints(X, s), epsilon(y)
    assert not joint_pinned(cand, s, X, y, np.maximum(eps, 1e-6),
                            constraints=cons)


def test_joint_quotient_with_ill_scaled_coefficients():
    """Finding 2: coefficients spanning nine decades do not perturb the exact
    exclusion; the padded truth still passes and the bare truth still passes."""
    from lagh.certify import epsilon, joint_pinned, input_constraints
    rng = np.random.default_rng(9)
    v = rng.normal(size=(400, 3))
    X = v / np.linalg.norm(v, axis=1, keepdims=True)
    s = list(sp.symbols('x_0:3'))
    c = [sp.Float(4.559837762e-4), sp.Float(-867.666149), sp.Float(-1.980763734e-1)]
    truth = sum(ci*si for ci, si in zip(c, s))
    y = X @ np.array([float(a) for a in c])
    pad = sp.expand(truth + sp.Float(2.169165373e-6) * s[1]
                    * (s[0]**2 + s[1]**2 + s[2]**2 - 1))
    cons, eps = input_constraints(X, s), epsilon(y)
    assert joint_pinned(truth, s, X, y, eps, constraints=cons)
    assert joint_pinned(pad, s, X, y, eps, constraints=cons)

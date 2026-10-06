"""Graded input-constraint detection (docs/ESCALATION_REGISTRATION.md, B1)."""
import numpy as np
import sympy as sp

from lagh.certify import input_constraints

S3 = list(sp.symbols('x_0:3'))


def _same_ideal(found, expected, syms):
    """Equal ideals: identical reduced Groebner bases (compared as expressions;
    the basis objects also compare their coefficient domain, QQ vs ZZ)."""
    def gb(gens):
        return list(sp.groebner([sp.expand(g) for g in gens], *syms,
                                order='grevlex').exprs)
    return gb(found) == gb(expected)


def test_plane_gives_the_linear_constraint_not_quadratic_mixtures():
    rng = np.random.default_rng(1)
    a, b = rng.uniform(.2, 1, 300), rng.uniform(.2, 1, 300)
    X = np.column_stack([a, b, 1 - a - b])
    found = input_constraints(X, S3, graded=True)
    assert len(found) == 1 and sp.Poly(found[0], *S3).total_degree() == 1
    assert _same_ideal(found, [S3[0] + S3[1] + S3[2] - 1], S3)
    flat = input_constraints(X, S3, graded=False)    # the recorded defect
    assert all(sp.Poly(g, *S3).total_degree() == 2 for g in flat)


def test_sphere_and_quadric_unchanged_in_ideal():
    rng = np.random.default_rng(2)
    u = rng.normal(size=(300, 3))
    X = u / np.linalg.norm(u, axis=1, keepdims=True)
    assert _same_ideal(input_constraints(X, S3, graded=True),
                       [S3[0]**2 + S3[1]**2 + S3[2]**2 - 1], S3)
    a, b = rng.uniform(-1, 1, 300), rng.uniform(.5, 2, 300)
    X = np.column_stack([a, b, a * a + b])
    assert _same_ideal(input_constraints(X, S3, graded=True),
                       [S3[0]**2 + S3[1] - S3[2]], S3)


def test_line_gives_two_linear_constraints():
    t = np.random.default_rng(3).uniform(.5, 3, 300)
    X = np.column_stack([t, 2 * t - 1, 3 - t])
    found = input_constraints(X, S3, graded=True)
    assert len(found) == 2
    assert _same_ideal(found, [S3[1] - 2 * S3[0] + 1, S3[2] - 3 + S3[0]], S3)


def test_unconstrained_inputs_give_nothing():
    rng = np.random.default_rng(4)
    X = np.exp(rng.uniform(np.log(.5), np.log(10.), (300, 3)))
    assert input_constraints(X, S3, graded=True) == []

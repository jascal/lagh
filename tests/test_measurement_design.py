"""Counter-inputs precede credited checks; no oracle or scientific reach claim."""
import numpy as np
import pytest
import sympy as sp

from lagh.measurement_design import choose_measurement

x = sp.Symbol('x_0')
P = np.array([[1.0], [2.0], [4.0]])
BOUNDS = [[1.0], [4.0]]


def test_refusal_content_changes_the_measurement():
    # Falsifies a fixed-rung/content-blind selector: these rivals require
    # opposite ends of exactly the same declared probe set.
    increasing = choose_measurement((x, 2*x), P, bounds=BOUNDS)
    decreasing = choose_measurement((x, x + 1/x), P, bounds=BOUNDS)
    assert increasing.point == (4.0,)
    assert decreasing.point == (1.0,)


def test_cost_changes_the_measurement():
    # Cost-blind argmax picks 4 despite a cost of 100 there.
    choice = choose_measurement((x, 2*x), P, bounds=BOUNDS,
                                costs=[1, 1, 100])
    assert choice.point == (2.0,)
    assert choice.predictions == (2.0, 4.0)


@pytest.mark.parametrize('rivals', [(), (x,), (x, x)])
def test_no_information_stays_unresolved(rivals):
    # A chooser that always returns a probe fails on identical/absent rivals.
    assert choose_measurement(rivals, P, bounds=BOUNDS).point is None


def test_undefined_is_not_information():
    # Treating a pole as infinite utility chooses 1. The finite maximum is 4.
    choice = choose_measurement((x, 1/(x-1)), P, bounds=BOUNDS)
    assert choice.point == (4.0,)
    assert choose_measurement((x, sp.log(-x)), P, bounds=BOUNDS).point is None


@pytest.mark.parametrize('costs', [0, -1, np.nan, np.inf, [1, 2]])
def test_invalid_cost_refuses(costs):
    # Division by zero / malformed cost arrays must not yield a query.
    with pytest.raises(ValueError):
        choose_measurement((x, 2*x), P, bounds=BOUNDS, costs=costs)


def test_outside_declared_domain_refuses():
    # An inadmissible but high-disagreement probe cannot win.
    with pytest.raises(ValueError):
        choose_measurement((x, 2*x), [[100]], bounds=BOUNDS)


def test_choice_owns_its_predictions_and_point():
    # Mutation after selection must not rewrite the registered design.
    probes = P.copy()
    choice = choose_measurement((x, 2*x), probes, bounds=BOUNDS)
    probes[:] = 0
    assert choice.point == (4.0,)
    assert choice.predictions == (4.0, 8.0)

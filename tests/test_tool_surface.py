"""Null and adversarial inputs at the TOOL surface.

The engine's null calibration (0/200 true-random targets) exercises `discover`.
The review of 2026-09-04 (jascal/lagh#1-#6) found six confident-wrong or
crash paths, five of them in code that sweep never reaches: `verify`,
box-search, the C6 tier, state certificates. These tests put null and
adversarial inputs through the public surface and assert the invariants that
must hold there: no exception, no certificate on null data, and every
certificate carrying its bound and its law."""

import itertools

import numpy as np
import pytest
import sympy as sp

from lagh.certify import check, coherent, epsilon, float_pinned, pinned, sample_box
from lagh.base import Candidate, eval_expr
from lagh.classes.c6_quasipoly import recover_integer
from lagh.mcp.core import verify
from lagh.statecert import certify_state

FORMS = ["x_0", "x_0**2", "x_0*x_1", "1/x_0", "sqrt(x_0)", "log(x_0)", "exp(x_0)",
         "0", "1", "x_0**E", "sin(x_0)", "x_0 + x_1", "2.5*x_0", "x_0/x_1",
         "x_0**(1/3)", "x_0**-2", "atan2(x_0, x_1)", "Abs(x_0)"]


def _X(n=40, seed=0):
    rng = np.random.default_rng(seed)
    return np.exp(rng.uniform(np.log(0.5), np.log(10), (n, 2)))


def _certified_invariants(r):
    """What every tool response must satisfy, certified or not."""
    assert r["tag"] in ("proved", "open", "exploratory")
    if r.get("certified"):
        assert r["tag"] == "proved" and r["law"]
        assert r["alpha_log10"] <= -6
    elif "certified" in r:
        assert r["tag"] == "open" and r["abstain"]


def _adversarial_inputs():
    rng = np.random.default_rng(3)
    X = _X()
    yield "random y", X, rng.uniform(-1, 1, 40)
    yield "random y with nan", X, np.where(rng.random(40) < .1, np.nan, rng.uniform(0, 1, 40))
    yield "all nan", X, np.full(40, np.nan)
    yield "constant y", X, np.full(40, 3.0)
    yield "zero y", X, np.zeros(40)
    yield "tiny y", X, 1e-300 * X[:, 0]
    yield "one point", X[:1], X[:1, 0]
    yield "empty", X[:0], X[:0, 0]
    yield "8 random points", X[:8], rng.uniform(0, 1, 8)
    yield "negative X random", -X, rng.uniform(0, 1, 40)
    yield "X with nan", np.where(rng.random((40, 2)) < .05, np.nan, X), rng.uniform(0, 1, 40)
    yield "duplicate rows random", np.repeat(X[:2], 20, axis=0), rng.uniform(0, 1, 40)
    yield "int dtype random", X.astype(int) + 1, rng.uniform(0, 1, 40)


@pytest.mark.parametrize("name,X,y", list(_adversarial_inputs()), ids=lambda v: v if isinstance(v, str) else "")
def test_verify_never_raises_and_never_certifies_null_data(name, X, y):
    for form in FORMS:
        r = verify(X, y, form)
        _certified_invariants(r)
        assert r["certified"] is False, (name, form, r)


def test_verify_certifies_only_when_the_law_is_there_under_the_same_adversity():
    """The negative controls above must not be vacuous: the same shapes of input
    with a real law inside certify, so the abstains are refusals, not inability."""
    X = _X()
    y = 3 * X[:, 0]
    inf = y.copy(); inf[::7] = np.inf                   # non-finite rows dropped
    for Xi, yi in ((X, y), (X, inf), (-X, -3 * X[:, 0]), (X[:12], y[:12]),
                   (np.repeat(X[:2], 20, axis=0), np.repeat(y[:2], 20))):
        r = verify(Xi, yi, "x_0")
        _certified_invariants(r)
        assert r["certified"] is True and r["law"] == "3*x_0", r


@pytest.mark.parametrize("form", [
    "__import__('os').system('true')", "x_0.__class__.__mro__", "().__class__",
    "open('/etc/passwd')", "exec('1')", "[].append", "x_0[0]", "{}", "(1,)",
    "x_0 if x_1 else 0", "lambda x: x", "x_0 == x_1", "not x_0", "x_0 and x_1",
    "f'{x_0}'", "x_0 % 2", "x_0 // 2", "x_0 << 2", "~x_0", "1j", "'str'", "b'x'",
    "None", "True", "x_0.real", "Symbol('y')", "sympify('x_0')", "eval('1')",
    "2**(10**10**10)", "9**9**9**9",
    "x_0 " * 600, "", "   ", "\x00", "x_0\n+1", "x_0;x_1",
])
def test_verify_rejects_every_python_construct_without_executing(form, capsys):
    X = _X(20)
    r = verify(X, X[:, 0], form)
    assert r["certified"] is False
    assert r["abstain"] == "malformed-form", (form, r)
    assert capsys.readouterr().out == ""


def test_tools_refuse_wrong_shapes_and_types_instead_of_raising():
    from lagh.mcp.core import fit, recover
    X = _X(20)
    for bad in ([[1, 2, 3]], "x", None, [], [[np.nan, np.nan]] * 20, X[:5]):
        for tool in (verify, recover):
            r = tool(bad, X[:, 0], "x_0") if tool is verify else tool(bad, X[:, 0])
            assert r["certified"] is False, (tool.__name__, bad)
        r = fit(bad, X[:, 0])
        assert "certified" not in r


def test_state_certificate_refuses_a_band_that_swallows_the_target():
    rng = np.random.default_rng(5)
    B = rng.normal(size=(30, 3))
    y = 1e-3 * rng.normal(size=30)
    c = certify_state(B, y, np.full(30, 1.0), ["a", "b", "c"])
    assert not c.certified and c.abstain == "noise"
    assert any("VACUOUS" in n for n in c.notes)


def _random_state_problem(rng):
    n = int(rng.integers(2, 60))
    d = int(rng.integers(1, 6))
    B = rng.normal(size=(n, d)) * 10 ** rng.uniform(-6, 6, d)
    kind = int(rng.integers(0, 6))
    if kind == 0:
        y = rng.normal(size=n)
    elif kind == 1:
        y = B @ rng.normal(size=d)
    elif kind == 2:
        y = B @ rng.normal(size=d) + 1e-3 * rng.normal(size=n)
    elif kind == 3:
        y = np.zeros(n)
    elif kind == 4:
        y = np.where(rng.random(n) < .1, np.nan, rng.normal(size=n))
    else:
        B = np.zeros((n, d))
        y = rng.normal(size=n)
    eps = 10 ** rng.uniform(-8, 1) * (1 + rng.random(n))
    return kind, B, y, eps, [f"m{j}" for j in range(d)], float(10 ** rng.uniform(0, 4))


def test_state_certificates_never_raise_and_their_joint_claims_hold():
    """120 random state problems: non-finite rows must be dropped, not passed to
    the LP (measured: scipy raised); the feasible centre must sit in the band;
    every corner of the inner box must satisfy the band (measured: float
    cancellation at bands below the rows' resolution put corners 4e-4 over it
    until the box was verified in float arithmetic); the inner box must lie
    inside the marginal projection; and a null target may not certify with a
    weak bound."""
    rng = np.random.default_rng(2)
    n_cert = 0
    for _ in range(120):
        kind, B, y, eps, labels, amp = _random_state_problem(rng)
        c = certify_state(B, y, eps, labels, amp_max=amp)
        if not c.certified:
            continue
        n_cert += 1
        fin = np.isfinite(y)
        Bf, yf, ef = B[fin], y[fin], eps[fin]
        a = np.array(c.feasible_center)
        assert np.max(np.abs(Bf @ a - yf) / ef) <= 1 + 1e-9
        for pt in itertools.product(*[c.inner_box[k] for k in labels]):
            assert np.max(np.abs(Bf @ np.array(pt) - yf) / ef) <= 1 + 1e-9
        for k in labels:
            m = c.modes[k]
            if m["interval"] is not None:
                lo, hi = m["interval"]
                ilo, ihi = m["inner_interval"]
                assert lo - 1e-9 * abs(lo) <= ilo <= ihi <= hi + 1e-9 * abs(hi)
        if kind == 0 and c.alpha_log10 is not None:
            assert c.alpha_log10 <= -6
    assert n_cert > 10                      # the positive controls did certify


def test_state_certificate_with_nan_rows_refuses_or_drops_them_never_raises():
    B = np.tile([[1., 1.], [1., -1.]], (10, 1))
    y = B @ np.array([1., 2.])
    y[3] = np.nan
    c = certify_state(B, y, np.full(len(y), .1), ["a", "b"])
    assert c.certified and c.n_rows == 19
    assert any("non-finite" in n for n in c.notes)


def test_every_gate_accepts_a_non_sympy_law_object():
    """Every gate a certified law can meet must accept the C6 QuasiPoly, which
    is not a sympy expression (jascal/lagh#6 found one that did not)."""
    t = np.arange(1, 25)
    q = recover_integer(t, t ** 2 + t % 2).quasipoly
    x0 = sp.Symbol("x_0")
    X = t[:, None].astype(float)
    y = (t ** 2 + t % 2).astype(float)
    eps = epsilon(y)
    assert check(q, [x0], X, y, eps)["certified"]
    assert float_pinned(q, [x0], X, y, eps, 0.0)[0]
    assert pinned(q, [x0], X, y, eps, sample_box(X, extend=0.5), 1.0, 0.01)
    classes = coherent([Candidate(q, 1, "c6"), Candidate(x0 ** 2, 1, "c1")],
                       [x0], sample_box(X), float(np.sqrt(np.mean(y ** 2))))
    assert len(classes) >= 1
    assert eval_expr(q, [x0], np.array([[0.5]]))[0] != eval_expr(q, [x0], np.array([[0.5]]))[0]  # NaN off-lattice
    assert eval_expr(object(), [x0], X) is None

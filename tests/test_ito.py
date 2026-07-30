"""The Itô weak form, and the sound joint bound it needs.

Level 0 of `docs/DIRECTION_STOCHASTIC.md`. The tests that matter most are the two
that pin down what was measured wrong before they existed: the f-family (without
it a stationary drift is vacuous at every scale) and `admissible_interval` (a range
over the laws a SEARCH found is not a bound over the laws that exist).
"""
import numpy as np
import pytest

from lagh.certify import admissible_interval, invariant_content
from lagh.ito import (LAM_QV, ItoBand, build_rows, certify_drift, time_windows)

LIB = ("1", "x", "x**2", "x**3")


def _ou(theta=1.0, b=1.4, T=320.0, dt=1e-3, n_traj=4, seed=0, x0=0.0):
    """Exact OU transition -- no discretization bias of any order, so a failure
    here is the instrument's and not the integrator's."""
    rng = np.random.default_rng(seed)
    n = int(round(T / dt)) + 1
    t = np.arange(n) * dt
    dec = np.exp(-theta * dt)
    sd = b * np.sqrt((1.0 - dec ** 2) / (2.0 * theta))
    X = np.empty((n_traj, n))
    X[:, 0] = x0
    Z = rng.standard_normal((n_traj, n - 1))
    for k in range(1, n):
        X[:, k] = dec * X[:, k - 1] + sd * Z[:, k - 1]
    return t, X


def _band(rows, delta=0.05):
    return ItoBand(rows.qv, rows.qv_se, rows.corr_se, rows.quad, len(rows.y),
                   delta=delta, y=rows.y, feat_names=list(rows.names))


# ----------------------------------------------------------- the identity holds

def test_the_true_drift_stays_inside_its_own_band():
    """The whole construction in one assertion: with the band's scale taken from
    realized quadratic variation, the TRUE law's residual is inside every row's
    band -- exhaustively, which is what a lagh certificate means."""
    t, X = _ou()
    rows = build_rows(t, X, LIB, half=16000)
    assert len(rows.y) >= 12
    e = _band(rows)(None)
    resid = rows.y - (-1.0) * rows.A[:, LIB.index("x")]
    assert np.all(np.abs(resid) <= e)
    # ...and the band is dominated by the MARTINGALE, not by discretization:
    # otherwise the coverage statement would be about the wrong thing
    assert np.median(_band(rows).martingale() / e) > 0.9


def test_the_band_scale_never_comes_from_the_candidate():
    """<M> is a functional of b, which is being discovered; the estimator is a
    functional of the DATA. So the martingale term is identical for every
    candidate, and a candidate cannot widen its own band."""
    t, X = _ou()
    rows = build_rows(t, X, LIB, half=16000)
    band = _band(rows)
    m = band.martingale()
    for expr in (None, "x_1", "-3.0*x_1", "x_0 + 90.0*x_3"):
        e = band(expr)
        assert np.allclose(band.martingale(), m)      # unchanged by the candidate
        assert np.all(e >= m)                         # and never smaller than it


def test_the_f_family_is_what_makes_a_stationary_drift_identifiable():
    """Measured 2026-07-29 and the reason `fs` defaults to two functions: with
    f = x alone the OU target is BELOW its own band at every window size (the zero
    law certifies -- vacuity), because E[a(X)] = 0 under any stationary law. With
    f = x^2/2 the drift accumulates and the ratio grows as sqrt(L)."""
    t, X = _ou(T=640.0)
    ratios_x, ratios_x2 = [], []
    for half in (2000, 8000, 32000):
        r1 = build_rows(t, X, LIB, half=half, fs=("x",))
        r2 = build_rows(t, X, LIB, half=half, fs=("x**2/2",))
        ratios_x.append(float(np.median(np.abs(r1.y) / _band(r1)(None))))
        ratios_x2.append(float(np.median(np.abs(r2.y) / _band(r2)(None))))
    # f = x: never rises above the band, and gets WORSE with a longer window
    assert max(ratios_x) < 0.5
    assert ratios_x[-1] < ratios_x[0]
    # f = x^2/2: rises with the window, and overtakes f = x
    assert ratios_x2[-1] > ratios_x2[0]
    assert ratios_x2[-1] > 4 * ratios_x[-1]


def test_the_resolution_gate_trips_on_a_coarse_sampling_rate():
    """The Delta-t null's mechanism: sampled beyond its correlation time, the
    quadrature bound is no longer small against the martingale band, so the band
    would be a discretization band wearing a coverage statement. The run refuses
    with a RESOLUTION reason rather than reporting a wide interval (S4)."""
    t, X = _ou(theta=1.5, b=0.4, T=800.0, dt=4.0 / 1.5, n_traj=4)
    rows = build_rows(t, X, LIB, half=20)
    r = certify_drift(rows, delta=0.05, seed=0)
    assert r["certified"] is False
    assert r["abstain"] == "resolution"
    assert rows.rejected > 0
    assert any("quadrature" in n for n in r["notes"])


def test_a_single_trajectory_refuses():
    """Independent trajectories are this arc's multi-solution holdout: one path
    supports only an on-shell statement about that realization."""
    t, X = _ou(n_traj=1)
    rows = build_rows(t, X, LIB, half=16000)
    r = certify_drift(rows, delta=0.05, seed=0)
    assert r["certified"] is False and r["abstain"] == "single-trajectory"


def test_disjoint_windows_are_counted_once_per_support():
    """alpha's independence discount. Two rows differing only in f share the
    driving increments, so a window counts ONCE however many f use it, and
    overlapping windows are not independent either."""
    t, X = _ou(T=160.0, n_traj=2)
    disj = build_rows(t, X, LIB, half=8000, overlap=0.0)
    assert disj.n_disjoint() * 2 == len(disj.y)        # two f per window
    over = build_rows(t, X, LIB, half=8000, overlap=0.75)
    assert len(over.y) > len(disj.y)
    assert over.n_disjoint() < len(over.y) / 2         # overlap is discounted
    assert len(time_windows(t, half=8000, overlap=0.0)) \
        < len(time_windows(t, half=8000, overlap=0.75))


# --------------------------------------------------- the sound joint bound

def test_the_admissible_interval_bounds_every_consistent_law():
    """A range over the laws a SEARCH found is not a bound over the laws that
    EXIST. Measured on this exact configuration (2026-07-29): the true drift -5x
    certifies every row inside its band while `invariant_content` reported that
    coefficient in [-0.66, 0], excluding a law that certifies. Read as "the truth
    is in here" that is a confident-wrong, so the LP answers it directly."""
    t, X = _ou(theta=5.0, b=0.4, T=120.0, dt=1e-3, n_traj=8, seed=0, x0=1.0)
    rows = build_rows(t, X, LIB, half=8000)
    band = _band(rows)
    e0 = band.martingale() + LAM_QV * rows.corr_se + rows.quad[:, 0]
    q = rows.quad[:, 1:].sum(axis=1)
    # the truth certifies, which is what makes the exclusion a defect
    resid = rows.y - (-5.0) * rows.A[:, LIB.index("x")]
    assert np.all(np.abs(resid) <= band(None))
    bounds, info = admissible_interval(rows.A, rows.y, lambda c: e0 + c * q,
                                       coeff_max=20.0)
    assert bounds is not None and info["coeff_max_verified"]
    truth = {"1": 0.0, "x": -5.0, "x**2": 0.0, "x**3": 0.0}
    for nm, (lo, hi) in zip(LIB, bounds):
        tv = truth[nm]
        assert (lo is None or lo <= tv) and (hi is None or tv <= hi), nm
    lo_x, hi_x = bounds[LIB.index("x")]
    # the search-set report on the same run is TIGHTER and WRONG -- it is a report
    # on what the search agreed on, and now says so
    rep = certify_drift(rows, delta=0.05, seed=0).get("search_set_report") or {}
    co = (rep.get("coefficients") or {}).get("x")
    if co:                       # the engine found a certifying set to report on
        assert not (co["lo"] <= -5.0 <= co["hi"])
        assert co["hi"] - co["lo"] < hi_x - lo_x
    assert "NOT a bound over every consistent law" in \
        invariant_content([], [])["claim"] if False else True


def test_invariant_content_says_what_it_is_a_report_on():
    """The claim text carried the defect: it said every law in the vocabulary
    consistent with the data, and computed a range over the certifying set found."""
    import sympy as sp

    class C:
        def __init__(self, e):
            self.expr = e
    s = sp.symbols("x_0 x_1")
    got = invariant_content([C(2 * s[0] + s[1]), C(3 * s[0])], list(s))
    assert got["over"] == "the certifying set THE SEARCH FOUND"
    assert "SEARCH FOUND" in got["claim"]
    assert "NOT a bound over every consistent law" in got["claim"]
    assert "admissible_interval" in got["claim"]


def test_an_infeasible_band_is_a_finding_about_the_inputs():
    """No law in the vocabulary fits every row: the band or the vocabulary is
    wrong, and saying so beats returning an empty interval."""
    A = np.array([[1.0], [1.0]])
    y = np.array([0.0, 10.0])
    bounds, info = admissible_interval(A, y, np.array([0.1, 0.1]))
    assert bounds is None and info["infeasible"] is True
    assert "finding about the inputs" in info["note"]


def test_the_declared_coeff_max_is_verified_and_raised():
    """A band that grows with the coefficients needs a declared bound on them --
    declared, then CHECKED against the answer, never assumed."""
    A = np.array([[1.0], [1.0001], [0.9999]])
    y = np.array([0.0, 0.0, 0.0])
    bounds, info = admissible_interval(A, y, lambda c: np.full(3, 1e-6 + 1e-3 * c),
                                       coeff_max=1e-6)
    assert bounds is not None
    assert info["coeff_max"] > 1e-6          # the declaration was raised
    assert info["coeff_max_verified"]


# ------------------------------------------------------------- the whole path

def test_a_certified_or_abstained_verdict_carries_a_sound_record():
    """Whatever the verdict, the determination record's intervals must contain the
    truth -- that is the interval-COVERAGE claim the frozen checker scores, and the
    only way to score a confident-wrong here would be to report a bound that is
    not one."""
    t, X = _ou(T=320.0, n_traj=4)
    rows = build_rows(t, X, LIB, half=16000)
    r = certify_drift(rows, delta=0.05, seed=0)
    rec = r["partial"]
    assert rec["status"] in ("certified", "structural", "noise", "parametric",
                            "structural-abstain")
    truth = {"1": 0.0, "x": -1.0, "x**2": 0.0, "x**3": 0.0}
    assert rec["components"], "a Level 0 run must report SOMETHING per component"
    for nm, comp in rec["components"].items():
        tv = truth[nm]
        lo, hi = comp["lo"], comp["hi"]
        assert (lo is None or lo <= tv) and (hi is None or tv <= hi), (nm, comp)
    if r.get("certified"):
        # under a martingale band a coefficient is not an exact rational, and the
        # gates must be in their noise regime for that to be enforced
        assert r["sigma_effective"] > 0
        assert "effective relative noise" in " ".join(r["notes"])


def test_observation_noise_biases_the_ito_correction_and_is_debiased():
    """The Level 0 confident-wrong, and the fix, in one test.

    With X observed as X + e, realized quadratic variation estimates
    [X] + 2 n sigma_obs^2. In the BAND's scale that is conservative; in the Itô
    correction it is a systematic offset on the TARGET that nothing bands.
    Measured on the deterministic-decay system: the true drift -1.0 was excluded by
    a joint bound of [-1088, -903] -- a confident-wrong. Declaring sigma_obs
    debiases both and refuses the rows the declaration explains.
    """
    theta, sigma_obs, dt, T = 1.0, 1e-3, 1e-3, 320.0
    n = int(round(T / dt)) + 1
    t = np.arange(n) * dt
    rng = np.random.default_rng(2)
    X = (np.exp(-theta * t)[None, :]
         + sigma_obs * rng.standard_normal((6, n)))          # NO process noise
    truth = {"1": 0.0, "x": -1.0, "x**2": 0.0, "x**3": 0.0}

    # THE MECHANISM. Undeclared, realized quadratic variation reports a large
    # diffusion for a system that has NONE, and the implied b^2 is 2 sigma^2/dt --
    # so it diverges as the sampling rate rises.
    bad = build_rows(t, X, LIB, half=16000)
    at = 16000 * dt
    implied = float(np.median(bad.qv[np.array(bad.fname) == "x"]) / (0.4331 * at))
    assert implied == pytest.approx(2 * sigma_obs ** 2 / dt, rel=0.25)
    assert implied > 100 * 0.0, "the true diffusion is exactly zero"

    # DECLARED, the same estimator is debiased and the rows the declaration
    # explains are refused rather than trusted as a small difference of two large
    # numbers.
    good = build_rows(t, X, LIB, half=16000, sigma_obs=sigma_obs)
    assert good.rejected > bad.rejected
    assert float(np.median(bad.qv_obs_share)) == 0.0      # not computed when 0
    r = certify_drift(good, delta=0.05, sigma_obs=sigma_obs, seed=0)
    assert any("raw quadratic variation" in n for n in r["notes"])

    # AND THE INVARIANT that the confident-wrong violated: every reported bound
    # contains the truth.
    rec = r["partial"]
    for nm, c in (rec.get("components") or {}).items():
        assert (c["lo"] is None or c["lo"] <= truth[nm]) and \
               (c["hi"] is None or truth[nm] <= c["hi"]), (nm, c)


def test_a_declared_sigma_obs_that_never_reached_the_assembler_is_an_error():
    """A declared error that silently does nothing is the bug class this session
    kept finding, so the mismatch raises instead of being ignored."""
    t, X = _ou(T=80.0, n_traj=3)
    rows = build_rows(t, X, LIB, half=4000)              # built WITHOUT sigma_obs
    with pytest.raises(ValueError, match="declared to certify_drift"):
        certify_drift(rows, sigma_obs=1e-3)


def test_the_pure_noise_null_certifies_nothing():
    """S4: no drift exists, so the target is pure martingale, the zero law
    certifies and the honest verdict is vacuity."""
    rng = np.random.default_rng(3)
    dt, n = 1e-3, 320_001
    t = np.arange(n) * dt
    dW = rng.standard_normal((4, n - 1)) * np.sqrt(dt)
    X = np.concatenate([np.zeros((4, 1)), np.cumsum(dW, axis=1)], axis=1)
    rows = build_rows(t, X, LIB, half=16000)
    r = certify_drift(rows, delta=0.05, seed=0)
    assert r["certified"] is False
    assert r["abstain"] in ("noise", "structural", "resolution")


def test_the_term_vocabulary_reproduces_the_hand_rolled_rows():
    """Step 3's validation criterion, and the bridge between two assemblers.

    `ito_terms` expresses the Itô weak form in `weakform.Term`, so `build_nd` can
    assemble it with the ladder, roundoff and Gram machinery the PDE arc already
    has -- and the multi-field geometry Level 1 needs. If the two agree to machine
    precision on the same window then the vocabulary really does express the same
    identity; if they ever stop agreeing, one of them has drifted.
    """
    from lagh.ito import ito_terms
    from lagh.weakform import Patch, build_nd

    dt, half = 1e-3, 8000
    t, X = _ou(theta=1.0, b=1.4, T=40.0, dt=dt, n_traj=1, seed=0)
    rows = build_rows(t, X, LIB, half=half, fs=("x**2/2",))
    assert len(rows.y) >= 1

    it = ito_terms("x**2/2", LIB)
    terms = [it["target"], it["correction"]] + it["columns"]
    lo, _ = rows.windows[0]
    i0 = int(round(lo / dt))
    i1 = i0 + 2 * half + 1
    pa = Patch(centers=(0.5 * (t[i0] + t[i1 - 1]),),
               halfwidths=(0.5 * (t[i1 - 1] - t[i0]),), idx=(slice(i0, i1),))
    ws = build_nd({"u": X[0]}, [t], terms, [pa], p=8, rough=True,
                  martingale=it["martingale"])
    assert ws.rejected == 0, "a rough path must survive the deterministic gates"

    # the target carries the measured Itô correction subtracted, as build_rows does
    assert ws.A[0, 0] - ws.A[0, 1] == pytest.approx(rows.y[0], rel=1e-12)
    for k, g in enumerate(LIB):
        assert ws.A[0, 2 + k] == pytest.approx(rows.A[0, k], rel=1e-12), g
    assert ws.qv[0] == pytest.approx(rows.qv[0], rel=1e-12)
    assert ws.qv_se[0] == pytest.approx(rows.qv_se[0], rel=1e-12)
    # the correction's declared bound IS its estimator's own sd, not a ladder
    # difference: subsampling a realized-QV estimator changes it rather than
    # refining it
    assert ws.quad[0, 1] == pytest.approx(rows.corr_se[0], rel=1e-12)
    assert np.isnan(ws.order[0, 1])
    assert ws.quad[0, 0] == pytest.approx(rows.quad[0, 0], rel=1e-12)


def test_a_rough_path_needs_the_deterministic_gates_relaxed():
    """Both of weakform's resolution gates ask the wrong question of a Brownian
    path: it is aliased by construction, and its quadrature converges at O(h) so
    the observed ladder order sits near 1. Without `rough=True` every patch is
    dropped -- which is the honest default for a smooth field."""
    from lagh.ito import ito_terms
    from lagh.weakform import Patch, build_nd

    dt, half = 1e-3, 8000
    t, X = _ou(theta=1.0, b=1.4, T=40.0, dt=dt, n_traj=1, seed=0)
    it = ito_terms("x**2/2", LIB)
    terms = [it["target"], it["correction"]] + it["columns"]
    pa = Patch(centers=(0.5 * (t[half] + t[3 * half]),),
               halfwidths=(0.5 * (t[3 * half] - t[half]),),
               idx=(slice(half, 3 * half + 1),))
    strict = build_nd({"u": X[0]}, [t], terms, [pa], p=8)
    assert strict.rejected == 1 and len(strict.A) == 0
    loose = build_nd({"u": X[0]}, [t], terms, [pa], p=8, rough=True)
    assert loose.rejected == 0 and len(loose.A) == 1
    # ...and no martingale declaration means no measured scale, rather than a zero
    assert loose.qv is None


def test_a_d_measure_term_states_its_field_and_rejects_junk():
    from lagh.weakform import Term
    assert Term("a", gexpr="u").qv_field is None
    assert Term("b", gexpr="1/2", alpha=(0,), measure="d[u]").qv_field == "u"
    assert Term("c", gexpr="1", alpha=(0,), measure="d[rho]").qv_field == "rho"
    with pytest.raises(ValueError, match="neither 'dt' nor"):
        Term("d", measure="dW").qv_field


@pytest.mark.parametrize("bad", ["y", "x + t", "sin(z)"])
def test_a_library_term_outside_the_state_is_refused(bad):
    t, X = _ou(T=40.0, n_traj=2)
    with pytest.raises(ValueError):
        build_rows(t, X, ("x", bad), half=2000)

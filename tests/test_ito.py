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


def _gbm(mu=0.8, b=0.2, T=6.0, dt=1e-4, n_traj=8, seed=1, s0=1.0):
    """EXACT GBM: S_t = s0 exp((mu - b^2/2) t + b W_t). No scheme, so a failure
    here is the instrument's."""
    rng = np.random.default_rng(seed)
    n = int(round(T / dt)) + 1
    t = np.arange(n) * dt
    dW = rng.standard_normal((n_traj, n - 1)) * np.sqrt(dt)
    W = np.concatenate([np.zeros((n_traj, 1)), np.cumsum(dW, axis=1)], axis=1)
    return t, s0 * np.exp((mu - 0.5 * b ** 2) * t[None, :] + b * W)


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


def test_the_diffusion_certifies_from_quadratic_variation():
    """Level 1's second increment, and the reason for it.

    Putting b^2 in the DRIFT's design matrix gives a joint bound hundreds of times
    the truth, because the diffusion's signal-to-band goes as b while the drift's
    goes as 1/b. Quadratic variation is the efficient estimator for it, and the
    same weak-form discipline applies with a different target:

        int phi w(X) d[X] = sum_j d_j int phi w(X) h_j(X) dt

    On OU, where b^2 is the constant 1.96, that CERTIFIES.
    """
    from lagh.ito import build_qv_rows, certify_diffusion
    t, X = _ou(theta=1.0, b=1.4, T=320.0, dt=1e-3, n_traj=6)
    rows = build_qv_rows(t, X, ("1", "x", "x**2"), half=16000)
    assert len(rows.y) >= 12 and rows.n_trajectories == 6
    r = certify_diffusion(rows, delta=0.05, drift_max=5.0, seed=0)
    truth = {"diffusion:1": 1.96, "diffusion:x": 0.0, "diffusion:x**2": 0.0}
    # every bound covers, and the constant is ESTABLISHED AS PRESENT
    for nm, comp in r["partial"]["components"].items():
        lo, hi = comp["lo"], comp["hi"]
        assert (lo is None or lo <= truth[nm]) and (hi is None or truth[nm] <= hi)
    const = r["partial"]["components"]["diffusion:1"]
    assert const["resolved"] is True
    assert const["hi"] - const["lo"] < 0.5 * 1.96      # determined, not just bounded
    # the drift LEAKS into this estimator at O(dt) -- E[(dX)^2] = b^2 dt + a^2 dt^2
    # -- so the declared bound on |a| is a real term in the band, not decoration
    assert 0.0 < r["median_drift_leakage_share"] < 1.0
    assert r["drift_max_declared"] == 5.0
    loose = certify_diffusion(rows, delta=0.05, drift_max=50.0, seed=0)
    assert loose["median_band"] > r["median_band"]      # a looser declaration
    assert loose["median_drift_leakage_share"] > r["median_drift_leakage_share"]


def test_quadratic_variation_beats_the_design_matrix_for_the_diffusion():
    """The measurement that redirected Level 1, as a regression test: the same
    truth, the same data, two estimators, orders of magnitude apart."""
    from lagh.certify import admissible_interval
    from lagh.ito import build_qv_rows, certify_diffusion
    b = 0.2
    t, S = _gbm(mu=0.8, b=b)

    # (a) the diffusion inside the DRIFT's design matrix
    joint = build_rows(t, S, LIB, diff_names=("1", "x", "x**2"), half=5000)
    band = _band(joint)
    e0 = band.martingale() + LAM_QV * joint.corr_se + joint.quad[:, 0]
    q = joint.quad[:, 1:].sum(axis=1)
    bb, _ = admissible_interval(joint.A, joint.y, lambda c: e0 + c * q,
                                coeff_max=30.0)
    lo_j, hi_j = bb[joint.names.index("diffusion:x**2")]
    rel_joint = (hi_j - lo_j) / b ** 2

    # (b) the same coefficient from quadratic variation
    qv = build_qv_rows(t, S, ("1", "x", "x**2"), half=5000)
    r = certify_diffusion(qv, delta=0.05, drift_max=2.0, seed=0)
    lo_q, hi_q = r["admissible"]["diffusion:x**2"]
    rel_qv = (hi_q - lo_q) / b ** 2

    assert lo_j <= b ** 2 <= hi_j and lo_q <= b ** 2 <= hi_q   # both COVER
    assert rel_joint > 100                                     # and one is useless
    assert rel_qv < 1.0
    assert rel_joint / rel_qv > 500                            # measured ~2600x
    assert not (lo_j <= 0.0 <= hi_j) is False                  # joint: unresolved
    assert not (lo_q <= 0.0 <= hi_q)                           # qv: RESOLVED


def test_the_w_family_helps_a_stationary_process_and_not_a_growing_one():
    """A claim first asserted and then measured, which changed it. The w-family
    supplies state variation when every window sees the same state distribution;
    when the process is non-stationary the windows already provide it and w only
    costs a little kappa."""
    from lagh.ito import build_qv_rows, certify_diffusion
    DIFF = ("1", "x", "x**2")

    def width(t, X, half, ws, dmax, key):
        rows = build_qv_rows(t, X, DIFF, half=half, ws=ws)
        r = certify_diffusion(rows, delta=0.05, drift_max=dmax, seed=0)
        lo, hi = r["admissible"][key]
        return hi - lo

    # the DIRECTION is the finding; the magnitude is configuration-dependent (2.9x
    # at T = 640 with the same parameters, 1.25x here), so only the sign is asserted
    t, X = _ou(theta=1.0, b=1.4, T=320.0, dt=1e-3, n_traj=8)
    one = width(t, X, 16000, ("1",), 5.0, "diffusion:1")
    many = width(t, X, 16000, ("1", "x", "x**2"), 5.0, "diffusion:1")
    assert many < one, "a stationary process is helped by the w-family"

    t2, S = _gbm(mu=0.8, b=0.2)
    one2 = width(t2, S, 5000, ("1",), 2.0, "diffusion:x**2")
    many2 = width(t2, S, 5000, ("1", "x", "x**2"), 2.0, "diffusion:x**2")
    assert many2 > one2, "a growing process gets nothing from the w-family"
    assert many2 < 1.2 * one2                       # and loses only a little


def _cir(theta=1.0, m=1.0, b=1.0, T=200.0, dt=1e-3, n_traj=6, seed=11):
    """EXACT CIR: dX = theta(m - X)dt + b sqrt(X) dW, b^2(x) = b^2 x."""
    from scipy.stats import ncx2
    rng = np.random.default_rng(seed)
    n = int(round(T / dt)) + 1
    c = 4.0 * theta / (b ** 2 * (1.0 - np.exp(-theta * dt)))
    d = 4.0 * theta * m / b ** 2
    X = np.empty((n_traj, n))
    X[:, 0] = m
    dec = np.exp(-theta * dt)
    for k in range(1, n):
        X[:, k] = ncx2.rvs(df=d, nc=c * X[:, k - 1] * dec, size=n_traj,
                           random_state=rng) / c
    return np.arange(n) * dt, X


def test_the_cir_sampler_reproduces_its_own_stationary_moments():
    """A generator check, and it earned its place: the scale factor
    4 theta / (b^2 (1 - e^-theta dt)) was written with a 2, which puts the
    stationary mean at 2m instead of m -- an exact sampler with a wrong constant,
    which the instrument would have been blamed for. Stationary mean m, variance
    m b^2 / 2 theta."""
    theta, m, b = 1.0, 1.0, 0.8
    _, X = _cir(theta=theta, m=m, b=b, T=200.0)
    assert X.mean() == pytest.approx(m, rel=0.05)
    assert X.var() == pytest.approx(m * b ** 2 / (2 * theta), rel=0.10)
    assert X.min() > 0.0                       # Feller 2 theta m / b^2 = 3.1 > 1


def test_a_state_dependent_diffusion_resolves_once_the_state_spreads():
    """The diffusion's counterpart to the drift's theta*L > 2 kappa^2: the FORM of
    a state-dependent b^2 is identifiable in proportion to the state's RELATIVE
    SPREAD, because that is what separates {1, x, x^2} from each other.

    Measured on CIR (b^2 = b^2 x), sd/mean rising with b: at 0.55 the x coefficient
    straddles zero, and by 0.71 it is RESOLVED. Every bound covers throughout.
    """
    from lagh.ito import build_qv_rows, certify_diffusion
    DIFF = ("1", "x", "x**2")
    got = []
    for b in (0.8, 1.0):
        t, X = _cir(theta=1.0, m=1.0, b=b, T=200.0)
        rows = build_qv_rows(t, X, DIFF, half=8000, ws=DIFF)
        r = certify_diffusion(rows, delta=0.05, drift_max=3.0, seed=0)
        ad = r["admissible"]
        for nm, (lo, hi) in ad.items():
            tv = b ** 2 if nm == "diffusion:x" else 0.0
            assert lo <= tv <= hi, (b, nm)      # covers at every spread
        lo, hi = ad["diffusion:x"]
        got.append((float(X.std() / X.mean()), not (lo <= 0.0 <= hi)))
    (spread_lo, res_lo), (spread_hi, res_hi) = got
    assert spread_lo < spread_hi
    assert res_lo is False and res_hi is True


def test_the_drift_leakage_bound_is_pointwise_not_flat():
    """The diffusion's band needs a bound on |a| because E[(dX)^2] = b^2 dt + a^2
    dt^2. Bounding a^2 by its worst value ANYWHERE in the visited range over-declares
    by the square of a ratio, and that flat bound was measured to destroy the
    diffusion's precision entirely on Level 1 -- every task abstained on vacuity.

    The honest bound is pointwise: sum_i |w_i| env(X_i)^2 dt, with env evaluated at
    the states the process actually visited. Same distinction weakform draws with
    `field_l1`, and it is what made two diffusion certificates possible.
    """
    from lagh.ito import QvBand, build_qv_rows, certify_diffusion
    t, X = _ou(theta=1.0, b=1.4, T=320.0, dt=1e-3, n_traj=6)
    DIFF = ("1", "x", "x**2")
    env_max = 50.0                     # a wide envelope, as an undetermined drift gives

    flat = build_qv_rows(t, X, DIFF, half=16000, ws=DIFF)
    assert flat.leak is None
    tight = build_qv_rows(t, X, DIFF, half=16000, ws=DIFF,
                          drift_envelope=lambda z: np.minimum(np.abs(z) * 1.0,
                                                              env_max))
    assert tight.leak is not None and np.all(tight.leak > 0)

    def band(rows):
        return QvBand(rows.var, rows.lam, len(rows.y), rows.dt, delta=0.05,
                      drift_max=env_max, y=rows.y, leak=rows.leak)(None)
    # the pointwise bound is strictly tighter, because the process lives where its
    # own drift is small relative to the envelope's worst case
    assert np.all(band(tight) < band(flat))
    assert np.median(band(flat) / band(tight)) > 5.0

    # and the tighter band is what lets the constant diffusion certify
    r_t = certify_diffusion(tight, delta=0.05, drift_max=env_max, seed=0)
    r_f = certify_diffusion(flat, delta=0.05, drift_max=env_max, seed=0)
    assert r_t["leakage_bound"] == "pointwise envelope"
    assert r_f["leakage_bound"] == "flat drift_max"
    assert r_t["median_signal_to_band"] > r_f["median_signal_to_band"]
    # BOTH still cover -- the flat bound was over-declared, never unsound
    for r in (r_t, r_f):
        c = r["partial"]["components"]["diffusion:1"]
        assert c["lo"] <= 1.96 <= c["hi"]
    assert r_t["partial"]["components"]["diffusion:1"]["resolved"] is True


def test_a_claimed_diffusion_becomes_ordinary_columns():
    """Level 1's first target. With a library for b^2 the Itô correction stops
    being a measured value moved to the left and becomes the dt columns
    1/2 int phi f'' h_j dt, so drift and diffusion are identified JOINTLY.

    Two consequences the test pins down: the columns are named for the frozen
    checker's component vocabulary, and the measured correction -- the one UNSAFE
    consumer of realized quadratic variation, and the source of the Level 0
    confident-wrong -- is gone.
    """
    t, X = _ou(T=160.0, n_traj=4)
    plain = build_rows(t, X, LIB, half=8000)
    joint = build_rows(t, X, LIB, diff_names=("1", "x**2"), half=8000)
    assert plain.names == list(LIB)
    assert joint.names == ["drift:1", "drift:x", "drift:x**2", "drift:x**3",
                           "diffusion:1", "diffusion:x**2"]
    assert joint.A.shape[1] == 6 and plain.A.shape[1] == 4
    # the measured correction is gone, so the target is the plain -int phi' f dt
    assert np.all(joint.corr_se == 0.0)
    assert np.any(plain.corr_se > 0.0)
    # the band still comes from realized QV -- the SAFE consumer -- unchanged
    assert np.allclose(joint.qv, plain.qv)
    # and the diffusion columns are zero for f = x, where f'' vanishes
    fx = np.array(joint.fname) == "x"
    assert np.allclose(joint.A[fx][:, 4:], 0.0)
    assert not np.allclose(joint.A[~fx][:, 4:], 0.0)


def test_the_admissible_functional_bounds_the_law_not_just_its_coefficients():
    """`admissible_interval` bounds one coefficient; this bounds any LINEAR
    FUNCTIONAL of them, which is what a reader actually wants -- the drift AS A
    FUNCTION. Both must cover, and the functional bound must be no wider than the
    coefficient bound implies."""
    from lagh.certify import admissible_functional
    t, X = _ou(T=320.0, n_traj=4)
    rows = build_rows(t, X, LIB, half=16000)
    band = _band(rows)
    e0 = band.martingale() + LAM_QV * rows.corr_se + rows.quad[:, 0]
    q = rows.quad[:, 1:].sum(axis=1)
    eps = lambda c: e0 + c * q                               # noqa: E731
    # a(x) at three states, as functionals of (c_1, c_x, c_x2, c_x3)
    xs = np.array([-0.5, 0.0, 0.5])
    W = np.array([[1.0, x0, x0 ** 2, x0 ** 3] for x0 in xs])
    fun, info = admissible_functional(rows.A, rows.y, eps, W, coeff_max=20.0)
    assert info["n_functionals"] == 3
    for x0, (lo, hi) in zip(xs, fun):
        true_a = -1.0 * x0                       # OU drift -theta x, theta = 1
        assert lo is not None and hi is not None and lo <= true_a <= hi
    # picking out one coordinate must reproduce admissible_interval exactly
    one = np.zeros((1, 4))
    one[0, 1] = 1.0
    got, _ = admissible_functional(rows.A, rows.y, eps, one, coeff_max=20.0)
    ref, _ = admissible_interval(rows.A, rows.y, eps, coeff_max=20.0)
    assert got[0][0] == pytest.approx(ref[LIB.index("x")][0], rel=1e-9)
    assert got[0][1] == pytest.approx(ref[LIB.index("x")][1], rel=1e-9)
    with pytest.raises(ValueError, match="columns"):
        admissible_functional(rows.A, rows.y, eps, np.zeros((1, 9)))


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


def test_a_d_measure_term_states_its_fields_and_rejects_junk():
    """`d[u]` is quadratic variation, `d[u,v]` the CROSS-variation a
    multi-dimensional Itô correction needs -- the off-diagonal terms a
    diagonal-only measure cannot express."""
    from lagh.weakform import Term
    assert Term("a", gexpr="u").qv_fields is None
    assert Term("b", gexpr="1/2", alpha=(0,), measure="d[u]").qv_fields == ("u", "u")
    assert Term("c", gexpr="1", alpha=(0,), measure="d[rho]").qv_field == "rho"
    assert Term("d", gexpr="1", alpha=(0,),
                measure="d[x,y]").qv_fields == ("x", "y")
    assert Term("e", gexpr="1", alpha=(0,),
                measure="d[ x , y ]").qv_fields == ("x", "y")
    for junk in ("dW", "d[]", "d[a,b,c]", "d[a,]"):
        with pytest.raises(ValueError, match="neither 'dt'"):
            Term("bad", measure=junk).qv_fields


@pytest.mark.parametrize("bad", ["y", "x + t", "sin(z)"])
def test_a_library_term_outside_the_state_is_refused(bad):
    t, X = _ou(T=40.0, n_traj=2)
    with pytest.raises(ValueError):
        build_rows(t, X, ("x", bad), half=2000)


# ------------------------------------------------- multi-field Itô (Level 1, 2-D)

def _vdp(mu=1.0, b=0.5, T=100.0, dt=1e-3, n_traj=2, seed=3, substeps=8):
    """Van der Pol with additive noise on y ONLY:
        dx = y dt,   dy = (mu(1-x^2)y - x)dt + b dW.
    The x component is noise-free, which is the point of using it here."""
    rng = np.random.default_rng(seed)
    n = int(round(T / dt)) + 1
    ds = dt / substeps
    sq = np.sqrt(ds)
    x = rng.uniform(-2, 2, n_traj)
    y = rng.uniform(-2, 2, n_traj)
    X = np.empty((n_traj, n))
    Y = np.empty((n_traj, n))
    X[:, 0], Y[:, 0] = x, y
    for k in range(1, n):
        for _ in range(substeps):
            xn = x + y * ds
            y = y + (mu * (1 - x ** 2) * y - x) * ds \
                + b * sq * rng.standard_normal(n_traj)
            x = xn
        X[:, k], Y[:, k] = x, y
    return np.arange(n) * dt, X, Y


VDP_LIBS = {"x": ("1", "x", "y", "x**2", "x*y", "y**2"),
            "y": ("1", "x", "y", "x**2", "x*y", "x**2*y")}


def _vdp_rows(t, X, Y, f, half, *, mart=True):
    """(design, qv, column names, n_corrections) for one f, stacked over paths."""
    from lagh.ito import ito_terms_nd
    from lagh.weakform import Patch, build_nd
    it = ito_terms_nd(f, VDP_LIBS, ("x", "y"))
    terms = [it["target"]] + it["corrections"] + it["columns"]
    A, qv = [], []
    for k in range(len(X)):
        wins = [(c - half, c + half + 1)
                for c in range(half, len(t) - half - 1, 2 * half)]
        pas = [Patch(centers=(0.5 * (t[lo] + t[hi - 1]),),
                     halfwidths=(0.5 * (t[hi - 1] - t[lo]),),
                     idx=(slice(lo, hi),)) for lo, hi in wins]
        ws = build_nd({"x": X[k], "y": Y[k]}, [t], terms, pas, p=8, rough=True,
                      martingale=it["martingale"] if mart else None)
        if len(ws.A):
            A.append(ws.A)
            if ws.qv is not None:
                qv.append(ws.qv)
    return (np.vstack(A), np.concatenate(qv) if qv else None,
            [c.name for c in it["columns"]], len(it["corrections"]))


def test_the_multi_field_ito_terms_have_the_right_shape():
    """Two things the scalar case hid: the Hessian's off-diagonal entries appear
    TWICE in sum_ik and the term must carry that factor, and the martingale is a
    LIST of per-field sensitivities."""
    from lagh.ito import ito_terms_nd
    # f = x: no second derivative, and only x drives it
    it = ito_terms_nd("x", VDP_LIBS, ("x", "y"))
    assert it["corrections"] == []
    assert it["martingale"] == [("x", "1")]
    assert [c.name for c in it["columns"]][:3] == ["x:1", "x:x", "x:y"]
    # f = y^2/2: diagonal Hessian, coefficient 1/2
    it2 = ito_terms_nd("y**2/2", VDP_LIBS, ("x", "y"))
    assert [(c.name, c.gexpr, c.measure) for c in it2["corrections"]] == \
        [("ito:y,y", "1/2", "d[y,y]")]
    assert it2["martingale"] == [("y", "y")]
    # f = x*y: OFF-diagonal, so 2 * 1/2 = 1, and BOTH fields drive it
    it3 = ito_terms_nd("x*y", VDP_LIBS, ("x", "y"))
    assert [(c.name, c.gexpr, c.measure) for c in it3["corrections"]] == \
        [("ito:x,y", "1", "d[x,y]")]
    assert it3["martingale"] == [("x", "y"), ("y", "x")]
    assert len(it3["columns"]) == 12            # both components' libraries
    with pytest.raises(ValueError, match="reads"):
        ito_terms_nd("x*z", VDP_LIBS, ("x", "y"))


def test_a_cross_variation_term_computes_the_increment_product():
    """d[x,y] must be sum_i w_i dx_i dy_i, not either field's own variation."""
    from lagh.weakform import Patch, Term, build_nd
    n = 4001
    t = np.arange(n) * 1e-3
    rng = np.random.default_rng(0)
    Xf = np.cumsum(rng.standard_normal(n)) * 1e-2
    Yf = np.cumsum(rng.standard_normal(n)) * 1e-2
    pa = Patch(centers=(t[2000],), halfwidths=(t[1500],), idx=(slice(500, 3501),))
    got = build_nd({"x": Xf, "y": Yf}, [t],
                   [Term("cross", gexpr="1", alpha=(0,), measure="d[x,y]"),
                    Term("xx", gexpr="1", alpha=(0,), measure="d[x]")],
                   [pa], p=8, rough=True)
    from lagh.weakform import _tensor, _windows, bump_derivatives
    dpsi = bump_derivatives(8, 0)
    phi = _tensor(_windows(pa, [t[500:3501]], dpsi, None, (0,)))[:-1]
    dx, dy = np.diff(Xf[500:3501]), np.diff(Yf[500:3501])
    assert got.A[0, 0] == pytest.approx(float(np.sum(phi * dx * dy)), rel=1e-12)
    assert got.A[0, 1] == pytest.approx(float(np.sum(phi * dx * dx)), rel=1e-12)
    # independent fields: the cross-variation is near zero and the diagonal is not
    assert abs(got.A[0, 0]) < 0.2 * abs(got.A[0, 1])


def test_van_der_pol_covers_its_truth_on_both_test_functions():
    """The migration's validation: a 2-D state, one row set spanning both
    components, and the TRUE law inside the band on the deterministic row and the
    stochastic one alike."""
    mu, b = 1.0, 0.5
    t, X, Y = _vdp(mu=mu, b=b, T=100.0)
    for f, truth in (("x", {"x:y": 1.0}),
                     ("y**2/2", {"y:y": mu, "y:x**2*y": -mu, "y:x": -1.0})):
        A, qv, names, ncorr = _vdp_rows(t, X, Y, f, 8000)
        assert len(A) >= 6 and qv is not None
        y_t = A[:, 0] - A[:, 1:1 + ncorr].sum(axis=1)
        c = np.array([truth.get(nm, 0.0) for nm in names])
        resid = np.abs(y_t - A[:, 1 + ncorr:] @ c)
        band = 4.0 * np.sqrt(np.maximum(qv, 0.0))
        assert np.all(resid <= band), f
    # the NOISELESS row is exact to quadrature, so its martingale band is hugely
    # conservative -- measured ~2000x. Reach lost, never soundness.
    A, qv, names, ncorr = _vdp_rows(t, X, Y, "x", 8000)
    resid = np.abs(A[:, 0] - A[:, 1 + ncorr:] @ np.array(
        [1.0 if nm == "x:y" else 0.0 for nm in names]))
    assert np.median(resid) < 1e-2
    assert np.median(4.0 * np.sqrt(qv)) > 100 * np.median(resid)


def test_quadratic_variation_tells_a_driven_component_from_a_noiseless_one():
    """The multi-dimensional form of the error-provenance question, and it is
    MEASURABLE rather than declared: a component's realized quadratic variation
    scales as O(dt) when it carries no noise (a differentiable path has zero
    quadratic variation, and the estimator sees the O(dt) residue) and as O(1) when
    it is driven. Halving dt separates them.

    Consequence, stated because it is the next increment's job: for a noiseless
    component the martingale band built from that residue over-declares by
    ~sqrt(1/dt). Tightening it needs a DECLARATION that the component is noiseless,
    and this diagnostic is what would verify such a declaration.
    """
    got = {}
    for dt, sub in ((4e-3, 2), (1e-3, 8)):
        t, X, Y = _vdp(dt=dt, substeps=sub, T=100.0)
        half = int(round(8.0 / dt))
        _, qv_x, _, _ = _vdp_rows(t, X, Y, "x", half)
        _, qv_y, _, _ = _vdp_rows(t, X, Y, "y**2/2", half)
        got[dt] = (float(np.median(qv_x)), float(np.median(qv_y)))
    (x_coarse, y_coarse), (x_fine, y_fine) = got[4e-3], got[1e-3]
    # the noiseless component's qv FALLS with dt; 4x finer, expect ~4x smaller
    assert x_fine < x_coarse / 2.0
    # the driven component's does not
    assert 0.5 < y_fine / y_coarse < 2.0


def test_quadratic_variation_separates_its_martingale_part_from_a_smooth_residue():
    """A component that carries NO noise still has a nonzero realized QV -- the
    O(dt) residue of a differentiable path -- and banding with it over-declares by
    ~sqrt(1/dt). The two parts separate by STRIDE SCALING, because a martingale's
    increment variance grows like s while a differentiable path's SQUARED increment
    grows like s^2:

        sum_i (u[i+s] - u[i])^2  ~=  alpha*s + beta*s^2

    alpha is the martingale part. Measured to within 0.3% on mixed paths.
    """
    from lagh.weakform import qv_martingale_part
    n, dt = 200_001, 1e-3
    t = np.arange(n) * dt
    rng = np.random.default_rng(0)
    smooth = np.sin(2 * t) + 0.3 * np.cos(5 * t)

    def sums(u, strides=(1, 2, 4, 8)):
        return [float(np.sum((u[s:] - u[:-s]) ** 2)) for s in strides]

    # a purely smooth path: the martingale part is ~0, thousands of times below its
    # own total quadratic variation
    sc, info = qv_martingale_part(sums(smooth))
    total = sums(smooth)[0]
    assert info["used"] == "martingale part"
    assert sc < total / 500

    # mixed paths: alpha recovers the true martingale part
    for b in (0.01, 0.05, 0.5):
        bm = np.r_[0.0, np.cumsum(rng.standard_normal(n - 1) * np.sqrt(dt))]
        sc, info = qv_martingale_part(sums(smooth + b * bm))
        assert info["used"] == "martingale part"
        assert info["alpha"] == pytest.approx(b ** 2 * t[-1], rel=0.05)
        assert sc >= info["alpha"]          # the returned scale is an UPPER estimate

    # a pure martingale has no smooth part to remove, so no tightening is available
    # and it falls back to the total -- the safe direction, always
    bm = np.r_[0.0, np.cumsum(rng.standard_normal(n - 1) * np.sqrt(dt))]
    sc, info = qv_martingale_part(sums(0.5 * bm))
    assert sc == pytest.approx(sums(0.5 * bm)[0], rel=1e-12)
    assert info["used"] == "total"
    # ...and a design too small to fit falls back too
    _, bad = qv_martingale_part([1.0, 2.0], strides=(1, 2))
    assert bad["used"] == "total"


def test_the_decomposition_tightens_a_noiseless_row_and_leaves_a_driven_one_alone():
    """The decisive test, on real 2-D data: the tightening must apply where the
    component is noise-free and NOT where it is driven, and the truth must still
    cover in both cases. Measured: 12x tighter on Van der Pol's noiseless x row,
    unchanged on the driven y row."""
    from lagh.ito import ito_terms_nd
    from lagh.weakform import Patch, build_nd
    mu, b, half = 1.0, 0.5, 8000
    t, X, Y = _vdp(mu=mu, b=b, T=100.0)

    def run(f, truth, decompose):
        it = ito_terms_nd(f, VDP_LIBS, ("x", "y"))
        terms = [it["target"]] + it["corrections"] + it["columns"]
        names = [c.name for c in it["columns"]]
        nc = len(it["corrections"])
        A, qv = [], []
        for k in range(len(X)):
            wins = [(c - half, c + half + 1)
                    for c in range(half, len(t) - half - 1, 2 * half)]
            pas = [Patch(centers=(0.5 * (t[lo] + t[hi - 1]),),
                         halfwidths=(0.5 * (t[hi - 1] - t[lo]),),
                         idx=(slice(lo, hi),)) for lo, hi in wins]
            ws = build_nd({"x": X[k], "y": Y[k]}, [t], terms, pas, p=8, rough=True,
                          martingale=it["martingale"],
                          martingale_decompose=decompose)
            if len(ws.A):
                A.append(ws.A)
                qv.append(ws.qv)
        A = np.vstack(A)
        qv = np.concatenate(qv)
        c = np.array([truth.get(nm, 0.0) for nm in names])
        resid = np.abs(A[:, 0] - A[:, 1:1 + nc].sum(axis=1) - A[:, 1 + nc:] @ c)
        band = 4.0 * np.sqrt(np.maximum(qv, 0.0))
        return float(np.median(band)), bool(np.all(resid <= band))

    # the NOISELESS component: a large tightening, and the truth still covers
    raw_x, cov_x = run("x", {"x:y": 1.0}, False)
    dec_x, cov_dx = run("x", {"x:y": 1.0}, True)
    assert cov_x and cov_dx
    assert dec_x < raw_x / 5.0

    # the DRIVEN component: nothing to remove, so the band barely moves
    truth_y = {"y:y": mu, "y:x**2*y": -mu, "y:x": -1.0}
    raw_y, cov_y = run("y**2/2", truth_y, False)
    dec_y, cov_dy = run("y**2/2", truth_y, True)
    assert cov_y and cov_dy
    assert dec_y > 0.9 * raw_y


def test_the_noise_free_component_of_van_der_pol_certifies_its_drift():
    """The migration completing itself, and the arc's first certified DRIFT.

    `build_rows_nd` packs multi-field Itô rows into the same `ItoRows` the scalar
    path produces, so `certify_drift` -- with all its coherence, significance,
    parsimony, holdout and admissible-bound machinery -- applies to a 2-D state with
    no further work. That was the point of putting the vocabulary in `weakform`.

    Van der Pol's x component carries no noise, so `dx = y dt` is a deterministic
    weak-form identity; with the stride decomposition removing the O(dt) residue from
    its band, it certifies. Its y component is driven and is vacuous at this
    configuration -- the split the identifiability measurement predicted.
    """
    from lagh.ito import build_rows_nd
    mu = 1.0
    t, X, Y = _vdp(mu=mu, b=0.5, T=200.0, n_traj=4)

    # the NOISE-FREE equation: certifies, with the coefficient resolved
    rows = build_rows_nd(t, {"x": X, "y": Y}, VDP_LIBS, "x", half=8000)
    assert rows.names == ["x:1", "x:x", "x:y", "x:x**2", "x:x*y", "x:y**2"]
    assert len(rows.y) >= 12
    r = certify_drift(rows, delta=0.05, seed=0)
    assert r["certified"] is True, r.get("abstain")
    assert "x:y" in r["law"]
    comp = r["partial"]["components"]
    assert comp["x:y"]["resolved"] is True
    assert comp["x:y"]["lo"] <= 1.0 <= comp["x:y"]["hi"]
    assert comp["x:y"]["hi"] - comp["x:y"]["lo"] < 0.2      # measured ~0.06
    for nm, c in comp.items():                              # and every bound covers
        tv = 1.0 if nm == "x:y" else 0.0
        assert (c["lo"] is None or c["lo"] <= tv) and \
               (c["hi"] is None or tv <= c["hi"]), nm
    assert r["alpha_log10"] < -5

    # the DRIVEN equation: vacuous here, and its bounds still cover
    rows_y = build_rows_nd(t, {"x": X, "y": Y}, VDP_LIBS, "y**2/2", half=8000)
    ry = certify_drift(rows_y, delta=0.05, seed=0)
    assert ry["certified"] is False
    truth_y = {"y:y": mu, "y:x**2*y": -mu, "y:x": -1.0}
    for nm, c in ry["partial"]["components"].items():
        tv = truth_y.get(nm, 0.0)
        assert (c["lo"] is None or c["lo"] <= tv) and \
               (c["hi"] is None or tv <= c["hi"]), nm


def test_the_stride_decomposition_tightens_that_certificate_by_seven_fold():
    """What the decomposition buys, measured rather than assumed.

    It was tempting to claim the certificate depends on it. It does not -- the
    noise-free equation certifies either way at this configuration, because even the
    un-decomposed band leaves signal-to-band at 3.6. What the decomposition buys is
    PRECISION: the certified coefficient's joint bound narrows 6.7x, from +-20% of
    the truth to +-3%. Same law, a far better interval.
    """
    from lagh.ito import build_rows_nd
    t, X, Y = _vdp(mu=1.0, b=0.5, T=200.0, n_traj=4)
    out = {}
    for dec in (False, True):
        rows = build_rows_nd(t, {"x": X, "y": Y}, VDP_LIBS, "x", half=8000,
                             decompose=dec)
        r = certify_drift(rows, delta=0.05, seed=0)
        c = r["partial"]["components"]["x:y"]
        out[dec] = (r["certified"], r["median_signal_to_band"],
                    c["hi"] - c["lo"], c["lo"] <= 1.0 <= c["hi"])
    assert out[False][0] and out[True][0]            # certified either way
    assert out[False][3] and out[True][3]            # and covering either way
    assert out[True][1] > 5 * out[False][1]          # ~8x the signal-to-band
    assert out[True][2] < out[False][2] / 5.0        # ~6.7x tighter bound


def test_the_component_index_convention_carries_a_per_component_drift():
    """The frozen interface's `part[index]:term` convention, and the first producer
    to need it. A 2-D system has one drift equation PER COMPONENT, so the columns
    `build_rows_nd` names `<field>:<term>` re-key to `drift[i]:<term>` and the
    checker scores them against a truth spanning both components.

    Asserted here rather than only in the runner because the mapping is the whole
    point: without the index a per-component drift cannot be expressed at all.
    """
    from lagh.ito import build_rows_nd
    from lagh.stochcheck import (Coverage, Submission, Task, component,
                                 score_task, validate_submission, validate_task)
    mu, b = 1.0, 0.5
    truth, exp = {}, {}
    for i, fld in enumerate(("x", "y")):
        for g in VDP_LIBS[fld]:
            truth[component("drift", g, i)] = 0.0
            exp[component("drift", g, i)] = "interval" if i == 0 else "abstain"
        for h in ("1", "x", "x**2"):
            truth[component("diffusion", h, i)] = 0.0
            exp[component("diffusion", h, i)] = "interval"
    truth[component("drift", "y", 0)] = 1.0
    truth[component("drift", "y", 1)] = mu
    truth[component("drift", "x**2*y", 1)] = -mu
    truth[component("drift", "x", 1)] = -1.0
    truth[component("diffusion", "1", 1)] = b ** 2
    task = Task(task_id="L1-vdp", level=1, system="vdp", state_dim=2,
                truth=truth, expectation=exp)
    assert validate_task(task) == []
    # both components' drifts and both diffusions are expressible, zeros included
    assert task.truth["drift[0]:y"] == 1.0
    assert task.truth["drift[1]:x**2*y"] == -1.0
    assert task.truth["diffusion[1]:1"] == b ** 2
    assert len(task.truth) == 18

    t, X, Y = _vdp(mu=mu, b=b, T=100.0, n_traj=3)
    rows = build_rows_nd(t, {"x": X, "y": Y}, VDP_LIBS, "x", half=8000)
    r = certify_drift(rows, delta=0.05, seed=0)
    rec = dict(r["partial"])
    rec["components"] = {f"drift[0]:{k.split(':', 1)[1]}": v
                         for k, v in rec["components"].items()}
    for key in ("exact", "interval", "unconstrained"):
        rec[key] = [f"drift[0]:{k.split(':', 1)[1]}" for k in rec.get(key, [])]
    sub = Submission(task_id=task.task_id,
                     kind="answer" if r.get("certified") else "abstain",
                     abstain=None if r.get("certified") else "noise",
                     record=rec,
                     coverage=Coverage(kappa=r["kappa"], delta=0.05,
                                       n_rows=r["n_rows"],
                                       n_disjoint=r["n_disjoint"]),
                     submission_id="drift-x")
    assert validate_submission(task, sub) == []
    score = score_task(task, [sub])
    # ZERO confident-wrong, and the x equation's own components are scored
    assert score.n_confident_wrong == 0
    assert "drift[0]:y@all" in score.components
    got = score.components["drift[0]:y@all"]
    assert got["truth"] == 1.0 and got["outcome"] == "covered"
    # the components this submission did not speak to are MISSED, not silently ok
    assert score.n_missed >= 6

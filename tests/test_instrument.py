"""What a real instrument does to a record, declared and measured
(lagh/instrument.py; docs/CASE_STUDY_TWEEZERS_C1.md §7). Every test simulates
the C-Trap bead of C1 -- theta 3297 /s, b^2 2.007e5 nm^2/s, 78125 Hz -- and
applies the distortions ourselves, so the truth is known."""

import numpy as np
import pytest

from lagh.instrument import (COMMON_MODE_SPREAD, acf_timescale, axis_gate,
                             band_loss, diode_response,
                             drive_scale, equipartition_ratio, faxen_drag_ratio,
                             faxen_height, gate_record, local_drag, lockin,
                             attribute_deviation, realized_diffusion, retention)
from lagh.ito import build_qv_rows, certify_diffusion

THETA, B2, FS = 3297.0, 2.007e5, 78125.0
DT = 1.0 / FS
ALPHA, F_DIODE, CUTOFF = 0.42, 8400.0, 23000.0


OVER = 8            # the instrument filters the CONTINUOUS signal, then samples


def ou(n=400_000, ntraj=1, seed=3, theta=THETA, b2=B2, over=1):
    """Discrete OU at `over` x the sampling rate (exact transition density).

    The recursion x[k] = e^{-theta dt} x[k-1] + sd z[k] IS a one-pole filter, so
    `lfilter` runs it in C: a 6-million-step Python loop was the test suite's
    bottleneck and kept the simulated records shorter than the real ones."""
    from scipy.signal import lfilter
    rng = np.random.default_rng(seed)
    dt = DT / over
    dec = np.exp(-theta * dt)
    sd = np.sqrt(b2 * (1 - dec ** 2) / (2 * theta))
    m = n * over
    x0 = rng.standard_normal((ntraj, 1)) * np.sqrt(b2 / (2 * theta))
    z = rng.standard_normal((ntraj, m - 1))
    rest = lfilter([sd], [1.0, -dec], z, axis=-1, zi=dec * x0)[0]
    return np.hstack([x0, rest])


def diode(X, alpha=ALPHA, f_diode=F_DIODE, over=OVER):
    """The instrument's own detector model, applied forwards on the FINE grid:
    alpha of the position instantly, the rest through a first-order low-pass."""
    from scipy.signal import lfilter
    a = np.exp(-2 * np.pi * f_diode * DT / over)
    y = lfilter([1 - a], [1.0, -a], X, axis=-1)
    return alpha * X + (1 - alpha) * y


def anti_alias(X, cutoff=CUTOFF, over=OVER):
    from scipy.signal import butter, filtfilt
    b, a = butter(4, cutoff / (0.5 * FS * over))
    return filtfilt(b, a, X, axis=-1)


def instrument(X_fine, over=OVER):
    """Detector, anti-alias filter, then SAMPLING -- in the instrument's order.
    Filtering the sampled series instead leaves the aliased above-Nyquist part
    of the process unfiltered, a measured 13% error in the loss."""
    return anti_alias(diode(X_fine, over=over), over=over)[..., ::over]


def test_the_axis_gate_passes_a_clean_record_and_refuses_a_contaminated_one():
    x = ou()[0]
    g = axis_gate(x, DT, THETA)
    assert g["passed"] and abs(g["ratio"] - 1) < 0.1
    # a slow contaminant below any PSD fit floor: the timescale collapses
    t = np.arange(len(x)) * DT
    bad = x + 3 * x.std() * np.sin(2 * np.pi * 7.0 * t)
    g = axis_gate(bad, DT, THETA)
    assert not g["passed"] and g["ratio"] < 0.5
    # and the gate is SCALE-FREE: a wrong displacement scale changes nothing
    assert axis_gate(1e3 * x, DT, THETA)["ratio"] == pytest.approx(g["ratio"] * 0 + axis_gate(x, DT, THETA)["ratio"])
    # a detector filter does not fool it either (it reads lags beyond the filter)
    assert axis_gate(instrument(ou(n=100_000, over=OVER))[0], DT, THETA)["passed"]


def test_acf_timescale_refuses_a_record_it_cannot_read():
    assert acf_timescale(np.zeros(1000), DT)["theta"] is None
    assert acf_timescale(np.random.default_rng(0).standard_normal(5000), DT)["theta"] is None


def test_the_band_loss_is_the_measured_qv_attenuation_and_its_factors_are_named():
    """Apply the diode (declared) and an anti-alias filter (measured) to a
    clean OU record; the declared band loss must equal the ratio of recorded to
    true quadratic variation, and its Nyquist factor must be the process's
    above-band share."""
    Xf = ou(n=300_000, over=OVER)
    X = Xf[..., ::OVER]
    Xr = instrument(Xf)
    qv_true = float(np.mean(np.diff(X, axis=-1) ** 2)) / DT          # ~ b^2
    qv_rec = float(np.mean(np.diff(Xr, axis=-1) ** 2)) / DT
    bl = band_loss(Xr[0], FS, fit_range=(100.0, CUTOFF), alpha=ALPHA,
                   f_diode=F_DIODE)
    # the loss is the ratio of recorded to SAMPLED-PROCESS quadratic variation
    # (the model's denominator is the exact discrete increment variance, which
    # sits at (1 - exp(-theta dt))/(theta dt) = 0.98 of b^2 dt for this bead);
    # what remains is the fine grid's own aliasing, ~3%
    assert abs(bl["r_total"] / (qv_rec / qv_true) - 1) < 0.05
    assert abs(qv_rec / bl["r_total"] / qv_true - 1) < 0.05         # corrected: truth
    assert abs(qv_true / B2 - 1) < 0.06                              # theta dt = 4%
    assert 0.7 < bl["r_nyquist"] < 0.85                              # 23% above Nyquist
    assert 0.35 < bl["r_diode"] < 0.55                               # alpha^2 plateau
    assert 0.75 < bl["r_antialias"] < 1.0        # filtfilt: 8th order at the ceiling
    assert bl["r_se"] < 0.02
    assert abs(bl["fc_fit_hz"] - THETA / (2 * np.pi)) / (THETA / (2 * np.pi)) < 0.2
    # the true record has no loss beyond Nyquist -- but the diode is DECLARED,
    # so a clean record analysed with a declared diode still reports its factor
    # (the declaration is the instrument's, not the data's)
    bl0 = band_loss(X[0], FS, fit_range=(100.0, CUTOFF), alpha=1.0, f_diode=F_DIODE)
    assert abs(bl0["r_diode"] - 1) < 1e-9 and bl0["r_antialias"] > 0.97


def test_the_declared_band_loss_closes_the_diffusion_coverage_miss():
    """The C1 failure in miniature: the diffusion's interval on a filtered
    record does not cover the truth; with the loss declared it does, and the
    loss's own uncertainty sits in the band at coefficient 1."""
    Xf = ou(n=150_000, ntraj=3, seed=5, over=OVER)
    X = Xf[..., ::OVER]
    Xr = instrument(Xf)
    t = np.arange(X.shape[1]) * DT
    bl = band_loss(Xr[0], FS, fit_range=(100.0, CUTOFF), alpha=ALPHA, f_diode=F_DIODE)
    kw = dict(half=2000, ws=("1",), drift_envelope=None)
    raw = build_qv_rows(t, Xr, ("1",), **kw)
    dec = build_qv_rows(t, Xr, ("1",), band_loss=(bl["r_total"], bl["r_se"]), **kw)
    assert dec.band_loss_built == bl["r_total"] and dec.loss_se is not None
    assert np.allclose(dec.y, raw.y / bl["r_total"])
    assert np.allclose(dec.var, raw.var / bl["r_total"] ** 2)
    assert any("declared band loss" in n for n in dec.notes)
    dmax = THETA * X.std() * 3
    r_raw = certify_diffusion(raw, delta=0.05, drift_max=dmax, seed=0)
    r_dec = certify_diffusion(dec, delta=0.05, drift_max=dmax, seed=0)
    lo, hi = r_raw["admissible"]["diffusion:1"]
    assert hi < B2                                     # the miss, reproduced
    lo, hi = r_dec["admissible"]["diffusion:1"]
    assert lo <= B2 <= hi                              # closed by the declaration
    assert np.all(certify_diffusion.__globals__["QvBand"](
        dec.var, dec.lam, len(dec.y), dec.dt, y=dec.y, loss_se=dec.loss_se)(None)
        >= dec.loss_se)
    with pytest.raises(ValueError, match="band_loss"):
        build_qv_rows(t, Xr, ("1",), band_loss=(1.5, 0.0), **kw)


def test_lockin_reads_a_drive_and_drive_scale_is_the_fluid_drag_response():
    """An active calibration in miniature: the stage moves at 38 Hz with a known
    amplitude in nm; the bead follows with the fluid-drag transfer f_d /
    sqrt(f_d^2 + f_c^2) and is recorded in VOLTS through an unknown scale. The
    scale must come back without any thermal quantity entering."""
    fd, A_stage, Rd = 38.12, 463.0, 4.57e3              # nm, nm/V
    fc = THETA / (2 * np.pi)
    x = ou(n=300_000, seed=8)[0]                        # thermal motion, nm
    t = np.arange(len(x)) * DT
    stage = A_stage * np.sin(2 * np.pi * fd * t)
    resp_nm = A_stage * fd / np.sqrt(fd ** 2 + fc ** 2) * np.sin(2 * np.pi * fd * t - 0.3)
    volts = (x + resp_nm) / Rd
    li = lockin(stage, FS, fd)
    assert abs(li["amplitude"] / A_stage - 1) < 1e-4 and li["power_share"] > 0.999
    ds = drive_scale(stage, volts, FS, fc)
    assert abs(ds["f_drive_hz"] - fd) < 0.05
    assert abs(ds["scale"] / Rd - 1) < 0.03            # the non-circular Rd
    # the drive is removed from the record before any stochastic claim
    resid_nm = ds["residual"] * ds["scale"]
    assert axis_gate(resid_nm, DT, THETA)["passed"]
    assert abs(resid_nm.var() / x.var() - 1) < 0.05


def test_diode_response_is_the_declared_model():
    f = np.array([0.0, F_DIODE, 1e9])
    g = diode_response(f, ALPHA, F_DIODE)
    assert g[0] == pytest.approx(1.0) and g[1] == pytest.approx(ALPHA ** 2 + (1 - ALPHA ** 2) / 2)
    assert g[2] == pytest.approx(ALPHA ** 2, rel=1e-6)


# ---------------------------------------------------------------------------
# C3: what the record determines about kappa and gamma separately, and the
# degeneracy that decides it (docs/CASE_STUDY_TWEEZERS_C3.md).

ETA, A_BEAD_UM = 0.001002, 0.52                 # water at 20 C, a 1.04 um bead
N_REAL = 780_000                                # the C-Trap active record's length:
# COMMON_MODE_SPREAD must sit above the ACF ratio's OWN sampling scatter, and that
# scatter is length-dependent (measured: 4.8% of the ratio at 400k samples, 2.3% at
# 780k), so a two-axis agreement test is only meaningful at a realistic length
G_BULK = 6 * np.pi * ETA * A_BEAD_UM * 1e-6
KBT = 1.380649e-23 * 293.15
RD_UM_PER_V = 4.568                             # truth: the displacement scale
F_DRIVE, A_STAGE_UM = 38.153, 0.527             # an active calibration's stage


def bead(drag_over_bulk=1.0, kappa_pN_nm=0.1194, n=400_000, seed=4,
         driven=True, alpha=0.451, f_diode=11684.0, over=OVER):
    """A trapped bead whose drag is `drag_over_bulk` x bulk Stokes, observed
    through the declared detector and an anti-alias filter, in VOLTS.

    Returns (stage, volts, truth). The trap is linear, so the driven and
    thermal parts superpose exactly and the drive is the analytic steady state
    x_bead = A f_d / sqrt(f_d^2 + f_c^2) with f_c the TRUE corner."""
    gamma = drag_over_bulk * G_BULK
    kappa = kappa_pN_nm * 1e-3
    theta = kappa / gamma
    b2 = 2 * KBT / gamma * 1e18                              # nm^2/s
    x = ou(n=n, seed=seed, theta=theta, b2=b2, over=over)[0]
    tt = np.arange(len(x)) * (DT / over)
    w = 2 * np.pi * F_DRIVE
    stage = A_STAGE_UM * np.sin(w * tt)                      # um, as the file records it
    if driven:
        x = x + (A_STAGE_UM * 1e3 * w / np.hypot(w, theta)
                 * np.sin(w * tt + np.arctan2(theta, w)))
    volts = anti_alias(diode(x, alpha=alpha, f_diode=f_diode, over=over),
                       over=over)[::over] / (RD_UM_PER_V * 1e3)
    truth = {"theta": theta, "gamma": gamma, "kappa": kappa, "b2": b2,
             "fc_true": theta / (2 * np.pi),
             "fc_bulk_referenced": kappa / (2 * np.pi * G_BULK),
             "drag_over_bulk": drag_over_bulk}
    return stage[::over], volts, truth


def _fit(volts, alpha=0.451, f_diode=11684.0):
    return retention(volts, FS, fit_range=(10.0, CUTOFF), alpha=alpha,
                     f_diode=f_diode)


def test_retention_serves_two_consumers_from_one_measurement():
    """Variance and quadratic variation weight frequency differently, so one
    measured response must give two different fractions -- and an estimator
    that used one where it needed the other would be wrong by their ratio."""
    _, volts, _ = bead(driven=False)
    r = _fit(volts)
    assert r.variance > 2.5 * r.fraction(1)        # ~0.87 against ~0.32
    f = r.fraction(1, factors=True)
    assert abs(f["r_nyquist"] * f["r_diode"] * f["r_antialias"] - f["r_total"]) < 1e-12
    # the exact Lorentzian moments must agree with brute-force integration
    F = np.linspace(0.0, 400.0 * FS, 2_000_001)
    lor = r.amplitude / (r.fc_fit ** 2 + F ** 2)
    assert abs(r._full(None) / np.trapezoid(lor, F) - 1) < 2e-3
    assert abs(r._full(1.0) / np.trapezoid(4 * np.sin(np.pi * F / FS) ** 2 * lor, F)
               - 1) < 2e-3
    bl = band_loss(volts, FS, fit_range=(10.0, CUTOFF), alpha=0.451, f_diode=11684.0)
    assert bl["r_var"] == pytest.approx(r.variance)
    assert bl["r_total"] == pytest.approx(r.fraction(1))


def test_the_diffusion_is_flat_in_stride_once_each_stride_carries_its_own_weight():
    """'Step past the filter' becomes a measured correction rather than a hope:
    the lag-s increment has weight 4 sin^2(pi f s / fs), so each stride has its
    OWN retention. Corrected that way, b^2 is flat in the stride, and its
    flatness is the self-check."""
    for drag in (1.0, 2.3):
        _, volts, tr = bead(drag_over_bulk=drag, driven=False)
        r = _fit(volts)
        th = acf_timescale(volts, DT, min_lag=8, max_lag=200)["theta"]
        got = []
        for s in (4, 8, 16, 32):
            m = float(np.mean((volts[s:] - volts[:-s]) ** 2)) * (RD_UM_PER_V * 1e3) ** 2
            got.append(m * th / (1 - np.exp(-th * s * DT)) / r.fraction(s))
        assert max(got) / min(got) - 1 < 0.06            # flat
        assert abs(np.mean(got) / tr["b2"] - 1) < 0.06   # and right


@pytest.mark.parametrize("drag", [1.0, 1.44, 2.30])
def test_the_drive_scale_is_degenerate_and_a_declared_corner_returns_the_assumed_drag(drag):
    """ONE drive frequency determines only Rd * f_c. Supplying the corner from a
    calibration whose reference drag is the quantity under test returns that
    assumption: the scale comes back wrong by exactly the drag ratio, which is
    the error C2 shipped. Taken from the record, it is right."""
    stage, volts, tr = bead(drag_over_bulk=drag)
    ds = drive_scale(stage, volts, FS, declared_f_c=tr["fc_bulk_referenced"],
                     acf_kw={"min_lag": 8, "max_lag": 200})
    assert ds["f_c_source"] == "record"
    assert abs(ds["f_c_record_hz"] / tr["fc_true"] - 1) < 0.05
    assert abs(ds["scale"] / RD_UM_PER_V - 1) < 0.05            # the record's corner
    assert abs(ds["degeneracy_factor"] - 1 / drag) < 0.06       # the declared one
    assert abs(ds["scale_at_declared_f_c"] / (RD_UM_PER_V / drag) - 1) < 0.06


def test_local_drag_recovers_a_known_drag_and_stiffness():
    for drag in (1.0, 1.44, 2.30):
        stage, volts, tr = bead(drag_over_bulk=drag)
        ds = drive_scale(stage, volts, FS, acf_kw={"min_lag": 8, "max_lag": 200})
        r = _fit(ds["residual"])
        th = acf_timescale(ds["residual"], DT, min_lag=8, max_lag=200)["theta"]
        out = local_drag(th, ds["scale"], float(ds["residual"].var()), 20.0,
                         r_var=r.variance, radius_um=A_BEAD_UM, viscosity=ETA)
        assert abs(out["drag_over_bulk"] / drag - 1) < 0.10
        assert abs(out["kappa_N_per_m"] / tr["kappa"] - 1) < 0.10
        assert abs(out["gamma_kg_per_s"] / tr["gamma"] - 1) < 0.10
        # the variance retention is not optional: using the QV one instead is
        # wrong by their ratio, which is a factor of three
        bad = local_drag(th, ds["scale"], float(ds["residual"].var()), 20.0,
                         r_var=r.fraction(1), radius_um=A_BEAD_UM, viscosity=ETA)
        assert bad["drag_over_bulk"] / out["drag_over_bulk"] < 0.5


def test_faxen_inverts_inside_its_domain_and_refuses_outside_it():
    for lam in (0.2, 0.4, 0.55):
        F = faxen_drag_ratio(lam)
        got = faxen_height(F, A_BEAD_UM)
        assert abs(got["lambda"] - lam) < 1e-6
        assert abs(got["height_um"] - A_BEAD_UM / lam) < 1e-6
    assert faxen_height(faxen_drag_ratio(0.9), A_BEAD_UM)["height_um"] is None
    assert faxen_height(0.8, A_BEAD_UM)["height_um"] is None


def test_the_record_gate_tells_a_contaminant_from_a_different_drag():
    """The distinction one axis cannot draw. Equipartition does not involve the
    drag, so a record whose drag is not the calibration's reference has healthy
    variances on EVERY axis and one common timescale factor; a contaminated axis
    has neither."""
    kappa_pN = 0.1194
    axes = {}
    for name, drag in (("x", 2.30), ("y", 2.30)):
        _, volts, tr = bead(drag_over_bulk=drag, driven=False, n=N_REAL,
                            seed={"x": 41, "y": 42}[name])
        r = _fit(volts)
        var_nm2 = float(volts.var()) * (RD_UM_PER_V * 1e3) ** 2
        axes[name] = axis_gate(volts, DT, 2 * np.pi * tr["fc_bulk_referenced"],
                               min_lag=8, max_lag=200,
                               var_ratio=equipartition_ratio(
                                   var_nm2, kappa_pN * 1e-3, 20.0, r.variance))
    v = gate_record(axes)
    assert v["verdict"] == "common-mode"
    assert abs(v["common_ratio"] * 2.30 - 1) < 0.10
    assert v["spread"] < COMMON_MODE_SPREAD
    assert all(not a["passed"] and a["var_ok"] for a in axes.values())
    # and the gate stops exactly there: kappa/gamma is ONE number and they are
    # two, so naming the drag is `attribute_deviation`'s job, not the gate's
    assert "drag_over_reference" not in v and "kappa/gamma" in v["note"]

    # now contaminate one axis with a slow line the calibration never sees
    _, volts, tr = bead(driven=False, seed=17)
    t = np.arange(len(volts)) * DT
    dirty = volts + 2 * volts.std() * np.sin(2 * np.pi * 5.0 * t)
    r = _fit(dirty)
    var_nm2 = float(dirty.var()) * (RD_UM_PER_V * 1e3) ** 2
    axes["x"] = axis_gate(dirty, DT, 2 * np.pi * tr["fc_bulk_referenced"],
                          min_lag=8, max_lag=200,
                          var_ratio=equipartition_ratio(var_nm2, kappa_pN * 1e-3,
                                                        20.0, r.variance))
    v2 = gate_record(axes)
    assert v2["verdict"] == "contaminated" and v2["contaminated"] == ["x"]


def test_a_clean_record_gates_clean_and_one_axis_is_never_common_mode():
    kappa_pN = 0.1194
    axes = {}
    for name in ("x", "y"):
        _, volts, tr = bead(driven=False, n=N_REAL, seed={"x": 43, "y": 44}[name])
        r = _fit(volts)
        var_nm2 = float(volts.var()) * (RD_UM_PER_V * 1e3) ** 2
        axes[name] = axis_gate(volts, DT, 2 * np.pi * tr["fc_bulk_referenced"],
                               min_lag=8, max_lag=200,
                               var_ratio=equipartition_ratio(
                                   var_nm2, kappa_pN * 1e-3, 20.0, r.variance))
    assert gate_record(axes)["verdict"] == "clean"
    # a single off axis is never promoted to a drag statement
    _, volts, tr = bead(drag_over_bulk=2.30, driven=False, seed=5)
    r = _fit(volts)
    var_nm2 = float(volts.var()) * (RD_UM_PER_V * 1e3) ** 2
    one = {"x": axis_gate(volts, DT, 2 * np.pi * tr["fc_bulk_referenced"],
                          min_lag=8, max_lag=200,
                          var_ratio=equipartition_ratio(var_nm2, kappa_pN * 1e-3,
                                                        20.0, r.variance))}
    assert gate_record(one)["verdict"] == "unresolved"


KAPPA_REF_PN = 0.1194
THETA_REF = (KAPPA_REF_PN * 1e-3) / G_BULK
B2_REF_V = 2 * KBT / (G_BULK * (RD_UM_PER_V * 1e3 * 1e-9) ** 2)   # V^2/s: the
# SDE coefficient b^2 = 2 D, which is what `realized_diffusion` returns
VAR_REF_V = KBT / (KAPPA_REF_PN * 1e-3 * (RD_UM_PER_V * 1e3 * 1e-9) ** 2)


def _ratios(volts, theta_ref=THETA_REF):
    """The three scale-free ratios a record offers against its calibration."""
    r = _fit(volts)
    th = acf_timescale(volts, DT, min_lag=8, max_lag=200)["theta"]
    d = realized_diffusion(volts, DT, th, r)
    return (th / theta_ref, d["plateau"] / B2_REF_V,
            float(volts.var()) / (VAR_REF_V * r.variance), d)


def test_realized_diffusion_is_flat_and_needs_no_displacement_scale():
    _, volts, tr = bead(drag_over_bulk=1.44, driven=False)
    r = _fit(volts)
    th = acf_timescale(volts, DT, min_lag=8, max_lag=200)["theta"]
    d = realized_diffusion(volts, DT, th, r)
    assert d["spread"] < 0.06                                   # the self-check
    assert abs(d["plateau"] / (tr["b2"] / (RD_UM_PER_V * 1e3) ** 2) - 1) < 0.06
    assert max(d["strides"]) * DT * th <= 1.0 + 1e-9            # capped, not extrapolated
    # scaling the signal scales b^2 by the square and the RATIO not at all
    d2 = realized_diffusion(3.0 * volts, DT, th, _fit(3.0 * volts))
    assert abs(d2["plateau"] / d["plateau"] / 9.0 - 1) < 0.02


@pytest.mark.parametrize("drag", [1.44, 2.30])
def test_the_three_ratios_name_the_drag_and_only_the_drag(drag):
    _, volts, _ = bead(drag_over_bulk=drag, driven=False)
    th_r, d_r, v_r, _ = _ratios(volts)
    a = attribute_deviation(th_r, d_r, v_r)
    assert a["verdict"] == "drag"
    assert abs(a["hypotheses"]["drag"]["value"] / drag - 1) < 0.10
    assert a["hypotheses"]["stiffness"]["max_residual"] > 5 * \
        a["hypotheses"]["drag"]["max_residual"]


def test_the_three_ratios_tell_a_stiffness_change_from_a_drag_change():
    """The same timescale deficit, two different causes: only the diffusion and
    the variance tell them apart, and both are read in the detector's own units."""
    c = 0.70
    _, volts, _ = bead(drag_over_bulk=1.0, kappa_pN_nm=c * KAPPA_REF_PN, driven=False)
    th_r, d_r, v_r, _ = _ratios(volts)
    a = attribute_deviation(th_r, d_r, v_r)
    assert a["verdict"] == "stiffness"
    assert abs(a["hypotheses"]["stiffness"]["value"] / c - 1) < 0.10
    assert abs(th_r / c - 1) < 0.10 and abs(d_r - 1) < 0.10      # D did NOT move


def test_attribution_abstains_when_no_single_change_explains_the_record():
    assert attribute_deviation(0.7, 1.4, 0.3)["verdict"] == "unattributed"
    assert attribute_deviation(0.0, 1.0, 1.0)["verdict"] == "unreadable"


def test_faxen_refuses_a_height_when_the_drag_is_consistent_with_bulk():
    """A drag ratio of 1.02 +- 0.10 inverts to a 12 um height and determines no
    such thing -- the wall's effect dies with distance, so the interval is
    bounded below only. Measured: this is what the drive route returned on the
    C-Trap before the uncertainty was carried."""
    got = faxen_height(1.025, A_BEAD_UM, 0.10)
    assert got["height_um"] is None and got["height_lower_bound_um"] is not None
    assert "consistent with bulk" in got["note"]
    tight = faxen_height(1.42, A_BEAD_UM, 0.02)
    assert tight["height_um"] is not None and 0.8 < tight["height_um"] < 1.1

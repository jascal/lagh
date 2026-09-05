"""What a real instrument does to a record, declared and measured
(lagh/instrument.py; docs/CASE_STUDY_TWEEZERS_C1.md §7). Every test simulates
the C-Trap bead of C1 -- theta 3297 /s, b^2 2.007e5 nm^2/s, 78125 Hz -- and
applies the distortions ourselves, so the truth is known."""

import numpy as np
import pytest

from lagh.instrument import (acf_timescale, axis_gate, band_loss, diode_response,
                             drive_scale, lockin)
from lagh.ito import build_qv_rows, certify_diffusion

THETA, B2, FS = 3297.0, 2.007e5, 78125.0
DT = 1.0 / FS
ALPHA, F_DIODE, CUTOFF = 0.42, 8400.0, 23000.0


OVER = 8            # the instrument filters the CONTINUOUS signal, then samples


def ou(n=400_000, ntraj=1, seed=3, theta=THETA, b2=B2, over=1):
    """Discrete OU at `over` x the sampling rate (exact transition density)."""
    rng = np.random.default_rng(seed)
    dt = DT / over
    dec = np.exp(-theta * dt)
    sd = np.sqrt(b2 * (1 - dec ** 2) / (2 * theta))
    X = np.empty((ntraj, n * over))
    X[:, 0] = rng.standard_normal(ntraj) * np.sqrt(b2 / (2 * theta))
    z = rng.standard_normal((ntraj, n * over - 1))
    for k in range(1, n * over):
        X[:, k] = dec * X[:, k - 1] + sd * z[:, k - 1]
    return X


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

"""What a real instrument does to a stochastic record, DECLARED and MEASURED
before any claim is made (docs/CASE_STUDY_TWEEZERS_C1.md §7, the C2 list).

Every simulated system in the stochastic arc was observed exactly. A real
trapped bead is not: the detector low-passes the position (the photodiode's
parasitic filtering, with parameters the instrument itself fits and declares),
the acquisition anti-alias-filters before storage, the sampling band ends at
Nyquist while the process's quadratic variation does not, and one axis may carry
a slow contaminant the calibration's own fit range never looks at. C1 measured
each of these on the LUMICKS C-Trap and found the arc's diffusion reading 32% of
truth and its drift 46% low, with zero confident-wrong only because everything
abstained.

Three tools, in the order C1 says to apply them:

* `axis_gate` -- the autocorrelation TIMESCALE of the record against the
  calibration's own corner frequency. A timescale is scale-free (no Rd, no b^2)
  and is read at lags beyond the detector's filter, so it is immune to both of
  the instrument's distortions; a slow contaminant below the calibration's fit
  floor shows up here and nowhere else (measured: one axis at 4% of its
  nominal timescale, 1.9x thermal variance, invisible to a 100 Hz-floored PSD
  fit). Gate first; nothing below is claimed on an axis that fails.

* `band_loss` -- the fraction of the process's quadratic variation the RECORD
  can contain, as a product of three declared or measured factors:
      r_nyquist   geometric: the sampled increments see the process through
                  w(f) = 4 sin^2(pi f dt), whose full-band integral is the
                  true QV; the record holds only |f| < f_N
      r_diode     the instrument's DECLARED detector model
                  g(f) = alpha^2 + (1 - alpha^2) / (1 + (f/f_diode)^2)
      r_antialias MEASURED: the recorded PSD over the declared model above the
                  calibration's fit ceiling, where the model is not fitted
  Measured on the C-Trap: 0.77 x 0.44 x 0.94 = 0.32, and the raw QV read 0.325
  of truth -- the corrected diffusion is 1.01-1.02 of truth on both axes. This
  is a DECLARED INPUT to the diffusion rows (`ito.build_qv_rows(band_loss=)`),
  the way sigma_obs is: it is measured from the record and the instrument's own
  declaration, never fitted alongside the claim.

* `lockin` / `drive_scale` -- an ACTIVE calibration drives the stage at a known
  amplitude; the bead's response at the drive frequency is the fluid-drag
  response A_stage * f_d / sqrt(f_d^2 + f_c^2), read by lock-in against the
  recorded stage position. Its ratio to the response in volts is a displacement
  scale Rd that does NOT pass through the thermal motion -- the one number on
  this instrument that breaks the passive calibration's circularity (C1 §1).
  The same lock-in removes the drive from the record before any stochastic
  claim is made on it.
"""
from __future__ import annotations

import numpy as np


# ------------------------------------------------------------------ the axis gate

def acf_timescale(x, dt: float, *, max_lag: int = 120, min_lag: int = 10,
                  floor: float = 0.05) -> dict:
    """theta from the autocorrelation decay, exp(-theta * lag * dt), fitted on
    lags in [min_lag, max_lag) while the ACF is above `floor`. Scale-free."""
    v = np.asarray(x, float) - np.mean(x)
    den = float(v @ v)
    if den <= 0 or len(v) < max_lag + 2:
        return {"theta": None, "note": "no variance or too short"}
    lags = np.arange(1, max_lag)
    acf = np.array([float(v[:-L] @ v[L:]) / den for L in lags])
    m = (lags >= min_lag) & (acf > floor)
    if m.sum() < 5:
        return {"theta": None, "acf_lag1": float(acf[0]),
                "note": "ACF decays too fast to fit above the floor"}
    slope, icpt = np.polyfit(lags[m] * dt, np.log(acf[m]), 1)
    return {"theta": float(-slope), "acf_lag1": float(acf[0]),
            "n_lags": int(m.sum()), "intercept": float(np.exp(icpt))}


def axis_gate(x, dt: float, theta_ref: float, *, tol: float = 1.25,
              **kw) -> dict:
    """PASS iff the record's autocorrelation timescale agrees with the declared
    one within a factor `tol`. `theta_ref` is the calibration's 2 pi f_c, a
    timescale the instrument fitted from the PSD within its own band; a record
    whose decay disagrees with it carries something the calibration did not see
    (a slow contaminant, an un-removed drive, a second bead)."""
    a = acf_timescale(x, dt, **kw)
    th = a.get("theta")
    if th is None or not np.isfinite(th) or th <= 0:
        return {"passed": False, "theta_acf": th, "theta_ref": float(theta_ref),
                "ratio": None, "tol": tol,
                "note": a.get("note", "no timescale could be read")}
    ratio = th / float(theta_ref)
    ok = bool(1.0 / tol <= ratio <= tol)
    return {"passed": ok, "theta_acf": th, "theta_ref": float(theta_ref),
            "ratio": float(ratio), "tol": tol, "acf_lag1": a["acf_lag1"],
            "note": ("timescale agrees with the calibration" if ok else
                     f"timescale is {ratio:.3g}x the calibration's: the axis "
                     "carries something the calibration's fit range never saw")}


# ---------------------------------------------------------------- the band loss

def diode_response(f, alpha: float, f_diode: float):
    """|H|^2 of the photodiode: the instrument's own declared model."""
    f = np.asarray(f, float)
    return alpha ** 2 + (1 - alpha ** 2) / (1 + (f / f_diode) ** 2)


def _lorentzian(f, A, fc):
    return A / (fc ** 2 + np.asarray(f, float) ** 2)


def band_loss(x, fs: float, *, fit_range=(100.0, 23000.0), alpha: float,
              f_diode: float, nperseg: int = 2 ** 16, halves: bool = True) -> dict:
    """The fraction of the process's quadratic variation the record contains.

    Fits the Lorentzian x diode SHAPE to the recorded PSD inside the
    calibration's fit range (free amplitude and corner: no calibration scale
    enters), takes the recorded PSD over that model ABOVE the range as the
    measured anti-alias response (clipped to [0, 1]: a filter removes, it does
    not add), and integrates all three factors against the discrete-increment
    weight w(f) = 4 sin^2(pi f / fs). Returns the factors and their product
    `r_total`, with `r_se` from the two halves of the record when `halves`.

    The denominator is the FULL-BAND integral of w x Lorentzian: the sampled
    increments see every frequency (w is periodic in f), so the process's true
    QV includes what lies above Nyquist, which the record can never hold.
    Truncating there was a measured 23% error."""
    from scipy.optimize import curve_fit
    from scipy.signal import welch
    x = np.asarray(x, float)
    dt = 1.0 / fs
    lo, hi = float(fit_range[0]), float(fit_range[1])

    def one(seg):
        fr, P = welch(seg - np.mean(seg), fs=fs, nperseg=min(nperseg, len(seg)))
        inb = (fr > lo) & (fr < hi)
        g = diode_response(fr, alpha, f_diode)
        p0 = [float(P[inb][0] * fr[inb][0] ** 2), float(0.1 * hi)]
        (A, fc), _ = curve_fit(lambda fq, A, fc: _lorentzian(fq, A, fc)
                               * diode_response(fq, alpha, f_diode),
                               fr[inb], P[inb], p0=p0, maxfev=20000)
        fc = abs(float(fc))
        model = _lorentzian(fr, A, fc)
        H2 = np.where(fr >= hi, P / np.maximum(model * g, 1e-300), 1.0)
        H2 = np.clip(H2, 0.0, 1.0)
        w = 4.0 * np.sin(np.pi * fr * dt) ** 2
        F = np.linspace(0.0, 400.0 * fs, 2_000_001)
        full = 2.0 * np.trapezoid(4.0 * np.sin(np.pi * F * dt) ** 2
                                  * _lorentzian(F, A, fc), F)
        trunc = 2.0 * np.trapezoid(w * model, fr)
        with_diode = 2.0 * np.trapezoid(w * model * g, fr)
        with_all = 2.0 * np.trapezoid(w * model * g * H2, fr)
        return {"fc_fit_hz": fc, "r_nyquist": trunc / full,
                "r_diode": with_diode / trunc, "r_antialias": with_all / with_diode,
                "r_total": with_all / full,
                "psd_ratio_above_range": float(np.mean(H2[fr >= hi]))}

    out = one(x)
    if halves and len(x) >= 4 * nperseg:
        n2 = len(x) // 2
        a, b = one(x[:n2]), one(x[n2:])
        out["r_se"] = float(abs(a["r_total"] - b["r_total"]) / 2.0)
        out["r_halves"] = [a["r_total"], b["r_total"]]
    else:
        out["r_se"] = float(0.05 * out["r_total"])      # declared, not measured
        out["r_se_note"] = "record too short for a split estimate; 5% declared"
    out.update(alpha=float(alpha), f_diode_hz=float(f_diode),
               fit_range_hz=[lo, hi], fs_hz=float(fs))
    return out


# --------------------------------------------------------- the active calibration

def drive_frequency(ref, fs: float) -> float:
    """The dominant line of a recorded drive, from its spectrum, refined below
    the bin width by a quadratic fit of the log-magnitude at the peak (a
    10 s record has 0.1 Hz bins; the lock-in needs better than that)."""
    s = np.asarray(ref, float) - np.mean(ref)
    spec = np.abs(np.fft.rfft(s * np.hanning(len(s))))
    fr = np.fft.rfftfreq(len(s), 1.0 / fs)
    k = int(np.argmax(spec[1:])) + 1
    if 1 <= k < len(spec) - 1 and spec[k - 1] > 0 and spec[k + 1] > 0:
        a, b, c = np.log(spec[k - 1]), np.log(spec[k]), np.log(spec[k + 1])
        den = a - 2 * b + c
        shift = 0.5 * (a - c) / den if den != 0 else 0.0
        return float(fr[k] + shift * (fr[1] - fr[0]))
    return float(fr[k])


def lockin(x, fs: float, f_drive: float) -> dict:
    """Least-squares lock-in of `x` on sin/cos at `f_drive`: amplitude, phase,
    and the record with that line removed. The removal subtracts ONE sinusoid;
    what it leaves (harmonics, stage jitter) is for the axis gate to judge."""
    x = np.asarray(x, float)
    t = np.arange(len(x)) / fs
    M = np.column_stack([np.sin(2 * np.pi * f_drive * t),
                         np.cos(2 * np.pi * f_drive * t)])
    c, *_ = np.linalg.lstsq(M, x - np.mean(x), rcond=None)
    fit = M @ c
    return {"amplitude": float(np.hypot(c[0], c[1])),
            "phase": float(np.arctan2(c[1], c[0])),
            "residual": x - fit,
            "power_share": float(1.0 - np.var(x - fit) / max(np.var(x), 1e-300))}


def drive_scale(stage, response, fs: float, f_c: float) -> dict:
    """Displacement scale from an active calibration: stage amplitude over the
    bead's lock-in response, corrected for the trap's fluid-drag transfer at the
    drive frequency, |x_bead / x_stage| = f_d / sqrt(f_d^2 + f_c^2). The units
    are those of `stage` per unit of `response`; `f_c` is a timescale, so no
    displacement scale enters -- this is the non-circular Rd."""
    fd = drive_frequency(stage, fs)
    n = min(len(stage), len(response))
    ls, lr = lockin(stage[:n], fs, fd), lockin(response[:n], fs, fd)
    transfer = fd / np.sqrt(fd ** 2 + float(f_c) ** 2)
    # x_bead = transfer * x_stage, so Rd = (A_stage * transfer) / A_response[V]
    return {"f_drive_hz": fd, "stage_amplitude": ls["amplitude"],
            "response_amplitude": lr["amplitude"], "transfer": float(transfer),
            "scale": float(ls["amplitude"] * transfer / lr["amplitude"]),
            "response_power_share": lr["power_share"],
            "residual": lr["residual"], "n": int(n)}

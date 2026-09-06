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
  fit). Gate first; stochastic law claims still refuse on a failed axis.
  Separate C3/C4 diagnostics retain the gate and report deviations, never
  upgrading its failure to a law certificate.

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

from dataclasses import dataclass, field

import numpy as np

KB = 1.380649e-23


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
              var_ratio: float | None = None, var_tol: float = 1.35,
              **kw) -> dict:
    """PASS iff the record's autocorrelation timescale agrees with the declared
    one within a factor `tol`. `theta_ref` is the calibration's 2 pi f_c, a
    timescale the instrument fitted from the PSD within its own band; a record
    whose decay disagrees with it carries something the calibration did not see
    (a slow contaminant, an un-removed drive, a second bead) -- OR a drag that
    is not the calibration's reference drag, which is a MEASUREMENT and not a
    fault. `var_ratio` (see `equipartition_ratio`) is what tells the two apart
    on a single axis, and `gate_record` is what tells them apart across a
    record; both are reported here, neither is guessed at when absent."""
    a = acf_timescale(x, dt, **kw)
    th = a.get("theta")
    if th is None or not np.isfinite(th) or th <= 0:
        return {"passed": False, "theta_acf": th, "theta_ref": float(theta_ref),
                "ratio": None, "tol": tol, "var_ratio": var_ratio,
                "note": a.get("note", "no timescale could be read")}
    ratio = th / float(theta_ref)
    ok = bool(1.0 / tol <= ratio <= tol)
    var_ok = None if var_ratio is None else \
        bool(1.0 / var_tol <= float(var_ratio) <= var_tol)
    if ok:
        note = "timescale agrees with the calibration"
    elif var_ok is False:
        note = (f"timescale is {ratio:.3g}x the calibration's AND the "
                f"equipartition variance is {var_ratio:.3g}x its own expectation: "
                "the axis carries something the calibration's fit range never saw")
    elif var_ok is True:
        note = (f"timescale is {ratio:.3g}x the calibration's while the "
                "equipartition variance is as expected: kappa/gamma differs from "
                "the calibration's REFERENCE ratio -- see gate_record")
    else:
        note = (f"timescale is {ratio:.3g}x the calibration's; no equipartition "
                "variance was supplied, so contamination and a different drag "
                "are not separated")
    return {"passed": ok, "theta_acf": th, "theta_ref": float(theta_ref),
            "ratio": float(ratio), "tol": tol, "acf_lag1": a["acf_lag1"],
            "var_ratio": None if var_ratio is None else float(var_ratio),
            "var_ok": var_ok, "note": note}


def equipartition_ratio(var_measured, kappa, temperature_c: float,
                        r_var: float = 1.0) -> float:
    """Measured position variance over its equipartition expectation kB T / kappa,
    with the record's own variance RETENTION divided out (`Retention.variance`).

    The discriminator's whole content is that this ratio does NOT depend on the
    drag: equipartition fixes the variance from kappa and temperature alone, so a
    record whose drag differs from the calibration's reference still has a ratio
    of 1 while its TIMESCALE moves by the drag factor. Contamination moves both.

    `var_measured` [nm^2] and `kappa` [N/m] must be in the same length units the
    displacement scale provides; the ratio is dimensionless."""
    expected = KB * (float(temperature_c) + 273.15) / float(kappa) * 1e18
    return float(var_measured) / (expected * float(r_var))


COMMON_MODE_SPREAD = 0.10          # axes agreeing this closely are common-mode


def gate_record(axes: dict, *, tol: float = 1.25,
                spread: float = COMMON_MODE_SPREAD) -> dict:
    """One verdict for a RECORD from its per-axis gates -- the distinction a
    single axis cannot draw.

    `axes` maps an axis name to its `axis_gate` result (with `var_ratio` set).
    Three outcomes:

      clean          every axis's timescale agrees with the calibration.
      contaminated   an axis whose equipartition variance is also wrong. That
                     axis carries something the calibration never saw, and it is
                     named; the others are judged on their own.
      common-mode    EVERY axis's timescale is off by the SAME factor while
                     every equipartition variance is as expected. Equipartition
                     does not involve the drag and the timescale is kappa/gamma,
                     so a common factor with healthy variances says the record's
                     kappa/gamma is not the calibration's REFERENCE ratio. That
                     is a measurement of 1/(gamma/gamma_reference) and must not
                     be thrown away as a refusal (measured: the C-Trap's
                     near-surface record, 0.70 on both parallel axes, against a
                     contaminated axis on the passive record at 0.04 with 3.6x
                     the equipartition variance).

    Two axes are the minimum for `common-mode`: one axis off by a factor is
    indistinguishable from one axis being wrong, and calling it a measurement
    would be the unsound direction."""
    named = {k: v for k, v in axes.items()
             if v.get("ratio") is not None and np.isfinite(v["ratio"])
             and v["ratio"] > 0}
    if not named:
        return {"verdict": "unreadable", "axes": list(axes),
                "note": "no axis yielded a timescale"}
    bad_var = sorted(k for k, v in named.items() if v.get("var_ok") is False)
    if bad_var:
        return {"verdict": "contaminated", "contaminated": bad_var,
                "ratios": {k: v["ratio"] for k, v in named.items()},
                "note": (f"{', '.join(bad_var)}: the equipartition variance is "
                         "wrong too, so the timescale deficit is not a drag "
                         "statement")}
    ratios = {k: float(v["ratio"]) for k, v in named.items()}
    lo, hi = min(ratios.values()), max(ratios.values())
    if len(named) != len(axes):
        return {"verdict": "unresolved", "ratios": ratios,
                "note": "an unreadable axis prevents a whole-record claim"}
    if not all(v.get("var_ok") is True for v in named.values()):
        return {"verdict": "unresolved", "ratios": ratios,
                "note": "missing equipartition evidence prevents a clean or common-mode claim"}
    if all(1.0 / tol <= r <= tol for r in ratios.values()):
        return {"verdict": "clean", "ratios": ratios,
                "note": "every axis agrees with the calibration's timescale"}
    if (len(ratios) >= 2 and hi / lo - 1.0 <= spread
            and all(v.get("var_ok") is True for v in named.values())):
        common = float(np.exp(np.mean(np.log(list(ratios.values())))))
        return {"verdict": "common-mode", "ratios": ratios,
                "common_ratio": common, "spread": float(hi / lo - 1.0),
                "note": ("every axis is off by the same factor with healthy "
                         f"equipartition variances: the record's kappa/gamma is "
                         f"{common:.3f} of the calibration's reference. WHICH of "
                         "kappa and gamma moved is not decided here -- the ratio "
                         "is one number and they are two; `attribute_deviation` "
                         "decides it from a second scale-free ratio")}
    return {"verdict": "unresolved", "ratios": ratios,
            "spread": float(hi / lo - 1.0),
            "note": ("axes disagree with the calibration by DIFFERENT factors "
                     "with healthy variances: not one drag, and not one "
                     "contaminant this test can name")}


# ---------------------------------------------------------------- the band loss

def diode_response(f, alpha: float, f_diode: float):
    """|H|^2 of the photodiode: the instrument's own declared model."""
    f = np.asarray(f, float)
    return alpha ** 2 + (1 - alpha ** 2) / (1 + (f / f_diode) ** 2)


def _lorentzian(f, A, fc):
    return A / (fc ** 2 + np.asarray(f, float) ** 2)


@dataclass
class Retention:
    """One measured instrument response, integrated against ANY weight.

    The record's spectrum is the process's spectrum times the sampling band,
    times the DECLARED detector model, times whatever else the acquisition did
    above the calibration's fit ceiling (MEASURED there, where the instrument's
    own model is not fitted and cannot be trusted). What fraction of a
    process quantity survives that depends on the weight the quantity puts on
    frequency, and different consumers weight it differently:

        variance          weight 1
        increment at lag s   weight 4 sin^2(pi f s / fs)

    so ONE measurement serves both, and their ratio is not free to be chosen
    per consumer. The denominators are the Lorentzian's own moments and are
    EXACT rather than numerically integrated to some multiple of the sampling
    rate: for S(f) = A/(fc^2 + f^2),

        int_0^inf S df                    = A pi / (2 fc)          = Var
        int_0^inf 4 sin^2(pi f tau) S df  = (A pi / fc)(1 - e^{-2 pi fc tau})

    which are exactly the Ornstein-Uhlenbeck variance and increment variance
    with theta = 2 pi fc -- so the denominator is the process's own moment, and
    the truncation that made an earlier version of this 23% wrong cannot occur.
    """
    freqs: np.ndarray
    model: np.ndarray
    diode: np.ndarray
    measured: np.ndarray
    fs: float
    amplitude: float
    fc_fit: float
    fit_range: tuple
    alpha: float
    f_diode: float
    notes: list = field(default_factory=list)

    def _full(self, stride: float | None) -> float:
        a, fc = self.amplitude, self.fc_fit
        if stride is None:
            return a * np.pi / (2.0 * fc)
        tau = stride / self.fs
        return (a * np.pi / fc) * (1.0 - np.exp(-2.0 * np.pi * fc * tau))

    def _weight(self, stride: float | None) -> np.ndarray:
        if stride is None:
            return np.ones_like(self.freqs)
        return 4.0 * np.sin(np.pi * self.freqs * stride / self.fs) ** 2

    def fraction(self, stride: float | None = 1.0, *, factors: bool = False):
        """The fraction of the process's weighted spectral content the record
        retains: `stride=None` for the variance, `stride=s` for the lag-s
        increment. With `factors`, the three named causes as well."""
        w = self._weight(stride)
        fr = self.freqs
        band = float(np.trapezoid(w * self.model, fr))
        diode = float(np.trapezoid(w * self.model * self.diode, fr))
        total = float(np.trapezoid(w * self.model * self.diode * self.measured, fr))
        full = self._full(stride)
        if not factors:
            return total / full
        return {"r_nyquist": band / full,
                "r_diode": diode / max(band, 1e-300),
                "r_antialias": total / max(diode, 1e-300),
                "r_total": total / full}

    @property
    def variance(self) -> float:
        return self.fraction(None)


def _fit_retention(x, fs, fit_range, alpha, f_diode, nperseg) -> Retention:
    from scipy.optimize import curve_fit
    from scipy.signal import welch
    x = np.asarray(x, float)
    lo, hi = float(fit_range[0]), float(fit_range[1])
    fr, P = welch(x - np.mean(x), fs=fs, nperseg=min(nperseg, len(x)))
    inb = (fr > lo) & (fr < hi)
    if inb.sum() < 8:
        raise ValueError(f"only {int(inb.sum())} spectral bins inside {fit_range}")
    g = diode_response(fr, alpha, f_diode)
    p0 = [float(P[inb][0] * fr[inb][0] ** 2), float(0.1 * hi)]
    (A, fc), _ = curve_fit(lambda fq, A, fc: _lorentzian(fq, A, fc)
                           * diode_response(fq, alpha, f_diode),
                           fr[inb], P[inb], p0=p0, maxfev=20000)
    fc = abs(float(fc))
    model = _lorentzian(fr, A, fc)
    # ABOVE the fit ceiling the instrument's model is not fitted, so the response
    # is measured there rather than trusted. Clipped to [0, 1]: a filter removes.
    H2 = np.clip(np.where(fr >= hi, P / np.maximum(model * g, 1e-300), 1.0), 0.0, 1.0)
    return Retention(fr, model, g, H2, float(fs), float(A), fc,
                     (lo, hi), float(alpha), float(f_diode))


def retention(x, fs: float, *, fit_range=(100.0, 23000.0), alpha: float,
              f_diode: float, nperseg: int = 2 ** 16) -> Retention:
    """Measure the instrument's response on this record; see `Retention`."""
    return _fit_retention(x, fs, fit_range, alpha, f_diode, nperseg)


def band_loss(x, fs: float, *, fit_range=(100.0, 23000.0), alpha: float,
              f_diode: float, nperseg: int = 2 ** 16, halves: bool = True) -> dict:
    """The fraction of the process's QUADRATIC VARIATION the record contains,
    with its three named factors and the VARIANCE retention alongside.

    `r_total` is `Retention.fraction(1)` -- the lag-1 increment, which is what a
    realized quadratic variation sums -- and `r_var` is `Retention.variance`,
    which is what an equipartition reading needs. They differ by a factor of
    three on a real record (measured on the C-Trap: 0.32 against 0.94), and an
    estimator that used one where it needed the other would be wrong by that.

    `r_se` is the spread over the two halves of the record, or a declared 5%
    when the record is too short to split."""
    x = np.asarray(x, float)
    r = _fit_retention(x, fs, fit_range, alpha, f_diode, nperseg)
    out = r.fraction(1.0, factors=True)
    out["r_var"] = r.variance
    out["fc_fit_hz"] = r.fc_fit
    out["psd_ratio_above_range"] = float(np.mean(r.measured[r.freqs >= r.fit_range[1]]))
    if halves and len(x) >= 4 * nperseg:
        n2 = len(x) // 2
        a = _fit_retention(x[:n2], fs, fit_range, alpha, f_diode, nperseg).fraction(1.0)
        b = _fit_retention(x[n2:], fs, fit_range, alpha, f_diode, nperseg).fraction(1.0)
        out["r_se"] = float(abs(a - b) / 2.0)
        out["r_halves"] = [a, b]
    else:
        out["r_se"] = float(0.05 * out["r_total"])
        out["r_se_note"] = "record too short for a split estimate; 5% declared"
    out.update(alpha=float(alpha), f_diode_hz=float(f_diode),
               fit_range_hz=list(r.fit_range), fs_hz=float(fs))
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


def drive_scale(stage, response, fs: float, f_c=None, *,
                declared_f_c: float | None = None, acf_kw: dict | None = None) -> dict:
    """Displacement scale from an active calibration, and the degeneracy that
    decides which corner frequency may be used to compute it.

    An active calibration drives the stage at a known amplitude; the bead
    follows through the fluid with the over-damped transfer

        |x_bead / x_stage|  =  f_d / sqrt(f_d^2 + f_c^2)

    so the lock-in response in volts determines the product Rd * transfer, and
    at f_d << f_c the transfer is f_d / f_c: ONE drive frequency determines only
    the PRODUCT Rd * f_c, never either alone. The corner must therefore come
    from somewhere else, and WHICH somewhere is not a detail:

    * from the RECORD (default, `f_c=None`): the bead's own relaxation, read off
      the drive-removed residual by `acf_timescale`. Scale-free, so it cannot
      import a calibration's assumptions.
    * from the CALIBRATION: only sound when the calibration's corner is already
      known to describe this record. A calibration whose corner is referenced to
      a BULK drag does not describe a record taken near a surface, and using it
      there returns the drag you assumed. Measured (jascal/lagh, tweezers C3):
      on a simulated bead whose drag is 2.30x bulk, the declared-corner scale is
      wrong by exactly 2.30x while the record's-corner scale is right to 1.4%.
      This is why the C2 reading of 1.48x was superseded.

    `declared_f_c` is reported alongside for exactly that comparison and never
    used for the returned `scale`. The units of `scale` are those of `stage` per
    unit of `response`; `f_c` is a timescale, so no displacement scale enters."""
    fd = drive_frequency(stage, fs)
    n = min(len(stage), len(response))
    ls, lr = lockin(stage[:n], fs, fd), lockin(response[:n], fs, fd)
    ac = acf_timescale(lr["residual"], 1.0 / fs, **(acf_kw or {}))
    fc_record = None if not ac.get("theta") else ac["theta"] / (2.0 * np.pi)
    if f_c is None:
        f_c, source = fc_record, "record"
    else:
        f_c, source = float(f_c), "declared"
    out = {"f_drive_hz": fd, "stage_amplitude": ls["amplitude"],
           "response_amplitude": lr["amplitude"],
           "response_power_share": lr["power_share"],
           "f_c_record_hz": fc_record, "f_c_used_hz": f_c, "f_c_source": source,
           "residual": lr["residual"], "n": int(n)}
    if f_c is None or not np.isfinite(f_c) or f_c <= 0:
        out.update(scale=None, transfer=None,
                   note="no corner frequency: the record yielded no timescale "
                        "and none was declared, so the scale is not determined")
        return out
    def _scale(fc):
        return float(ls["amplitude"] * (fd / np.hypot(fd, fc)) / lr["amplitude"])
    out["transfer"] = float(fd / np.hypot(fd, f_c))
    out["scale"] = _scale(f_c)
    if declared_f_c:
        out["scale_at_declared_f_c"] = _scale(float(declared_f_c))
        out["degeneracy_factor"] = float(out["scale_at_declared_f_c"] / out["scale"])
        out["declared_f_c_hz"] = float(declared_f_c)
    return out


# ------------------------------------------------------- what the record determines

def faxen_drag_ratio(lam: float) -> float:
    """Faxen's parallel-to-wall drag enhancement gamma/gamma_bulk at lam = a/h
    (bead radius over centre height). The series is a near-field expansion and
    is not used outside lam <= 0.6 by `faxen_height`."""
    return 1.0 / (1.0 - 9 / 16 * lam + 1 / 8 * lam ** 3
                  - 45 / 256 * lam ** 4 - 1 / 16 * lam ** 5)


FAXEN_LAMBDA_MAX = 0.6


def faxen_height(drag_ratio: float, radius_um: float,
                 drag_ratio_se: float | None = None) -> dict:
    """Conditional Faxen inversion and, only with uncertainty, a height reading.

    None means UNKNOWN uncertainty, not zero. `conditional_height_um` is the
    mathematical inversion, never a measured position on its own. Exact
    synthetic ratios can explicitly declare se=0. A supplied finite uncertainty
    reaching bulk gives only a lower bound; invalid uncertainty refuses.
    The five-term series is restricted to lambda <= FAXEN_LAMBDA_MAX.
    """
    from scipy.optimize import brentq
    F, radius = float(drag_ratio), float(radius_um)
    out = {"height_um": None, "lambda": None, "height_determined": False,
           "conditional_height_um": None, "evidence": "empirical",
           "uncertainty_supplied": drag_ratio_se is not None}
    if not np.isfinite(F) or not np.isfinite(radius) or radius <= 0:
        return dict(out, note="finite drag ratio and positive finite radius required")
    se = None if drag_ratio_se is None else float(drag_ratio_se)
    if se is not None and (not np.isfinite(se) or se < 0):
        return dict(out, note="drag uncertainty must be finite and nonnegative")
    max_F = faxen_drag_ratio(FAXEN_LAMBDA_MAX)

    def invert(ratio):
        if not 1 < ratio <= max_F:
            return None
        # Solve for z=lambda/delta, delta=1-1/F. (F-1)/F avoids
        # subtracting two nearly equal reciprocals. The scaled polynomial
        # stays O(1) even for nextafter(1,2), unlike the old lambda>=1e-9
        # bracket. delta/lambda >= 9/16 - .6**2/8 > 1/2, so z<2.
        delta = (ratio - 1) / ratio
        def equation(z):
            lam = delta*z
            return z*(9/16 - lam**2/8 + 45*lam**3/256 + lam**4/16) - 1
        z = brentq(equation, 0., 2., xtol=1e-14)
        lam = delta*z
        return radius/lam, lam

    central = invert(F)
    if central is not None:
        out.update(conditional_height_um=float(central[0]), **{"lambda": float(central[1])})
    if se is None:
        return dict(out, note="uncertainty not supplied: conditional series inversion only; height not determined")
    if F-se <= 1:
        high = invert(F+se)
        return dict(out, height_lower_bound_um=None if high is None else float(high[0]),
                    note="drag interval is consistent with bulk: no finite upper height bound")
    if central is None:
        return dict(out, note="drag ratio outside the Faxen series domain")
    low_height, high_height = invert(F+se), invert(F-se)
    if low_height is None or high_height is None:
        return dict(out, note="drag uncertainty extends outside the Faxen series domain")
    out.update(height_um=float(central[0]), gap_um=float(central[0]-radius),
               height_interval_um=[float(low_height[0]), float(high_height[0])],
               height_determined=True,
               note="Faxen five-term series, conditional on the supplied drag uncertainty; no added coverage claim")
    return out


def local_drag(theta: float, scale_um_per_v: float, var_volts: float,
               temperature_c: float, *, r_var: float = 1.0,
               radius_um: float | None = None,
               viscosity: float | None = None) -> dict:
    """What the record determines about kappa and gamma SEPARATELY.

    The record alone gives only the ratio: theta = kappa/gamma, scale-free. A
    displacement scale adds the missing dimension, and then equipartition and
    the timescale split the ratio:

        kappa = kB T / (Rd^2 Var_volts / r_var)        gamma = kappa / theta

    `r_var` is the record's VARIANCE retention (`Retention.variance`), not its
    quadratic-variation retention: the two differ by a factor of three on a
    real instrument. Both kappa and gamma scale as 1/Rd^2, so both inherit the
    displacement scale's units and any uncertainty in them EXACTLY -- their
    ratio does not, which is why theta is the quantity the record can state
    without a scale."""
    var_nm2 = float(var_volts) * (float(scale_um_per_v) * 1e3) ** 2 / float(r_var)
    kt = KB * (float(temperature_c) + 273.15)
    kappa = kt / (var_nm2 * 1e-18)
    gamma = kappa / float(theta)
    out = {"kappa_N_per_m": kappa, "kappa_pN_per_nm": kappa * 1e3,
           "gamma_kg_per_s": gamma, "theta_per_s": float(theta),
           "variance_nm2": var_nm2, "r_var": float(r_var)}
    if radius_um and viscosity:
        bulk = 6 * np.pi * float(viscosity) * float(radius_um) * 1e-6
        out["gamma_bulk_kg_per_s"] = bulk
        out["drag_over_bulk"] = gamma / bulk
        out["faxen"] = faxen_height(gamma / bulk, float(radius_um))
    return out


# ------------------------------------------- what moved: three scale-free ratios

def realized_diffusion(x, dt: float, theta: float, ret: Retention, *,
                       strides=None, theta_tau_max: float = 1.0) -> dict:
    """The record's OWN diffusion, in the signal's units squared per second.

    Each lag-s increment is corrected by the retention at ITS OWN weight
    (`Retention.fraction(s)`), and the Ornstein-Uhlenbeck relation
    E[(X_{t+tau} - X_t)^2] = (b^2/theta)(1 - e^{-theta tau}) is inverted. The
    answer must then be FLAT in the stride, and its flatness is the self-check:
    'step past the filter' becomes a measured correction with a residual rather
    than a hope. Strides are capped at theta*tau <= `theta_tau_max`, beyond
    which the exponential saturates and the inversion is ill-conditioned.

    No displacement scale enters, so this ratio against a calibration's own
    fitted diffusion is scale-free -- which is what makes it usable beside the
    timescale ratio to say WHICH of kappa and gamma moved."""
    x = np.asarray(x, float)
    if strides is None:
        cap = max(2, int(theta_tau_max / (float(theta) * dt)))
        strides = [s for s in (2, 3, 4, 6, 8, 12, 16, 24, 32) if s <= cap]
    if not strides:
        return {"plateau": None, "note": "no stride is both past the filter and "
                                         "inside the exponential's usable range"}
    per = {}
    for s in strides:
        m = float(np.mean((x[s:] - x[:-s]) ** 2))
        r = ret.fraction(s)
        per[int(s)] = m * float(theta) / (1 - np.exp(-float(theta) * s * dt)) / r
    vals = np.array(list(per.values()), float)
    return {"per_stride": per, "strides": list(per),
            "plateau": float(np.median(vals)),
            "spread": float(vals.max() / vals.min() - 1.0)}


ATTRIBUTION_TOL = 0.08  # historical compatibility bar, not calibrated joint coverage
# Synthetic theta precision is regime-dependent: PRECISION_AND_TEST_SELECTION.md.
CONTAMINANT_VARIANCE_FLOOR = 0.08  # separate materiality declaration, C4 control failure


def attribute_deviation(theta_ratio: float, diffusion_ratio: float,
                        variance_ratio: float, *, tol: float = ATTRIBUTION_TOL,
                        contaminant_variance_floor: float = CONTAMINANT_VARIANCE_FLOOR) -> dict:
    """WHICH single quantity moved, from three ratios that need no length scale.

    A record deviating from its calibration deviates in a signature, because the
    three observables respond differently to the candidate causes. Each ratio
    is (record / calibration's expectation):

        drag       gamma -> F gamma      (theta, D, Var) = (1/F, 1/F, 1)
        stiffness  kappa -> c kappa      (theta, D, Var) = (c,   1,   1/c)
        scale      Rd    -> s Rd         (theta, D, Var) = (1,   1/s^2, 1/s^2)
        slow-contaminant  added variance v  (free, 1, 1+v), v > materiality floor

    -- because theta = kappa/gamma, D = kB T/(gamma Rd^2) and Var = kB T/(kappa
    Rd^2), all read in the DETECTOR's units so that no length is assumed. Each
    drag/stiffness/scale hypothesis carries one free parameter fitted in log
    space. Slow contamination fits only excess variance; apparent theta is
    unpredicted and short-lag diffusion must agree independently. Its excess must
    exceed contaminant_variance_floor (C4: a 3.8% control excess otherwise
    falsely attributed). Theta is unpredicted and cannot gate admission. More
    than one passing signature refuses, including the intrinsic overlap between
    a stiffness decrease and slow contamination. A timescale for the contaminant
    must come from its covariance, not this triple.

    Measured (tweezers C3): the C-Trap's near-surface record reads
    (0.701, 0.785, 1.012) -- the drag signature at F = 1.348, and neither of the
    others explains it within tolerance."""
    r = np.array([float(theta_ratio), float(diffusion_ratio), float(variance_ratio)])
    if np.any(r <= 0) or not np.all(np.isfinite(r)):
        return {"verdict": "unreadable", "ratios": r.tolist()}
    admission = {
        "variance_excess": float(r[2] - 1),
        "variance_floor": float(contaminant_variance_floor),
        "distance_to_floor": float(r[2] - 1 - contaminant_variance_floor),
        "admitted": bool(r[2] > 1 + contaminant_variance_floor),
        "theta_used_for_admission": False,
        "note": "declared materiality boundary, no uncertainty or hysteresis model"}
    lr = np.log(r)
    if np.max(np.abs(np.expm1(np.abs(lr)))) <= tol:
        return {"verdict": "consistent", "ratios": r.tolist(), "tol": tol,
                "slow_contaminant_admission": admission,
                "note": ("all three ratios are 1 within tolerance: the record "
                         "agrees with its calibration and there is nothing to "
                         "attribute")}
    out = {}
    # each hypothesis: log-ratios = coef * log(parameter); fit the parameter by
    # least squares on the coefficient vector, then measure the residual
    for name, coef, param in (("drag", np.array([-1.0, -1.0, 0.0]), "F"),
                              ("stiffness", np.array([1.0, 0.0, -1.0]), "c"),
                              ("scale", np.array([0.0, -2.0, -2.0]), "s")):
        lp = float(coef @ lr / (coef @ coef))
        resid = lr - coef * lp
        out[name] = {"parameter": param, "value": float(np.exp(lp)),
                     "max_residual": float(np.max(np.abs(np.expm1(resid))))}
    # The apparent theta is deliberately NOT predicted by this signature. The
    # variance estimates v; short-lag diffusion is the remaining prediction.
    # It cannot identify a contaminant timescale, or distinguish a stiffness
    # change that also passes. Ambiguity keeps the existing refusal rule.
    if admission["admitted"]:
        out["slow-contaminant"] = {
            "parameter": "v", "value": float(r[2] - 1),
            "max_residual": float(abs(r[1] - 1)),
            "unpredicted": ["apparent_theta"],
            "note": "excess variance in reference-variance units; not a timescale"}
    ranked = sorted(out, key=lambda k: out[k]["max_residual"])
    best, second = ranked[0], ranked[1]
    ok = out[best]["max_residual"] <= tol and out[second]["max_residual"] > tol
    return {"verdict": best if ok else "unattributed", "hypotheses": out,
            "ranked": ranked, "tol": tol, "ratios": r.tolist(),
            "slow_contaminant_admission": admission,
            "note": (f"{best}: {out[best]['parameter']} = {out[best]['value']:.3f} "
                     f"matches its predicted ratios to {out[best]['max_residual']:.1%}, "
                     f"while {second} leaves {out[second]['max_residual']:.1%}"
                     if ok else
                     "no unique signature passes: " + ", ".join(
                         f"{k} leaves {out[k]['max_residual']:.1%}" for k in ranked))}


def covariance_timescale(x, dt: float, *, min_lag: int = 8,
                         max_lag: int = 200, block_size: int = 2048,
                         shrinkage: float = 0.2) -> dict:
    """Empirical GLS of A exp(-theta*t)+C on a fixed covariance window.

    C absorbs a locally slow covariance component; it is not a contaminant
    amplitude measurement. Block product covariance supplies correlated weights,
    regularized toward its diagonal. No reference rate or ACF sign selects lags.
    Precision must be calibrated on the acquisition regime, not assumed from
    optimizer curvature. This opt-in estimator leaves historical ACF readings
    unchanged. See FALSIFIABILITY_REGISTRATION.md.
    """
    from scipy.linalg import cholesky, solve_triangular
    from scipy.optimize import least_squares
    v = np.asarray(x, float)
    if (v.ndim != 1 or not np.all(np.isfinite(v)) or not np.isfinite(dt) or dt <= 0
            or not 0 < shrinkage <= 1 or min_lag < 1
            or max_lag <= min_lag + 3 or block_size <= max_lag):
        raise ValueError("invalid covariance estimator inputs")
    nblocks = len(v) // block_size
    metadata = {"evidence": "empirical", "n_blocks": nblocks,
                "n_lags": max_lag-min_lag, "shrinkage": shrinkage,
                "method": "fixed-window block-covariance GLS",
                "uncertainty": "requires regime calibration"}
    if nblocks < 4 or np.var(v) <= 0:
        return {**metadata, "theta": None, "amplitude_normalized": None,
                "offset_normalized": None,
                "note": "requires variance and four covariance blocks"}
    v = (v - v.mean()) / v.std()
    blocks = v[:nblocks*block_size].reshape(nblocks, block_size)
    lags = np.arange(min_lag, max_lag)
    products = np.array([np.mean(blocks[:, :-k]*blocks[:, k:], axis=1) for k in lags]).T
    cov = np.cov(products, rowvar=False) / nblocks
    diagonal = np.maximum(np.diag(cov), 1e-12)
    cov = (1-shrinkage)*cov + shrinkage*np.diag(diagonal)
    root = cholesky(cov, lower=True)
    target = products.mean(axis=0)

    def residual(p):
        model = p[0]*np.exp(-p[1]*lags)+p[2]
        return solve_triangular(root, model-target, lower=True)

    # A fixed initial dimensionless rate, not an input calibration or truth.
    fit = least_squares(residual, [1., .03, 0.],
                        bounds=([0, 1e-6, -2], [10, 2., 2]), max_nfev=1000)
    rate = float(fit.x[1])
    readable = fit.success and not np.any(fit.active_mask) and np.isfinite(rate)
    return {**metadata, "theta": rate/dt if readable else None,
            "amplitude_normalized": float(fit.x[0]),
            "offset_normalized": float(fit.x[2]), "n_blocks": nblocks,
            "n_lags": len(lags), "shrinkage": shrinkage,
            "method": "fixed-window block-covariance GLS",
            "evidence": "empirical", "uncertainty": "requires regime calibration",
            "note": "fit converged" if readable else "fit failed or reached a bound"}

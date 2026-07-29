"""The PDE verify track (registered: docs/CASE_STUDY_PDE_C2.md).

A weak-form certificate is a claim about patch integrals. This asks whether the
certified law, integrated FORWARD from an initial condition it never saw,
reproduces the field -- the strongest falsification available to the arc, and
the same shape as the RRab/Poisson verify-the-declared-form track.

Everything the forecast adds to the claim is declared, not assumed:
  * solver error -- MEASURED on a tolerance ladder (rtol, rtol/10); a run that
    does not converge on the ladder is refused, not banded,
  * parameter uncertainty -- the certificate claims an INTERVAL, so the forecast
    is the envelope over the interval endpoints, never a single trajectory,
  * the field's own declared noise sigma.
"""
from __future__ import annotations

import numpy as np
from scipy.integrate import solve_ivp


def _spectral_derivs(u, k):
    """u_x, u_xx, u_xxx on a periodic grid by FFT (exact for band-limited u)."""
    U = np.fft.rfft(u)
    return (np.fft.irfft(1j * k * U, len(u)),
            np.fft.irfft(-(k ** 2) * U, len(u)),
            np.fft.irfft(-1j * k ** 3 * U, len(u)))


def integrate(u0, x, t_eval, coeffs, *, rtol=1e-10, atol=1e-12):
    """Method of lines for u_t = c_uxx u_xx + c_uux u u_x + c_ux u_x + c_uxxx u_xxx.

    `coeffs` maps the weak-form term names to values, so the SAME vocabulary
    that was certified is what gets integrated -- no hand translation."""
    L = float(x[-1] - x[0]) + float(x[1] - x[0])
    k = 2 * np.pi * np.fft.rfftfreq(len(x), d=L / len(x))
    c_xx = coeffs.get("u_xx", 0.0)
    c_uux = coeffs.get("u*u_x", 0.0)
    c_x = coeffs.get("u_x", 0.0)
    c_xxx = coeffs.get("u_xxx", 0.0)
    c_u = coeffs.get("u", 0.0)

    def rhs(_, u):
        ux, uxx, uxxx = _spectral_derivs(u, k)
        return c_xx * uxx + c_uux * u * ux + c_x * ux + c_xxx * uxxx + c_u * u

    sol = solve_ivp(rhs, (t_eval[0], t_eval[-1]), u0, t_eval=t_eval,
                    method="RK45", rtol=rtol, atol=atol)
    return sol.y if sol.success else None


MACHINE_REL = 1e3 * np.finfo(float).eps          # the program's float floor


def solver_bound(u0, x, t_eval, coeffs, *, rtol=1e-10):
    """Declared solver error: the tolerance ladder's own disagreement, declared
    UNIFORMLY over the domain and floored at the float term.

    Measured on exact heat: the pointwise ladder difference under-covers the
    fine run's actual error at ~6% of points (by up to 1e-12) -- adaptive
    stepping makes the per-point error fluctuate while the difference does not,
    and at this level part of the residual is float roundoff the ladder cannot
    see at all. A trajectory's integration error is a global property; declaring
    it per-point from a noisy difference claims more resolution than the
    measurement has. Returns (forecast at the finer tolerance, scalar bound)."""
    a = integrate(u0, x, t_eval, coeffs, rtol=rtol, atol=rtol * 1e-2)
    b = integrate(u0, x, t_eval, coeffs, rtol=rtol * 0.1, atol=rtol * 1e-3)
    if a is None or b is None:
        return None, None
    return b, float(np.max(np.abs(a - b))) + MACHINE_REL * float(
        np.max(np.abs(b)))


def forecast_envelope(u0, x, t_eval, intervals, *, rtol=1e-10):
    """Integrate at every combination of interval ENDPOINTS plus the centre and
    return (lo, hi, centre, solver_bound). The certificate claims a family of
    laws, so the forecast is a family too."""
    names = list(intervals)
    corners = [{}]
    for n in names:
        lo, hi = intervals[n]
        mid = 0.5 * (lo + hi)
        corners = [dict(c, **{n: v}) for c in corners for v in (lo, mid, hi)]
    runs, sb = [], None
    for c in corners:
        y, bound = solver_bound(u0, x, t_eval, c, rtol=rtol)
        if y is None:
            return None
        runs.append(y)
        sb = bound if sb is None else max(sb, bound)
    R = np.stack(runs)
    centre = {n: 0.5 * (intervals[n][0] + intervals[n][1]) for n in names}
    y_c, _ = solver_bound(u0, x, t_eval, centre, rtol=rtol)
    return R.min(axis=0), R.max(axis=0), y_c, sb


def ic_noise_bound(u0, x, t_eval, coeffs, sigma, *, draws=3, rtol=1e-10,
                   seed=0):
    """How far the DECLARED measurement noise on the initial condition moves the
    forecast, measured by re-integrating from perturbed initial conditions.

    This term is not optional and was measured the hard way: a forecast started
    from measured data carries that data's noise for the whole window, and for a
    NON-DISSIPATIVE equation it never damps. Advection at sigma = 1e-5 moved the
    forecast by 2.8 sigma, so a band carrying only the target's 4 sigma failed at
    the tail of every noisy rung while the law itself was exact to 3e-11.
    Computable from measured data alone: sigma is declared, so the perturbation
    can be simulated without knowing the clean field."""
    if sigma <= 0:
        return 0.0
    base = integrate(u0, x, t_eval, coeffs, rtol=rtol, atol=rtol * 1e-2)
    if base is None:
        return None
    rng = np.random.default_rng(seed)
    worst = 0.0
    for _ in range(draws):
        pert = integrate(u0 + rng.normal(0, sigma, u0.shape), x, t_eval, coeffs,
                         rtol=rtol, atol=rtol * 1e-2)
        if pert is None:
            return None
        worst = max(worst, float(np.max(np.abs(pert - base))))
    return worst


def verify(u_true, u0, x, t_eval, intervals, *, sigma=0.0, rtol=1e-10,
           u_clean=None):
    """FORECAST-VERIFIED iff the field lies inside the declared band everywhere.

    Two claims, kept apart because they are different:
      * LAW-VERIFY (when `u_clean` is given, i.e. a dev campaign where the truth
        is known): forecast from the clean initial condition against the clean
        field, band = solver + parameter-interval envelope. This tests the LAW.
      * DATA-VERIFY: forecast from the MEASURED initial condition against the
        measured field, band additionally carrying the target's own 4-sigma and
        the propagated initial-condition noise. This tests the law as a
        predictor of measurements, which is a strictly harder claim.
    """
    got = forecast_envelope(u0, x, t_eval, intervals, rtol=rtol)
    if got is None:
        return {"verified": False, "refusal": "solver-ladder-did-not-converge"}
    lo, hi, centre, sb = got
    centre_c = {n: 0.5 * (v[0] + v[1]) for n, v in intervals.items()}
    icb = ic_noise_bound(u0, x, t_eval, centre_c, sigma, rtol=rtol)
    if icb is None:
        return {"verified": False, "refusal": "ic-noise-probe-did-not-converge"}
    band_lo = lo - sb - 4.0 * sigma - icb
    band_hi = hi + sb + 4.0 * sigma + icb
    inside = (u_true >= band_lo) & (u_true <= band_hi)
    out = {"data_verified": bool(np.all(inside)),
           "n_outside": int((~inside).sum()), "n_points": int(inside.size),
           "worst_excess": float(np.max(np.maximum(band_lo - u_true,
                                                   u_true - band_hi))),
           "envelope_width_med": float(np.median(hi - lo)),
           "solver_bound": float(sb), "noise_term": 4.0 * sigma,
           "ic_noise_bound": float(icb)}
    if u_clean is not None:
        gc = forecast_envelope(u_clean[:, 0], x, t_eval, intervals, rtol=rtol)
        if gc is None:
            out["verified"] = False
            out["refusal"] = "solver-ladder-did-not-converge (clean)"
            return out
        clo, chi, cc, csb = gc
        ins = (u_clean >= clo - csb) & (u_clean <= chi + csb)
        out["verified"] = bool(np.all(ins))              # the LAW claim
        out["law_n_outside"] = int((~ins).sum())
        out["law_max_err"] = (None if cc is None
                              else float(np.max(np.abs(cc - u_clean))))
    else:
        out["verified"] = out["data_verified"]
    return out

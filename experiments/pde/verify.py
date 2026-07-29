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


def _phi(z):
    """(e^z, phi_1(z), phi_2(z)) with the small-|z| series where the closed
    forms cancel catastrophically. phi_1 = (e^z-1)/z, phi_2 = (e^z-1-z)/z^2."""
    E = np.exp(z)
    small = np.abs(z) < 1e-3
    zs = np.where(small, 1.0, z)
    p1 = np.where(small, 1 + z / 2 + z ** 2 / 6 + z ** 3 / 24, (E - 1) / zs)
    p2 = np.where(small, 0.5 + z / 6 + z ** 2 / 24 + z ** 3 / 120,
                  (E - 1 - z) / zs ** 2)
    return E, p1, p2


def integrate(u0, x, t_eval, coeffs, *, rtol=1e-10, atol=1e-12,
              scheme="direct", nsub=64):
    """Method of lines for u_t = c_uxx u_xx + c_uux u u_x + c_ux u_x + c_uxxx u_xxx.

    `coeffs` maps the weak-form term names to values, so the SAME vocabulary
    that was certified is what gets integrated -- no hand translation.

    `scheme="direct"` is explicit RK45 in real space: what C0-C3 ran, and the
    default so those results stand. It is DIFFUSION-LIMITED -- the stable step
    goes like 1/(nu k_max^2), so on a 512-point grid at nu = 0.2 it needs ~250k
    steps per trajectory and a verify run stops being feasible. That is not a
    tolerance problem and no tolerance ladder can rescue it.

    `scheme="etd"` is an exponential integrator (ETD-RK2): every LINEAR term is
    diagonal in Fourier space and is solved EXACTLY over each substep, so only
    the nonlinearity is stepped and the diffusive stiffness is gone. Its error
    is declared the same way everything else here is -- by a ladder, in `nsub`
    rather than in `rtol` (see solver_bound). This is the one to use on
    PDEBench-scale grids.

    (An integrating factor unwrapped from t0 rather than per substep was tried
    first and is WRONG for diffusion: it carries e^{+nu k^2 t}, which overflows
    within a fraction of a time unit at any real resolution.)
    """
    L = float(x[-1] - x[0]) + float(x[1] - x[0])
    k = 2 * np.pi * np.fft.rfftfreq(len(x), d=L / len(x))
    c_xx = coeffs.get("u_xx", 0.0)
    c_uux = coeffs.get("u*u_x", 0.0)
    c_x = coeffs.get("u_x", 0.0)
    c_xxx = coeffs.get("u_xxx", 0.0)
    c_u = coeffs.get("u", 0.0)
    # POINTWISE reaction terms: a reaction-diffusion law is u_t = nu u_xx +
    # rho u (1 - u), whose quadratic part is a plain function of u, not a
    # divergence. Omitting them does not make the reference solve a slightly
    # different equation, it makes it solve a different problem (measured: the
    # deviation came back at 0.875 of the field scale, which is the missing
    # reaction, not the shipped solver's error).
    c_u2 = coeffs.get("u^2", 0.0)
    c_u3 = coeffs.get("u^3", 0.0)

    if scheme == "direct":
        def rhs(_, u):
            ux, uxx, uxxx = _spectral_derivs(u, k)
            return (c_xx * uxx + c_uux * u * ux + c_x * ux + c_xxx * uxxx
                    + c_u * u + c_u2 * u ** 2 + c_u3 * u ** 3)

        sol = solve_ivp(rhs, (t_eval[0], t_eval[-1]), u0, t_eval=t_eval,
                        method="RK45", rtol=rtol, atol=atol)
        return sol.y if sol.success else None

    # d_x -> i k, so the linear symbol is exact per mode
    lam = (-c_xx * k ** 2 + 1j * c_x * k - 1j * c_xxx * k ** 3 + c_u)
    n = len(x)
    ik = 1j * k
    if n % 2 == 0:
        # the Nyquist mode of a real field has no odd-derivative image: keeping
        # one would rotate its phase for no physical reason
        ik = ik.copy()
        ik[-1] = 0.0
        lam = lam.copy()
        lam[-1] = -c_xx * k[-1] ** 2 + c_u

    def nonlinear(U):
        if not (c_uux or c_u2 or c_u3):
            return np.zeros_like(U)
        u = np.fft.irfft(U, n)
        out = np.zeros_like(U)
        if c_uux:
            out = out + c_uux * ik * np.fft.rfft(0.5 * u * u)
        if c_u2 or c_u3:
            out = out + np.fft.rfft(c_u2 * u ** 2 + c_u3 * u ** 3)
        return out

    U = np.fft.rfft(np.asarray(u0, float))
    out = np.empty((n, len(t_eval)))
    out[:, 0] = np.fft.irfft(U, n)
    for j in range(1, len(t_eval)):
        h = (float(t_eval[j]) - float(t_eval[j - 1])) / nsub
        E, p1, p2 = _phi(lam * h)
        for _ in range(nsub):
            N0 = nonlinear(U)
            a = E * U + h * p1 * N0
            U = a + h * p2 * (nonlinear(a) - N0)      # ETD-RK2
        out[:, j] = np.fft.irfft(U, n)
    return out


MACHINE_REL = 1e3 * np.finfo(float).eps          # the program's float floor


def solver_bound(u0, x, t_eval, coeffs, *, rtol=1e-10, scheme="direct",
                 nsub=64):
    """Declared solver error: the tolerance ladder's own disagreement, declared
    UNIFORMLY over the domain and floored at the float term.

    Measured on exact heat: the pointwise ladder difference under-covers the
    fine run's actual error at ~6% of points (by up to 1e-12) -- adaptive
    stepping makes the per-point error fluctuate while the difference does not,
    and at this level part of the residual is float roundoff the ladder cannot
    see at all. A trajectory's integration error is a global property; declaring
    it per-point from a noisy difference claims more resolution than the
    measurement has. Returns (forecast at the finer tolerance, scalar bound)."""
    if scheme == "etd":
        # the ladder is in SUBSTEPS for an exponential integrator: its accuracy
        # knob is the step, not a solver tolerance
        a = integrate(u0, x, t_eval, coeffs, scheme=scheme, nsub=nsub)
        b = integrate(u0, x, t_eval, coeffs, scheme=scheme, nsub=2 * nsub)
    else:
        a = integrate(u0, x, t_eval, coeffs, rtol=rtol, atol=rtol * 1e-2,
                      scheme=scheme)
        b = integrate(u0, x, t_eval, coeffs, rtol=rtol * 0.1, atol=rtol * 1e-3,
                      scheme=scheme)
    if a is None or b is None:
        return None, None
    return b, float(np.max(np.abs(a - b))) + MACHINE_REL * float(
        np.max(np.abs(b)))


def forecast_envelope(u0, x, t_eval, intervals, *, rtol=1e-10,
                      scheme="direct", nsub=64, n_samples=3):
    """Integrate at every combination of interval ENDPOINTS plus the centre and
    return (lo, hi, centre, solver_bound). The certificate claims a family of
    laws, so the forecast is a family too."""
    names = list(intervals)
    corners = [{}]
    for n in names:
        lo, hi = intervals[n]
        # SAMPLE the interval, do not just take its ends. For a translation the
        # forecast is not monotone in the parameter, so the pointwise min/max
        # over {lo, mid, hi} is NOT an envelope of the family: measured on
        # PDEBench advection, where a wide beta interval left 45309 of 205824
        # points outside a "envelope" built from three corners while every
        # member of the family was fine. `n_samples` is reported.
        vs = np.linspace(lo, hi, max(3, n_samples)) if hi > lo else [lo]
        corners = [dict(c, **{n: float(v)}) for c in corners for v in vs]
    runs, sb = [], None
    for c in corners:
        y, bound = solver_bound(u0, x, t_eval, c, rtol=rtol, scheme=scheme,
                                nsub=nsub)
        if y is None:
            return None
        runs.append(y)
        sb = bound if sb is None else max(sb, bound)
    R = np.stack(runs)
    centre = {n: 0.5 * (intervals[n][0] + intervals[n][1]) for n in names}
    y_c, _ = solver_bound(u0, x, t_eval, centre, rtol=rtol, scheme=scheme,
                          nsub=nsub)
    return R.min(axis=0), R.max(axis=0), y_c, sb


def ic_noise_bound(u0, x, t_eval, coeffs, sigma, *, draws=3, rtol=1e-10,
                   seed=0, scheme="direct", nsub=64):
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
    base = integrate(u0, x, t_eval, coeffs, rtol=rtol, atol=rtol * 1e-2,
                     scheme=scheme, nsub=nsub)
    if base is None:
        return None
    rng = np.random.default_rng(seed)
    worst = 0.0
    for _ in range(draws):
        pert = integrate(u0 + rng.normal(0, sigma, u0.shape), x, t_eval, coeffs,
                         rtol=rtol, atol=rtol * 1e-2, scheme=scheme, nsub=nsub)
        if pert is None:
            return None
        worst = max(worst, float(np.max(np.abs(pert - base))))
    return worst


PERIODIC_SEAM_MAX = 2.0


def periodicity_seam(u) -> float:
    """|u(x_first) - u(x_last)| measured in units of an ORDINARY interior step.

    The raw wrap gap does not discriminate: on an endpoint-excluded grid the
    first and last cells are one dx apart, so a perfectly periodic field has a
    gap of dx*|u_x| -- measured at 2.6e-2 to 3.0e-2 on this program's own
    periodic C2 fields, which is larger than PDEBench advection's. Dividing by
    the field's own 99th-percentile interior step makes the test scale-free and
    shock-tolerant: a periodic field's seam is an ordinary step.

    MEASURED separation: every periodic field tested sits at 0.14-0.91 (our C2
    heat and Burgers, PDEBench advection, Burgers and periodic CFD) while
    PDEBench's TRANSMISSIVE-boundary CFD sits at 4.63.
    """
    u = np.asarray(u, float)
    wrap = np.abs(u[0] - u[-1])
    step = np.abs(np.diff(u, axis=0))
    return float(np.median(wrap) / max(np.quantile(step, 0.99), 1e-300))


def _dx(f, k, order):
    """d^order f / dx^order on a periodic grid by FFT."""
    if order == 0:
        return f
    F = np.fft.rfft(f)
    m = (1j * k) ** order
    if order % 2 == 1 and len(f) % 2 == 0:
        m = m.copy()
        m[-1] = 0.0             # the Nyquist mode has no odd-derivative image
    return np.fft.irfft(m * F, len(f))


def integrate_system(q0, x, t_eval, equations, terms, to_fields, *,
                     rtol=1e-10, atol=None):
    """Forward-integrate a SYSTEM in the weak form's own vocabulary.

    `equations` maps each target term name to {feature term name: coefficient}
    -- exactly what a system certificate reports. Every library term is
    d^alpha(g(fields)), so the evolution of the target's own g is
    d/dt g_target = sum_k c_k d_x^{alpha_k} g_k(fields), evaluated spectrally.
    Nothing is hand-translated: the vocabulary that was certified is the
    vocabulary that gets integrated.

    `to_fields(q)` maps the evolved quantities back to the physical fields --
    a DECLARATION the caller makes (shallow water evolves (h, hu) while the
    fields are (h, u)), not something this function can guess.
    """
    n = len(x)
    Lx = float(x[-1] - x[0]) + float(x[1] - x[0])
    k = 2 * np.pi * np.fft.rfftfreq(n, d=Lx / n)
    tgts = list(equations)

    def rhs(_, q):
        Q = q.reshape(len(tgts), n)
        fields = to_fields(Q)
        out = np.zeros_like(Q)
        for i, tg in enumerate(tgts):
            acc = np.zeros(n)
            for name, c in equations[tg].items():
                if c == 0.0:
                    continue
                tm = terms[name]
                g = tm.g(fields) if name != "1" else np.ones(n)
                acc = acc + c * _dx(np.asarray(g, float).ravel(), k, tm.ax)
            out[i] = acc
        return out.ravel()

    s = solve_ivp(rhs, (t_eval[0], t_eval[-1]), np.asarray(q0, float).ravel(),
                  t_eval=t_eval, method="RK45", rtol=rtol,
                  atol=rtol * 1e-2 if atol is None else atol)
    return s.y.reshape(len(tgts), n, len(t_eval)) if s.success else None


def system_solver_bound(q0, x, t_eval, equations, terms, to_fields, *,
                        rtol=1e-10):
    """The declared solver error for a system forecast: the tolerance ladder's
    own disagreement, uniform over the domain, floored at the float term."""
    a = integrate_system(q0, x, t_eval, equations, terms, to_fields, rtol=rtol)
    b = integrate_system(q0, x, t_eval, equations, terms, to_fields,
                         rtol=rtol * 0.1)
    if a is None or b is None:
        return None, None
    return b, float(np.max(np.abs(a - b))) + MACHINE_REL * float(
        np.max(np.abs(b)))


def verify_system(fields_true, q0, x, t_eval, intervals, terms, to_fields,
                  field_of=None, *, sigma=0.0, rtol=1e-10, fields_clean=None):
    """FORECAST-VERIFY a certified SYSTEM from an initial condition it never saw.

    `intervals` maps target -> {term: (lo, hi)}: the certificate claims a family
    of systems, so the forecast is a family too. The envelope is built by
    perturbing ONE parameter at a time to each end of its interval and SUMMING
    the resulting deviations -- the triangle inequality, which upper-bounds the
    corner envelope exactly for a linear-in-parameters response and conservatively
    otherwise. That is a deliberate choice over the scalar track's full corner
    product: a 2-equation system with 3 parameters each has 3^6 = 729 corners,
    and a bound that costs 2P integrations and can only be too WIDE is the right
    trade. Reported as `envelope_method`, never silently substituted.

    `field_of` maps a target term name to the physical field it evolves, so the
    forecast can be compared against the measured fields.
    """
    centre = {tg: {n: 0.5 * (lo + hi) for n, (lo, hi) in iv.items()}
              for tg, iv in intervals.items()}
    tgts = list(intervals)

    def forecast(q_start):
        """(centre trajectory, solver bound, interval envelope) from one start."""
        yc, s_b = system_solver_bound(q_start, x, t_eval, centre, terms,
                                      to_fields, rtol=rtol)
        if yc is None:
            return None
        d = np.zeros_like(yc)
        runs = 0
        for tg, iv in intervals.items():
            for name, (lo, hi) in iv.items():
                if hi <= lo:
                    continue
                worst = np.zeros_like(yc)
                for v in (lo, hi):
                    pert = {t2: dict(c2) for t2, c2 in centre.items()}
                    pert[tg][name] = v
                    yp = integrate_system(q_start, x, t_eval, pert, terms,
                                          to_fields, rtol=rtol)
                    runs += 1
                    if yp is None:
                        return None
                    worst = np.maximum(worst, np.abs(yp - yc))
                d = d + worst
        return yc, s_b, d, runs

    got = forecast(q0)
    if got is None:
        return {"verified": False, "refusal": "solver-ladder-did-not-converge"}
    y_c, sb, dev, n_runs = got
    out = {"envelope_method": "per-parameter triangle bound",
           "n_endpoint_runs": n_runs, "solver_bound": float(sb),
           "envelope_width_med": float(np.median(2 * dev)),
           "noise_term": 4.0 * sigma, "targets": tgts}

    # A forecast started from MEASURED data carries that data's noise for the
    # whole window -- the term C2 had to add after advection failed at the tail
    # of every noisy rung with the law exact to 3e-11. Measured the same way
    # here: re-integrate from perturbed initial states.
    icb = 0.0
    if sigma > 0:
        rng = np.random.default_rng(0)
        for _ in range(3):
            q0p = np.asarray(q0, float) + rng.normal(0, sigma, np.shape(q0))
            yp = integrate_system(q0p, x, t_eval, centre, terms, to_fields,
                                  rtol=rtol)
            if yp is None:
                return {"verified": False,
                        "refusal": "ic-noise-probe-did-not-converge"}
            icb = max(icb, float(np.max(np.abs(yp - y_c))))
    out["ic_noise_bound"] = icb

    def _compare(truth_fields, band_extra, yc, s_b, d):
        n_out, worst, npts, per = 0, -np.inf, 0, {}
        for i, tg in enumerate(tgts):
            g = np.asarray(terms[tg].g(truth_fields), float)   # the evolved q
            lo = yc[i] - d[i] - s_b - band_extra
            hi = yc[i] + d[i] + s_b + band_extra
            inside = (g >= lo) & (g <= hi)
            n_out += int((~inside).sum())
            npts += int(inside.size)
            worst = max(worst, float(np.max(np.maximum(lo - g, g - hi))))
            per[tg] = int((~inside).sum())
        return n_out, worst, npts, per

    if field_of:
        out["fields_evolved"] = {tg: field_of.get(tg) for tg in tgts}
    if fields_clean is not None:
        # LAW-VERIFY forecasts from the CLEAN initial state against the clean
        # fields: it tests the LAW, so it must not be charged for the initial
        # condition's noise -- and must not be credited with a band for it
        # either. Reusing the measured start here (the first version of this
        # function) makes the law claim inherit a noise it does not band, which
        # is the same asymmetry C2 found in the scalar track.
        q0c = np.stack([np.asarray(terms[tg].g(fields_clean), float)[..., 0]
                        for tg in tgts])
        gc = forecast(q0c)
        if gc is None:
            out["verified"] = False
            out["refusal"] = "solver-ladder-did-not-converge (clean)"
            return out
        yc, sbc, devc, _ = gc
        n_out, worst, npts, per = _compare(fields_clean, 0.0, yc, sbc, devc)
        out["verified"] = bool(n_out == 0)
        out["law_n_outside"] = n_out
        out["law_worst_excess"] = worst
        out["law_outside_per_target"] = per
        out["n_points"] = npts
    # the measured comparison carries the field's own declared noise, propagated
    # through the target's g by its sensitivity at first order (g may be a
    # product of fields: (h u) moves by ~ |u| dh + |h| du), plus the IC term
    n_out, worst, npts, per = _compare(
        fields_true, 4.0 * sigma * _g_sensitivity(terms, tgts, fields_true)
        + icb, y_c, sb, dev)
    out["data_verified"] = bool(n_out == 0)
    out["n_outside"] = n_out
    out["worst_excess"] = worst
    out["outside_per_target"] = per
    out["n_points"] = npts
    if fields_clean is None:
        out["verified"] = out["data_verified"]
    return out


def _g_sensitivity(terms, tgts, fields):
    """max_i sum_f |dg/df| over the targets: how much a per-field noise of size
    sigma can move the evolved quantity (first order, declared)."""
    s = 0.0
    for tg in tgts:
        tm = terms[tg]
        for f in fields:
            s = max(s, float(np.max(np.abs(tm.dg(fields, f)))))
    return max(s, 1.0)


def verify_state(u_true, x, t_eval, modes, basis_fns, labels, law, *,
                 sigma=0.0, rtol=1e-10, u_clean=None):
    """Forward-integrate the RECOVERED INITIAL CONDITION under the known law and
    ask whether it reproduces the observed field.

    This is the state certificate's falsification track, and it is nearly free:
    the law is known, so the only new error is the integrator's, and the
    certificate's per-mode intervals propagate the same way parameter intervals
    do in the law track -- one perturbation per mode, summed (the triangle
    bound).

    Modes reported UNDETERMINED enter the forecast at their fitted value with NO
    claimed bound, and are named in the result: the certificate does not claim
    them, so a reader must not read the forecast as evidence about them.
    """
    x = np.asarray(x, float)
    centre = {lab: (0.5 * (modes[lab][0] + modes[lab][1])
                    if isinstance(modes[lab], (tuple, list)) else modes[lab])
              for lab in modes}
    idx = {lab: i for i, lab in enumerate(labels)}

    def build_u0(vals):
        u0 = np.zeros_like(x)
        for lab, v in vals.items():
            u0 = u0 + v * np.asarray(basis_fns[idx[lab]](x), float)
        return u0

    u0_c = build_u0(centre)
    y_c, sb = solver_bound(u0_c, x, t_eval, law, rtol=rtol)
    if y_c is None:
        return {"verified": False, "refusal": "solver-ladder-did-not-converge"}
    dev = np.zeros_like(y_c)
    for lab, iv in modes.items():
        if not isinstance(iv, (tuple, list)) or iv[1] <= iv[0]:
            continue
        worst = np.zeros_like(y_c)
        for v in iv:
            vals = dict(centre)
            vals[lab] = v
            yp = integrate(build_u0(vals), x, t_eval, law, rtol=rtol,
                           atol=rtol * 1e-2)
            if yp is None:
                return {"verified": False,
                        "refusal": "interval-endpoint-did-not-converge"}
            worst = np.maximum(worst, np.abs(yp - y_c))
        dev = dev + worst
    out = {"solver_bound": float(sb), "noise_term": 4.0 * sigma,
           "envelope_width_med": float(np.median(2 * dev)),
           "envelope_method": "per-mode triangle bound",
           "n_modes_in_forecast": len(centre)}
    lo, hi = y_c - dev - sb, y_c + dev + sb
    if u_clean is not None:
        inside = (u_clean >= lo) & (u_clean <= hi)
        out["verified"] = bool(np.all(inside))
        out["law_n_outside"] = int((~inside).sum())
    inside = (u_true >= lo - 4.0 * sigma) & (u_true <= hi + 4.0 * sigma)
    out["data_verified"] = bool(np.all(inside))
    out["n_outside"] = int((~inside).sum())
    out["n_points"] = int(inside.size)
    out["worst_excess"] = float(np.max(np.maximum(lo - 4.0 * sigma - u_true,
                                                  u_true - hi - 4.0 * sigma)))
    if u_clean is None:
        out["verified"] = out["data_verified"]
    return out


def verify(u_true, u0, x, t_eval, intervals, *, sigma=0.0, rtol=1e-10,
           u_clean=None, scheme="direct", nsub=64, field_err=0.0,
           n_samples=3, periodic_tol=PERIODIC_SEAM_MAX):
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
    # DOMAIN GUARD. Every forecast here differentiates spectrally, which assumes
    # periodicity; on a non-periodic field the FFT sees a jump at the seam and
    # rings. Measured on PDEBench's transmissive-boundary CFD: the interior
    # derivative error tracks the seam (8e-2), and NOTHING in the pipeline
    # noticed -- check_geometry returned ok with no notes. A capability that
    # applies itself outside its domain without saying so is the failure this
    # program exists to prevent, so the track refuses instead of forecasting.
    seam = periodicity_seam(u_true if u_clean is None else u_clean)
    if seam > periodic_tol:
        return {"verified": False, "data_verified": False,
                "refusal": "field is not periodic",
                "periodicity_seam": seam, "periodic_tol": periodic_tol,
                "note": ("the spectral forecast assumes a periodic domain; this "
                         "field's wrap seam is " f"{seam:.2f}x an ordinary "
                         "interior step. Weak-form CERTIFICATION is unaffected "
                         "(its test functions vanish inside the domain) -- it is "
                         "the forecast track that does not apply here")}
    got = forecast_envelope(u0, x, t_eval, intervals, rtol=rtol,
                            scheme=scheme, nsub=nsub, n_samples=n_samples)
    if got is None:
        return {"verified": False, "refusal": "solver-ladder-did-not-converge"}
    lo, hi, centre, sb = got
    centre_c = {n: 0.5 * (v[0] + v[1]) for n, v in intervals.items()}
    icb = ic_noise_bound(u0, x, t_eval, centre_c, sigma, rtol=rtol,
                         scheme=scheme, nsub=nsub)
    if icb is None:
        return {"verified": False, "refusal": "ic-noise-probe-did-not-converge"}
    # A DECLARED field error belongs in the forecast band as well as in the
    # certification band: the measured field is only that accurate, so a
    # forecast compared against it cannot be held to a tighter standard.
    # Coefficient 1, not 4 -- it is a computed/declared bound, not a stochastic
    # scale (certify.epsilon's `hard` channel convention). Measured on PDEBench
    # advection, where omitting it failed 58736 of 205824 points while the
    # certified law was right to 3e-5.
    band_lo = lo - sb - 4.0 * sigma - icb - field_err
    band_hi = hi + sb + 4.0 * sigma + icb + field_err
    inside = (u_true >= band_lo) & (u_true <= band_hi)
    out = {"envelope_samples_per_parameter": int(max(3, n_samples)),
           "data_verified": bool(np.all(inside)),
           "n_outside": int((~inside).sum()), "n_points": int(inside.size),
           "worst_excess": float(np.max(np.maximum(band_lo - u_true,
                                                   u_true - band_hi))),
           "envelope_width_med": float(np.median(hi - lo)),
           "solver_bound": float(sb), "noise_term": 4.0 * sigma,
           "ic_noise_bound": float(icb)}
    if u_clean is not None:
        gc = forecast_envelope(u_clean[:, 0], x, t_eval, intervals, rtol=rtol,
                               scheme=scheme, nsub=nsub, n_samples=n_samples)
        if gc is None:
            out["verified"] = False
            out["refusal"] = "solver-ladder-did-not-converge (clean)"
            return out
        clo, chi, cc, csb = gc
        ins = ((u_clean >= clo - csb - field_err)
               & (u_clean <= chi + csb + field_err))
        out["verified"] = bool(np.all(ins))              # the LAW claim
        out["law_n_outside"] = int((~ins).sum())
        out["law_max_err"] = (None if cc is None
                              else float(np.max(np.abs(cc - u_clean))))
    else:
        out["verified"] = out["data_verified"]
    return out

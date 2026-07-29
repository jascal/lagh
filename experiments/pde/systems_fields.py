"""Fields for the PDE system curriculum (registration docs/CASE_STUDY_PDE_C3.md).

Stage 1 is EXACTLY solvable, so a stage-1 failure is the instrument's and never
an integrator's -- the same reason C0 was entirely analytic. Every later stage
needs a reference solver, and a solver's error is a new undeclared quantity
entering the claim, so each solved field arrives with a MEASURED bound:

    declared field error = (tolerance-ladder disagreement)
                         + (spatial-resolution-ladder disagreement)

both measured, both uniform over the field, refused rather than banded when the
ladder does not converge. Downstream it enters the weak-form band through the L1
sensitivity (weakform.WeakSystem.det), because a solver error is one fixed
function and nothing cancels across it.
"""
from __future__ import annotations

import numpy as np
from scipy.integrate import solve_ivp
from scipy.linalg import expm

L = 2 * np.pi


# --------------------------------------------------------------------------
# Stage 1 -- linear coupled, exact (2x2 propagator per Fourier mode)
# --------------------------------------------------------------------------

def linear_pair(seed, *, a=0.1, b=0.5, c=0.05, d=-0.3, nx=257, nt=81, tmax=1.0,
                modes=(1, 2, 3), periodic=False):
    """Exact (u, v) for  u_t = a u_xx + b v,  v_t = c v_xx + d u.

    ONE phase per mode, shared by both fields: the 2x2 propagator acts on the
    amplitude pair of the SAME basis function, so independent phases would make
    the coupling term b*v spatially orthogonal to u and the system would not hold
    at all (measured in the scoping probe: the truth then misses its own band by
    0.85, which looks exactly like an identifiability finding until you check the
    truth against its own band first).
    """
    # `periodic` drops the duplicated right endpoint: the weak-form patches do
    # not care, but the spectral verify track does -- an endpoint-duplicated
    # grid makes every FFT wavenumber slightly wrong and the forecast fails for
    # a reason that has nothing to do with the law (measured: 18244 of 41634
    # points outside the band at sigma = 0, with the law exact).
    x = np.linspace(0.0, L, nx, endpoint=not periodic)
    t = np.linspace(0.0, tmax, nt)
    rng = np.random.default_rng(seed)
    amps = rng.uniform(0.2, 1.0, (2, len(modes)))
    ph = rng.uniform(0, 2 * np.pi, len(modes))
    u = np.zeros((nx, nt))
    v = np.zeros((nx, nt))
    for j, k in enumerate(modes):
        M = np.array([[-a * k ** 2, b], [d, -c * k ** 2]])
        for i, ti in enumerate(t):
            at = expm(M * ti) @ np.array([amps[0, j], amps[1, j]])
            u[:, i] += at[0] * np.cos(k * x + ph[j])
            v[:, i] += at[1] * np.cos(k * x + ph[j])
    return {"u": u, "v": v}, (x, t)


# --------------------------------------------------------------------------
# The reference solver + its declared error
# --------------------------------------------------------------------------

def _spectral(nx, Lx=L):
    return 2 * np.pi * np.fft.rfftfreq(nx, d=Lx / nx)


def _solve(rhs, y0, t, rtol):
    s = solve_ivp(rhs, (t[0], t[-1]), y0, t_eval=t, method="RK45", rtol=rtol,
                  atol=rtol * 1e-2)
    return s.y if s.success else None


def solve_declared(make_rhs, ic, x, t, *, rtol=1e-11, n_fields=2):
    """Integrate, and MEASURE the error of the field that comes out.

    Two ladders, because there are two discretizations: the time tolerance
    (rtol vs rtol/10) and the spatial resolution (nx vs 2nx, compared on the
    coarse grid). The declared bound is their sum, uniform over the field -- a
    trajectory's numerical error is a global property, and declaring it
    per-point from a fluctuating difference claims more resolution than the
    measurement has (the C2 lesson). Returns (fields, err) or (None, None).
    """
    nx = len(x)
    y0 = np.concatenate(ic(x))
    a = _solve(make_rhs(x), y0, t, rtol)
    b = _solve(make_rhs(x), y0, t, rtol * 0.1)
    if a is None or b is None:
        return None, None
    x2 = np.linspace(x[0], x[0] + L, 2 * nx, endpoint=False)
    c = _solve(make_rhs(x2), np.concatenate(ic(x2)), t, rtol * 0.1)
    if c is None:
        return None, None
    A = a.reshape(n_fields, nx, len(t))
    B = b.reshape(n_fields, nx, len(t))
    C = c.reshape(n_fields, 2 * nx, len(t))[:, ::2, :]
    err = float(np.max(np.abs(A - B))) + float(np.max(np.abs(B - C)))
    return B, err


# --------------------------------------------------------------------------
# Stage 2 -- reaction-diffusion
# --------------------------------------------------------------------------

BRU_A, BRU_B, BRU_DU, BRU_DV = 1.0, 1.7, 0.02, 0.01


def _modes(rng, n_modes):
    """n distinct wavenumbers, amplitudes and phases -- the spectral richness of
    the initial data, which Y4 says the nonlinear stages need more of."""
    ks = rng.choice([1, 2, 3, 4], size=n_modes, replace=False)
    return ks, rng.uniform(0, 2 * np.pi, n_modes)


def brusselator(seed, *, nx=256, nt=41, tmax=1.0, rtol=1e-11, n_modes=2):
    """u_t = Du u_xx + A - (B+1) u + u^2 v,  v_t = Dv v_xx + B u - u^2 v.

    The first CROSS term (u^2 v): a feature neither field can express alone, and
    the reason the multi-field library exists at all. B < 1 + A^2 keeps the
    homogeneous state stable, so the fields stay smooth and well resolved over
    the window instead of running to a limit cycle the grid cannot follow.
    """
    x = np.linspace(0.0, L, nx, endpoint=False)
    t = np.linspace(0.0, tmax, nt)
    rng = np.random.default_rng(seed)
    ks, phs = _modes(rng, n_modes)
    amps = rng.uniform(0.05, 0.2, (2, n_modes))

    def ic(xx):
        u0 = BRU_A + sum(amps[0, j] * np.cos(ks[j] * xx + phs[j])
                         for j in range(n_modes))
        v0 = BRU_B / BRU_A + sum(amps[1, j] * np.cos(ks[j] * xx - phs[j])
                                 for j in range(n_modes))
        return u0, v0

    def make_rhs(xx):
        k = _spectral(len(xx))

        def rhs(_, y):
            u, v = y[:len(xx)], y[len(xx):]
            uxx = np.fft.irfft(-(k ** 2) * np.fft.rfft(u), len(xx))
            vxx = np.fft.irfft(-(k ** 2) * np.fft.rfft(v), len(xx))
            r = u * u * v
            return np.concatenate([BRU_DU * uxx + BRU_A - (BRU_B + 1) * u + r,
                                   BRU_DV * vxx + BRU_B * u - r])
        return rhs

    y, err = solve_declared(make_rhs, ic, x, t, rtol=rtol)
    if y is None:
        return None, None, None
    return {"u": y[0], "v": y[1]}, (x, t), err


FHN_D, FHN_EPS, FHN_A, FHN_B, FHN_I = 0.05, 0.3, 0.7, 0.8, 0.5


def fitzhugh_nagumo(seed, *, nx=256, nt=41, tmax=1.0, rtol=1e-11, n_modes=2):
    """u_t = D u_xx + u - u^3/3 - v + I,  v_t = eps (u + a - b v).

    Nonlinear in ONE field (the cubic) with linear coupling: the stage that
    separates 'nonlinear' from 'cross-coupled', so a stage-2 failure can be
    attributed to the right one.
    """
    x = np.linspace(0.0, L, nx, endpoint=False)
    t = np.linspace(0.0, tmax, nt)
    rng = np.random.default_rng(seed)
    ks, phs = _modes(rng, n_modes)
    amps = rng.uniform(0.2, 0.6, (2, n_modes))

    def ic(xx):
        return (sum(amps[0, j] * np.cos(ks[j] * xx + phs[j])
                    for j in range(n_modes)),
                sum(amps[1, j] * np.cos(ks[j] * xx - phs[j])
                    for j in range(n_modes)))

    def make_rhs(xx):
        k = _spectral(len(xx))

        def rhs(_, y):
            u, v = y[:len(xx)], y[len(xx):]
            uxx = np.fft.irfft(-(k ** 2) * np.fft.rfft(u), len(xx))
            return np.concatenate([
                FHN_D * uxx + u - u ** 3 / 3.0 - v + FHN_I,
                FHN_EPS * (u + FHN_A - FHN_B * v)])
        return rhs

    y, err = solve_declared(make_rhs, ic, x, t, rtol=rtol)
    if y is None:
        return None, None, None
    return {"u": y[0], "v": y[1]}, (x, t), err


# --------------------------------------------------------------------------
# Stage 3 -- shallow water (a conservation-law system with a real flux)
# --------------------------------------------------------------------------

SW_G = 9.81


def shallow_water(seed, *, nx=256, nt=41, tmax=0.5, rtol=1e-11, amp=0.05,
                  n_modes=2):
    """h_t + (h u)_x = 0,  (h u)_t + (h u^2 + g h^2 / 2)_x = 0.

    Solved in CONSERVATIVE variables (h, m = h u) so the conservation form holds
    to solver error rather than to a discretization's idea of it; the fields
    handed to the factory are the primitive (h, u), which is what a measurement
    would give. Amplitudes stay small so nothing steepens inside the window --
    shock formation is a declared out-of-scope, not an accident.
    """
    x = np.linspace(0.0, L, nx, endpoint=False)
    t = np.linspace(0.0, tmax, nt)
    rng = np.random.default_rng(seed)
    ks, phs = _modes(rng, n_modes)
    a2 = rng.uniform(0.5, 1.0, (2, n_modes)) * amp

    def ic(xx):
        h0 = 1.0 + sum(a2[0, j] * np.cos(ks[j] * xx + phs[j])
                       for j in range(n_modes))
        u0 = sum(a2[1, j] * np.cos(ks[j] * xx - phs[j])
                 for j in range(n_modes))
        return h0, h0 * u0

    def make_rhs(xx):
        k = _spectral(len(xx))

        def d1(f):
            F = np.fft.rfft(f)
            m = 1j * k
            if len(xx) % 2 == 0:
                m = m.copy()
                m[-1] = 0.0
            return np.fft.irfft(m * F, len(xx))

        def rhs(_, y):
            h, m = y[:len(xx)], y[len(xx):]
            u = m / h
            return np.concatenate([-d1(m),
                                   -d1(m * u + 0.5 * SW_G * h ** 2)])
        return rhs

    y, err = solve_declared(make_rhs, ic, x, t, rtol=rtol)
    if y is None:
        return None, None, None
    h, m = y[0], y[1]
    return {"h": h, "u": m / h}, (x, t), err


# --------------------------------------------------------------------------
# Stage 4 -- Navier-Stokes, vorticity-streamfunction (2-D space)
# --------------------------------------------------------------------------

NS_NU = 0.02


def _ns_rhs_factory(n, nu):
    kx = 2 * np.pi * np.fft.fftfreq(n, d=L / n)
    ky = 2 * np.pi * np.fft.fftfreq(n, d=L / n)
    KX, KY = np.meshgrid(kx, ky, indexing="ij")
    K2 = KX ** 2 + KY ** 2
    K2i = np.where(K2 == 0, 1.0, K2)
    mask = (np.abs(KX) < (2 / 3) * np.abs(kx).max()) & \
           (np.abs(KY) < (2 / 3) * np.abs(ky).max())      # 2/3 dealiasing

    def parts(W):
        P = W / K2i                       # psi_hat, with the mean mode killed
        P[0, 0] = 0.0
        U = np.real(np.fft.ifft2(1j * KY * P))
        V = np.real(np.fft.ifft2(-1j * KX * P))
        return P, U, V

    def rhs(_, w):
        w = w.reshape(n, n)
        W = np.fft.fft2(w)
        _, U, V = parts(W)
        adv = np.fft.fft2(U * w) * 1j * KX + np.fft.fft2(V * w) * 1j * KY
        return np.real(np.fft.ifft2(mask * (-adv - nu * K2 * W))).ravel()

    return rhs, parts, (KX, KY, K2i)


def ns_vorticity(seed, *, n=64, nt=33, tmax=0.5, rtol=1e-10, nu=NS_NU):
    """omega_t = -(u omega)_x - (v omega)_y + nu (omega_xx + omega_yy), with
    u = psi_y, v = -psi_x and the MACHINE-EXACT input constraint
    psi_xx + psi_yy = -omega.

    Four fields are handed to the factory (omega, psi, u, v) precisely because
    three of them are related to the fourth by exact linear relations: in the
    weak-form column space the `u` column EQUALS the `psi_y` column, which is
    the constrained-input situation the engine closed for algebraic constraints
    (CASE_STUDY_GAIA_P3.md) appearing here as physics rather than as an artifact.

    Declared error: the time-tolerance ladder only. A resolution ladder in 2-D
    costs 4x the state and is reported separately by the runner's spectral-tail
    check instead -- stated here rather than folded in silently.
    """
    x = np.linspace(0.0, L, n, endpoint=False)
    y = np.linspace(0.0, L, n, endpoint=False)
    t = np.linspace(0.0, tmax, nt)
    rng = np.random.default_rng(seed)
    X, Y = np.meshgrid(x, y, indexing="ij")
    w0 = np.zeros((n, n))
    for _ in range(3):
        kx, ky = rng.choice([1, 2, 3], size=2)
        a = rng.uniform(0.4, 1.0)
        ph = rng.uniform(0, 2 * np.pi, 2)
        w0 += a * np.sin(kx * X + ph[0]) * np.cos(ky * Y + ph[1])
    rhs, parts, _ = _ns_rhs_factory(n, nu)
    out = []
    for rt in (rtol, rtol * 0.1):
        s = solve_ivp(rhs, (t[0], t[-1]), w0.ravel(), t_eval=t, method="RK45",
                      rtol=rt, atol=rt * 1e-2)
        if not s.success:
            return None, None, None
        out.append(s.y.reshape(n, n, len(t)))
    err = float(np.max(np.abs(out[0] - out[1])))
    w = out[1]
    psi = np.zeros_like(w)
    u = np.zeros_like(w)
    v = np.zeros_like(w)
    for i in range(len(t)):
        W = np.fft.fft2(w[:, :, i])
        P, U, V = parts(W)
        psi[:, :, i] = np.real(np.fft.ifft2(P))
        u[:, :, i] = U
        v[:, :, i] = V
    return {"w": w, "psi": psi, "u": u, "v": v}, (x, y, t), err


def burgers_shock(nu=0.005, *, nx=256, nt=161, tmax=1.2, rtol=1e-11):
    """u_t = nu u_xx - u u_x from u0 = -sin(x): the classic steepening solution,
    which forms a near-discontinuity at t ~ 1 for small nu.

    The Cole-Hopf family used elsewhere CANNOT do this -- with a fixed positive
    potential its initial profile is nu times a FIXED shape, so lowering nu
    scales the field down instead of steepening it (measured: the basis
    truncation error was identical at nu = 0.02, 0.01 and 0.005). Testing a
    shock prediction on that family would have tested nothing.

    Solved with 2/3 dealiasing and a declared tolerance-ladder error. Once the
    front is thinner than the grid, the solution is no longer resolved -- which
    is the point: that is what "the information is destroyed" looks like, and
    the instrument should refuse rather than report a wide interval.
    """
    x = np.linspace(0.0, L, nx, endpoint=False)
    t = np.linspace(0.0, tmax, nt)
    k = _spectral(nx)
    kmask = k <= (2.0 / 3.0) * k.max()

    def rhs(_, u):
        U = np.fft.rfft(u)
        m = 1j * k
        m[-1] = 0.0
        ux = np.fft.irfft(m * U * kmask, nx)
        uxx = np.fft.irfft(-(k ** 2) * U * kmask, nx)
        return nu * uxx - u * ux

    out = []
    for rt in (rtol, rtol * 0.1):
        s = solve_ivp(rhs, (t[0], t[-1]), -np.sin(x), t_eval=t, method="RK45",
                      rtol=rt, atol=rt * 1e-2)
        if not s.success:
            return None, None, None
        out.append(s.y)
    return out[1], (x, t), float(np.max(np.abs(out[0] - out[1])))


def smooth_random_pair(x, t, seed=0, kmax=3):
    """Two smooth fields solving no library system: the hard null. 'Smooth and
    coupled-looking' must not be enough to certify."""
    rng = np.random.default_rng(seed)
    out = {}
    for name in ("u", "v"):
        f = np.zeros((len(x), len(t)))
        for k in range(1, kmax + 1):
            a, b = rng.normal(0, 1, 2)
            w = rng.uniform(0.5, 2.0)
            f += a * np.sin(k * x[:, None] + b) * (
                1.0 + 0.5 * np.sin(w * t[None, :]))
        out[name] = f
    return out

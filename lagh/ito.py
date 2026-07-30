"""The Itô weak form: certified drift discovery from SDE trajectories.

Level 0 of `docs/DIRECTION_STOCHASTIC.md`. Step 3 landed the VOCABULARY in
`weakform.py` -- `Term(measure="d[u]")` plus `build_nd(rough=, martingale=)`, with
`ito_terms` below emitting the identity in it and a machine-precision equality test
binding the two. The assembler here stays scalar-only on purpose: delegating it to
`build_nd` buys nothing until multi-field state is needed, which is Level 1's
requirement, so the migration is sequenced there.

The enabling identity is the one the PDE arc already built, one dimension down. For
dX = a(X)dt + b(X)dW, a compactly supported test function phi on a time window and
a smooth state function f,

    -int phi' f(X) dt  =  int phi f'(X) a(X) dt  +  1/2 int phi f''(X) b^2 dt + M
    M = int phi f'(X) b dW,        <M> = int phi^2 f'(X)^2 b^2 dt

and BOTH b^2 integrals are realized quadratic variation -- b^2 dt = d[X] -- so
they are computed from the path and need no knowledge of b at all. Moving the Itô
correction to the left,

    y  =  -int phi' f(X) dt - 1/2 int phi f''(X) d[X]  =  int phi f'(X) a(X) dt + M

is a linear system in the declared drift library, with a band whose scale is
MEASURED. Three properties carry over from the PDE arc unchanged:

* **No differentiated data.** The derivative lives on phi. A Brownian-driven path
  is nowhere differentiable, so this is not a convenience -- it is the only way
  the question is well posed.
* **The band's scale never comes from the candidate.** <M> is a functional of b,
  which is being discovered; realized quadratic variation is a functional of the
  DATA. That is what lets `certify.coverage_factor` be applied without letting a
  candidate widen its own band (docs/STOCHASTIC_CHECKER.md §1c).
* **Independent trajectories are the multi-solution holdout,** for the reason
  `pdesystem.discover_equation` gives: one path satisfies on-shell relations of
  its own realization.

## Why the f-family is not optional (measured, 2026-07-29)

With f = x alone the OU drift is NOT certifiable at any patch size, and this is
structural rather than a tuning failure: for ANY 1-D diffusion with a stationary
law, E[a(X)] = 0 (take f = x in the stationary generator identity), so the drift's
phi-weighted integral is pure FLUCTUATION of the same scale as the martingale.
Measured on OU (theta = 1.5, b = 0.4): |y| / band falls from 0.14 to 0.025 as the
patch grows from 0.1 to 40 -- the zero law certifies everywhere, so the honest
verdict is vacuity, at every scale.

f = x^2/2 changes the question. Then the drift column is int phi X^2 dt, and
E[-theta X^2] = -theta Var(X) != 0: the drift now ACCUMULATES while the martingale
grows as the root of the window. Measured on the same paths, |y| / band grows as
sqrt(L) and crosses 1 near L = 50, matching sqrt(theta L / 2) / kappa. So the
requirement is

    theta * L  >  2 kappa(n, delta)^2

which is the sharpest identifiability statement in the arc so far: a mean-reverting
drift needs a window long against its own relaxation time, and no amount of
sampling rate substitutes for it.

## Two differences from the deterministic weak form

1. **The aliasing gate does not apply and must not be run.**
   `weakform.build_nd` drops a patch whose windowed field has energy near Nyquist,
   because a smooth field with such energy is under-resolved. A Brownian path has
   energy at every frequency BY CONSTRUCTION -- the gate would reject every patch,
   and it would be answering the wrong question. What still applies is the
   refinement ladder (the path is fixed data, so subsampling it is a genuine
   coarser rule), plus a stochastic-specific gate: the quadrature bound must stay
   small against the martingale band, because a band the discretization dominates
   is not a martingale band and the coverage statement would be about something
   else. That gate is what the coarse-dt null trips.
2. **The discretization error is a CONFIDENT-WRONG risk, not just lost reach.**
   The O(dt) bias of the drift quadrature does not average away over patches while
   the martingale part does, so a certifying interval could end up centred off the
   truth and narrow enough to exclude it. It therefore enters the band as a
   DETERMINISTIC term with coefficient 1 (`certify.epsilon`'s `hard` discipline),
   never as a scale multiplied by a coverage factor.
3. **Realized quadratic variation has TWO consumers here and the same contamination
   is safe for one and unsafe for the other** -- measured, as the one
   confident-wrong Level 0 produced. It sets the band's scale AND supplies the Itô
   correction on the target. When the state is observed with error, the estimator
   picks up 2 n sigma_obs^2: in the band that is conservative, in the correction it
   is a systematic offset nothing bands. On the deterministic-decay system the true
   drift -1.0 was excluded by a joint bound of [-1088, -903]. Hence `sigma_obs`
   debiases both inside `_one`, rows the declaration explains are refused, and
   declaring it to `certify_drift` without passing it to `build_rows` RAISES: a
   declared error that does not reach the assembler is an error, not a default.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
import sympy as sp

from .certify import (KAPPA, MACHINE_REL, Abstain, admissible_interval,
                      coverage_factor, determination, free_atoms,
                      parameter_interval)
from .engine import discover
from .weakform import (MACHINE_EPS, MAX_TRUSTED_ORDER, MIN_LADDER_ORDER,
                       bump_derivatives)

# The quadrature bound may take at most this fraction of the martingale band.
# Above it the band is dominated by discretization rather than by the martingale,
# so the coverage statement kappa(n, delta) licenses would not be the statement
# being made. This is the gate the coarse-dt null is designed to trip.
QUAD_MAX_FRAC = 0.25
# When the resolution gate rejects at least this fraction of the attempted rows,
# the SAMPLING RATE is the binding constraint and the verdict says so, whatever
# the downstream engine reports. Measured on the coarse-dt null: 80% of windows
# dropped, and the survivors then abstained on vacuity -- a true statement about
# the survivors that names the wrong cause. S4 asks for the resolution reason.
RESOLUTION_BINDING_FRAC = 0.5
# Coverage factor for the MEASURED quantities' own estimation error (the realized
# quadratic variation in the band's scale, and the Ito correction on the target).
# Same 4-sigma convention as certify.LAM_B, and stated separately because it is a
# different consumer of the same physical b (docs/STOCHASTIC_CHECKER.md §4).
LAM_QV = 4.0


def _fderivs(fs):
    """State test functions f -> (name, and f through its third derivative).

    The third derivative is carried because the DIFFUSION columns need it: the
    sensitivity of 1/2 int phi f'' h dt to a per-sample state error is
    1/2 phi (f_3 h + f'' h'), and a band that dropped the f_3 part would be tight
    for a reason nobody stated. It is zero for every f in the default family --
    which is exactly why it had to be computed rather than assumed.
    """
    x = sp.Symbol("x")
    out = []
    for nm in fs:
        e = sp.sympify(nm)
        if e.free_symbols - {x}:
            raise ValueError(f"f {nm!r} reads a symbol other than x")
        out.append((str(nm),) + tuple(
            sp.lambdify(x, sp.diff(e, x, k), "numpy") for k in (0, 1, 2, 3)))
    return out


def _term_fns(names):
    """Drift library term strings -> (g, g') numpy callables of the state.

    g' is ANALYTIC, not a finite difference along the path: the path is a
    stochastic trajectory and is not monotone in time, so differencing g(X) with
    respect to X is meaningless there. It is needed for the observational Gram.
    """
    x = sp.Symbol("x")
    out = []
    for nm in names:
        e = sp.sympify(nm)
        if e.free_symbols - {x}:
            raise ValueError(f"term {nm!r} reads a symbol other than x")
        out.append((sp.lambdify(x, e, "numpy"),
                    sp.lambdify(x, sp.diff(e, x), "numpy")))
    return out


def ito_terms(f: str, library, *, field: str = "u") -> dict:
    """The Itô weak form for one state test function f, as `weakform.Term`s.

    Step 3 of `docs/DIRECTION_STOCHASTIC.md`: the identity this module assembles by
    hand, expressed in the vocabulary of record. Returns

        {"target": Term, "correction": Term, "columns": [Term...],
         "martingale": (field, gexpr)}

    where, with u the state field,

        target      (-1)^1 int (d_t phi) f(u)         == -int phi' f(u) dt
        correction  int phi (f''(u)/2) d[u]           == 1/2 int phi f'' b^2 dt
        columns     int phi f'(u) g_k(u) dt           the drift library
        martingale  the declared structure of the row's own noise, f'(u)

    and the row is `y = target - correction`. Two of the three are ordinary dt
    terms; the correction is the one lagh could not previously express, and it is a
    `measure="d[u]"` term precisely because b^2 dt IS d[u] -- so the diffusion never
    has to be modelled to write the drift's equation down.

    The martingale scale int phi^2 f'^2 d[u] is NOT among them, and cannot be: it is
    quadratic in the test function while the weak form is linear in it. It is
    returned as a declaration for `build_nd(martingale=...)`, which measures it and
    puts it beside `noise_l2` -- the same kind of object, differing only in that
    `noise_l2` multiplies a declared sigma and this one is measured.

    `tests/test_ito.py::test_the_term_vocabulary_reproduces_the_hand_rolled_rows`
    checks this against `build_rows` to machine precision. That test is the bridge
    between the two assemblers and is what keeps them from drifting apart.
    """
    from .weakform import Term
    u, x = sp.Symbol(field), sp.Symbol("x")
    fe = sp.sympify(f).xreplace({x: u})
    if fe.free_symbols - {u}:
        raise ValueError(f"f {f!r} reads a symbol other than the state")
    f1, f2 = sp.diff(fe, u), sp.diff(fe, u, 2)
    cols = []
    for g in library:
        ge = sp.sympify(g).xreplace({x: u})
        if ge.free_symbols - {u}:
            raise ValueError(f"term {g!r} reads a symbol other than the state")
        cols.append(Term(str(g), gexpr=str(sp.expand(f1 * ge)), alpha=(0,)))
    return {"target": Term(f"target:{f}", gexpr=str(fe), alpha=(1,)),
            "correction": Term(f"ito:{f}", gexpr=str(f2 / 2), alpha=(0,),
                               measure=f"d[{field}]"),
            "columns": cols,
            "martingale": (field, str(f1))}


def _at(fn, v):
    """Evaluate a lambdified scalar function on an array, broadcasting constants
    (sympy lambdifies `1` to a function returning a scalar)."""
    return np.broadcast_to(np.asarray(fn(v), float), np.shape(v)).astype(float)


@dataclass
class ItoRows:
    """Weak-form rows from one or more SDE trajectories.

    One row per (trajectory, window, f). `y` carries the target with the MEASURED
    Itô correction already subtracted; `A` the drift-library columns
    int phi f'(X) g_k(X) dt; `qv` the realized int phi^2 f'^2 d[X], which is the
    martingale variance of that row's residual.
    """
    y: np.ndarray
    A: np.ndarray
    names: list
    qv: np.ndarray
    qv_se: np.ndarray                 # sd of the qv ESTIMATOR (not of M)
    corr_se: np.ndarray               # sd of the measured Itô correction on y
    quad: np.ndarray                  # (n_rows, 1 + n_terms), target column first
    traj: np.ndarray
    fname: list = field(default_factory=list)
    windows: list = field(default_factory=list)
    gram_obs: np.ndarray | None = None   # (n_rows, 1+n_terms, 1+n_terms)
    qv_obs_share: np.ndarray | None = None
    sigma_obs_built: float = 0.0      # the sigma_obs the ROWS were assembled with
    dt: float = 0.0
    rejected: int = 0
    attempted: int = 0
    notes: list = field(default_factory=list)

    @property
    def n_trajectories(self) -> int:
        return int(len(np.unique(self.traj)))

    def n_disjoint(self) -> int:
        """Windows with DISJOINT support -- the honest held-out count for alpha
        (docs/STOCHASTIC_CHECKER.md §1d). Two rows differing only in f share the
        driving increments, so a window is counted ONCE however many f use it;
        overlapping windows on one path share increments too. Greedy interval
        scheduling per trajectory, which is exact for this.
        """
        total = 0
        for tr in np.unique(self.traj):
            wins = sorted({self.windows[i]
                           for i in np.where(self.traj == tr)[0]},
                          key=lambda w: w[1])
            end = -np.inf
            for lo, hi in wins:
                if lo >= end:
                    total += 1
                    end = hi
        return total

    def subset(self, idx) -> "ItoRows":
        idx = np.asarray(idx, int)
        return ItoRows(self.y[idx], self.A[idx], list(self.names), self.qv[idx],
                       self.qv_se[idx], self.corr_se[idx], self.quad[idx],
                       self.traj[idx], [self.fname[i] for i in idx],
                       [self.windows[i] for i in idx],
                       None if self.gram_obs is None else self.gram_obs[idx],
                       None if self.qv_obs_share is None
                       else self.qv_obs_share[idx], self.sigma_obs_built,
                       self.dt, self.rejected, self.attempted, list(self.notes))


@dataclass
class ItoBand:
    """The per-candidate band for Itô rows, one named channel per consumer.

    Four terms, and the interface's whole point is that the first three are
    different QUANTITIES in different units (`stochcheck.CONSUMER_QUANTITY`):

      kappa * sqrt(<M> + LAM_QV*se)   the INTRINSIC martingale term. kappa comes
                                      from `certify.coverage_factor(n, delta)` and
                                      <M> from realized quadratic variation --
                                      from the DATA, never from the candidate. Its
                                      estimator's own error goes INSIDE the root,
                                      because <M> is what the inequality
                                      conditions on.
      + LAM_QV * se(correction)       the measured Itô correction 1/2 int phi f''
                                      d[X] sits on the TARGET, so its estimation
                                      error is a target error, not a band scale.
      + quad . (1, |c|)               the DETERMINISTIC discretization/roundoff
                                      bound, coefficient 1 -- never multiplied by
                                      a coverage factor -- weighted by the
                                      candidate's own coefficients because the
                                      residual carries the feature columns' errors.
      + kappa * sigma_obs*sqrt(a'Ga)  the OBSERVATIONAL channel, if the state is
                                      measured with error. Zero when the path is
                                      observed exactly.

    `delta` is the false-abstain budget the certificate claims and kappa is derived
    from it and the declared row count, so the two cannot disagree.
    """
    qv: np.ndarray
    qv_se: np.ndarray
    corr_se: np.ndarray
    quad: np.ndarray
    n_rows_declared: int
    delta: float = 0.05
    sigma_obs: float = 0.0
    gram_obs: np.ndarray | None = None
    y: np.ndarray | None = None
    feat_names: list | None = None

    def __post_init__(self):
        self.kappa = coverage_factor(self.n_rows_declared, self.delta)
        self.syms = [sp.Symbol(f"x_{i}")
                     for i in range(len(self.feat_names or []))]

    def coefficients(self, expr):
        """(n_rows, n_feat) of df/dX_k, or None when the law is not linear in the
        columns -- the only case a declared library produces."""
        if expr is None:
            return None
        try:
            e = sp.sympify(expr)
            if e.free_symbols - set(self.syms):
                return None
            poly = e.as_poly(*self.syms) if self.syms else None
            if poly is None or poly.total_degree() > 1:
                return None
            c = np.array([float(e.coeff(s)) for s in self.syms])
        except Exception:                                      # noqa: BLE001
            return None
        return np.broadcast_to(c, (len(self.qv), len(self.syms))) \
            if np.all(np.isfinite(c)) else None

    def martingale(self) -> np.ndarray:
        v = np.maximum(self.qv + LAM_QV * self.qv_se, 0.0)
        return self.kappa * np.sqrt(v)

    def __call__(self, expr=None) -> np.ndarray:
        eps = self.martingale() + LAM_QV * self.corr_se
        if self.y is not None:
            eps = eps + MACHINE_REL * np.abs(self.y)
        C = self.coefficients(expr)
        if C is None:
            eps = eps + self.quad[:, 0] + self.quad[:, 1:].sum(axis=1)
        else:
            eps = eps + self.quad[:, 0] + np.einsum(
                "nk,nk->n", np.abs(C), self.quad[:, 1:])
        if self.sigma_obs > 0 and self.gram_obs is not None:
            n = len(self.qv)
            a = np.zeros((n, self.gram_obs.shape[1]))
            a[:, 0] = 1.0
            if C is not None:
                a[:, 1:] = -C
            q = np.einsum("nk,nkl,nl->n", a, self.gram_obs, a)
            eps = eps + self.kappa * self.sigma_obs * np.sqrt(np.maximum(q, 0.0))
        return eps


def time_windows(t: np.ndarray, *, half: int, n_windows: int = 0,
                 overlap: float = 0.0) -> list:
    """Index windows tiling the grid: `half` steps either side of each centre.

    `overlap` 0 gives DISJOINT supports, so the honest alpha count equals the
    window count; a positive fraction packs in more rows at the cost of the
    independence discount, which `ItoRows.n_disjoint` then reports.
    """
    n = len(t)
    step = max(1, int(round(2 * half * (1.0 - overlap))))
    centres = list(range(half, n - half - 1, step))
    if n_windows and len(centres) > n_windows:
        pick = np.linspace(0, len(centres) - 1, n_windows).round().astype(int)
        centres = [centres[i] for i in dict.fromkeys(pick)]
    return [(c - half, c + half + 1) for c in centres]


def _one(t, X, win, fns, fam, dpsi, step: int, sigma_obs: float = 0.0,
         dfns=None):
    """(y, features, corr, qv, qv_se, corr_se, absum) for one (window, f) at grid
    stride `step`, or None when the window cannot be read at that stride."""
    lo, hi = win
    tt, xx = t[lo:hi][::step], X[lo:hi][::step]
    if len(tt) < 5:
        return None
    tc = 0.5 * (tt[0] + tt[-1])
    at = 0.5 * (tt[-1] - tt[0])
    s = (tt - tc) / at
    h = float(tt[1] - tt[0])
    phi, dphi = dpsi[0](s), dpsi[1](s) / at
    _, f0, f1, f2, f3 = fam
    fx, f1x, f2x = _at(f0, xx), _at(f1, xx), _at(f2, xx)
    f3x = _at(f3, xx)
    dX = np.diff(xx)
    dX2, dX4 = dX ** 2, dX ** 4
    # the measured Ito correction, and the band's scale: both realized QV
    w_corr = 0.5 * phi[:-1] * f2x[:-1]
    w_qv = phi[:-1] ** 2 * f1x[:-1] ** 2
    corr = float(np.sum(w_corr * dX2))
    qv = float(np.sum(w_qv * dX2))
    # sd of a weighted realized-QV estimator: Var(sum w (dX)^2) = sum w^2
    # Var((dX)^2) = 2 sum w^2 (b^2 dt)^2, and E[(dX)^4] = 3 (b^2 dt)^2 gives the
    # self-normalizing estimate below -- no knowledge of b, same as the value.
    corr_se = float(np.sqrt(2.0 / 3.0 * np.sum(w_corr ** 2 * dX4)))
    qv_se = float(np.sqrt(2.0 / 3.0 * np.sum(w_qv ** 2 * dX4)))
    # OBSERVATION NOISE CONTAMINATES REALIZED QV, and not symmetrically.
    #
    # With X observed as X + e, every increment carries Var 2 sigma_obs^2 on top of
    # b^2 dt, so E[sum w (dX_obs)^2] = sum w (b^2 dt + 2 sigma_obs^2): a bias that
    # DIVERGES as dt -> 0. In the band's scale that is conservative (a wider band
    # loses laws, never admits them). In the Itô CORRECTION it is a systematic
    # offset on the TARGET, and nothing bands it -- measured as a confident-wrong on
    # the L0-ode-obs task, where the true drift -1.0 was excluded by a joint bound
    # of [-1088, -903]. Two consumers of one measured quantity, safe for one and
    # unsafe for the other: exactly the rule the checker's declarations encode.
    #
    # So both are debiased with the DECLARED sigma_obs, and the share of the raw
    # value the declaration explains is returned, so a caller can refuse a row the
    # observation noise dominates rather than trust a small difference of two large
    # numbers.
    share = 0.0
    if sigma_obs > 0:
        two_s2 = 2.0 * sigma_obs ** 2
        bias_c = float(np.sum(w_corr) * two_s2)
        bias_q = float(np.sum(w_qv) * two_s2)
        share = abs(bias_q) / max(abs(qv), 1e-300)
        corr = corr - bias_c
        qv = qv - bias_q
        # the debias has its own error: Var((de)^2) for de ~ N(0, 2 sigma^2) is
        # 2(2 sigma^2)^2, and adjacent increments share a sample, so the factor 3
        # below is deliberately conservative rather than exact
        corr_se = float(np.hypot(corr_se,
                                 np.sqrt(3.0 * np.sum(w_corr ** 2)) * two_s2))
        qv_se = float(np.hypot(qv_se,
                               np.sqrt(3.0 * np.sum(w_qv ** 2)) * two_s2))
    if dfns:
        # THE DIFFUSION IS CLAIMED, NOT MEASURED. b^2 = sum_j d_j h_j(x) goes into
        # the design matrix as the ordinary dt columns 1/2 int phi f'' h_j dt, so
        # drift and diffusion are identified JOINTLY from the same rows, and the
        # target keeps the plain -int phi' f dt with no correction subtracted.
        #
        # This also removes the one UNSAFE consumer of realized quadratic
        # variation: the measured correction was the term where observation-noise
        # contamination became a systematic offset on the target (the Level 0
        # confident-wrong). Here realized QV is used only to set the BAND, where
        # contamination is conservative.
        corr, corr_se = 0.0, 0.0
    y = -float(np.sum(dphi * fx) * h) - corr
    feats = np.array([float(np.sum(phi * f1x * _at(g, xx)) * h) for g, _ in fns]
                     + [0.5 * float(np.sum(phi * f2x * _at(hh, xx)) * h)
                        for hh, _ in (dfns or [])])
    absum = float(np.sum(np.abs(dphi * fx)) * h)
    # FIRST-ORDER sensitivity of each functional to a per-sample state error, for
    # the observational channel's Gram. The target reads -phi' f(X) so its
    # sensitivity is -phi' f'(X) h; column k reads phi f'(X) g_k(X) so its own is
    # phi (f'' g_k + f' g_k')(X) h. Independent per-sample error makes the band one
    # quadratic form a'Ga with a = (1, -c), the same shape weakform._noise_gram
    # assembles for a field.
    nu = [-dphi * f1x * h]
    for g, dgf in fns:
        nu.append(phi * (f2x * _at(g, xx) + f1x * _at(dgf, xx)) * h)
    for hh, dhf in (dfns or []):
        nu.append(0.5 * phi * (f3x * _at(hh, xx) + f2x * _at(dhf, xx)) * h)
    NU = np.vstack(nu)
    gram = NU @ NU.T
    return y, feats, corr, qv, qv_se, corr_se, absum, share, gram


def build_rows(t, paths, names, *, fs=("x", "x**2/2"), diff_names=None,
               windows=None, half: int = 200, n_windows: int = 0,
               overlap: float = 0.0, p: int = 8, delta: float = 0.05,
               sigma_obs: float = 0.0, qv_obs_max: float = 0.5) -> ItoRows:
    """Assemble Itô weak-form rows from `paths` (n_traj, n_steps) on grid `t`.

    `fs` is the family of state test functions. `("x",)` alone reproduces the
    plain dX form and, per the module docstring, cannot certify a stationary
    drift at any patch size -- `x**2/2` is what makes one identifiable, so it is
    in the default.

    `diff_names` turns the DIFFUSION into a claim instead of a measurement. Given a
    library {h_j} for b^2, the Itô correction becomes the ordinary dt columns
    1/2 int phi f''(X) h_j(X) dt and drift and diffusion are identified JOINTLY from
    the same rows -- the Level 1 target. Two things follow and neither is
    incidental:

      * **The Delta-t requirement disappears from the diffusion.** Nothing here
        estimates b from quadratic variation, so the classic dt -> 0 demand applies
        to the QV estimator and not to this. S2 was registered in the QV picture and
        is measured against this one.
      * **The unsafe consumer of realized QV is gone.** The measured correction was
        the one place where observation-noise contamination became a systematic
        offset on the TARGET (the Level 0 confident-wrong). With the diffusion in
        the design matrix, realized QV only sets the BAND, where contamination is
        conservative.

    Columns are named `drift:<g>` and `diffusion:<h>`, which is what
    `stochcheck.component` expects.

    `sigma_obs` is the DECLARED per-sample measurement error on the state. It is
    not cosmetic: it debiases realized quadratic variation (see `_one`) and it
    builds the observational channel's Gram. A row whose raw quadratic variation
    is more than `qv_obs_max` explained by that declaration is REFUSED -- past
    that point the process diffusion is a small difference of two large numbers
    and separating process noise from measurement noise is what Level 2 exists to
    supply ground truth for.
    """
    t = np.asarray(t, float)
    P = np.atleast_2d(np.asarray(paths, float))
    fns = _term_fns(names)
    dfns = _term_fns(diff_names) if diff_names else None
    colnames = ([f"drift:{g}" for g in names]
                + [f"diffusion:{h}" for h in (diff_names or [])]) \
        if diff_names else list(names)
    fam = _fderivs(fs)
    dpsi = bump_derivatives(p, 1)
    wins = windows if windows is not None else time_windows(
        t, half=half, n_windows=n_windows, overlap=overlap)
    kappa = coverage_factor(max(len(wins) * len(P) * len(fam), 1), delta)
    ys, As, qvs, qses, cses, quads, trs, fnm, wl = [], [], [], [], [], [], [], [], []
    grams, shares = [], []
    rejected, gate_hits, obs_hits = 0, [], []
    for k in range(len(P)):
        X = P[k]
        for win in wins:
            for fa in fam:
                base = _one(t, X, win, fns, fa, dpsi, 1, sigma_obs, dfns)
                lad = [_one(t, X, win, fns, fa, dpsi, s, sigma_obs, dfns)
                       for s in (2, 4)]
                if base is None or any(z is None for z in lad):
                    rejected += 1
                    continue
                y, feats, corr, qv, qv_se, corr_se, absum, share, gram = base
                if sigma_obs > 0 and (share > qv_obs_max or qv <= 0.0):
                    rejected += 1
                    obs_hits.append(float(share))
                    continue
                round_h = MACHINE_EPS * absum
                # THE LADDER, in the one axis there is. Only the phi-quadrature
                # part is laddered: the realized-QV terms are not quadratures of a
                # smooth integrand and subsampling them changes the estimator
                # rather than refining it, so their error is the statistical one
                # above, declared separately.
                d1 = abs((lad[0][0] + lad[0][2]) - (y + corr))
                d2 = abs((lad[1][0] + lad[1][2]) - (lad[0][0] + lad[0][2]))
                if d1 <= round_h:
                    q_t = round_h
                else:
                    r = np.log2(d2 / d1) if d1 > 0 and d2 > 0 else 0.0
                    if np.isfinite(r) and r >= MIN_LADDER_ORDER:
                        q_t = d1 / (2.0 ** float(np.clip(
                            r, MIN_LADDER_ORDER, MAX_TRUSTED_ORDER)) - 1.0)
                    else:
                        q_t = d1     # not asymptotic: the raw difference stands
                q_f = np.abs(lad[0][1] - feats)
                band_m = kappa * np.sqrt(max(qv + LAM_QV * qv_se, 0.0))
                tot_q = max(q_t, round_h) + float(q_f.sum())
                if not np.isfinite(band_m) or band_m <= 0 \
                        or tot_q > QUAD_MAX_FRAC * band_m:
                    rejected += 1
                    gate_hits.append(tot_q / band_m if band_m > 0 else np.inf)
                    continue
                ys.append(y)
                As.append(feats)
                qvs.append(qv)
                qses.append(qv_se)
                cses.append(corr_se)
                quads.append(np.concatenate([[max(q_t, round_h)], q_f]))
                trs.append(k)
                fnm.append(fa[0])
                wl.append((float(t[win[0]]), float(t[win[1] - 1])))
                grams.append(gram)
                shares.append(share)
    notes = []
    if obs_hits:
        notes.append(
            f"{len(obs_hits)} rows refused because the DECLARED sigma_obs "
            f"explains a median {np.median(obs_hits):.0%} of their raw quadratic "
            f"variation (bar {qv_obs_max:.0%}): past that the process diffusion "
            "is a small difference of two large numbers, and separating process "
            "from measurement noise is Level 2's job, not a band's")
    if rejected:
        med = float(np.median(gate_hits)) if gate_hits else float("nan")
        notes.append(
            f"{rejected} of {len(wins) * len(P) * len(fam)} (window, f) rows "
            f"dropped by the resolution gate: the quadrature bound was a median "
            f"{med:.3g} of the martingale band against a {QUAD_MAX_FRAC:g} bar, "
            "so the band would have been a discretization band wearing a "
            "coverage statement (Abstain.RESOLUTION territory)")
    K = 1 + len(colnames)
    return ItoRows(np.array(ys),
                   np.array(As) if As else np.zeros((0, len(colnames))),
                   list(colnames), np.array(qvs), np.array(qses), np.array(cses),
                   np.array(quads) if quads else np.zeros((0, K)),
                   np.array(trs, int), fnm, wl,
                   np.array(grams) if grams else np.zeros((0, K, K)),
                   np.array(shares) if shares else np.zeros(0), float(sigma_obs),
                   float(t[1] - t[0]), rejected,
                   len(wins) * len(P) * len(fam), notes)


# ---------------------------------------------------------------------------
# The DIFFUSION, from quadratic variation (docs/CASE_STUDY_STOCHASTIC_L1.md)
# ---------------------------------------------------------------------------
# Measured 2026-07-29: putting b^2 in the drift's design matrix gives a joint bound
# 2000x the truth, because the diffusion's signal-to-band goes as b while the
# drift's goes as 1/b -- the thing the diffusion measures IS the noise. Realized
# quadratic variation determines the same b^2 to ~1.4% on the same data. So the
# diffusion gets its own weak form, with the SAME discipline and a different target:
#
#     int phi w(X) d[X]  =  sum_j d_j int phi w(X) h_j(X) dt
#
# the left side observable as sum_i (phi w)_i (dX_i)^2 and the right side ordinary
# dt columns. `w` is a family of state weights and plays exactly the role f plays
# for the drift: one row per (window, w), all sharing the coefficients d_j.

# The declared bound on |a(x)| used for the drift-leakage term below. Declared,
# then VERIFIED against a drift certificate when one exists -- the same
# declared-then-checked discipline weakform.declared_epsilon uses for coeff_max.
DRIFT_MAX_DEFAULT = 10.0


@dataclass
class QvRows:
    """Quadratic-variation rows: one per (trajectory, window, w)."""
    y: np.ndarray                     # int phi w d[X], debiased for sigma_obs
    A: np.ndarray                      # int phi w h_j dt
    names: list
    var: np.ndarray                    # Var of the y estimator, MEASURED
    lam: np.ndarray                    # int |phi w| dt, for the flat leakage bound
    leak: np.ndarray | None = None     # int |phi w| env(X)^2 dt, when an envelope
                                       # for |a(x)| was supplied (the tight bound)
    traj: np.ndarray = None
    wname: list = field(default_factory=list)
    windows: list = field(default_factory=list)
    dt: float = 0.0
    sigma_obs_built: float = 0.0
    notes: list = field(default_factory=list)

    @property
    def n_trajectories(self) -> int:
        return int(len(np.unique(self.traj)))

    def n_disjoint(self) -> int:
        total = 0
        for tr in np.unique(self.traj):
            wins = sorted({self.windows[i]
                           for i in np.where(self.traj == tr)[0]},
                          key=lambda z: z[1])
            end = -np.inf
            for lo, hi in wins:
                if lo >= end:
                    total += 1
                    end = hi
        return total


@dataclass
class QvBand:
    """The band for a quadratic-variation row.

        kappa * sqrt(Var)              the estimator's own fluctuation. Var is
                                       MEASURED from the FOURTH moment of the
                                       increments: E[(dX)^4] = 3 (b^2 dt)^2 while
                                       Var((dX)^2) = 2 (b^2 dt)^2, so
                                       Var = (2/3) sum (phi w)^2 (dX)^4 needs no
                                       knowledge of b -- the same property that
                                       made the martingale scale measurable.
        + dt * drift_max^2 * lam       DRIFT LEAKAGE, deterministic, coefficient 1.
                                       E[(dX)^2] = b^2 dt + a^2 dt^2, so the drift
                                       biases this estimator at O(dt) exactly as the
                                       diffusion biased the drift's -- the same
                                       coupling, the other way round. `drift_max` is
                                       declared and must be verified.
        + machine

    kappa comes from `certify.coverage_factor` as everywhere else, and the
    justification differs in one stated way: the increments (phi w)((dX)^2 - E) are
    martingale differences with CHI-SQUARE tails rather than a continuous
    martingale, so the applicable inequality is Bernstein's, whose sub-exponential
    correction is a factor 1/(1 + c*kappa/sqrt(Var)) inside the exponent with
    c ~ max_i (phi w)_i b^2 dt. That correction is O(sqrt(dt/L)) -- negligible when
    a window holds many increments -- and `report` returns its measured size rather
    than asserting it is small.
    """
    var: np.ndarray
    lam: np.ndarray
    n_rows_declared: int
    dt: float
    delta: float = 0.05
    drift_max: float = DRIFT_MAX_DEFAULT
    y: np.ndarray | None = None
    max_term: np.ndarray | None = None      # max_i |phi w|_i (dX_i)^2, for the
                                            # Bernstein correction's scale c
    leak: np.ndarray | None = None          # POINTWISE leakage, when an envelope
                                            # for |a(x)| was supplied

    def __post_init__(self):
        self.kappa = coverage_factor(self.n_rows_declared, self.delta)

    def bernstein_correction(self) -> np.ndarray:
        """c*kappa/sqrt(Var) per row: how far the chi-square tails move the bound
        away from the Gaussian-like reading. Reported, never assumed away."""
        if self.max_term is None:
            return np.zeros(len(self.var))
        return (self.max_term * self.kappa
                / np.maximum(np.sqrt(np.maximum(self.var, 0.0)), 1e-300))

    def __call__(self, expr=None) -> np.ndarray:
        eps = self.kappa * np.sqrt(np.maximum(self.var, 0.0))
        # POINTWISE leakage where an envelope was supplied, the flat bound
        # otherwise. The difference is large and it is the same distinction
        # weakform draws with `field_l1`: sum_i |w_i| a(X_i)^2 dt^2 bounded by the
        # envelope AT X_i, rather than by its worst value anywhere in the visited
        # range. A process spends its time where its own drift is small -- most of
        # all a bistable one, whose wells sit at a(x) = 0 -- so the flat bound
        # over-declares by the square of a ratio that can be an order of magnitude.
        if self.leak is not None:
            eps = eps + self.dt * self.leak
        else:
            eps = eps + self.dt * self.drift_max ** 2 * self.lam
        if self.y is not None:
            eps = eps + MACHINE_REL * np.abs(self.y)
        return eps


def build_qv_rows(t, paths, diff_names, *, ws=("1", "x", "x**2"), windows=None,
                  half: int = 200, n_windows: int = 0, overlap: float = 0.0,
                  p: int = 8, delta: float = 0.05, sigma_obs: float = 0.0,
                  qv_obs_max: float = 0.5, drift_envelope=None) -> QvRows:
    """Assemble quadratic-variation rows for the DIFFUSION.

    `ws` is the family of state weights, playing the role f plays for the drift.
    What it buys depends on the process, and this was MEASURED after being asserted
    wrongly: the w-family supplies state variation when the process is STATIONARY,
    where every window sees the same state distribution -- on OU the joint bounds
    tighten 2.9-5.6x going from `("1",)` to `("1", "x", "x**2")`. When the process
    is NON-STATIONARY the windows already sample different state regions and w adds
    nothing: on GBM, whose state grows exponentially, the same change left the
    bounds marginally WIDER (0.0079 -> 0.0083 on the x^2 coefficient), because the
    extra rows cost a little kappa and bought no new information.
    """
    t = np.asarray(t, float)
    P = np.atleast_2d(np.asarray(paths, float))
    hfns = _term_fns(diff_names)
    wfns = _term_fns(ws)
    dpsi = bump_derivatives(p, 0)
    wins = windows if windows is not None else time_windows(
        t, half=half, n_windows=n_windows, overlap=overlap)
    ys, As, vs, lams, lks, trs, wn, wl = [], [], [], [], [], [], [], []
    rejected, obs_hits = 0, []
    dt = float(t[1] - t[0])
    for k in range(len(P)):
        X = P[k]
        for lo, hi in wins:
            tt, xx = t[lo:hi], X[lo:hi]
            if len(tt) < 5:
                rejected += 1
                continue
            tc = 0.5 * (tt[0] + tt[-1])
            at = 0.5 * (tt[-1] - tt[0])
            phi = dpsi[0]((tt - tc) / at)
            dX2 = np.diff(xx) ** 2
            dX4 = dX2 ** 2
            for wnm, (wf, _) in zip(ws, wfns):
                wx = _at(wf, xx)
                wt = (phi * wx)[:-1]
                q = float(np.sum(wt * dX2))
                var = float(2.0 / 3.0 * np.sum(wt ** 2 * dX4))
                share = 0.0
                if sigma_obs > 0:
                    bias = float(np.sum(wt) * 2.0 * sigma_obs ** 2)
                    share = abs(bias) / max(abs(q), 1e-300)
                    q = q - bias
                    var = var + 3.0 * float(np.sum(wt ** 2)) \
                        * (2.0 * sigma_obs ** 2) ** 2
                    if share > qv_obs_max:
                        rejected += 1
                        obs_hits.append(share)
                        continue
                cols = np.array([float(np.sum(phi * wx * _at(hf, xx)) * dt)
                                 for hf, _ in hfns])
                ys.append(q)
                As.append(cols)
                vs.append(var)
                aw = np.abs(phi * wx)
                lams.append(float(np.sum(aw) * dt))
                if drift_envelope is not None:
                    env = np.abs(np.asarray(drift_envelope(xx), float))
                    lks.append(float(np.sum(aw * env ** 2) * dt))
                trs.append(k)
                wn.append(wnm)
                wl.append((float(tt[0]), float(tt[-1])))
    notes = []
    if obs_hits:
        notes.append(
            f"{len(obs_hits)} rows refused: the declared sigma_obs explains a "
            f"median {np.median(obs_hits):.0%} of their raw quadratic variation "
            f"(bar {qv_obs_max:.0%}), so b^2 would be a small difference of two "
            "large numbers")
    return QvRows(np.array(ys),
                  np.array(As) if As else np.zeros((0, len(diff_names))),
                  [f"diffusion:{h}" for h in diff_names], np.array(vs),
                  np.array(lams),
                  np.array(lks) if lks else None,
                  np.array(trs, int), wn, wl, dt, float(sigma_obs), notes)


def certify_diffusion(rows: QvRows, *, delta: float = 0.05,
                      drift_max: float = DRIFT_MAX_DEFAULT,
                      sigma_obs: float = 0.0, seed: int = 0,
                      max_tier: int = 3, holdout: bool = True,
                      qualifier: dict | None = None) -> dict:
    """Certify b^2 as a linear law over the declared diffusion library.

    Mirrors `certify_drift`: same engine, same holdout-by-trajectory, and the
    determination record comes from the joint LP so interval coverage means "the
    truth is in here".
    """
    if abs(float(sigma_obs) - rows.sigma_obs_built) > 0:
        raise ValueError(
            f"sigma_obs={sigma_obs!r} declared to certify_diffusion but the rows "
            f"were built with {rows.sigma_obs_built!r}; it debiases the estimator "
            "inside build_qv_rows and cannot be applied afterwards")
    n = len(rows.y)
    out = {"n_rows": n, "n_trajectories": rows.n_trajectories,
           "n_disjoint": rows.n_disjoint(), "features": list(rows.names),
           "ws": sorted(set(rows.wname)), "notes": list(rows.notes)}
    if n < 6:
        out.update(certified=False, abstain=Abstain.RANGE.value)
        out["partial"] = determination([], status=Abstain.RANGE.value,
                                       note="too few rows", qualifier=qualifier)
        return out
    if holdout and rows.n_trajectories < 2:
        out.update(certified=False, abstain="single-trajectory")
        out["partial"] = determination([], status="single-trajectory",
                                       note="nothing was determined",
                                       qualifier=qualifier)
        return out
    mx = None
    band = QvBand(rows.var, rows.lam, n, rows.dt, delta=delta,
                  drift_max=drift_max, y=rows.y, max_term=mx,
                  leak=rows.leak)
    e_all = band(None)
    out.update(kappa=band.kappa, delta=delta,
               median_band=float(np.median(e_all)),
               median_signal_to_band=float(np.median(
                   np.abs(rows.y) / np.maximum(e_all, 1e-300))),
               median_drift_leakage_share=float(np.median(
                   band.dt * (rows.leak if rows.leak is not None
                              else band.drift_max ** 2 * rows.lam)
                   / np.maximum(e_all, 1e-300))),
               leakage_bound=("pointwise envelope" if rows.leak is not None
                              else "flat drift_max"),
               drift_max_declared=drift_max)
    rng = np.random.default_rng(seed)
    if holdout:
        held = np.unique(rows.traj)[-1]
        tr = rng.permutation(np.where(rows.traj != held)[0])
        ce = np.where(rows.traj == held)[0]
    else:
        idx = rng.permutation(n)
        b = int(0.8 * n)
        tr, ce = idx[:b], idx[b:]
    if len(ce) < 3 or len(tr) < 4:
        out.update(certified=False, abstain=Abstain.RANGE.value)
        out["partial"] = determination([], status=Abstain.RANGE.value,
                                       note="splits too small",
                                       qualifier=qualifier)
        return out
    a = int(0.75 * len(tr))
    eps_ce = QvBand(rows.var[ce], rows.lam[ce], n, rows.dt, delta=delta,
                    drift_max=drift_max, y=rows.y[ce],
                    leak=None if rows.leak is None else rows.leak[ce])
    sigma_eff = float(np.median(e_all)
                      / max(KAPPA * np.median(np.abs(rows.y)), 1e-300))
    out["sigma_effective"] = sigma_eff
    r = discover(rows.A[tr[:a]], rows.y[tr[:a]], rows.A[tr[a:]], rows.y[tr[a:]],
                 rows.A[ce], rows.y[ce], sigma=sigma_eff, eps_model=eps_ce,
                 max_tier=max_tier, declared_basis=True, linear_basis=True,
                 band_sel=QvBand(rows.var[tr[a:]], rows.lam[tr[a:]], n, rows.dt,
                                 delta=delta, drift_max=drift_max,
                                 y=rows.y[tr[a:]],
                                 leak=None if rows.leak is None
                                 else rows.leak[tr[a:]])(None))
    c = r.certificate
    out.update(certified=bool(c.certified), abstain=c.abstain,
               alpha_log10=c.alpha_log10, n_cert_rows=int(len(ce)), tier=r.tier)
    out["notes"] += [str(x)[:220] for x in c.notes][:3]
    e0 = (band.kappa * np.sqrt(np.maximum(rows.var[ce], 0.0))
          + MACHINE_REL * np.abs(rows.y[ce]))
    lp, info = admissible_interval(
        rows.A[ce], rows.y[ce],
        lambda cm: e0 + rows.dt * (rows.leak[ce] if rows.leak is not None
                                   else drift_max ** 2 * rows.lam[ce]),
        coeff_max=10.0)
    out["admissible_info"] = info
    if lp is not None:
        out["admissible"] = {nm: [v[0], v[1]] for nm, v in zip(rows.names, lp)}
        out["partial"] = determination(
            [(nm, v[0], v[1]) for nm, v in zip(rows.names, lp)],
            status=("certified" if c.certified
                    else (c.abstain or "structural-abstain")),
            note=("JOINT bound by LP over every b^2 in the declared library "
                  "consistent with the held-out quadratic-variation rows at the "
                  f"declared band, coverage 1-delta={1 - delta:g}"),
            qualifier=qualifier)
    if c.certified:
        out["law"] = _readable(str(r.expr), list(rows.names))
        out["expr"] = str(r.expr)
    return out


def certify_drift(rows: ItoRows, *, delta: float = 0.05, sigma_obs: float = 0.0,
                  seed: int = 0, max_tier: int = 3, holdout: bool = True,
                  qualifier: dict | None = None) -> dict:
    """Certify the drift as a linear law over the declared library.

    Same shape as `pdesystem.discover_equation`, and for the same reason: the
    engine's coherence, significance, parsimony and parameter-interval machinery
    IS the discipline, and re-deriving it here would be re-deriving it wrongly.
    What is Itô-specific is the band (`ItoBand`) and the holdout unit -- a
    TRAJECTORY, this arc's "several solutions".

    Returns a dict carrying the verdict, the coverage statement, the disjoint row
    count alpha must be read at, and a `certify.determination` record: exactly
    what `stochcheck.Submission` consumes.
    """
    # A DECLARED ERROR THAT DOES NOT REACH THE ASSEMBLER IS AN ERROR, not a
    # no-op. sigma_obs debiases realized quadratic variation inside `build_rows`;
    # declaring it only here would leave the Ito correction biased and the band
    # unaware, which is exactly how the L0-ode-obs confident-wrong was produced.
    if abs(float(sigma_obs) - rows.sigma_obs_built) > 0:
        raise ValueError(
            f"sigma_obs={sigma_obs!r} was declared to certify_drift but the rows "
            f"were assembled with sigma_obs={rows.sigma_obs_built!r}. Pass it to "
            "build_rows as well: it debiases realized quadratic variation, and "
            "without that the Ito correction carries a systematic offset nothing "
            "bands (measured as a confident-wrong on L0-ode-obs)")
    n = len(rows.y)
    out = {"n_rows": n, "n_trajectories": rows.n_trajectories,
           "n_disjoint": rows.n_disjoint(), "rejected": rows.rejected,
           "features": list(rows.names), "fs": sorted(set(rows.fname)),
           "notes": list(rows.notes)}
    if n < 6:
        out.update(certified=False, abstain=Abstain.RESOLUTION.value)
        out["notes"].append(
            f"{n} rows survived the resolution gate: at this sampling rate the "
            "quadrature error is not small against the martingale band, so no "
            "coverage statement can be made. A resolution refusal, not a wide "
            "interval.")
        out["partial"] = determination([], status=Abstain.RESOLUTION.value,
                                       note="nothing was determined",
                                       qualifier=qualifier)
        return out
    if holdout and rows.n_trajectories < 2:
        out.update(certified=False, abstain="single-trajectory")
        out["notes"].append(
            "a drift claim needs rows from a HELD-OUT trajectory: one path "
            "supports only an on-shell statement about that realization")
        out["partial"] = determination([], status="single-trajectory",
                                       note="nothing was determined",
                                       qualifier=qualifier)
        return out
    band = ItoBand(rows.qv, rows.qv_se, rows.corr_se, rows.quad, n, delta=delta,
                   sigma_obs=sigma_obs, gram_obs=rows.gram_obs, y=rows.y,
                   feat_names=list(rows.names))
    e_all = band(None)
    out["kappa"] = band.kappa
    out["delta"] = delta
    out["median_band"] = float(np.median(e_all))
    out["median_martingale_share"] = float(np.median(
        band.martingale() / np.maximum(e_all, 1e-300)))
    out["median_signal_to_band"] = float(np.median(
        np.abs(rows.y) / np.maximum(e_all, 1e-300)))
    rng = np.random.default_rng(seed)
    if holdout:
        held = np.unique(rows.traj)[-1]
        tr = rng.permutation(np.where(rows.traj != held)[0])
        ce = np.where(rows.traj == held)[0]
    else:
        idx = rng.permutation(n)
        b = int(0.8 * n)
        tr, ce = idx[:b], idx[b:]
    if len(ce) < 3 or len(tr) < 4:
        out.update(certified=False, abstain=Abstain.RANGE.value)
        out["notes"].append("splits too small after the trajectory holdout")
        # every return path carries a record: a consumer must never have to tell
        # "nothing was determined" apart from "the field is missing"
        out["partial"] = determination([], status=Abstain.RANGE.value,
                                       note="nothing was determined",
                                       qualifier=qualifier)
        return out
    a = int(0.75 * len(tr))
    eps_ce = _sub_band(band, ce, rows, delta, sigma_obs)
    # THE EFFECTIVE RELATIVE NOISE, declared. `discover`'s `sigma` does not widen
    # the band when an eps_model is supplied -- it drives the GATES: at sigma <= 0
    # `float_pinned` snaps a coefficient to an exact rational and `pinned` stands
    # down entirely. Passing 0 therefore tells the engine this data is CLEAN, and
    # it is not: measured on GBM, the drift came back as 495334*x/606799, an exact
    # rational claim about a diffusion-scale coefficient -- precisely the false
    # exactness RNOISE_STUDY.md's parametric gate exists to refuse. So the band is
    # translated into the relative scale it corresponds to and declared.
    sigma_eff = float(np.median(e_all) / max(KAPPA * np.median(np.abs(rows.y)),
                                             1e-300))
    out["sigma_effective"] = sigma_eff
    out["notes"].append(
        f"effective relative noise declared to the gates: {sigma_eff:.3g} "
        f"(median band / (KAPPA * median |y|)). The certification band is the "
        "martingale band, not this number; this only puts the exactness and "
        "parametric gates in their noise regime, where they belong")
    r = discover(rows.A[tr[:a]], rows.y[tr[:a]], rows.A[tr[a:]], rows.y[tr[a:]],
                 rows.A[ce], rows.y[ce], sigma=sigma_eff, eps_model=eps_ce,
                 max_tier=max_tier, declared_basis=True, linear_basis=True,
                 band_sel=_sub_band(band, tr[a:], rows, delta, sigma_obs)(None))
    c = r.certificate
    out.update(certified=bool(c.certified), abstain=c.abstain,
               alpha_log10=c.alpha_log10, n_cert_rows=int(len(ce)),
               tier=r.tier)
    out["notes"] += [str(x)[:220] for x in c.notes][:3]
    # NAME THE BINDING CONSTRAINT. When the resolution gate took most of the
    # evidence, the sampling rate is what refused -- the survivors' own verdict
    # (usually vacuity) is true about the survivors and points at the wrong cause.
    frac = rows.rejected / max(rows.attempted, 1)
    out["rejected_fraction"] = float(frac)
    if not c.certified and frac >= RESOLUTION_BINDING_FRAC:
        out["downstream_abstain"] = c.abstain
        out["abstain"] = Abstain.RESOLUTION.value
        out["notes"].append(
            f"reason RESOLUTION: the gate took {frac:.0%} of the attempted rows, "
            f"so the sampling rate is the binding constraint; the "
            f"{len(ce)}-row survivors then abstained on {c.abstain!r}, which is "
            "true of them and names the wrong cause")
    out["n_disjoint_cert"] = rows.subset(ce).n_disjoint()
    # THE INTERVAL RECORD COMES FROM THE LP, on both paths.
    #
    # `parameter_interval` bisects the certification predicate with every OTHER
    # parameter HELD -- a conditional statement about one law, which is not the
    # claim interval-coverage scoring reads ("the truth is in here"). And
    # `invariant_content` ranges over the certifying laws the SEARCH FOUND, which
    # was measured to exclude a law that certifies (see admissible_interval).
    # The LP over the consistent polytope is the bound that actually holds.
    lp, lp_info = _admissible(rows, ce, band, delta, sigma_obs)
    out["admissible_info"] = lp_info
    if lp is not None:
        out["admissible"] = {nm: [None if v[0] is None else float(v[0]),
                                  None if v[1] is None else float(v[1])]
                             for nm, v in zip(rows.names, lp)}
        out["partial"] = determination(
            [(nm, v[0], v[1]) for nm, v in zip(rows.names, lp)],
            status=("certified" if c.certified
                    else (c.abstain or "structural-abstain")),
            note=("JOINT bound by LP over every law in the declared vocabulary "
                  "consistent with the held-out rows at the martingale band: if "
                  "the truth is in the vocabulary and the band covers, the truth "
                  "is in these ranges"),
            qualifier=qualifier)
    out["alpha_kind"] = (
        f"stochastic: coverage 1-delta={1 - delta:g} at kappa={band.kappa:.3f} "
        f"over {n} rows. NOT comparable with a deterministic alpha, and to be "
        f"read at n_disjoint_cert={out['n_disjoint_cert']} rather than "
        f"n_cert_rows={len(ce)} -- q^h assumes independent held-out rows and "
        "overlapping windows share the driving path")
    if not c.certified:
        # the search-set report travels too, clearly labelled as such -- it says
        # what the SEARCH agreed on, which is a different and weaker thing
        out["search_set_report"] = _relabel(c.partial, list(rows.names),
                                            qualifier, c.abstain
                                            or "structural-abstain")
        return out
    syms = [sp.Symbol(f"x_{i}") for i in range(len(rows.names))]
    out["law"] = _readable(str(r.expr), list(rows.names))
    out["expr"] = str(r.expr)
    e = sp.expand(sp.sympify(r.expr))
    coeffs = {nm: float(e.coeff(s)) for nm, s in zip(rows.names, syms)}
    out["coefficients"] = coeffs
    ivs = {}
    atoms = free_atoms(r.expr)
    for nm, s in zip(rows.names, syms):
        atom = next((a for a in atoms
                     if abs(float(a) - coeffs[nm])
                     <= 1e-12 * max(1.0, abs(coeffs[nm]))), None)
        if coeffs[nm] == 0.0 or atom is None:
            ivs[nm] = None
            continue
        ivs[nm] = parameter_interval(r.expr, syms, rows.A[ce], rows.y[ce],
                                     eps_ce, atom)
    # CONDITIONAL, and labelled as such: the range over which THIS law still
    # certifies with every OTHER coefficient held. Informative about the certified
    # law, and deliberately NOT the determination record, because it is not a bound
    # on the true coefficient -- the joint LP above is.
    out["conditional_intervals"] = {
        k: (None if v is None else [float(v[0]), float(v[1])])
        for k, v in ivs.items()}
    out["conditional_note"] = (
        "one coefficient varied with the others HELD: a statement about this law, "
        "not a bound on the truth. The determination record carries the joint LP "
        "bound instead")
    return out


def _admissible(rows: ItoRows, idx, band: ItoBand, delta, sigma_obs):
    """The joint LP bound on every coefficient, over the CERTIFICATION rows.

    The band handed to the LP must not depend on the candidate, or the polytope
    would not be one. The candidate-independent part (martingale + measured-
    correction error + machine + the target's own quadrature) is taken as is, and
    the feature columns' quadrature enters as `coeff_max` times their sum -- the
    declared-then-verified discipline, with `admissible_interval` raising the
    declaration until the answer respects it.
    """
    idx = np.asarray(idx, int)
    sub = _sub_band(band, idx, rows, delta, sigma_obs)
    e0 = (sub.martingale() + LAM_QV * rows.corr_se[idx]
          + MACHINE_REL * np.abs(rows.y[idx]) + rows.quad[idx, 0])
    qsum = rows.quad[idx, 1:].sum(axis=1)
    if sigma_obs > 0 and sub.gram_obs is not None and len(sub.gram_obs):
        # the observational channel is candidate-dependent too; bound it at the
        # declared coeff_max the same way rather than dropping it
        d = np.sqrt(np.maximum(np.diagonal(sub.gram_obs, axis1=1, axis2=2), 0.0))
        e0 = e0 + sub.kappa * sigma_obs * d[:, 0]
        qsum = qsum + sub.kappa * sigma_obs * d[:, 1:].sum(axis=1)
    return admissible_interval(rows.A[idx], rows.y[idx],
                               lambda c: e0 + c * qsum, coeff_max=10.0)


def _sub_band(band: ItoBand, idx, rows: ItoRows, delta, sigma_obs) -> ItoBand:
    """A band on a SUBSET of rows that keeps the parent's coverage factor.

    kappa is a property of the exhaustive check being CLAIMED -- over all n rows
    the run assembled -- not of whichever split the engine happens to be looking
    at. Recomputing it per split would quietly shrink the band on the
    certification subset, which is the impostor-admitting direction.
    """
    idx = np.asarray(idx, int)
    g = None if rows.gram_obs is None or not len(rows.gram_obs) \
        else rows.gram_obs[idx]
    return ItoBand(rows.qv[idx], rows.qv_se[idx], rows.corr_se[idx],
                   rows.quad[idx], band.n_rows_declared, delta=delta,
                   sigma_obs=sigma_obs, gram_obs=g, y=rows.y[idx],
                   feat_names=list(rows.names))


def _readable(expr: str, names) -> str:
    out = expr
    for i, nm in reversed(list(enumerate(names))):
        out = out.replace(f"x_{i}", f"[{nm}]")
    return out


def _relabel(partial, names, qualifier, status):
    """Re-key an engine-level partial record by library term name and attach the
    run's domain -- the same job as `pdesystem._relabel` / `_qualify`."""
    if not partial:
        return determination([], status=status,
                             note="nothing was determined over this domain",
                             qualifier=qualifier)
    out = dict(partial)
    co = out.get("coefficients")
    if isinstance(co, dict):
        ren = {}
        for k, v in co.items():
            bit = k.split("_")[-1]
            ren[names[int(bit)] if k.startswith("x_") and bit.isdigit()
                and int(bit) < len(names) else k] = v
        out["coefficients"] = ren
    if qualifier is not None:
        out["qualifier"] = qualifier
    return out

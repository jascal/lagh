"""Weak-form feature factory for PDE claims (docs/DIRECTION_PDE.md,
registration docs/CASE_STUDY_PDE_DEV.md).

lagh certifies |pred - y| <= eps per point under a DECLARED error model.
Derivative fields estimated by finite differences carry neither a measured nor
a declarable error (sigma/h^k amplification plus truncation), so a PDE claim
made against them would be exactly the undeclared-error claim the float32
lesson forbids. The weak form removes the problem instead of estimating it:
integrate every candidate term against a smooth compactly-supported test
function phi and move all derivatives onto phi by parts,

    int phi * d^alpha( g(u) )  =  (-1)^|alpha| int (d^alpha phi) * g(u)

so the data enters only through a POINTWISE g of the raw field, weighted by an
ANALYTIC known function. Each patch integral is one "data point"; the
certification domain is the patch family.

The library is divergence-form terms {d^alpha(g(u))} only -- u*u_xx and friends
cannot be moved onto phi and are out of reach by construction (stated in the
registration, not hidden).

Error accounting per patch, all parts computed or declared:
  * float summation R -- rigorous: eps_mach * sum|contributions| (the quadrature
    sum cancels heavily; this bounds its roundoff),
  * quadrature truncation Q -- MEASURED on a three-level (h, 2h, 4h) ladder with
    the Richardson estimate at h as the declared value; a patch whose ladder does
    not converge is REJECTED as unresolved rather than silently kept,
  * declared noise -- sigma * ||w||_2 through the known linear functional
    (returned per patch; C0 runs sigma = 0, see the registration on why the
    noisy case needs an errors-in-variables eps model first).
"""
from __future__ import annotations

from dataclasses import dataclass, field
from functools import lru_cache

import numpy as np
import sympy as sp

from .certify import KAPPA, MACHINE_REL   # one epsilon vocabulary

MACHINE_EPS = float(np.finfo(float).eps)

# A patch is unresolved when the declared quadrature bound exceeds this fraction
# of the patch's own term scale -- the band would be a sizable part of the
# signal. This is a coarse sanity filter only: whether a patch family is
# INFORMATIVE is a separate, reported quantity (WeakSystem.signal_to_band),
# because a wide band does not make a certificate wrong, only weak.
UNRESOLVED_REL = 1e-3
# ...and when the (h, 2h, 4h) ladder is not in the asymptotic regime: below this
# observed order refinement is not buying accuracy and the Richardson estimate
# of the h-level error means nothing.
MIN_LADDER_ORDER = 2.0
# The most convergence the declared bound will credit from a coarse ladder pair
# (see build(): trusting a reported order of 10 under-declares the band).
MAX_TRUSTED_ORDER = 4.0
# Spectral-energy fraction above half-Nyquist that marks a patch as aliased (see
# aliasing_ratio): the ladder alone cannot detect an under-resolved field.
# MEASURED separation on 33x17-point patches: resolved fields (modes 1-16) sit at
# 7e-7 - 1.2e-6 -- the bump window's own spectral tail, which is the floor --
# while a 4-points-per-wavelength field and white noise sit at 0.3 - 0.8. The bar
# sits in the middle of that five-decade gap.
ALIAS_MAX = 1e-3


@dataclass(frozen=True)
class Term:
    """One divergence-form library term: (-1)^|alpha| int (d^alpha phi) g(u)."""
    name: str
    ax: int                      # x-derivatives moved onto phi
    at: int                      # t-derivatives moved onto phi
    gexpr: str                   # pointwise function of u, as a sympy string

    @property
    def sign(self) -> int:
        return (-1) ** (self.ax + self.at)

    def g(self, u: np.ndarray) -> np.ndarray:
        f = _gfun(self.gexpr)
        return np.broadcast_to(np.asarray(f(u), float), u.shape)

    def dg(self, u: np.ndarray) -> np.ndarray:
        """g'(u): how a perturbation of the field at a point moves this term's
        integrand -- the sensitivity the errors-in-variables band is built from."""
        f = _gfun(str(sp.diff(sp.sympify(self.gexpr), sp.Symbol("u"))))
        return np.broadcast_to(np.asarray(f(u), float), u.shape)


@lru_cache(maxsize=None)
def _gfun(gexpr: str):
    return sp.lambdify(sp.Symbol("u"), sp.sympify(gexpr), "numpy")


# The C0 curriculum library: heat, advection, wave, Burgers, KdV, and the
# constant/source terms that let a law be wrong in a detectable way.
LIBRARY: dict[str, Term] = {
    t.name: t for t in [
        Term("1", 0, 0, "1"),
        Term("u", 0, 0, "u"),
        Term("u^2", 0, 0, "u**2"),
        Term("u_x", 1, 0, "u"),
        Term("u_xx", 2, 0, "u"),
        Term("u_xxx", 3, 0, "u"),
        Term("u_xxxx", 4, 0, "u"),
        Term("u*u_x", 1, 0, "u**2/2"),
        Term("u^2*u_x", 1, 0, "u**3/3"),
        Term("u_t", 0, 1, "u"),
        Term("u_tt", 0, 2, "u"),
    ]
}


def bump_derivatives(p: int, max_order: int) -> list:
    """psi(s) = (1 - s^2)^p on |s| < 1, zero outside, and its derivatives as
    numpy callables. psi is C^(p-1) with every derivative up to p-1 vanishing at
    the support boundary -- which is what makes the by-parts identity exact (no
    boundary terms) and the uniform-grid quadrature converge fast."""
    s = sp.Symbol("s")
    psi = (1 - s ** 2) ** p
    out = []
    for k in range(max_order + 1):
        d = sp.diff(psi, s, k)
        f = sp.lambdify(s, d, "numpy")

        def make(f=f):
            def call(v):
                v = np.asarray(v, float)
                y = np.zeros_like(v)
                m = np.abs(v) < 1.0
                if np.any(m):
                    y[m] = np.asarray(f(v[m]), float)
                return y
            return call
        out.append(make())
    return out


@dataclass
class Patch:
    """A test-function support: center and half-widths in (x, t), plus the grid
    index window it covers."""
    xc: float
    tc: float
    ax: float
    at: float
    ix: slice
    it: slice


@dataclass
class WeakSystem:
    """The design matrix over patches. `A[:, k]` is term `names[k]`; the runner
    picks which column is the target."""
    A: np.ndarray
    names: list[str]
    patches: list[Patch] = field(default_factory=list)
    quad: np.ndarray | None = None      # (n_patches, n_terms) declared quadrature bound
    roundoff: np.ndarray | None = None  # (n_patches, n_terms) rigorous float bound
    noise_l2: np.ndarray | None = None  # (n_patches, n_terms) ||w||_2: sigma multiplier
    order: np.ndarray | None = None     # (n_patches, n_terms) observed ladder order
    rejected: int = 0                   # patches dropped as unresolved
    gram: np.ndarray | None = None      # (n_patches, n_terms, n_terms) noise Gram

    def signal_to_band(self, target: str, **kw) -> float:
        """median |target| / declared band: how much evidence a patch family
        actually carries. A family near 1 certifies nothing worth having; the
        runner reports this rather than hiding it inside a filter."""
        y = self.A[:, self.names.index(target)]
        return float(np.median(np.abs(y) / self.declared_epsilon(target, **kw)))

    def declared_epsilon(self, target: str, *, coeff_max: float = 10.0,
                         sigma: float = 0.0) -> np.ndarray:
        """Per-patch declared band for the TARGET residual.

        The certified relation is y = sum_k c_k X_k, so the residual carries the
        target column's own error plus the feature columns' errors weighted by
        the coefficients. The coefficients are not known when eps is assembled,
        so they are bounded by a DECLARED `coeff_max` -- which the runner must
        check against the certified law and re-assemble if it was too small
        (an assumption stated and verified, never assumed silently)."""
        j = self.names.index(target)
        err = (self.quad + self.roundoff)
        if sigma > 0:
            err = err + sigma * self.noise_l2
        others = np.delete(err, j, axis=1)
        return err[:, j] + coeff_max * others.sum(axis=1)


@dataclass
class PatchEpsilon:
    """The errors-in-variables band for weak-form rows: a CALLABLE epsilon
    (`certify.band`) assembled per candidate law.

    Why it has to be per candidate: the residual of y = sum_k c_k X_k is
        (delta_y - sum_k c_k delta_k)  +  (deterministic quadrature/roundoff),
    and every delta is a functional of the SAME field noise, so the stochastic
    part is `KAPPA * sigma * ||nu_y - sum_k c_k nu_k||_2` -- one norm over the
    combined sensitivity vector, computed here as the quadratic form a'Ga with
    a = (1, -c). Bounding the coefficients instead (the C0 `coeff_max` band)
    both loosens epsilon and smuggles in an assumption; loosening epsilon is the
    direction that admits impostors.

    For a candidate that is not linear in the columns the same first-order
    formula holds with c_k = df/dX_k evaluated PER ROW, which is what
    `coefficients` returns. If the gradient cannot be evaluated, the band falls
    back to the conservative `coeff_max` form rather than guessing.
    """
    names: list[str]                 # term order of the Gram / quad columns
    target: str
    y: np.ndarray                    # (n_rows,) target column
    X: np.ndarray                    # (n_rows, n_feat) feature columns
    det: np.ndarray                  # (n_rows, n_terms) quadrature + roundoff
    gram: np.ndarray                 # (n_rows, n_terms, n_terms)
    sigma: float = 0.0
    floor_abs: float = 0.0
    coeff_max: float = 2.0

    def __post_init__(self):
        self.j = self.names.index(self.target)
        self.feat = [n for n in self.names if n != self.target]
        self.cols = [self.names.index(n) for n in self.feat]
        self.syms = [sp.Symbol(f"x_{i}") for i in range(len(self.feat))]

    def subset(self, idx) -> "PatchEpsilon":
        idx = np.asarray(idx, int)
        return PatchEpsilon(self.names, self.target, self.y[idx], self.X[idx],
                            self.det[idx], self.gram[idx], self.sigma,
                            self.floor_abs, self.coeff_max)

    def coefficients(self, expr):
        """(n_rows, n_feat) of df/dX_k, or None when the law cannot be read that
        way (non-finite, or not a function of the columns)."""
        if expr is None:
            return None
        try:
            e = sp.sympify(expr)
            if e.free_symbols - set(self.syms):
                return None
            grads = [sp.lambdify(self.syms, sp.diff(e, s), "numpy")
                     for s in self.syms]
            cols = [np.broadcast_to(np.asarray(
                g(*[self.X[:, i] for i in range(len(self.feat))]), float),
                (len(self.y),)) for g in grads]
        except Exception:                                      # noqa: BLE001
            return None
        C = np.column_stack(cols) if cols else np.zeros((len(self.y), 0))
        return C if np.all(np.isfinite(C)) else None

    def __call__(self, expr=None) -> np.ndarray:
        C = self.coefficients(expr)
        n, m = len(self.y), len(self.names)
        eps = MACHINE_REL * np.abs(self.y) + self.floor_abs
        if C is None:
            # conservative fallback: coefficients unknown, so bound them, and
            # add the stochastic parts by triangle inequality (no cancellation)
            d = self.det[:, self.j] + self.coeff_max * np.delete(
                self.det, self.j, axis=1).sum(axis=1)
            if self.sigma > 0:
                sd = np.sqrt(np.maximum(np.diagonal(self.gram, axis1=1,
                                                    axis2=2), 0.0))
                d = d + KAPPA * self.sigma * (
                    sd[:, self.j] + self.coeff_max
                    * np.delete(sd, self.j, axis=1).sum(axis=1))
            return eps + d
        a = np.zeros((n, m))
        a[:, self.j] = 1.0
        for k, col in enumerate(self.cols):
            a[:, col] = -C[:, k]
        eps = eps + self.det[:, self.j] + np.einsum(
            "nk,nk->n", np.abs(C), self.det[:, self.cols])
        if self.sigma > 0:
            q = np.einsum("nk,nkl,nl->n", a, self.gram, a)
            eps = eps + KAPPA * self.sigma * np.sqrt(np.maximum(q, 0.0))
        return eps


def make_patches(x: np.ndarray, t: np.ndarray, *, nx_half: int, nt_half: int,
                 n_x: int, n_t: int) -> list[Patch]:
    """A grid of patch centers whose supports lie strictly inside the domain.
    `nx_half`/`nt_half` are half-widths in GRID CELLS (the ladder needs the
    window to stay a multiple of 4)."""
    x = np.asarray(x, float)
    t = np.asarray(t, float)
    hx = float(x[1] - x[0])
    ht = float(t[1] - t[0])
    xs = np.linspace(nx_half, len(x) - 1 - nx_half, n_x).astype(int)
    ts = np.linspace(nt_half, len(t) - 1 - nt_half, n_t).astype(int)
    out = []
    for i in xs:
        for j in ts:
            out.append(Patch(float(x[i]), float(t[j]), nx_half * hx,
                             nt_half * ht,
                             slice(i - nx_half, i + nx_half + 1),
                             slice(j - nt_half, j + nt_half + 1)))
    return out


def aliasing_ratio(u_patch: np.ndarray, patch: Patch, xs: np.ndarray,
                   ts: np.ndarray, dpsi, sigma: float = 0.0) -> float:
    """Fraction of the windowed patch field's spectral energy sitting in the top
    half of the representable band.

    The refinement ladder alone can be FOOLED (measured): a field at ~4 points
    per wavelength aliases to a different low frequency at every coarse level,
    so the three integrals agree smoothly and report a high convergence order
    while all three are garbage. This is the direct test -- energy near Nyquist
    means the grid cannot represent the field, let alone its derivatives.
    Windowing by the test function first is what makes the DFT meaningful on a
    non-periodic patch (psi*u is smooth and compactly supported)."""
    wx = dpsi[0]((xs - patch.xc) / patch.ax)
    wt = dpsi[0]((ts - patch.tc) / patch.at)
    W = np.outer(wx, wt)
    F = np.abs(np.fft.fft2(W * u_patch)) ** 2
    kx = np.fft.fftfreq(F.shape[0])                      # cycles/sample
    kt = np.fft.fftfreq(F.shape[1])
    mask = (np.abs(kx) > 0.25)[:, None] | (np.abs(kt) > 0.25)[None, :]
    # DECLARED NOISE IS NOT ALIASING: white noise puts energy in every mode, so
    # an honest resolution test subtracts what the declared sigma is expected to
    # contribute (E|F_k|^2 = sigma^2 * sum(W^2) per mode, uniformly) before
    # asking whether the FIELD has unrepresentable structure. Without this the
    # gate rejects patches for carrying the noise the eps model already bands.
    per_mode = sigma ** 2 * float((W ** 2).sum())
    hi = float(F[mask].sum()) - per_mode * int(mask.sum())
    tot = float(F.sum()) - per_mode * F.size
    if tot <= 0:
        return 1.0 if hi > 0 else 0.0
    return max(0.0, hi) / tot


def _weights(term: Term, xs: np.ndarray, ts: np.ndarray, patch: Patch,
             dpsi) -> np.ndarray:
    """The analytic weight array (-1)^|alpha| (d^alpha phi)(z_i) * hx * ht."""
    hx = float(xs[1] - xs[0])
    ht = float(ts[1] - ts[0])
    wx = dpsi[term.ax]((xs - patch.xc) / patch.ax) / patch.ax ** term.ax
    wt = dpsi[term.at]((ts - patch.tc) / patch.at) / patch.at ** term.at
    return np.outer(wx, wt) * (hx * ht) * term.sign


def _integrate(term: Term, u_patch: np.ndarray, xs: np.ndarray, ts: np.ndarray,
               patch: Patch, dpsi) -> tuple:
    """One patch integral at one resolution: value, sum|contributions|,
    ||w||_2 (the noise multiplier for a linear-in-u term)."""
    W = _weights(term, xs, ts, patch, dpsi)
    C = W * term.g(u_patch)
    return (float(C.sum()), float(np.abs(C).sum()),
            float(np.sqrt(np.sum(W ** 2))))


def _ladder_noise(term: Term, u_patch: np.ndarray, xs: np.ndarray,
                  ts: np.ndarray, patch: Patch, dpsi, sa: int, sb: int) -> float:
    """||d||_2 for the difference functional I_{sa*h} - I_{sb*h}: the ladder
    difference is ALSO a linear functional of the field noise, with a computable
    weight vector, so sigma * this is exactly how much of an observed ladder
    difference the declared noise explains."""
    d = np.zeros(u_patch.shape)
    for s, sign in ((sa, 1.0), (sb, -1.0)):
        sub = (slice(None, None, s), slice(None, None, s))
        d[sub] += sign * (_weights(term, xs[::s], ts[::s], patch, dpsi)
                          * term.dg(u_patch[sub]))
    return float(np.sqrt((d ** 2).sum()))


def _noise_gram(terms, u_patch: np.ndarray, xs: np.ndarray, ts: np.ndarray,
                patch: Patch, dpsi) -> np.ndarray:
    """Gram matrix of the per-term noise SENSITIVITY vectors on this patch:
    G[k,l] = sum_i (W_k g_k'(u))_i (W_l g_l'(u))_i.

    To first order a field perturbation e moves term k's integral by
    <W_k g_k'(u), e>, so the residual of y = sum_k c_k X_k moves by <v, e> with
    v = nu_target - sum_k c_k nu_k -- ONE vector, because every column is a
    functional of the SAME realization. Its length is a'Ga with
    a = (1, -c_1, ..., -c_K), which is why the Gram is what gets stored: the
    band for any candidate is then a quadratic form costing K^2 flops per patch
    instead of a pass over the grid."""
    nu = np.stack([(_weights(tm, xs, ts, patch, dpsi)
                    * tm.dg(u_patch)).ravel() for tm in terms])
    return nu @ nu.T


def build(u: np.ndarray, x: np.ndarray, t: np.ndarray, terms: list[str],
          patches: list[Patch], *, p: int = 8, sigma: float = 0.0) -> WeakSystem:
    """The factory: field u[i, j] on grid (x[i], t[j]) -> patch design matrix
    with a declared error bound per entry.

    Each entry is computed at three resolutions (h, 2h, 4h) on the SAME patch;
    the declared quadrature bound comes from that ladder and the observed order
    is recorded. A patch is DROPPED, never silently kept, when it is aliased
    (`aliasing_ratio` -- the ladder alone cannot see an under-resolved field),
    when the ladder is not in the asymptotic regime, or when the resulting bound
    is a sizable fraction of the patch's own term scale.

    `p` is the bump exponent: it sets how many derivatives of phi vanish at the
    support edge and so how fast the quadrature converges for the HIGH-derivative
    weights. Measured on the heat field, the u_xxx bound falls 1.5e-6 -> 7.5e-11
    from p = 8 to p = 16 at fixed patch size; since the declared band is set by
    the roughest term in the library, a library with third derivatives wants
    p >= 12.
    """
    u = np.asarray(u, float)
    x = np.asarray(x, float)
    t = np.asarray(t, float)
    T = [LIBRARY[n] if isinstance(n, str) else n for n in terms]
    maxo = max(max(tm.ax, tm.at) for tm in T)
    dpsi = bump_derivatives(p, maxo)
    rows, quad, rnd, l2, orders, kept, grams = [], [], [], [], [], [], []
    rejected = 0
    for pa in patches:
        xs_f, ts_f = x[pa.ix], t[pa.it]
        up = u[pa.ix, pa.it]
        vals, qs, rs, ls, ords = [], [], [], [], []
        ok = aliasing_ratio(up, pa, xs_f, ts_f, dpsi, sigma) <= ALIAS_MAX
        if not ok:
            rejected += 1
            continue
        for tm in T:
            ladder = []
            for step in (1, 2, 4):
                v, absum, w2 = _integrate(tm, up[::step, ::step], xs_f[::step],
                                          ts_f[::step], pa, dpsi)
                ladder.append((v, absum, w2))
            v_h, absum, w2 = ladder[0]
            d1 = abs(ladder[1][0] - ladder[0][0])          # |I_2h - I_h|
            d2 = abs(ladder[2][0] - ladder[1][0])          # |I_4h - I_2h|
            scale = max(absum, MACHINE_EPS)
            round_h = MACHINE_EPS * absum
            # Under declared noise the ladder differences are dominated by the
            # NOISE, not by truncation, and the raw convergence test then rejects
            # patches for carrying noise the eps model already bands (measured:
            # 16 of 24 heat patches lost at sigma = 1e-4). Subtract what sigma
            # explains before asking whether refinement bought accuracy.
            ns1 = ns2 = 0.0
            if sigma > 0:
                ns1 = sigma * _ladder_noise(tm, up, xs_f, ts_f, pa, dpsi, 1, 2)
                ns2 = sigma * _ladder_noise(tm, up, xs_f, ts_f, pa, dpsi, 2, 4)
                d1 = max(0.0, d1 - KAPPA * ns1)
                d2 = max(0.0, d2 - KAPPA * ns2)
            if d1 <= max(round_h, ns1):
                # truncation is not detectable above the roundoff/noise floor:
                # declare it AT that floor rather than pretending to resolve it
                q, r = max(round_h, ns1), np.inf
            else:
                # ASYMPTOTIC-REGIME TEST: the Richardson estimate is only
                # meaningful once refinement actually buys accuracy. Observed
                # order r = log2(|I_4h - I_2h| / |I_2h - I_h|); below order 2 the
                # ladder is not converging and no declared bound would be honest.
                r = np.log2(d2 / d1) if d2 > 0 and d1 > 0 else 0.0
                if not np.isfinite(r) or r < MIN_LADDER_ORDER:
                    ok = False
                    break
                # |e_h| <= d1 / (2^r - 1) holds when the error really shrinks by
                # 2^r per refinement. The ORDER USED IS CAPPED (measured reason:
                # for the exponentially-convergent bump quadrature the coarse
                # (4h, 2h) pair reports r ~ 10 while the true h-level ratio is
                # ~2^7 -- an overstated order under-declares the band and
                # rejects true laws). The observed order is recorded separately.
                r_used = float(np.clip(r, MIN_LADDER_ORDER, MAX_TRUSTED_ORDER))
                q = d1 / (2.0 ** r_used - 1.0)
                if q > UNRESOLVED_REL * scale:
                    ok = False                 # resolved, but too coarsely to use
                    break
            vals.append(v_h)
            qs.append(max(q, round_h))
            rs.append(round_h)
            ls.append(w2)
            ords.append(float(r))
        if not ok:
            rejected += 1
            continue
        rows.append(vals)
        quad.append(qs)
        rnd.append(rs)
        l2.append(ls)
        orders.append(ords)
        kept.append(pa)
        grams.append(_noise_gram(T, up, xs_f, ts_f, pa, dpsi))
    return WeakSystem(A=np.array(rows, float), names=[tm.name for tm in T],
                      patches=kept, quad=np.array(quad, float),
                      roundoff=np.array(rnd, float),
                      noise_l2=np.array(l2, float),
                      order=np.array(orders, float), rejected=rejected,
                      gram=np.array(grams, float) if grams else None)

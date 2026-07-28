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
                   ts: np.ndarray, dpsi) -> float:
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
    f = np.outer(wx, wt) * u_patch
    F = np.abs(np.fft.rfft2(f)) ** 2
    tot = float(F.sum())
    if tot <= 0:
        return 0.0
    kx = np.fft.fftfreq(f.shape[0])           # cycles/sample, |kx| <= 1/2
    kt = np.fft.rfftfreq(f.shape[1])
    mask = (np.abs(kx) > 0.25)[:, None] | (kt > 0.25)[None, :]
    return float(F[mask].sum()) / tot


def _integrate(term: Term, u_patch: np.ndarray, xs: np.ndarray, ts: np.ndarray,
               patch: Patch, dpsi) -> tuple:
    """One patch integral at one resolution: value, sum|contributions|,
    ||w||_2 (the noise multiplier for a linear-in-u term)."""
    hx = float(xs[1] - xs[0])
    ht = float(ts[1] - ts[0])
    sx = (xs - patch.xc) / patch.ax
    st = (ts - patch.tc) / patch.at
    wx = dpsi[term.ax](sx) / patch.ax ** term.ax
    wt = dpsi[term.at](st) / patch.at ** term.at
    W = np.outer(wx, wt) * (hx * ht) * term.sign
    C = W * term.g(u_patch)
    return (float(C.sum()), float(np.abs(C).sum()),
            float(np.sqrt(np.sum(W ** 2))))


def build(u: np.ndarray, x: np.ndarray, t: np.ndarray, terms: list[str],
          patches: list[Patch], *, p: int = 8) -> WeakSystem:
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
    rows, quad, rnd, l2, orders, kept = [], [], [], [], [], []
    rejected = 0
    for pa in patches:
        xs_f, ts_f = x[pa.ix], t[pa.it]
        up = u[pa.ix, pa.it]
        vals, qs, rs, ls, ords = [], [], [], [], []
        ok = aliasing_ratio(up, pa, xs_f, ts_f, dpsi) <= ALIAS_MAX
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
            if d1 <= round_h:
                q, r = round_h, np.inf        # converged into the roundoff floor
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
    return WeakSystem(A=np.array(rows, float), names=[tm.name for tm in T],
                      patches=kept, quad=np.array(quad, float),
                      roundoff=np.array(rnd, float),
                      noise_l2=np.array(l2, float),
                      order=np.array(orders, float), rejected=rejected)

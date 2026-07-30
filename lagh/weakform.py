"""Weak-form feature factory for PDE claims (docs/DIRECTION_PDE.md,
registration docs/CASE_STUDY_PDE_DEV.md).

lagh certifies |pred - y| <= eps per point under a DECLARED error model.
Derivative fields estimated by finite differences carry neither a measured nor
a declarable error (sigma/h^k amplification plus truncation), so a PDE claim
made against them would be exactly the undeclared-error claim the float32
lesson forbids. The weak form removes the problem instead of estimating it:
integrate every candidate term against a smooth compactly-supported test
function phi and move all derivatives onto phi by parts,

    int phi * d^alpha( g(fields) )  =  (-1)^|alpha| int (d^alpha phi) * g(fields)

so the data enters only through a POINTWISE g of the raw field(s), weighted by
an ANALYTIC known function. Each patch integral is one "data point"; the
certification domain is the patch family.

The library is divergence-form terms {d^alpha(g)} only -- u*u_xx and friends
cannot be moved onto phi and are out of reach by construction (stated in the
registration, not hidden). `g` may depend on SEVERAL fields (u*v, u**2*v, h*u),
which is what makes systems of PDEs reachable with no change to the engine
(docs/DIRECTION_PDE_SYSTEMS.md): a coupled system is one row set with one
target per equation and features spanning every field.

Error accounting per patch, all parts computed or declared:
  * float summation R -- rigorous: eps_mach * sum|contributions| (the quadrature
    sum cancels heavily; this bounds its roundoff),
  * quadrature truncation Q -- MEASURED on a three-level (h, 2h, 4h) ladder with
    the Richardson estimate at h as the declared value; a patch whose ladder does
    not converge is REJECTED as unresolved rather than silently kept,
  * declared field noise -- sigma * ||w||_2 through the known linear functional
    (per patch, and as a Gram over terms so the errors-in-variables band can be
    assembled per candidate: see PatchEpsilon),
  * declared FIELD error -- a field produced by a reference solver carries the
    solver's own declared error, which is deterministic (no cancellation across
    a realization) and so enters through the L1 sensitivity ||w g'||_1 rather
    than the L2 one. `WeakSystem.field_l1` carries it; `WeakSystem.det()`
    multiplies it by the declared bound.

Geometry is n-dimensional: coordinates are (x, ..., t) with TIME LAST. The 1-D+t
entry points (`build`, `make_patches`) are the special case every campaign
through C2 used; `build_nd` / `make_patches_nd` take an arbitrary coordinate
list, which is what 2-D fields (Navier-Stokes vorticity, the 2-D half of
PDEBench) need.
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

FIELD_SEP = ":"          # system term names are "<field>:<term>" (u:u_xx, v:u)


@dataclass(frozen=True)
class Term:
    """One divergence-form library term: (-1)^|alpha| int (d^alpha phi) g dmu.

    `gexpr` is a pointwise function of ONE OR MORE named fields; the symbols it
    is written in ARE the fields it reads (`u*v` reads u and v). `ax`/`at` are
    the 1-D+t convenience form of the derivative multi-index; `alpha` states it
    in full for n-D geometry (space axes first, TIME LAST).

    `measure` is the dmu above and is `"dt"` for every deterministic term -- the
    Lebesgue quadrature every campaign through the PDE arc uses. `"d[<field>]"`
    integrates against that field's REALIZED QUADRATIC VARIATION along the last
    axis instead, summing (d^alpha phi) g (Delta field)^2. It exists for the Itô
    weak form (`lagh/ito.py`, docs/DIRECTION_STOCHASTIC.md step 3): the Itô
    correction 1/2 int phi f'' b^2 dt IS 1/2 int phi f'' d[u], so the one term
    lagh could not previously express becomes a library term like any other, and
    the field's diffusion never has to be modelled to write it down.

    Two properties of a d[] term that follow from what it is, not from choice:
    it carries NO grid-spacing factor (the squared increment is the measure), and
    it does not refine under the quadrature ladder (subsampling changes the
    estimator rather than resolving it), so its error is statistical and declared
    separately. Both are handled in `build_nd`.
    """
    name: str
    ax: int = 0
    at: int = 0
    gexpr: str = "u"
    alpha: tuple | None = None
    measure: str = "dt"

    @property
    def qv_fields(self) -> tuple | None:
        """The field(s) whose (co)variation this term integrates against, or None
        for an ordinary dt term.

        `d[u]` is the quadratic variation of u and `d[u,v]` the CROSS-variation of
        the pair. The cross form is what multi-dimensional Itô needs: for
        dX_i = a_i dt + sum_j b_ij dW_j the correction is
        1/2 sum_ik int phi (d2f/dx_i dx_k) d[X_i, X_k], and the off-diagonal terms
        are exactly the ones a diagonal-only measure cannot express.
        """
        m = str(self.measure)
        if m == "dt":
            return None
        if m.startswith("d[") and m.endswith("]") and len(m) > 3:
            parts = tuple(p.strip() for p in m[2:-1].split(","))
            if len(parts) in (1, 2) and all(parts):
                return parts * 2 if len(parts) == 1 else parts
        raise ValueError(f"term {self.name!r}: measure {m!r} is neither 'dt', "
                         "'d[<field>]' nor 'd[<field>,<field>]'")

    @property
    def qv_field(self) -> str | None:
        """The first (co)variation field, kept for the scalar callers."""
        f = self.qv_fields
        return None if f is None else f[0]

    def multi(self, ndim: int) -> tuple:
        """The derivative multi-index over `ndim` axes (space..., time)."""
        if self.alpha is not None:
            if len(self.alpha) != ndim:
                raise ValueError(f"term {self.name!r}: alpha has "
                                 f"{len(self.alpha)} axes, geometry has {ndim}")
            return tuple(int(a) for a in self.alpha)
        if ndim < 2:
            raise ValueError("geometry needs at least one space axis and time")
        return (self.ax,) + (0,) * (ndim - 2) + (self.at,)

    def order(self, ndim: int) -> int:
        return sum(self.multi(ndim))

    @property
    def sign(self) -> int:
        """(-1)^|alpha| -- the by-parts sign, independent of the geometry it is
        evaluated in (the multi-index total is the same either way)."""
        return (-1) ** (sum(self.alpha) if self.alpha is not None
                        else self.ax + self.at)

    @property
    def fields(self) -> tuple:
        """The field names g reads, in sorted order (empty for a constant)."""
        return _gsyms(self.gexpr)

    def g(self, fields) -> np.ndarray:
        f, shape = _as_fields(fields)
        syms = self.fields
        out = _gfun(self.gexpr, syms)(*[f[s] for s in syms])
        return np.broadcast_to(np.asarray(out, float), shape)

    def dg(self, fields, wrt: str | None = None) -> np.ndarray:
        """dg/d(field): how a perturbation of that field at a point moves this
        term's integrand -- the sensitivity every band is built from. `wrt=None`
        means the single-field default, kept for the scalar campaigns."""
        f, shape = _as_fields(fields)
        if wrt is None:
            wrt = next(iter(f))
        syms = self.fields
        if wrt not in syms:
            return np.zeros(shape)
        d = str(sp.diff(sp.sympify(self.gexpr), sp.Symbol(wrt)))
        out = _gfun(d, syms)(*[f[s] for s in syms])
        return np.broadcast_to(np.asarray(out, float), shape)


def _as_fields(fields):
    """Accept a bare array (the single-field scalar case) or a name->array dict."""
    if isinstance(fields, dict):
        shape = np.shape(next(iter(fields.values())))
        return {k: np.asarray(v, float) for k, v in fields.items()}, shape
    a = np.asarray(fields, float)
    return {"u": a}, a.shape


@lru_cache(maxsize=None)
def _gsyms(gexpr: str) -> tuple:
    return tuple(sorted(str(s) for s in sp.sympify(gexpr).free_symbols))


@lru_cache(maxsize=None)
def _gfun(gexpr: str, syms: tuple):
    return sp.lambdify([sp.Symbol(s) for s in syms], sp.sympify(gexpr), "numpy")


# The C0 curriculum library: heat, advection, wave, Burgers, KdV, and the
# constant/source terms that let a law be wrong in a detectable way.
LIBRARY: dict[str, Term] = {
    t.name: t for t in [
        Term("1", 0, 0, "1"),
        Term("u", 0, 0, "u"),
        Term("u^2", 0, 0, "u**2"),
        Term("u^3", 0, 0, "u**3"),
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


def field_terms(fieldname: str, names) -> list:
    """The scalar LIBRARY rewritten in one field's symbol, named `<field>:<term>`.

    For a LINEAR system every feature is a single-field term of one field or the
    other, which is why the system probe needed no factory extension at all; this
    just makes that vocabulary explicit and unambiguous when several fields share
    one row set."""
    out = []
    for n in names:
        base = LIBRARY[n] if isinstance(n, str) else n
        gx = str(sp.sympify(base.gexpr).xreplace(
            {sp.Symbol("u"): sp.Symbol(fieldname)}))
        out.append(Term(f"{fieldname}{FIELD_SEP}{base.name}", base.ax, base.at,
                        gx, base.alpha))
    return out


def bump_derivatives(p: int, max_order: int, onesided: bool = False) -> list:
    """psi(s) = (1 - s^2)^p and its derivatives as numpy callables.

    Two-sided (the default): support |s| < 1. psi is C^(p-1) with every
    derivative up to p-1 vanishing at BOTH support boundaries -- which is what
    makes the by-parts identity exact (no boundary terms) and the uniform-grid
    quadrature converge fast.

    One-sided: support 0 <= s < 1, i.e. the SAME psi restricted to the right
    half. It still vanishes to order p at s = 1, but psi(0) = 1, so the by-parts
    identity in time keeps its boundary term at the window's left edge -- which
    is exactly the term a state certificate reads the initial condition out of
    (docs/DIRECTION_PDE_STATE.md). psi'(0) = 0, so the window is still smooth
    where it meets the initial time.
    """
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
                m = (v >= 0.0) & (v < 1.0) if onesided else np.abs(v) < 1.0
                if np.any(m):
                    y[m] = np.asarray(f(v[m]), float)
                return y
            return call
        out.append(make())
    return out


@dataclass
class Patch:
    """A test-function support: centre and half-width per axis (space..., time
    last), plus the grid index window it covers.

    `onesided` marks a time window running from its LEFT edge (the centre, in
    that case) to centre + halfwidth, with the test function equal to 1 at the
    left edge: the shape a state certificate needs.
    """
    centers: tuple
    halfwidths: tuple
    idx: tuple
    onesided: bool = False

    # 1-D+t names, kept because every campaign through C2 speaks them
    @property
    def xc(self) -> float:
        return self.centers[0]

    @property
    def tc(self) -> float:
        return self.centers[-1]

    @property
    def ax(self) -> float:
        return self.halfwidths[0]

    @property
    def at(self) -> float:
        return self.halfwidths[-1]

    @property
    def ix(self):
        return self.idx[0]

    @property
    def it(self):
        return self.idx[-1]

    @property
    def ndim(self) -> int:
        return len(self.centers)


@dataclass
class WeakSystem:
    """The design matrix over patches. `A[:, k]` is term `names[k]`; the runner
    picks which column is the target."""
    A: np.ndarray
    names: list[str]
    patches: list[Patch] = field(default_factory=list)
    quad: np.ndarray | None = None      # (n_patches, n_terms) declared quadrature bound
    roundoff: np.ndarray | None = None  # (n_patches, n_terms) rigorous float bound
    noise_l2: np.ndarray | None = None  # (n_patches, n_terms) ||w g'||_2: sigma multiplier
    field_l1: np.ndarray | None = None  # (n_patches, n_terms) ||w g'||_1: solver-error mult.
    order: np.ndarray | None = None     # (n_patches, n_terms) observed ladder order
    rejected: int = 0                   # patches dropped as unresolved
    gram: np.ndarray | None = None      # (n_patches, n_terms, n_terms) noise Gram
    # THE MARTINGALE SCALE, when `build_nd(martingale=...)` declared one: the
    # measured int phi^2 nu^2 d[field] per patch, its estimator's own sd, and the
    # share of it a declared sigma_obs explains. It sits beside `noise_l2` because
    # it is the same KIND of object -- a per-patch scale a coverage factor
    # multiplies -- and apart from it because that one multiplies a DECLARED sigma
    # while this one is measured (docs/STOCHASTIC_CHECKER.md §1c).
    qv: np.ndarray | None = None        # (n_patches,)
    qv_se: np.ndarray | None = None     # (n_patches,)
    qv_obs_share: np.ndarray | None = None
    # per-patch record of the stride decomposition, when it ran: which scale was
    # used (martingale part or the total fallback) and by how much it tightened
    qv_decompose: list | None = None

    def normalize(self, by: str = "1", drop: bool = True) -> "WeakSystem":
        """Divide every row by its own `by` integral (default ∫φ), and drop that
        column.

        Rows from different patch SCALES are not commensurable: ∫φ differs, so a
        source term contributes a patch-dependent amount and the `1` column is a
        real, varying input. Handing it to a general library is a modelling
        error, and a measured one -- the library promptly built laws like
        u_xx·[1]^(3/2) and u_xx/[1], and 21 materially different classes
        certified on clean heat data, because `[1]` takes one value per scale and
        any function hitting the right value at three points fits.

        Normalizing removes the scale from the feature space entirely: every
        column becomes a patch AVERAGE, the linear relation is unchanged (row
        scaling cancels), a source term becomes a plain intercept the engine
        already proposes, and the bands scale with the rows they band (the Gram
        by the square, since the noise sensitivity vector scales linearly).
        """
        j = self.names.index(by)
        w = self.A[:, j].astype(float)
        if not np.all(np.isfinite(w)) or np.any(w == 0):
            raise ValueError(f"cannot normalize by {by!r}: zero/non-finite")
        keep = [k for k in range(len(self.names)) if not (drop and k == j)]
        return WeakSystem(
            A=(self.A / w[:, None])[:, keep],
            names=[self.names[k] for k in keep],
            patches=self.patches,
            quad=(self.quad / np.abs(w)[:, None])[:, keep],
            roundoff=(self.roundoff / np.abs(w)[:, None])[:, keep],
            noise_l2=(self.noise_l2 / np.abs(w)[:, None])[:, keep],
            field_l1=(self.field_l1 / np.abs(w)[:, None])[:, keep]
            if self.field_l1 is not None else None,
            order=self.order[:, keep],
            rejected=self.rejected,
            gram=(self.gram / (w ** 2)[:, None, None])[:, keep][:, :, keep]
            if self.gram is not None else None,
            # the martingale VARIANCE scales by the square of the row scaling,
            # exactly as the Gram does -- it is a variance, not an amplitude
            qv=self.qv / w ** 2 if self.qv is not None else None,
            qv_se=self.qv_se / w ** 2 if self.qv_se is not None else None,
            qv_obs_share=self.qv_obs_share)

    def det(self, field_err: float = 0.0) -> np.ndarray:
        """The DETERMINISTIC per-entry error: quadrature + roundoff, plus a
        declared field error (a reference solver's own tolerance-ladder bound)
        propagated through the L1 sensitivity.

        L1, not L2: a solver error is one fixed function, not a realization, so
        nothing cancels across the patch and the honest propagation is
        sum_i |w_i g'(u)_i| * delta. First order in delta, and declared as such."""
        d = self.quad + self.roundoff
        if field_err > 0:
            if self.field_l1 is None:
                raise ValueError("no field_l1: rebuild with the factory that "
                                 "records it before declaring a field error")
            d = d + field_err * self.field_l1
        return d

    def signal_to_band(self, target: str, **kw) -> float:
        """median |target| / declared band: how much evidence a patch family
        actually carries. A family near 1 certifies nothing worth having; the
        runner reports this rather than hiding it inside a filter."""
        y = self.A[:, self.names.index(target)]
        return float(np.median(np.abs(y) / self.declared_epsilon(target, **kw)))

    def declared_epsilon(self, target: str, *, coeff_max: float = 10.0,
                         sigma: float = 0.0,
                         field_err: float = 0.0) -> np.ndarray:
        """Per-patch declared band for the TARGET residual.

        The certified relation is y = sum_k c_k X_k, so the residual carries the
        target column's own error plus the feature columns' errors weighted by
        the coefficients. The coefficients are not known when eps is assembled,
        so they are bounded by a DECLARED `coeff_max` -- which the runner must
        check against the certified law and re-assemble if it was too small
        (an assumption stated and verified, never assumed silently)."""
        j = self.names.index(target)
        err = self.det(field_err)
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

    For a SYSTEM the vector nu spans every field (independent per-field noise
    makes the Gram block-diagonal over fields, and the blocks are summed when the
    Gram is built), so nothing here changes: it is still one quadratic form.
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
    feat_names: list | None = None   # feature order, when it is not "names minus target"

    def __post_init__(self):
        self.j = self.names.index(self.target)
        self.feat = ([n for n in self.names if n != self.target]
                     if self.feat_names is None else list(self.feat_names))
        self.cols = [self.names.index(n) for n in self.feat]
        self.syms = [sp.Symbol(f"x_{i}") for i in range(len(self.feat))]

    def subset(self, idx) -> "PatchEpsilon":
        idx = np.asarray(idx, int)
        return PatchEpsilon(self.names, self.target, self.y[idx], self.X[idx],
                            self.det[idx], self.gram[idx], self.sigma,
                            self.floor_abs, self.coeff_max, self.feat)

    def coefficients(self, expr):
        """(n_rows, n_feat) of df/dX_k, or None when the law cannot be read that
        way (non-finite, or not a function of the columns)."""
        if expr is None:
            return None
        try:
            e = sp.sympify(expr)
            if e.free_symbols - set(self.syms):
                return None
            # FAST PATH for a law that is LINEAR in the columns, which every
            # weak-form PDE law is: the gradient is the coefficient, a constant,
            # so there is nothing to differentiate or lambdify. The general path
            # costs one sympy.diff + one lambdify PER FEATURE PER CANDIDATE, and
            # with an exhaustive linear-basis search over 13 declared columns
            # that is ~8000 candidates x 13 -- measured as the dominant cost of
            # the PDEBench CFD run, above the weak-form integration itself.
            poly = e.as_poly(*self.syms) if self.syms else None
            if poly is not None and poly.total_degree() <= 1:
                c = np.array([float(e.coeff(s)) for s in self.syms])
                C = np.broadcast_to(c, (len(self.y), len(self.syms)))
                return np.array(C) if np.all(np.isfinite(C)) else None
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
    return make_patches_nd((x, t), (nx_half, nt_half), (n_x, n_t))


def make_patches_nd(coords, halves, counts) -> list[Patch]:
    """n-D patch grid: `coords` are the axis coordinate vectors (space..., time
    last), `halves` the half-widths in grid cells per axis, `counts` how many
    centres per axis. Supports lie strictly inside the domain on every axis."""
    coords = [np.asarray(c, float) for c in coords]
    halves = [int(h) for h in halves]
    centres = [np.linspace(h, len(c) - 1 - h, n).astype(int)
               for c, h, n in zip(coords, halves, counts)]
    steps = [float(c[1] - c[0]) for c in coords]
    out = []
    for combo in _product(centres):
        out.append(Patch(
            centers=tuple(float(c[i]) for c, i in zip(coords, combo)),
            halfwidths=tuple(h * s for h, s in zip(halves, steps)),
            idx=tuple(slice(i - h, i + h + 1) for i, h in zip(combo, halves))))
    return out


def onesided_patches(x: np.ndarray, t: np.ndarray, *, nx_half: int,
                     nt_cells: int, n_x: int, t0_index: int = 0) -> list[Patch]:
    """Patches whose TIME window starts at t[t0_index] and runs forward, with a
    test function equal to 1 there.

    The by-parts identity in time then keeps its boundary term, and that term is
    a linear functional of the INITIAL CONDITION with an analytic weight -- the
    whole mechanism of a state certificate (docs/DIRECTION_PDE_STATE.md). Space
    windows are the usual two-sided bumps, so nothing leaks at the spatial
    boundary."""
    x = np.asarray(x, float)
    t = np.asarray(t, float)
    hx = float(x[1] - x[0])
    ht = float(t[1] - t[0])
    xs = np.linspace(nx_half, len(x) - 1 - nx_half, n_x).astype(int)
    out = []
    if t0_index + nt_cells > len(t) - 1:
        # the window would run past the observations. A patch whose support is
        # not covered by data is not a patch -- silently truncating the slice
        # would leave the weight array claiming a window the field never fills.
        return out
    for i in xs:
        out.append(Patch(centers=(float(x[i]), float(t[t0_index])),
                         halfwidths=(nx_half * hx, nt_cells * ht),
                         idx=(slice(i - nx_half, i + nx_half + 1),
                              slice(t0_index, t0_index + nt_cells + 1)),
                         onesided=True))
    return out


def multiscale_patches(x: np.ndarray, t: np.ndarray, scales, *, n_x: int,
                       n_t: int) -> list[Patch]:
    """Patches at several supports pooled into one family.

    `scales` is a list of (nx_half, nt_half) in grid cells. Two reasons this is
    the default rather than a refinement:

    * a SINGLE-scale family makes any constant-integrand column exactly constant
      across rows -- the `1` term integrates to the same number on every patch --
      and the constrained-input detector then correctly reads that as a
      machine-exact input constraint and switches the engine to its
      domain-restricted path. The degeneracy is an artifact of the patch family,
      not of the physics (measured in the C1 ladder).
    * the test-function scale is a knob the law must be INDEPENDENT of; pooling
      scales puts that independence inside the certificate's own domain rather
      than leaving it to a separate patch-family sweep.
    """
    out = []
    for nxh, nth in scales:
        out += make_patches(x, t, nx_half=int(nxh), nt_half=int(nth),
                            n_x=n_x, n_t=n_t)
    return out


def multiscale_patches_nd(coords, scales, counts) -> list[Patch]:
    """`multiscale_patches` for n-D geometry: each entry of `scales` is a full
    per-axis tuple of half-widths in grid cells."""
    out = []
    for half in scales:
        out += make_patches_nd(coords, half, counts)
    return out


def _product(vectors):
    """Cartesian product of index vectors, as tuples."""
    out = [()]
    for v in vectors:
        out = [c + (int(i),) for c in out for i in v]
    return out


def _windows(patch: Patch, coords, dpsi, dpsi_one, alpha):
    """Per-axis analytic weight vectors (d^a psi)(s)/halfwidth^a, with the time
    axis one-sided (and trapezoid-weighted) when the patch says so."""
    ws = []
    for a, (c, cen, hw) in enumerate(zip(coords, patch.centers,
                                         patch.halfwidths)):
        k = alpha[a]
        last = a == len(coords) - 1
        d = dpsi_one if (patch.onesided and last) else dpsi
        w = d[k]((np.asarray(c, float) - cen) / hw) / hw ** k
        if patch.onesided and last:
            # A window that does NOT vanish at its left edge loses the
            # Euler-Maclaurin cancellation that makes the two-sided rule
            # spectrally accurate, so the rule there is composite trapezoid
            # (half weight at the endpoints) and the value is Richardson-
            # extrapolated across the same ladder that declares its error
            # (see _romberg_weights). Straight rectangle summation is O(h) here
            # and would put a quadrature floor two decades above the noise.
            w = w.copy()
            w[0] *= 0.5
            w[-1] *= 0.5
        ws.append(w)
    return ws


def _tensor(vs):
    out = vs[0]
    for v in vs[1:]:
        out = out[..., None] * v
    return out


def _weights(term: Term, coords, patch: Patch, dpsi, dpsi_one=None
             ) -> np.ndarray:
    """The analytic weight array (-1)^|alpha| (d^alpha phi)(z_i) * prod(h_a)."""
    ndim = len(coords)
    alpha = term.multi(ndim)
    ws = _windows(patch, coords, dpsi, dpsi_one, alpha)
    h = 1.0
    for c in coords:
        c = np.asarray(c, float)
        h *= float(c[1] - c[0])
    return _tensor(ws) * h * term.sign


def _embed(term: Term, coords, patch: Patch, dpsi, dpsi_one,
           steps) -> np.ndarray:
    """The weight array of the rule that subsamples axis a by steps[a],
    embedded (zero-padded) on the FULL patch grid."""
    shape = tuple(len(np.asarray(c, float)) for c in coords)
    W = np.zeros(shape)
    sub = tuple(slice(None, None, s) for s in steps)
    W[sub] = _weights(term, [np.asarray(c)[::s] for c, s in zip(coords, steps)],
                      patch, dpsi, dpsi_one)
    return W


def _romberg_weights(term: Term, coords, patch: Patch, dpsi, dpsi_one,
                     step: int, xstep: int = 1) -> np.ndarray:
    """Weights of the Richardson-extrapolated (Romberg) rule in the TIME axis,
    embedded on the full patch grid.

    E = I_h + (I_h - I_2h) / 3 cancels the h^2 endpoint term the one-sided
    window leaves behind. Returning it as one weight ARRAY (rather than
    combining computed integrals) is what keeps the noise vector, the Gram, the
    L1 sensitivity and the roundoff bound all referring to the SAME functional
    as the reported value -- the band must band the number that was actually
    reported.

    ONLY the time axis is refined, and that is the whole point: the space
    windows vanish at both ends, so their rule is spectrally accurate at h and
    garbage at 4h (a 21-point window subsampled to 6 points is not a
    quadrature). Extrapolating across a joint refinement mixes the two and
    CORRECTS with the junk -- measured: it moved the u_t integral by 1.3% and
    every one-sided patch was rejected as unresolved. The space direction gets
    its own separate, non-extrapolated bound in build_nd.
    """
    nd = len(coords)
    sx = (xstep,) * (nd - 1)
    Wa = _embed(term, coords, patch, dpsi, dpsi_one, sx + (step,))
    Wb = _embed(term, coords, patch, dpsi, dpsi_one, sx + (2 * step,))
    return Wa + (Wa - Wb) / 3.0


def aliasing_ratio(u_patch: np.ndarray, patch: Patch, coords, dpsi,
                   sigma: float = 0.0, dpsi_one=None) -> float:
    """Fraction of the windowed patch field's spectral energy sitting in the top
    half of the representable band.

    The refinement ladder alone can be FOOLED (measured): a field at ~4 points
    per wavelength aliases to a different low frequency at every coarse level,
    so the three integrals agree smoothly and report a high convergence order
    while all three are garbage. This is the direct test -- energy near Nyquist
    means the grid cannot represent the field, let alone its derivatives.
    Windowing by the test function first is what makes the DFT meaningful on a
    non-periodic patch (psi*u is smooth and compactly supported).

    A ONE-SIDED patch is diagnosed with its own two-sided window instead: that
    window does not taper at t0, so the windowed signal has a step there and its
    spectrum is broad no matter how well resolved the FIELD is (measured: 8e-3
    against a 1e-3 bar, rejecting every state patch). The resolution test is a
    statement about the field, so it gets a tapered window centred in the same
    time span -- the test function's own shape is not evidence about the grid.
    """
    if patch.onesided:
        patch = Patch(centers=patch.centers[:-1]
                      + (patch.centers[-1] + 0.5 * patch.halfwidths[-1],),
                      halfwidths=patch.halfwidths[:-1]
                      + (0.5 * patch.halfwidths[-1],),
                      idx=patch.idx, onesided=False)
        dpsi_one = None
    zero = Term("w", 0, 0, "u", alpha=(0,) * len(coords))
    W = _tensor(_windows(patch, coords, dpsi, dpsi_one, zero.multi(len(coords))))
    F = np.abs(np.fft.fftn(W * u_patch)) ** 2
    mask = None
    for a in range(len(coords)):
        k = np.abs(np.fft.fftfreq(F.shape[a])) > 0.25
        shp = [1] * len(coords)
        shp[a] = -1
        m = k.reshape(shp)
        mask = m if mask is None else (mask | m)
    mask = np.broadcast_to(mask, F.shape)
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


# A weighted realized-QV estimator's own standard deviation, self-normalized.
# Var(sum w (du)^2) = sum w^2 Var((du)^2) = 2 sum w^2 (b^2 dt)^2, and
# E[(du)^4] = 3 (b^2 dt)^2, so the estimate below needs no knowledge of b -- the
# same property that makes the VALUE computable without modelling the diffusion.
def _qv_se(w: np.ndarray, du: np.ndarray) -> float:
    return float(np.sqrt(2.0 / 3.0 * np.sum(w ** 2 * du ** 4)))


def _qv_entry(W: np.ndarray, term: Term, fp: dict, names, duprod: np.ndarray,
              sigma_obs: float = 0.0, diagonal: bool = True) -> tuple:
    """(value, |contributions|, se, observation-explained share) for a d[] term.

    `duprod` is the increment PRODUCT the measure integrates against: (du)^2 for
    d[u], and du_i * du_k for the cross measure d[u,v]. The weight array W carries
    NO grid spacing (see `Term.measure`), and both it and g are truncated to the
    increments' length.

    The estimator's variance differs between the two cases and both are stated
    rather than shared. For jointly Gaussian increments with variances v_i, v_k and
    covariance c, Var(du_i du_k) = v_i v_k + c^2 while
    E[(du_i du_k)^2] = v_i v_k + 2 c^2 -- so the fourth-moment sum is an UPPER bound
    on the variance in the cross case and is used as is (conservative), while on the
    diagonal c = v and the 2/3 factor makes it exact.

    `sigma_obs` DEBIASES the estimator, and it is not optional bookkeeping: with a
    field observed as u + e every increment carries Var 2 sigma_obs^2 on top of
    b^2 dt, so E[sum w (du_obs)^2] = sum w (b^2 dt + 2 sigma_obs^2) -- a bias that
    DIVERGES as dt -> 0. Measured as a confident-wrong when it was left in
    (docs/CASE_STUDY_STOCHASTIC_L0.md): the contamination is conservative wherever
    this quantity sets a BAND and a systematic offset wherever it sits on a
    TARGET, which is the same one-quantity-two-consumers rule the error-provenance
    direction states.
    """
    w = (W * term.g(fp))[..., :-1]
    val = float(np.sum(w * duprod))
    fac = 2.0 / 3.0 if diagonal else 1.0
    se = float(np.sqrt(fac * np.sum(w ** 2 * duprod ** 2)))
    share = 0.0
    if sigma_obs > 0 and diagonal:
        # observation noise adds 2 sigma^2 to each DIAGONAL increment variance; on
        # the cross measure independent per-field errors add nothing in expectation,
        # so there is no bias to remove there (only extra variance, already in `se`)
        bias = float(np.sum(w) * 2.0 * sigma_obs ** 2)
        share = abs(bias) / max(abs(val), 1e-300)
        val = val - bias
        se = float(np.hypot(se, np.sqrt(3.0 * np.sum(w ** 2))
                            * 2.0 * sigma_obs ** 2))
    return val, float(np.abs(w * duprod).sum()), se, share


# The martingale part of a realized quadratic variation is separable from the
# SMOOTH part, and measurably so -- which matters because a component that carries
# no noise still has a nonzero realized QV (the O(dt) residue of a differentiable
# path) and banding with it over-declares by ~sqrt(1/dt). Measured on Van der Pol's
# noise-free x component: residual 1.8e-4 against a band of 0.35.
#
# The separation is the STRIDE SCALING. Summing over all offsets at lag s,
#     sum_i (u[i+s] - u[i])^2  ~=  alpha * s  +  beta * s^2
# because a martingale increment's variance grows like s while a differentiable
# path's increment grows like s (so its square like s^2). alpha is therefore the
# martingale part and beta the smooth residue. Measured to within 0.3% of the truth
# on mixed paths, and 2000x below the total for a purely smooth one.
#
# Two safety rules, because using this SHRINKS a band:
#   * the returned scale is alpha + LAM_QV * se(alpha), an UPPER estimate, never the
#     point estimate;
#   * a poor two-term fit FALLS BACK to the total quadratic variation, which is
#     always sound. The two-term model is an approximation (a curved path
#     contributes O(s^3), a state-dependent b contributes O(s^1.5)), so the fallback
#     is what keeps the tightening from ever being a claim the data cannot support.
QV_STRIDES = (1, 2, 4, 8)
QV_FIT_MAX_REL = 0.05          # fit residual, relative to the stride-1 value


def qv_martingale_part(vals: np.ndarray, strides=QV_STRIDES) -> tuple:
    """(scale, info) -- the MARTINGALE part of a stride-indexed QV sum.

    `vals[j]` is the weighted increment-square sum at `strides[j]`. Returns an upper
    estimate of the martingale part and an info dict recording whether the two-term
    fit was good enough to use; when it was not, `scale` is the stride-1 total.
    """
    S = np.asarray(strides, float)
    q = np.asarray(vals, float)
    total = float(q[0])
    if len(S) < 3 or total <= 0:
        return total, {"used": "total", "reason": "too few strides"}
    A = np.vstack([S, S ** 2]).T
    coef, *_ = np.linalg.lstsq(A, q, rcond=None)
    resid = float(np.max(np.abs(q - A @ coef)))
    alpha, beta = float(coef[0]), float(coef[1])
    # se(alpha) from the fit's own residual through the normal equations
    try:
        cov = np.linalg.inv(A.T @ A)
        dof = max(len(S) - 2, 1)
        s2 = float(np.sum((q - A @ coef) ** 2)) / dof
        se = float(np.sqrt(max(s2 * cov[0, 0], 0.0)))
    except np.linalg.LinAlgError:                              # pragma: no cover
        return total, {"used": "total", "reason": "singular stride design"}
    if resid > QV_FIT_MAX_REL * total or alpha < 0:
        return total, {"used": "total", "reason": "two-term fit poor",
                       "fit_residual_rel": resid / total, "alpha": alpha}
    up = alpha + KAPPA * se
    if up >= total:
        return total, {"used": "total", "reason": "no tightening available"}
    return up, {"used": "martingale part", "alpha": alpha, "beta": beta,
                "se_alpha": se, "scale": up, "total": total,
                "tightened_by": total / max(up, 1e-300),
                "fit_residual_rel": resid / total}


def _martingale_scale(W0: np.ndarray, nus, dus, sigma_obs: float = 0.0,
                      fields_for_strides=None) -> tuple:
    """(<M>, se, observation-explained share) for a row whose target functional is
    the stochastic integral sum_i int phi nu_i dX_i^mart.

    MULTI-FIELD, and the algebra collapses pleasingly. In d dimensions
    <M> = int phi^2 sum_ik nu_i nu_k d[X_i, X_k], which as a sum over increments is
    sum_n (phi_n sum_i nu_i,n du_i,n)^2 -- the SQUARE OF THE SUMMED increment, one
    scalar per sample. No cross terms to enumerate and no d^2 loop; the scalar case
    is the same formula with one field.

    Its variance is int phi^2 nu^2 d[u] -- QUADRATIC in the test function, which
    is why it is not a Term and never could be: the weak form is linear in phi.
    It belongs exactly where `noise_l2` belongs, and is the same kind of object:
    a per-patch scale the band multiplies by a coverage factor. The difference is
    that `noise_l2` multiplies a DECLARED sigma and this one is MEASURED, which is
    what lets `certify.coverage_factor` be applied to an intrinsic term without
    letting a candidate widen its own band.
    """
    xi = np.zeros_like(np.asarray(dus[0], float))
    for nu, du in zip(nus, dus):
        xi = xi + (W0 * nu)[..., :-1] * du
    xi2 = xi ** 2
    val = float(np.sum(xi2))
    se = float(np.sqrt(2.0 / 3.0 * np.sum(xi2 ** 2)))
    info = None
    if fields_for_strides is not None:
        # STRIDE DECOMPOSITION (opt-in): separate the martingale part of this scale
        # from the smooth residue a differentiable component contributes. See
        # qv_martingale_part -- it falls back to the total whenever the two-term fit
        # cannot support the tightening, so the direction is always safe.
        vals = []
        for sd in QV_STRIDES:
            xs = None
            for nu, u in zip(nus, fields_for_strides):
                inc = u[..., sd:] - u[..., :-sd]
                w = (W0 * nu)[..., :-sd]
                xs = w * inc if xs is None else xs + w * inc
            vals.append(float(np.sum(xs ** 2)))
        val, info = qv_martingale_part(vals)
        # The returned `val` is ALREADY alpha + kappa*se(alpha) -- an upper estimate
        # of the martingale part that carries its own estimation error. So the
        # downstream LAM_QV * qv_se inflation, which exists to cover exactly that
        # error, would double-count it, and `se` is zero here.
        #
        # The first version scaled the total's `se` down in proportion to the
        # tightening. That was arbitrary AND in the band-shrinking direction, which
        # is the pair of properties this arc keeps finding at the bottom of its
        # defects. When the fit falls back to the total, so does the error term.
        se = 0.0 if info.get("used") == "martingale part" else se
    share = 0.0
    if sigma_obs > 0:
        # each field's observation error inflates its own increment variance by
        # 2 sigma^2, and the sensitivities enter squared
        w2 = sum(float(np.sum(((W0 * nu)[..., :-1]) ** 2)) for nu in nus)
        bias = w2 * 2.0 * sigma_obs ** 2
        share = abs(bias) / max(abs(val), 1e-300)
        val = val - bias
        se = float(np.hypot(se, np.sqrt(3.0 * w2) * 2.0 * sigma_obs ** 2))
    return val, se, share, info


def _entry(W: np.ndarray, term: Term, fp: dict, names) -> tuple:
    """(value, sum|contributions|, ||w g'||_2, ||w g'||_1) for one weight array.

    The sensitivity norms sum over EVERY field the term reads: with independent
    per-field noise the sensitivity vectors simply concatenate, so the L2 norm
    is the root of the summed per-field squares (docs/DIRECTION_PDE_SYSTEMS.md).
    """
    C = W * term.g(fp)
    l2sq, l1 = 0.0, 0.0
    for f in names:
        d = W * term.dg(fp, f)
        l2sq += float(np.sum(d ** 2))
        l1 += float(np.sum(np.abs(d)))
    return (float(C.sum()), float(np.abs(C).sum()), float(np.sqrt(l2sq)), l1)


def _wnorm(W: np.ndarray, term: Term, fp: dict, names) -> float:
    """||W g'||_2 summed over fields: how much declared noise a difference
    functional with weight array W explains."""
    return float(np.sqrt(sum(np.sum((W * term.dg(fp, f)) ** 2)
                             for f in names)))


def _ladder_noise(term: Term, fp: dict, coords, patch: Patch, dpsi, dpsi_one,
                  sa: int, sb: int, names) -> float:
    """||d||_2 for the difference functional I_{sa*h} - I_{sb*h}: the ladder
    difference is ALSO a linear functional of the field noise, with a computable
    weight vector, so sigma * this is exactly how much of an observed ladder
    difference the declared noise explains."""
    shape = np.shape(next(iter(fp.values())))
    total = 0.0
    for f in names:
        d = np.zeros(shape)
        for s, sign in ((sa, 1.0), (sb, -1.0)):
            sub = tuple(slice(None, None, s) for _ in coords)
            sc = [np.asarray(c)[::s] for c in coords]
            fsub = {k: v[sub] for k, v in fp.items()}
            d[sub] += sign * (_weights(term, sc, patch, dpsi, dpsi_one)
                              * term.dg(fsub, f))
        total += float((d ** 2).sum())
    return float(np.sqrt(total))


def _noise_gram(terms, Ws, fp: dict, names) -> np.ndarray:
    """Gram matrix of the per-term noise SENSITIVITY vectors on this patch:
    G[k,l] = sum_fields sum_i (W_k dg_k/df)_i (W_l dg_l/df)_i.

    To first order a field perturbation e moves term k's integral by
    <W_k g_k'(u), e>, so the residual of y = sum_k c_k X_k moves by <v, e> with
    v = nu_target - sum_k c_k nu_k -- ONE vector, because every column is a
    functional of the SAME realization. Its length is a'Ga with
    a = (1, -c_1, ..., -c_K), which is why the Gram is what gets stored: the
    band for any candidate is then a quadratic form costing K^2 flops per patch
    instead of a pass over the grid. With several fields the vector simply
    concatenates over them (independent per-field noise), which is this sum.
    """
    G = np.zeros((len(terms), len(terms)))
    for f in names:
        nu = np.stack([(W * tm.dg(fp, f)).ravel() for tm, W in zip(terms, Ws)])
        G += nu @ nu.T
    return G


def build(u, x: np.ndarray, t: np.ndarray, terms: list, patches: list[Patch],
          *, p: int = 8, sigma: float = 0.0) -> WeakSystem:
    """The 1-D+t factory: field(s) on grid (x[i], t[j]) -> patch design matrix
    with a declared error bound per entry.

    `u` is either an array (the single field "u") or a dict of named fields
    sharing the grid -- a system is one row set over several fields.
    """
    fields, _ = _as_fields(u)
    return build_nd(fields, (x, t), terms, patches, p=p, sigma=sigma)


def build_nd(fields, coords, terms: list, patches: list[Patch], *,
             p: int = 8, sigma: float = 0.0, rough: bool = False,
             martingale: tuple | None = None, sigma_obs: float = 0.0,
             martingale_decompose: bool = False) -> WeakSystem:
    """The factory over arbitrary geometry: `coords` are the axis coordinate
    vectors (space..., TIME LAST) and every field array is shaped accordingly.

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

    `sigma` is ONE declared field-noise scale for every field. Per-field scales
    that differ must be declared at the largest (conservative) -- an unequal
    declaration would need per-field Gram blocks and is not claimed here.

    Three arguments serve the Itô weak form (docs/DIRECTION_STOCHASTIC.md step 3)
    and are inert at their defaults, so no campaign predating them changes:

    `rough=True` SKIPS THE ALIASING GATE, and it has to. That gate drops a patch
    whose windowed field has energy near Nyquist, because a smooth field with such
    energy is under-resolved. A Brownian-driven path has energy at every frequency
    BY CONSTRUCTION -- the gate would reject every patch, and it would be answering
    the wrong question. The refinement ladder still applies (the path is fixed
    data, so subsampling it is a genuine coarser rule) and it is what decides
    whether the quadrature converged.

    `martingale=(field, gexpr)` declares that the target functional of these rows
    is the stochastic integral int phi (gexpr) dW driven by `field` -- so its
    variance int phi^2 gexpr^2 d[field] can be MEASURED per patch and returned as
    `qv`/`qv_se`. It is a declaration about the STRUCTURE of the row, not about any
    magnitude: the magnitude comes from the data.

    `sigma_obs` is the declared per-sample measurement error on the fields, and it
    debiases every quadratic-variation quantity (see `_qv_entry`).
    """
    fields, shape = _as_fields(fields)
    coords = [np.asarray(c, float) for c in coords]
    ndim = len(coords)
    T = [LIBRARY[n] if isinstance(n, str) else n for n in terms]
    fnames = sorted(fields)
    for tm in T:
        missing = set(tm.fields) - set(fnames)
        if missing:
            raise ValueError(f"term {tm.name!r} reads unknown field(s) "
                             f"{sorted(missing)}; have {fnames}")
    maxo = max(max(tm.multi(ndim)) for tm in T)
    dpsi = bump_derivatives(p, maxo)
    dpsi_one = (bump_derivatives(p, maxo, onesided=True)
                if any(pa.onesided for pa in patches) else None)
    rows, quad, rnd, l2, l1s, orders, kept, grams = [], [], [], [], [], [], [], []
    qvs, qv_ses, qv_shares, qv_decomp = [], [], [], []
    rejected = 0
    for pa in patches:
        sub = tuple(pa.idx)
        cs = [c[s] for c, s in zip(coords, sub)]
        fp = {k: v[sub] for k, v in fields.items()}
        vals, qs, rs, ls, l1v, ords, Ws = [], [], [], [], [], [], []
        ok = rough or all(aliasing_ratio(fp[f], pa, cs, dpsi, sigma, dpsi_one)
                          <= ALIAS_MAX for f in fnames)
        if not ok:
            rejected += 1
            continue
        for tm in T:
            if pa.onesided:
                # The reported value, its weight array and its declared bound
                # all come from the Romberg-in-time rule (see
                # _romberg_weights). TWO separate truncation bounds are declared
                # and summed, because the two directions converge for different
                # reasons: the time bound is the disagreement between the two
                # extrapolants, the space bound is the plain h-vs-2h difference
                # in x alone (over-covering, since that direction converges
                # spectrally).
                W1 = _romberg_weights(tm, cs, pa, dpsi, dpsi_one, 1)
                W2 = _romberg_weights(tm, cs, pa, dpsi, dpsi_one, 2)
                W4 = _romberg_weights(tm, cs, pa, dpsi, dpsi_one, 4)
                Wx = _romberg_weights(tm, cs, pa, dpsi, dpsi_one, 1, xstep=2)
                v_h, absum, w2, wl1 = _entry(W1, tm, fp, fnames)
                v_2, _, _, _ = _entry(W2, tm, fp, fnames)
                v_4, _, _, _ = _entry(W4, tm, fp, fnames)
                v_x, _, _, _ = _entry(Wx, tm, fp, fnames)
                round_h = MACHINE_EPS * absum
                d1, d2, dx_ = abs(v_h - v_2), abs(v_2 - v_4), abs(v_h - v_x)
                ns1 = ns2 = 0.0
                if sigma > 0:
                    ns1 = sigma * _wnorm(W1 - W2, tm, fp, fnames)
                    ns2 = sigma * _wnorm(W1 - Wx, tm, fp, fnames)
                    d1 = max(0.0, d1 - KAPPA * ns1)
                    d2 = max(0.0, d2 - KAPPA * sigma
                             * _wnorm(W2 - W4, tm, fp, fnames))
                    dx_ = max(0.0, dx_ - KAPPA * ns2)
                # SAME asymptotic-regime discipline as the two-sided ladder: the
                # extrapolant is 4th order, so |e_h| <= d1 / (2^r - 1) once
                # refinement really buys accuracy at the observed order r; below
                # order 2 no Richardson estimate would be honest and the raw
                # difference stands.
                r_t = (np.log2(d2 / d1) if d1 > 0 and d2 > 0 else 0.0)
                if np.isfinite(r_t) and r_t >= MIN_LADDER_ORDER:
                    r_used = float(np.clip(r_t, MIN_LADDER_ORDER,
                                           MAX_TRUSTED_ORDER))
                    q_t = d1 / (2.0 ** r_used - 1.0)
                else:
                    q_t = d1
                q = max(q_t + dx_, round_h, ns1, ns2)
                if q > UNRESOLVED_REL * max(absum, MACHINE_EPS):
                    ok = False
                    break
                vals.append(v_h); qs.append(q); rs.append(round_h)
                ls.append(w2); l1v.append(wl1)
                ords.append(float(r_t))
                Ws.append(W1)
                continue
            if tm.qv_fields is not None:
                # A d[] TERM DOES NOT REFINE. Subsampling a realized-quadratic-
                # variation estimator does not evaluate the same integral more
                # accurately -- it computes a different estimator, over half the
                # increments -- so the ladder would report a "truncation error"
                # that is really the estimator's own statistical spread. That
                # spread is what `_qv_entry` returns, and it is declared as the
                # term's bound directly.
                W0 = _tensor(_windows(pa, cs, dpsi, dpsi_one,
                                      tm.multi(ndim))) * tm.sign
                fa, fb = tm.qv_fields
                duprod = np.diff(fp[fa], axis=-1) * np.diff(fp[fb], axis=-1)
                v_h, absum, se, share = _qv_entry(W0, tm, fp, fnames, duprod,
                                                  sigma_obs, diagonal=fa == fb)
                vals.append(v_h)
                qs.append(max(se, MACHINE_EPS * absum))
                rs.append(MACHINE_EPS * absum)
                # a d[] term carries no first-order sensitivity to a declared
                # field noise the way a dt term does: its contamination is the
                # 2 sigma_obs^2 bias `_qv_entry` has already removed, and what is
                # left sits in `se`. Zero here keeps the L1/L2 channels from
                # double-counting it.
                ls.append(0.0)
                l1v.append(0.0)
                ords.append(float("nan"))
                Ws.append(np.zeros_like(W0))
                continue
            ladder = []
            for step in (1, 2, 4):
                sc = [c[::step] for c in cs]
                fsub = {k: v[tuple(slice(None, None, step) for _ in coords)]
                        for k, v in fp.items()}
                W = _weights(tm, sc, pa, dpsi, dpsi_one)
                ladder.append((_entry(W, tm, fsub, fnames), W))
            (v_h, absum, w2, wl1), W_h = ladder[0]
            d1 = abs(ladder[1][0][0] - ladder[0][0][0])    # |I_2h - I_h|
            d2 = abs(ladder[2][0][0] - ladder[1][0][0])    # |I_4h - I_2h|
            scale = max(absum, MACHINE_EPS)
            round_h = MACHINE_EPS * absum
            # Under declared noise the ladder differences are dominated by the
            # NOISE, not by truncation, and the raw convergence test then rejects
            # patches for carrying noise the eps model already bands (measured:
            # 16 of 24 heat patches lost at sigma = 1e-4). Subtract what sigma
            # explains before asking whether refinement bought accuracy.
            ns1 = ns2 = 0.0
            if sigma > 0:
                ns1 = sigma * _ladder_noise(tm, fp, cs, pa, dpsi, dpsi_one,
                                            1, 2, fnames)
                ns2 = sigma * _ladder_noise(tm, fp, cs, pa, dpsi, dpsi_one,
                                            2, 4, fnames)
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
                    if not rough:
                        ok = False
                        break
                    # A ROUGH PATH IS NOT EXPECTED TO BE IN THE ASYMPTOTIC REGIME:
                    # the integrand is Holder-1/2, the quadrature converges at
                    # O(h) rather than spectrally, and the observed order sits near
                    # 1. Dropping the patch would drop every patch. The RAW ladder
                    # difference stands as the declared bound instead -- strictly
                    # larger than any Richardson estimate, so the direction is
                    # conservative -- and whether it is small enough to use is the
                    # CALLER's decision, because only the caller knows the coverage
                    # factor the martingale band will carry.
                    vals.append(v_h)
                    qs.append(max(d1, round_h))
                    rs.append(round_h)
                    ls.append(w2)
                    l1v.append(wl1)
                    ords.append(float(r))
                    Ws.append(W_h)
                    continue
                # |e_h| <= d1 / (2^r - 1) holds when the error really shrinks by
                # 2^r per refinement. The ORDER USED IS CAPPED (measured reason:
                # for the exponentially-convergent bump quadrature the coarse
                # (4h, 2h) pair reports r ~ 10 while the true h-level ratio is
                # ~2^7 -- an overstated order under-declares the band and
                # rejects true laws). The observed order is recorded separately.
                r_used = float(np.clip(r, MIN_LADDER_ORDER, MAX_TRUSTED_ORDER))
                q = d1 / (2.0 ** r_used - 1.0)
                if q > UNRESOLVED_REL * scale and not rough:
                    ok = False                 # resolved, but too coarsely to use
                    break
                # rough: the bound stands and the caller gates on it. A relative
                # bar against the term's OWN scale is the wrong test here -- the
                # martingale, not the quadrature, is what the band is made of, so
                # the comparison that matters is against the martingale scale and
                # the caller is the one holding it.
            vals.append(v_h)
            qs.append(max(q, round_h))
            rs.append(round_h)
            ls.append(w2)
            l1v.append(wl1)
            ords.append(float(r))
            Ws.append(W_h)
        if not ok:
            rejected += 1
            continue
        if martingale is not None:
            # one (field, gexpr) pair, or a LIST of them for a multi-field state:
            # the target's sensitivity to each driving field in turn
            pairs = ([martingale] if isinstance(martingale[0], str)
                     else list(martingale))
            W0 = _tensor(_windows(pa, cs, dpsi, dpsi_one, (0,) * ndim))
            nus = [Term(f"nu{i}", gexpr=mg, alpha=(0,) * ndim).g(fp)
                   for i, (_, mg) in enumerate(pairs)]
            dus = [np.diff(fp[mf], axis=-1) for mf, _ in pairs]
            v, se, share, dec = _martingale_scale(
                W0, nus, dus, sigma_obs,
                fields_for_strides=[fp[mf] for mf, _ in pairs]
                if martingale_decompose else None)
            qvs.append(v)
            qv_ses.append(se)
            qv_shares.append(share)
            if dec is not None:
                qv_decomp.append(dec)
        rows.append(vals)
        quad.append(qs)
        rnd.append(rs)
        l2.append(ls)
        l1s.append(l1v)
        orders.append(ords)
        kept.append(pa)
        grams.append(_noise_gram(T, Ws, fp, fnames))
    return WeakSystem(A=np.array(rows, float), names=[tm.name for tm in T],
                      patches=kept, quad=np.array(quad, float),
                      roundoff=np.array(rnd, float),
                      noise_l2=np.array(l2, float),
                      field_l1=np.array(l1s, float),
                      order=np.array(orders, float), rejected=rejected,
                      gram=np.array(grams, float) if grams else None,
                      qv=np.array(qvs, float) if qvs else None,
                      qv_se=np.array(qv_ses, float) if qv_ses else None,
                      qv_obs_share=np.array(qv_shares, float) if qv_shares
                      else None,
                      qv_decompose=qv_decomp or None)


# --------------------------------------------------------------------------
# Initial-condition columns (state certificates, docs/DIRECTION_PDE_STATE.md)
# --------------------------------------------------------------------------

def ic_columns(basis, x: np.ndarray, patches: list[Patch], *, p: int = 16,
               nodes: int = 400) -> tuple:
    """The boundary columns B[i, j] = int phi_i(x, t0) b_j(x) dx, and their
    declared quadrature error.

    These are the pleasant asymmetry of a state certificate: phi is analytic and
    the basis is declared, so the columns carry NO DATA ERROR AT ALL -- only the
    quadrature error of an integral of two known smooth functions, which
    Gauss-Legendre drives to the float floor. They are computed OFF the data
    grid on purpose: nothing about the measurement enters them.

    Returns (B, err) with err[j] the per-column bound, measured by doubling the
    node count (a convergence check, never an assumed exactness).
    """
    psi = bump_derivatives(p, 0, onesided=False)[0]

    def quad(n):
        out = np.zeros((len(patches), len(basis)))
        z, w = np.polynomial.legendre.leggauss(n)
        for i, pa in enumerate(patches):
            xs = pa.centers[0] + pa.halfwidths[0] * z
            ww = w * pa.halfwidths[0] * psi(z)
            for j, b in enumerate(basis):
                out[i, j] = float(np.sum(ww * np.asarray(b(xs), float)))
        return out

    B = quad(nodes)
    Bc = quad(nodes // 2)
    err = np.max(np.abs(B - Bc), axis=0) + MACHINE_REL * np.max(np.abs(B), axis=0)
    return B, err

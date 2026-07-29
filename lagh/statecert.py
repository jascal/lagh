"""State certificates: inverting for an initial condition
(docs/DIRECTION_PDE_STATE.md, registration docs/CASE_STUDY_PDE_C4.md).

Everything this program has certified so far is a LAW over a stated domain.
Recovering an initial condition is a claim about one system's particular
history -- a different KIND of claim, kept in its own category deliberately
rather than slid into.

The mechanism is the boundary term the weak form has been throwing away. Every
test function used through C3 vanishes on the whole patch boundary, which is
what kills the boundary terms in the by-parts identity. Drop that IN TIME ONLY
-- a phi that is nonzero at the window's initial time -- and

    int int phi u_t  =  [int phi u dx]_{t0}^{T}  -  int int phi_t u

so with phi(., T) = 0 the initial condition enters as -int phi(x, t0) u(x, t0) dx:
a linear functional of the unknown, with an ANALYTIC weight, inside the same
linear system, under the same declared band. Rearranged, with the law KNOWN,

    int phi(., t0) u0 dx  =  [the u_t column]  -  sum_k c_k [the law's columns]
    <=>  sum_j a_j B_ij   =  y_i

where B_ij = int phi_i(x, t0) b_j(x) dx are known analytic constants carrying NO
data error at all, and y_i carries the ordinary weak-form band. The pleasant
asymmetry: the band is SIMPLER here than in the law case, not harder -- the
design matrix is exact and only the target is noisy.

What a state certificate says, and how it differs from a law certificate:

* CLAIM. Over the stated observation window and patch family, every initial
  condition whose mode coefficients lie in the reported intervals reproduces the
  observations within the declared band. Modes outside the reported set are NOT
  determined -- they are named, not silently dropped. An ill-posed inversion's
  honest output is a resolution statement.
* DOMAIN. The observation window and the basis. It says nothing about other
  times, other solutions, or modes above the reported cut.
* ALPHA. With a fixed basis and a known law there is no search: |H| = 1, so
  alpha is a pure chance-agreement bound rather than a multiple-comparisons
  corrected one. Defensible, but a DIFFERENT quantity from a law certificate's
  alpha, and labelled as such (`alpha_kind`).
* RESOLUTION BOUND, stated up front: dof = number of basis modes, so
  h = n_rows - dof must stay positive with margin. At most as many modes can be
  certified as there are independent patch equations -- the same mechanism as
  the 35-term-interpolation lesson, and it belongs in the registration rather
  than in the results.
"""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
import sympy as sp

from .certify import (band, check, determination, parameter_interval,
                      significance_log10)
from .weakform import LIBRARY, PatchEpsilon, build, ic_columns

# h = n_rows - dof must clear this: a certificate resting on fewer independent
# patch equations than it has unknowns is an interpolation, and one resting on
# barely more is not much better. Stated up front, in the registration.
MIN_HELDOUT = 8
# The interval search calls a parameter UNDETERMINED when no bound is found
# within this relative range -- deliberately huge, so "undetermined" means
# genuinely unbounded rather than "wider than we liked".
MAX_REL = 1e6


def _scaled(B, y, eps):
    """Row-scale by the band and column-scale by column norm.

    Both scalings are exact changes of variable, and both are needed: the rows
    span a decade of band across patch scales and the columns span many decades
    (a wide test function's spectrum is tiny at high k), so the unscaled program
    is badly conditioned and the solver returns failure on exactly the hard
    cases -- measured on the shock rungs, where the LP failed, the code fell
    back to least squares, and a state whose TRUTH sat at 0.57 of its band came
    back as 'no state explains the observations'. A spurious refusal is still a
    wrong answer."""
    r = 1.0 / np.asarray(eps, float)
    Br = np.asarray(B, float) * r[:, None]
    d = np.sqrt(np.einsum("ij,ij->j", Br, Br))
    d = np.where(d > 0, d, 1.0)
    return Br / d, np.asarray(y, float) * r, d


def _feasible_center(B, y, eps):
    """The most comfortably certifying amplitude vector: maximize the slack s in
    |B a - y| <= eps - s. A linear program, because that is exactly the
    certificate's question -- IS there an initial condition in the declared
    basis consistent with the observations? -- while least squares answers a
    different one and answers it badly here.

    Measured on advection at sigma = 0: the plain lstsq fit violated the band at
    9.8x on some rows while the TRUE amplitudes sat at 0.97x. The band is
    per-row and spans a decade across patch scales, so an L2 fit happily trades
    a loose row's slack for a tight row's violation. Returns (a, slack), or
    (None, None) when the LP fails to solve (the caller falls back to lstsq and
    the certification check still decides)."""
    from scipy.optimize import linprog
    Bs, ys, dscale = _scaled(B, y, eps)          # constraints become |Bs a - ys| <= 1
    n, d = Bs.shape
    A_ub = np.vstack([np.hstack([Bs, np.ones((n, 1))]),
                      np.hstack([-Bs, np.ones((n, 1))])])
    b_ub = np.concatenate([ys + 1.0, -ys + 1.0])
    c = np.zeros(d + 1)
    c[-1] = -1.0
    r = linprog(c, A_ub=A_ub, b_ub=b_ub,
                bounds=[(None, None)] * d + [(None, None)], method="highs")
    if not r.success:
        return None, None
    return np.asarray(r.x[:d], float) / dscale, float(r.x[-1])


def _mode_bounds(B, y, eps, j):
    """The exact JOINT interval for mode j: min and max of a_j over the whole
    feasible set {a : |B a - y| <= eps}, by two linear programs.

    This is the honest resolution statement, and it differs from bisecting one
    coefficient with the others HELD FIXED (certify.parameter_interval), which
    answers a conditional question and is necessarily narrower. Both are
    reported; where they differ the joint one is the claim.

    Returns (lo, hi), or None when the projection is unbounded -- which is
    exactly what UNDETERMINED means."""
    from scipy.optimize import linprog
    Bs, ys, dscale = _scaled(B, y, eps)
    n, d = Bs.shape
    A_ub = np.vstack([Bs, -Bs])
    b_ub = np.concatenate([ys + 1.0, -ys + 1.0])
    out = []
    for sgn in (1.0, -1.0):
        c = np.zeros(d)
        c[j] = sgn
        r = linprog(c, A_ub=A_ub, b_ub=b_ub, bounds=[(None, None)] * d,
                    method="highs")
        if r.status == 3 or not r.success:      # 3 = unbounded
            return None
        out.append(float(r.x[j]) / dscale[j])
    return (min(out), max(out))


@dataclass
class StateCertificate:
    certified: bool
    modes: dict = field(default_factory=dict)   # label -> value/interval/determined
    undetermined: list = field(default_factory=list)
    abstain: str | None = None
    alpha_log10: float | None = None
    alpha_kind: str = "state (|H| = 1: fixed basis, known law, no search)"
    n_rows: int = 0
    dof: int = 0
    heldout: int = 0
    window: tuple = ()
    basis: list = field(default_factory=list)
    notes: list = field(default_factory=list)
    residual_ratio: float | None = None         # max |resid| / band
    amp_max_declared: float | None = None
    amp_max_certified: float | None = None
    slack: float | None = None                  # LP margin inside the band
    # PARTIAL DETERMINATION in the shared vocabulary (`certify.determination`),
    # keyed by MODE LABEL. `modes` above stays as it is -- it carries the
    # conditional interval and the fitted value, which the shared record does
    # not model -- but the determined/undetermined/resolved statement now reads
    # the same here as it does on a law certificate and a structural abstain,
    # which is the whole point of having one vocabulary. This is the MODE
    # dimension of the five.
    partial: dict = field(default_factory=dict)

    def one_line(self) -> str:
        if not self.certified:
            return f"ABSTAIN[{self.abstain}] rows={self.n_rows} dof={self.dof}"
        det = len(self.modes) - len(self.undetermined)
        return (f"STATE CERTIFIED: {det}/{len(self.modes)} modes determined, "
                f"{len(self.undetermined)} reported UNDETERMINED, "
                f"alpha <= 1e{self.alpha_log10:.0f}")


def fourier_basis(kmax: int, *, sine: bool = True) -> tuple:
    """The declared basis: {1, cos kx, sin kx} up to kmax, as (labels, callables).

    Fourier on a periodic domain is not an aesthetic choice: the interval search
    has to happen in a basis where the forward operator is DIAGONAL, or the
    per-mode question is not well posed. Anything non-periodic is a later
    problem, not a first one -- stated, not discovered."""
    labels = ["1"]
    fns = [lambda x: np.ones_like(np.asarray(x, float))]
    for k in range(1, kmax + 1):
        labels.append(f"cos{k}")
        fns.append((lambda k: lambda x: np.cos(k * np.asarray(x, float)))(k))
        if sine:
            labels.append(f"sin{k}")
            fns.append((lambda k: lambda x: np.sin(k * np.asarray(x, float)))(k))
    return labels, fns


def truth_check_state(B, y, eps, a_true) -> dict:
    """Does the TRUE state, expressed in the declared basis, sit inside its own
    band? The state-certificate twin of pdesystem.truth_check, and it separates
    two refusals that look identical from the outside:

      * the DECLARED BASIS cannot represent the state (a truncated Fourier basis
        against a Cole-Hopf profile, whose spectrum is infinite) -- a statement
        about the basis, fixable by declaring a bigger one,
      * the information really is not in the observations.

    Reading the first as the second is exactly the mistake the system probe's
    construction bug taught: check the truth against its own band first.
    """
    r = np.abs(B @ np.asarray(a_true, float) - y) / eps
    return {"truth_max_ratio": float(np.max(r)),
            "truth_median_ratio": float(np.median(r)),
            "truth_certifies": bool(np.max(r) <= 1.0)}


def backpropagate(cert: StateCertificate, kind: str, t0: float, *, nu=0.1,
                  c=0.7) -> dict:
    """Carry a certificate of the state at t0 back to the INITIAL condition at
    t = 0 through the known law, and report what that costs per mode.

    A one-sided window reads the state at the window's own start, so a window
    beginning at t = 0 sees every mode at full amplitude and its resolution is
    flat in k -- measured, and the reason the registered exp(+nu k^2 T) curve
    does NOT appear there. The exponential ill-posedness of the inverse problem
    is entirely in the propagation: for the heat equation the forward map is
    diagonal with gain d_k = exp(-nu k^2 t0), so the certified interval at t0
    maps to interval / d_k at time zero -- an exact transformation of the claim,
    since the map is diagonal and known.

    For advection the map is a rotation in each (cos k, sin k) pair with unit
    gain, so the widths are unchanged; that is reported rather than recomputed,
    and no per-mode box claim is made for it (a rotation does not map a box to a
    box).

    A mode is RESOLVED when its back-propagated interval EXCLUDES ZERO -- a
    threshold-free criterion: below it the certificate cannot even determine
    that the mode is present.
    """
    out = {"t0": float(t0), "kind": kind, "modes": {}, "k_cut": None}
    kmax_res = 0
    for lab, m in cert.modes.items():
        if not lab.startswith(("cos", "sin")):
            continue
        k = int(lab[3:])
        gain = np.exp(-nu * k ** 2 * t0) if kind == "heat" else 1.0
        iv = m["interval"]
        if iv is None:
            out["modes"][lab] = {"interval": None, "resolved": False,
                                 "gain": gain}
            continue
        lo, hi = iv[0] / gain, iv[1] / gain
        resolved = not (lo <= 0.0 <= hi)
        out["modes"][lab] = {"interval": [lo, hi],
                             "half_width": 0.5 * (hi - lo),
                             "value": m["value"] / gain,
                             "gain": gain, "resolved": bool(resolved)}
        if resolved:
            kmax_res = max(kmax_res, k)
    out["k_cut"] = kmax_res
    return out


def assemble_state(u, x, t, law: dict, patches, basis_fns, *, p: int = 16,
                   sigma: float = 0.0, field_err: float = 0.0,
                   amp_max: float = 10.0, target: str = "u_t"):
    """(B, y, eps, info): the state linear system and its declared band.

    `law` maps weak-form term names to their KNOWN coefficients (the law is an
    input here, not the unknown). `patches` must be one-sided in time
    (weakform.onesided_patches) or the boundary term this rests on is zero.
    """
    names = [target] + [k for k in law if k != target]
    terms = [LIBRARY[n] for n in names]
    s = build(u, x, t, terms, patches, p=p, sigma=sigma)
    info = {"rejected": int(s.rejected), "n_patches": int(len(s.A))}
    if len(s.A) == 0:
        return None, None, None, info
    j = s.names.index(target)
    feat = [n for n in s.names if n != target]
    cols = [s.names.index(n) for n in feat]
    det = s.det(field_err)
    m = PatchEpsilon(s.names, target, s.A[:, j], s.A[:, cols], det, s.gram,
                     sigma=sigma, floor_abs=0.0, coeff_max=2.0, feat_names=feat)
    syms = [sp.Symbol(f"x_{i}") for i in range(len(feat))]
    law_expr = sum((sp.Float(law[n]) * syms[i] for i, n in enumerate(feat)),
                   sp.S.Zero)
    y = s.A[:, j].copy()
    for i, n in enumerate(feat):
        y = y - law[n] * s.A[:, cols[i]]
    eps = band(m, law_expr)
    B, berr = ic_columns(basis_fns, x, s.patches, p=p)
    # the IC columns carry no DATA error; what they do carry is their own
    # quadrature error, which enters weighted by the (unknown) amplitudes and is
    # therefore bounded by a DECLARED amp_max -- audited afterwards against the
    # amplitudes that actually certified, never assumed silently
    eps = eps + amp_max * float(np.sum(berr))
    info.update({"ic_column_err": float(np.max(berr)),
                 "median_band": float(np.median(eps)),
                 "median_target": float(np.median(np.abs(y)))})
    return B, y, eps, info


def certify_state(B, y, eps, labels, *, window=(), amp_max: float = 10.0,
                  info: dict | None = None) -> StateCertificate:
    """Fit the amplitudes, certify the reconstruction, and report per-mode
    intervals -- with UNDETERMINED as a first-class outcome."""
    info = dict(info or {})
    n, dof = len(y), B.shape[1]
    h = n - dof
    cert = StateCertificate(certified=False, n_rows=int(n), dof=int(dof),
                            heldout=int(h), window=tuple(window),
                            basis=list(labels), amp_max_declared=amp_max)
    if h < MIN_HELDOUT:
        cert.abstain = "resolution"
        cert.notes.append(
            f"h = n - dof = {h} < {MIN_HELDOUT}: at most as many modes can be "
            "certified as there are independent patch equations, with margin; "
            "reduce the basis or add patches")
        return cert
    a, slack = _feasible_center(B, y, eps)
    if a is None:
        a, *_ = np.linalg.lstsq(B, y, rcond=None)
    cert.slack = None if slack is None else float(slack)
    syms = [sp.Symbol(f"x_{i}") for i in range(dof)]
    expr = sum((sp.Float(v) * s for v, s in zip(a, syms)), sp.S.Zero)
    resid = np.abs(B @ a - y)
    cert.residual_ratio = float(np.max(resid / eps))
    cert.amp_max_certified = float(np.max(np.abs(a)))
    if cert.amp_max_certified > amp_max:
        cert.abstain = "amplitude-bound-violated"
        cert.notes.append(
            f"declared amp_max={amp_max:g} but the fit needs "
            f"{cert.amp_max_certified:.3g}: re-assemble the band")
        return cert
    if not check(expr, syms, B, y, eps)["certified"]:
        cert.abstain = "no-state-explains-the-observations"
        cert.notes.append(
            "no initial condition in the declared basis reproduces the "
            f"observations within the band (worst residual "
            f"{cert.residual_ratio:.3g} x eps): the information is not there, "
            "or the declared law/basis is wrong for this data")
        return cert
    for jj, lab in enumerate(labels):
        iv = _mode_bounds(B, y, eps, jj)            # the JOINT projection
        others = y - (B @ a - B[:, jj] * a[jj])     # ...and the conditional one
        cond = parameter_interval(sp.Float(a[jj]) * syms[0], [syms[0]],
                                  B[:, [jj]], others, eps, sp.Float(a[jj]),
                                  max_rel=MAX_REL)
        cert.modes[lab] = {"value": float(a[jj]),
                           "interval": None if iv is None else
                           [float(iv[0]), float(iv[1])],
                           "conditional_interval": None if cond is None else
                           [float(cond[0]), float(cond[1])],
                           "determined": iv is not None,
                           # DETERMINED means the projection is bounded;
                           # RESOLVED means it excludes zero, i.e. the
                           # certificate can tell the mode is there at all. A
                           # bounded-but-zero-straddling interval is a real
                           # result and must not read as a recovered mode.
                           "resolved": bool(iv is not None
                                            and not (iv[0] <= 0.0 <= iv[1])),
                           "half_width": None if iv is None else
                           float(0.5 * (iv[1] - iv[0]))}
        if iv is None:
            cert.undetermined.append(lab)
    cert.partial = determination(
        [(lab, None if m["interval"] is None else m["interval"][0],
          None if m["interval"] is None else m["interval"][1])
         for lab, m in cert.modes.items()],
        status="state",
        note="joint projections over the feasible set of initial conditions: "
             "every state whose amplitudes lie in these intervals reproduces "
             "the observations within the declared band")
    cert.certified = True
    # |H| = 1: fixed basis, known law, no search over forms. This is a pure
    # chance-agreement bound and is NOT comparable with a law certificate's
    # alpha, which is corrected for the candidate space it searched.
    cert.alpha_log10 = significance_log10(expr, y, eps, 1)
    cert.notes.append(
        f"state certificate over window {tuple(window)}: every initial "
        f"condition whose amplitudes lie in the reported intervals reproduces "
        f"the observations within the declared band; "
        f"{len(cert.undetermined)} of {dof} modes are UNDETERMINED "
        f"({', '.join(cert.undetermined) or 'none'}) and are not claimed")
    if info:
        cert.notes.append(f"assembly: {info}")
    return cert

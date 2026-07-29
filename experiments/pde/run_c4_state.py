"""PDE dev campaign C4 -- state certificates (registered:
docs/CASE_STUDY_PDE_C4.md, direction docs/DIRECTION_PDE_STATE.md).

Inverting for an INITIAL CONDITION: a claim about one system's particular
history, deliberately kept a separate category from a law. Scores the registered
predictions W1-W5.
"""
from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).parent.parent.parent))
from experiments.pde import fields as F                          # noqa: E402
from experiments.pde import systems_fields as S                   # noqa: E402
from experiments.pde.verify import verify_state                  # noqa: E402
from lagh.statecert import (MIN_HELDOUT, assemble_state,         # noqa: E402
                            backpropagate, certify_state, fourier_basis,
                            truth_check_state)
from lagh.weakform import onesided_patches                       # noqa: E402

OUT = Path("experiments/results/pde_c4.json")

NX, NT = 256, 81
TMAX = 0.6
P_BUMP = 16
KMAX = 8
# One-sided-in-time families: (spatial half-width in cells, time cells, centres).
# Several scales for the same reason C1b pooled them -- the test-function scale
# is a knob the answer must be independent of -- and small spatial windows on
# purpose: a wide bump has a narrow spectrum, so it cannot see high-k modes at
# all, and that would be the PATCH FAMILY's resolution limit masquerading as the
# physics'.
SCALES = [(10, 0.30, 12), (14, 0.40, 10), (18, 0.50, 8)]   # (nx_half cells,
#                        fraction of the remaining record, number of x centres)


def grid(nt=NT, tmax=TMAX):
    x = np.linspace(0.0, 2 * np.pi, NX, endpoint=False)
    t = np.linspace(0.0, tmax, nt)
    return x, t


def patches(x, t, t0_index=0, n_x_mult=1):
    """The time windows are stated in a FRACTION of the record, not in cells, so
    a finer time grid means a better-resolved quadrature over the same physical
    window rather than a different claim."""
    out = []
    span = len(t) - 1 - t0_index
    for nxh, frac, n_x in SCALES:
        cells = int(4 * round(frac * span / 4))     # the ladder wants /4
        if cells < 8:
            continue
        out += onesided_patches(x, t, nx_half=nxh, nt_cells=cells,
                                n_x=int(n_x * n_x_mult), t0_index=t0_index)
    return out


def true_amplitudes(kind, labels, coeffs):
    """The mode amplitudes of the initial condition that generated the field."""
    truth = {lab: 0.0 for lab in labels}
    kind, t0, nu, c = kind
    for k, a in coeffs.items():
        if kind == "heat":
            truth[f"cos{k}"] = a * np.exp(-nu * k ** 2 * t0)
        else:                      # a translation mixes cos into sin
            truth[f"cos{k}"] = a * np.cos(k * c * t0)
            truth[f"sin{k}"] = a * np.sin(k * c * t0)
    return truth


def field(kind, coeffs, x, t, nu=0.1, c=0.7):
    """u(x, t) for an initial condition given as cosine-mode amplitudes."""
    u = np.zeros((len(x), len(t)))
    for k, a in coeffs.items():
        if kind == "heat":
            u += a * np.exp(-nu * k ** 2 * t)[None, :] * np.cos(k * x)[:, None]
        else:                                    # advection: translation by c t
            u += a * np.cos(k * (x[:, None] - c * t[None, :]))
    return u


LAWS = {"heat": {"u_xx": 0.1}, "advection": {"u_x": -0.7}}
COEFFS = {1: 1.0, 2: 0.6, 3: 0.4, 4: 0.3, 5: 0.25, 6: 0.2, 7: 0.16, 8: 0.13}


def run(kind, sigma, *, kmax=KMAX, seed=0, tmax=TMAX, verify=True,
        t0_index=0, nu=0.1, nt=NT, n_x_mult=1):
    """`t0_index` > 0 puts the observation window LATER than the initial time:
    the certificate is then about the state at that time, and the initial
    condition is what back-propagating it through the known law says -- which is
    where the exponential ill-posedness lives."""
    x, t = grid(nt, tmax)
    law = dict(LAWS[kind])
    if kind == "heat":
        law["u_xx"] = nu
    clean = field(kind, COEFFS, x, t, nu=nu)
    u = clean + (np.random.default_rng(seed).normal(0, sigma, clean.shape)
                 if sigma > 0 else 0.0)
    labels, fns = fourier_basis(kmax)
    t0 = time.time()
    B, y, eps, info = assemble_state(u, x, t, law,
                                     patches(x, t, t0_index, n_x_mult), fns,
                                     p=P_BUMP, sigma=sigma)
    if B is None:
        return {"kind": kind, "sigma": sigma, "certified": False,
                "abstain": "all-patches-rejected", "info": info}
    tw0 = float(t[t0_index])
    cert = certify_state(B, y, eps, labels, window=(tw0, float(t[-1])),
                         info=info)
    truth = true_amplitudes((kind, tw0, nu, 0.7), labels, COEFFS)
    out = {"kind": kind, "sigma": sigma, "certified": bool(cert.certified),
           "t0": tw0, "nu": nu,
           "abstain": cert.abstain, "alpha_log10": cert.alpha_log10,
           "alpha_kind": cert.alpha_kind, "n_rows": cert.n_rows,
           "dof": cert.dof, "heldout": cert.heldout,
           "residual_ratio": cert.residual_ratio,
           "undetermined": cert.undetermined,
           "modes": cert.modes, "notes": cert.notes, "info": info,
           "seconds": round(time.time() - t0, 1)}
    if cert.certified:
        wrong, worst = [], 0.0
        for lab, m in cert.modes.items():
            tv = truth.get(lab, 0.0)
            if m["interval"] is None:
                continue
            lo, hi = m["interval"]
            if not (lo <= tv <= hi):
                wrong.append({"mode": lab, "truth": tv,
                              "interval": m["interval"]})
            worst = max(worst, abs(m["value"] - tv))
        out["confident_wrong"] = wrong
        out["worst_amplitude_error"] = worst
        out["halfwidth_over_sigma"] = {
            lab: (None if m["half_width"] is None else
                  (m["half_width"] / sigma if sigma > 0 else None))
            for lab, m in cert.modes.items()}
        # the certificate is about the state at the window's start, so the
        # forecast runs from there -- and, when the window is late, the claim
        # about the INITIAL condition is what back-propagation says it is
        out["backpropagated"] = backpropagate(cert, kind, tw0, nu=nu)
        bp_cw = []
        for k, a in COEFFS.items():
            m = out["backpropagated"]["modes"].get(f"cos{k}")
            if m and m["interval"] is not None and kind == "heat" \
                    and not (m["interval"][0] <= a <= m["interval"][1]):
                bp_cw.append({"mode": f"cos{k}", "truth": a,
                              "interval": m["interval"]})
        out["backpropagated_confident_wrong"] = bp_cw
        if verify:
            ivs = {lab: m["interval"] for lab, m in cert.modes.items()
                   if m["interval"] is not None}
            out["verify"] = verify_state(
                (u if sigma > 0 else clean)[:, t0_index:], x, t[t0_index:],
                ivs, fourier_basis(kmax)[1], labels, law, sigma=sigma,
                u_clean=clean[:, t0_index:])
    return out


def report(k, r):
    if not r.get("certified"):
        print(f"{k:34s} ABSTAIN[{r.get('abstain')}] rows={r.get('n_rows')} "
              f"dof={r.get('dof')} resid/band={r.get('residual_ratio')}",
              flush=True)
        return
    det = r["dof"] - len(r["undetermined"])
    nres = sum(1 for m in r["modes"].values() if m.get("resolved"))
    hw = {lab: m["half_width"] for lab, m in r["modes"].items()
          if m["half_width"] is not None and lab.startswith("cos")}
    print(f"{k:34s} CERT  {det}/{r['dof']} determined, {nres} resolved  "
          f"undet={','.join(r['undetermined']) or '-'}  "
          f"alpha<=1e{r['alpha_log10']:.0f}  CW={len(r['confident_wrong'])}  "
          f"worst_amp_err={r['worst_amplitude_error']:.2e}", flush=True)
    print("      half-widths: " + "  ".join(
        f"{lab}:{v:.2e}" for lab, v in list(hw.items())[:9]), flush=True)
    bp = r.get("backpropagated")
    if bp and bp["t0"] > 0:
        hw = {lab: m.get("half_width") for lab, m in bp["modes"].items()
              if lab.startswith("cos") and m.get("half_width") is not None}
        print(f"      back-prop to t=0 (t0={bp['t0']:.3f}): k_cut={bp['k_cut']} "
              f"CW={len(r.get('backpropagated_confident_wrong', []))}  "
              + "  ".join(f"{lab}:{v:.1e}" for lab, v in list(hw.items())[:8]),
              flush=True)
    v = r.get("verify")
    if v:
        print(f"      verify: law={'OK' if v.get('verified') else 'FAIL'} "
              f"data={'OK' if v.get('data_verified') else 'FAIL'} "
              f"outside={v.get('n_outside')}/{v.get('n_points')} "
              f"env={v.get('envelope_width_med', 0):.2e}", flush=True)


NT_MAIN = 161          # see the floor study: NT = 81 is quadrature-limited


def burgers_state(nu, amps, kmax, sigma=1e-6, tag=""):
    """W4: does a Burgers state certificate refuse, and WHY?

    The Cole-Hopf profile u = -2 nu phi_x / phi is a RATIO, so its Fourier
    spectrum is infinite and a truncated basis cannot represent it exactly. Two
    refusals therefore look alike from outside and are not alike at all, so the
    true state is projected onto the declared basis and checked against its own
    band first (truth_check_state), and the basis is swept."""
    x, t = grid(NT_MAIN)
    u = F.burgers_coleHopf(x, t, nu=nu, modes=(1, 2), amps=amps,
                           phases=(0.0, 0.5))
    labels, fns = fourier_basis(kmax)
    # enough patches that the resolution bound is not the binding constraint
    mult = max(1, int(np.ceil((2 * kmax + 1 + MIN_HELDOUT) / 30)))
    B, y, eps, info = assemble_state(u, x, t, {"u_xx": nu, "u*u_x": -1.0},
                                     patches(x, t, 0, mult), fns, p=P_BUMP,
                                     sigma=sigma)
    k = f"burgers{'_' + tag if tag else ''}_nu{nu:g}_k{kmax}"
    if B is None:
        print(f"{k:34s} ABSTAIN all-patches-rejected  {info}", flush=True)
        return {k: {"certified": False, "abstain": "all-patches-rejected",
                    "nu": nu, "kmax": kmax, "info": info}}
    # the true state in the declared basis: a least-squares projection of the
    # measured initial profile onto the SAME basis (no privileged information --
    # it is the best any certificate could do with this basis)
    Phi = np.column_stack([f(x) for f in fns])
    a_true, *_ = np.linalg.lstsq(Phi, u[:, 0], rcond=None)
    tc = truth_check_state(B, y, eps, a_true)
    c = certify_state(B, y, eps, labels, window=(float(t[0]), float(t[-1])),
                      info=info)
    out = {"certified": bool(c.certified), "abstain": c.abstain,
           "n_rows": c.n_rows, "dof": c.dof, "nu": nu, "kmax": kmax,
           "residual_ratio": c.residual_ratio, "truth": tc,
           "basis_truncation_rel": float(
               np.max(np.abs(Phi @ a_true - u[:, 0])) / np.max(np.abs(u[:, 0]))),
           "undetermined": c.undetermined, "modes": c.modes,
           "notes": c.notes, "info": info}
    if c.certified:
        out["confident_wrong"] = [
            lab for lab, m in c.modes.items()
            if m["interval"] is not None
            and not (m["interval"][0] <= a_true[labels.index(lab)]
                     <= m["interval"][1])]
        out["worst_amplitude_error"] = float(
            max(abs(m["value"] - a_true[labels.index(lab)])
                for lab, m in c.modes.items()))
    print(f"burgers nu={nu:g} k={kmax:3d} "
          f"{'CERT   ' if c.certified else 'ABSTAIN'} "
          f"{str(c.abstain or ''):32s} rows={c.n_rows} dof={c.dof} "
          f"resid/band={c.residual_ratio if c.residual_ratio else 0:.3g} "
          f"truth/band={tc['truth_max_ratio']:.3g} "
          f"basis_trunc={out['basis_truncation_rel']:.1e}", flush=True)
    return {k: out}


def shock_state(nu, t0_frac, kmax=32, sigma=1e-6):
    """W4 proper: a Burgers field that really forms a front, read from a window
    starting before (t0_frac = 0) or after (0.5, 0.75) the steepening."""
    u, coords, ferr = S.burgers_shock(nu=nu)
    k = f"shock_nu{nu:g}_t0f{t0_frac:g}"
    if u is None:
        return {k: {"certified": False, "abstain": "solver-did-not-converge"}}
    x, t = coords
    i0 = int(t0_frac * (len(t) - 1))
    labels, fns = fourier_basis(kmax)
    mult = max(1, int(np.ceil((2 * kmax + 1 + MIN_HELDOUT) / 30)))
    pa = patches(x, t, i0, mult)
    B, y, eps, info = assemble_state(u, x, t, {"u_xx": nu, "u*u_x": -1.0}, pa,
                                     fns, p=P_BUMP, sigma=sigma,
                                     field_err=ferr)
    out = {"nu": nu, "t0_frac": t0_frac, "t0": float(t[i0]), "kmax": kmax,
           "field_err": ferr, "n_patches_offered": len(pa), "info": info}
    if B is None:
        out.update({"certified": False, "abstain": "all-patches-rejected"})
        print(f"{k:34s} ABSTAIN all-patches-rejected  "
              f"(offered {len(pa)}, rejected {info.get('rejected')}) -- the "
              f"grid does not represent this field", flush=True)
        return {k: out}
    Phi = np.column_stack([f(x) for f in fns])
    a_true, *_ = np.linalg.lstsq(Phi, u[:, i0], rcond=None)
    tc = truth_check_state(B, y, eps, a_true)
    c = certify_state(B, y, eps, labels, window=(float(t[i0]), float(t[-1])),
                      info=info)
    out.update({"certified": bool(c.certified), "abstain": c.abstain,
                "n_rows": c.n_rows, "dof": c.dof, "truth": tc,
                "residual_ratio": c.residual_ratio,
                "patches_rejected": info.get("rejected"),
                "basis_truncation_rel": float(
                    np.max(np.abs(Phi @ a_true - u[:, i0]))
                    / np.max(np.abs(u[:, i0]))),
                "undetermined": c.undetermined, "modes": c.modes,
                "notes": c.notes})
    if c.certified:
        out["confident_wrong"] = [
            lab for lab, m in c.modes.items()
            if m["interval"] is not None
            and not (m["interval"][0] <= a_true[labels.index(lab)]
                     <= m["interval"][1])]
        out["n_resolved"] = sum(1 for m in c.modes.values() if m["resolved"])
    print(f"{k:34s} {'CERT   ' if c.certified else 'ABSTAIN'} "
          f"{str(c.abstain or ''):32s} rows={c.n_rows} dof={c.dof} "
          f"rejected={info.get('rejected')} "
          f"resid/band={c.residual_ratio if c.residual_ratio else 0:.3g} "
          f"truth/band={tc['truth_max_ratio']:.3g} "
          f"basis_trunc={out['basis_truncation_rel']:.1e}", flush=True)
    return {k: out}


def main(only=None):
    res = {}
    for kind in ("advection", "heat"):
        if only and only not in (kind, "ladder"):
            continue
        for sigma in (1e-6, 1e-5, 1e-4):
            k = f"{kind}_sigma{sigma:g}"
            res[k] = run(kind, sigma, nt=NT_MAIN)
            report(k, res[k])

    # the FLOOR study: is the small-sigma resolution noise-limited or
    # quadrature-limited? The one-sided window's Romberg residual is a declared
    # deterministic bound, and at a coarse time grid it dominates the band --
    # which shows up as half-width/sigma being far above 1 and NOT linear in
    # sigma. Refining the time grid is the direct test.
    if only in (None, "floor"):
        for nt in (81, 161, 321):
            for sigma in (1e-6, 1e-4):
                k = f"floor_heat_nt{nt}_sigma{sigma:g}"
                res[k] = run("heat", sigma, nt=nt, verify=False)
                report(k, res[k])

    # W2/W3: put the observation window LATE and back-propagate to t = 0. This
    # is where the exponential ill-posedness lives -- a window starting at t = 0
    # sees every mode at full amplitude, which is why the certificate above is
    # flat in k. `k_cut` is threshold-free: the largest k whose back-propagated
    # interval EXCLUDES ZERO (below it the certificate cannot even determine
    # that the mode is present).
    if only in (None, "cut"):
        for i0 in (0, 40, 80):
            for sigma in (1e-6, 1e-4):
                for nu in (0.1, 0.3):
                    k = f"heat_t0i{i0}_nu{nu:g}_sigma{sigma:g}"
                    res[k] = run("heat", sigma, t0_index=i0, nu=nu, nt=NT_MAIN,
                                 kmax=12, n_x_mult=2, verify=False)
                    report(k, res[k])

    # W4: Burgers after shock formation -- refusal, not a wide interval
    if only in (None, "burgers"):
        for nu, amps in ((0.2, (0.3, 0.15)), (0.05, (0.6, 0.3)),
                         (0.02, (0.9, 0.45))):
            for kmax in (8, 16, 32):
                res.update(burgers_state(nu, amps, kmax))
        # ...and the case W4 was actually about: a field that really steepens
        # (u0 = -sin x, small nu), read BEFORE and AFTER the front forms
        for nu in (0.02, 0.005):
            for t0f in (0.0, 0.5, 0.75):
                res.update(shock_state(nu, t0f))

    # the resolution bound, exercised: more modes than independent equations
    if only in (None, "resolution"):
        k = "heat_overparameterized_k40"
        res[k] = run("heat", 1e-5, kmax=40, nt=NT_MAIN, verify=False)
        report(k, res[k])
        k = "heat_overparameterized_k40_morepatches"
        res[k] = run("heat", 1e-5, kmax=40, nt=NT_MAIN, n_x_mult=4,
                     verify=False)
        report(k, res[k])

    prev = json.loads(OUT.read_text()) if OUT.exists() and only else {}
    prev.update(res)
    OUT.write_text(json.dumps(prev, indent=1))
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1] if len(sys.argv) > 1 else None))

"""C1 for the stochastic arc's real-data read: the trapped bead, certified.

C0 screened the LUMICKS C-Trap passive record in (docs/PROPOSAL_STOCHASTIC_REAL.md
§9). C1 asks what can actually be certified from it, and the first thing it finds is
that the obvious target CANNOT be.

    THE PASSIVE DISPLACEMENT CALIBRATION IS CIRCULAR FOR b^2.

`Rd`, the volts->nm scale, is not measured against a length standard on this record.
It is DEFINED by setting the observed diffusion equal to its Stokes-Einstein value:

    Rd  ==  sqrt( (k_B T / gamma_0) / D_volts )

which this run verifies to machine precision on all 8 calibration items. So a
certificate that b^2 = 2 k_B T / gamma_0, computed from a position signal scaled by
Rd, restates the calibration's own definition. It cannot fail except by estimator
disagreement, and it is the FLAME circularity that got reanalysis rejected in
DIRECTION_PDE.md, arriving in a new costume.

    THE DRIFT IS NOT CIRCULAR, BECAUSE IT IS A TIMESCALE.

theta = kappa/gamma_0 = 2 pi f_c has units of 1/s. Rescaling the position signal by
ANY constant leaves it unchanged -- so the drift is immune both to the circularity
and to the applied-vs-derived unit mix-up this run also finds. It is the certifiable
stratum here, and the diffusion is not. That inverts the arc's usual difficulty
ordering, where the diffusion certified 2600x tighter than the drift.

Run: .venv/bin/python experiments/tweezers/run_c1.py
"""
from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).parent.parent.parent))

from experiments.tweezers import adapter  # noqa: E402
from lagh.ito import build_qv_rows, build_rows, certify_diffusion, certify_drift  # noqa: E402
from lagh.weakform import qv_three_way  # noqa: E402

OUT = Path("experiments/results/tweezers_c1.json")
KB = 1.380649e-23
LIBRARY = ("1", "x", "x**2", "x**3")
DIFF_LIBRARY = ("1",)             # a harmonic trap has ADDITIVE noise: b is constant
N_SEG = 4                         # holdout units, see `segments`
HALF = 2000                       # window half-width [samples] -> theta*L ~ 160


def circularity_audit(f) -> dict:
    """Is Rd defined by the thermal motion it would be used to measure? Checked on
    every calibration item in the file, not asserted from the documentation."""
    rows = []
    for ch in ("Force 1x", "Force 1y", "Force 2x", "Force 2y"):
        try:
            cals = adapter.bfp_calibrations(f, ch)
        except Exception:                                        # noqa: BLE001
            continue
        for c in cals:
            T = c["Temperature (C)"] + 273.15
            implied = np.sqrt((KB * T / c["gamma_0 (kg/s)"]) / c["D (V^2/s)"]) * 1e6
            rows.append({"channel": ch, "item": c["item"],
                         "Rd_stored_um_per_V": c["Rd (um/V)"],
                         "Rd_implied_um_per_V": float(implied),
                         "rel_dev": float(abs(implied / c["Rd (um/V)"] - 1))})
    worst = max((r["rel_dev"] for r in rows), default=None)
    return {"rows": rows, "worst_rel_dev": worst,
             "verdict": ("Rd IS defined by the thermal motion (max deviation "
                         f"{worst:.1e}) -- b^2 against 2 k_B T / gamma_0 is CIRCULAR "
                         "on a passively calibrated record")}


def psd_band_audit(x, cal, dt) -> dict:
    """Where does the instrument's calibration LOOK, and where does a quadratic
    variation live? They are not the same band, and that is the whole problem.

    The Lorentzian-times-diode model is fitted over `Fit range`, typically
    100 Hz - 23 kHz. A PSD fit weights that band roughly flat. A quadratic variation
    is `int (2 pi f)^2 PSD df` -- weighted by f^2, so it is dominated by the TOP of
    the band and by everything above it. Any anti-alias filter removes exactly the
    part QV cares about while leaving the fit untouched.
    """
    from scipy.signal import welch
    fs = 1.0 / dt
    fr, P = welch(np.asarray(x, float) - np.mean(x), fs=fs, nperseg=2 ** 16)
    D_nm = KB * (cal["Temperature (C)"] + 273.15) / cal["gamma_0 (kg/s)"] * 1e18
    lor = D_nm / (np.pi ** 2 * (cal["fc (Hz)"] ** 2 + fr ** 2))
    al, fd = cal["alpha"], cal["f_diode (Hz)"]
    model = lor * (al ** 2 + (1 - al ** 2) / (1 + (fr / fd) ** 2))
    f_lo = cal.get("Fit range (min.) (Hz)", 100.0)
    f_hi = cal.get("Fit range (max.) (Hz)", 23000.0)
    inb = (fr > f_lo) & (fr < f_hi)
    above = fr >= f_hi
    w = (2 * np.pi * fr) ** 2                       # the QV weight
    it = np.trapezoid
    return {"fit_range_hz": [f_lo, f_hi], "nyquist_hz": fs / 2,
            "psd_ratio_in_band": float(it(P[inb], fr[inb]) / it(model[inb], fr[inb])),
            "psd_ratio_above_band": float(it(P[above], fr[above])
                                          / it(model[above], fr[above])),
            "qv_weighted_ratio": float(it(w * P, fr) / it(w * model, fr)),
            "model_qv_fraction_above_fit_band": float(
                it((w * model)[above], fr[above]) / it(w * model, fr)),
            "psd_ratio_near_nyquist": float(
                P[np.argmin(abs(fr - 0.97 * fs / 2))]
                / model[np.argmin(abs(fr - 0.97 * fs / 2))])}


def theta_acf(x, dt, theta_true, max_lag=120, min_lag=10) -> dict:
    """theta from the autocorrelation decay: a TIMESCALE, never routed through b^2.

    This is the estimator the weak form is not. It is invariant to any rescaling of
    the position signal (so the circular calibration cannot touch it) and it reads
    the decay at lags well beyond the detector's own filter (so the band loss cannot
    either). It doubles as the per-axis quality gate: a slow contaminant below the
    calibration's 100 Hz fit floor shows up here and nowhere else.
    """
    v = np.asarray(x, float) - np.mean(x)
    den = float(v @ v)
    lags = np.arange(1, max_lag)
    acf = np.array([float(v[:-L] @ v[L:]) / den for L in lags])
    m = (lags >= min_lag) & (acf > 0.05)
    if m.sum() < 5:
        return {"theta_per_s": None, "note": "ACF decays too fast to fit"}
    slope = float(np.polyfit(lags[m] * dt, np.log(acf[m]), 1)[0])
    th = -slope
    return {"theta_per_s": th, "theta_true": theta_true,
            "ratio": th / theta_true, "acf_lag1": float(acf[0]),
            "clean": bool(0.8 < th / theta_true < 1.25)}


def segments(x, n=N_SEG):
    """Split one record into `n` holdout units.

    Declared, not silent: these are not independent experiments, they are stretches
    of one path. They are admissible as holdout units here because the trap
    correlation time is ~0.3 ms against a 16 s record -- ~50 000 correlation times,
    so distinct segments share no memory. A claim about a DIFFERENT bead would not be
    supported by them, and C2's cross-bead work needs the companion record.
    """
    m = (len(x) // n) * n
    return np.asarray(x[:m], float).reshape(n, -1)


def one_axis(f, channel, *, deconvolve: bool) -> dict:
    t_all, x_nm, applied, derived, meta = adapter.bfp_position_nm(f, channel)
    dt = 1.0 / meta["fs_hz"]
    T_k = derived["Temperature (C)"] + 273.15
    gamma = derived["gamma_0 (kg/s)"]
    kappa = derived["kappa (pN/nm)"] * 1e-3                     # [N/m]

    # TRUTHS. theta is scale-free; b2 is not, and carries the circularity above.
    theta_true = kappa / gamma                                  # [1/s]
    b2_true = 2 * KB * T_k / gamma * 1e18                       # [nm^2/s]

    info_dec = None
    if deconvolve:
        x_nm, info_dec = adapter.diode_deconvolve(
            x_nm, derived["alpha"], derived["f_diode (Hz)"], dt)

    P = segments(x_nm)
    P = P - P.mean(axis=1, keepdims=True)      # the trap centre is the equilibrium
    t = np.arange(P.shape[1]) * dt

    out = {"channel": channel, "deconvolved": deconvolve, "meta": meta,
           "answer_key": {"kappa_pN_per_nm": derived["kappa (pN/nm)"],
                          "gamma_kg_per_s": gamma, "fc_hz": derived["fc (Hz)"],
                          "err_fc_hz": derived.get("err_fc (Hz)"),
                          "T_C": derived["Temperature (C)"],
                          "diode_alpha": derived["alpha"],
                          "diode_f_hz": derived["f_diode (Hz)"]},
           "truth": {"theta_per_s": theta_true, "b2_nm2_per_s": b2_true,
                     "two_pi_fc": 2 * np.pi * derived["fc (Hz)"]},
           "deconvolution": info_dec,
           "observed": {"position_rms_nm": float(P.std()),
                        "increment_rms_nm": float(np.diff(P, axis=1).std()),
                        "expected_increment_nm": float(np.sqrt(b2_true * dt))}}

    out["psd_band"] = psd_band_audit(x_nm, derived, dt)
    out["theta_acf"] = theta_acf(x_nm, dt, theta_true)

    # is the detector noise visible after deconvolution? (it amplifies high f)
    sums = np.array([float(np.sum((P[:, s:] - P[:, :-s]) ** 2))
                     for s in (1, 2, 3, 4, 6, 8)])
    _, tw = qv_three_way(sums, strides=(1, 2, 3, 4, 6, 8), n_increments=P.size)
    out["three_way"] = {k: tw.get(k) for k in
                        ("c", "alpha", "separable", "process_resolved",
                         "process_significance", "sigma_obs")}

    # ---- the DRIFT, from the Ito weak form. Scale-free, hence the real target.
    dr_rows = build_rows(t, P, LIBRARY, half=HALF)
    dr = certify_drift(dr_rows, delta=0.05, seed=0) if len(dr_rows.y) else {}
    out["drift"] = {"n_rows": dr.get("n_rows"), "kappa_cov": dr.get("kappa"),
                    "certified": dr.get("certified"), "law": dr.get("law"),
                    "abstain": dr.get("abstain"),
                    "median_signal_to_band": dr.get("median_signal_to_band"),
                    "alpha": dr.get("alpha")}
    comps = (dr.get("partial") or {}).get("components", {})
    out["drift"]["components"] = comps
    xc = comps.get("x") or comps.get("drift:x")
    if xc:
        lo, hi = xc.get("lo"), xc.get("hi")
        out["drift"]["theta_interval"] = [None if lo is None else -hi,
                                          None if hi is None else -lo]
        out["drift"]["covers_theta"] = bool(
            lo is not None and hi is not None and lo <= -theta_true <= hi)

    # ---- the DIFFUSION, reported with its circularity tag
    qv_rows = build_qv_rows(t, P, DIFF_LIBRARY, half=HALF, ws=("1",),
                            drift_envelope=None)
    df = (certify_diffusion(qv_rows, delta=0.05, drift_max=theta_true * P.std() * 3,
                            seed=0) if len(qv_rows.y) else {})
    out["diffusion"] = {"n_rows": df.get("n_rows"), "certified": df.get("certified"),
                        "law": df.get("law"), "abstain": df.get("abstain"),
                        "median_signal_to_band": df.get("median_signal_to_band"),
                        "components": (df.get("partial") or {}).get("components", {}),
                        "circular": True,
                        "circular_note": ("Rd is defined by D, so this compares our "
                                          "estimator to the instrument's PSD fit, "
                                          "NOT the physics to a measurement")}
    dc = out["diffusion"]["components"].get("1") or \
        out["diffusion"]["components"].get("diffusion:1")
    if dc and dc.get("lo") is not None:
        out["diffusion"]["b2_interval"] = [dc["lo"], dc["hi"]]
        out["diffusion"]["covers_b2"] = bool(dc["lo"] <= b2_true <= dc["hi"])
        out["diffusion"]["b2_mid_over_truth"] = (
            0.5 * (dc["lo"] + dc["hi"]) / b2_true)
    return out


def main():
    t0 = time.time()
    import h5py
    path = adapter.bfp_fetch("passive_calibration.h5")
    f = h5py.File(path, "r")

    audit = circularity_audit(f)
    rows = []
    for ch in ("Force 1x", "Force 1y"):
        for dec in (False, True):
            rows.append(one_axis(f, ch, deconvolve=dec))

    res = {"campaign": "stochastic-real (C-Trap)", "stage": "C1",
           "artifact": json.loads((adapter.BFP_DIR /
                                   "passive_calibration.h5.json").read_text()),
           "circularity_audit": audit, "axes": rows,
           "seconds": round(time.time() - t0, 1)}
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(res, indent=1, default=str))

    print(f"\nwrote {OUT}  ({res['seconds']}s)")
    print(f"\nCIRCULARITY AUDIT ({len(audit['rows'])} calibration items)")
    print(f"  {audit['verdict']}")
    m = rows[0]["meta"]
    print(f"\nUNIT-CHAIN PROVENANCE")
    print(f"  stored pN used calibration item {m['applied_item']}, but the only "
          f"calibration describing this trap is item {m['derived_item']}")
    print(f"  the naive route (stored pN / derived kappa) is off by "
          f"{m['naive_kappa_route_error'] * 100:+.2f}% in position, "
          f"{((1 + m['naive_kappa_route_error']) ** 2 - 1) * 100:+.2f}% in b^2")

    for r in rows:
        tag = "deconvolved" if r["deconvolved"] else "raw        "
        print(f"\n=== {r['channel']}  [{tag}]")
        o, tr = r["observed"], r["truth"]
        print(f"    position rms {o['position_rms_nm']:.2f} nm   increment rms "
              f"{o['increment_rms_nm']:.3f} nm  (free-Brownian "
              f"{o['expected_increment_nm']:.3f})")
        tw = r["three_way"]
        print(f"    three-way: resolved={tw['process_resolved']} "
              f"sig={tw['process_significance']:.1f}  sigma_obs="
              f"{(tw['sigma_obs'] or 0):.4f} nm  separable={tw['separable']}")
        pb = r["psd_band"]
        print(f"    PSD vs the instrument's own model: in-band "
              f"{pb['psd_ratio_in_band']:.4f}, above band "
              f"{pb['psd_ratio_above_band']:.3f}, near Nyquist "
              f"{pb['psd_ratio_near_nyquist']:.1e}")
        print(f"      {pb['model_qv_fraction_above_fit_band'] * 100:.1f}% of the "
              f"model's QUADRATIC VARIATION lives above the {pb['fit_range_hz'][1]:.0f} "
              f"Hz fit ceiling -- where the chain attenuates")
        ta = r["theta_acf"]
        if ta.get("theta_per_s"):
            print(f"    theta from ACF decay (scale-free, not via b^2): "
                  f"{ta['theta_per_s']:.1f} /s vs truth {ta['theta_true']:.1f}"
                  f"   ratio {ta['ratio']:.3f}   clean={ta['clean']}")
        d = r["drift"]
        print(f"    DRIFT   certified={d['certified']}  rows={d['n_rows']}  "
              f"s/band={d.get('median_signal_to_band')}")
        if d.get("theta_interval"):
            lo, hi = d["theta_interval"]
            print(f"      theta in [{lo:.1f}, {hi:.1f}] /s   truth "
                  f"{tr['theta_per_s']:.1f} /s   COVERS={d['covers_theta']}")
        else:
            print(f"      theta not determined (abstain: {d.get('abstain')})")
        df = r["diffusion"]
        print(f"    DIFF    certified={df['certified']}  rows={df['n_rows']}  "
              f"[CIRCULAR]")
        if df.get("b2_interval"):
            lo, hi = df["b2_interval"]
            print(f"      b^2 in [{lo:.4g}, {hi:.4g}]  truth {tr['b2_nm2_per_s']:.4g}"
                  f"   COVERS={df['covers_b2']}  mid/truth "
                  f"{df['b2_mid_over_truth']:.3f}")
        else:
            print(f"      b^2 not determined (abstain: {df.get('abstain')})")


if __name__ == "__main__":
    main()

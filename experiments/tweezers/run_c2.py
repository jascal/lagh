"""C2 for the stochastic arc's real-data read: the instrument, declared.

C1 (docs/CASE_STUDY_TWEEZERS_C1.md) measured what a real trapped bead does to the
arc's estimators -- the diffusion read 32% of truth, the drift 46% low -- and
listed what C2 must do, in order: gate axes on the ACF timescale; declare the
band loss as a measured input; use the b^2-free drift form; and move to the
ACTIVE record, where the nanostage gives a displacement scale that does not pass
through the thermal motion. This run does all four, on both records.

    passive_calibration.h5             the C1 record, 16 s, both axes
    near_surface_active_calibration.h5 the stage driven at ~38 Hz along X

Nothing here uses `answer_key`-style external truth: the calibration's own
2 pi f_c is the drift reference (a timescale, scale-free), 2 k_B T / gamma_0 the
diffusion reference (bulk Stokes; near a surface the true drag is higher, and
that is what the active scale measures rather than assumes).

Run: .venv/bin/python experiments/tweezers/run_c2.py
"""
from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).parent.parent.parent))

from experiments.tweezers import adapter                                 # noqa: E402
from lagh.instrument import axis_gate, band_loss, drive_scale           # noqa: E402
from lagh.ito import build_lag_rows, build_qv_rows, certify_diffusion   # noqa: E402

OUT = Path("experiments/results/tweezers_c2.json")
KB = 1.380649e-23
N_SEG = 4
HALF = 2000
LAGS = (16, 32)
LIB = ("1", "x")


def segments(x, n=N_SEG):
    m = (len(x) // n) * n
    P = np.asarray(x[:m], float).reshape(n, -1)
    return P - P.mean(axis=1, keepdims=True)


def theta_lagform(P, dt, theta_ref):
    """The b^2-free drift, scored: lagged form at two lags, inverted analytically
    (run_lagform.py). An estimator, not a certificate -- C1's ridge stands."""
    t = np.arange(P.shape[1]) * dt
    out = {}
    for lag in LAGS:
        r = build_lag_rows(t, P, LIB, lags=(lag,), ws=("x",), half=HALF)
        c, *_ = np.linalg.lstsq(r.A, r.y, rcond=None)
        a = -float(c[list(r.names).index("x")])
        z = a * lag * dt
        inv = float(-np.log(1 - z) / (lag * dt)) if 0 < z < 1 else float("nan")
        out[f"lag{lag}"] = {"theta_apparent": a, "theta_inverted": inv,
                            "ratio": inv / theta_ref}
    return out


def one_axis(f, channel, *, record: str, stage=None) -> dict:
    t_all, x_nm, applied, de, meta = adapter.bfp_position_nm(f, channel)
    fs = meta["fs_hz"]
    dt = 1.0 / fs
    T_k = de["Temperature (C)"] + 273.15
    gamma0 = de["gamma_0 (kg/s)"]
    theta_ref = 2 * np.pi * de["fc (Hz)"]
    b2_bulk = 2 * KB * T_k / gamma0 * 1e18                     # nm^2/s, bulk Stokes
    out = {"record": record, "channel": channel, "fs_hz": fs, "n": int(len(x_nm)),
           "calibration": {"fc_hz": de["fc (Hz)"], "err_fc_hz": de.get("err_fc (Hz)"),
                           "alpha": de["alpha"], "f_diode_hz": de["f_diode (Hz)"],
                           "Rd_um_per_V": de["Rd (um/V)"], "gamma0_kg_per_s": gamma0,
                           "T_C": de["Temperature (C)"],
                           "fit_range_hz": [de.get("Fit range (min.) (Hz)", 100.0),
                                            de.get("Fit range (max.) (Hz)", 23000.0)]},
           "reference": {"theta_per_s": theta_ref, "b2_bulk_nm2_per_s": b2_bulk}}

    # ---- the ACTIVE scale: lock in on the stage, remove the drive, read Rd
    scale_note = "stored Rd (defined by the thermal motion: circular for b^2)"
    x_used = x_nm
    if stage is not None:
        volts = x_nm / (de["Rd (um/V)"] * 1e3)
        ds = drive_scale(stage, volts, fs, de["fc (Hz)"])
        Rd_active = ds["scale"]                                 # stage units / V
        out["active_scale"] = {k: ds[k] for k in ("f_drive_hz", "stage_amplitude",
                                                  "response_amplitude", "transfer",
                                                  "response_power_share", "n")}
        out["active_scale"]["Rd_active_stage_units_per_V"] = Rd_active
        out["active_scale"]["Rd_active_over_stored"] = Rd_active / de["Rd (um/V)"]
        driven = ds["response_power_share"] > 0.05
        out["active_scale"]["driven_axis"] = bool(driven)
        if driven:
            # the stage reads in um (Bluelake convention): Rd_active in um/V
            x_used = ds["residual"] * Rd_active * 1e3
            scale_note = ("ACTIVE Rd from the nanostage drive: no thermal quantity "
                          "enters, so b^2 against 2 k_B T / gamma is a physics claim "
                          "and its ratio to bulk Stokes is the near-surface drag")
        else:
            x_used = ds["residual"] * de["Rd (um/V)"] * 1e3
            scale_note += "; drive line removed (undriven axis)"
    out["scale_note"] = scale_note

    # ---- 1. the axis gate, before anything is claimed
    gate = axis_gate(x_used, dt, theta_ref)
    out["axis_gate"] = gate

    # ---- 2. the band loss, declared
    bl = band_loss(x_used, fs, fit_range=tuple(out["calibration"]["fit_range_hz"]),
                   alpha=de["alpha"], f_diode=de["f_diode (Hz)"])
    out["band_loss"] = bl
    qv_raw = float(np.mean(np.diff(x_used) ** 2) / dt)
    out["qv"] = {"raw_over_bulk": qv_raw / b2_bulk,
                 "declared_over_bulk": qv_raw / bl["r_total"] / b2_bulk}
    if not gate["passed"]:
        out["verdict"] = "axis refused by the ACF gate; nothing claimed"
        return out

    P = segments(x_used)
    t = np.arange(P.shape[1]) * dt
    # ---- 3. the diffusion, with the loss declared
    dmax = theta_ref * P.std() * 3
    res = {}
    for tag, loss in (("raw", None), ("declared", (bl["r_total"], bl["r_se"]))):
        rows = build_qv_rows(t, P, ("1",), half=HALF, ws=("1",), band_loss=loss)
        d = certify_diffusion(rows, delta=0.05, drift_max=dmax, seed=0)
        iv = (d.get("admissible") or {}).get("diffusion:1")
        res[tag] = {"n_rows": d.get("n_rows"), "certified": d.get("certified"),
                    "abstain": d.get("abstain"),
                    "median_signal_to_band": d.get("median_signal_to_band"),
                    "b2_interval": iv,
                    "covers_bulk": (None if iv is None else
                                    bool(iv[0] <= b2_bulk <= iv[1])),
                    "mid_over_bulk": (None if iv is None else
                                      0.5 * (iv[0] + iv[1]) / b2_bulk)}
    out["diffusion"] = res
    iv = res["declared"]["b2_interval"]
    if iv is not None and stage is not None and out["active_scale"]["driven_axis"]:
        g_eff = 2 * KB * T_k / (0.5 * (iv[0] + iv[1]) * 1e-18)
        out["drag"] = {"gamma_eff_kg_per_s": g_eff, "gamma_eff_over_bulk": g_eff / gamma0,
                       "interval_over_bulk": [2 * KB * T_k / (iv[1] * 1e-18) / gamma0,
                                              2 * KB * T_k / (iv[0] * 1e-18) / gamma0],
                       "note": "near-surface drag from an independent scale, "
                               "tag empirical: Faxen's height is not in the file"}
    # ---- 4. the drift, b^2-free, scored
    out["drift"] = {"theta_acf_over_ref": gate["ratio"],
                    "lagform": theta_lagform(P, dt, theta_ref)}
    out["verdict"] = "axis gated in; diffusion declared; drift scored"
    return out


def main():
    t0 = time.time()
    import h5py
    rows = []
    fp = h5py.File(adapter.BFP_DIR / "passive_calibration.h5", "r")
    for ch in ("Force 1x", "Force 1y"):
        rows.append(one_axis(fp, ch, record="passive"))
    fa = h5py.File(adapter.BFP_DIR / "near_surface_active_calibration.h5", "r")
    stage = fa["Nanostage position/X"][()]
    for ch in ("Force 1x", "Force 1y"):
        rows.append(one_axis(fa, ch, record="active", stage=stage))
    res = {"campaign": "stochastic-real (C-Trap)", "stage": "C2", "axes": rows,
           "seconds": round(time.time() - t0, 1)}
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(res, indent=1, default=str))
    print(f"\nwrote {OUT}  ({res['seconds']}s)")
    for r in rows:
        print(f"\n=== {r['record']:<8} {r['channel']}")
        if "active_scale" in r:
            a = r["active_scale"]
            print(f"    drive {a['f_drive_hz']:.2f} Hz, response share {a['response_power_share']:.3f}, "
                  f"driven={a['driven_axis']}; Rd_active {a['Rd_active_stage_units_per_V']:.4f} "
                  f"[stage/V] vs stored {r['calibration']['Rd_um_per_V']:.4f} um/V "
                  f"(ratio {a['Rd_active_over_stored']:.4f})")
        g = r["axis_gate"]
        print(f"    ACF GATE  passed={g['passed']}  theta_acf/2pi fc = "
              f"{(g['ratio'] or float('nan')):.3f}")
        b = r["band_loss"]
        print(f"    BAND LOSS r_total={b['r_total']:.3f} +- {b['r_se']:.1e}  "
              f"(nyquist {b['r_nyquist']:.3f} x diode {b['r_diode']:.3f} x anti-alias "
              f"{b['r_antialias']:.3f}); raw QV/bulk {r['qv']['raw_over_bulk']:.3f} -> "
              f"declared {r['qv']['declared_over_bulk']:.3f}")
        if "diffusion" not in r:
            print(f"    {r['verdict']}")
            continue
        for tag in ("raw", "declared"):
            d = r["diffusion"][tag]
            print(f"    DIFF {tag:<8} rows={d['n_rows']} s/band={d['median_signal_to_band']:.2f} "
                  f"mid/bulk={d['mid_over_bulk']:.3f} covers_bulk={d['covers_bulk']}")
        if "drag" in r:
            dg = r["drag"]
            print(f"    DRAG gamma_eff/bulk = {dg['gamma_eff_over_bulk']:.3f} "
                  f"[{dg['interval_over_bulk'][0]:.3f}, {dg['interval_over_bulk'][1]:.3f}]")
        for k, v in r["drift"]["lagform"].items():
            print(f"    DRIFT {k}: inverted theta / 2pi fc = {v['ratio']:.3f}")


if __name__ == "__main__":
    main()

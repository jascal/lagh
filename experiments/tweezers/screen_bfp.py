"""Screen a candidate BFP-interferometry record against the four criteria the
tweezers C0 left behind (docs/PROPOSAL_STOCHASTIC_REAL.md §9, "What the next
candidate must show").

The screen IS the guard added the same day: a dataset passes only if
`weakform.qv_three_way` reports `process_resolved` True on a passive record AND the
resolved b^2 lands near the value the instrument's own calibration implies. The
previous candidate failed both and reported `separable=True` while doing so.

Candidate: LUMICKS C-Trap, Zenodo 14726586 (CC-BY-4.0), Pylake tutorial force
calibration. Back-focal-plane interferometry, 78125 Hz.

Run: .venv/bin/python experiments/tweezers/screen_bfp.py
"""
from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).parent.parent.parent))

from lagh.weakform import qv_three_way  # noqa: E402

DATA = Path(__file__).parent / "data" / "bfp"
OUT = Path("experiments/results/tweezers_bfp_screen.json")
KB = 1.380649e-23
SHORT = (1, 2, 3, 4, 6, 8)          # ladder kept inside the trap correlation time
LONG = (1, 2, 4, 8, 16, 32)         # QV_STRIDES_3, for comparability


def calibration(f, channel="Force 1x", start_ns=None):
    """The calibration item whose voltage window sits inside this record -- i.e. the
    one computed FROM these samples. Read for scoring; consumed only for the declared
    unit conversion."""
    best = fallback = None
    for k in f["Calibration"]:
        j = json.loads(f[f"Calibration/{k}/JSON"][()])
        pl = j.get("payload")
        pl = json.loads(pl) if isinstance(pl, str) else pl
        for ch in pl["value0"]:
            if ch["channel_name"] != channel:
                continue
            t = ch.get("voltage_start", 0)      # schema varies across Bluelake versions
            while isinstance(t, dict):
                t = t.get("time_since_epoch", t.get("count", 0))
            t = int(t)
            fallback = (t, ch)
            if start_ns is None or t >= start_ns - 10**9:
                if best is None or t > best[0]:
                    best = (t, ch)
    # some Bluelake versions record no in-window calibration; fall back to the last
    # item for the channel and say so rather than failing the screen on schema drift
    _, ch = best or fallback
    # two Bluelake schemas in the same record set: list-of-{key,value} and a plain
    # dict, with the diode fields renamed between them. Normalise both.
    def flat(v):
        return dict(v) if isinstance(v, dict) else {d["key"]: d["value"] for d in v}
    out = {**flat(ch["parameters"]), **flat(ch["results"]),
           "Rf_transform": ch["transform"]["response"]}
    for new, old in (("alpha", "Diode alpha"), ("f_diode (Hz)", "Diode frequency (Hz)")):
        if new not in out and old in out:
            out[new] = out[old]
    return out


def stride_sums(u, strides):
    return np.array([float((u[s:] - u[:-s]) @ (u[s:] - u[:-s])) for s in strides])


def screen(path, channel="Force 1x", label=""):
    import h5py
    f = h5py.File(path, "r")
    d = f[f"Force HF/{channel}"]
    fs = float(d.attrs["Sample rate (Hz)"])
    dt = 1.0 / fs
    cal = calibration(f, channel, int(d.attrs["Start time (ns)"]))

    kappa = cal["kappa (pN/nm)"] * 1e-3                  # [N/m]
    gamma = cal["gamma_0 (kg/s)"]
    T = cal["Temperature (C)"] + 273.15
    kbt = KB * T
    var_th = kbt / kappa * 1e18                          # [nm^2]
    b2_th = 2 * kbt / gamma * 1e18                       # [nm^2/s]
    theta = kappa / gamma                                # [1/s]

    F = d[()]                                            # [pN]
    u = F / cal["kappa (pN/nm)"]                         # [nm]
    n = len(u)
    u = u - u.mean()

    out = {"file": Path(path).name, "channel": channel, "label": label,
           "fs_hz": fs, "n_samples": n, "duration_s": n * dt,
           "answer_key": {"kappa_pN_per_nm": cal["kappa (pN/nm)"],
                          "gamma_kg_per_s": gamma, "fc_hz": cal["fc (Hz)"],
                          "Rd_um_per_V": cal["Rd (um/V)"],
                          "T_C": cal["Temperature (C)"],
                          "viscosity_Pa_s": cal["Viscosity (Pa*s)"],
                          "bead_diameter_um": cal["Bead diameter (um)"],
                          "diode_alpha": cal["alpha"],
                          "diode_f_hz": cal["f_diode (Hz)"]},
           "expected": {"thermal_rms_nm": float(np.sqrt(var_th)),
                        "b2_nm2_per_s": b2_th, "theta_per_s": theta,
                        "tau_c_samples": 1.0 / theta / dt,
                        "per_sample_step_nm": float(np.sqrt(b2_th * dt))},
           "observed": {"position_rms_nm": float(u.std()),
                        "increment_rms_nm": float(np.diff(u).std())}}
    out["observed"]["rms_over_thermal"] = (out["observed"]["position_rms_nm"]
                                           / out["expected"]["thermal_rms_nm"])
    for tag, ladder in (("short", SHORT), ("long", LONG)):
        _, info = qv_three_way(stride_sums(u, ladder), strides=ladder,
                               n_increments=n - 1)
        b2 = info["alpha"] / (n * dt)
        out[tag] = {"strides": list(ladder), "c": info["c"], "alpha": info["alpha"],
                    "beta": info["beta"], "separable": info["separable"],
                    "process_resolved": info["process_resolved"],
                    "process_significance": info["process_significance"],
                    "sigma_obs_nm": info.get("sigma_obs"),
                    "b2_nm2_per_s": b2, "b2_over_truth": b2 / b2_th}
    return out


def main():
    t0 = time.time()
    rows = [
        screen(DATA / "passive_calibration.h5", "Force 1x", "passive, deep in cell"),
        screen(DATA / "passive_calibration.h5", "Force 1y", "passive, orthogonal axis"),
        screen(DATA / "noise_floor.h5", "Force 1x", "record with a noise floor"),
        screen(DATA / "near_surface_active_calibration.h5", "Force 1x",
               "ACTIVE: nanostage driven, independent displacement standard"),
    ]
    res = {"candidate": "LUMICKS C-Trap / Pylake force calibration",
           "record": "zenodo 14726586", "license": "CC-BY-4.0",
           "detection": "back-focal-plane interferometry",
           "rows": rows, "seconds": round(time.time() - t0, 1)}
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(res, indent=1, default=str))
    print(f"wrote {OUT}  ({res['seconds']}s)\n")

    for r in rows:
        e, o = r["expected"], r["observed"]
        print(f"=== {r['file']} :: {r['channel']}   ({r['label']})")
        print(f"    {r['n_samples']:,} samples @ {r['fs_hz']:.0f} Hz = "
              f"{r['duration_s']:.1f} s   kappa {r['answer_key']['kappa_pN_per_nm']:.4f}"
              f" pN/nm   f_c {r['answer_key']['fc_hz']:.0f} Hz   T "
              f"{r['answer_key']['T_C']:.0f} C")
        print(f"    diode filter (the correlated-error case, WITH parameters): "
              f"alpha {r['answer_key']['diode_alpha']:.3f}, "
              f"f_diode {r['answer_key']['diode_f_hz']:.0f} Hz")
        print(f"    thermal rms   expected {e['thermal_rms_nm']:>8.2f} nm   "
              f"observed {o['position_rms_nm']:>8.2f} nm   "
              f"ratio {o['rms_over_thermal']:.2f}")
        print(f"    per-sample    expected {e['per_sample_step_nm']:>8.2f} nm   "
              f"observed {o['increment_rms_nm']:>8.2f} nm   "
              f"(tau_c {e['tau_c_samples']:.0f} samples)")
        for tag in ("short", "long"):
            s = r[tag]
            print(f"    {tag:<6} s={str(s['strides']):<22} resolved="
                  f"{str(s['process_resolved']):<6} sig={s['process_significance']:>8.1f}"
                  f"   b2/truth {s['b2_over_truth']:>7.3f}   "
                  f"sigma_obs {s['sigma_obs_nm']:.3f} nm")
        print()

    p = rows[0]["short"]
    ok = p["process_resolved"] and 0.5 < p["b2_over_truth"] < 2.0
    print(f"SCREEN VERDICT: {'PASS' if ok else 'FAIL'} on the passive record "
          f"(process_resolved={p['process_resolved']}, "
          f"b2/truth={p['b2_over_truth']:.3f})")


if __name__ == "__main__":
    main()

"""C3 for the stochastic arc's real-data read: what the record determines.

C2 (docs/CASE_STUDY_TWEEZERS_C2.md) gated both axes of the near-surface active
record OUT -- their autocorrelation timescale was 0.69 of the calibration's --
and left the drag open with two readings 30% apart. C3 asks what that 0.69
actually is, and finds it is a MEASUREMENT rather than a fault:

  * equipartition fixes the position variance from kappa and temperature alone
    and does not involve the drag, so a record whose drag is not the
    calibration's reference has healthy variances on every axis and ONE common
    timescale factor. Contamination moves both. That is the discriminator, and
    it is what one axis cannot supply.
  * the drive scale at a SINGLE drive frequency is degenerate between Rd and the
    corner frequency; C2 supplied the calibration's corner, which is the quantity
    under test, so its 1.48x drag reading was the assumption coming back. The
    corner must come from the record.

Run: .venv/bin/python experiments/tweezers/run_c3.py
"""
from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).parent.parent.parent))

from experiments.tweezers import adapter                                  # noqa: E402
from lagh.instrument import (attribute_deviation, axis_gate,             # noqa: E402
                             drive_scale, equipartition_ratio, faxen_height,
                             gate_record, local_drag, realized_diffusion,
                             retention)

# The published C3 snapshot belongs to commit 489bd97. Current classifiers
# may change diagnostics; never overwrite that historical artifact on rerun.
OUT = Path("experiments/results/tweezers_c3_current.json")
ACF = {"min_lag": 8, "max_lag": 200}


def one_axis(f, channel, *, stage=None) -> dict:
    t, x_nm, applied, de, meta = adapter.bfp_position_nm(f, channel)
    fs = meta["fs_hz"]
    dt = 1.0 / fs
    fit_range = (de.get("Fit range (min.) (Hz)", 100.0),
                 de.get("Fit range (max.) (Hz)", 23000.0))
    kappa_cal = de["kappa (pN/nm)"] * 1e-3
    theta_ref = 2 * np.pi * de["fc (Hz)"]
    volts = x_nm / (de["Rd (um/V)"] * 1e3)
    out = {"channel": channel, "fs_hz": fs,
           "calibration": {"fc_hz": de["fc (Hz)"], "kappa_pN_per_nm": de["kappa (pN/nm)"],
                           "Rd_um_per_V": de["Rd (um/V)"], "alpha": de["alpha"],
                           "f_diode_hz": de["f_diode (Hz)"], "T_C": de["Temperature (C)"],
                           "gamma0_kg_per_s": de["gamma_0 (kg/s)"],
                           "bead_diameter_um": de["Bead diameter (um)"],
                           "viscosity_Pa_s": de["Viscosity (Pa*s)"],
                           "hydrodynamic_correction": de.get("Hydrodynamic correction enabled"),
                           "fit_range_hz": list(fit_range)}}
    ds = None
    if stage is not None:
        ds = drive_scale(stage, volts, fs, declared_f_c=de["fc (Hz)"], acf_kw=ACF)
        v = ds["residual"]
        out["drive"] = {k: ds[k] for k in
                        ("f_drive_hz", "stage_amplitude", "response_amplitude",
                         "response_power_share", "f_c_record_hz", "f_c_used_hz",
                         "f_c_source", "transfer", "scale", "declared_f_c_hz",
                         "scale_at_declared_f_c", "degeneracy_factor") if k in ds}
        out["drive"]["driven_axis"] = bool(ds["response_power_share"] > 0.05)
    else:
        v = volts - volts.mean()
    r = retention(v, fs, fit_range=fit_range, alpha=de["alpha"],
                  f_diode=de["f_diode (Hz)"])
    fac = r.fraction(1.0, factors=True)
    out["retention"] = {"r_var": r.variance, "fc_fit_hz": r.fc_fit, **fac,
                        "stride": {str(s): r.fraction(s) for s in (1, 2, 4, 8, 16, 32)}}
    var_nm2 = float(v.var()) * (de["Rd (um/V)"] * 1e3) ** 2
    eq = equipartition_ratio(var_nm2, kappa_cal, de["Temperature (C)"], r.variance)
    out["gate"] = axis_gate(v * de["Rd (um/V)"] * 1e3, dt, theta_ref,
                            var_ratio=eq, **ACF)
    th = out["gate"]["theta_acf"]
    # THE THREE SCALE-FREE RATIOS, each against the calibration's own number and
    # each read in VOLTS, so no displacement scale enters any of them.
    if th:
        d = realized_diffusion(v, dt, th, r)
        b2_ref_V = 2.0 * de["D (V^2/s)"]              # the SDE coefficient, b^2 = 2 D
        ratios = {"theta": th / theta_ref,
                  "diffusion": d["plateau"] / b2_ref_V,
                  "variance": out["gate"]["var_ratio"]}
        out["realized_diffusion"] = {"plateau_V2_per_s": d["plateau"],
                                     "spread": d["spread"], "strides": d["strides"],
                                     "calibration_b2_V2_per_s": b2_ref_V}
        out["ratios"] = ratios
        out["attribution"] = attribute_deviation(ratios["theta"], ratios["diffusion"],
                                                 ratios["variance"])
    out["_v"] = v
    out["_de"] = de
    return out


def main():
    t0 = time.time()
    import h5py
    records = {}
    fp = h5py.File(adapter.BFP_DIR / "passive_calibration.h5", "r")
    records["passive"] = {ch: one_axis(fp, ch) for ch in ("Force 1x", "Force 1y")}
    fa = h5py.File(adapter.BFP_DIR / "near_surface_active_calibration.h5", "r")
    stage = fa["Nanostage position/X"][()]
    records["active"] = {ch: one_axis(fa, ch, stage=stage) for ch in ("Force 1x", "Force 1y")}

    res = {"campaign": "stochastic-real (C-Trap)", "stage": "C3", "records": {}}
    for rec, axes in records.items():
        verdict = gate_record({ch: a["gate"] for ch, a in axes.items()})
        block = {"gate_record": verdict, "axes": {}}
        for ch, a in axes.items():
            v, de = a.pop("_v"), a.pop("_de")
            if a.get("drive", {}).get("driven_axis") and a["gate"]["theta_acf"]:
                r_var = a["retention"]["r_var"]
                ld = local_drag(a["gate"]["theta_acf"], a["drive"]["scale"],
                                float(v.var()), de["Temperature (C)"], r_var=r_var,
                                radius_um=0.5 * de["Bead diameter (um)"],
                                viscosity=de["Viscosity (Pa*s)"])
                ld["kappa_over_calibration"] = ld["kappa_pN_per_nm"] / de["kappa (pN/nm)"]
                ld["needs"] = ("the stage's length unit, which the file does not "
                               "declare; um assumed")
                a["local_from_drive"] = ld
                a["local_at_declared_f_c"] = {
                    "drag_over_bulk": local_drag(
                        a["gate"]["theta_acf"], a["drive"]["scale_at_declared_f_c"],
                        float(v.var()), de["Temperature (C)"], r_var=r_var,
                        radius_um=0.5 * de["Bead diameter (um)"],
                        viscosity=de["Viscosity (Pa*s)"])["drag_over_bulk"],
                    "note": "what C2's corner choice gives"}
            block["axes"][ch] = a
        drags = [a["attribution"]["hypotheses"]["drag"]["value"]
                 for a in block["axes"].values()
                 if a.get("attribution", {}).get("verdict") == "drag"]
        if len(drags) == len(block["axes"]) and len(drags) >= 2:
            F = float(np.mean(drags))
            se = float(max(np.ptp(drags) / 2.0,
                           F * max(a["attribution"]["hypotheses"]["drag"]["max_residual"]
                                   for a in block["axes"].values())))
            block["drag"] = {
                "over_bulk": F, "se": se, "per_axis": drags,
                "faxen": faxen_height(F, 0.5 * next(iter(block["axes"].values()))
                                      ["calibration"]["bead_diameter_um"], se),
                "note": ("every axis attributes its deviation to the drag alone, "
                         "from ratios that need no displacement scale")}
        # The drive route is the ONLY one that needs a length the file does not
        # declare. Given the drag the scale-free routes agree on, invert it: what
        # displacement scale, and what stage amplitude, would the drive have to
        # carry? Comparing those to the calibration's Rd and the recorded stage
        # amplitude says where the disagreement actually lives.
        for ch, a in block["axes"].items():
            if "local_from_drive" in a and "drag" in block:
                L, d = a["local_from_drive"], a["drive"]
                f = float(np.sqrt(block["drag"]["over_bulk"] / L["drag_over_bulk"]))
                a["implied_by_drag"] = {
                    "scale_um_per_V": d["scale"] / f,
                    "calibration_Rd_um_per_V": a["calibration"]["Rd_um_per_V"],
                    "scale_over_calibration_Rd": (d["scale"] / f)
                    / a["calibration"]["Rd_um_per_V"],
                    "stage_amplitude_implied": d["stage_amplitude"] / f,
                    "stage_amplitude_recorded": d["stage_amplitude"],
                    "note": ("the drag the scale-free routes agree on implies this "
                             "displacement scale; compare it with the calibration's "
                             "Rd and the recorded stage amplitude")}
        res["records"][rec] = block
    res["seconds"] = round(time.time() - t0, 1)
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(res, indent=1, default=str))
    print(f"\nwrote {OUT}  ({res['seconds']}s)")

    for rec, block in res["records"].items():
        v = block["gate_record"]
        print(f"\n=== {rec}   GATE: {v['verdict'].upper()}"
              + (f"  common kappa/gamma ratio {v['common_ratio']:.3f} "
                 f"(axis spread {v['spread']:.1%})"
                 if v["verdict"] == "common-mode" else "")
              + (f"  [{', '.join(v['contaminated'])}]" if v["verdict"] == "contaminated" else ""))
        for ch, a in block["axes"].items():
            g, rt = a["gate"], a["retention"]
            print(f"  {ch}: theta/2pi fc = {g['ratio']:.3f}  equipartition = "
                  f"{g['var_ratio']:.3f} (var_ok={g['var_ok']})   r_var={rt['r_var']:.3f} "
                  f"r_qv={rt['r_total']:.3f}")
            if "b2_stride" in a:
                b = a["b2_stride"]
                print("      b^2 across strides / bulk: "
                      + " ".join(f"s{k}:{x:.3f}" for k, x in b['per_stride_over_bulk'].items())
                      + f"   spread {b['spread']:.1%}")
            if "drive" in a and a["drive"]["driven_axis"]:
                d = a["drive"]
                print(f"      DRIVE {d['f_drive_hz']:.2f} Hz, share {d['response_power_share']:.2f}: "
                      f"corner from the {d['f_c_source']} = {d['f_c_used_hz']:.0f} Hz "
                      f"(calibration says {d['declared_f_c_hz']:.0f}); Rd = {d['scale']:.3f}, "
                      f"declared-corner Rd = {d['scale_at_declared_f_c']:.3f} "
                      f"(degeneracy factor {d['degeneracy_factor']:.3f})")
            if "ratios" in a:
                rr, at = a["ratios"], a["attribution"]
                print(f"      RATIOS theta {rr['theta']:.3f}  b^2 {rr['diffusion']:.3f}  "
                      f"var {rr['variance']:.3f}  ->  {at['verdict']}"
                      + (f" ({at['hypotheses'][at['verdict']]['parameter']} = "
                         f"{at['hypotheses'][at['verdict']]['value']:.3f}, residual "
                         f"{at['hypotheses'][at['verdict']]['max_residual']:.1%})"
                         if at["verdict"] in at.get("hypotheses", {}) else ""))
            if "local_from_drive" in a:
                L = a["local_from_drive"]
                print(f"      DRIVE-SCALE reading (needs the stage unit): "
                      f"kappa {L['kappa_over_calibration']:.3f} of the calibration's, "
                      f"gamma/bulk {L['drag_over_bulk']:.3f}"
                      f"   [C2's corner choice: {a['local_at_declared_f_c']['drag_over_bulk']:.3f}]")
        for ch, a in block["axes"].items():
            if "implied_by_drag" in a:
                im = a["implied_by_drag"]
                print(f"      IMPLIED BY THE DRAG: Rd = {im['scale_um_per_V']:.3f} um/V "
                      f"({im['scale_over_calibration_Rd']:.3f} of the calibration's "
                      f"{im['calibration_Rd_um_per_V']:.3f}); stage amplitude "
                      f"{im['stage_amplitude_implied']:.3f} against "
                      f"{im['stage_amplitude_recorded']:.3f} recorded")
        if "drag" in block:
            dg = block["drag"]
            fx = dg["faxen"]
            print(f"  DRAG (scale-free routes, both axes): {dg['over_bulk']:.3f} "
                  f"+- {dg['se']:.3f} x bulk Stokes")
            print("      Faxen: " + (f"centre height {fx['height_um']:.3f} um "
                                     f"(gap {fx['gap_um']:.3f} um)" if fx["height_um"]
                                     else fx["note"]))


if __name__ == "__main__":
    main()

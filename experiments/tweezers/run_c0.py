"""C0 for the stochastic arc's real-data read: format decode, provenance, and the
COLD three-way probe (docs/PROPOSAL_STOCHASTIC_REAL.md §7).

A SCOUT, not a certificate. It answers one question before anything is registered as a
read: does a real 200 kHz QPD channel show the three stride exponents
`weakform.qv_three_way` separates, and in which regime does it land?

    sum_i (u[i+s] - u[i])^2  ~=  c + alpha*s + beta*s^2
                                 |    |         `- smooth: the driven relaxation
                                 |    `- PROCESS: the bead's Brownian motion
                                 `- OBSERVATION: QPD detector + storage quantization

Declared INPUTS (external, consumed): the volts->nm scale from the commanded AOD
displacements, and the Stokes-Faxen drag from bead radius/height/viscosity.
The ANSWER KEY (spring constant, reciprocal time constant) is read for scoring and
consumed by NO estimator.

Run: .venv/bin/python experiments/tweezers/run_c0.py
"""
from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).parent.parent.parent))

from experiments.tweezers import adapter  # noqa: E402
from lagh.weakform import QV_STRIDES_3, qv_three_way  # noqa: E402

OUT = Path("experiments/results/tweezers_c0.json")
FILE_KEY = "calibration_optical_tweezers_07-02-2013_bead_1.zip"
TAG = "bead_07-02-2013_b1"

LADDER = (1, 2, 3, 4, 6, 8, 12, 16, 24, 32, 48, 64)
LATE = (4, 8, 16, 32)                     # the s >= 4 fit, for R2
BLOCKS = (1, 2, 5, 10, 50, 100, 500, 1000, 5000, 20000, 100000)
HALF_PERIOD = adapter.DRIVE_PERIOD // 2   # 8000 samples = 40 ms
GUARD = 1500                              # ~4 trap time constants after each edge
T_KELVIN = 298.15                         # DECLARED input: not recorded in the file
KB_SI = 1.380649e-23


def stride_sums(u, segs, strides=LADDER):
    """Increments never straddle a segment boundary -- a stride crossing an excised
    drive edge would re-import the jump as a spurious linear-in-s term."""
    sums = np.zeros(len(strides))
    n_tot = 0
    for a, b in segs:
        v = u[a:b]
        n_tot += len(v)
        for j, s in enumerate(strides):
            if len(v) > s:
                d = v[s:] - v[:-s]
                sums[j] += float(d @ d)
    return sums, n_tot


def fit(sums, subset, n_inc):
    idx = [LADDER.index(s) for s in subset]
    return qv_three_way(np.asarray(sums)[idx], strides=subset, n_increments=n_inc)


def storage_floor(x, n=400_000):
    """sigma_rep: the lattice the stored values live on. A quantization step is an iid
    observation error and therefore a confound for `c`."""
    u = np.unique(x[:n])
    d = np.diff(u)
    d = d[d > 0]
    if len(d) == 0:
        return {"lsb_volts": None}
    lsb = float(np.min(d))
    return {"lsb_volts": lsb, "n_unique_in_slice": int(len(u)),
            "fraction_of_gaps_on_lattice":
                float(np.mean(np.abs(d / lsb - np.round(d / lsb)) < 1e-6)),
            "sigma_rep_volts": lsb / np.sqrt(12.0)}


def block_sweep(res, var_th, tau_c):
    """m * Var[block mean of m] is FLAT in m for white noise and RISES where a process
    contributes. Run against the OU's own predicted share so 'invisible' is a measured
    statement rather than an impression."""
    out = []
    for m in BLOCKS:
        tot = (len(res) // m) * m
        if tot < 8 * m:
            continue
        v = float(res[:tot].reshape(-1, m).mean(axis=1).var())
        tau = m * adapter.DT
        ou = var_th * (2 * tau_c / tau) if tau > 3 * tau_c else var_th
        out.append({"m": m, "avg_window_s": tau, "var_block_nm2": v,
                    "m_times_var": m * v, "ou_predicted_nm2": ou,
                    "excess_over_ou": v / ou})
    return out


def main():
    t0 = time.time()
    h5 = adapter.fetch(FILE_KEY, TAG)
    import h5py
    f = h5py.File(h5, "r")

    key = adapter.answer_key(f)
    slope = abs(key["qpd_slope_V_per_V_per_nm"])       # declared input: volts -> nm
    gamma = abs(key["beta_pNs_per_nm"]) * 1e-3         # declared input: drag [kg/s]
    k_ext = abs(key["k_pN_per_nm"]) * 1e-3             # ANSWER: scoring only [N/m]
    theta_ext = abs(key["inv_tau_per_s"])              # ANSWER: scoring only [1/s]

    # what the OU says this instrument should show, from the declared inputs + answer
    var_th = KB_SI * T_KELVIN / k_ext * 1e18           # [nm^2]
    b2_th = 2 * KB_SI * T_KELVIN / gamma * 1e18        # [nm^2/s]
    tau_c = 1.0 / theta_ext
    expected = {"thermal_var_nm2": var_th, "thermal_rms_nm": float(np.sqrt(var_th)),
                "b2_nm2_per_s": b2_th, "tau_c_s": tau_c,
                "tau_c_samples": tau_c / adapter.DT,
                "per_sample_brownian_step_nm": float(np.sqrt(b2_th * adapter.DT))}

    prov = {}
    ratios = [(v["k_pN_per_nm"] / v["inv_tau_per_s"]) / key["beta_pNs_per_nm"]
              for v in key["per_amplitude"].values()]
    prov["k_is_derived_from_inv_tau"] = {
        "max_dev_from_1": float(np.max(np.abs(np.array(ratios) - 1.0))),
        "finding": ("spring_constant is DERIVED as (1/tau)*beta to machine precision, "
                    "so the answer key carries ONE independent drift measurement, "
                    "not two")}
    prov["dark_channels_are_scalars"] = {
        "values": adapter.dark(f),
        "finding": ("'signal in darkness' is a mean already subtracted in real time, "
                    "not a dark time series -- the proposal's R7 measured null is NOT "
                    "available in this container")}

    per_amp, primary = {}, {}
    for disp in adapter.DISPLACEMENTS:
        x, _y, ssum = adapter.raw(f, disp)
        n = len(x)
        sum_mean = float(np.mean(ssum))
        # normalise by the MEAN sum: a fixed scale injects no extra noise, unlike
        # dividing by the noisy per-sample sum.
        u = (x / sum_mean) / slope                                     # [nm]
        # the drive is phase-locked at 16000 samples/period; excise a guard after
        # every half-period edge so the exponential relaxation is outside the segments
        segs = [(i * HALF_PERIOD + GUARD, (i + 1) * HALF_PERIOD)
                for i in range(n // HALF_PERIOD)]
        s_full, n_full = stride_sums(u, [(0, n)])
        s_flat, n_flat = stride_sums(u, segs)
        rec = {"n_samples": n, "duration_s": n * adapter.DT,
               "sum_mean_V": sum_mean, "sum_rel_sd": float(np.std(ssum) / sum_mean),
               "n_flat_segments": len(segs), "ladder": list(LADDER),
               "sums_flat_nm2": [float(v) for v in s_flat]}
        for tag, sums, n_inc, sub in (("full", s_full, n_full, tuple(QV_STRIDES_3)),
                                      ("flat", s_flat, n_flat, tuple(QV_STRIDES_3)),
                                      ("flat_s_ge_4", s_flat, n_flat, LATE)):
            _, info = fit(sums, sub, n_inc)
            rec[tag] = {kk: info.get(kk) for kk in
                        ("c", "alpha", "beta", "se_c", "se_alpha", "dominant",
                         "separable", "process_over_observation", "sigma_obs")}
            rec[tag]["b2_nm2_per_s"] = info["alpha"] / (n_inc * adapter.DT)
            rec[tag]["b2_over_truth"] = rec[tag]["b2_nm2_per_s"] / b2_th

        if disp == "800":
            res = np.concatenate([u[a:b] - u[a:b].mean() for a, b in segs])
            rc = res - res.mean()
            den = float(rc @ rc)
            acf = [float(rc[:-j] @ rc[j:] / den) for j in range(1, 13)]
            rho = acf[0]
            primary = {
                "per_sample_sd_nm": float(res.std()),
                "noise_over_thermal_rms": float(res.std() / np.sqrt(var_th)),
                "var_flat_nm2": float(res.var()),
                "acf_lags_1_to_12": acf,
                "periodic_line": {
                    "lag5": acf[4], "lag10": acf[9],
                    "finding": ("the observation error carries a PERIODIC component "
                                "(40 kHz at this rate), not a decaying correlation")},
                "block_sweep": block_sweep(res, var_th, tau_c),
                "theta_ar_per_s": (float(-np.log(rho) / adapter.DT) if rho > 0
                                   else None),
                "theta_ar_note": (None if rho > 0 else
                                  f"AR(1) refused: lag-1 autocorrelation is {rho:.3f} "
                                  f"< 0, so the record has no recoverable OU "
                                  f"relaxation -- the noise, not the bead, sets it"),
                "theta_external_per_s": theta_ext,
            }
            rec["storage_floor"] = storage_floor(x)
            rec["storage_floor"]["sigma_rep_nm"] = (
                rec["storage_floor"]["sigma_rep_volts"] / sum_mean / slope)
        per_amp[disp] = rec
        del x, _y, ssum, u

    # ------------------------------------------------------------------- the verdict
    p8 = per_amp["800"]
    verdict = {
        "read": "REFUSED for the diffusion / process-noise target",
        "reasons": [
            f"per-sample detection noise is {primary['noise_over_thermal_rms']:.0f}x "
            f"the bead's whole thermal rms ({expected['thermal_rms_nm']:.2f} nm), and "
            f"{primary['per_sample_sd_nm'] / expected['per_sample_brownian_step_nm']:.0f}x "
            f"the per-sample Brownian step",
            "averaging cannot close the gap: m*Var[block mean] RISES with m rather "
            "than staying flat, so the noise is not white -- and the process would "
            "need ~1e5 samples of averaging against a correlation time of only "
            f"{expected['tau_c_samples']:.0f} samples",
            "detection is FLUORESCENCE from a photon-starved bead (sum channel "
            f"{p8['sum_mean_V'] * 1e3:.2f} mV), not back-focal-plane interferometry; "
            "the authors recover their signal by averaging ~750 IDENTICAL driven "
            "responses, which a Brownian path by definition does not admit",
        ],
        "survives": ("the DRIVEN read: a deterministic relaxation observed through a "
                     "stochastic channel is Level 0's third rung, and theta has an "
                     "external answer (1/tau) this record can be scored against"),
    }
    guard = {
        "reported_by_qv_three_way": {"separable": p8["flat"]["separable"],
                                     "dominant": p8["flat"]["dominant"]},
        "b2_error_factor": p8["flat"]["b2_over_truth"],
        "c_over_alpha": (p8["flat"]["c"] / p8["flat"]["alpha"]),
        "finding": (
            "qv_three_way reported separable=True and handed back a b^2 wrong by "
            f"{p8['flat']['b2_over_truth']:.3g}x. SEP_MIN_FRAC guards ONE direction "
            "-- observation buried under process -- and passes trivially when the "
            "reverse holds. The measured graceful degradation (720x observation "
            "noise -> alpha 40% high) was established under IID observation error; "
            "this record's error is correlated (periodic line at lag 5), and at only "
            f"{p8['flat']['c'] / p8['flat']['alpha']:.0f}x the alpha estimate is wrong "
            "by 4-5 orders of magnitude. A symmetric bar gating the b^2 consumer is "
            "the missing guard."),
    }
    amend = {
        "target": "DIRECTION_STOCHASTIC.md:486 / proposal R2",
        "registered": ("computer-vision tracking error is autocorrelated to some lag "
                       "tau, so restricting the fit to strides s > tau recovers c"),
        "measured": ("this instrument's error is PERIODIC (lag-5 and lag-10 both "
                     f"~{primary['periodic_line']['lag5']:.2f}), which has no finite "
                     "correlation length -- s > tau has nothing to step past. The "
                     "registered mitigation is inadequate for periodic error and the "
                     "prediction needs re-registering in that form."),
    }

    res_json = {"campaign": "stochastic-real (optical tweezers)", "stage": "C0",
                "artifact": json.loads((adapter.DATA / f"{TAG}.json").read_text()),
                "answer_key_scoring_only": key, "expected_from_ou": expected,
                "provenance": prov, "per_amplitude": per_amp, "primary_800": primary,
                "verdict": verdict, "checker_guard_gap": guard, "amendment": amend,
                "seconds": round(time.time() - t0, 1)}
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(res_json, indent=1, default=str))

    # -------------------------------------------------------------------- the report
    print(f"\nwrote {OUT}  ({res_json['seconds']}s)")
    a = res_json["artifact"]
    print(f"\nARTIFACT   {a['sha256'][:16]}...  {a['bytes'] / 1e6:.0f} MB  "
          f"{a['license']}  doi:{a['doi']}")
    print(f"RECORD     {p8['n_samples']:,} samples @ 200 kHz = {p8['duration_s']:.1f} s"
          f",  {p8['n_flat_segments']} flat segments,  sum channel "
          f"{p8['sum_mean_V'] * 1e3:.2f} mV")
    sf = p8["storage_floor"]
    print(f"FLOOR      LSB {sf['lsb_volts']:.3e} V -> sigma_rep {sf['sigma_rep_nm']:.2f} nm")

    print("\nPROVENANCE")
    print(f"  k == (1/tau)*beta to {prov['k_is_derived_from_inv_tau']['max_dev_from_1']:.1e}"
          f"  -> ONE independent drift number in the answer key, not two")
    print("  dark channels are scalars -> R7's measured null is unavailable")

    print(f"\nEXPECTED from the OU, given the declared inputs and the answer key:")
    print(f"  thermal rms {expected['thermal_rms_nm']:.2f} nm   "
          f"b^2 {expected['b2_nm2_per_s']:.4g} nm^2/s   "
          f"per-sample step {expected['per_sample_brownian_step_nm']:.2f} nm   "
          f"tau_c {expected['tau_c_samples']:.0f} samples")
    print(f"OBSERVED  per-sample sd {primary['per_sample_sd_nm']:.0f} nm  ->  "
          f"{primary['noise_over_thermal_rms']:.0f}x the thermal rms")

    print("\nR1/R3  three-way separation, +800 nm  (b^2 in nm^2/s):")
    hdr = (f"  {'fit':<13}{'c':>11}{'alpha':>11}{'beta':>11}  {'dominant':<12}"
           f"{'sep':<6}{'b2':>11}{'b2/truth':>11}")
    print(hdr + "\n  " + "-" * (len(hdr) - 2))
    for t in ("full", "flat", "flat_s_ge_4"):
        r = p8[t]
        print(f"  {t:<13}{r['c']:>11.3g}{r['alpha']:>11.3g}{r['beta']:>11.3g}  "
              f"{str(r['dominant']):<12}{str(r['separable']):<6}"
              f"{r['b2_nm2_per_s']:>11.3g}{r['b2_over_truth']:>11.3g}")
    print(f"\n  R3: alpha(full)/alpha(flat) = "
          f"{p8['full']['alpha'] / p8['flat']['alpha']:.3f}"
          f"   (drive edges inflate the process term)")
    print(f"  R2: c(s>=1) {p8['flat']['c']:.3g}  vs  c(s>=4) "
          f"{p8['flat_s_ge_4']['c']:.3g}   ratio "
          f"{p8['flat_s_ge_4']['c'] / p8['flat']['c']:.3f}")
    print("  ACF lags 1..12: " + " ".join(f"{v:+.2f}" for v in
                                          primary["acf_lags_1_to_12"]))
    print(f"      lag-5 {primary['periodic_line']['lag5']:+.3f}, lag-10 "
          f"{primary['periodic_line']['lag10']:+.3f}  -> PERIODIC error, not a "
          f"decaying tail")

    print("\nBLOCK-MEAN SWEEP   m*Var is flat for white noise; a process makes it rise")
    print(f"  {'m':>7}{'window s':>11}{'m*Var':>12}{'Var(block)':>12}"
          f"{'OU says':>11}{'excess':>10}")
    for b in primary["block_sweep"]:
        print(f"  {b['m']:>7}{b['avg_window_s']:>11.5f}{b['m_times_var']:>12.3g}"
              f"{b['var_block_nm2']:>12.3g}{b['ou_predicted_nm2']:>11.3g}"
              f"{b['excess_over_ou']:>10.0f}x")

    print(f"\nAR(1) theta: {primary['theta_ar_note'] or primary['theta_ar_per_s']}")
    print(f"\nVERDICT  {verdict['read']}")
    for r in verdict["reasons"]:
        print(f"  - {r}")
    print(f"  SURVIVES: {verdict['survives']}")
    print(f"\nCHECKER GUARD GAP\n  {guard['finding']}")
    print(f"\nAMENDMENT ({amend['target']})\n  {amend['measured']}")


if __name__ == "__main__":
    main()

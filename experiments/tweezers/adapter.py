"""Zenodo optical-tweezers adapter (docs/PROPOSAL_STOCHASTIC_REAL.md): frozen-artifact
acquisition. Every fetch records the exact URL, DOI, byte count and SHA-256 before any
analysis touches the bytes -- the registration discipline applied to a remote archive,
per the Gaia and Materials adapter pattern.

NO credential of any kind. Zenodo record 3333905 is CC-0 over anonymous HTTP; that is
the reason this dataset replaced highD, whose access is granted by manual review
against an institutional address.

The HDF5 layout is decoded from the PRODUCER'S OWN script (`calibrateopticaltweezers.py`,
distributed in the same record) rather than inferred from the file, so the field
semantics below are the authors' and not ours.
"""
from __future__ import annotations

import hashlib
import json
import urllib.request
import zipfile

import numpy as np
from pathlib import Path

DATA = Path(__file__).parent / "data"
RECORD = "3333905"
DOI = "10.5281/zenodo.3333905"
CITE = ("Skidmore, B., Rajasekharan, V., & Farrell, B. (2019). Calibration of the "
        "temporal and displacement characteristics of an optical tweezers "
        "[Data set]. Zenodo. https://doi.org/10.5281/zenodo.3333905  (CC-0)")

# The instrument's fixed acquisition parameters, from the record's own description.
FS = 200_000.0                    # QPD sampling rate [Hz]
DT = 1.0 / FS                     # [s]
DRIVE_PERIOD = 16_000             # samples per AOD square-wave period (80 ms @ 200 kHz)
DISPLACEMENTS = ("800", "-800", "500", "-500")   # commanded amplitudes [nm]


def file_url(key: str) -> str:
    return f"https://zenodo.org/records/{RECORD}/files/{key}?download=1"


def fetch(key: str, tag: str) -> Path:
    """Download one archive file, freeze its provenance, unzip, return the .h5 path.
    Idempotent by tag: an existing artifact is reused rather than re-fetched."""
    DATA.mkdir(parents=True, exist_ok=True)
    zpath = DATA / f"{tag}.zip"
    meta = DATA / f"{tag}.json"
    outdir = DATA / tag
    if not zpath.exists():
        url = file_url(key)
        with urllib.request.urlopen(url, timeout=1800) as r, open(zpath, "wb") as fh:
            while chunk := r.read(1 << 20):
                fh.write(chunk)
    if not meta.exists():
        h = hashlib.sha256()
        with open(zpath, "rb") as fh:
            while chunk := fh.read(1 << 20):
                h.update(chunk)
        meta.write_text(json.dumps({
            "record": RECORD, "doi": DOI, "file": key, "url": file_url(key),
            "bytes": zpath.stat().st_size, "sha256": h.hexdigest(),
            "license": "CC-0", "cite": CITE}, indent=1))
    if not outdir.exists():
        with zipfile.ZipFile(zpath) as z:
            z.extractall(outdir)
    h5 = next(outdir.glob("*.h5"))
    return h5


# --------------------------------------------------------------------------- BFP
# The replacement candidate (docs/PROPOSAL_STOCHASTIC_REAL.md §9): LUMICKS C-Trap,
# back-focal-plane interferometry, where the fluorescence record's detection noise
# problem does not arise.
BFP_RECORD = "14726586"
BFP_DOI = "10.5281/zenodo.14726586"
BFP_DIR = DATA / "bfp"


def bfp_fetch(key: str) -> Path:
    """One Pylake .h5, with provenance frozen alongside it. Idempotent."""
    BFP_DIR.mkdir(parents=True, exist_ok=True)
    out = BFP_DIR / key
    if not out.exists():
        url = f"https://zenodo.org/records/{BFP_RECORD}/files/{key}?download=1"
        with urllib.request.urlopen(url, timeout=900) as r, open(out, "wb") as fh:
            while chunk := r.read(1 << 20):
                fh.write(chunk)
    meta = BFP_DIR / f"{key}.json"
    if not meta.exists():
        h = hashlib.sha256()
        with open(out, "rb") as fh:
            while chunk := fh.read(1 << 20):
                h.update(chunk)
        meta.write_text(json.dumps({
            "record": BFP_RECORD, "doi": BFP_DOI, "file": key,
            "bytes": out.stat().st_size, "sha256": h.hexdigest(),
            "license": "CC-BY-4.0",
            "cite": "LUMICKS. Pylake tutorial dataset: Force calibration. Zenodo. "
                    "https://doi.org/10.5281/zenodo.14726586"}, indent=1))
    return out


def _ts(v) -> int:
    """Timestamps appear as a bare int or nested under time_since_epoch/count,
    depending on the Bluelake version that wrote the file."""
    while isinstance(v, dict):
        v = v.get("time_since_epoch", v.get("count", 0))
    return int(v or 0)


def bfp_calibrations(f, channel: str) -> list:
    """Every calibration item for `channel`, normalised across the two Bluelake
    schemas (list-of-{key,value} and plain dict) and the diode-field rename."""
    def flat(v):
        return dict(v) if isinstance(v, dict) else {d["key"]: d["value"] for d in v}
    out = []
    for k in f["Calibration"]:
        j = json.loads(f[f"Calibration/{k}/JSON"][()])
        pl = j.get("payload")
        pl = json.loads(pl) if isinstance(pl, str) else pl
        for ch in pl["value0"]:
            if ch["channel_name"] != channel:
                continue
            rec = {**flat(ch["parameters"]), **flat(ch["results"]),
                   "item": k, "Rf_transform": abs(ch["transform"]["response"]),
                   "voltage_start": _ts(ch.get("voltage_start")),
                   "voltage_stop": _ts(ch.get("voltage_stop")),
                   "conversion_start": _ts(ch.get("conversion_start"))}
            for new, old in (("alpha", "Diode alpha"),
                             ("f_diode (Hz)", "Diode frequency (Hz)")):
                if new not in rec and old in rec:
                    rec[new] = rec[old]
            out.append(rec)
    return sorted(out, key=lambda r: r["conversion_start"])


def bfp_position_nm(f, channel: str):
    """(t, x_nm, applied, derived, meta) for one high-frequency force channel.

    THE UNIT CHAIN IS NOT ONE CALIBRATION (C1 provenance finding). The stored pN
    values were produced with whichever calibration was ACTIVE at acquisition
    (`conversion_start` <= record start), while the calibration DERIVED from this
    record -- the one whose voltage window lies inside it, and the only one whose
    kappa describes this trap -- is a later item. Dividing the stored pN by the
    derived kappa silently mixes the two and costs ~5% in position, ~9% in b^2.
    The honest chain is  pN -> V with the APPLIED Rf,  V -> nm with the DERIVED Rd.
    """
    d = f[f"Force HF/{channel}"]
    fs = float(d.attrs["Sample rate (Hz)"])
    t0, t1 = int(d.attrs["Start time (ns)"]), int(d.attrs["Stop time (ns)"])
    cals = bfp_calibrations(f, channel)
    applied = [c for c in cals if c["conversion_start"] <= t0]
    applied = applied[-1] if applied else cals[0]
    inside = [c for c in cals if t0 <= c["voltage_start"] <= t1]
    derived = inside[-1] if inside else applied
    volts = d[()] / applied["Rf_transform"]
    x_nm = volts * derived["Rd (um/V)"] * 1e3
    t = np.arange(len(x_nm)) / fs
    meta = {"fs_hz": fs, "n": len(x_nm), "duration_s": len(x_nm) / fs,
            "applied_item": applied["item"], "derived_item": derived["item"],
            "same_item": applied["item"] == derived["item"],
            "naive_kappa_route_error": float(
                applied["Rf_transform"] / derived["Rf_transform"] - 1.0)}
    return t, x_nm, applied, derived, meta


def diode_deconvolve(u, alpha: float, f_diode: float, dt: float):
    """Invert the photodiode's parasitic filtering.

    The standard model (Berg-Sorensen & Flyvbjerg) states the measured PSD as
    `g(f) = alpha^2 + (1 - alpha^2)/(1 + (f/f_diode)^2)`, and that is exactly the
    squared modulus of

        H(f) = alpha + (1 - alpha) / (1 + i f / f_diode)

    i.e. in the time domain the detector reports `alpha` of the true position
    instantly and passes the rest through a first-order low-pass. So the distortion
    is invertible in closed form rather than merely characterisable, and the
    correction consumes only DECLARED inputs (alpha, f_diode) that the instrument
    fitted independently of anything certified here.

    Discretising the low-pass with a = exp(-2 pi f_diode dt):

        y[n] = a y[n-1] + (1-a) x[n]              (the filtered part)
        u[n] = alpha x[n] + (1-alpha) y[n]        (what is recorded)

    eliminating x gives a stable first-order recursion for y, and x follows.
    """
    from scipy.signal import lfilter
    a = float(np.exp(-2 * np.pi * f_diode * dt))
    den = alpha + (1 - alpha) * (1 - a)
    c = (1 - alpha) * a
    A = a - (1 - a) * c / den
    B = (1 - a) / den
    y = lfilter([B], [1.0, -A], np.asarray(u, float))
    y_prev = np.concatenate([[0.0], y[:-1]])
    return (np.asarray(u, float) - c * y_prev) / den, {"a": a, "pole": A}


def answer_key(f) -> dict:
    """The EXTERNAL values, produced by measurements that never touch the passive
    fluctuations we certify against them (proposal §3). Read for reporting only --
    no estimator in this campaign may consume this dict.

      beta  Stokes-Faxen drag  [pN*s/nm]  from bead radius, height above dish, viscosity
      k     spring constant    [pN/nm]    from the DRIVEN step response
      inv_tau  reciprocal time constant [1/s], ditto, with its across-amplitude sd
      qpd_slope  [V/V/nm]                 from the COMMANDED AOD displacements
    """
    cal = f["calibration"]
    p = f["mean_background_corrected_datum"][
        "optical_tweezer_calibration_parameters_for_experiment"]
    return {
        "beta_pNs_per_nm": float(cal["Stokes_Faxen_coeffiecent"][()]),
        "gain": float(cal["gain"][()]),
        "power_mW": float(cal["power at microscope objective"][()]),
        "k_pN_per_nm": float(p["mean_spring_constant"][()]),
        "inv_tau_per_s": float(p["mean_reciprocal_time_constant"][()]),
        "inv_tau_sd_per_s": float(p["mean_reciprocal_time_constant_sd"][()]),
        "qpd_slope_V_per_V_per_nm": float(p["QPD_slope_of_line"][()]),
        "per_amplitude": {d: {
            "k_pN_per_nm": float(f["mean_background_corrected_datum"][d]
                                 ["spring_constant"][()]),
            "inv_tau_per_s": float(f["mean_background_corrected_datum"][d]
                                   ["reciprocal_time_constant"][()]),
        } for d in DISPLACEMENTS if d in f["mean_background_corrected_datum"]},
    }


def raw(f, displacement: str):
    """(x, y, sum) raw QPD channels [V] at 200 kHz for one commanded amplitude."""
    g = f["measurement_datum"]["planned_displacement_of_trapped_bead"][displacement]
    return (g["displacement_x"][()], g["displacement_y"][()],
            g["sum_signal_in_light"][()])


def dark(f) -> dict:
    """The 'signal in darkness' fields. NOTE (C0 finding): these are SCALARS -- the
    mean dark level, already subtracted in real time -- not a dark time series. The
    proposal's R7 null is therefore not available in this container as written."""
    g = f["measurement_datum"]["mean_background_signal_in_darkness"]
    return {k: float(g[k][()]) for k in g}

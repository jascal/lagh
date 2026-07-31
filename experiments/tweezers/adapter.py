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

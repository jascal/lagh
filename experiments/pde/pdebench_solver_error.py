"""Measure PDEBench's own numerical error where an exact solution exists.

`docs/PDEBENCH_READINESS.md` declaration 2 says the solver error of a shipped
field cannot be measured after the fact, because the coarse levels a tolerance
ladder would need were never saved. That is true in general -- but not for
1-D ADVECTION, where the stated law has an exact solution: u(x, t) = u0(x - beta t),
a circular shift, computable spectrally to machine precision from the file's own
first time slice.

So for that one family the undeclared term can be MEASURED rather than assumed,
and the number is worth having, because it sets the floor for every certificate
against these files: a band assembled from float32 storage noise alone would be
certifying against numbers whose error was never declared -- the exact failure
mode the weak-form arc exists to prevent.

Reports the deviation as a function of time, so its GROWTH is visible: storage
rounding is flat in t, a solver's error is not.
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).parent.parent.parent))


def exact_shift(u0, x, distance):
    """u0 translated by `distance` on a periodic grid, spectrally exact."""
    n = len(x)
    L = n * float(x[1] - x[0])
    F = np.fft.rfft(u0)
    k = np.fft.rfftfreq(n, d=1.0 / n)
    return np.fft.irfft(F * np.exp(-2j * np.pi * k * distance / L), n)


def advection_error(path, beta, *, samples=None):
    import h5py
    with h5py.File(path, "r") as h:
        ten = h["tensor"]
        n_s = ten.shape[0] if samples is None else min(samples, ten.shape[0])
        x = np.asarray(h["x-coordinate"][()], float)
        t = np.asarray(h["t-coordinate"][()], float)
        data = np.asarray(ten[:n_s], float)
    n_t = data.shape[1]
    t = t[:n_t]
    per_t = np.zeros(n_t)
    scale = 0.0
    for s in range(len(data)):
        u = data[s]
        scale = max(scale, float(np.max(np.abs(u))))
        for j in range(n_t):
            pred = exact_shift(u[0], x, beta * (t[j] - t[0]))
            per_t[j] = max(per_t[j], float(np.max(np.abs(u[j] - pred))))
    return {"path": str(path), "beta": beta, "n_samples": int(len(data)),
            "field_scale": scale,
            "max_abs_error": float(per_t.max()),
            "max_rel_error": float(per_t.max() / scale),
            "error_vs_time": [[float(t[j]), float(per_t[j])]
                              for j in range(0, n_t, max(1, n_t // 10))],
            "note": ("deviation from the EXACT solution of the file's own "
                     "stated law, computed spectrally from its first time "
                     "slice; it grows with t, which storage rounding cannot")}


def reference_error(path, law, *, samples=None, nsub=64):
    """The same measurement for any family whose stated law this program can
    integrate: re-solve from the file's OWN first time slice with the spectral
    exponential integrator, at two substep counts so the reference's own error
    is bounded, and report the deviation.

    The reference is independent of the shipped solver, and its error is
    declared rather than assumed -- if the reference ladder does not converge
    well below the deviation, the measurement is refused instead of reported."""
    import h5py

    from experiments.pde.verify import integrate
    with h5py.File(path, "r") as h:
        ten = h["tensor"]
        n_s = ten.shape[0] if samples is None else min(samples, ten.shape[0])
        x = np.asarray(h["x-coordinate"][()], float)
        t = np.asarray(h["t-coordinate"][()], float)
        data = np.asarray(ten[:n_s], float)
    n_t = data.shape[1]
    t = t[:n_t]
    per_t = np.zeros(n_t)
    scale, ref_err = 0.0, 0.0
    for s in range(len(data)):
        u = data[s]
        scale = max(scale, float(np.max(np.abs(u))))
        a = integrate(u[0], x, t, law, scheme="etd", nsub=nsub)
        b = integrate(u[0], x, t, law, scheme="etd", nsub=2 * nsub)
        if a is None or b is None:
            return {"refusal": "reference solver did not run"}
        ref_err = max(ref_err, float(np.max(np.abs(a - b))))
        per_t = np.maximum(per_t, np.max(np.abs(u.T - b), axis=0))
    if ref_err > 0.1 * per_t.max():
        return {"refusal": "the reference's own ladder error is not small "
                           "against the deviation being measured",
                "reference_ladder_error": ref_err,
                "deviation": float(per_t.max())}
    return {"path": str(path), "law": law, "n_samples": int(len(data)),
            "field_scale": scale, "reference_ladder_error": ref_err,
            "max_abs_error": float(per_t.max()),
            "max_rel_error": float(per_t.max() / scale),
            "error_vs_time": [[float(t[j]), float(per_t[j])]
                              for j in range(0, n_t, max(1, n_t // 10))],
            "note": ("deviation from an INDEPENDENT high-accuracy solve of the "
                     "file's own stated law, started from the file's first "
                     "time slice; the reference's own error is bounded by its "
                     "substep ladder and reported")}


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("path")
    ap.add_argument("--beta", type=float, default=None,
                    help="advection speed: uses the EXACT shift solution")
    ap.add_argument("--law", default=None,
                    help="JSON law, e.g. '{\"u*u_x\": -1, \"u_xx\": 0.0032}': "
                         "uses an independent high-accuracy solve as reference")
    ap.add_argument("--nsub", type=int, default=64)
    ap.add_argument("--samples", type=int, default=None)
    a = ap.parse_args(argv)
    if a.beta is not None:
        r = advection_error(a.path, a.beta, samples=a.samples)
    elif a.law:
        r = reference_error(a.path, json.loads(a.law), samples=a.samples,
                            nsub=a.nsub)
    else:
        ap.error("one of --beta or --law is required")
    print(json.dumps(r, indent=1))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

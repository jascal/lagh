"""NewtonBench vanilla-equation oracle for lagh's active loop.

Queried offline (no API key) via each module's run_experiment_for_module at
noise_level=0, system='vanilla_equation'. Law bodies are NOT read by the learner;
the source is used only for the dense-grid reference at scoring time.

Input names are read from the interface signatures (the wyly N1 work), NOT the law
bodies. Domain boxes are each function's declared positive-real safe range.
"""

from __future__ import annotations

import importlib
import sys
from pathlib import Path

import numpy as np

NB = Path("/home/allans/code/NewtonBench")
if str(NB) not in sys.path:
    sys.path.insert(0, str(NB))

# module -> (vanilla input names, box_lo, box_hi). Names from signatures only.
MODULES = {
    "m0_gravity": (["mass1", "mass2", "distance"], [0.5, 0.5, 0.5], [10, 10, 10]),
    "m1_coulomb_force": (["q1", "q2", "distance"], [0.5, 0.5, 0.5], [10, 10, 10]),
    "m2_magnetic_force": (["current1", "current2", "distance"],
                          [0.5, 0.5, 0.5], [10, 10, 10]),
    "m3_fourier_law": (["k", "A", "delta_T", "d"], [0.5]*4, [10]*4),
    "m4_snell_law": (["refractive_index_1", "refractive_index_2", "incidence_angle"],
                     [1.0, 1.0, 10.0], [2.0, 2.0, 80.0]),   # angle in DEGREES
    "m5_radioactive_decay": (["N0", "lambda_constant", "t"], [1.0, 0.05, 0.1],
                             [100, 1.0, 10]),
    "m6_underdamped_harmonic": (["k_constant", "mass", "b_constant"],
                                [2.0, 1.0, 0.1], [10, 5, 1.0]),
    "m7_malus_law": (["I_0", "theta"], [1.0, 0.05], [10, 1.5]),
    "m8_sound_speed": (["adiabatic_index", "temperature", "molar_mass"],
                       [1.1, 100, 0.01], [1.7, 500, 0.05]),
    "m9_hooke_law": (["x"], [0.1], [10]),
    "m10_be_distribution": (["omega", "temperature"], [50, 50], [500, 500]),
    "m11_heat_transfer": (["m", "c", "delta_T"], [0.5, 0.5, 1.0], [10, 10, 100]),
}


def _runner(module: str):
    return importlib.import_module(f"modules.{module}.core").run_experiment_for_module


def make_oracle(module: str, law_version: str, difficulty: str = "easy"):
    inputs = MODULES[module][0]
    run = _runner(module)

    def oracle(X: np.ndarray) -> np.ndarray:
        X = np.atleast_2d(np.asarray(X, float))
        out = []
        for row in X:
            kw = {k: float(v) for k, v in zip(inputs, row)}
            try:
                y = float(run(noise_level=0.0, difficulty=difficulty,
                              system="vanilla_equation", law_version=law_version, **kw))
            except Exception:                                # noqa: BLE001
                y = float("nan")
            # values at/beyond the float64 exact-integer ceiling (2**52) are
            # overflow/precision saturation artifacts, not signal -> invalid
            if abs(y) >= 2.0**52:
                y = float("nan")
            out.append(y)
        return np.array(out)

    return oracle


def available_versions(module: str, difficulty: str = "easy"):
    return importlib.import_module(f"modules.{module}.laws") \
        .get_available_law_versions(difficulty)

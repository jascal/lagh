"""System discoverer (H3b): equations, invariants, shared constants."""
import numpy as np

from lagh.systems import discover_invariants, discover_system


def _chain(n=300, k1=0.7, k2=1.3, seed=5):
    rng = np.random.default_rng(seed)
    t = np.sort(rng.uniform(0, 4, n))
    A0 = 3.0
    A = A0 * np.exp(-k1 * t)
    B = A0 * k1 / (k2 - k1) * (np.exp(-k1 * t) - np.exp(-k2 * t)) + 0.5 * np.exp(-k2 * t)
    C = A0 + 0.5 + 0.2 - A - B + 0.0 * t + 0.2 * 0
    return {"A": A, "B": B, "C": C,
            "dA_dt": -k1 * A, "dB_dt": k1 * A - k2 * B, "dC_dt": k2 * B}


def test_invariant_conservation_found():
    data = _chain()
    invs = discover_invariants(data)
    joined = " | ".join(iv["expr"] for iv in invs)
    assert any(iv["n_terms"] == 3 and "A" in iv["expr"] and "B" in iv["expr"]
               and "C" in iv["expr"] and "**" not in iv["expr"]
               for iv in invs), joined


def test_equations_and_shared_constants():
    data = _chain()
    cert = discover_system(data, per_target_s=30)
    assert "dA_dt" in cert.equations or "A" in cert.equations
    assert len(cert.equations) >= 2
    assert cert.alpha_log10_total is not None


def test_invariant_certified_constant_only():
    rng = np.random.default_rng(0)
    data = {"a": rng.uniform(1, 2, 200), "b": rng.uniform(1, 2, 200)}
    invs = discover_invariants(data)
    assert invs == []          # independent random columns conserve nothing

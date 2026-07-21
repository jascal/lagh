"""Lévy-exponent oracle: recover ψ(u) from the empirical characteristic function.

For a SYMMETRIC infinitely divisible law, ψ(u) is real and even, and the real
observable L(u) = log|φ̂(u)| = ψ(u) maps onto lagh's real-function engine. The novel
piece is the STATISTICAL epsilon: φ̂ carries O(1/√n) sampling error, so certification
is at a significance level, via a bootstrap se fed to the engine's `se` argument.
"""

from __future__ import annotations

import numpy as np

# --- symmetric infinitely-divisible samplers (the generators) ---------------

def sample_gaussian(n, sigma, rng):
    return rng.normal(0, sigma, n)


def sample_sym_compound_poisson(n, rate, jump, rng):
    """N ~ Poisson(rate) jumps of +-jump each (symmetric)."""
    counts = rng.poisson(rate, n)
    out = np.zeros(n)
    for i, c in enumerate(counts):
        if c:
            out[i] = jump * (2 * rng.integers(0, 2, c) - 1).sum()
    return out


def sample_variance_gamma(n, b, rng):
    """Symmetric variance-gamma-like: Gaussian with Gamma-mixed variance.
    L(u) = -c*log(1 + u^2/b^2) shape."""
    g = rng.gamma(1.0, 1.0, n)
    return rng.normal(0, 1, n) * np.sqrt(g) / b


def sample_sym_stable(n, alpha, c, rng):
    """Symmetric alpha-stable via CMS. L(u) = -c^alpha * |u|^alpha."""
    U = rng.uniform(-np.pi / 2, np.pi / 2, n)
    W = rng.exponential(1.0, n)
    if abs(alpha - 1.0) < 1e-9:
        X = np.tan(U)
    else:
        X = (np.sin(alpha * U) / np.cos(U) ** (1 / alpha)
             * (np.cos(U - alpha * U) / W) ** ((1 - alpha) / alpha))
    return c * X


GENERATORS = {
    "L1_gaussian": lambda n, rng: sample_gaussian(n, 1.3, rng),
    "L2_compound_poisson": lambda n, rng: sample_sym_compound_poisson(n, 2.0, 1.0, rng),
    "L3_variance_gamma": lambda n, rng: sample_variance_gamma(n, 1.5, rng),
    "L4_stable_rational": lambda n, rng: sample_sym_stable(n, 1.5, 1.0, rng),
    "L5_stable_irrational": lambda n, rng: sample_sym_stable(n, np.sqrt(2), 1.0, rng),
    "L6_mixed": lambda n, rng: (sample_gaussian(n, 1.0, rng)
                                + sample_sym_compound_poisson(n, 1.5, 1.2, rng)),
}

# true L(u) = log|phi(u)| for scoring (computed, never fed to recovery)
TRUE_L = {
    "L1_gaussian": lambda u: -0.5 * 1.3**2 * u**2,
    "L2_compound_poisson": lambda u: 2.0 * (np.cos(1.0 * u) - 1.0),
    "L3_variance_gamma": lambda u: -np.log(1.0 + u**2 / (2 * 1.5**2)),  # shape ref
    "L4_stable_rational": lambda u: -np.abs(u) ** 1.5,
    "L5_stable_irrational": lambda u: -np.abs(u) ** np.sqrt(2),
    "L6_mixed": lambda u: -0.5 * 1.0**2 * u**2 + 1.5 * (np.cos(1.2 * u) - 1.0),
}


def empirical_L(X, u):
    """L̂(u) = log|φ̂(u)|, φ̂(u) = mean(exp(i u X))."""
    phi = np.exp(1j * np.outer(u, X)).mean(axis=1)
    return np.log(np.abs(phi))


def bootstrap_se(X, u, B=40, seed=0):
    """Statistical se of L̂(u) by resampling the data -- the significance epsilon."""
    rng = np.random.default_rng(seed)
    n = len(X)
    vals = np.empty((B, len(u)))
    for b in range(B):
        vals[b] = empirical_L(X[rng.integers(0, n, n)], u)
    return vals.std(axis=0)


def reliable_u_max(X, U, n, floor_mult=8.0, n_probe=200):
    """Largest u where |φ̂(u)| stays above the sampling floor ~1/√n. Beyond it the
    empirical CF is pure noise (|φ̂| of an average of unit-modulus numbers cannot go
    below O(1/√n)), and log|φ̂| saturates/biases. This is adaptive ranging in the CF
    domain -- the u-analog of the signal-floor box contraction."""
    floor = floor_mult / np.sqrt(n)
    ug = np.linspace(U / n_probe, U, n_probe)
    phi = np.abs(np.exp(1j * np.outer(ug, X)).mean(axis=1))
    below = np.where(phi < floor)[0]
    return float(ug[below[0]]) if len(below) else float(U)


def adaptive_cf_grid(X, n, *, U_wide=16.0, floor_mult=8.0, n_probe=600, n_keep=70):
    """Sample u over the UNION of reliable regions {u : |φ̂(u)| > floor}, not just the
    first interval. For an oscillatory CF (compound Poisson) |φ̂| dips at u=π then
    RECOVERS at u=2π; capturing those recovery bumps reveals the oscillation and
    distinguishes it from a monotonic exponent. For a monotonic CF this reduces to the
    single reliable interval. Adaptive ranging that unions instead of truncating."""
    floor = floor_mult / np.sqrt(n)
    ug = np.linspace(U_wide / n_probe, U_wide, n_probe)
    phi = np.abs(np.exp(1j * np.outer(ug, X)).mean(axis=1))
    reliable = ug[phi > floor]
    if len(reliable) == 0:
        reliable = ug[:5]
    if len(reliable) > n_keep:
        reliable = reliable[np.linspace(0, len(reliable) - 1, n_keep).astype(int)]
    return reliable


def make_cf_dataset(name, *, n=20000, U=6.0, n_grid=60, seed=0):
    """Sample the process, build (u, L̂(u), se(u)) over the RELIABLE u-range (where the
    empirical CF is above the sampling floor). u is the 1-D input; L̂(u) the output;
    se the statistical epsilon. Returns (u, L, se, u_max) for disclosure."""
    rng = np.random.default_rng(seed)
    X = GENERATORS[name](n, rng)
    u = adaptive_cf_grid(X, n)                       # union of reliable regions
    L = empirical_L(X, u)
    se = bootstrap_se(X, u, seed=seed + 1)
    return u[:, None], L, se, float(u.max())

# Exoplanet Archive C2 — conditioning on stellar properties

**Registered 2026-07-27, predictions before the run.** Governing proposal:
`docs/PROPOSAL_EXOPLANET.md`. Same frozen snapshot as C1
(`c1_d2e83168c0d1.csv` — the fetch is idempotent); conditioning variables
st_met, st_mass, st_teff were in the frozen query. All astrophysical:
conjecture track, `recover` must abstain everywhere.

## Registered predictions

- **P1 (giant-planet metallicity preference).** Hosts of giants
  (M > 100 M⊕, true-mass window) are metal-rich relative to hosts of small
  planets (R < 4 R⊕, same quality window): mean offset Δ[Fe/H] in
  **[0.05, 0.25] dex** (Fischer–Valenti/Johnson core–accretion signature),
  bootstrap 95% CI excluding 0. Caveat registered: the sample mixes
  detection methods; the offset is reported, not decomposed.
- **P2 (the valley moves with stellar mass).** On the C1 Fulton-like sample
  split into stellar-mass bins: d log R_v / d log M★ in **[0.0, 0.5]**
  (Berger et al. ≈ +0.26). Estimator: per-bin KDE minimum on a widened grid
  [1.2, 2.6] R⊕ (the C1 grid-edge caveat addressed); slope over ≥ 3 bins.
- **P3 (metallicity does NOT reshape the volatile M–R slope at current
  precision).** Volatile-window (1.8–4 R⊕) log-log slope computed for
  above/below-median [Fe/H] halves: |Δslope| < **0.8**. A larger split
  would be a band miss flagging real conditional structure.
- **P4.** `recover` structurally abstains on every conditioned window; zero
  confident-wrong; no certificate anywhere in C2.

## Results (2026-07-27, `experiments/results/exoplanet_c2.json`; frozen
snapshot `c1_d2e83168c0d1.csv`)

- **P1 MET:** giant hosts are metal-rich by **Δ[Fe/H] = 0.124 dex**,
  bootstrap 95% CI [0.102, 0.147] — excludes 0, inside the registered
  [0.05, 0.25]: the core-accretion signature, with the registered
  detection-method caveat standing.
- **P2 ESTIMATOR-INVALID — reported as a miss, not a discovery.** The
  measured slope (−0.44) is outside the band with inverted sign, but it is
  not credible: the two highest-stellar-mass bins returned identical valley
  locations of 1.4034 R⊕ — the first grid point inside the [1.4, 2.4]
  search window. The registered edge check tested the KDE grid boundary
  (1.2 R⊕), not the search-window boundary, so the pinning went unflagged:
  in those bins the estimator found no interior minimum (sparse counts near
  1.4–1.5 R⊕), and the fitted slope is a window artifact. The honest
  verdict: this estimator, on this sample, does not measure the
  valley–stellar-mass slope; the prediction is neither confirmed nor
  refuted. (A per-bin minimum-interiority requirement is the fix for any
  future valley conditioning.)
- **P3 MET:** volatile-window slopes 1.14 (low-Z) vs 1.30 (high-Z),
  |Δ| = 0.15 < 0.8 — no metallicity restructuring of the mass–radius
  relation at current precision; `recover` structurally abstained on both
  halves.
- **P4 MET:** zero certificates anywhere in C2.

**The sentence this stage earns:** the one clean stellar-conditioning
signal in the sample — giant planets prefer metal-rich hosts, at 0.12 dex —
came through inside its registered band, and the stage's failed estimator
was reported as a failed estimator rather than dressed up as an
inverted-sign discovery.

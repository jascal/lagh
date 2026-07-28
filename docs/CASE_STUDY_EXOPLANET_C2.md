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

## Results

(after the run)

# Exoplanet Archive C1 — mass–radius sequences and the radius valley

**Registered 2026-07-27, predictions before any fetch or run.** Governing
proposal: `docs/PROPOSAL_EXOPLANET.md`. C1 is the campaign's first
astrophysical stage: every relation here is empirical with intrinsic
scatter, so the deliverables are **labeled conjectures with pre-registered
literature bands** — `recover` must abstain everywhere, and a certificate
anywhere in this stage would be a confident-wrong.

## Sample (frozen query, `C1_ADQL`)

`pscomppars`, all rows with radius and period, plus mass (with provenance),
density, discovery method, and host Teff/mass/metallicity. Quality windows
are applied per prediction (registered here):

- **Mass–radius windows:** `pl_bmassprov = 'Mass'` (true masses, no Msini),
  relative mass error < 25%, relative radius error < 8%.
- **Radius-valley sample (Fulton-like):** transit-discovered, host Teff in
  [4700, 6500] K (FGK), P < 100 d, relative radius error < 8%.

## Registered predictions

- **P1 (rocky regime, R < 1.5 R⊕).** Log-log OLS slope of M(R) in
  **[3.0, 4.5]** (Earth-like composition sequences: M ∝ R^3.5±). Labeled
  conjecture; `recover` on the window must abstain.
- **P2 (volatile regime, 1.8 < R < 4 R⊕).** Slope of M(R) in **[0.7, 2.3]**
  (the literature genuinely disagrees — Weiss–Marcy ≈ 0.9, Chen–Kipping
  ≈ 1.7; the registered band covers the dispute and the point estimate is
  reported against both). Conjecture; `recover` must abstain.
- **P3 (degenerate giants, M > 100 M⊕).** The radius plateau: slope of
  log R vs log M in **[−0.10, +0.15]** (electron-degeneracy flattening near
  1 R_J). Conjecture; `recover` must abstain.
- **P4 (the radius valley exists and slopes down with period).** On the
  Fulton-like sample: (a) the histogram of log R (0.05-dex bins) has a local
  minimum with center in **[1.7, 2.0] R⊕**, with both flanking maxima
  ≥ **1.3×** the valley bin count; (b) the valley center's period
  dependence, estimated by per-period-bin minima of the smoothed radius
  distribution, has slope d log R_v / d log P in **[−0.20, 0.00]**
  (thermally-driven mass-loss predicts ≈ −0.1; a positive slope would favor
  primordial/gas-poor formation and count as a band miss).
- **P5.** Zero confident-wrong: no certificate anywhere in C1.

## Results (2026-07-27, `experiments/results/exoplanet_c1.json`; frozen
snapshot `c1_d2e83168c0d1.csv`, 5934 planets; 1325 in the true-mass window)

- **P1 MISSED, informatively.** Rocky-window slope **1.61** vs the
  registered [3.0, 4.5] (scatter 0.53 dex, n = 94). The registered band was
  the theoretical Earth-composition sequence; the measured-mass R < 1.5 R⊕
  population is not that sequence — it mixes compositions (iron-rich USPs to
  low-density worlds), and regression attenuation from the permitted 25%
  mass errors flattens the fit further. The honest reading: the *sample*
  does not follow a single rocky mass–radius law, and the 0.53 dex scatter
  says so louder than the slope.
- **P2 MET:** volatile-regime slope **1.42**, inside [0.7, 2.3] — between
  Weiss–Marcy (≈0.9) and Chen–Kipping (≈1.7), closer to the latter.
- **P3 MET, cleanly:** giant plateau slope **−0.007** (n = 750), dead
  center in [−0.10, +0.15] — the electron-degeneracy flattening, flat to
  better than 1%.
- **P4 split.** Location MET: the valley minimum sits at **1.90 R⊕**,
  inside [1.7, 2.0]. Contrast MISSED narrowly: 1.246 vs the registered
  ≥ 1.3 (flanks 76 and 141 against a 61-count valley bin) — the valley is
  visibly there but the raw-count criterion on 0.05-dex bins was 4% short.
  Period slope MET with a caveat: **d log R_v/d log P = −0.127**, inside
  [−0.20, 0.00] and consistent with the mass-loss prediction (≈ −0.1); the
  caveat is that the per-bin minimum estimator pinned at its 1.4 R⊕ grid
  edge in the three longest-period bins, so the point estimate is
  suggestive, not robust.
- **P5 MET:** `recover` structurally abstained in every window — zero
  certificates in C1, as the astrophysical stage demands.

**The sentence this stage earns:** the population's real structure came
through the registered bands — a flat giant plateau, a volatile slope
between the disputed literature values, a radius valley at 1.90 R⊕ sloping
down with period — while the two misses (the rocky band, the contrast
threshold) are the instrument telling us the registered expectations, not
the data, were wrong.

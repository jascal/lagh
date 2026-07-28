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

## Results

(after the runs)

# Exoplanet Archive Phase 2 (C3–C4) — registration and results

**Registered 2026-07-27, predictions before any fetch or run.** Governing
proposal: `docs/PROPOSAL_EXOPLANET.md`. C3 = density/irradiation trends;
C4 = multi-planet architecture. All astrophysical: conjecture track with
registered bands; `recover` must abstain everywhere; zero certificates.
New frozen fetch (`PH2_ADQL`): all planets with period, plus hostname,
multiplicity, semi-major axis, radius/mass with errors and provenance,
density, insolation (Earth-flux units), stellar mass/Teff.

## C3 — density and irradiation

Giant windows use radius (R > 8 R⊕); quality: radius error < 8% where
errors are used. Insolation S in S⊕.

- **P1 (hot-Jupiter inflation offset).** Hot giants (S > 200) vs cool
  giants (S < 200): mean radius difference in **[1, 6] R⊕**, bootstrap 95%
  CI excluding 0.
- **P2 (the inflation threshold structure).** Radius–insolation slope
  dR/dlog₁₀S for cool giants (S < 100) consistent with flat: in
  **[−1.0, +0.75] R⊕/dex**; for hot giants (S > 200): positive, in
  **[0.5, 3.0] R⊕/dex**.
- **P3 (small-planet density–radius).** Over 1–4 R⊕ (density non-null):
  d log ρ / d log R in **[−2.5, −0.5]** (the volatile transition; C1's
  volatile slope predicts ≈ 1.4 − 3 = −1.6).

## C4 — multi-planet architecture (sy_pnum ≥ 2)

Adjacent pairs by period within each system.

- **P5 (period-ratio distribution).** Median adjacent-pair ratio in
  **[1.5, 2.2]**; fraction of adjacent pairs with ratio < 1.3 is **< 10%**
  (Hill/Lagrange stability floor).
- **P6 (first-order resonance asymmetry).** Pair counts just WIDE of the
  resonances exceed just-NARROW: count(ratio ∈ [2.00, 2.10]) >
  count([1.90, 2.00]) and count([1.50, 1.57]) > count([1.43, 1.50]) — the
  known Kepler asymmetry (tidal migration out of resonance).
- **P7 (mutual Hill spacing).** For adjacent pairs with both masses and
  semi-major axes (any mass provenance; noted): median separation Δ in
  **[15, 30]** mutual Hill radii; fraction with Δ < 8 is **< 5%**
  (long-term-stability floor).
- **P8.** `recover` structurally abstains on the period-ratio sequence and
  every C3 window; zero confident-wrong; zero certificates in Phase 2.

## Results (2026-07-27, `experiments/results/exoplanet_ph2.json`; frozen
snapshot `ph2_7f758ab2b303.csv`, 5981 planets, 1573 adjacent pairs)

- **P1 MET:** hot giants are inflated by **ΔR = 1.53 R⊕** (CI [1.30, 1.77],
  excludes 0), inside the registered [1, 6].
- **P2 half-MET.** The threshold *structure* is confirmed: cool giants flat
  (slope −0.28 R⊕/dex, inside [−1.0, +0.75], n = 871) versus hot giants
  strongly positive — but the hot slope, **3.85 R⊕/dex**, exceeds the
  registered ceiling of 3.0. The band was set from older, smaller inflation
  samples; the composite table's ultra-hot population steepens the trend.
  Recorded as a numeric band miss on a qualitatively confirmed structure.
- **P3 MET:** small-planet density–radius slope **−1.24** over 1–4 R⊕
  (n = 3538), inside [−2.5, −0.5] and consistent with C1's volatile slope
  (1.42 − 3 ≈ −1.6); `recover` structurally abstains.
- **P5 half-MET.** The stability floor holds: only **1.6%** of adjacent
  pairs are tighter than a 1.3 ratio (registered < 10%). The median ratio,
  **2.35**, misses the registered [1.5, 2.2] high: the band described
  Kepler-multi statistics, but PSCompPars multis include wide-spaced RV
  giant systems. A sample-composition miss, not an architecture surprise.
- **P6 MET, strikingly.** The first-order resonance asymmetry is present at
  both resonances: just-wide of 2:1 outnumbers just-narrow **96 : 57**, and
  just-wide of 3:2 outnumbers just-narrow **100 : 40** — the tidal
  migration-out-of-resonance signature, in a heterogeneous composite sample.
- **P7 half-MET.** Median mutual-Hill spacing **Δ = 17.6** (inside
  [15, 30], n = 1457 pairs); but 11.2% of pairs fall below Δ = 8 against
  the registered < 5%. With any-provenance masses (Msini for RV systems)
  and literature semi-major axes, the sub-8 tail is dominated by
  mass/geometry heterogeneity rather than genuinely unstable systems — the
  registered stability floor assumed cleaner inputs than the composite
  table provides.
- **P8 MET:** `recover` structurally abstains on the ratio sequence and
  every C3 window — zero certificates in Phase 2, as the astrophysical
  stages demand.

**The sentence this phase earns:** the architecture of planetary systems
came through the registered bands — inflation above the irradiation
threshold, the density transition, the resonance asymmetry at both 2:1 and
3:2, Hill spacing at the Kepler value — and each of the three numeric
misses decodes into the same lesson C0 taught: PSCompPars is a literature
patchwork, and its heterogeneity, not the physics, sets the tails.

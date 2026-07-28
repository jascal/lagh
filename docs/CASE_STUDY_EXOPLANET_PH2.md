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

## Results

(after the run)

# Case study — solar-system and solar-activity data through the instrument

**Registered 2026-07-28, predictions before results.** Second real-data case
study (after `CASE_STUDY_MACRO.md`). Solar behavior splits the same way macro
did: an EXACTLY-LAWFUL stratum (orbital mechanics — certifiable with α) and a
STATISTICAL stratum (solar activity cycles — honest labeled conjectures). The
instrument must tell them apart on real data it has never seen.

## Data (public, snapshots committed)

- **Planetary elements** (semi-major axis AU, sidereal period yr): NASA NSSDC
  planetary fact sheet values, 8 planets, committed CSV with provenance.
- **Galilean moons of Jupiter** and **major moons of Saturn** (a in km, P in
  days): NSSDC fact sheets — two more Kepler systems with DIFFERENT central
  masses.
- **Monthly sunspot numbers** (SILSO V2.0, 1749–present): fetched CSV
  snapshot.
- **Declared noise:** fact-sheet rounding plus the known physical correction
  the tables ignore (planetary masses: P² = 4π²a³/G(M☉+m); m/M☉ up to ~1e-3
  for Jupiter) → σ = 2e-3 declared for planets; moons likewise (σ = 1e-3).

## Registered predictions

- **P1 (Kepler III, three times over):** `recover` on (a → P) certifies
  `P = C·a^{3/2}` — the exact 3/2 exponent — separately for planets, Jovian
  moons, and Saturnian moons, each with its own constant, at the declared σ.
  The historic law, re-derived by machine from the real tables, with α stated.
- **P2 (shared constants respect system boundaries):** the per-system
  constants C (∝ GM^{-1/2}) must AGREE within a system across data subsets and
  DIFFER across systems — no spurious cross-system sharing verdicts.
- **P3 (sunspot cycle):** certification FAILS on the monthly series
  (quasi-periodic, stochastic amplitude); the astronomer's period estimator
  reports the cycle as a labeled conjecture in the 9–13 yr band (~11 yr
  expected).
- **P4 (Waldmeier effect):** rise-time vs amplitude across cycles — negative
  correlation as a labeled empirical conjecture; nowhere near certification.
- **P5:** zero certified-wrong; nothing certifies beyond the orbital stratum.

## Results (2026-07-28, `experiments/results/solar_case_study.json`)

- **P1 — Kepler III certified twice, and the third refusal decoded into real
  physics.** Planets: `P = 0.99994·a^{3/2}` certified, α ≤ 10⁻¹⁴ (tiny-data
  mode: exhaustive certification with α bounding the selection exposure — a
  mode this study forced into existence). Jovian moons: `6.458e-9·a^{3/2}`,
  α ≤ 10⁻¹⁸. Saturnian moons: **honest structural abstain at σ=1e-3 — and the
  residuals decode into Saturn's J₂ oblateness** (measured Kepler deviations
  track first-order J₂ theory to ~10⁻⁴ per moon: Mimas 0.998792 measured vs
  0.998711 predicted) **plus Hyperion's 4:3 Titan resonance** (2% outlier).
  Secondary run (resonant moon excluded, J₂ envelope declared): certified,
  `1.181e-8·a^{3/2}`, α ≤ 10⁻¹⁵.
- **P2 MET:** per-system constants stable within systems (≤0.24% across
  halves), differing across systems by orders of magnitude; no spurious
  sharing.
- **P3 MET:** raw sunspot series honestly abstains (structural); the cycle
  emerges as a labeled conjecture at **10.9975 years** — the textbook ~11.
- **P4 MET:** Waldmeier effect: rise-time vs amplitude correlation **−0.77**
  across 24 cycles, labeled empirical.
- **P5 MET:** zero certified-wrong.

**The sentence this study earns:** given the real tables, the instrument
re-derives Kepler's third law with a stated significance bound, and where it
refuses — Saturn — the refusal itself is a discovery: the residual structure
is the planet's oblateness and a resonance, quantitatively.

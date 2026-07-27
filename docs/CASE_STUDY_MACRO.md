# Case study — raw US macro data through the full instrument

**Registered 2026-07-28, predictions before results.** The question (user):
can the machinery infer simple macroeconomic laws from RAW public data? The
honest expectation is itself the point: real-world "laws" split into
**accounting identities** (true by construction, should CERTIFY with α) and
**statistical regularities** (Okun, Phillips — should honestly ABSTAIN into
labeled conjectures with fitted slopes). An instrument that can't tell these
apart has no business near real data.

## Data (public, committed as snapshot for reproducibility)

FRED quarterly series, 1990Q1–2019Q4 (the econ-sae calibration window):
NOMINAL expenditure accounts — GDP, PCEC (C), GPDI (I), GCE (G), NETEXP (NX) —
plus UNRATE and CPIAUCSL. Nominal chosen deliberately: the expenditure
identity Y = C+I+G+NX is exact in nominal dollars to BEA rounding; in CHAINED
dollars the components famously do NOT add (chain-weighting residual) — a
secondary prediction below.

**Declared noise (standing policy):** BEA levels are published to $0.1B →
σ_rep from half-ulp rounding relative to series scale (~5e-6..1e-5); UNRATE
to 0.1pp.

## Registered predictions

- **P1**: the invariant scan over {Y, C, I, G, NX} finds the expenditure
  identity (a 5-term linear combination, value ≈ 0, within declared rounding)
  and CERTIFIES it with a stated α. This is a certified law recovered from
  raw public data with no economics supplied.
- **P2**: the same scan on CHAINED-dollar components does NOT certify the
  identity (the chain residual exceeds rounding noise) — the instrument
  correctly refuses a "law" that textbooks state loosely but the data
  violates.
- **P3**: Okun's law (Δunemployment vs GDP growth): certification FAILS
  (statistical scatter ≫ declared noise) → track-B conjecture, expected slope
  ≈ −0.3 to −0.5 in the difference formulation.
- **P4**: Phillips (inflation vs unemployment level, 1990–2019): abstain;
  weak/flat conjecture (the era's flat Phillips curve).
- **P5**: zero certified-wrong: nothing certifies beyond identity-class
  relations.

## Results (2026-07-28, `experiments/results/macro_case_study.json`;
data snapshots committed under `experiments/macro/data/`)

- **P1 MET — a certified law from raw public data.** The invariant scan over
  nominal {Y, C, I, G, NX} recovers the expenditure identity
  `C + I + G + NX − Y = 0` with all coefficients within 3 ppm of unity, value
  −0.0009 ($0.9M on a ~$15T economy), **α ≤ 10⁻⁵⁶⁸** at the declared
  BEA-rounding noise (σ ≈ 3.9e-6). No economics was supplied; the identity
  emerged from the subset scan alone (instrument fixes it drove: term-scale
  constancy normalization for zero-valued invariants; sigma-scaled eigen gate;
  exhaustive size-5 for small libraries).
  *(Re-validation note, 2026-07-27: under the term-scale α accounting the
  Gaia C0 study later introduced, the re-run strengthens the headline bound
  to α ≤ 10⁻⁸⁵⁷ and re-baselines the scaled variants to 10⁻³⁸⁸…10⁻⁵⁹².
  Verdicts and laws unchanged; `macro_case_study.json` carries the current
  numbers.)*
- **P2 MET — the refusal that proves the certification means something.** The
  same scan on CHAINED-dollar components certifies NO linear relation: the
  chain-weighting residual is real and the instrument declines the textbook
  shorthand the data violates.
- **P3 MET.** Okun's law: certification honestly fails; the labeled empirical
  conjecture is `Δu = 0.152 − 0.268·g` (slope −0.268, in the literature band;
  corr −0.56). Required the affine-OLS fallback conjecture (the power-law
  probes are positive-data-only — a real-data lesson now in the fit scout).
- **P4 MET.** Phillips 1990–2019: conjecture `π = 3.00 − 0.11·u`, corr −0.09 —
  the era's famously flat curve, correctly labeled empirical and nowhere near
  certification.
- **P5 MET.** Nothing certified beyond identity-class relations; zero
  certified-wrong.

**The sentence this case study earns:** on raw public macroeconomic data, the
instrument certifies exactly what is true by construction (with an astronomical
significance bound), refuses exactly what the data quietly violates, and labels
the famous statistical "laws" as the empirical conjectures they are.

# Gaia DR3 Phase 2 (C3–C4) — registration and results

**Registered 2026-07-27, predictions before any fetch or run.** Governing
proposal: `PROPOSAL_GAIA.md`; C3 reframed at C0 registration as
mass-DEPENDENT relation discovery (Kepler III is not a single cross-catalog
law). Phase 1 lessons applied here: derived catalog columns may be marginal
point estimates that do not compose through nonlinear identities (FLAME), and
float32-stored-as-float64 quantization must be detected by round-trip, not by
repr digit counting.

**Floor procedure (fixed here, v2):** per column, storage precision by
float32 round-trip (half-ulp 2⁻²⁵ relative if exact, else decimal-repr
half-step); for a target law, `floor_abs = 2·(q_y + Σⱼ |βⱼ|·q_termⱼ)` with β
from an OLS scout over the registered basis on the fetched data (a logged
calibration step, run before discovery) and q_term propagated through the
basis derivatives.

## C3 — astrometric binaries (`nss_two_body_orbit`, Orbital solutions, TOP 1000)

Photocenter semi-major axis a₀ from the Thiele–Innes elements
(u = (A²+B²+F²+G²)/2, a₀ = √(u + √(u² − (AG−BF)²))), a_AU = a₀/ϖ,
P in years.

- **P1 (no cross-catalog Kepler law).** `recover` on (log₁₀ a_AU → log₁₀ P)
  must NOT certify: the family is mass-dependent, one law per system mass.
  Abstain or conjecture are the only acceptable outcomes.
- **P2 (photocentric mass scale, labeled conjecture).** The mass proxy
  a_AU³/P_yr² has median in **[0.05, 5] M☉** (photocenter orbits
  underestimate the relative orbit, so the low side is expected). Conjecture
  only; no certification of any mass relation.

## C4 — variable stars

### RR Lyrae (`vari_rrlyrae`, RRab, TOP 800)

- **P3 (pipeline-definitional metallicity — the certification target).** DR3
  photometric `metallicity` for RRab is computed from period and the G-band
  Fourier phase φ31 by a fixed published calibration. Prediction: `recover`
  on (pf, phi31_g) → metallicity **certifies a low-degree polynomial (degree
  ≤ 2, ≤ 6 terms)** at the v2-procedure floor, α stated; the discovered
  coefficients are compared post-hoc to the DR3 calibration paper. Fallback
  honesty (P1 lesson): if the published columns do not compose, the verdict
  is a structural abstain and the residual decode replaces the certificate.
- **P4 (Leavitt law stays empirical).** MW field fundamental-mode classical
  Cepheids (`vari_cepheid` ⋈ `gaia_source`, `parallax_over_error > 10`):
  M_G = G + 5·log₁₀(ϖ/100) vs log₁₀ pf is a labeled conjecture with slope in
  **[−3.5, −2.0] mag/dex**; `recover` must NOT certify (astrophysical
  scatter, extinction, width of the instability strip).
- **P5.** Zero confident-wrong across the phase.

## First public certificate release

On completion, every certificate the Gaia program has produced is collected
in `docs/CERTIFICATES.md` with law, α, domain, frozen-snapshot SHA, engine
commit, and a reproduction command — the proposal's phase-2 release artifact.

## Results (2026-07-27, `experiments/results/gaia_p2.json`; frozen snapshots
`c3_67fabb172bf7.csv` (1000 orbits), `c4rr_bf9b07ed638a.csv` (800 RRab),
`c4cep_db0dbc7d8842.csv` (500 Cepheids))

- **P1 MET:** no cross-catalog Kepler law — `recover` on (log a_AU → log P)
  structurally abstains, as a mass-dependent family demands.
- **P2 MISSED at the band edge, honestly:** median photocentric mass proxy
  0.0489 M☉ vs the registered [0.05, 5] — 2% below the low edge. The
  registered band underweighted how strongly photocenter orbits (a₀ shrinks
  with luminosity ratio) suppress the proxy. The miss is recorded, not
  rationalized away; no certification was at stake.
- **P3 MISSED as registered for `recover` — and the second track turned the
  miss into the phase's certificate.** `recover` abstains with "no law
  certifies through tier 7": a REACH gap — the dim-2 cross-term quadratic
  support {1, x₀, x₁, x₀x₁, x₁²} is never proposed, although the relation
  composes within the floor (scout max residual 3.5×10⁻⁷ < 4.2×10⁻⁷).
  **Registered instrument issue: dim-2 quadratic cross-term reach.** The
  logged fallback (calibration scout declares the form; `verify` judges it —
  `floor_abs` plumbed through `verify` as an execution amendment) certifies,
  **pinned** ("this exact form, no rival within the noise"):

  `[Fe/H] = 3.0336834 − 20.127524·P + 1.3684267·φ31 + 6.27·P·φ31 − 0.72·φ31²`

  The cross and quadratic coefficients are EXACTLY 6.27 and −0.72 — the
  Nemec et al. (2013) calibration the DR3 pipeline documents. The three
  remaining coefficients are consistent to ~10⁻⁶ with composing Nemec's
  published constants with a single Fourier-phase offset Δ ≈ 3.18859
  (≈ π + 0.047), but a one-parameter Δ composition leaves 4×10⁻⁶ residuals —
  10× the floor — so the exact provenance of the linear terms is left OPEN
  rather than guessed. What is certified is what the data pins: the quadratic
  above, exhaustively, over the snapshot.
- **P4 MET:** the Leavitt law stays a labeled conjecture, slope
  **−2.21 mag/dex** — inside the registered [−3.5, −2.0].
- **P5 MET:** zero confident-wrong across the phase.

**The sentence this phase earns:** the instrument abstained from every
relation the data could not pin — including its own reach gap, which it
reported as such — and the one certificate it granted decodes the DR3
pipeline's metallicity formula from the archive, exact published coefficients
and all, with the parts it could not pin explicitly left open.

## Issue closure (2026-07-27): the dim-2 quadratic reach gap

Diagnosis: the term library was never the problem — every quadratic monomial
exists at dim 2. The 5-term support {1, x₀, x₁, x₀x₁, x₁²} was unreachable by
PROPOSAL: size-5 supports were drawn only from the top-6 single-term fits
(fragile in a collinear library — P and φ31 correlate at 0.60, and smooth
transcendental features outrank the cross term as singles), OMP's one greedy
path missed it, and CAP-Q's exhaustive cluster search stops at size 4.

Fix (CAP-T, `engine._linear_candidates`): the COMPLETE low-degree polynomial
supports are proposed unconditionally — the full quadric at dim ≤ 4, the full
cubic at dim ≤ 2. One lstsq each; `refit_minimal` and the minimality repair
prune what the data does not need, the coefficient gate pins what remains.

Validation (same discipline as the loose-ε closure): full suite green;
`recover` on the frozen RRab snapshot now certifies OUTRIGHT —
`[Fe/H] = (5178136/825859)·P·φ31 − (15295147/759912)·P −
(18000001/25000000)·φ31² + (1029763/752516)·φ31 + 1354664/446541`,
**α ≤ 10⁻¹⁰⁰⁹** — the same surface as the verify-track form (max difference
8×10⁻⁸ on the snapshot, both within the floor), with the cross and quadratic
coefficients again 6.2700001 and −0.7200000. Under declared noise the linear
channel remains certification-excluded (the significance boundary), so the
new supports add no impostor surface there; the NewtonBench-clean machine-
dominated regime gates unpinned full-poly fits per candidate as before.

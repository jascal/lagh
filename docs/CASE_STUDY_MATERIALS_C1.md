# Materials Project C1 — elastic scaling: certificates and banded conjectures

**Registered 2026-07-28, predictions before any fetch or run.** Governing
proposal: `docs/PROPOSAL_MATERIALS.md`. C1 mixes the two strata explicitly:
the VRH averaging convention is pipeline-definitional (certifiable with
exact rational coefficients), while every physical scaling relation belongs
on the banded conjecture track.

## Sample (frozen query)

`summary` endpoint, stable materials, fields: material_id, composition,
nsites, volume, density, bulk_modulus, shear_modulus; first 2000 by
material_id, filtered client-side to rows with elastic moduli (the API
serves the VRH/Voigt/Reuss dictionaries where elasticity was computed).
Rows with non-positive or non-finite moduli are dropped (logged count).

## Registered predictions

- **P1 (the VRH identities — definitional).** The pipeline computes the
  Voigt–Reuss–Hill moduli as the arithmetic mean of the Voigt and Reuss
  bounds. `recover` on (K_voigt, K_reuss) → K_VRH certifies
  `x₀/2 + x₁/2` — exact coefficients **1/2, 1/2** — at the v2/σ_rep floor;
  same for the shear moduli (G). α stated for both.
- **P2 (Voigt–Reuss ordering).** K_reuss ≤ K_vrh ≤ K_voigt holds for every
  row (bound-property sanity; violations counted and, if any, decoded —
  they would indicate ingestion defects, not physics).
- **P3 (Birch/Anderson-type scaling within a family, banded conjecture).**
  Within the binary oxide family (materials whose composition is {metal, O}):
  log K_vrh vs log(V/nsites) has OLS slope in **[−2.0, −0.8]**
  (Anderson–Nafe-type volume scaling). `recover` on the window must
  structurally abstain (real scatter).
- **P4 (Pugh ratio, banded conjecture).** The distribution of G_vrh/K_vrh
  has median in **[0.35, 0.70]** (metals cluster near 0.4, ceramics higher;
  the Pugh ductile/brittle boundary sits at 0.57). Conjecture only.
- **P5 (cross-family scaling must NOT certify).** `recover` on the FULL
  mixed sample (log V/atom → log K) structurally abstains — family mixtures
  are not one law; the scout slope is reported as a labeled conjecture.
- **P6.** Zero confident-wrong; certificates clear α ≤ 10⁻⁶.

**Execution amendments (logged, 2026-07-28, after the first run):** (1) the
API caps `_limit` at 1000/page — adapter paginates (no data change). (2) The
moduli are stored with FIXED 3-decimal rounding in GPa — ABSOLUTE precision,
the Gaia absolute-floor lesson inverted: σ_rep under-covers small moduli;
`floor_abs = 3×10⁻³ GPa` declared from the composed half-steps. (3) The
bound-ordering census runs at rounding tolerance, violations decoded by ID.

## Results (2026-07-28, `experiments/results/materials_c1.json`; frozen
snapshot `c1_ecbea094bbdd.json`, 2000 fetched, 1368 clean elastic rows)

- **P1 split by track, both halves sound.** K: **certified outright,
  `K_VRH = K_voigt/2 + K_reuss/2` — exact 1/2 coefficients,
  α ≤ 10⁻²⁵²⁸** (the program's strongest bound to date). G: `recover`
  honestly abstains (parametric — with shear moduli mostly < 400 GPa, the
  ±10⁻⁵-relative perturbation shifts predictions below the 3×10⁻³ floor, so
  the coefficient is identified to 10⁻⁴ but not 10⁻⁵; K passes because its
  larger values break the same perturbation); the declared exact form
  `x₀/2 + x₁/2` then **verifies pinned** (no rational rival) — the RRab
  second-track precedent.
- **P2: 12 genuine violations, decoded.** Beyond rounding tolerance, 12
  bound-ordering violations across 9 materials — every one **Reuss > Voigt,
  impossible for a valid elastic tensor** (one shear modulus of 4×10⁷ GPa).
  All nine are ICSD-matched compounds with `deprecated: false` and empty
  `warnings` (checked against the live API): defective elastic-tensor fits
  served unflagged. The VRH identity still holds on every defective row —
  the pipeline averaged the defective bounds faithfully — which is exactly
  what lets the arithmetic certificate and the physical-validity census be
  independent findings. Reported upstream:
  [materialsproject/emmet#1499](https://github.com/materialsproject/emmet/issues/1499).
- **P3 MET:** binary-oxide volume scaling slope **−1.17** (n = 60), inside
  the registered [−2.0, −0.8]; `recover` structurally abstains.
- **P4 MET:** Pugh ratio median **0.522** (q10 0.34, q90 0.73), inside the
  registered [0.35, 0.70].
- **P5 MET:** the cross-family mixture does not certify (scout slope −1.38,
  labeled conjecture; structural abstain).
- **P6 MET:** zero confident-wrong.

**The sentence this stage earns:** the instrument certified the pipeline's
arithmetic at the strongest significance in the program's history, refused
the same claim where a smaller dynamic range honestly cannot pin the fifth
decimal, and its validity census pulled nine unflagged impossible entries
out of a database of computed materials — arithmetic, physics, and data
quality, cleanly separated.

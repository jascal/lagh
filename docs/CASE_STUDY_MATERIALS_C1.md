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

## Results

(after the run)

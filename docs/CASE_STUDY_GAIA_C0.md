# Gaia DR3 — Phase 0 / C0 registration (proposal: docs/PROPOSAL_GAIA.md)

**Registered 2026-07-28, predictions before any discovery run.** Phase 0 per
the proposal: the ADQL adapter + one known relation end-to-end under full
discipline. Adjustments to the proposal recorded here: C0 anchors on an
EXACTLY-LAWFUL definitional identity (catalog magnitudes are computed from
fluxes with fixed zeropoints), not on empirical physics; the binaries stage
(C3) is reframed as mass-dependent-relation discovery since Kepler III is not
a single cross-catalog law.

## C0 sample (frozen query)

`gaiadr3.gaia_source`, TOP 400 by random_index, quality cuts `ruwe < 1.4`,
`parallax_over_error > 20`, `phot_g_mean_flux_over_error > 100`, columns:
G/BP/RP magnitudes + fluxes + errors, parallax, bp_rp. Snapshot committed
(frozen artifact).

## Registered predictions

- **P1 (the definitional identity):** `recover` on (log10 flux -> mag)
  certifies `G = zp_G − 2.5·log10(flux_G)` — structure exact, zp a fitted
  constant — at declared representation noise (mags published to ~1e-4).
  Same for BP and RP with their own zeropoints. α stated.
- **P2 (color definition):** the invariant scan over {bp_rp, mag_BP, mag_RP}
  certifies `bp_rp − mag_BP + mag_RP = 0` (exact by construction).
- **P3 (real physics stays honest):** the main-sequence color–magnitude
  relation (abs G vs bp_rp on the parallax-selected sample) does NOT certify —
  labeled conjecture only (astrophysical scatter, binaries, giants).
- **P4:** zero certified-wrong.

## Amendment (logged, 2026-07-28, during the run)

The registration's "~1e-4 publication precision" assumption was WRONG: catalog
columns carry full float precision (measured true-law residuals ~1e-6). At the
loose 2e-4 floor, 14-term approximants CERTIFIED in all three bands — the
loose-epsilon parsimony exposure (same family as the sigma>0 approximant
boundary; the sigma=0 floor-dominated regime needs the minimality repair too —
**registered instrument issue**, not patched mid-study). Floor amended to the
measured 5e-6; the loose-floor artifacts are preserved in git history.

## Results (2026-07-28, `experiments/results/gaia_c0.json`; frozen snapshot
`experiments/gaia/data/c0_*.csv`, 400 sources, live TAP)

- **P1 MET (post-amendment):** all three bands certify the definitional law
  with the slope as the EXACT rational 5/2:
  `G = 25.68735 − (5/2)·log₁₀(flux)` (α ≤ 10⁻⁴⁶¹), BP (α ≤ 10⁻⁴⁵⁹),
  RP (α ≤ 10⁻⁴⁵⁸). The fitted zeropoints match the published DR3 Vega
  zeropoints (G: 25.6874).
- **P2 MET:** the color identity `bp_rp − mag_BP + mag_RP = 0` recovered with
  exact unit coefficients, **α ≤ 10⁻⁴⁵⁸⁶** (after the term-scale α fix this
  study forced — a zero-valued identity's chance-range is its constituent
  terms').
- **P3 MET:** the color–magnitude relation stays a labeled conjecture
  (`M_G ≈ 2.89·(BP−RP) + 1.29`), never a certificate.
- **P4:** met post-amendment; the loose-floor event is the registered
  exception and is documented above rather than hidden.

Phase 0 of the proposal is complete: adapter live, frozen-artifact fetches,
one known relation end-to-end under full discipline — plus two instrument
hardenings (absolute-precision declaration in `recover`; term-scale invariant
α) and one registered issue (loose-ε parsimony at σ=0) that Gaia forced on
day one, exactly as the proposal's methodological goals predicted.

## Issue closure (2026-07-27): the loose-ε parsimony fix

Diagnosis on the frozen C0 snapshot, engine internals exposed. Three
compounding mechanisms, none of them the assumed one:

1. **The exact-coefficient gate INVERTS at a floor-dominated ε.** At floor
   2e-4, 28 tier-1 candidates certify, including the 2-term truth. But
   `float_pinned` (perturb each coefficient ±1e-5/±1e-4 relative; cert must
   break) REJECTS the truth — its slope is honestly not identified to 1e-5 at
   that ε — while PASSING the multi-term whales, whose large mutually
   canceling coefficients are hyper-sensitive to any perturbation. The gate
   built to block dyadic garbage removed the true rival and admitted the
   approximants; coherence then saw a whale-only class.
2. **No floor separates a superset fit from the truth.** The whale bases
   CONTAIN {1, x}: lstsq spreads catalog-rounding noise over junk terms that
   mutually cancel, so whale residuals ≈ truth residuals (~1e-6). Even at the
   amended 5e-6 floor all 28 candidates certify; term-by-term dropping
   (`reduce_to_minimal`) cannot prune delicately canceling terms.
3. **Re-splits dissolved genuine ambiguity.** The amended-floor "clean" result
   was itself split luck: split 0 sees 3 rival classes (honest structural
   abstain), a later re-split's coherence collapses to one class, and the
   full-data gate — which checks only fit, not rivalry — passes it.

Fix (validated against the full suite, this snapshot at both floors, and the
NewtonBench-default 1e-12 floor which must not trip it):

- **Floor-dominated regime flag** (`sigma = 0` and
  `floor_abs > max(1e-9, 100·MACHINE_REL·median|y|)`): the per-candidate
  coefficient gate stands down (same lesson as the σ>0 PO12 move — coherence
  must see the full certifying set) and moves to the winner; the CAP-S cheap
  pre-pass stands down.
- **Refit-parsimony collapse** (`certify.refit_minimal`): each certifying
  candidate is greedily forward-selected over its OWN basis with refitting;
  a candidate whose sub-support already certifies collapses to it. Every C0
  whale collapses to the 2-term line, so coherence sees the one class that is
  actually there.
- **Sticky re-split ambiguity** (`passive`): a split that witnesses materially
  different rival classes vetoes certification by any other split.

Post-fix verdicts on the frozen snapshot: floor 2e-4 → **parametric abstain**
(whales collapse to the line; the slope is honestly unpinnable at a loose
floor); floor 5e-6 → **certified `25.68735 − (5/2)·log₁₀(flux)`, α ≤ 10⁻⁴⁶¹**
— the published result, now reached through a sound path. The full C0 study
re-run reproduces every verdict and law; the BP zeropoint snaps to an
equivalent rational (Δ ≈ 1.5×10⁻⁷, well inside the floor) with a marginally
sharper α (10⁻⁴⁶⁶ vs 10⁻⁴⁵⁹).

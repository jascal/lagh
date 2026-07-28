# Proposal: Certified Law Discovery on the Materials Project with lagh

**Drafted 2026-07-28 as the campaign's governing proposal** (third archive
campaign, after `PROPOSAL_GAIA.md` and `PROPOSAL_EXOPLANET.md`; domain
selected in-session from a shortlist of Open Catalyst / Materials Project /
ERA5, with Open Catalyst designated the follow-up campaign and ERA5 deferred
to the PDE capability arc).

## 1. Summary

Apply lagh to the Materials Project — a large database of computed
(DFT/VASP) materials properties with a documented, fixed pipeline — under
the established archive discipline: query-driven acquisition with frozen
artifacts, registration before every run, certificates with significance
bounds or explicit abstentions, zero confident-wrong.

The domain fit mirrors the exoplanet campaign: a **definitional stratum**
(columns the pipeline computes from other columns — density from structure,
per-atom normalizations, hull constructions) that should certify with
constants decoding the pipeline's unit conventions; an **empirical stratum**
(elastic/thermodynamic scaling relations) that belongs on the banded
conjecture track; and a **must-not-certify stratum** (DFT band gaps —
famously method-biased and scattered) where a certificate would be a
confident-wrong. Because every value comes from one documented pipeline
rather than a literature patchwork, the composite-table heterogeneity that
set the exoplanet tails should be largely absent — itself a registered,
testable expectation.

## 2. Data access

REST API (`api.materialsproject.org`, `X-API-KEY` from `machine/.env`,
key never committed), `summary` and `elasticity` endpoints, selected fields
only, TOP-N samples with deterministic sort. Every fetch is a frozen
artifact (exact endpoint + params + SHA-keyed snapshot) per the Gaia
adapter pattern.

## 3. Curriculum

- **C0 — calibration & identity checks.** The density identity (cell mass /
  volume, with the amu→g/cm³ constant 1.66054 decoding the unit convention);
  registered non-certification of band-gap relations; baseline certificate
  quality and floor procedure on API-served float columns.
- **C1 — elastic and thermodynamic scaling (banded conjectures).**
  Bulk-modulus–volume scaling within chemical families (Birch-type
  power laws, registered exponent bands); shear/bulk ratios; cohesive
  trends vs melting where joinable. `recover` must abstain on scattered
  cross-family relations.
- **C2 — open sweep under heavy abstention.** The Gaia-C6/exoplanet-C5
  protocol over the summary columns: scout-gated pairs, registered
  definitional triples, census; any unregistered certificate must decode
  into the pipeline documentation.
- **C3+ (with Open Catalyst as the successor campaign):** adsorption-energy
  scaling relations with theory-predicted RATIONAL slopes
  (Abild-Pedersen (4−x)/4 bond-counting) — the exact-rational machinery's
  best scientific target — once the OC20 relaxation tables are engineered
  into joinable form.

## 4. Discipline

As established: predictions frozen in git before any fetch; the v2 floor
procedure (storage precision on RAW columns, σ_rep for relative rounding);
provenance-flag subsampling where the archive marks computed values;
amendment protocol for floor mis-declarations; every certificate behind the
α ≤ 10⁻⁶ significance gate; misses reported as misses.

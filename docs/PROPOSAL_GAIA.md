# Proposal: Certified Law Discovery on Gaia DR3 with lagh
### Adaptive, Certificate-Driven Exploration of Stellar and Galactic Regularities
**Author: James (jascal). Date: July 2026.** (Committed verbatim as the
project's governing proposal; execution deltas are recorded in the per-phase
registration docs, starting with `CASE_STUDY_GAIA_C0.md`.)

lagh — a certified law discoverer with a zero-confident-wrong invariant,
machine-checked certificates or explicit abstentions, curriculum registry,
adaptive ranging acquisition, and registration discipline — applied to
selected subsets of Gaia DR3 (~1.81B sources; ~33-34M with 6D phase space;
~470M with astrophysical parameters). Adaptive acquisition maps onto Gaia's
ADQL/TAP + HEALPix-partitioned bulk access: fetch only what each curriculum
stage demands; never ingest the full catalog.

**Goals:** (1) methodological — stress-test adaptive ranging, error-model
assembly (per-source heteroscedastic uncertainties exercise the epsilon
model's `se` term for the first time on real data), certificate generation,
and curriculum progression at observational scale; (2) scientific —
independently certified rediscoveries of known relations, and search for new
or refined regularities, with abstention where the data do not support
confident claims.

**Curriculum (frozen per stage):** C0 calibration on definitional identities;
C1 photometric/astrophysical-parameter relations; C2 solar-neighborhood 6D
kinematics; C3 binary/orbital solutions (reframed at registration:
mass-DEPENDENT relation discovery — Kepler III is not a single cross-catalog
law); C4 variable-star empirical laws; C5 galactic-scale structure; C6 open
discovery under heavy abstention.

**Discipline:** every fetch is a frozen artifact with its exact ADQL;
predictions registered before every discovery run; certificates (with α) or
abstention reports close each stage before the next ranging decision.

**Phases:** 0 — adapter + C0 end-to-end (this repo, `experiments/gaia/`);
1 — C0-C2; 2 — C3-C4 + first public certificate release; 3 — C5-C6.
Compute modest; storage well under a few hundred GB at peak by construction.

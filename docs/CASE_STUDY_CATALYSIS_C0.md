# Catalysis C0 — acquisition and census (registration)

**Registered 2026-07-28, before any fetch.** Governing proposal:
`docs/PROPOSAL_CATALYSIS.md`. C0 makes NO law predictions — it freezes the
acquisition protocol and the census questions, so that C1's rational-slope
registration can be written against known sample sizes without having seen
a single energy–energy correlation.

## Protocol

- **Adapter:** GraphQL POST to `api.catalysis-hub.org/graphql`,
  frozen-artifact snapshots keyed by the exact query (the established
  pattern). Cursor pagination.
- **Fetch:** all reactions for `pubId: MamunHighT2019`, fields:
  chemical composition, facet, sites, reactants, products, reactionEnergy.
- **Census (the C0 deliverables):** (1) adsorbate species inventory and
  counts; (2) joinable-surface counts per scaling pair (surfaces where both
  the hydrogenated species and its central atom have energies, keyed by
  composition + facet); (3) duplicate-entry count per (surface, adsorbate)
  key — the σ_rep estimator (median absolute spread of duplicates), with
  the literature 0.1 eV as fallback if duplicates are < 50.
- **Firewall:** the census computes NO energy–energy regressions; scaling
  slopes are not looked at until C1's bands are frozen.

## Census results

(after the fetch)

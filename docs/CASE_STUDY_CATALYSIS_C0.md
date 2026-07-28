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

**Acquisition amendments (logged, 2026-07-28):** the Catalysis-Hub GraphQL
API now requires an API key whose registration flow appears to require
institutionally-affiliated ORCID — unavailable to this program. (The API
answers requests carrying the public website's Origin header; that
workaround is REJECTED without explicit authorization as a circumvention of
the stated access policy.) Fallback attempted: the public Ulissi-group
research mirror (`ulissigroup/uncertainty_benchmarking`,
`cathub.pkl` @ sha256 861a4569…, itself a pull of the CC-BY Mamun data).
Extraction required two logged fixes (supercell-dependent raw formulas →
reduced compositions; the systems lists carry no clean slab → strip
adsorbate elements from the largest system).

## Census results (2026-07-28, frozen table
`experiments/catalysis/data/mamun_mirror_energies.csv`, 30,420 rows)

- Mirror coverage by adsorbate: H 10,074 · N 9,000 · C 6,664 · O 3,534 ·
  OH 1,148 — **no CHₓ, NH, or SH** (the Ulissi pull fetched only its
  paper's adsorbate set).
- **The mirror is a per-adsorbate patchwork of DISJOINT surface slices**:
  13 distinct alloy compositions carry O, 11 carry OH, and the
  intersection is EMPTY (O on AgTa/Cr₃Sn/CuAg…; OH on Ag/Ag₃Bi/Ag₃Pt…).
  Zero of the six registered scaling pairs are joinable from this mirror.
- Site-multiplicity: 92 multi-entry (surface, adsorbate) keys, median
  spread 0.22 eV — confirming the most-stable-site convention matters and
  informing the eventual σ declaration.
- **C0 verdict: the mirror cannot serve C1.** The firewall held throughout
  (no energy–energy regression was computed). Remaining acquisition paths,
  in order of preference: (a) a public CatHub GitHub issue requesting an
  independent-researcher API path (user-mediated); (b) the Materials Cloud
  raw deposit (25 GiB of QE outputs, CC-BY, fully open — the heavy but
  autonomous path); (c) the Origin-header route only under explicit
  authorization.

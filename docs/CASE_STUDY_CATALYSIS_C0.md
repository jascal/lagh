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

## The path-2 rebuild (2026-07-28; user directed: raw deposit only, no
CatHub engagement)

All 25.3 GiB fetched from Materials Cloud (CC-BY); ~92k QE logs
stream-parsed with ASE (`parse_raw.py`): references 3,648 ok/414 fail —
**no parseable gas-phase molecules** (so the builder works in
ΔE* = E(slab+species) − E(slab) space; gas terms are species constants
that cannot touch scaling slopes); adsorbate archives 76,738 adsorbed rows
(per-archive failure rates 20–56%, consistent with restart fragments and
unconverged runs — a census fact, gated below). Builder (`build_table.py`):
clean slabs by vacuum height (1,614 compositions), pairing on EXACT metal
counts, most-stable site per (surface, species).

**The rebuilt table serves all six pairs**: CH/CH₂/CH₃ vs C at 255/249/249
joinable surfaces, NH vs N at 255, OH vs O at 239, SH vs S at 255
(`data/mamun_rebuilt_energies.csv`, committed).

**Validation gate (mirror overlap, per-species offset calibrated):**
median |residual| after offset — **O 0.0007 eV, OH 0.0000 eV, C 0.010 eV**
(the rebuild reproduces the canonical processing exactly where coverage
aligns; the large offsets are the gas-reference constants, as designed);
H 0.22 and N 0.15 miss the 0.05 gate with a diagnosis: min-over-different-
site-sets asymmetry (the mirror slices are tiny and our parse drops
unconverged sites), not a systematic energy shift. Gate: ACCEPTED for C1
with the H/N caveat recorded; σ for C1's conjecture track declared from
the measured site-multiplicity spread (median 0.22 eV).

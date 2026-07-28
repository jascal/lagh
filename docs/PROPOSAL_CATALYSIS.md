# Proposal: Certified Law Discovery in Catalysis — the Rational-Slope Campaign

**Drafted 2026-07-28 as the campaign's governing proposal** (fourth
campaign). Successor to `PROPOSAL_MATERIALS.md` per its own closing note.

## 1. The scientific target

Adsorption-energy **scaling relations** (Abild-Pedersen et al. 2007): across
surfaces, the adsorption energy of a partially hydrogenated species scales
linearly with that of its central atom, with slope predicted by BOND
COUNTING to be a **simple rational**: γ = (x_max − x)/x_max, where x_max is
the central atom's maximum valency. Concretely: CH vs C → 3/4; CH₂ vs C →
1/2; CH₃ vs C → 1/4; NH vs N → 2/3; OH vs O → 1/2; SH vs S → 1/2.

This is the best match in any proposed domain between an open scientific
question and this instrument's exact-rational machinery: theory predicts
specific small rationals, data carries real DFT scatter (~0.1–0.3 eV), and
the honest deliverable is a REGISTERED BAND TEST of the rational values —
conjecture-track throughout, with certification explicitly forbidden (a
certificate on scattered chemistry would be confident-wrong).

## 2. Data (amendment from the materials proposal, logged)

The materials proposal named Open Catalyst (OC20/OC22) as this campaign's
substrate. OC20 ships as LMDB relaxation trajectories requiring substantial
engineering before a joinable (adsorbate × surface) energy table exists.
The same scientific target is served directly by **Catalysis-Hub.org**
(GraphQL API, query-driven — the established adapter pattern) and
specifically the **Mamun et al. 2019 high-throughput dataset**
(`MamunHighT2019`): adsorption energies of C, CH, CH₂, CH₃, N, NH, O, OH,
S, SH across thousands of alloy surfaces — the scaling-relation table,
pre-built. OC20 remains the designated later stage (transfer of the same
registrations to a second, larger substrate).

## 3. Curriculum

- **C0 — acquisition & census (no law predictions).** GraphQL adapter with
  frozen artifacts; dataset shape census (species coverage, per-surface
  join counts); declared-σ procedure fixed (DFT replication scatter, from
  duplicate surface entries where present, else the literature 0.1 eV).
- **C1 — the rational-slope registration (the centerpiece).** Frozen bands
  around the six bond-counting rationals; consistency tests of the exact
  values; pooled recover must abstain; zero certificates permitted.
- **C2 — BEP and beyond (stretch).** Brønsted–Evans–Polanyi relations if
  activation energies are joinable; transfer to OC20.

## 4. Discipline

As established across three closed campaigns: registrations frozen in git
before fetches/runs; logged amendments; census before claims; every miss
reported; zero confident-wrong.

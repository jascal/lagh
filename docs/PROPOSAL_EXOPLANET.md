# Proposal: Certified Empirical Law Discovery on the NASA Exoplanet Archive with lagh

**Author: James (jascal). Date: July 2026.** Instrument: `jascal/lagh`
(certified law discoverer). Data source: NASA Exoplanet Archive (Planetary
Systems Composite Parameters + Stellar Hosts). (Committed verbatim as the
campaign's governing proposal; execution deltas are recorded in the per-phase
registration docs, following the `PROPOSAL_GAIA.md` precedent.)

## 1. Executive Summary

We propose to apply lagh to the NASA Exoplanet Archive, focusing on the
Planetary Systems Composite Parameters (PSCompPars) table and associated
stellar properties. The goal is to discover, refine, or cleanly abstain on
empirical relations that govern exoplanet populations — mass–radius,
radius–period, density trends, the radius valley, and multiplicity/spacing
laws — while strictly preserving lagh's zero-confident-wrong invariant,
adaptive ranging, curriculum discipline, and formal certificate (or
abstention) output.

The archive is the right scale and character for the next lagh campaign:
large enough for meaningful adaptive acquisition and conditioning, small
enough to remain fully controllable, and scientifically active enough that
certified results (or well-documented refusals) would constitute new
contributions. Lessons from the Gaia Phase 1–2 work — especially the
diagnosis of catalogue-level inconsistencies in published parameters —
transfer directly.

## 2. Motivation and Fit

lagh has already demonstrated, on Gaia, the ability to:

- Recover classical relations under registration discipline,
- Detect when published columns fail to satisfy expected physical identities
  (FLAME medians vs Stefan–Boltzmann),
- Keep conjectures properly labeled and abstain rather than force
  certificates,
- Operate via adaptive, query-driven acquisition against a remote archive.

The exoplanet domain offers a natural continuation with higher potential for
novel scientific output. Empirical relations in exoplanet demographics remain
actively debated (radius valley location and slope, mass–radius breaks,
density trends with irradiation or stellar type, period-ratio preferences).
A system that produces machine-checked certificates or explicit,
quantitative abstentions would be a distinctive contribution.

## 3. Relevant Data Products

Primary tables (NASA Exoplanet Archive):

| Table | Content | Approximate scale | lagh use |
|---|---|---|---|
| Planetary Systems Composite Parameters (PSCompPars) | Homogenized planet parameters (mass, radius, period, semi-major axis, density, Teq, etc.) with uncertainties and flags | ~6,300+ confirmed planets | Core relations |
| Planetary Systems (PS) | Default parameter set per planet | Similar | Cross-check / alternative homogenization |
| Stellar Hosts | Host-star mass, radius, Teff, metallicity, etc. | Thousands of hosts | Conditioning variables |

Key measurable quantities for relations: planetary mass (M_p), radius (R_p),
orbital period (P), semi-major axis (a), equilibrium temperature /
insolation, bulk density, stellar mass (M_star), metallicity [Fe/H], and
multiplicity indicators.

Access is via the Archive's TAP service (programmatic ADQL-style queries) and
bulk downloads of selected columns — ideal for lagh's adaptive ranging.

## 4. Technical Approach — Adaptive Acquisition

Exactly as with Gaia:

- All acquisition is query-driven or selective bulk. No full-archive
  materialization is required.
- lagh's ranging decisions (period cuts, mass/radius precision thresholds,
  stellar-type windows, multiplicity filters) become TAP queries or filtered
  table subsets.
- Each curriculum stage freezes an artifact (selected rows + columns +
  quality cuts) before any discovery step, following the
  predictions-before-run harness.
- Uncertainties and literature flags are retained so the assembled error
  model can operate.

This keeps the working set modest while allowing controlled expansion.

## 5. Curriculum Design

A concrete, pre-registered curriculum (C0–C5 style):

- **C0 – Calibration & identity checks.** High-precision subsample with both
  mass and radius. Test basic physical identities and catalogue consistency
  (analogous to the Gaia Stefan–Boltzmann diagnosis). Establish baseline
  certificate quality.
- **C1 – Global mass–radius and radius–period relations.** Recover or refine
  the primary empirical sequences; locate and characterize the radius valley
  under registration discipline.
- **C2 – Conditioning on stellar properties.** Split or regress on stellar
  mass, metallicity, or Teff. Seek certified conditional slopes or cleanly
  abstain where the data do not support them.
- **C3 – Density and irradiation trends.** Bulk density versus insolation /
  equilibrium temperature; look for breaks or scaling relations.
- **C4 – Multi-planet systems.** Period ratios, spacing statistics, and
  mutual Hill-radius relations in systems with N_p >= 2.
- **C5 – Open / stress exploration.** Broader or lower-precision samples;
  lagh must rely more heavily on abstention and the error model.

Each stage ends with a registered certificate package or explicit abstention
report before the next ranging decision.

## 6. Experimental Plan & Phases

**Phase 0 (1–2 weeks)**

- Implement TAP adapter and selective-column loader for PSCompPars + Stellar
  Hosts.
- Reproduce a small set of known literature relations end-to-end under full
  registration discipline.
- Document any immediate catalogue inconsistencies.

**Phase 1 (3–5 weeks)**

- Execute C0–C2.
- Quantify certificate quality, abstention rate, and the effect of precision
  cuts and stellar conditioning.

**Phase 2 (4–6 weeks)**

- C3–C4 (density/irradiation + multi-planet).
- Compare lagh outputs against current literature consensus; release first
  public certificate set.

**Phase 3 (ongoing)**

- C5 expansion and iterative refinement.
- Integration of new Archive updates as they appear.
- Public release of curriculum registry, frozen artifacts, and certificates
  (companion to the Gaia CERTIFICATES.md).

## 7. Expected Outcomes

**Methodological**

- Demonstration that lagh's adaptive, certificate-driven workflow transfers
  cleanly from stellar catalogues to exoplanet demographics.
- Further hardening of the acquisition and error-model layers on
  heterogeneous, literature-homogenized data.
- Reusable TAP / tabular adapter.

**Scientific**

- Independently certified (or cleanly refused) statements on mass–radius,
  radius-valley, density, and multiplicity relations.
- Quantitative diagnosis of any internal inconsistencies in the composite
  parameters, analogous to the FLAME median result.
- Public, machine-readable certificate packages that the exoplanet community
  can inspect, reproduce, or challenge.

## 8. Resources & Feasibility

- **Data volume:** Tables are modest (megabytes to low tens of megabytes for
  selected columns). Adaptive queries keep the working set small.
- **Access:** Public TAP service + bulk CSV/VOTable downloads. No special
  authorization required.
- **Compute:** Laptop-scale or modest cloud instance; far lighter than Gaia
  HEALPix or simulation workloads.
- **Software:** lagh core + Archive TAP client (or simple pandas/polars +
  requests) + standard scientific Python.

## 9. Risks & Mitigations

| Risk | Mitigation |
|---|---|
| Heterogeneous literature origins and homogenization artifacts | Explicit quality/precision cuts; retain provenance flags; treat catalogue consistency itself as a discovery target (Gaia lesson) |
| Small-number statistics in some conditioned bins | Pre-registered minimum sample sizes; forced abstention below thresholds |
| Rapid Archive updates | Version-pin frozen artifacts; re-run certificates on new releases under the same registration discipline |
| Selection biases (detection method, stellar type) | Curriculum stages that isolate or explicitly condition on detection technique and stellar parameters |

## 10. Timeline (indicative)

- Weeks 1–2: Adapter + C0 validation
- Weeks 3–7: C1–C2 + first certificate release
- Weeks 8–14: C3–C4
- Month 4+: C5 and ongoing Archive monitoring

## 11. Conclusion

The NASA Exoplanet Archive offers a high-value, tractable domain in which
lagh can move from rediscovery and consistency diagnosis (Gaia) toward
certified empirical statements on still-open questions in planetary
demographics. The combination of adaptive TAP acquisition, strict
registration discipline, and zero-confident-wrong output is well-matched to
the character of the data and to the scientific needs of the field.

This proposal is ready for immediate implementation as the next lagh
campaign after the current Gaia Phase 3 run.

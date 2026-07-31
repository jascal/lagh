# Case study: C1 on a trapped bead — what a real instrument will and will not certify

**Run 2026-07-30.** `experiments/tweezers/run_c1.py`,
`experiments/results/tweezers_c1.json`. LUMICKS C-Trap passive calibration record
(Zenodo 14726586, CC-BY-4.0), back-focal-plane interferometry, 78125 Hz, 16.0 s,
Force 1x and Force 1y of one 4.89 µm bead.

Registered in `PROPOSAL_STOCHASTIC_REAL.md` §9. C0 screened this record IN after the
fluorescence record was refused. C1 asks what can be certified from it.

**Headline: the diffusion cannot be certified from this record, for two independent
reasons, and the drift can be recovered to 0.8% — but not by the estimator the arc
currently uses.** Zero confident-wrong: nothing certified, every verdict an
abstention.

## 1. The passive displacement calibration is circular for b²

`Rd`, the volts→nm scale, is not measured against a length standard on this record.
It is **defined** by setting the observed diffusion equal to its Stokes–Einstein
value:

    Rd  ==  sqrt( (k_B T / gamma_0) / D_volts )

Verified on **all 8 calibration items in the file, maximum deviation 0.0e+00** — not
approximately, identically. So a certificate that b² = 2k_BT/γ₀, computed from a
position signal scaled by Rd, restates the calibration's own definition. It cannot
fail except through estimator disagreement.

This is the FLAME circularity that got reanalysis rejected in `DIRECTION_PDE.md`,
arriving in a new costume, and `PROPOSAL_STOCHASTIC_REAL.md` §3 walked straight into
it: the proposal argued the C-Trap "breaks the circularity" on the strength of the
*active* record's nanostage, then ran C1 on the *passive* one. The break is real but
it lives in a different file.

## 2. Anti-alias filtering destroys quadratic variation

The instrument fits its Lorentzian-times-diode model over **100 Hz – 23 kHz**. Our
position signal agrees with that model **to 0.06% in-band** (ratio 0.9994 on the clean
axis), so the unit chain is right and the calibration is excellent at what it claims.

Outside the band it is a different story:

| | measured / model |
|---|---|
| in fit band, 100 Hz – 23 kHz | **0.9994** |
| above 23 kHz | **0.777** |
| at 0.97 × Nyquist (37.9 kHz) | **6.0e-05** |

And a quadratic variation is `∫ (2πf)² PSD df` — weighted by **f²**, so it is
dominated by the top of the band and above. **24.2% of the model's quadratic
variation lives above the 23 kHz fit ceiling**, exactly where the chain attenuates.

A PSD fit weights its band roughly flat and is blind to this. A QV estimator is not.
So the two disagree by the attenuation the fit never saw — and no amount of data
fixes it, because the information was removed before storage.

**The general lesson, which is what this record was worth acquiring for: realized
quadratic variation is destroyed by anti-alias filtering, and every real instrument
has one.** The arc's simulated systems have none, so every diffusion result to date
was obtained in the one regime where this cannot bite.

Measured b², mid-interval over truth:

| axis | raw | after diode deconvolution |
|---|---|---|
| Force 1x | 0.344 | **0.580** |
| Force 1y | 0.328 | **0.554** |

Deconvolving the *modelled* part of the attenuation (the photodiode, using the
instrument's own published α ≈ 0.42 and f_diode ≈ 8.4 kHz) recovers about half the
gap, which confirms the mechanism. The rest is the unmodelled anti-alias filter.

The diffusion's partial determination therefore **fails to cover** on all four runs.
Not a confident-wrong — nothing was certified — but a coverage miss with an
identified cause.

## 3. The f = x²/2 family makes the drift inherit the diffusion's bias

`DIRECTION_STOCHASTIC.md` records Level 0's finding that a stationary drift is vacuous
with `f = x` alone, because `E[a(X)] = 0` under any stationary law, and that `f = x²/2`
is what makes it accumulate. That fix has a consequence nobody noticed, because it
is invisible without a band-limited instrument.

With `f = x²/2` the Itô identity is

    ∫φ d(x²/2)  =  −θ ∫φ x² dt  +  (b²/2) ∫φ dt  +  martingale

and for a stationary path the left side is a small boundary term, so the fit is
essentially `θ ≈ b² / (2⟨x²⟩)`. **The drift is determined through b².** If b² is
attenuated by a factor r, so is θ:

| Force 1y, deconvolved | value | ratio to truth |
|---|---|---|
| b², mid-interval | 1.11e5 nm²/s | **0.554** |
| θ, weak form (`("1","x")`, half 20000) | 1778 /s | **0.539** |

The two ratios agree to 3%. That is the coupling, measured.

**And it is the estimator's fault, not the data's.** θ from the autocorrelation decay
— a genuine timescale, never routed through b², invariant to any rescaling of the
signal and read at lags beyond the detector's filter — recovers the truth from the
**raw** data:

| axis | θ from ACF | truth | ratio |
|---|---|---|---|
| Force 1y, raw | 3323.4 /s | 3297.0 /s | **1.008** |
| Force 1y, deconvolved | 3333.7 /s | 3297.0 /s | 1.011 |

0.8% from data the weak form reads 46% low on. The quantity was always recoverable;
the arc's estimator routes it through the one thing this instrument destroys.

**This is the actionable engine finding.** `f = x²/2` was adopted for identifiability
and is correct in simulation. On a band-limited instrument it transfers the entire
diffusion bias into the drift. Any real-data drift claim from `build_rows` needs
either a band-loss declaration or an estimator that does not pass through b².

## 4. Force 1x carries a contaminant the calibration cannot see

θ from the ACF is **125.2 /s against a truth of 3138.2 — a factor of 25**. The x axis
is dominated by a slow contaminant (stage drift or laser pointing), which inflates its
position variance to 1.9× thermal while leaving the increments normal.

The instrument's calibration does not catch it because **its fit range starts at
100 Hz** and the contaminant lives below that. The ACF timescale is the gate that
does catch it, and it costs one line. Force 1y is clean by the same test.

Per-axis quality gating is therefore not optional on this dataset, and "the
manufacturer calibrated it" is not a substitute — the calibration's own fit range
excludes the failure.

## 5. Provenance: the unit chain spans two calibrations

The stored pN values were produced with the calibration **active at acquisition**
(item 5), while the only calibration describing this trap — the one whose voltage
window lies inside the record — is **item 6**. Dividing stored pN by item 6's κ
silently mixes them: **−4.75% in position, −9.28% in b²**.

The honest chain is pN → V with the *applied* Rf, V → nm with the *derived* Rd, and
it is what `adapter.bfp_position_nm` implements. This is the fourth campaign in a row
to find that an archive's columns come from more than one pipeline state.

## 6. Verdicts

| target | verdict | why |
|---|---|---|
| b² (diffusion) | **abstain**, interval does not cover | circular calibration; 24% of QV above the fit band, attenuated |
| θ (drift), weak form, 4-term library | **abstain**, interval covers but spans zero | signal/band 0.11–0.31; honestly undetermined |
| θ (drift), weak form, narrow library | **abstain**, interval tight and does NOT cover | the b² coupling of §3 |
| θ (drift), ACF timescale | **0.8% on the clean axis** | not a certificate — an estimator the arc does not yet have |

Zero confident-wrong across all four runs, which is the ranking's first key. The
instrument declined every claim it could not support, and the two coverage misses are
both on partial determinations rather than certificates.

## 7. What C2 must do

1. **Move to the active record.** `near_surface_active_calibration.h5` drives the
   nanostage at a known amplitude and records `Nanostage position`, giving Rd from a
   length standard rather than from the thermal motion. That is the only way b²
   becomes a physics claim on this instrument. §1 is not fixable on a passive record.
2. **Declare the band loss.** The QV estimator needs the anti-alias response as a
   declared input, the way `sigma_obs` is — measured from the PSD's high-frequency
   ratio rather than assumed. Without it every real-instrument b² is biased low by a
   factor no row count reduces.
3. **Give the arc a b²-free drift estimator**, and register `f = x²/2`'s coupling in
   `DIRECTION_STOCHASTIC.md` alongside the identifiability reason it was adopted for.
4. **Gate axes on the ACF timescale** before any of the above.

Tag: **empirical** throughout. Nothing here is `proved`; §1's identity is exact
arithmetic on the file's own fields and is the one claim that could be promoted with
a stated domain.

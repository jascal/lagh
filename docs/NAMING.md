# Naming decision (2026-07-21)

**Keep `lagh`. Not `sr-orca`.** Considered and declined, for two independent reasons:

1. **Not Orca-shaped.** The Orca family is spec-DSLs you author *in*, statically verified then
   compiled to a heavy backend (`i-orca:Isabelle :: n-orca:PyTorch :: orca-lang:XState`). lagh
   authors nothing — it is an autonomous instrument that *emits* the law, with the certifier as the
   product, not a pre-check before a backend. By shape it is a sibling of `rosetta`/`ergo`
   (checkers), which share the untrusted-proposer/sound-checker discipline but are not Orca-named
   for the same reason.
2. **"sr" mislabels it.** lagh's identity is being *not* symbolic regression — exact/certified/
   abstaining vs approximate/fitted/always-answering (see `HARD_CLASS_ANALYSIS.md` on refusing
   `x^e`). Branding it "sr" stamps it with the name of the thing it is defined against.

By lineage lagh descends from `wyly ← pil` (the SAE-interpretability line), not the Orca line.

**Actual family:** the Gaelic-named instrument/serving cluster — `sgiandubh`, `claymore`, `lagh` —
which is exactly where the tool-shape places it (a certified bounded expert federating under
claymore). Cross-project kinship is real but is *the certified-tool federation* (i-orca certifies
proofs, lagh laws, rosetta equivalences, ergo deductions — united by discipline), not the Orca
language family. Subtitle for legibility: *"lagh — certified law discovery (a checker in the
i-orca/rosetta discipline)."*

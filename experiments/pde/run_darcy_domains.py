"""Two Darcy phases, two domain-qualified verdicts, and the refusal to merge them.

    .venv/bin/python experiments/pde/run_darcy_domains.py

The variable-coefficient route registered 2026-07-29 (`DIRECTION_PDE.md` (c)) is
the HONEST DOMAIN RESTRICTION: certify where the coefficient field is locally
constant and report where. What makes that a claim rather than a filter is the
composition rule -- a verdict established in the high-conductivity phase must not
silently merge with one established in the low phase, because their conjunction
is defined only where both were established, and nothing in this program can
decide that two differently-worded regions are the same set.

`certify.conjoin_determination` enforces it. This script checks the enforcement
on REAL RUN OUTPUT rather than on constructed records: it reads the two Darcy
verdicts from `pdebench.json` and conjoins them three ways.

  1. high phase with itself      -> composes (same predicate)
  2. high phase with low phase   -> REFUSES, naming both domains
  3. either with an unqualified record -> composes, and the result inherits the
     qualifier, because a record without one claims the whole field

Both phases abstain structurally at the declared band, which is the verdict the
first Darcy pass already reported; nothing here changes it. What changed is that
the abstain now CARRIES its domain in the shared vocabulary instead of in a
free-text sibling field, so a consumer reading `partial` can no longer mistake a
statement about 16% of the field for a statement about the field.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent.parent))
from lagh.certify import conjoin_determination                    # noqa: E402

RESULTS = Path("experiments/results/pdebench.json")
OUT = Path("experiments/results/darcy_domains.json")
HIGH, LOW = "darcy_beta0.1_qualified", "darcy_beta0.1_lowphase_qualified"


def main():
    d = json.loads(RESULTS.read_text())
    missing = [k for k in (HIGH, LOW) if k not in d]
    if missing:
        print(f"REFUSED: run run_pdebench_darcy.py first; missing {missing}")
        return 1
    hi = d[HIGH]["equation"]["partial"]
    lo = d[LOW]["equation"]["partial"]
    plain = {"status": "certified", "components":
             {"1/a": {"kind": "interval", "lo": -0.11, "hi": -0.09,
                      "resolved": True}}}
    for r, name in ((hi, HIGH), (lo, LOW)):
        print(f"{name}: status={r['status']} components={len(r['components'])} "
              f"coverage={r['qualifier']['coverage']:.3f}")
        print(f"   domain: {r['qualifier']['predicate']}")

    same = conjoin_determination([hi, hi])
    cross = conjoin_determination([hi, lo])
    with_plain = conjoin_determination([hi, plain])
    print(f"\nhigh + high  -> {same['status']}")
    print(f"high + low   -> {cross['status']}: {cross.get('refusal')} "
          f"{cross.get('domains')}")
    print(f"high + unqualified -> {with_plain['status']}, qualifier kept: "
          f"{with_plain.get('qualifier', {}).get('predicate')!r}, "
          f"components {sorted(with_plain['components'])}")

    ok = (same["status"] != "refused" and cross["status"] == "refused"
          and with_plain.get("qualifier", {}).get("predicate")
          == hi["qualifier"]["predicate"])
    print(f"\n{'OK' if ok else 'FAILED'}: the domain rule holds on real output")
    OUT.write_text(json.dumps(
        {"high": hi, "low": lo, "same_domain": same, "cross_domain": cross,
         "with_unqualified": with_plain, "rule_holds": bool(ok)}, indent=1))
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())

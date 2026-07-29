"""C3's verify track: does a certified SYSTEM forecast?

The C2 question, asked of systems: integrate the certified system forward, from
an initial condition no stage of the pipeline has seen, in the weak form's own
vocabulary, and ask whether every field stays inside the declared band. The
certificate claims a family of systems (interval parameters per equation), so
the forecast is a family too.

Reads the certificates C3 actually produced (experiments/results/pde_c3.json).
"""
from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).parent.parent.parent))
from experiments.pde import systems_fields as S                    # noqa: E402
from experiments.pde.run_c3_systems import (BRU_TERMS, FHN_TERMS,  # noqa: E402
                                            LINEAR_TERMS, SW_TERMS)
from experiments.pde.verify import verify_system                   # noqa: E402

IN = Path("experiments/results/pde_c3.json")
OUT = Path("experiments/results/pde_c3_verify.json")
NX, NT, TMAX = 256, 41, 0.4
WRONG_FACTOR = 10.0            # a system this far outside its own intervals


def two_field(q):
    return {"u": q[0], "v": q[1]}


def shallow_fields(q):
    """Shallow water evolves (h, hu); the FIELDS are (h, u) -- the caller
    declares the map, because nothing in the vocabulary can guess it."""
    return {"h": q[0], "u": q[1] / q[0]}


STAGES = {
    "S1_linear": {
        "terms": LINEAR_TERMS, "to_fields": two_field,
        "field_of": {"u:u_t": "u", "v:u_t": "v"},
        "fresh": lambda: S.linear_pair(101, periodic=True, nx=256),
    },
    "S2a_fhn": {
        "terms": FHN_TERMS, "to_fields": two_field,
        "field_of": {"u:u_t": "u", "v:u_t": "v"},
        "fresh": lambda: S.fitzhugh_nagumo(101)[:2],
    },
    "S2b_brusselator": {
        "terms": BRU_TERMS, "to_fields": two_field,
        "field_of": {"u:u_t": "u", "v:u_t": "v"},
        "fresh": lambda: S.brusselator(101)[:2],
    },
    "S3_shallow": {
        "terms": SW_TERMS, "to_fields": shallow_fields,
        "field_of": {"h_t": "h", "(hu)_t": "h"},
        "fresh": lambda: S.shallow_water(101)[:2],
    },
}


def intervals_of(entry):
    """{target: {term: (lo, hi)}} from the certificate, degenerate where the
    rung reported an exact value."""
    out = {}
    for eq in entry["equations"]:
        if not eq.get("certified"):
            return None
        ivs = {}
        for k, v in (eq.get("intervals") or {}).items():
            c = (eq.get("coefficients") or {}).get(k)
            ivs[k] = tuple(v) if v is not None else (c, c)
        out[eq["target"]] = ivs
    return out


def main():
    src = json.loads(IN.read_text())
    res = {}
    for key, entry in src.items():
        stage = next((s for s in STAGES if key.startswith(s)), None)
        if stage is None or "sigma" not in key or not entry.get("equations"):
            continue
        ivs = intervals_of(entry)
        if ivs is None:
            continue
        sigma = float(key.split("_sigma")[1])
        spec = STAGES[stage]
        fields, coords = spec["fresh"]()
        x, t = coords[0], coords[1]
        terms = {tm.name: tm for tm in spec["terms"]}
        clean = {k: v.copy() for k, v in fields.items()}
        if sigma > 0:
            rng = np.random.default_rng(7)
            fields = {k: v + rng.normal(0, sigma, v.shape)
                      for k, v in fields.items()}
        q0 = np.stack([np.asarray(terms[tg].g(fields), float)[:, 0]
                       for tg in ivs])
        t0 = time.time()
        r = verify_system(fields, q0, x, t, ivs, terms, spec["to_fields"],
                          spec["field_of"], sigma=sigma, rtol=1e-9,
                          fields_clean=clean)
        r["seconds"] = round(time.time() - t0, 1)
        r["intervals"] = {tg: {k: list(v) for k, v in d.items()}
                          for tg, d in ivs.items()}
        res[key] = r
        print(f"{key:30s} law={'OK  ' if r.get('verified') else 'FAIL'} "
              f"data={'OK  ' if r.get('data_verified') else 'FAIL'} "
              f"outside={r.get('n_outside')}/{r.get('n_points')} "
              f"env={r.get('envelope_width_med', 0):.2e} "
              f"solver={r.get('solver_bound', 0):.1e} "
              f"ic={r.get('ic_noise_bound', 0):.1e} {r['seconds']}s",
              flush=True)

        # the control: every coefficient pushed WRONG_FACTOR half-widths off
        # centre must FAIL, or the track carries no information
        wrong = {}
        for tg, d in ivs.items():
            wrong[tg] = {}
            for k, (lo, hi) in d.items():
                c = 0.5 * (lo + hi)
                hw = max(0.5 * (hi - lo), abs(c) * 1e-6)
                wrong[tg][k] = (c + WRONG_FACTOR * hw, c + WRONG_FACTOR * hw)
        rw = verify_system(fields, q0, x, t, wrong, terms, spec["to_fields"],
                           spec["field_of"], sigma=sigma, rtol=1e-9,
                           fields_clean=clean)
        res[key + "_WRONG"] = {"verified": rw.get("verified"),
                               "law_n_outside": rw.get("law_n_outside"),
                               "n_points": rw.get("n_points"),
                               "refusal": rw.get("refusal")}
        print(f"{'  ^ 10x-wrong control':30s} "
              f"{'VERIFIED (BAD)' if rw.get('verified') else 'FAILED (good)':15s}"
              f" law_outside={rw.get('law_n_outside')}/{rw.get('n_points')} "
              f"{rw.get('refusal') or ''}", flush=True)
    OUT.write_text(json.dumps(res, indent=1))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

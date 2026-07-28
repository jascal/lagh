"""Materials Project REST adapter (docs/PROPOSAL_MATERIALS.md): query-driven
acquisition with frozen-artifact snapshots. API key from machine/.env
(MP_API_KEY), never committed."""
from __future__ import annotations

import hashlib
import json
from pathlib import Path

DATA = Path(__file__).parent / "data"
API = "https://api.materialsproject.org"


def _key() -> str:
    env = Path(__file__).parent.parent.parent / "machine" / ".env"
    for line in env.read_text().splitlines():
        if line.startswith("MP_API_KEY="):
            return line.split("=", 1)[1].strip()
    raise RuntimeError("MP_API_KEY not found in machine/.env")


def fetch(path: str, params: dict, tag: str):
    """GET an endpoint; snapshot the JSON docs + exact request to a frozen
    artifact. Idempotent by (path, params) hash."""
    DATA.mkdir(parents=True, exist_ok=True)
    canon = json.dumps({"path": path, "params": params}, sort_keys=True)
    h = hashlib.sha256(canon.encode()).hexdigest()[:12]
    base = DATA / f"{tag}_{h}"
    out = base.with_suffix(".json")
    if out.exists():
        return json.loads(out.read_text())["data"], out
    import requests
    docs, page_params = [], dict(params)
    while True:
        r = requests.get(f"{API}{path}", params=page_params,
                         headers={"X-API-KEY": _key()}, timeout=120)
        r.raise_for_status()
        j = r.json()
        docs.extend(j.get("data", []))
        if len(docs) >= int(params.get("_limit", 0)) or not j.get("data"):
            break
        page_params["_skip"] = len(docs)
    out.write_text(json.dumps({"request": json.loads(canon),
                               "n": len(docs), "data": docs}, indent=1))
    return docs, out


C0_PARAMS = {
    "_fields": "material_id,composition,nsites,volume,density,"
               "formation_energy_per_atom,energy_per_atom,band_gap",
    "energy_above_hull_max": 0.0,
    "_limit": 500,
    "_sort_fields": "material_id",
}

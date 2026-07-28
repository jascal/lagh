"""Catalysis-Hub GraphQL adapter (docs/PROPOSAL_CATALYSIS.md): query-driven
acquisition with frozen-artifact snapshots."""
from __future__ import annotations

import hashlib
import json
from pathlib import Path

DATA = Path(__file__).parent / "data"
API = "https://api.catalysis-hub.org/graphql"


def _key() -> str:
    env = Path(__file__).parent.parent.parent / "machine" / ".env"
    for line in env.read_text().splitlines():
        if line.startswith("CATALYSIS_HUB_API_KEY="):
            return line.split("=", 1)[1].strip()
    raise RuntimeError("CATALYSIS_HUB_API_KEY not found in machine/.env "
                       "(obtain at https://api.catalysis-hub.org/auth/login)")


def fetch_reactions(pub_id: str, tag: str, page_size: int = 500):
    """Fetch all reactions for a publication via cursor pagination; snapshot
    to a frozen artifact. Idempotent by (pub_id, page_size)."""
    DATA.mkdir(parents=True, exist_ok=True)
    canon = json.dumps({"pubId": pub_id, "page": page_size}, sort_keys=True)
    h = hashlib.sha256(canon.encode()).hexdigest()[:12]
    out = DATA / f"{tag}_{h}.json"
    if out.exists():
        return json.loads(out.read_text())["data"], out
    import requests
    q = """
    query($pubId: String!, $first: Int!, $after: String) {
      reactions(pubId: $pubId, first: $first, after: $after) {
        totalCount
        pageInfo { hasNextPage endCursor }
        edges { node {
          chemicalComposition surfaceComposition facet sites coverages
          reactants products reactionEnergy
        } }
      }
    }"""
    docs, after = [], None
    while True:
        r = requests.post(API, json={"query": q, "variables": {
            "pubId": pub_id, "first": page_size, "after": after}},
            headers={"X-API-Key": _key()}, timeout=180)
        r.raise_for_status()
        j = r.json()
        if "errors" in j:
            raise RuntimeError(str(j["errors"])[:300])
        blk = j["data"]["reactions"]
        docs.extend(e["node"] for e in blk["edges"])
        if not blk["pageInfo"]["hasNextPage"]:
            break
        after = blk["pageInfo"]["endCursor"]
    out.write_text(json.dumps({"request": json.loads(canon),
                               "n": len(docs), "data": docs}, indent=1))
    return docs, out

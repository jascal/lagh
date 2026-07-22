"""Provider-agnostic LLM config for the scientist harness.

The harness talks to an OpenAI-compatible chat API, so one config -- `base_url`,
`api_key`, `model` -- serves MiniMax M3, Nous Hermes, or any compatible endpoint.
Credentials load from the environment, backfilled from a gitignored `.env` at the
repo root. No external dependency (a tiny `.env` parser; `python-dotenv` not required).
"""

from __future__ import annotations

import os
from dataclasses import dataclass
from pathlib import Path

_REPO_ROOT = Path(__file__).resolve().parents[2]


def load_dotenv(path: str | Path | None = None) -> None:
    """Backfill os.environ from a `.env` file (existing env vars win)."""
    p = Path(path) if path else _REPO_ROOT / ".env"
    if not p.exists():
        return
    for line in p.read_text().splitlines():
        line = line.strip()
        if not line or line.startswith("#") or "=" not in line:
            continue
        k, v = line.split("=", 1)
        os.environ.setdefault(k.strip(), v.strip().strip('"').strip("'"))


@dataclass(frozen=True)
class LLMConfig:
    provider: str
    base_url: str
    api_key: str
    model: str
    group_id: str | None = None

    @classmethod
    def from_env(cls, *, dotenv: bool = True) -> "LLMConfig":
        if dotenv:
            load_dotenv()
        key = os.environ.get("LLM_API_KEY", "").strip()
        if not key or key.startswith("put-your"):
            raise RuntimeError(
                "LLM_API_KEY is not set. Copy .env.example to .env and fill in your "
                "provider key (see .env.example). `.env` is gitignored.")
        return cls(
            provider=os.environ.get("LLM_PROVIDER", "minimax").strip(),
            base_url=os.environ.get("LLM_BASE_URL", "https://api.minimax.io/v1").strip(),
            api_key=key,
            model=os.environ.get("LLM_MODEL", "MiniMax-M3").strip(),
            group_id=(os.environ.get("LLM_GROUP_ID") or None))

    def client(self):
        """An OpenAI-compatible client pointed at this provider. Requires `openai`
        (`uv pip install openai`); imported lazily so config loading needs no SDK."""
        try:
            from openai import OpenAI
        except ImportError as e:                              # noqa: BLE001
            raise RuntimeError("the `openai` client is not installed "
                               "(`uv pip install openai`)") from e
        return OpenAI(api_key=self.api_key, base_url=self.base_url)

    def redacted(self) -> dict:
        """Safe-to-log view -- never exposes the key."""
        tail = self.api_key[-4:] if len(self.api_key) >= 4 else "?"
        return {"provider": self.provider, "base_url": self.base_url,
                "model": self.model, "api_key": f"…{tail}",
                "group_id": self.group_id}

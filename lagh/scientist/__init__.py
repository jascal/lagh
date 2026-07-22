"""The lagh 'scientist' harness: the LLM that orchestrates recover/verify/fit.

Provider-agnostic (OpenAI-compatible), so MiniMax M3 / Hermes / etc. are one config
away. Credentials live in a gitignored `.env` at the repo root (see `.env.example`).
"""

from .config import LLMConfig, load_dotenv

__all__ = ["LLMConfig", "load_dotenv"]

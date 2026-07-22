"""The LLM config loads provider-agnostically from env and never leaks the key."""

import pytest

from lagh.scientist import LLMConfig


def test_config_from_env(monkeypatch):
    monkeypatch.setenv("LLM_PROVIDER", "minimax")
    monkeypatch.setenv("LLM_BASE_URL", "https://api.minimax.io/v1")
    monkeypatch.setenv("LLM_API_KEY", "sk-secret-1234")
    monkeypatch.setenv("LLM_MODEL", "MiniMax-M3")
    cfg = LLMConfig.from_env(dotenv=False)
    assert cfg.provider == "minimax" and cfg.model == "MiniMax-M3"
    assert cfg.api_key == "sk-secret-1234"
    # redacted view must NOT contain the key
    r = cfg.redacted()
    assert "sk-secret-1234" not in str(r) and r["api_key"].endswith("1234")


def test_config_errors_clearly_when_key_missing(monkeypatch):
    monkeypatch.delenv("LLM_API_KEY", raising=False)
    with pytest.raises(RuntimeError, match="LLM_API_KEY is not set"):
        LLMConfig.from_env(dotenv=False)


def test_config_rejects_the_placeholder(monkeypatch):
    monkeypatch.setenv("LLM_API_KEY", "put-your-minimax-key-here")
    with pytest.raises(RuntimeError, match="LLM_API_KEY is not set"):
        LLMConfig.from_env(dotenv=False)

"""Provider-agnostic LLM proposer for the LawResearch machine's `proposing` state.

The LLM is a BOUNDED proposer: given the samples + lagh's (uncertified) characterization,
it proposes ONE sympy form in x_0.. for lagh to `verify` -- or NONE. A wrong guess is simply
refuted by the sound checker, so the proposer is told to give its single best structural
guess. This is the ONLY place the composite calls a model.

Reads the SAME LLM_* env as deploy/hermes (LLM_MODEL / LLM_BASE_URL / LLM_API_KEY /
LLM_API_MODE); loads machine/.env if present. Supports the OpenAI chat dialect
(Grok/GPT, LLM_API_MODE empty) and the Anthropic messages dialect (MiniMax /anthropic,
LLM_API_MODE=anthropic_messages). If nothing is configured, returns None -> the machine
routes NO_FORM and the loop runs on lagh alone (never worse than the tool).
"""
from __future__ import annotations

import os
from pathlib import Path

import numpy as np


def _load_env() -> None:
    try:
        from dotenv import load_dotenv
        p = Path(__file__).with_name(".env")
        if p.exists():
            load_dotenv(p)
    except Exception:                                          # noqa: BLE001
        pass


_PROMPT = """You are proposing candidate closed-form laws to be CHECKED by a sound certifier. You are an UNTRUSTED proposer: a wrong guess is simply rejected, never accepted -- so propose your best STRUCTURALLY DISTINCT guesses, most likely first, and do not hedge.

The inputs are x_0..x_{last}. A certifier could not find an exact law; its uncertified read of the samples:
  class: {cls}
  note:  {why}

Sample rows (x_0..x_{last} -> y):
{rows}

Propose up to {k} sympy expressions in x_0.. that could fit these samples -- for a "non-algebraic" class try trig / inverse-trig forms (sin, cos, asin, atan), for "additive-or-mixed" try sums of monomials, and make the guesses structurally DIFFERENT from each other (not the same shape with different constants). Use ** for powers. Reply with one line per guess and nothing else:
FORM: <a single sympy expression in x_0..>
FORM: <another, structurally different>
(or a single line "FORM: NONE" if you have no plausible guess)
"""


def propose_forms(info: dict, k: int = 3) -> list[str]:
    """The propose_fn the driver injects. `info` = {class, why, X, y, box, problem}.
    Returns up to k structurally-distinct sympy form strings, best first -- possibly
    empty (not configured / no proposal / call failed). Verification is sound, so
    every extra proposal is free upside: a wrong form is refuted, never accepted."""
    _load_env()
    key = os.environ.get("LLM_API_KEY")
    model = os.environ.get("LLM_MODEL")
    base = os.environ.get("LLM_BASE_URL")
    mode = os.environ.get("LLM_API_MODE", "").strip()
    if not (key and model and base):
        return []                                             # unconfigured -> lagh alone
    X = np.asarray(info.get("X", []), float)
    y = np.asarray(info.get("y", []), float)
    if X.ndim != 2 or len(X) == 0:
        return []
    last = X.shape[1] - 1
    rows = "\n".join(
        "  " + ", ".join(f"{v:.4g}" for v in X[i]) + f" -> {y[i]:.4g}"
        for i in range(min(12, len(X))))
    prompt = _PROMPT.format(last=last, k=k, cls=info.get("class", "unresolved"),
                            why=info.get("why", ""), rows=rows)
    try:
        text = _call(base, key, model, mode, prompt)
    except Exception:                                         # noqa: BLE001
        return []                                             # any API error -> no proposal
    forms: list[str] = []
    for ln in text.splitlines():
        if ln.strip().upper().startswith("FORM:"):
            f = ln.split(":", 1)[1].strip().strip("`")
            if f and f.upper() != "NONE" and f not in forms:
                forms.append(f)
    return forms[:k]


def propose_form(info: dict) -> str | None:
    """Single-form compatibility wrapper over propose_forms."""
    forms = propose_forms(info, k=1)
    return forms[0] if forms else None


def _call(base: str, key: str, model: str, mode: str, prompt: str, max_tokens: int = 300) -> str:
    if mode == "anthropic_messages":
        import requests
        url = base.rstrip("/") + "/v1/messages"
        r = requests.post(url, timeout=60, headers={
            "x-api-key": key, "anthropic-version": "2023-06-01", "content-type": "application/json"},
            json={"model": model, "max_tokens": max_tokens,
                  "messages": [{"role": "user", "content": prompt}]})
        r.raise_for_status()
        return "".join(b.get("text", "") for b in r.json().get("content", []))
    from openai import OpenAI                                 # OpenAI chat dialect (Grok/GPT)
    client = OpenAI(base_url=base, api_key=key)
    resp = client.chat.completions.create(
        model=model, max_tokens=max_tokens,
        messages=[{"role": "user", "content": prompt}])
    return resp.choices[0].message.content or ""

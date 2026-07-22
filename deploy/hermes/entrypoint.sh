#!/usr/bin/env bash
# Runtime setup: install lagh into its venv (editable, from the mounted source), then
# write the Hermes config (MiniMax model + lagh MCP server) if it isn't there yet.
# Idempotent -- safe to re-run on every container start.
set -euo pipefail

LAGH_SRC="${LAGH_SRC:-/opt/lagh}"
LAGH_PY="/opt/lagh-venv/bin/python"
HERMES_CFG="${HERMES_CFG:-/root/.hermes/config.yaml}"

# 1) lagh (+ the mcp extra) into its isolated venv, editable so host edits reflect live
if [ -f "${LAGH_SRC}/pyproject.toml" ]; then
    # reinstall if either the mcp OR the lab subpackage is missing (a warm venv from an
    # older image may predate lagh.lab)
    if ! "${LAGH_PY}" -c "import lagh.mcp, lagh.lab" 2>/dev/null; then
        echo "[entrypoint] installing lagh[mcp] (editable) from ${LAGH_SRC} ..."
        "${LAGH_PY}" -m pip install -q --upgrade pip
        "${LAGH_PY}" -m pip install -q -e "${LAGH_SRC}[mcp]"
    fi
    # the `newtonbench` lab source imports NewtonBench's eval chain (openai, requests, ...)
    # at module load; install its deps when NewtonBench is mounted so the source can load
    NB_DIR="${NEWTONBENCH_DIR:-/opt/NewtonBench}"
    if [ -d "${NB_DIR}" ] && ! "${LAGH_PY}" -c "import openai, requests" 2>/dev/null; then
        echo "[entrypoint] installing newtonbench source deps (openai, requests, its requirements) ..."
        [ -f "${NB_DIR}/requirements.txt" ] && "${LAGH_PY}" -m pip install -q -r "${NB_DIR}/requirements.txt" || true
        "${LAGH_PY}" -m pip install -q openai requests
    fi
    echo "[entrypoint] lagh ready: $("${LAGH_PY}" -c 'import lagh; print(lagh.__version__)')"
else
    echo "[entrypoint] WARNING: no lagh source at ${LAGH_SRC} (mount it in docker-compose)."
fi

# 2) Hermes config: MiniMax as a custom OpenAI-compatible provider + lagh MCP server.
#    The API key is NOT written here -- `key_env` tells Hermes to read it from the
#    container env (compose injects LLM_API_KEY from your .env).
mkdir -p "$(dirname "${HERMES_CFG}")"
# ALWAYS regenerate from the current env, so a `.env` change + `up -d --build` takes
# effect (the old skip-if-exists guard cached a stale config in the volume). If you
# hand-edit the config, set HERMES_KEEP_CONFIG=1 to preserve it.
if [ "${HERMES_KEEP_CONFIG:-0}" != "1" ] || [ ! -f "${HERMES_CFG}" ]; then
    # provider-agnostic config gen (MiniMax anthropic / Grok / GPT openai) -- one place,
    # so swapping the backing model is just an .env change.
    LAGH_PY="${LAGH_PY}" HERMES_CFG="${HERMES_CFG}" \
        "${LAGH_PY}" "${LAGH_SRC}/deploy/hermes/gen_config.py"
else
    echo "[entrypoint] HERMES_KEEP_CONFIG=1: keeping existing config (${HERMES_CFG})."
fi

if [ -z "${LLM_API_KEY:-}" ]; then
    echo "[entrypoint] NOTE: LLM_API_KEY is empty -- set it in deploy/hermes/.env before running Hermes."
fi

# install the repo's Hermes skills (auto-discovered from ~/.hermes/skills on startup).
# The law-discovery skill gives weaker models the explicit sampling/delegation recipe;
# it doesn't constrain strong models.
if [ -d "${LAGH_SRC}/deploy/hermes/skills" ]; then
    mkdir -p /root/.hermes/skills
    cp -r "${LAGH_SRC}/deploy/hermes/skills/." /root/.hermes/skills/
    echo "[entrypoint] skills installed: $(ls /root/.hermes/skills 2>/dev/null | tr '\n' ' ')"
fi

echo "[entrypoint] ready. Run:  docker compose exec -it hermes hermes"
exec "$@"

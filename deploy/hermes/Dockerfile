# Hermes Agent + lagh in one container.
#   - Hermes Agent (NousResearch) is the MCP host + reasoning loop, backed by MiniMax.
#   - lagh is installed in its own venv and registered as an MCP server (recover/verify/fit).
# UNTESTED by the author (built on the "you run it" path) -- the likely adjustment points
# are called out in deploy/hermes/README.md. Base is Debian-slim Python 3.12.
FROM python:3.12-slim

ENV DEBIAN_FRONTEND=noninteractive PIP_DISABLE_PIP_VERSION_CHECK=1

# System deps: git/curl for installers, ripgrep+ffmpeg (Hermes tools), node+npm for
# npx-launched MCP servers, build-essential for any wheels that need compiling.
RUN apt-get update && apt-get install -y --no-install-recommends \
        git curl ca-certificates build-essential ripgrep ffmpeg nodejs npm \
    && rm -rf /var/lib/apt/lists/*

# uv (Hermes uses it for its own venv)
ADD https://astral.sh/uv/install.sh /tmp/uv-install.sh
RUN sh /tmp/uv-install.sh && rm /tmp/uv-install.sh
ENV PATH="/root/.local/bin:${PATH}"

# --- Hermes Agent (non-interactive: skip the setup wizard and the browser toolchain) ---
RUN curl -fsSL https://hermes-agent.nousresearch.com/install.sh -o /tmp/hermes-install.sh \
    && bash /tmp/hermes-install.sh --skip-setup --skip-browser \
    && rm /tmp/hermes-install.sh
ENV PATH="/root/.hermes/bin:${PATH}"

# --- lagh in its OWN venv, isolated from Hermes' deps ---
# The lagh source is mounted at /opt/lagh at runtime; the entrypoint installs it editable
# into this venv on first boot so `/opt/lagh-venv/bin/python -m lagh.mcp.server` is the
# command Hermes spawns. PRE-INSTALL the deps here (baked into the image) so a fresh venv
# on every recreate/model-swap already has them -- the entrypoint then only re-links lagh
# (fast) instead of reinstalling scipy/sklearn/pandas each time. Includes NewtonBench's
# eval-chain deps (openai/requests/scipy/sklearn/pandas/...) so the newtonbench lab source
# imports without a runtime install.
RUN python -m venv /opt/lagh-venv \
    && /opt/lagh-venv/bin/pip install --no-cache-dir -q --upgrade pip \
    && /opt/lagh-venv/bin/pip install --no-cache-dir -q \
       numpy sympy "mcp>=1.0" \
       openai requests scipy scikit-learn pandas python-dotenv fix-busted-json

COPY entrypoint.sh /usr/local/bin/entrypoint.sh
RUN chmod +x /usr/local/bin/entrypoint.sh

WORKDIR /work
ENTRYPOINT ["/usr/local/bin/entrypoint.sh"]
# keep the container alive; interact with:  docker compose exec -it hermes hermes
CMD ["sleep", "infinity"]

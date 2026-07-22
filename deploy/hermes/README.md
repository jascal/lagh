# Hermes Agent + lagh in Docker

Runs **Hermes Agent** (Nous Research — the MCP host + reasoning loop) backed by **MiniMax**,
with **lagh** wired in as an MCP server (`recover` / `verify` / `fit`). One container, clean
start/stop, nothing installed on your host.

> **Authored but not run by me** — I don't have Docker access in this session, so this is the
> "you run it" path. It's built to the documented Hermes config formats, but expect one or two
> tweaks; the likely spots are flagged under **If something's off** below. Paste me any error
> and I'll fix it.

## Prerequisites
- Docker with access for your user (the daemon is up on this box, but your login isn't in the
  `docker` group — run these with `sudo docker …`, or `sudo usermod -aG docker $USER` then
  re-login to drop the `sudo`).

## Setup & run
```bash
cd deploy/hermes
cp .env.example .env
$EDITOR .env                     # set LLM_API_KEY, LLM_MODEL, LLM_BASE_URL (MiniMax)
docker compose up -d --build     # build image, start container (prefix sudo if needed)
```

## Use it
```bash
docker compose exec -it hermes hermes            # interactive agent; lagh's tools are available
docker compose exec -it hermes hermes mcp test lagh   # verify lagh connected
```
Inside the agent, lagh appears as three tools:
- `recover(X, y[, sigma])` — discover an exact law → certificate or reasoned abstention
- `verify(X, y, form[, sigma])` — check a declared form → certificate (`pinned`/`consistent`) or abstain
- `fit(X, y[, sigma])` — scout: conjectures + identifiability diagnosis (NOT a certificate)

## Stop / restart
```bash
docker compose down      # stop; config + sessions persist in the `hermes-home` volume
docker compose up -d     # back up, state intact
```

## What's wired
- `Dockerfile` — Python 3.12 + uv + Hermes (installed `--skip-setup --skip-browser`) + a
  separate `/opt/lagh-venv`.
- `entrypoint.sh` — on boot, `pip install -e /opt/lagh[mcp]` into that venv, then writes
  `~/.hermes/config.yaml` (MiniMax as a `provider: custom` model + lagh under `mcp_servers`,
  key read from env via `key_env`). Idempotent.
- `docker-compose.yml` — mounts the lagh repo at `/opt/lagh` (live edits reflect) and a named
  volume for `~/.hermes`; injects the MiniMax creds from `.env`.

## If something's off (the untested seams)
1. **Config path.** Docs show both `~/.hermes/config.yaml` and `~/.hermes/hermes-agent/config.yaml`.
   If Hermes ignores the config, set `HERMES_CFG` in `docker-compose.yml` to the other path and
   `docker compose up -d` again (or run `hermes config path` inside the container to find it).
2. **Model provider form.** If `provider: custom` isn't picked up, Hermes also supports a named
   `custom_providers:` block selected via `/model custom:<name>:<model-id>`. Or configure it
   interactively once: `docker compose exec -it hermes hermes model`.
3. **MiniMax base_url / model id.** Set to your region's endpoint (ends in `/v1`; Hermes appends
   `/chat/completions`) and the exact model string from your console.
4. **`key_env` support.** If your Hermes build wants the key inline instead, put it in
   `~/.hermes/.env` inside the volume, or set `api_key:` directly (less clean).
5. **lagh MCP command.** The entrypoint points Hermes at `/opt/lagh-venv/bin/python -m
   lagh.mcp.server`. Confirm with `docker compose exec -it hermes /opt/lagh-venv/bin/python -m
   lagh.mcp.server </dev/null` (should start, then block waiting for a client).

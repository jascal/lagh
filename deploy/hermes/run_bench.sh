#!/usr/bin/env bash
# Full-run benchmark: drives MiniMax+Hermes+lagh over NewtonBench-dev cells and scores.
# Runs the driver INSIDE the container (one exec, not per-cell), so it needs no sudo per
# cell and scores in-process against the real oracle. Results land in bench_results/.
#
#   ./run_bench.sh --subset easy            # all 36 easy cells
#   ./run_bench.sh --subset easy --limit 6  # a quick first pass
#   ./run_bench.sh --cells m0_gravity/easy/v0,m4_snell_law/easy/v0
#
# Requires the container up with LAB_SOURCE=newtonbench (see README). Defaults to
# `sudo docker`; override with DOCKER=docker if you're in the docker group.
set -euo pipefail
cd "$(dirname "$0")"
DOCKER="${DOCKER:-sudo docker}"

# sanity: the lab must be on the newtonbench source (else you'd score proxy problems)
if ! $DOCKER exec lagh-hermes printenv LAB_SOURCE 2>/dev/null | grep -qx newtonbench; then
  echo "WARNING: container LAB_SOURCE is not 'newtonbench'." >&2
  echo "  set LAB_SOURCE=newtonbench in .env and: $DOCKER compose up -d --build" >&2
  exit 1
fi

# DETERMINISTIC skill environment: refresh our law-discovery skill from the repo and drop
# any Hermes `/learn`-captured loop (e.g. law-discovery-loop) that would route the model
# back into slow manual sampling. Reproducibility -- the score must not depend on whatever
# got auto-learned. Targeted: touches only these, not the bundled skills.
echo "[run_bench] resetting skill env (fresh law-discovery, drop learned loops)" >&2
$DOCKER exec lagh-hermes bash -lc '
  rm -rf /root/.hermes/skills/law-discovery-loop /root/.hermes/skills/*law-discovery*loop* 2>/dev/null
  mkdir -p /root/.hermes/skills
  cp -r /opt/lagh/deploy/hermes/skills/. /root/.hermes/skills/' >/dev/null 2>&1 || true

exec $DOCKER exec -it lagh-hermes \
  /opt/lagh-venv/bin/python /opt/lagh/deploy/hermes/bench.py "$@"

#!/usr/bin/env bash
set -euo pipefail
MODE="${1:-perf}"
CONFIG="${2:-configs/nonint3body_general_doublecover_config.in}"
mkdir -p output bin
export PYTHONUNBUFFERED=1
printf '[run] started at %s\n' "$(date)"
printf '[run] MODE=%s CONFIG=%s\n' "$MODE" "$CONFIG"
printf '[run] [--------------------] 0%% starting python executable\n'
if [[ "$MODE" == "smoke" ]]; then
  TMP="output/nonint3body_general_doublecover_smoke_config.in"
  cp "$CONFIG" "$TMP"
  sed -i 's/^Ecm_max.*/Ecm_max = 0.33/' "$TMP" || true
  sed -i 's/^nmax.*/nmax = 3/' "$TMP" || true
  python3 source/nonint3body_general_doublecover.py "$TMP" 2>&1 | tee output/nonint3body_general_doublecover_smoke.log
else
  python3 source/nonint3body_general_doublecover.py "$CONFIG" 2>&1 | tee output/nonint3body_general_doublecover_run.log
fi
printf '[run] [####################] 100%% complete\n'

#!/usr/bin/env bash
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
# CPU runs use the package-local INI-style config by default.
CONFIG="${1:-$ROOT/configs/config_cpu_component_cachegen.in}"
BIN="$ROOT/bin/generate_cpu_component_cache"
if [[ ! -x "$BIN" ]]; then
  # Build on demand so the package is runnable immediately after unzip.
  echo "[run-cpu] executable missing; compiling first"
  "$ROOT/scripts/compile_cpu_component_cachegen.sh"
fi
mkdir -p "$ROOT/logs"
# The CPU path emits a single consolidated log.
"$BIN" "$CONFIG" 2>&1 | tee "$ROOT/logs/run_cpu_component_cachegen.log"

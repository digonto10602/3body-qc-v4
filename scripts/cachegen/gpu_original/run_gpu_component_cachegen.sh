#!/usr/bin/env bash
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
# Default to the package-local config and binary. The config is sourced so its
# shell variables become the runtime settings for the generator.
CONFIG="${V32ZU_CONFIG:-$ROOT/configs/config_gpu_component_cachegen.sh}"
BIN="${V32ZU_BIN:-$ROOT/bin/cachegen_gpu_v32zu_allobj_qc}"

if [[ ! -f "$CONFIG" ]]; then
  echo "[run ERROR] config not found: $CONFIG" >&2
  exit 1
fi
# shellcheck source=/dev/null
source "$CONFIG"

usage() {
  cat >&2 <<USAGE
Usage:
  bash scripts/run_gpu_component_cachegen.sh <Lbyas>
  bash scripts/run_gpu_component_cachegen.sh <Lbyas_min> <Lbyas_max>
  bash scripts/run_gpu_component_cachegen.sh <Lbyas_min> <Lbyas_max> <step_size>

Examples:
  bash scripts/run_gpu_component_cachegen.sh 20
  bash scripts/run_gpu_component_cachegen.sh 19 22.4
  bash scripts/run_gpu_component_cachegen.sh 19 22.4 0.2

Notes:
  - 1 argument  => single-L run.
  - 2 arguments => inclusive L range using default_step_size from config, default 0.2.
  - 3 arguments => inclusive L range using the provided step_size.
  - 0 or 4+ arguments are invalid.
USAGE
}

if [[ $# -eq 0 || $# -gt 3 ]]; then
  usage
  exit 2
fi

if [[ ! -x "$BIN" ]]; then
  # Keep the package runnable from a fresh unzip: build the binary on demand
  # instead of requiring the user to remember a separate compile step.
  echo "[run-gpu-components] executable missing; compiling component-cache binary first..."
  V32ZU_ALLOBJ_QC=1 bash "$ROOT/scripts/compile_gpu_component_cachegen.sh"
fi

# Create the configured output root before launching any runs.
mkdir -p "$output_root"

# Build the Lbyas list without floating-point drift.
if [[ $# -eq 1 ]]; then
  RUNMODE="single-L mode"
  LVALUES=("$1")
else
  RUNMODE="L-range mode"
  LMIN="$1"
  LMAX="$2"
  STEP="${3:-${default_step_size:-0.2}}"
  # Python Decimal is the simplest safe way to step through lattice sizes
  # like 19, 19.2, 19.4 without floating-point accumulation errors.
  mapfile -t LVALUES < <(python3 - "$LMIN" "$LMAX" "$STEP" <<'PY'
import sys, decimal
D = decimal.Decimal
lm = D(sys.argv[1]); lx = D(sys.argv[2]); st = D(sys.argv[3])
if st <= 0:
    raise SystemExit("step_size must be positive")
if lx < lm:
    raise SystemExit("Lbyas_max must be >= Lbyas_min")
# scale enough for decimal inputs like 0.2, 0.05, etc.
vals = []
x = lm
# include max with tiny tolerance by integer-like decimal loop
while x <= lx + D('1e-30'):
    vals.append(x)
    x = x + st
for v in vals:
    s = format(v, 'f')
    if '.' in s:
        s = s.rstrip('0').rstrip('.')
    print(s)
PY
  )
fi

tagL() {
  python3 - "$1" <<'PY'
import sys
s = sys.argv[1]
if '.' in s:
    s = s.rstrip('0').rstrip('.')
print(s.replace('.', 'p'))
PY
}

format_seconds() {
  python3 - "$1" <<'PY'
import sys
s=int(round(float(sys.argv[1])))
print(f"{s//3600}:{(s%3600)//60:02d}:{s%60:02d}")
PY
}

RANGE_START=$(date +%s)
echo "============================================================"
echo "[RUNMODE] $RUNMODE"
echo "[CONFIG] $CONFIG"
# The next lines are the values that the sourced config injected.
echo "[CONFIG] Ecm_min=$Ecm_min Ecm_max=$Ecm_max coarseN=$coarseN"
echo "[CONFIG] output_root=$output_root"
echo "[CONFIG] use_fused_f2k2=${use_fused_f2k2:-true} validate_fused_f2k2=${validate_fused_f2k2:-false}"
echo "[CONFIG] use_kinematic_precompute=${use_kinematic_precompute:-true} use_stream_pipeline=${use_stream_pipeline:-true}"
echo "[CONFIG] max_batch_energies=${max_batch_energies:-0} gpu_memory_safety_fraction=${gpu_memory_safety_fraction:-0.85}"
echo "[CONFIG] use_multi_gpu=${use_multi_gpu:-false} use_mixed_precision_f2k2=${use_mixed_precision_f2k2:-false}"
echo "[CONFIG] use_concurrent_streams=${use_concurrent_streams:-true} concurrent_streams_safety_fraction=${concurrent_streams_safety_fraction:-0.85}"
echo "[CONFIG] irreps=${irreps:-000_A1m,100_A2,110_A2,111_A2,200_A2}"
echo "[LGRID] values: ${LVALUES[*]}"
echo "[LGRID] total L values=${#LVALUES[@]}"
echo "============================================================"

idx=0
for L in "${LVALUES[@]}"; do
  idx=$((idx+1))
  LTAG=$(tagL "$L")
  # Every L gets its own cache directory, so a partial run does not overwrite
  # a different ensemble or irrep sweep.
  OUTDIR="$output_root/v32zu_Lbyas${LTAG}_gpu_cache"
  mkdir -p "$OUTDIR"
  echo "============================================================"
  echo "[L ${idx}/${#LVALUES[@]}] Lbyas=$L"
  echo "[OUTPUT] $OUTDIR"
  echo "[EXEC] $BIN $L $xi $Ecm_min $Ecm_max $coarseN $OUTDIR ${debug:-n} ${max_batch_energies:-0} ${gpu_memory_safety_fraction:-0.85} ${use_fused_f2k2:-true} ${validate_fused_f2k2:-false} ${use_kinematic_precompute:-true} ${use_stream_pipeline:-true} ${use_multi_gpu:-false} ${use_mixed_precision_f2k2:-false} ${use_concurrent_streams:-true} ${concurrent_streams_safety_fraction:-0.85} ${irreps:-000_A1m,100_A2,110_A2,111_A2,200_A2}"
  echo "============================================================"
  LSTART=$(date +%s)
  # The executable writes the component caches; tee preserves a per-L log.
  "$BIN" "$L" "$xi" "$Ecm_min" "$Ecm_max" "$coarseN" "$OUTDIR" "${debug:-n}" "${max_batch_energies:-0}" "${gpu_memory_safety_fraction:-0.85}" "${use_fused_f2k2:-true}" "${validate_fused_f2k2:-false}" "${use_kinematic_precompute:-true}" "${use_stream_pipeline:-true}" "${use_multi_gpu:-false}" "${use_mixed_precision_f2k2:-false}" "${use_concurrent_streams:-true}" "${concurrent_streams_safety_fraction:-0.85}" "${irreps:-000_A1m,100_A2,110_A2,111_A2,200_A2}" | tee "$OUTDIR/run_v32zu_Lbyas${LTAG}.log"
  LEND=$(date +%s)
  LSEC=$((LEND-LSTART))
  echo "[RUN SCRIPT L SUMMARY] Lbyas=$L wall_time_sec=$LSEC wall_time_hms=$(format_seconds "$LSEC") output=$OUTDIR"
done

RANGE_END=$(date +%s)
RANGE_SEC=$((RANGE_END-RANGE_START))
echo "============================================================"
echo "[RUN SCRIPT RANGE SUMMARY]"
echo "Lbyas values completed = ${#LVALUES[@]}"
IFS=',' read -r -a IRREP_ARRAY <<< "${irreps:-000_A1m,100_A2,110_A2,111_A2,200_A2}"
# This summary makes the run self-documenting for postmortem analysis.
echo "irreps per Lbyas       = ${#IRREP_ARRAY[@]}"
echo "coarseN                = $coarseN"
echo "total expected rows    = $((${#LVALUES[@]} * ${#IRREP_ARRAY[@]} * coarseN))"
echo "total range wall time  = ${RANGE_SEC} sec ($(format_seconds "$RANGE_SEC"))"
echo "output root            = $output_root"
echo "============================================================"

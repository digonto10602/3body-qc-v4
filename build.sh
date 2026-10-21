#!/usr/bin/env bash
# Build the CPU tools of 3body-qc-v4.
# Compile flags are taken unchanged from the original packages:
#   fitter + hot-window builder : KKpi-L20-user-validated-hotwindow-fitter/build.sh
#   v33h scanner, v33g builder   : KKpi-L20-.../scripts/compile_v33f_all.sh, compile_v33g_all.sh
#   CPU cache generator          : full_cache_generator_package/cpu_component_cachegen/scripts/compile_cpu_component_cachegen.sh
# Only the include paths changed (sources now live in src/<module>/).
# NOTE: this script has not yet been run in the new layout (see notes/STATUS.md).
#
# Env: CXX (default g++), EIGEN_INC (default /usr/include/eigen3), MINUIT2_ROOT (optional prefix).
# Usage: ./build.sh [all|core|pyext|fitter|hotwindow|scan|runtime|cachegen_cpu]
#   all    = hotwindow scan runtime cachegen_cpu core; the Minuit2 fitter is added only if Minuit2 headers are found
#   core   = pybind11 module qc3._core from src/qc3_core -> src/python/qc3/_core*.so (pyext is an alias)
#   fitter = on demand only (needs Minuit2)
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CXX="${CXX:-g++}"
EIGEN_INC="${EIGEN_INC:-/usr/include/eigen3}"
MINUIT2_ROOT="${MINUIT2_ROOT:-${CONDA_PREFIX:-}}"
TARGET="${1:-all}"

[[ -f "$EIGEN_INC/Eigen/Dense" ]] || { echo "[build] Eigen not found at EIGEN_INC=$EIGEN_INC" >&2; exit 2; }

INC=(-I"$EIGEN_INC" -I"$ROOT/src/qc3_core" -I"$ROOT/src/core" -I"$ROOT/src/fitter" -I"$ROOT/src/cache_io" -I"$ROOT/src/cache_build")
COMMON=(-std=c++17 -O3 -fopenmp)
NO_MINUIT=(-DV32F_HAS_MINUIT2=0 -DV32F_DISABLE_MINUIT)
MINUIT_INC=(); MINUIT_LIB=(-lMinuit2 -lMinuit2Math)
if [[ -n "$MINUIT2_ROOT" && -d "$MINUIT2_ROOT/include/Minuit2" ]]; then
  MINUIT_INC=(-I"$MINUIT2_ROOT/include" -I"$MINUIT2_ROOT/include/Minuit2")
  MINUIT_LIB=(-L"$MINUIT2_ROOT/lib" -Wl,-rpath,"$MINUIT2_ROOT/lib" -lMinuit2 -lMinuit2Math)
fi
FAD="$ROOT/src/core/Faddeeva.cc"
mkdir -p "$ROOT/bin"

build_fitter() {   # multi-L K3df fitter (needs Minuit2)
  "$CXX" "${COMMON[@]}" "${INC[@]}" "${MINUIT_INC[@]}" \
    "$ROOT/src/fitter/qc_fitter_norm_refine_v2_multiL.cpp" "$FAD" "${MINUIT_LIB[@]}" \
    -o "$ROOT/bin/k3df_fitter_multiL"
}
build_hotwindow() { # projected-basis hot-window cache builder (reads raw F3inv/Vsel cache)
  "$CXX" "${COMMON[@]}" "${NO_MINUIT[@]}" "${INC[@]}" \
    "$ROOT/src/cache_build/v33p_build_projected_basis_hotwindow_cache.cpp" "$FAD" \
    -o "$ROOT/bin/build_projected_basis_hotwindow_cache"
}
build_scan() {      # trusted CPU/OpenMP determinant scan of the raw cache (v33h)
  "$CXX" "${COMMON[@]}" "${NO_MINUIT[@]}" "${INC[@]}" \
    "$ROOT/src/cache_io/v33h_patched_gpu_cache_oldscale_det_scan.cpp" "$FAD" \
    -o "$ROOT/bin/det_scan_raw_cache_v33h"
}
build_runtime() {   # full-grid runtime K3-basis cache builder (v33g, v33n copy)
  "$CXX" "${COMMON[@]}" "${NO_MINUIT[@]}" "${INC[@]}" \
    "$ROOT/src/cache_build/v33g_build_runtime_k3basis_cache.cpp" "$FAD" \
    -o "$ROOT/bin/build_runtime_k3basis_cache_v33g"
}
build_cachegen_cpu() { # CPU reference component-cache generator (v33d)
  "$CXX" "${COMMON[@]}" "${NO_MINUIT[@]}" -DV31L_DISABLE_MINUIT=1 -DV31R_DISABLE_MINUIT=1 \
    -DV31S_DISABLE_MINUIT=1 -DV31T_DISABLE_MINUIT=1 "${INC[@]}" \
    "$ROOT/src/cachegen_cpu/generate_cpu_component_cache.cpp" "$FAD" \
    -o "$ROOT/bin/generate_cpu_component_cache"
}

build_core() {      # pybind11 module qc3._core (prompts/001 step 1)
  python3 -c "import pybind11" 2>/dev/null || { echo "[build] pybind11 missing: pip install pybind11" >&2; exit 2; }
  local PYINC EXT GIT
  PYINC="$(python3 -m pybind11 --includes)"
  EXT="$(python3 -c 'import sysconfig; print(sysconfig.get_config_var("EXT_SUFFIX"))')"
  GIT="$(git -C "$ROOT" rev-parse --short HEAD 2>/dev/null || echo unknown)"
  # shellcheck disable=SC2086
  "$CXX" "${COMMON[@]}" "${NO_MINUIT[@]}" -DV31L_DISABLE_MINUIT=1 -DV31R_DISABLE_MINUIT=1 \
    -DV31S_DISABLE_MINUIT=1 -DV31T_DISABLE_MINUIT=1 -shared -fPIC -DQC3_GIT_COMMIT="\"$GIT\"" "${INC[@]}" $PYINC \
    "$ROOT/src/qc3_core/module.cpp" "$FAD" \
    -o "$ROOT/src/python/qc3/_core$EXT"
}
have_minuit2() { [[ -n "$MINUIT2_ROOT" && -d "$MINUIT2_ROOT/include/Minuit2" ]] || [[ -f /usr/include/Minuit2/FCNBase.h ]]; }

case "$TARGET" in
  all) build_hotwindow; build_scan; build_runtime; build_cachegen_cpu; build_core
       if have_minuit2; then build_fitter; else echo "[build] Minuit2 not found: skipping fitter (use ./build.sh fitter)"; fi ;;
  core|pyext) build_core; echo "[build] done: core -> $ROOT/src/python/qc3/"; exit 0 ;;
  fitter) build_fitter ;; hotwindow) build_hotwindow ;; scan) build_scan ;;
  runtime) build_runtime ;; cachegen_cpu) build_cachegen_cpu ;;
  *) echo "unknown target $TARGET" >&2; exit 2 ;;
esac
echo "[build] done: $TARGET -> $ROOT/bin/"

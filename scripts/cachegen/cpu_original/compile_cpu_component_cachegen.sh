#!/usr/bin/env bash
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
PKG_ROOT="$(cd "$ROOT/.." && pwd)"
# Keep the output directories local to the package.
mkdir -p "$ROOT/bin" "$ROOT/logs" "$ROOT/output"
CXX="${CXX:-g++}"
if [[ -n "${EIGEN_INC:-}" ]]; then
  :
elif [[ -n "${CONDA_PREFIX:-}" && -d "$CONDA_PREFIX/include/eigen3/Eigen" ]]; then
  EIGEN_INC="$CONDA_PREFIX/include/eigen3"
elif [[ -d /usr/include/eigen3/Eigen ]]; then
  EIGEN_INC="/usr/include/eigen3"
elif [[ -d /usr/local/include/eigen3/Eigen ]]; then
  EIGEN_INC="/usr/local/include/eigen3"
else
  echo "[compile-error] Could not find Eigen/Dense. Set EIGEN_INC=/path/to/eigen3" >&2
  exit 2
fi
COMMON_FLAGS=(-O3 -std=c++17 -fopenmp -DV32F_DISABLE_MINUIT=1 -DV31L_DISABLE_MINUIT=1 -DV31R_DISABLE_MINUIT=1 -DV31S_DISABLE_MINUIT=1 -DV31T_DISABLE_MINUIT=1)
INCLUDES=(-I"$ROOT/source" -I"$PKG_ROOT/legacy_source" -I"$EIGEN_INC")
echo "[compile-cpu] CXX=$CXX"
echo "[compile-cpu] EIGEN_INC=$EIGEN_INC"
echo "[compile-cpu] PKG_ROOT=$PKG_ROOT"
set -x
# The CPU generator is a single translation unit plus the Faddeeva helper.
"$CXX" "${COMMON_FLAGS[@]}" "${INCLUDES[@]}" \
  "$ROOT/source/generate_cpu_component_cache.cpp" "$PKG_ROOT/legacy_source/Faddeeva.cc" \
  -o "$ROOT/bin/generate_cpu_component_cache"
set +x
echo "[compile-cpu] done: $ROOT/bin/generate_cpu_component_cache"

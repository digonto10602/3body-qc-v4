#!/usr/bin/env bash
# Build the GPU (CUDA) tools of 3body-qc-v4.
# The cache-generator command is copied from
#   full_cache_generator_package_20260629/gpu_component_cachegen/scripts/compile_gpu_component_cachegen.sh
# with include paths changed to the new layout. The original script is kept verbatim in
# scripts/cachegen/gpu_original/ together with the run script and configs/cachegen_gpu/.
# NOTE: not yet run in the new layout (see notes/STATUS.md).
#
# Env: EIGEN_INC, CUDA_ARCH_FLAGS (default -arch=native), V32ZU_ALLOBJ_QC (default 1 = all-object writer,
#      the production v33d build that writes the combined *_F3inv_Vsel_gpu.bin plus component caches).
# Usage: ./build_gpu.sh [cachegen|detlu]
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
EIGEN_INC="${EIGEN_INC:-/usr/include/eigen3}"
CUDA_ARCH_FLAGS="${CUDA_ARCH_FLAGS:--arch=native}"
ALLOBJ="${V32ZU_ALLOBJ_QC:-1}"
TARGET="${1:-cachegen}"
mkdir -p "$ROOT/bin"

if [[ "$TARGET" == cachegen ]]; then
  OUT="$ROOT/bin/F3_gpu_cachegen_v32zu"; EXTRA=""
  if [[ "$ALLOBJ" == 1 ]]; then OUT="$ROOT/bin/cachegen_gpu_v32zu_allobj_qc"; EXTRA="-DV32ZU_ALLOBJ_QC_WRITER=1"; fi
  nvcc -O3 -std=c++17 -DEIGEN_NO_CUDA --expt-relaxed-constexpr \
    -Xcompiler -fopenmp \
    -Xcompiler -Wno-unused-result -Xcompiler -Wno-unused-variable -Xcompiler -Wno-unused-but-set-variable \
    --diag-suppress=177 --diag-suppress=550 --diag-suppress=1650 --diag-suppress=20208 \
    -DUSE_CONCURRENT_STREAMS=1 $EXTRA $CUDA_ARCH_FLAGS \
    -I"$ROOT/src/cachegen_gpu" -I"$ROOT/src/core" -I"$EIGEN_INC" \
    "$ROOT/src/cachegen_gpu/F3_gpu_cachegen_v32zu.cu" "$ROOT/src/core/Faddeeva.cc" \
    -lcublas -lcusolver -o "$OUT"
  echo "[build_gpu] $OUT"
elif [[ "$TARGET" == detlu ]]; then
  # Experimental GPU determinant scan (v33p). Not production: see docs/OPEN_ISSUES.md.
  # No original compile script survived; this command is a reconstruction.
  nvcc -O3 -std=c++17 $CUDA_ARCH_FLAGS -c "$ROOT/src/gpu_det/v33p_gpu_det_lu.cu" -o "$ROOT/bin/v33p_gpu_det_lu.o"
  g++ -std=c++17 -O3 -fopenmp -DV32F_HAS_MINUIT2=0 -DV32F_DISABLE_MINUIT \
    -I"$EIGEN_INC" -I"$ROOT/src/gpu_det" -I"$ROOT/src/core" -I"$ROOT/src/fitter" -I"$ROOT/src/cache_io" \
    "$ROOT/src/gpu_det/v33p_cpuassembled_gpu_det_scan.cpp" "$ROOT/src/core/Faddeeva.cc" "$ROOT/bin/v33p_gpu_det_lu.o" \
    -lcudart -lcublas -lcusolver -o "$ROOT/bin/v33p_cpuassembled_gpu_det_scan"
  echo "[build_gpu] $ROOT/bin/v33p_cpuassembled_gpu_det_scan"
else
  echo "usage: $0 [cachegen|detlu]" >&2; exit 2
fi

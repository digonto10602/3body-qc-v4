#!/usr/bin/env bash
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
# Source and output live inside the package tree.
SRC="$ROOT/source/F3_gpu_cachegen_v32zu.cu"
OUT="$ROOT/bin/F3_gpu_cachegen_v32zu"
# Keep the compiler include path explicit so the package works on machines
# with either system Eigen or a custom Eigen install.
EIGEN_INC="${EIGEN_INC:-/usr/include/eigen3}"
LEGACY_INC="$ROOT/../legacy_source"
CUDA_ARCH_FLAGS="${CUDA_ARCH_FLAGS:--arch=native}"
NVCC_DIAG_SUPPRESS_FLAGS="${NVCC_DIAG_SUPPRESS_FLAGS:---diag-suppress=177 --diag-suppress=550 --diag-suppress=1650 --diag-suppress=20208}"
HOST_WARNING_SUPPRESS_FLAGS="${HOST_WARNING_SUPPRESS_FLAGS:- -Wno-unused-result -Wno-unused-variable -Wno-unused-but-set-variable}"
MIXED_PRECISION_FLAG="${MIXED_PRECISION_FLAG:-}"
CONCURRENT_STREAMS_FLAG="${CONCURRENT_STREAMS_FLAG:--DUSE_CONCURRENT_STREAMS=1}"
ALLOBJ_FLAG="${V32ZU_ALLOBJ_QC:-0}"
if [[ "$ALLOBJ_FLAG" == "1" || "${ALLOBJ_FLAG,,}" == "true" ]]; then
  # The component-cache build uses the all-object writer binary.
  OUT="$ROOT/bin/cachegen_gpu_v32zu_allobj_qc"
  EXTRA_DEFINES="-DV32ZU_ALLOBJ_QC_WRITER=1"
else
  EXTRA_DEFINES=""
fi
mkdir -p "$ROOT/bin"
echo "[compile] root=$ROOT"
echo "[compile] source=$SRC"
echo "[compile] output=$OUT"
echo "[compile] EIGEN_INC=$EIGEN_INC"
echo "[compile] LEGACY_INC=$LEGACY_INC"
echo "[compile] CUDA_ARCH_FLAGS=$CUDA_ARCH_FLAGS"
echo "[compile] NVCC_DIAG_SUPPRESS_FLAGS=$NVCC_DIAG_SUPPRESS_FLAGS"
echo "[compile] HOST_WARNING_SUPPRESS_FLAGS=$HOST_WARNING_SUPPRESS_FLAGS"
echo "[compile] MIXED_PRECISION_FLAG=$MIXED_PRECISION_FLAG"
echo "[compile] CONCURRENT_STREAMS_FLAG=$CONCURRENT_STREAMS_FLAG"
echo "[compile] ALLOBJ_FLAG=$ALLOBJ_FLAG"
set -x
# The package keeps the exact compile command visible so the build is
# reproducible without guessing flags from surrounding context.
nvcc -O3 -std=c++17 -DEIGEN_NO_CUDA --expt-relaxed-constexpr \
  -Xcompiler -fopenmp \
  -Xcompiler -Wno-unused-result -Xcompiler -Wno-unused-variable -Xcompiler -Wno-unused-but-set-variable \
  $NVCC_DIAG_SUPPRESS_FLAGS \
  $MIXED_PRECISION_FLAG \
  $CONCURRENT_STREAMS_FLAG \
  $EXTRA_DEFINES \
  $CUDA_ARCH_FLAGS \
  -I"$ROOT/source" -I"$ROOT/../include" -I"$LEGACY_INC" -I"$EIGEN_INC" \
  "$SRC" "$LEGACY_INC/Faddeeva.cc" -lcublas -lcusolver -o "$OUT"
set +x
echo "[compile] done: $OUT"
if [[ "$ALLOBJ_FLAG" != "1" && "${ALLOBJ_FLAG,,}" != "true" ]]; then
  # This is only a reminder; the package can still compile the default binary.
  echo "[compile] note: for component caches compile with V32ZU_ALLOBJ_QC=1"
fi

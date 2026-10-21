#!/usr/bin/env bash
set -euo pipefail
MODE=${1:-fit}
COARSEN=${2:-10000}
CONFIG=${3:-configs/config_v32w_QC_fitter_norm_refine_v2_fixedcoarse.in}
THREADS=${OMP_NUM_THREADS:-18}
export OMP_NUM_THREADS="$THREADS"
export EIGEN_INC=${EIGEN_INC:-/usr/include/eigen3}
export MINUIT2_ROOT=${MINUIT2_ROOT:-${CONDA_PREFIX:-}}
CXX=${CXX:-}
if [[ -z "$CXX" ]]; then
  if [[ -n "${CONDA_PREFIX:-}" && -x "$CONDA_PREFIX/bin/x86_64-conda-linux-gnu-c++" ]]; then
    CXX="$CONDA_PREFIX/bin/x86_64-conda-linux-gnu-c++"
  else
    CXX=g++
  fi
fi
if [[ -n "$MINUIT2_ROOT" ]]; then
  export LD_LIBRARY_PATH="$MINUIT2_ROOT/lib:${LD_LIBRARY_PATH:-}"
fi
mkdir -p bin logs cache output_v32w_QC_fitter_norm_refine_v2_fixedcoarse
CFG_RUN="configs/config_v32w_QC_fitter_norm_refine_v2_fixedcoarse_run.in"
cp "$CONFIG" "$CFG_RUN"
python3 - <<PY
from pathlib import Path
import re
p=Path('$CFG_RUN')
s=p.read_text()
s=re.sub(r'^coarseN\s*=.*$', 'coarseN = $COARSEN', s, flags=re.M)
p.write_text(s)
PY

echo "[compile] v32w QC fitter normalized determinant + verbose cache/refinement progress"
echo "[compile] mode=$MODE coarseN=$COARSEN threads=$THREADS"
echo "[compile] CXX=$CXX"
echo "[compile] EIGEN_INC=$EIGEN_INC"
echo "[compile] MINUIT2_ROOT=$MINUIT2_ROOT"
"$CXX" -std=c++17 -O3 -fopenmp \
  -I"$EIGEN_INC" -I. -Isource \
  ${MINUIT2_ROOT:+-I$MINUIT2_ROOT/include -I$MINUIT2_ROOT/include/Minuit2} \
  source/qc_fitter_norm_refine_v2.cpp \
  source/Faddeeva.cc \
  ${MINUIT2_ROOT:+-L$MINUIT2_ROOT/lib} -lMinuit2 -lMinuit2Math \
  -o bin/qc_fitter_norm_refine_v2_fixedcoarse

echo "[run] starting fit with config $CFG_RUN"
echo "[run] log will be written to logs/v32w_QC_fitter_norm_refine_v2_fixedcoarse.log"
./bin/qc_fitter_norm_refine_v2_fixedcoarse "$CFG_RUN" 2>&1 | tee logs/v32w_QC_fitter_norm_refine_v2_fixedcoarse.log

echo "[done] outputs in output_v32w_QC_fitter_norm_refine_v2_fixedcoarse"

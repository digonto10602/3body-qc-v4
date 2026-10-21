# v33d CPU/GPU Component Cache Generator Package

This package is the self-contained cache-generator bundle for the v33d KKpi
component caches. It keeps the GPU path and the CPU path side by side so you can
build, validate, and run both from one directory tree.

## What is inside

- `gpu_component_cachegen/`
- `cpu_component_cachegen/`
- `legacy_source/`
- `k3iso01_swave_cpu_cache_qc_omp_package/source/fv_projector_cartesian_l1_v30m.hpp`
- `codex_report.md`

The package intentionally omits run logs, previous report dumps, and unrelated
workspace artifacts.

## How the cache generator works

The generator produces one binary component cache per ingredient and irrep.
The saved objects are:

- `plm_config`
- `klm_config`
- `K2inv`
- `F2`
- `G`
- `F3`
- `F3inv`
- `Vsel`
- `projected_F3inv`

The production flow is:

```text
config -> basis construction -> component matrix build -> solve -> projector
-> cache write -> sidecar metadata
```

## Execution tree

```text
1. Choose backend
   - CPU: `cpu_component_cachegen/scripts/run_cpu_component_cachegen.sh`
   - GPU: `gpu_component_cachegen/scripts/run_gpu_component_cachegen.sh`

2. Load config
   - GPU: shell config is sourced
   - CPU: INI-style config is parsed

3. Build energy grid
   - Lbyas sweep
   - `coarseN` energy rows
   - irrep list

4. Build basis/configuration objects
   - `plm_config`
   - `klm_config`

5. Build physics ingredients
   - `K2inv`
   - `F2`
   - `G`
   - `F3`
   - `F3inv`

6. Build the projector
   - `Vsel`
   - `projected_F3inv`

7. Write cache files
   - one binary per ingredient
   - one sidecar metadata file per cache stream

8. Emit the run log
   - per-L summaries
   - per-irrep summaries
   - total runtime summary
```

## Key functions and where they live

### GPU path

- `v32zt_write_matrix_record`
  - `gpu_component_cachegen/source/F3_gpu_cachegen_v32zu.cu:849`
  - writes one matrix record into the binary cache stream
- `v32zt_process_swap_bytes`
  - `gpu_component_cachegen/source/F3_gpu_cachegen_v32zu.cu:1192`
  - reads host swap usage from `/proc/self/status`
- `v32zt_should_use_gpu_lane`
  - `gpu_component_cachegen/source/F3_gpu_cachegen_v32zu.cu:1267`
  - decides whether a task should remain on the GPU lane
- ordered flush loop
  - `gpu_component_cachegen/source/F3_gpu_cachegen_v32zu.cu:2062`
  - flushes completed chunks in order
- single-task wave fallback
  - `gpu_component_cachegen/source/F3_gpu_cachegen_v32zu.cu:2131`
  - forces a single-task wave when the writer or host memory says to back off
- `main`
  - `gpu_component_cachegen/source/F3_gpu_cachegen_v32zu.cu:2527`
  - parses CLI arguments and starts the generator

### CPU path

- `config_from_kv`
  - `cpu_component_cachegen/source/generate_cpu_component_cache.cpp:136`
  - parses the config file into runtime settings
- `write_exact`
  - `cpu_component_cachegen/source/generate_cpu_component_cache.cpp:185`
  - low-level binary writer helper
- `write_matrix_record`
  - `cpu_component_cachegen/source/generate_cpu_component_cache.cpp:195`
  - writes one matrix record into the output cache
- `build_component_entry`
  - `cpu_component_cachegen/source/generate_cpu_component_cache.cpp:223`
  - builds the physics objects for one `(Lbyas, irrep, Ecm)` row
- `main`
  - `cpu_component_cachegen/source/generate_cpu_component_cache.cpp:347`
  - drives the full CPU cache-generation loop

### Scripts

- `gpu_component_cachegen/scripts/compile_gpu_component_cachegen.sh`
  - compiles the all-object GPU binary with `nvcc`
- `gpu_component_cachegen/scripts/run_gpu_component_cachegen.sh`
  - loads the GPU config, builds the L grid, runs the binary, and tees logs
- `cpu_component_cachegen/scripts/compile_cpu_component_cachegen.sh`
  - compiles the CPU generator with `g++` and OpenMP
- `cpu_component_cachegen/scripts/run_cpu_component_cachegen.sh`
  - builds on demand and runs the CPU generator

## Decision tree used to get to the current generator

```text
Start from crash report
  -> confirm RAM + swap exhaustion
  -> confirm duplicate generator process symptom
  -> inspect scheduler and writer path
  -> add host-pressure aware fallback
     -> writer-bound single-task waves
     -> prefer the existing cache writer path, not a new abstraction
  -> verify on small smoke runs
  -> benchmark fused vs separate F2/K2
     -> choose the faster separate-kernel path for coarseN=1000
  -> run production coarseN=20000
     -> confirm one process only
     -> confirm zero swap
     -> confirm all requested irreps completed
```

## Current status at packaging time

- `Lbyas=24`, `coarseN=20000`
- `100_A2`, `110_A2`, `111_A2`, `200_A2` completed
- production cache files were written successfully
- the host-pressure guard and writer-bound single-task behavior are part of the current GPU path

## Build and run

### GPU

```bash
cd gpu_component_cachegen
V32ZU_ALLOBJ_QC=1 bash scripts/compile_gpu_component_cachegen.sh
bash scripts/run_gpu_component_cachegen.sh 24
```

### CPU

```bash
cd cpu_component_cachegen
bash scripts/compile_cpu_component_cachegen.sh
bash scripts/run_cpu_component_cachegen.sh
```

## Prerequisites

- CUDA toolkit for the GPU binary
- `g++` with OpenMP support for the CPU binary
- Eigen headers available via `EIGEN_INC` or a system install
- `Faddeeva.cc/.hh/.h` and the legacy physics headers in `legacy_source/`


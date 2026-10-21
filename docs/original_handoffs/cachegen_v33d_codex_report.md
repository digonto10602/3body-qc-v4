# Codex Report: Full Cache Generator Package

## Goal

Build a self-contained package for the v33d CPU/GPU cache generator, keep only
the files needed to run the generator successfully, add human-readable
documentation, and preserve the current validated state.

## Starting point

The work started from a crash investigation on the GPU cache generator. The
failure mode was:

- host RAM filled up
- swap filled up
- two generator instances appeared during the last failing run

The root cause was in the GPU scheduler / writer interaction in
`gpu_component_cachegen/source/F3_gpu_cachegen_v32zu.cu`.

## Changes made earlier in this session

Single source patch:

- `gpu_component_cachegen/source/F3_gpu_cachegen_v32zu.cu`

What changed:

- added `v32zt_process_swap_bytes()` at line 1192
- extended the lane-selection logic in `v32zt_should_use_gpu_lane()` at line 1267
- made the scheduler treat `rss + swap` as host pressure
- forced `single_task` waves when the writer is the bottleneck
- logged `swap_GB` and `host_pressure_GB`
- preserved the existing cache write path and cache format

The important behavior shift was simple:

```text
if host pressure is high or the writer is bound, stop trying to pack more work
and drain a single task cleanly.
```

## Validation performed before packaging

### Build validation

- rebuilt the GPU all-object binary successfully
- compiled and ran the GPU generator through the package scripts

### Smoke validation

- `coarseN=2` on `000_A1m`
- `coarseN=20` on `000_A1m`

Purpose:

- confirm the binary still runs
- confirm the new scheduler branch is reachable
- confirm cache writing still completes

### Benchmark validation

- `Lbyas=24`, `irrep=100_A2`, `coarseN=1000`
- compared separate-kernel and fused F2/K2 paths

Result:

- separate-kernel path was faster
- chosen path for the production coarse20000 run was the separate-kernel config

### Production validation

Completed production run:

- `Lbyas=24`
- `coarseN=20000`
- irreps: `100_A2`, `110_A2`, `111_A2`, `200_A2`

Observed:

- one generator process only
- swap remained `0`
- all requested irreps completed
- cache files were written successfully

## Current status

The current state at packaging time is a finished, validated GPU and CPU
package plus the production cache results already generated in the working
tree.

## Decision tree used

```text
1. Do not rewrite physics
   -> keep determinant / projector definitions unchanged

2. Keep CPU as correctness reference
   -> use the CPU path as the trusted baseline

3. Fix the actual crash root cause
   -> writer-bound single-task scheduling
   -> host-pressure detection from RSS + swap

4. Verify on small inputs first
   -> smoke runs

5. Benchmark before choosing the production GPU path
   -> fused vs separate F2/K2

6. Choose the faster verified path
   -> separate-kernel path for the production coarse20000 run

7. Confirm production behavior
   -> no duplicate process
   -> no swap growth
   -> successful cache output
```

## Files included in the package

- `gpu_component_cachegen/`
- `cpu_component_cachegen/`
- `legacy_source/`
- `k3iso01_swave_cpu_cache_qc_omp_package/source/fv_projector_cartesian_l1_v30m.hpp`
- `README.md`

## Files changed for the package

Package-local only:

- `gpu_component_cachegen/scripts/run_gpu_component_cachegen.sh`
- `gpu_component_cachegen/scripts/compile_gpu_component_cachegen.sh`
- `cpu_component_cachegen/scripts/run_cpu_component_cachegen.sh`
- `cpu_component_cachegen/scripts/compile_cpu_component_cachegen.sh`
- `README.md`
- `codex_report.md`

What changed in the scripts:

- added explanatory comments
- kept command behavior unchanged
- preserved the same defaults and fallback compile logic

## Tests and checks run while building the package

- log inspection for the completed `100_A2` / `110_A2` / `111_A2` / `200_A2` production run
- process table checks to confirm a single generator process
- `VmSwap` checks during the crash investigation
- benchmark timing comparison for `coarseN=1000`
- zip integrity check on the minimal package built earlier

## Output artifacts

Earlier minimal package:

- `v33d_gpu_component_cachegen_minimal_20260629.zip`

Current full package:

- `full_cache_generator_package_20260629/`

## Risks / limits

- The package still depends on system CUDA, OpenMP, and Eigen headers.
- The GPU binary remains architecture- and driver-sensitive.
- The package excludes logs and previous run outputs on purpose.


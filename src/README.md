# src

All C++ and CUDA sources are copied **verbatim** from the newest package that had them; see docs/PROVENANCE.tsv for each file's md5 and origin. Includes are flat (`#include "file.h"`) and resolved through `-I` for every `src/*` folder (see build.sh). The include closure of each tool was checked on 2026-10-21 and no files are missing.

| folder | entry points |
|---|---|
| `core/` | header-only physics (+ `Faddeeva.cc`) |
| `fitter/` | `qc_fitter_norm_refine_v2_multiL.cpp` (main). It `#include`s `qc_fitter_norm_refine_v2.cpp`, the `K3df_minuit_fit_v32f` header and the 2203-line `debug_v30r_*.cpp` (which provides helpers such as the irrep labels and the projector cache). Splitting these into proper headers is planned. |
| `cache_io/` | `v33h_patched_gpu_cache_oldscale_det_scan.cpp` (main unless `V33H_NO_MAIN`) |
| `cache_build/` | `v33p_build_projected_basis_hotwindow_cache.cpp`, `v33g_build_runtime_k3basis_cache.cpp` |
| `cachegen_gpu/` | `F3_gpu_cachegen_v32zu.cu` |
| `cachegen_cpu/` | `generate_cpu_component_cache.cpp` |
| `gpu_det/` | `v33p_cpuassembled_gpu_det_scan.cpp` + `v33p_gpu_det_lu.cu` (experimental) |

Version choices where copies differed:
- `K3df_minuit_fit_v32f_*.hpp`, `qc_fitter_norm_refine_v2.cpp`, `v33g_runtime_k3basis_cache.hpp`: the L20 versions, which are supersets of the older ones.
- `v33g_build_runtime_k3basis_cache.cpp`: the v33n copy (natural complex read; the L20 copy swaps once).
- The CPU cache generator originally used the older `K3df_minuit_fit_v32f` header. With the L20 header it also precomputes the projected K3 basis per row (extra work only). This is not yet compiled or tested.

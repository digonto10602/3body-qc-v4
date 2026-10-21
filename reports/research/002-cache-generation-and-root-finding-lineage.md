# Report B — F3inv/Vsel cache generation (CPU + GPU) and root-finding lineage for the projected QC determinant

> **Research report — consolidation investigation.** Date: 2026-10-21. Git commit: none (written before the repo existed; first commit of `3body-qc-v4`). Prompt in force: `prompts/000-consolidation-and-roadmap.md`. Related: `reports/research/001`–`006`.
> Absolute paths below refer to the unpacked original zips (`/home/digimonk/Projects/three-body-qc-codes/extracted/...`) and to a temporary scratchpad. For each file's new location in this repo, see `docs/PROVENANCE.tsv`.
> The system is two kaons + one pion treated as a 2+1 system (identical kaons). Where this report says "K K-bar pi", read "KKπ"; that wording came from the investigation brief and is not established. See `docs/OPEN_ISSUES.md` #7.


- Date: 2026-10-21
- Git commit: n/a. The extracted folders are not git repositories; every citation below is `path:line` relative to `/home/digimonk/Projects/three-body-qc-codes/extracted/`.
- Prompt in force: coordinator brief "Report B" (cache generation + root-finding evolution). Read-only; nothing was built or run except small read-only Python one-liners for byte counts and size estimates.
- Related reports: the other consolidation reports (A, C, …) written to the same `scratchpad/reports/` folder. This report does not depend on them.
- Short names used below:
  - `FCG` = `full_cache_generator_package_20260629-main/full_cache_generator_package_20260629-main`
  - `L20P` = `KKpi-L20-user-validated-hotwindow-fitter-main/KKpi-L20-user-validated-hotwindow-fitter-main`
  - `V33P` = `KKpi-v33p-cache-anchored-classifier-fitter-main/KKpi-v33p-cache-anchored-classifier-fitter-main`
  - `V33N` = `kkpi-v33n-corrected-gpu-cache-reader-main/kkpi-v33n-corrected-gpu-cache-reader-main`
  - `Q3V3` = `3body_quantization_v3-master/3body_quantization_v3-master`
  - `OLD` = `_nested_old_versions`

---

## 1. What was asked

1. The exact binary formats: the combined raw `*_F3inv_Vsel_gpu.bin` file written by v32zu, its `.meta.json` sidecar, the per-object "component" caches written by the CPU and GPU generators, and where the `variant_04_real_imag_swapped` complex-number convention comes from. Also which objects are stored for each energy, and how the energy grid is built.
2. How the GPU generator builds the matrices: kernels, batched LU solves, CUDA streams, the host-memory guard, and fused versus separate F2/K2 kernels. Also any performance numbers that were recorded.
3. The history of root-finding methods for the QC determinant. This includes a full search of all packages (including `_nested_old_versions`) for gradient-, slope- or derivative-based root detection.
4. Why the roots are hard to find, as a list of known pitfalls.
5. Which generator is the current one, and the smallest set of files needed to regenerate caches.

Terms used throughout:

- **QC**: the three-body quantization condition. In this code it is $\det\!\big[V^\dagger (F_3^{-1}+\mathcal K_{\rm df,3})V\big]=0$.
- **Projector / `Vsel`**: an $N\times N_{\rm proj}$ matrix whose columns span one lattice irrep (one symmetry sector).
- **Cache**: binary files of matrices precomputed on an energy grid, so the fitter does not rebuild $F_3^{-1}$ every time.
- **Batched LU**: an LU factorisation (Gaussian elimination with row pivoting) of many same-size matrices in one GPU call.
- **CUDA stream**: an independent GPU work queue. **cuBLAS/cuSOLVER**: NVIDIA's dense linear-algebra libraries.
- **Sign flip / bracket**: two neighbouring grid energies where the real part of the determinant changes sign.

---

## 2. Actions taken

- Read all of `FCG/README.md`, `FCG/codex_report.md`, `FCG/legacy_source/F3_GPU_pipeline_version_history.md`, the GPU generator `FCG/gpu_component_cachegen/source/F3_gpu_cachegen_v32zu.cu` (2699 lines, the relevant sections), and its included v16 file, Vsel projector header, configs and scripts.
- Read all of the CPU generator `FCG/cpu_component_cachegen/source/generate_cpu_component_cache.cpp` (428 lines) plus its config and scripts.
- Traced the `variant_04` convention through every reader that names it: `L20P/source/gpu_cache_reader_convention.hpp`, `v33h_patched_gpu_cache_oldscale_det_scan.cpp`, `qc_fitter_norm_refine_v2_multiL.cpp`, `v33i_single_row_gpu_reader_path_identity_test.cpp`, `debug_matrix_read_conventions_L20_100_A2.cpp` and `compare_gpu_cache_vs_cpu_openmp_det_oldscale.cpp`. I also read `V33N/GPU_CACHE_READER_FIX.md` and the v33n reader.
- Read the root-finding code in every legacy version from v7 to v21 (GPU) and v21 to v24 (CPU), the v29 to v31 debug zero scanners, `qc_fitter_norm_refine_v2.cpp` (v2/v3/v4 classifiers), the v33i/j/l classifier sweep scripts, `pole_searching.h`, and the READMEs in `OLD/KKpi_repro_context_transfer_v31u_package/.../docs/legacy/`, `OLD/v31y_*`, `OLD/v31zo`, `OLD/v32r`.
- Ran a case-insensitive grep over all text files in every package for `gradient|slope|derivative|dlogdet|d/dE|curvature|\bgrad\b|grad_`, then widened it to `ddet|dlog|finite diff|inflection|local min|kink|shoulder|crossing direction|monoton|np.diff|Hellmann`. I removed content duplicates with `md5sum`: 1043 matching paths collapsed to 29 unique files for the first pattern. Results are in §5.3.
- Resolved the full `#include` closure of both generators with a small read-only Python script (`scratchpad/incl.py`) to get the minimal file set (§7).
- Computed byte counts, grid spacing and cache-size estimates with `python3 -I -c` one-liners. The formulas are shown wherever a number is quoted.

---

## 3. Cache formats and stored objects (Question 1)

### 3.1 Physics objects stored per energy row

The algebra is the same in the CPU and GPU generators: `FCG/cpu_component_cachegen/source/generate_cpu_component_cache.cpp:251-273` and `F3_gpu_cachegen_v32zu.cu:1777-1843`.

$$
H = K_2^{-1} + F_2 + G,\qquad H X_1 = F_2,\qquad
F_3 = \tfrac13 F_2 - F_2 X_1 \;\big(= \tfrac13F_2 - F_2 H^{-1}F_2\big),
$$
$$
F_3^{-1}\ \text{from}\ F_3 X = \mathbb 1,\qquad
\texttt{projected\_F3inv} = V_{\rm sel}^\dagger F_3^{-1} V_{\rm sel}.
$$

| object | shape | meaning |
|---|---|---|
| `plm_config` | GPU: $N_1\times5$; CPU: $5\times N_1$ (see pitfall P-F4) | spectator basis for flavour-1 (KK̄π, spectator K): 5 rows = $p_x,p_y,p_z,\ell,m$ |
| `klm_config` | GPU $N_2\times5$ / CPU $5\times N_2$ | spectator basis for flavour 2 (spectator π) |
| `K2inv` | $N\times N$, $N=N_1+N_2$ | block-diagonal $\mathrm{diag}(K_{2,1}^{-1},\,2K_{2,2}^{-1})$ (factor 2 at `F3_gpu_cachegen_v32zu.cu:1122`) |
| `F2`, `G` | $N\times N$ | finite-volume F and G matrices |
| `F3`, `F3inv` | $N\times N$ | as above |
| `Vsel` | $N\times N_{\rm proj}$ | eigenvectors of the irrep projector $P_I$ with eigenvalue in $[1-0.05,\,1+0.05]$ (`projections_gpu_v32y_exact_device.cuh:264-277`) |
| `projected_F3inv` | $N_{\rm proj}\times N_{\rm proj}$ | $V^\dagger F_3^{-1}V$, computed on the host (`F3_gpu_cachegen_v32zu.cu:1906`) |

Physics constants are hard-coded in the GPU binary and copied into the sidecar: `F3_device_resident_true_batch_v16_custom_batched_getrs.cu` (gpu copy) `:558-592` and `F3_gpu_cachegen_v32zu.cu:2531-2580`.

- Masses: $a_t m_\pi=0.06906$, $a_t m_K=0.09698$.
- Anisotropy: $\xi=3.444$.
- Other parameters: $\eta_1=1$, $\eta_2=0.5$, $\alpha=0.5$, $\epsilon_h=0$, `max_shell_num`=20, `Q0norm`=true.
- Partial waves: `waves_vec_1`=$\{0,1\}$, `waves_vec_2`=$\{0\}$.
- Scattering parameters: `scatter1_00`=4.04, `scatter1_10`=−43.2, `scatter2_00`=4.12.
- Parity: forced to −1 at `:2535`.
- Projector tolerance: `eig_tol`=0.05, hard-coded at `:1643`.

The GPU binary cannot change the scattering parameters at run time. The config only documents them (`configs/config_gpu_component_cachegen.sh:11-16`).

### 3.2 Energy grid

The grid has equal spacing in the centre-of-mass energy $E_{\rm cm}$: `F3_gpu_cachegen_v32zu.cu:2352-2357` (GPU) and `generate_cpu_component_cache.cpp:367,389` (CPU).

$$
E_i = E_{\min} + i\,\frac{E_{\max}-E_{\min}}{N_c-1},\quad i=0,\dots,N_c-1,
\qquad E_{\rm lab}=\sqrt{E_i^2+|\vec P|^2},\ \ \vec P=\tfrac{2\pi}{L}\vec n_P,\ L=\xi\,L/a_s .
$$

- Production window: $E_{\min}=0.26310$, $E_{\max}=0.36$, in lattice units $a_t^{-1}$, with $N_c=20000$ (`configs/config_gpu_component_cachegen_sep_f2k2_coarse20000.sh:8`). The three-particle threshold is $2a_tm_K+a_tm_\pi=0.26302$, so the window starts $8\times10^{-5}$ above it. `V33P/README.md:52-55` warns that the nominal value 0.26301 must not be substituted for the sidecar value.
- Grid spacing: $\Delta E=0.0969/19999=4.84524\times10^{-6}\,a_t^{-1}$ (exact arithmetic). This matches the bracket widths seen in user-label CSVs, e.g. $0.340575423771-0.340570578529=4.84524\times10^{-6}$ (`V33N/CLASSIFIER_RESTART_NOTES.md:15`).
- Default $N_c$ is 10000 in the base GPU config and 1000 in the CPU config. Comments in both say "50000 for production", but 20000 is what the v33e caches actually use.
- Lbyas sweep: `scripts/run_gpu_component_cachegen.sh` takes one $L/a_s$, or `Lmin Lmax [step]` (default step 0.2), and builds the list with Python `Decimal` to avoid floating-point drift (lines 61-88). Each value of $L$ gets its own `.../v32zu_Lbyas<tag>_gpu_cache/` directory. Production caches named in packages: $L/a_s=20$ (v33e, one directory per irrep) and $L/a_s=24$ (100/110/111/200_A2 per `FCG/README.md:166-171`). Spatial $L=\xi\cdot L/a_s=68.88\,a_t$ for $L/a_s=20$.
- Default irreps: `000_A1m`, `100_A2`, `110_A2`, `111_A2`, `200_A2` (`F3_gpu_cachegen_v32zu.cu:662-671`).

### 3.3 Combined raw cache `v32zu_Lbyas<L>_xi<xi>_irrep<label>_coarse<N>_F3inv_Vsel_gpu.bin`

The name is built at `F3_gpu_cachegen_v32zu.cu:2329-2333`, and the tag `p` replaces `.` (`:687-696`). The writer is `v32zt_write_cache_record` (`:718-754`), called from the ordered flush `v32zu_write_staged_rows` (`:1953-1976`). All values are little-endian.

| offset | type | content |
|---|---|---|
| file start | 39 bytes ASCII | `"KKPI_F3INV_VSEL_GPU_V32ZT_EXACT_V32Y 1\n"` (`:2339`) |
| per record +0 | `uint64` | magic `0x5653325a4f524543` (`:58`). On disk this reads as ASCII `CEROZ2SV`. The source comment "VS2ZPREC" is a misspelling of the big-endian bytes `VS2ZOREC`. |
| +8 | `int32` ×5 | `grid_i, dim1(=N_1), dim2(=N_2), n(=N), vdim(=N_proj)` |
| +28 | `double` ×2 | `Ecm`, `En` (lab energy) |
| +44 | `16·n²` bytes | `F3inv`, **column-major**, each element a raw `cuDoubleComplex` `{x=Re, y=Im}` (memcpy of the GPU buffer, `:752`) |
| … | `16·n·vdim` bytes | `Vsel`, column-major, same element layout |

The record size is $44+16N^2+16N N_{\rm proj}$ bytes. There is no record count and no index, so readers scan sequentially and check the magic.

Record order is **not** energy order. Rows are grouped by the shape key $(N_1,N_2,N)$ in `std::map` order (key comparison at `v16:867-872`), and are in increasing `grid_i` only within a group (`:2391-2400`, ordered flush `:2047-2068`). Readers must sort by `Ecm` or `grid_i`. The trusted v33h reader does this at `L20P/source/v33h_patched_gpu_cache_oldscale_det_scan.cpp:331-347`.

Rows with an empty basis ($N=0$) are skipped with a log line and never written (`:2384-2387`). The sidecar then reports `rows_written < coarseN`, and the v33h validator rejects that (`v33h…:311`).

Size example for L20/100_A2: $N/N_{\rm proj}$ runs from 34/11 to 108/30 over 12 dimension segments (`V33P/output/v33p_classifier_restart_L20_100_A2_E026310_0360/cache_window_audit_E026310_0360.md:20-33`). Summing $44+16N^2+16NN_{\rm proj}$ over the 20000 rows gives about **2.60 GB** for the combined file, and about 10.8 GB for the F2/G/K2inv/F3/F3inv/Vsel/projected component files of the same irrep. This is my estimate from the formula, not a measured file size.

### 3.4 Combined-cache sidecar `<bin>.meta.json`

Written by `v32zt_write_sidecar_json` (`F3_gpu_cachegen_v32zu.cu:1301-1371`) as hand-assembled JSON with 17 significant digits. Fields:

- **Identity:** `version` (either `"v33d_gpu_component_cachegen_all_ingredients"` when compiled with `V32ZU_ALLOBJ_QC_WRITER`, or `"v32zu_gpu_cachegen_F3inv_Vsel_exact_v32y_opt9_concurrent_streams"`, `:59-63`), `cache_kind:"F3inv_plus_Vsel_gpu_full"`, `cache_file`, `projection_source` (names `generate_swave_cpu_cache_qc.cpp`, which is **not** in any package I read), `basis_ordering:"config_maker_4_momentum_first"`.
- **Sector:** `irrep_label`, `irrep`, `irrep_tag`, `nnP`, `Lbyas`, `xi`, `spatial_L`, `Ecm_min`, `Ecm_max`, `coarseN`.
- **Physics:** `waves_vec_1/2`, `scatter1_00/10`, `scatter2_00`, `masses{atmK,atmpi}`, `alpha`, `epsilon_h`, `max_shell_num`, `Q0norm`, `device_id`. Parity and `eig_tol` are **not** recorded.
- **Run flags:** `optimization_flags{use_fused_f2k2, validate_fused_f2k2, use_kinematic_precompute, use_stream_pipeline_requested/enabled, use_concurrent_streams, concurrent_streams_safety_fraction, use_multi_gpu, use_mixed_precision_f2k2}`, `writer_thread`.
- **Counts:** `rows_written`, `chunks`, `max_total_dim`, `max_vdim`, `bytes_written_binary_payload`.
- **Timing:** `timing_seconds{basis_construction_time, F2_K2_time, G_F3_matrix_build_time, matrix_inverse_solve_time, Vsel_construction_time, binary_write_time, total_irrep_time}`.

Readers treat the sidecar as authoritative for the window, row count and provenance (`v33h…:287-329`; `V33P/README.md:52`).

### 3.5 Component caches (format `KKPI_COMPONENT_CACHE_V001`), CPU and GPU

There is one file per object: `v32zu_Lbyas…_<object>_gpu.bin` for GPU (`F3_gpu_cachegen_v32zu.cu:804-819`) and `cpu_Lbyas…_<object>.bin` for CPU (`generate_cpu_component_cache.cpp:285-289`).

- **File header** (GPU `:912-942`, CPU `:299-317`): the 25-byte ASCII `KKPI_COMPONENT_CACHE_V001` (no newline), then `string object`, `string created_by` (`gpu_v32zu_component_cachegen` or `cpu_v33d_component_cachegen`), `double Lbyas`, `double xi`, `string irrep`, `int32 coarseN`, `int32 row_count`. Here `string` means a `uint32` length followed by the bytes.
- **Record** (GPU `:849-867`, CPU `:195-203`): `int32 row`, `double Ecm`, `int32 rows`, `int32 cols`, then **row-major** pairs `(double Re, double Im)`. The combined raw file is column-major, so the two formats differ in element order.
- **Sidecars:**
  - GPU, `v32zt_finalize_allobj_sidecars` (`:951-1027`): `component_cache_file_format`, `cache_kind`=object, physics block, `layout:"row-major records: energy_index,Ecm,rows,cols,re,im"`, `payload_bytes`, `payload_sha256` (from shelling out to `sha256sum`), host, CUDA and GPU info, timings. It also writes constant placeholders: `cuda_streams:1`, `cuda_batch_E:0`, `peak_gpu_bytes:0`, `validation:{"validated_against":"smoke_only","max_abs_diff":null}`. These are not measured values.
  - CPU (`:319-345`): version `v33d_cpu_component_cachegen_all_ingredients`, window, waves, scatter params, `attempted_rows`, `success_rows`, elapsed, timestamps.
- **Row-count inconsistency:** the GPU `row_count` is the number of non-empty rows. The CPU writes `row_count = coarseN` in the header (`:385`) but skips failed rows (`:393-406`), so its header can over-count.

### 3.6 The `variant_04_real_imag_swapped` convention: where it comes from

**Finding: the writer does not swap real and imaginary parts.** `v32zt_write_cache_record` memcpy's `cuDoubleComplex` buffers (`F3_gpu_cachegen_v32zu.cu:752-753`). That type is `double2{x,y}` with `x` = real part, and every kernel builds it as `make_cuDoubleComplex(re, im)` (e.g. `:1131`). The component writer writes `(Re, Im)` explicitly (`:861-864`). So the bytes on disk are in the natural order $(\mathrm{Re},\mathrm{Im})$.

The "swap" comes from the readers:

1. The diagnostic tools that found variant 04 read each complex entry as `comp(read_f64(), read_f64())`: `L20P/source/debug_matrix_read_conventions_L20_100_A2.cpp:358,363` and `L20P/source/v33i_single_row_gpu_reader_path_identity_test.cpp:248,253`. C++ does not specify the order in which function arguments are evaluated. With GCC on x86-64 they are typically evaluated right to left, so the second double read becomes the real part. *This compiler behaviour is my inference; it is not documented in the packages.* In that case the "current reader" baseline (`variant_00`) was itself swapped, and swapping again (`variant_04`, `:375`) restored the true matrix. The symptoms recorded before the fix (Hermiticity error ≈ 2, nonzero imaginary determinant; `V33N/GPU_CACHE_READER_FIX.md:7-11`) fit this explanation.
2. The shared "corrected" reader then hard-coded the swap with sequential reads, `return comp(b, a)` (`L20P/source/gpu_cache_reader_convention.hpp:179-186`). By itself that **is** a swap of the stored data.
3. The trusted v33h scanner applies `swap_real_imag` a **second** time after that reader (`v33h_patched_gpu_cache_oldscale_det_scan.cpp:240-242, 383-393`). The production fitter copies this double swap and explains why: "the L20 audit showed the one-swap path loses the accepted sign-change bracket" (`L20P/source/qc_fitter_norm_refine_v2_multiL.cpp:664-705`). Two swaps cancel, so the trusted path reads the natural $(\mathrm{Re},\mathrm{Im})$ order, which is exactly what the writer stored. The older v33n fitter reader also reads `comp(re,im)` directly with no swap (`V33N/source/qc_fitter_norm_refine_v2_multiL.cpp:627-646`).
4. At least one consumer still applies a **single** swap: `L20P/source/compare_gpu_cache_vs_cpu_openmp_det_oldscale.cpp:300-301` (and anything else that calls `read_gpu_cache_matrix_col_major` without the second swap). Those paths read $\tilde z = \mathrm{Im}\,z + i\,\mathrm{Re}\,z = i\,\overline{z}$. For a real `Vsel` (it comes from a real symmetric eigensolver, `cusolverDnDsyevd`, `projections_gpu_v32y_exact_device.cuh:264`), $(iV)^\dagger M (iV)=V^\dagger M V$, so the `Vsel` swap does no harm. The `F3inv` swap is harmful: it becomes $i\,\overline{F_3^{-1}}$.

Recommendation for the clean repo: read $(\mathrm{Re},\mathrm{Im})$ in sequence, drop the `variant_04` name, and add a unit test that round-trips a known non-Hermitian complex matrix through writer and reader.

---

## 4. GPU generator internals and performance (Question 2)

### 4.1 Pipeline per irrep (`v32zt_generate_one_irrep`, `F3_gpu_cachegen_v32zu.cu:2308-2519`)

1. **Prepass on CPU with OpenMP** (`:2351-2367`, `v32zt_energy_shape_prepass_momentum_first` `:1374-1426`). For every $E_i$ it builds the momentum-first spectator basis (`config_maker_4_momentum_first`) and records $(N_1,N_2,N)$.
2. **Group** rows by $(N_1,N_2,N)$ (`:2389-2398`). Each group is cut into chunks whose size is set by free VRAM: per-row estimate $18N^2\cdot16\,\mathrm{B}+2\,\mathrm{MB}$ (`:1153-1162`), cap $=\lfloor 0.85\cdot\text{free VRAM}/\text{per-row}\rfloor$ (`:1284-1299`). A finer estimate of the whole wave is at `:1219-1251`.
3. **Per chunk** (`v32zt_process_same_shape_chunk_cache`, `:1506-1951`):
   - CPU/OpenMP rebuild of the configs for contiguous upload.
   - **Vsel per energy on the GPU** (`:1634-1665`, `build_projection_Vsel_device_exact_v32y`). The projector $P_I$ is built element-wise with real Wigner-D matrices for $\ell\le1$ (Cartesian $\ell=1$ action, `projections_gpu_v32y_exact_device.cuh:1-60`). It is diagonalised with `cusolverDnDsyevd`, and eigenvectors with $|\lambda-1|\le0.05$ are kept. All rows in a chunk must have the same $N_{\rm proj}$, otherwise the chunk throws (`:1645-1648`).
   - **F2 and K2inv:** either the fused kernel `v32zu_fused_F2_K2inv_batch_same_shape_kernel` (one CUDA block per energy row, 256 threads covering the $N^2$ elements, `:1042-1133`, launch `:1688`), or the separate v14 kernels `v14_build_F2_batch_same_shape_kernel` and `v14_build_K2_batch_same_shape_kernel` (flat one-thread-per-element launch, `:1698-1715`). On the GPU, the imaginary error function `erfi` is summed as a series (`erfi_max_terms=100`, tol $10^{-15}$, `:1695`). The CPU uses `Faddeeva::erfi`.
   - **G** (`v14_build_G_batch_same_shape_kernel`, `:1777`) and $H=F_2+G+K_2^{-1}$ (`v14_combine_H_batch_kernel`, `:1781`, `v16:377-401`).
   - **Solve** $HX_1=F_2$ by batched LU (`:1803`). Then `cublasZgemmBatched` computes $F_2X_1$ (`:1819`), and `v14_build_F3_batch_kernel` forms $F_3=\tfrac13F_2-F_2X_1$ (`:1824`, `v16:403-415`).
   - **Solve** $F_3X=\mathbb 1$ by batched LU, giving the full $F_3^{-1}$ (`:1842`).
   - Copy device to host, then stage or write records (`:1855-1940`).
4. **Batched LU** (`v32zu_cublas_lu_solve_inplace_batched_same_n_stream`, `:1428-1504`):
   - `cublasZgetrfBatched` factorises.
   - A **custom triangular-solve kernel** `v16_getrs_from_getrf_batched_kernel` (`v16:990-1060`) replaces `cublasZgetrsBatched`, which crashed in v15 (`F3_GPU_pipeline_version_history.md`, v15/v16 sections). The kernel uses one CUDA block per (matrix, right-hand-side column) with a single thread doing pivoting, forward and back substitution serially.
   - Batch size 1 falls back to cuSOLVER `getrf/getrs` behind a mutex (`:1443-1451`).
   - Singular LU (`info>0`) is **printed and then ignored** (`:1473-1491`). No finiteness check runs before writing (grep for `isfinite`/`allFinite` in the file finds only unrelated lines). By contrast, the CPU generator drops non-finite rows (`generate_cpu_component_cache.cpp:275-277`).

### 4.2 Streams and the host-pressure guard

- `run_groups_concurrent` (`:1992-2306`) builds one task per chunk and orders tasks by estimated VRAM, largest first (`:2071-2073`). In each "wave" it launches several chunks, each on its own `cudaStreamNonBlocking` stream with its own cuBLAS handle (`:2185-2187`). Results are written in the original task order through `completed_by_order`/`flush_ready` (`:2047-2068`), so the output order does not depend on scheduling.
- Lane choice (`v32zt_should_use_gpu_lane`, `:1267-1282`) uses three tests:
  - **gpu_fit**: wave estimate ≤ free VRAM × safety (0.85).
  - **host_fit**: process RSS + swap < 0.8 × physical RAM. RSS and swap are read from `/proc/self/status` `VmRSS`/`VmSwap` (`:1173-1209`).
  - **worthwhile**: ≥2 rows and predicted time ≥0.25 s. Predicted time $t=\frac{\text{rows}}{3.525}\big(\frac{N}{189}\big)^2$ s (`:1253-1265`).
- Guard (`:2128-2150`): `force_single_task_wave = allobj_writer != nullptr || host_pressure ≥ 0.8·RAM`. **In the all-object production build this is always true**, so production runs one chunk per wave and does not overlap streams (the "writer_bound" mode). `V32ZU_ALLOBJ_QC_WRITER` also forces `use_concurrent_streams=true` (`:2555-2557`), but this only routes work through the wave scheduler.
- Why the guard exists (`FCG/codex_report.md:11-36`): a production run filled RAM and swap, and two generator processes appeared. The fix treats RSS+swap as host pressure and drains one task at a time when the writer is the bottleneck. Validation logged one process and zero swap.
- OOM fallback in the sequential path: on a `cudaMalloc`/out-of-memory exception the chunk is halved and retried (`:2431-2461`).

### 4.3 Fused vs separate F2/K2

`use_fused_f2k2` defaults to true in `main` (`:2546`) and in the base config. The production config `config_gpu_component_cachegen_sep_f2k2.sh:9` sets it to **false**. `FCG/codex_report.md:56-66` records the benchmark at $L/a_s=24$, `100_A2`, $N_c=1000$: "separate-kernel path was faster", with **no timing number given**. An optional validation mode compares fused and separate outputs (`:1720-1775`, stats printed `:2654-2660`). No result from it is packaged.

### 4.4 Performance numbers found (wall-clock seconds, GPU and host unspecified unless stated)

| quantity | value | source |
|---|---|---|
| v32zu scheduler baseline | $N_c=1000$ rows in 283.672794974 s, i.e. 3.525 rows/s at $N=189$ | code comment `F3_gpu_cachegen_v32zu.cu:1253-1257`. The irrep/L is not stated, and there is no uncertainty or log. |
| extrapolation (mine) | 20000 rows at 3.525 rows/s ≈ 5.7×10³ s ≈ 1.6 h per irrep at $N\approx189$ | arithmetic only, not measured |
| v16 true-batch small test ($N=28$, batch 5) | build 3.42 s per group (0.685 s per energy); solve 0.0173 s per group (3.46×10⁻³ s per energy) | `FCG/legacy_source/F3_GPU_pipeline_version_history.md`, v16 section (log excerpt) |
| production coarse20000 | $L/a_s=24$; 100/110/111/200_A2 completed; one process; swap 0 | `FCG/codex_report.md:68-84`. **No wall time packaged** because logs were deliberately excluded (`:166-169`). |
| fitter-side root classification (for context) | 1.37–1.43 s per FCN evaluation (10-FCN mean), 32/32 levels found | `V33N/reports/v33m_production_classifier_benchmark_and_fit_full_report.md:31-39` |

Earlier GPU history in brief (`F3_GPU_pipeline_version_history.md`):

| version | change |
|---|---|
| v3–v8 | CPU/OpenMP builds the matrices, GPU does batched LU and GEMM |
| v6 | recursive halving of chunks on OOM (`F3_gpu_omp_cublas_pipeline_v6_oom_safe_chunks.cu:972-1130`) |
| v10 | GPU-safe builders, but with a round trip through the CPU |
| v11 | device-resident data, one energy at a time, cuSOLVER LU |
| v12 | streams per thread |
| v13 | batched-LU prep |
| v14 | batch kernels added but not selected at run time |
| v15 | batch path forced on, which isolated the `getrsBatched` crash |
| v16 | custom GETRS kernel |
| v17 | `Vsel` reuse per shape group, chunk 512, refine 200 |
| v19 | moving-frame parity fix |

In `Q3V3/` the pre-pipeline experiments are:

- `gpu_varsize_batched_inverse.cu` (root, 382 lines) and `source/gpu_varsize_batched_inverse.cu` (1145 lines, adds RAM/VRAM-aware chunk sizing and cached buffers): CPU-built matrices grouped by size, then `cublasZgetrfBatched` + `cublasZgetrsBatched` to invert.
- `source/gpu_solvers_batched_streams*.cpp`: a QR pipeline over streams (cuSOLVER `Zgeqrf`/`Zunmqr` + `cublasZtrsm`).
- `source/test_gpu.cpp`: drivers (`test_F3iso_gpu*`, `test_projections_gpu`).

None of these write the cache format. They are historical only.

---

## 5. Root-finding lineage (Question 3)

Notation: $y(E)$ is the quantity being scanned. Early versions used $\mathrm{Re}\det[V^\dagger F_3^{-1}V]$ normalised by its maximum over the scan. Later versions used $\mathrm{Re}\det\big[(L\xi)^{-6}V^\dagger(F_3^{-1}+\mathcal K_{\rm df,3})V\big]$ (`V33N/docs/DIGONTO_CLASSIFIER_V3_ALGORITHM.md:9-13`). A "bracket" is $[E_i,E_{i+1}]$ with $y_iy_{i+1}<0$. "Zero-like"/"pole-like" are the code's own labels.

### 5.1 Table of methods (oldest → newest)

| # | method (code name) | where | what it does exactly | outcome recorded |
|---|---|---|---|---|
| 0 | **Secant iteration** on $\det F_3$ | `Q3V3/source/pole_searching.h:11-114` | $E_0=E_1-f_1\frac{E_1-E_2}{f_1-f_2}$ with a $10^{-10}i$ offset, stops when $|f_0|\le10^{-15}$. Needs two guesses near the root. | Used for K3=0 "poles" via `generate_pole.cpp`/`pole_maker.py`. No reports. |
| 1 | **Coarse sign flip + spike/slope-ratio classifier** (`classify_sign_flips_zero_vs_pole`) | `FCG/legacy_source/file_processor.hpp:59-330`; same pattern in pipelines v7–v10 and CPU v21–v22, e.g. `F3_gpu_omp_cublas_pipeline_v7_adaptive_zero_search.cu:1604-1700` | Linear zero estimate in each bracket. Window of ±5 points. `spike_ratio`=max\|y\|/min (or median)\|y\|. **`slope_ratio`=max/min (or median) of \|Δy/ΔE\|** (`file_processor.hpp:135-145, 258-264`). "likely_zero" if \|y\| is below `small_y_threshold`; "likely_pole" if spike or slope spike. | This is the only place where a **slope** decides the label. v9 dropped it as threshold-dependent (`F3_GPU_pipeline_version_history.md`, v9: "Removed dependence on small_y_threshold, spike_ratio_threshold, slope_ratio_threshold"). |
| 2 | **Adaptive zero search** (v7) | `…v7_adaptive_zero_search.cu:2288-2465` | Coarse 2000 points, classify (method 1, call-site thresholds `small_y`=1e-9, `spike`=1e1, `slope`=1e1, `:2297-2303`), re-scan each non-pole window [i−5, i+6] with 300 points, re-classify, keep refined likely_zero, otherwise keep the coarse estimate. | Superseded by v9. |
| 3 | **Inward-shape refinement** (v9 `classify_refined_file_by_inward_shape`) | `…v9_inward_shape_refinement.cu:2494-2735`; GPU v16–v17 `classify_inward_shape` `v16:1700-1735` | On the refined window, define $g_k=\|y_k\|+\|y_{N-1-k}\|$ moving inward. If ≥65 % of steps decrease (margin ≥0.20): zero. If ≥65 % increase: pole. Otherwise ambiguous, and the refinement N is doubled (up to 4 rounds, `v16:1758-1785`). | Removed the old thresholds. Problem: mixes the left and right sides. |
| 4 | **Zoom refine, constant N** (v18) | `F3_device_resident_true_batch_v18_zoom_refine_constantN.cu:1810-1940` | Same classifier as 3. When ambiguous, shrink the window to [flip−5, flip+1+5] of the refined grid and rescan with fixed `Nref`=200, up to 12 rounds, stopping if the width shrinks by <0.0001 %. | — |
| 5 | **Parity fix** (v19) | `…v19_parityfix_v9match.cu:1511-1544` | Moving-frame irreps without a g/u suffix get parity −1 (v18 used +1, which made 001/A2 empty). | Root cause found with `test_projection_v9_vs_v18_gpu_vsel_issue.cu`. |
| 6 | **Exact sign-flip refine** (v20) | `…v20_exact_signflip_refine.cu:1871-1922` | Initial and zoom windows are exactly $[E_i,E_{i+1}]$ (no ±5 padding). | — |
| 7 | **Side-|det| trend refine** (GPU v21, CPU v21–v24 "multi-signflip side_abs") | `…v21_side_abs_trend_refine.cu:1770-1890`; `F3_cpu_openmp_v23_…cpp:1-22`; `F3_cpu_openmp_v24_K3_QC_…cpp:880-1014` (max 12 rounds) | Treats the left side ($0\to j$) and right side ($N-1\to j+1$) separately. Zero if \|y\| decreases inward on **both** sides with fraction ≥0.65 and margin ≥0.20. Pole if it increases on both. Several flips inside one refined mesh are classified separately. v24 applies this to $\det V^\dagger(F_3^{-1}+\mathcal K)V$. | Became the basis of the shoulder classifiers. |
| 8 | **Determinant eigenbranch matching** (v29o–v29q) | `debug_v29h_F3inv_zero_compare_inputfile.cpp:1-19`, `debug_v29q…` | Per coarse interval and per eigenbranch, adjacent-interval eigenvector matching. Compares zeros of $1/\langle F_3\rangle_{\rm iso}$, $\det F_3^{-1}$ and $\det(V^\dagger F_3V)^{-1}$. | "fixes missed determinant candidates" (header). |
| 9 | **SVD-validated zero scan** (v30r) | `debug_v30r_svd_validated_cartesian_l1_100_A2_zero_scan.cpp:1-25`; README `OLD/…v31u…/docs/legacy/README_v30r_svd_validated_100_A2.md` | Branch sign flips, then accept only at a local minimum of the **minimum singular value** $\sigma_{\min}$ of projected QC. Caches `Vsel` per basis signature. | README: "For a blind irrep, prefer the SVD files." |
| 10 | **Eigen-min bisection** (v30s) | `debug_v30s_…:1-4` | Sign flips of the eigenvalue closest to zero (same $N_{\rm proj}$ segment), bisection to $10^{-9}$, tracking the branch by eigenvector overlap. | Failure mode named in v30t: small eigenvalue but QC not singular. |
| 11 | **Residual-checked eigen-min** (v30t) | `debug_v30t_…` | 10 plus acceptance only if the eigenvalue residual is small, the overlap is continuous, and $\sigma_{\min}\le$ tol. Rejection codes `REJECT_FINAL_EIGEN_RESIDUAL_TOO_LARGE`, `REJECT_SVD_CROSSCHECK_FAILED`. | — |
| 12 | **Raw bisection + pole reject** (v30y) and **local-min + pole reject** (v30z) | READMEs v30y, v30z; `debug_v30z_…:1-5` | Accept if $\min\|\lambda_{\rm closest}\|<10^{-6}$, not likely_pole, and stable under a finer local grid. Branch-index switches are recorded but not used to reject. | `CONVERSATION_CONTEXT.md:157-159`: "highlighted problems with closest-eigenvalue envelope switching". |
| 13 | **Closest-eigenvalue crossing direction** (v31a/v31b) | `debug_v31b_closest_eig_midpoint_A2_zero_pole_scan.cpp:1-12, 138-167` | Same $N_{\rm proj}$ only. $\lambda_{\rm closest}$: $+\to-$ is **ZERO**, $-\to+$ is **POLE**, midpoint estimate, no refinement. | **Gradient-sign rule** (§5.3). No result files found in any package. |
| 14 | **Hungarian eigenbranch tracking** (v31l) | `README_v31l_QC_tracked_eigenbranches_hungarian.md` | Overlap $O_{ab}=\|v_a^\dagger v_b\|$ between neighbouring energies; global assignment that maximises $\sum O$ (Hungarian algorithm). A new segment starts when $N_{\rm proj}$ changes. | Diagnostic only. |
| 15 | **Blind pole classifier for FG** (v31x/v31y) | `OLD/v31y_…/README.md` | Uses $\min\|\lambda\|$, $\max\|\lambda\|$ growth (≥1.05), norm growth and dimension-jump rejection. Classifies $\det(F_2+G)$ sign flips without the non-interacting list (that list is used only as an oracle afterwards). | — |
| 16 | **digonto_classifier_v1** (v31zo) | `OLD/v32_…/README_v31zo_…md` | Accept if $\det$ flips (either direction) **and** $\min\|\lambda\|$ and $\sigma_{\min}$ both show a "Gaussian-like local peak". | — |
| 17 | **v2_refine / "norm-refine"** (v32r–v32w) | `Q3V3/main_fitter/v33f_…/source/qc_fitter_norm_refine_v2.cpp:144, 280-383` | "Norm" = divide QC by $(L\xi)^6$ (`:85-101, 445`). Each bracket gets `v2_refine_points`=10 new energies; brackets with several flips are split (depth ≤8). Zero if near-flip \|y\| ≤ `zero_ratio`(0.8)·\|y\| at the bracket edge on either side, otherwise pole. | Old production default (`:34`). Builds new rows at every FCN evaluation. |
| 18 | **digonto_v3_window** (coarse only) | `qc_fitter_norm_refine_v2.cpp:147-226`; spec `V33N/docs/DIGONTO_CLASSIFIER_V3_ALGORITHM.md` | Six points $i-2..i+3$ rescaled by the window max. A shoulder is zero-like if each step toward the flip drops by ≥15 %, pole-like if it rises by ≥15 %. 3→2 trimming at a local extremum. True zero if one zero-like shoulder and not both pole-like. | **L20/100_A2: TP 2, FP 0, FN 0, TN 23** against user labels (`V33P/output/…/ACCEPTED_CLASSIFIER_DECISION_E026310_0360.md:9-14`) |
| 19 | **digonto_v4_window** | `qc_fitter_norm_refine_v2.cpp:237-278` | v3 plus merging of true zeros within 1e-3 (keeps the smallest \|y\| as representative). | Same TP/FP on L20/100_A2. In v33m, "COVERAGE_ONLY_NOT_EQUIVALENT", χ² 0.0848 vs 0.2544 (`v33m…report.md:19-25`). |
| 20 | **Sweep algorithms** raw_sign_only / raw_sign_clustered / v6 branch-count / v7 eigenbranch / v7 hybrid | `L20P/scripts/plot_classifier_algorithm_sweep_v33l.py:414-650` | Scores: `det_zero_score`/`pole_score` from log\|y\| at the edges vs the centre. Eigenbranch flag requires a crossing, overlap ≥0.35, condition number ≤1e14 and \|λ_c\|≤1e-6. Candidates are assigned to **lattice target levels** at minimum cost. | v33m chose `raw_sign_only` for production (32/32, χ²=0.2544). Assignment to lattice levels guides this choice, so it is **not** a blind root finder. |
| 21 | **Human labels** (current acceptance gate) | `L20P/input/v33p_*labels*.csv`; `V33P/docs/NEXT_VALIDATION_PLAN_…md:22` | The user marks every sign-change bracket true/false. A sector is accepted only if the classifier has FP=0 and FN=0. | 19 accepted windows for L20 up to 0.36 (`L20P/FULL_L20_FIT_AND_VALIDATION_REPORT.md:36-38`) |

### 5.2 How well things worked (numbers from reports)

- Fraction of determinant sign-change brackets (dimension jumps already excluded) that the user labelled true: counts read from the label CSVs in `L20P/input/` (`python3 -I` count of `user_label`):

| sector | brackets | true | false |
|---|---:|---:|---:|
| L20/000_A1m | 10 | 3 | 7 |
| L20/100_A2 | 25 | 2 | 23 |
| L20/110_A2 | 36 | 4 | 32 |
| L20/111_A2 | 30 | 5 | 25 |
| L20/200_A2 | 35 | 5 | 30 |
| L24/000_A1m | 20 | 3 | 17 |
| L24/100_A2 | 52 | 9 | 43 |

  Across these 208 brackets, 31 are true (15 %). **About 85 % of sign flips are not true zeros.**
- The v3/v4 window classifiers were scored against complete labels only for L20/100_A2 (TP 2/FP 0/FN 0/TN 23). The other sectors are labelled, but their scoring is still listed as pending in `V33P/docs/NEXT_VALIDATION_PLAN_ALL_L_IRREPS_CACHE_ANCHORED.md`.
- In the v33m multi-L benchmark all seven classifier modes reached 32/32 assigned levels. Only `raw_sign_only`, v7 eigenbranch, v7 hybrid and v3 gave the same assignment (χ²=0.25435). v4, clustered and v6 gave χ²=0.08475 with a different assignment (`v33m…report.md:17-25`). Per-FCN times differ by less than 5 % (1.37–1.43 s).

### 5.3 Search result: gradient-, slope- or derivative-based root identification

After removing content duplicates, these are **all** the places where a slope or derivative of the scanned function is used. Every other hit is unrelated: Gauss–Legendre Newton steps in `threebody_fullmatrix…/include/threebody_contour.hpp:69-84`, finite-difference χ² sensitivities in `L20_FULL_CUTOFF_DEPENDENCE_FINAL_REPORT.md:67-68` and `K3df_minuit_fit_v24.hpp:562-590`, and false matches on `-DEIGEN` or `output/debug`.

1. **Slope-ratio pole flag**, `file_processor.hpp:135-145, 258-276` (also `QC3_release/file_processor.hpp`, threshold 1e4) and the v7–v10/v21–v22 copies (`slope_ratio_threshold`=1e2 by default, 1e1 at the v7/v8 call sites). It computes $s_i=|\Delta y_i/\Delta E_i|$ and $\max s/\min s$ (or the median) in a ±5-point window. A large ratio means "likely_pole". It is not a zero detector. The v9 notes record it as dropped for depending on thresholds.
2. **Crossing-direction rule**, v31a/v31b (`debug_v31b…:158-163`). The sign of the discrete derivative of the closest-to-zero eigenvalue across the bracket decides: $\Delta\lambda<0$ is a zero, $\Delta\lambda>0$ is a pole. This is the most direct "gradient" root classifier in the code base. It assumes every zero is crossed with negative slope and every pole jumps from $-\infty$ to $+\infty$. Nothing in the packages tests or justifies that assumption, and no outputs are packaged.
3. **Monotonicity (gradient-sign) tests on |y|**, in methods 3, 7, 17, 18, 19. In discrete form they all test $\operatorname{sign}\frac{d|y|}{dE}$ on each side of the bracket: a zero has $|y|$ decreasing toward the flip from both sides (a V-shaped minimum), and a pole has $|y|$ increasing (a divergent peak). The best-scoring method (v3, TP 2/FP 0/FN 0/TN 23) is a three-point gradient-sign test with a 15 % minimum step change. **In effect the "gradient" idea is already the production classifier, applied to $|\det|$ on the coarse grid.**
4. **Least-squares side slopes** `slope_left`/`slope_right`/`slope_ratio` in `L20P/scripts/plot_classifier_algorithm_sweep_v33j.py:447-455, 507-520` (and v33i `:432-488`). These are computed and saved but **not used** in `zero_score`/`pole_score`, so they are diagnostic only. **Bug** in the same block, v33i `:489` and v33j `:521`: `y_est = y[i] - y[i]*(e[i+1]-e[i])/(y[i+1]-y[i])` uses `y[i]` where `e[i]` belongs, so the stored `E_estimate` is wrong in those two scripts. This bug is copied into the V33N, v33o and V33P packages.
5. **Secant method** (method 0). It uses a finite-difference slope to move toward the root, but it is not a classifier.
6. **Not found anywhere**: an analytic or finite-difference $dE$-derivative of $\det$, $\log|\det|$ or the eigenvalues; Hellmann–Feynman eigenvalue derivatives; any second derivative, curvature or inflection test; or the phrase "gradient change". The user's idea has not been written down in any package.

**Proposal (my own, not in the code):** a derivative test that does not depend on thresholds. By Jacobi's formula,
$$
g(E)\equiv\frac{d}{dE}\ln\det M(E)=\operatorname{Tr}\!\big[M^{-1}(E)\,M'(E)\big]=\sum_k\frac{\lambda_k'(E)}{\lambda_k(E)},\qquad M=(L\xi)^{-6}V^\dagger(F_3^{-1}+\mathcal K)V .
$$
Near a zero of order $m$, $g\approx +m/(E-E_0)$. Near a pole of order $m$, $g\approx -m/(E-E_p)$. So the **sign of the jump of $g$ (its "gradient change")** separates zeros ($g$ goes $-\infty\to+\infty$) from poles ($+\infty\to-\infty$). It also counts even-order zeros, which produce no sign flip. $M'$ can be taken by central differences between neighbouring cached rows inside one $N_{\rm proj}$ segment. Conditions: `Vsel` must be continuous inside the segment (eigenvectors returned by the eigensolver may rotate within degenerate subspaces), and the spacing $\Delta E=4.85\times10^{-6}$ must be small compared with the distance between a zero and the nearest pole. This needs validating against the human labels in §5.2 before anyone relies on it.

---

## 6. Why the roots are hard: known pitfalls (Question 4)

Labels P-x are mine, used only for reference.

- **P-1. Poles look like zeros in a sign scan.** $\det M$ changes sign at a simple zero **and** at a simple pole of one eigenvalue (a zero eigenvalue of $F_3$ gives a divergent eigenvalue of $F_3^{-1}$). In the labelled sectors ≈85 % of sign changes were not accepted as true zeros (§5.2). `VALIDATION_AND_DEBUGGING.md:84-95` (`OLD/…v31u…/`) lists the false-zero types: "true smooth zero, pole/asymptote, basis/shell dimension jump, closest-eigenvalue branch switch, numerical underflow or overflow artifact".
- **P-2. Dimension jumps.** As $E$ grows, more spectator momenta pass the cutoff, so $N$ and $N_{\rm proj}$ jump. L20/100_A2 has **11 jumps** in $[0.2631,0.36]$, from 34/11 to 108/30 (`cache_window_audit…md:20-33`). $\det$ of matrices of different size are not comparable, and the matrix-level $(L\xi)^{-6}$ scaling changes $|\det|$ by $(L\xi)^{-6\Delta N_{\rm proj}}$ with $(L\xi)^6=1.068\times10^{11}$ at $L/a_s=20$. All current tools drop brackets that straddle a jump: `V33N/CLASSIFIER_RESTART_NOTES.md:31-33`, v31b `:142`, v30s "within the same projected dimension segment". A true root at a jump edge would be missed.
- **P-3. Non-interacting (free) three-particle energies.** $F_2$ and $G$ have singular features at the free energies (v31r–v31y "FGiso pole" studies). Those energies also correlate with determinant spikes, dimension jumps and numerical instability (`VALIDATION_AND_DEBUGGING.md:39-41`). Free-shell positions are used only as an "oracle" overlay, never as input to the blind classifiers (v31y README).
- **P-4. Closest-to-zero eigenvalue switching.** Following "the eigenvalue nearest zero" jumps between branches and creates sign changes that are not zeros (v30z notes; `CONVERSATION_CONTEXT.md:159`). Overlap tracking (v30s) or global assignment (v31l Hungarian) reduce this but are diagnostic only.
- **P-5. Small eigenvalue without singularity.** v30s accepted points where $|\lambda|$ was small but QC was not truly singular. v30t added the residual and $\sigma_{\min}$ checks.
- **P-6. Underflow.** Raw determinants underflowed to exactly 0 for $\ell\le1$ (v31t, `CONVERSATION_CONTEXT.md:226-227`). Use LU-based slogdet (sign × log|det|), as v31u, v33h (`:244-269`) and the fitter `det_info` do.
- **P-7. Even-multiplicity or near-degenerate zeros.** Two eigenvalues crossing zero inside one grid step give no sign change, so they are invisible to every sign-flip method. A zero and a pole closer together than a few $\Delta E$ contaminate the shoulder tests. This is why v2 splits multi-flip brackets and v3 has 3→2 trimming.
- **P-8. Reader convention errors** produce a non-Hermitian $F_3^{-1}$ and a complex $\det$ (§3.6). After the fix the trusted scan reports $\max|\mathrm{Im}\det|=0$ (`cache_window_audit…md:14`).
- **P-9. Projection convention.** Legacy pipelines build $V^\dagger F_3^{-1}V$ (project after inverting). The v29–v30 debug tools build $(V^\dagger F_3V)^{-1}$. These agree only if $F_3$ commutes with $P_I$, which is why leakage checks exist (`max_fv_projection_leakage_rel_full_F3 ~1e-13`, README v30r). Wrong parity (+1 for moving-frame A2, v18) empties the irrep.
- **P-10. Only $\mathrm{Re}\det$ is scanned.** This is correct only if $M$ is Hermitian to rounding error. v33h enforces this with `--hard-hermiticity-check`.

GPU-generator specific:

- **P-F1.** Singular LU is logged and the row is still written (§4.1). A non-finite `F3inv` can enter the cache silently. The CPU generator drops such rows.
- **P-F2.** Record order follows shape groups, not energy (§3.3).
- **P-F3.** Every row in a batch must have the same $N_{\rm proj}$, otherwise the whole chunk throws (`:1645-1648`). There is no recovery by splitting in the concurrent path.
- **P-F4.** `plm_config`/`klm_config` matrices are transposed between CPU ($5\times N_1$, `generate_cpu_component_cache.cpp:205-211`) and GPU ($N_1\times5$, `F3_gpu_cachegen_v32zu.cu:772-784`).
- **P-F5.** Physics constants, parity and `eig_tol` are compiled in, and parity and `eig_tol` are missing from the sidecars.
- **P-F6.** The GPU `erfi` series is not Faddeeva, and the component sidecar says `validated_against:"smoke_only"`. The L20 hot-window table lists matrix max differences of ≤4.66×10⁻¹⁰ and determinant relative differences up to 4.69×10⁻⁴ (200_A2) with full sign agreement (`L20P/FULL_L20_FIT_AND_VALIDATION_REPORT.md:50-56`). That report does not state what reference those differences are measured against.

---

## 7. Canonical generator and minimal regeneration file set (Question 5)

**Canonical now:** the GPU v33d all-object build of `FCG/gpu_component_cachegen/source/F3_gpu_cachegen_v32zu.cu`.

- Compiled with `V32ZU_ALLOBJ_QC=1 bash scripts/compile_gpu_component_cachegen.sh`, which produces `bin/cachegen_gpu_v32zu_allobj_qc` with `-DV32ZU_ALLOBJ_QC_WRITER=1` and `-lcublas -lcusolver`, `-arch=native`.
- Run with `configs/config_gpu_component_cachegen_sep_f2k2_coarse20000.sh` (separate F2/K2 kernels, $N_c=20000$).
- It writes both the combined raw `*_F3inv_Vsel_gpu.bin` and the nine component caches. Every current reader (v33h scan, v33p/L20 fitter) validates the version string `v33d_gpu_component_cachegen_all_ingredients` (`v33h…:305-306`).
- The CPU generator `FCG/cpu_component_cachegen/source/generate_cpu_component_cache.cpp` is the **correctness reference** (`FCG/codex_report.md:96-101`). It writes only component caches, not the combined raw file.
- The legacy determinant scanners (v3–v24) are superseded for cache generation. Their root-finding logic now lives in the fitter (`qc_fitter_norm_refine_v2*.cpp`).

**Minimal GPU file set** (include closure resolved by script; none missing):

```text
gpu_component_cachegen/source/F3_gpu_cachegen_v32zu.cu
gpu_component_cachegen/source/F3_device_resident_true_batch_v16_custom_batched_getrs.cu   # GPU copy, adds K3 fields; differs from legacy_source copy by 4 lines
gpu_component_cachegen/source/{functions_gpu_config_maker4.cuh, F2_gpu_safe_builder.cuh, G_gpu_safe_builder.cuh,
                               K2_functions_gpu_safe.cuh, projections_gpu_safe.cuh, projections_gpu_v32y_exact_device.cuh,
                               functions.h, spherical_functions.h, real_wigner_d.hpp}   # last 3 identical to legacy_source copies
legacy_source/{K3_functions_2plus1.hpp, QC_functions_v2.h, F2_functions_v2.h, K2_functions_v2.h, G_functions_v2.h,
               Faddeeva.hh, Faddeeva.cc}
k3iso01_swave_cpu_cache_qc_omp_package/source/fv_projector_cartesian_l1_v30m.hpp       # identical to legacy_source copy
gpu_component_cachegen/scripts/{compile,run}_gpu_component_cachegen.sh
gpu_component_cachegen/configs/config_gpu_component_cachegen{,_sep_f2k2,_sep_f2k2_coarse20000}.sh
```

**Minimal CPU file set:** `cpu_component_cachegen/source/generate_cpu_component_cache.cpp` plus, from `legacy_source/`:

```text
K3df_minuit_fit_v32f_fullF3inv_QCfull_cached_classifier.hpp
debug_v30r_svd_validated_cartesian_l1_100_A2_zero_scan.cpp
F3_cpu_openmp_v25_K3QC_cached_core.hpp
functions.h, spherical_functions.h
F2_functions_v2.h, K2_functions_v2.h, G_functions_v2.h, QC_functions_v2.h
projections_v1.hpp, dig_tools.hpp, real_wigner_d.hpp
K3_functions_2plus1.hpp, fv_projector_cartesian_l1_v30m.hpp, lattice_data_covariance_cpp.hpp
Faddeeva.hh, Faddeeva.cc
```

Plus the CPU scripts and config. Note that `debug_v30r…cpp` is #included as a header, and Minuit is turned off with the `V3x*_DISABLE_MINUIT` macros. For the clean repo, the physics core should be split out of the fitter and debug files.

---

## 7b. Reference reader recipe and consolidation recommendations

### 7b.1 Reader for the combined raw file (format of §3.3, natural complex order)

This pseudocode describes the format; it does not come from any package.

```text
open file (binary)
line = read until '\n'; require line == "KKPI_F3INV_VSEL_GPU_V32ZT_EXACT_V32Y 1"
rows = []
while not EOF:
    magic = u64;  require magic == 0x5653325a4f524543
    grid_i, dim1, dim2, n, vdim = i32 x5
    Ecm, En = f64 x2
    F3inv = complex128[n*n]    reshaped column-major (Fortran order), element = (Re, Im) read in sequence
    Vsel  = complex128[n*vdim] reshaped column-major
    rows.append(...)
sort rows by grid_i (equivalently by Ecm); require len(rows) == sidecar.rows_written == coarseN
QC     = F3inv + K3df(Ecm)                 # K3df rebuilt on the same basis (plm/klm from config_maker_4_momentum_first)
Mproj  = Vsel^dagger QC Vsel / (Lbyas*xi)^6
(sign, logabs) = slogdet(Mproj) via partial-pivot LU; also record Nfull=n, Nproj=vdim
```

In NumPy the element array is `np.frombuffer(buf, dtype='<c16')`, which already interprets the bytes as (Re, Im). This matches the writer and avoids the `comp(read(), read())` evaluation-order problem.

Required sanity checks (all already present in v33h in some form):

1. $\|F_3^{-1}-F_3^{-1\dagger}\|_F/\|F_3^{-1}\|_F$ at rounding level. A value near 2 means the complex order is wrong.
2. $\|V^\dagger V-\mathbb 1\|_F$ at rounding level.
3. $\max|\mathrm{Im}\det M_{\rm proj}|=0$ up to rounding.
4. The dimension-jump list is identical to the sidecar or CPU prepass.

### 7b.2 Recommended structure for the clean repo (my suggestions)

- **Generator:** keep one physics core (F2, G, K2inv, F3 algebra, projector) shared by the CPU reference and the GPU path. The GPU path stays a thin batched layer. Remove the debug `.cpp` and fitter headers from the CPU generator's include closure (§7).
- **Writer:**
  - One format, either the component format or a new indexed format with a row table (offset, `grid_i`, `Ecm`, `n`, `vdim`) so readers can seek.
  - Rows written in energy order.
  - Parity, `eig_tol`, scattering parameters and the git commit recorded in the sidecar.
  - Non-finite rows and `info>0` from LU must fail loudly instead of being written.
- **Tests:**
  - (a) Writer→reader round trip on a random complex, non-Hermitian matrix.
  - (b) CPU vs GPU per-row agreement for F2, G, K2inv, F3inv and `Vsel`-projected objects. Compare projected objects, not `Vsel` itself, because eigenvectors are unique only up to a unitary rotation inside degenerate eigenspaces. For $V\to VU$ with unitary $U$, $\det(U^\dagger M U)=\det M$.
  - (c) Grid and dimension-jump table compared against the CPU prepass.
- **Root finding:** treat it as a separate module with a fixed contract:
  - input: (Ecm, sign, logabs, Nfull, Nproj) rows;
  - output: brackets with labels {true_zero, pole, uncertain, dimension_jump}.
  - Keep `digonto_v3_window` as the baseline, since it is the only automated mode with a labelled FP=0/FN=0 result.
  - Benchmark every new method, including the $g(E)$ log-derivative test of §5.3, against the 208 labelled brackets of §5.2 before using it in a fit.

### 7b.3 Compact lineage diagram

```text
secant (pole_searching.h)
  -> sign flip + spike/slope-ratio thresholds (file_processor.hpp, v7 adaptive)
  -> inward-shape |y_l|+|y_r| trend (v9, v16/v17)      -> zoom constant-N (v18) -> exact pair (v20)
  -> side-wise |det| trend, multi-flip (GPU v21, CPU v21-v24 with K3)
  -> eigenbranch matching (v29o-q) -> SVD-validated (v30r) -> eigen-min bisect (v30s/t/y) -> local-min (v30z)
  -> crossing direction of closest eigenvalue (v31a/b)   [gradient-sign rule, unvalidated]
  -> Hungarian tracking (v31l), FG blind pole classifier (v31x/y)  [diagnostics]
  -> digonto_classifier_v1 (det flip + min|lambda| + sigma_min peaks, v31zo)
  -> v2_refine "norm-refine" (v32r-v32w, refined points every FCN)
  -> digonto_v3_window / v4 merge (v33a+, coarse-only shoulder gradient-sign test)
  -> sweep modes raw_sign_only / eigenbranch / hybrid (v33i-l, target-assignment)
  -> human user_label gate (v33p, current acceptance)
```

## 8. Problems, assumptions and open items

1. Nothing was compiled or executed. All behaviour is inferred from source, and all performance numbers come from comments and reports, not my measurements.
2. The claim that GCC evaluates function arguments right to left (§3.6) is an inference that explains the observed swap. It is not stated in the packages. The conclusion that the writer stores natural $(\mathrm{Re},\mathrm{Im})$ comes directly from the writer code and the double-swap comments.
3. No production wall-clock time for coarse20000 is packaged. The 283.67 s / 1000 rows baseline has no stated irrep, $L$ or hardware.
4. The `projection_source` sidecar field points to `generate_swave_cpu_cache_qc.cpp`, which is absent from every package I read.
5. Cache sizes in §3.3 are estimates from the record formula and the L20/100_A2 dimension table.
6. The $g(E)$ log-derivative test in §5.3 is my proposal. Standard mathematics (Jacobi's formula, residues at simple zeros and poles) supports it, but the code base does not contain or validate it.
7. The `y_est` bug in the v33i/v33j sweep scripts affects only their diagnostic `E_estimate` column. It should not be carried into the clean repo.

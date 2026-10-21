# E — `threebody_fullmatrix_hpc_cuda_bench_v3_codex`: what it is and whether it belongs in 3body-qc-v4

> **Research report — consolidation investigation.** Date: 2026-10-21. Git commit: none (written before the repo existed; first commit of `3body-qc-v4`). Prompt in force: `prompts/000-consolidation-and-roadmap.md`. Related: `reports/research/001`–`006`.
> Absolute paths below refer to the unpacked original zips (`/home/digimonk/Projects/three-body-qc-codes/extracted/...`) and to a temporary scratchpad. For each file's new location in this repo, see `docs/PROVENANCE.tsv`.
> The system is two kaons + one pion treated as a 2+1 system (identical kaons). Where this report says "K K-bar pi", read "KKπ"; that wording came from the investigation brief and is not established. See `docs/OPEN_ISSUES.md` #7.


- **Date:** 2026-10-21
- **Package studied (read only):** `/home/digimonk/Projects/three-body-qc-codes/extracted/threebody_fullmatrix_hpc_cuda_bench_v3_codex-main/threebody_fullmatrix_hpc_cuda_bench_v3_codex-main/` (called `PKG/` below). Files dated 21 Aug 2026.
- **Git commit:** none. The extracted folder is not a git repository. Its own record says it was pushed to `https://github.com/digonto10602/threebody_fullmatrix_hpc_cuda_bench_v3_codex`, baseline commit `d965bc0` (`PKG/reports/GITHUB_REPOSITORY.md:3-5`).
- **Duplicate:** `extracted/threebody_fullmatrix_hpc_cuda_bench_v3_codex-main_1/` is identical. `diff -rq` finds only one difference: an extra `graphify-out/cache` folder in the copy studied here.
- **How I oriented:** I ran `graphify query`, `graphify explain "threebody_physics.hpp"` and `graphify god-nodes` from `PKG/` (graph: 886 nodes). After that I read the files cited below. Nothing was built or run.

---

## 0. Summary

1. **Physics.** The package solves the **infinite-volume** three-body ladder integral equation in momentum space for identical scalars of mass $m=1$. It works in the elastic window between the bound-state+particle threshold and the three-particle threshold, with scattering length $a=16$. The unknown is the **full matrix** $d(p_i,k_j)$, not one column of it. It is **not** the finite-volume quantization condition (QC): there is no box size $L$, no finite-volume $F$ or $F_3$, no $\mathcal K_{\rm df,3}$ and no spectrum anywhere in the benchmark code.
2. **Benchmark.** It times the dense complex-double solve $BD=-G$ at $N=6144$ ($B$, $G$, $D$ are $N\times N$). The only valid production numbers are:
   - one-core CPU: $554.7\ \mathrm{s}$ per energy;
   - 8-thread OpenMP: $426.3\ \mathrm{s}$ per energy;
   - GPU-native CUDA: $9.21\ \mathrm{s}$ per energy over 80 energies, on an RTX 3070.

   The report headlines a speed-up of $60.2\times$. That figure compares unlike quantities, and the like-for-like value is about $24\times$ (see §2.4). The adaptive hybrid workflow, the hypothesis the campaign was built to test, was **never benchmarked**: its CPU/GPU agreement check failed at $N=6144$.
3. **For the QC fitter.** No code should go into a finite-volume $\mathcal K_{\rm df,3}$ fitter. What carries over is a set of engineering lessons (§3) and two small validation ideas. Recommendation: a **document pointer only**. A possible later exception is the $F_3^\infty$ routine in `copied_reference/`, for the separate job of turning fitted $\mathcal K_{\rm df,3}$ values into infinite-volume amplitudes (§6).

---

## 1. What physics does it solve?

### 1.1 Ingredients (all in `PKG/include/threebody_physics.hpp`)

Conventions: $s=E^2$ is the squared total energy in the centre-of-momentum frame. $\omega_k=\sqrt{k^2+m^2}$ is the energy of a single particle. $\sigma_k$ is the squared invariant mass of the pair left over when the third particle (the *spectator*) carries momentum $k$.

- **Pair invariant mass**, at `:64-66`:
$$\sigma_k=(\sqrt s-\omega_k)^2-k^2 .$$

- **Two-body s-wave amplitude, scattering-length-only**, at `:77-82`. Here $\mathrm{mysqrt}(x)=i\sqrt{-x}$ (`:23-25`), which puts the square-root branch cut on the positive real axis:
$$\mathcal M_2(\sigma)=\frac{16\pi\sqrt{\sigma}}{-1/a-i\,\mathrm{mysqrt}\!\left(\tfrac{\sigma+i\epsilon}{4}-m^2\right)} .$$
  For $a>0$ this has a two-body bound state (a "dimer", the $b$ in $\varphi b$) at
$$\sigma_b=4\left(m^2-\frac{1}{a^2}\right)\quad(\text{`:68-70`}).$$

- **One-particle exchange (OPE) projected to s-wave**, at `:85-92`. $H$ is a smooth cutoff built from the function $J$ (`:35-49`):
$$G(p,k)=-\frac{H(\sigma_p,\sigma_k)}{4pk}\,\ln\frac{\alpha_{pk}-2pk+i\epsilon}{\alpha_{pk}+2pk+i\epsilon},\qquad \alpha_{pk}=(\sqrt s-\omega_p-\omega_k)^2-p^2-k^2-m^2 .$$
  A second-sheet version adds $-\frac{H}{4pk}\,2\pi i$ (`:95-102`).

- **Kernel**, at `:126-129`:
$$K(p,k)=\frac{k^2}{(2\pi)^2\,\omega_k}\,G(p,k)\,\mathcal M_2(\sigma_k).$$

### 1.2 The equation and how it is turned into a matrix problem

The code solves the ladder equation for the amputated amplitude $d$:
$$d(p,k)=-G(p,k)-\int_0^{q_{\max}}\!\frac{dq\,q^2}{(2\pi)^2\,\omega_q}\,G(p,q)\,\mathcal M_2(\sigma_q)\,d(q,k),\qquad q_{\max}=\frac{s-m^2}{2\sqrt s}.$$

It discretises the integral with Gauss–Legendre nodes $k_j$ and weights $w_j$ (the *Nyström* method: replace the integral by a weighted sum over quadrature nodes). That gives
$$B\,D=-G,\qquad B_{ij}=\delta_{ij}+w_j\,K(k_i,k_j),\qquad G_{ij}=G(k_i,k_j),\qquad D_{ij}=d(k_i,k_j).$$

Where this happens:
- assembly: `PKG/apps/threebody_fullmatrix_cpu_serial.cpp:82-125`;
- one LU factorization (Gaussian elimination with partial pivoting), then a solve with $N$ right-hand sides: `:135-141`.

`PKG/docs/CODEBASE_ANALYSIS.md:136-151` states the purpose of the package. The upstream production code solved for a single column, $d(p_i;q_b)$, where $q_b$ is the spectator momentum when the other two particles form the bound state. This package instead benchmarks the **full matrix** $D$.

### 1.3 Parameters and energy window

Defaults: $a=16$, $m=1$, $\epsilon=10^{-10}$, complex double precision (`PKG/benchmarks/N6144/build_info.txt`). The window runs between the $\varphi b$ threshold (particle plus dimer, `phibthreshold`, `:164-167`) and the three-particle threshold $9m^2$:
$$m_b=2\sqrt{1-1/a^2}=1.99609,\qquad s_{\varphi b}=(m+m_b)^2=8.97655,\qquad s_{3\varphi}=9 .$$
Source: `PKG/benchmarks/N6144/energy_convention.txt`. The production energies come from `generate_energy_grid` (`:194-208`), which places $E$ evenly between the two thresholds and squares it.

The earlier v2 runs used $s=80.78$. That value comes from the default `s = 8.988*8.988` (`cpu_serial.cpp:293`), which lies far above the three-particle threshold. The Codex audit therefore threw those runs out (`PKG/reports/CODEX_INITIAL_AUDIT.md:10`). **All timings in `PKG/docs/FINAL_REPORT.md` and `PKG/docs/BENCHMARK_REPORT.md` were taken at $s=80.78$ on a different machine, so they are historical only.**

### 1.4 Which literature this matches (my identification — the package cites no paper)

A search for arXiv/Jackura/Dawid/Hansen/Briceño/Romero/Islam finds no citation anywhere in `PKG/`. Judging from the equations and from the upstream folder names (`threebody_ac_efimov_production`, `Efimov_physics`), this is the integral-equation framework of:

- A. Jackura, R. Briceño, S. Dawid, Md H. E. Islam, C. McCarty, *Solving relativistic three-body integral equations in the presence of bound states*, PRD 104, 014507 (2021), arXiv:2010.09820;
- S. Dawid, Md H. E. Islam, R. Briceño, *Analytic continuation of the relativistic three-particle scattering amplitudes*, PRD 108, 034016 (2023), arXiv:2303.04394;
- the follow-up Efimov-states study by the same group, arXiv:2309.01732.

So the answer to "is it the infinite-volume ladder equation for the $d$ amplitude?" is **yes**. The benchmark covers only the $\mathcal K_{\rm df,3}=0$ ladder part. The one place a three-body K matrix appears is in `copied_reference/` (§4.2), which the benchmark never calls.

### 1.5 Caveats on the physics as benchmarked (problems I found)

- **No contour deformation in the benchmark.** The contour in `PKG/include/threebody_contour.hpp:95-138` is the straight real segment $[0,q_{\max}]$. Every node is tagged `OPESheet::First` (`:131`). The upstream code has deformed contours (`copied_reference/contour_planner.hpp`: `GlockleDeformed`, `PhiBSheetIIDimerBypass`), but they were not carried over. The parabolic option is a fixed bulge of 10% of $q_{\max}$ (`:117-119`).
  - Consequence: at real $s$ with $\epsilon=10^{-10}$, the dimer pole of $\mathcal M_2$ at $k=q_b$ and the logarithmic singularities of $G$ lie on or next to the integration path. Without an $\epsilon\to0$ extrapolation or a deformed contour, $D$ is **not a converged physical amplitude**. It is a realistic *workload*, with the right matrix structure and conditioning: a condition estimate of $2.85\times10^4$ at $N=6144$, $s=8.99$ (`benchmarks/N6144/validation_results.csv`).
- **Second-sheet branch is dead code.** In `cpu_serial.cpp:102-108` both branches call the same first-sheet `kernel_pk_2eps`. Because no node is ever tagged `Second`, the second-sheet logic is never exercised at all.
- **CPU and GPU disagree at $N=6144$.** At $s=8.99$ the relative difference is $7.8\times10^{-10}$ in the matrices and $3.4\times10^{-9}$ in the solution. The required matrix tolerance was $10^{-11}$, so the agreement gate is **BLOCKED**. At $N=64$ the same check passes ($6.8\times10^{-14}$ and $1.8\times10^{-13}$). Sources: `benchmarks/N6144/equivalence/physics_equivalence.csv`; `reports/workflow8_equivalence_N6144.log`.
  - My unverified hypothesis: rounding differences between the device and host `sqrt`/`log`/`exp` get amplified near the nearly singular dimer pole. A denser grid at larger $N$ puts nodes closer to that pole.
- **Dead variable.** `gpunative_fixed.cu:335` computes `q_max = sqrt(s-1)` and never uses it. The real contour comes from the host `generate_contour`.

---

## 2. What it benchmarks

### 2.1 Workflows (`PKG/reports/CODEX_INITIAL_AUDIT.md:14-23`)

| # | File (`PKG/apps/`) | What it does |
|---|---|---|
| 1 | `threebody_fullmatrix_cpu_serial.cpp` | single core; Eigen `PartialPivLU` |
| 2 | `threebody_fullmatrix_cpu_openmp.cpp` | OpenMP over rows during assembly (`:93`); the Eigen solve is not threaded |
| 3 | `threebody_fullmatrix_omp_energies.cpp` | OpenMP over energies (`:176`) plus OpenMP assembly (`:81`) |
| 4 | `threebody_fullmatrix_cuda_sequential.cu` | CPU assembly, H2D copy, cuSOLVER `Zgetrf`/`Zgetrs` (`:116,131`) |
| 5 | `threebody_fullmatrix_cuda_streams.cu` | as 4 with a pool of streams and handles (`:262-275`) |
| 6 | `threebody_fullmatrix_cuda_gpunative{,_fixed}.cu` | matrices built on the GPU by `construct_BG_kernel` (`_fixed.cu:195-218`), then `Zgetrf`, `Zscal`, `Zgetrs` on device-resident buffers (`:385-390`) |
| 7 | `threebody_fullmatrix_cuda_pipelined.cu` | OpenMP CPU assembly of a batch, then one OpenMP thread per GPU stream context (`:313`, `:322-331`) |
| 8 | `threebody_fullmatrix_cuda_hybrid.cu` | splits energies between CPU and GPU using a `cpu_fraction` heuristic (`:59`, `:380-412`) |

H2D and D2H mean host-to-device and device-to-host memory copies. `debug_gpu_native.cu` is a tool that compares GPU and CPU on single samples.

### 2.2 Matrix size and memory

$N=6144$ was required (`PKG/prompts/...HANDOFF.md:245-266`). One $N\times N$ complex-double matrix takes
$$N^2\times16\ \mathrm{B}=603{,}979{,}776\ \mathrm{B}=576\ \mathrm{MiB}$$
(`PKG/reports/N6144_GPU_MEMORY_MODEL.md`). Workflow 6 keeps two matrices on the device: $B$, which the LU overwrites, and $G$, which becomes $D$. Total reported: $1152.42\ \mathrm{MiB}$ including workspace (`benchmarks/N6144/memory_summary.csv`). The model predicts at most 5 such slots on the 8 GiB card.

### 2.3 Measured results (only rows that are valid)

Hardware: Intel i7-12700KF, RTX 3070 with 8 GiB, CUDA 13.0 from HPC SDK 25.9 (`benchmarks/N6144/system_info.txt`). Each figure comes from a single run; no repeats, so no uncertainty can be given. Sources: `benchmarks/N6144/raw_timings.csv`, `stage_breakdown.csv`, and the logs in `PKG/reports/`.

| Method | $N$ | energies | wall time | per energy | assembly | factor | solve with $N$ RHS | residual |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| CPU serial (1 thread) | 6144 | 1 | 554.74 s | 554.74 s | 9.57 s | 51.63 s | 158.07 s | $2.15\times10^{-14}$ |
| CPU OpenMP (8 threads) | 6144 | 1 | 426.35 s | 426.35 s | 1.49 s | 10.89 s | 161.59 s | $2.15\times10^{-14}$ |
| CUDA GPU-native, fixed (WF6) | 6144 | 80 | 737.18 s | 9.215 s | 0.486 s | 8.729 s (factor + solve) | (included) | not emitted |
| CUDA streams (WF5) | — | — | skipped | — | — | — | — | — |
| Adaptive hybrid (WF8) | — | — | **blocked** | — | — | — | — | — |

Headline claims as reported: WF6 is $60.2\times$ the serial throughput and $46.3\times$ the OpenMP throughput (`benchmarks/N6144/speedup_summary.csv`). `scripts/check_benchmark_arithmetic.py` confirms that throughput $=$ energies $/$ wall time.

### 2.4 My corrections to how the numbers read (derived from the CSVs above)

- **The headline speed-up includes extra CPU work.** The CPU "total" contains a residual product $B\cdot D$ and a condition estimate $\|B\|\,\|B^{-1}\|$ that forms $B^{-1}$ explicitly (`cpu_serial.cpp:148-164`). Together these cost $554.74-219.27=335.5$ s. The GPU number includes neither, and does not copy $D$ back to the host either.
  - Like-for-like serial (assembly + factor + solve): $9.57+51.63+158.07=219.27$ s, so speed-up $219.27/9.215\approx23.8$.
  - Like-for-like OpenMP: $173.97$ s, so speed-up $173.97/9.215\approx18.9$.
- **The OpenMP gain is small because the multi-RHS triangular solve stayed single-threaded.** The solve took 161.6 s versus 158.1 s serial.
- **The GPU runs close to the card's double-precision limit.** Counting operations:
  - complex LU: $\tfrac{8}{3}N^3\approx6.2\times10^{11}$ real floating-point operations (flop);
  - two triangular solves with $N$ right-hand sides: $8N^3\approx1.9\times10^{12}$ flop;
  - total $\approx2.5\times10^{12}$ flop in $8.73$ s, i.e. $\approx0.28$ TFLOP/s.

  The RTX 3070's FP64 peak is about $0.32$ TFLOP/s (1/64 of its FP32 rate). That is about 88% of peak (my estimate). A data-centre GPU with strong FP64 should therefore be one to two orders of magnitude faster. That is projected, not measured. Note also that the $N$-column solve costs about 3× the factorization; asking for one column instead of all $N$ would make the solve almost free.
- **The Nsight "memcpy-bound" reading is misleading.** Nsight is NVIDIA's profiler. It attributes 91.2% of CUDA API time to `cudaMemcpy`, but only 0.418 MB was moved (`benchmarks/N6144/profiling/nsys/...summary.txt`). The time is the blocking copy of the 4-byte `d_info` status flag (`gpunative_fixed.cu:386`), which waits for the LU to finish. It measures time spent waiting for the GPU, not data transfer.
- **The "streams" workflow (WF5) never overlaps anything.** It runs a serial host loop, shares one set of `d_B/d_G/d_D` across all streams (`streams.cu:280-284`, `:311-330`), calls `cudaStreamSynchronize` after every step (`:179-219`), and copies from pageable `std::vector` memory. Workflow 7 also copies from pageable memory (`pipelined.cu:228-240`). Neither uses pinned host memory or batched cuSOLVER calls.
- **The hybrid workflow (WF8) is not adaptive.** It splits work by a fixed table indexed on $N$, clamped to $[0.1,0.9]$ (`hybrid.cu:380-412`). A small-scale smoke test assigned a GPU job that reported zero solve time (`reports/CODEX_WORKFLOW8_N6144_DEN200_FINAL.md`; `workflow8_scheduler_trace.csv`, labelled `ACCOUNTING_FAILURE`). The measured-rate scheduler exists only as a plan (`reports/WORKFLOW8_ADAPTIVE_SCHEDULER.md`).

### 2.5 `release/`

`PKG/release/` contains only `SHA256SUMS.txt`: the checksum of a `threebody_fullmatrix_hpc_cuda_benchmark_results.zip` that is not in this extract. There is no further result here.

### 2.6 The package's own conclusion

From `reports/CODEX_N6144_FINAL_BENCHMARK_REPORT.md`: the GPU-native workflow is the fastest measured path. The hypothesis that the adaptive hybrid would win was not established. Release gates: G3 (numerical agreement) PARTIAL; G8 (Nsight Compute, NVIDIA's kernel-level profiler) BLOCKED by `ERR_NVGPUCTRPERM`, a permissions error; G9 (hybrid) PARTIAL.

---

## 3. Techniques worth reusing in a fast finite-volume QC fitter

The regime is different. A finite-volume $\mathcal K_{\rm df,3}$ fit needs **many small matrices**: one per energy and per set of $\mathcal K_{\rm df,3}$ parameters. Their size is set by the finite-volume momentum shells and the $\ell m$ (angular-momentum) indices, typically tens to a few hundred. This package has **one large matrix per energy**. Its own historical small-$N$ data, though invalid as physics, consistently show the GPU losing at $N\lesssim128$, mostly to setup cost. Examples: CUDA sequential at $N=64$ runs at $0.005\times$ serial; `cudaStreamCreate` is 77–83% of API time (`docs/FINAL_REPORT.md:114`; `docs/CHANGES_SUMMARY.md:278-285`).

| Technique | Where it is in PKG | Measured evidence | Use in a QC fitter |
|---|---|---|---|
| Matrices built on the GPU, one thread per entry $(i,j)$ | `gpunative_fixed.cu:195-218` | $0.486$ s per energy for $2N^2=7.5\times10^7$ entries, vs $9.57$ s serial and $1.49$ s on 8 cores (§2.3) | Directly reusable as a pattern: build the $\mathcal K_{\rm df,3}$-independent finite-volume matrices for every energy in one kernel |
| Device-resident buffers, persistent cuSOLVER/cuBLAS handles and workspace, no `cudaMalloc` inside the loop | `gpunative_fixed.cu:309-322` | Lets one buffer set run 80 energies with flat memory use ($1152$ MiB) | Essential: allocate once per fit, reuse every energy and iteration |
| Factor in place, overwrite the right-hand side with the solution (`Zgetrf`, then `Zscal(-1)`, then `Zgetrs` into $G$) | `gpunative_fixed.cu:385-390` | Halves device memory compared with a separate $D$ buffer | Same pattern for $\det$ or solves at small $N$, but via **batched** calls (`cublasZgetrfBatched` or `cusolverDn` batched routines). Batched calls are not used anywhere in PKG |
| Memory model computed before launch, stream slots limited by available VRAM | `reports/N6144_GPU_MEMORY_MODEL.md` | Predicts 5 slots | Reusable as a written check |
| Stream pool with one handle per stream | `streams.cu:262-275` | **No measured gain** (serial loop, shared buffers, pageable memory) | Lesson: streams help only with separate buffers, pinned host memory and no per-step synchronisation |
| CPU/GPU split ("hybrid") | `hybrid.cu` | **Never measured at production size**; static heuristic | No evidence it pays off; skip |
| Throughput-arithmetic checker | `scripts/check_benchmark_arithmetic.py` (19 lines) | Passes | Worth reproducing: it enforces throughput $=$ energies $/$ wall time in every benchmark CSV |
| Energy-convention test (print thresholds, reject points outside the window) | `scripts/check_energy_convention.py`, `benchmarks/N6144/energy_convention.txt` | Caught the invalid $s=80.78$ | Worth reproducing for the QC: record whether $s$ or $E$ is used, the frame, and the threshold values |
| Separate CPU/GPU agreement check on sampled $B$, $G$, $D$ (`--validate`) | `gpunative_fixed.cu:354-396` | Exposed the $10^{-9}$ disagreement at $N=6144$ | Required by the user's CPU/GPU agreement rule; copy the idea |

**Not reusable:** the measured speed-ups themselves. They concern one $6144\times6144$ dense LU per energy on a consumer GPU with weak FP64, which says nothing about many small finite-volume matrices. A general point not taken from this package: in a $\mathcal K_{\rm df,3}$ fit, everything that does not depend on $\mathcal K_{\rm df,3}$ (for example $F_3$ per energy) can be computed once and cached across fit iterations. That caching is likely to matter more than any of the GPU tricks above. The sibling packages `full_cache_generator_package_20260629-main` and `kkpi-v33n-corrected-gpu-cache-reader-main` appear, from their names, to address exactly that; I did not inspect them.

---

## 4. `copied_reference/`: where it comes from and what it does

### 4.1 Origin

According to `PKG/docs/REPOSITORY_ISOLATION_REPORT.md:30-41` (list of copied files) and `:9-29` (list of files read), these are copies from Digonto's own upstream repository on another machine. The upstream root is
`/home/digonto/Codes/Github_repos/3body_integral_analytic_package/threebody_ac_efimov_production_v2_20260731T054359Z/`.
That package is labelled "V14 production" in `docs/CODEBASE_ANALYSIS.md:5-6`. It is an analytic-continuation / Efimov production code.

- `functions_momrep_based.h`, `integralequation_momrep_based.h`: from `software/src/reference/reused/`.
- `solvers_gpu.h`, `d_secondsheet_printer.cpp`: from `reference_snapshots/Efimov_physics/`.
- `contours.hpp`, `continuation.hpp`, `contour_planner.hpp`, `three_phi_sheets.hpp`, `unitarity_domain.hpp`, `pole_certificate.hpp`, `gauss_legendre.hpp`: from `software/include/threebody/`.

None of these files (checked by name and MD5 checksum) appears in any other extracted package under `/home/digimonk/Projects/three-body-qc-codes/extracted/`. The upstream `3body_integral_analytic_package` is **not** among the extracted packages.

### 4.2 Role

`CMakeLists.txt:90` adds `copied_reference/` to the include path, but no benchmark app includes these files. The graphify `path` query also finds no link from `F3infvol()` to any app. They are **reference only**. The apps use rewritten copies in `PKG/include/`.

| File | Contents (from reading it) |
|---|---|
| `functions_momrep_based.h` (1272 lines) | Original physics functions: `M2kfunc` and its variants at `:283,300,316`, `GS_pk`, kernels, `gfuncConst` |
| `integralequation_momrep_based.h` (414 lines) | Builders for $G$ as a vector or matrix, and for $B$ with weights and sheet tags (`:10-306`); `rho3_func` (`:352-366`); `F3infvol` (`:368-412`) |
| `solvers_gpu.h` (498 lines) | Eigen LU wrappers `LinearSolver{,_2,_3}` (`:17-46`); a cuSOLVER **QR**-based complex solve `cusolverComplex` (`cusolverDnZgeqrf`, `Zunmqr`, `:47-200`) |
| `d_secondsheet_printer.cpp` (560 lines) | Older driver with $a=-8$ (`:28`), $N=500$, $s=8.8+10^{-4}i$; solves for a vector and a matrix (`:72-73`). Build lines in `:1-6` mention the cluster "wahab", CUDA 10 and HPC SDK 2022 |
| `contours.hpp`, `contour_planner.hpp`, `continuation.hpp` | Deformed contour families (`RealAxisRegulated`, `SourceSebaImspos`, `GlockleDeformed`, `PhiBSheetIIDimerBypass`), contour validation, labels for regions of analytic continuation |
| `three_phi_sheets.hpp`, `unitarity_domain.hpp` | Bookkeeping for the sheet of the three-particle cut (winding labels 0, +1); a gate checking that $s$ is in the elastic-unitarity window |
| `pole_certificate.hpp` | A capped-size diagnostic record for poles on unphysical sheets (residue samples, convergence rows) |

**`F3infvol`** is the one piece that touches $\mathcal K_{\rm df,3}$ physics. Write
$$\int_q\equiv\sum_i w_i\,\frac{q_i^2}{(2\pi)^2\,\omega_{q_i}},\qquad \rho_q=\texttt{rho3\_func}(s,q).$$
It returns
$$F_3^\infty(s)=\int_q\left[\frac{\rho_q}{3}-\rho_q\,\mathcal M_2(\sigma_q)\,\rho_q\right]-\int_q\int_{q'}\rho_q\,\mathcal M_2(\sigma_q)\,d(q,q')\,\mathcal M_2(\sigma_{q'})\,\rho_{q'} .$$
This needs the **full** matrix $d(q,q')$, which is why a full-matrix solve was worth benchmarking at all.

The structure matches the infinite-volume $F_3^\infty$ that links an isotropic (constant) $\mathcal K_{\rm df,3}$ to the infinite-volume three-body amplitude, in the form $\mathcal M_{\rm df,3}\sim[1/\mathcal K_{\rm df,3}+F_3^\infty]^{-1}$ "dressed" by ladder factors. That identification is mine (Hansen & Sharpe, arXiv:1504.04248, and the integral-equation papers above). **I have not checked the normalisation or $\rho$ conventions against any paper.**

---

## 5. `prompts/`: what was planned

There is one file: `PKG/prompts/CODEX_THREEBODY_N6144_HPC_CUDA_HANDOFF.md` (922 lines). It is a handoff to Codex (OpenAI's coding agent) to take the v2 repository (`/home/digonto/Codes/Github_repos/threebody_fullmatrix_hpc_cuda_bench_v2`) to v3. Plan, with section numbers:

1. Strict workspace boundary (§0). Correctness before speed; no single precision or TF32 (§1, `:35-57`).
2. Audit first (§2), and settle the $s=E^2$ issue that produced $s=80.78$ (§3, `:98-120`).
3. Publish a new GitHub repository (§4). Build a graphify code graph (§5).
4. Hypothesis (§6, `:226-228`): *an adaptive OpenMP + multi-stream CUDA hybrid gives the highest sustained throughput* for batches of $N=6144$ systems. Do not fake or cherry-pick results.
5. Mandatory $N=6144$ (§7). Batch sizes (number of energies): serial ≤ 3, OpenMP 50, CUDA/hybrid 100 (§8). Calibrate against a time budget and freeze the batch size beforehand (§9). Compare by throughput, labelling projected numbers (§10). Warm-ups (§11). Nsight Systems and Compute (§12-13).
6. Optimise the hybrid (§14, `:530-596`): persistent buffers, pinned memory, a stream pool driven by events, ring buffers, routing to the GPU-native path, a scheduler that is genuinely adaptive. Memory safety at $N=6144$ (§15).
7. Required outputs, plots, final report, README, release gates G1–G11, execution order (§16-24).

Outcome against the plan: steps 1–3 done. Step 5 partial (serial and OpenMP calibrated with only 1 energy each; GPU at 80 energies). Step 6 not done: the hybrid is blocked by the agreement gate. A later attempt at a 200-energy hybrid run also stopped at that gate (`reports/CODEX_WORKFLOW8_N6144_DEN200_FINAL.md`).

---

## 6. Recommendation for `3body-qc-v4` (finite-volume $\mathcal K_{\rm df,3}$ fits)

**Verdict: document pointer only. Copy no code into the core of the fitter.**

Reasons:
1. The physics is infinite-volume (§1). The QC fit itself does not need it.
2. The production code path is not validated. CPU and GPU disagree at $10^{-9}$ for $N=6144$. The benchmark contour is not physical. The second-sheet branch is dead code (§1.5).
3. The performance regime (one huge dense LU per energy) is the opposite of the QC regime (many small matrices) (§3).
4. The good code it contains (`copied_reference/`) is a copy. The authoritative source is the upstream `3body_integral_analytic_package` (V14), which is not in this collection.

**What to put in v4:**
- `docs/RELATED_CODES.md`, a pointer entry containing:
  - this package's GitHub URL and commit `d965bc0`;
  - a one-paragraph description from §1;
  - the measured table from §2.3 with the corrections from §2.4;
  - the lessons table from §3.
- *Re-implement* (not copy) two validation helpers in the v4 test or benchmark tooling:
  - a throughput-arithmetic checker (pattern from `scripts/check_benchmark_arithmetic.py`);
  - an energy-convention and threshold recorder (pattern from `scripts/check_energy_convention.py`).
- Keep the design pattern "assemble on the GPU, one thread per matrix entry, persistent handles, reused device buffers" (`gpunative_fixed.cu:195-218`, `:309-322`) as guidance for a GPU backend. For small finite-volume matrices, use batched LU or determinants instead.

**Possible later, optional module** (only if v4 grows a "from fitted $\mathcal K_{\rm df,3}$ to infinite-volume $\mathcal M_3$" stage):
- take `integralequation_momrep_based.h` (the $B$/$G$ builders and `F3infvol`) and the contour planner from the **upstream** V14 package, not from `copied_reference/`;
- check the $F_3^\infty$ normalisation against the literature before use;
- record its method in `docs/METHODS.md`.

---

## 7. Problems and assumptions

- I built and ran nothing. All numbers are copied from the package's CSVs and logs. The derived numbers in §2.4 are my arithmetic on those values.
- There are no repeat measurements, so no timing uncertainty is available. Serial and OpenMP each have a single energy at $s=8.99$.
- The paper identifications (§1.4, §4.2) are my inference. The package cites nothing.
- The explanation of the $N=6144$ CPU/GPU disagreement (§1.5) is a hypothesis.
- The RTX 3070 FP64 peak ($\approx0.32$ TFLOP/s) is a vendor-level figure I quoted from memory, not something measured here.
- The upstream paths (`/home/digonto/...`) belong to another machine. I did not check whether they exist on this one.

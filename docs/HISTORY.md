# HISTORY — where this code came from

`3body-qc-v4` consolidates eight packages (2026-06-15 → 2026-08-21) plus 71 archived older versions (v31u–v32y). This page is the map. Detail:
- [`reports/research/006-history-v31-to-v32.md`](../reports/research/006-history-v31-to-v32.md): v26–v32y, phase by phase, with lessons;
- [`reports/research/002-…`](../reports/research/002-cache-generation-and-root-finding-lineage.md) §5: the root-finding lineage;
- [`reports/research/003-…`](../reports/research/003-fitter-lineage-and-fast-fcn-design.md) §3: the fitter lineage v32w → L20;
- `docs/original_handoffs/`: the original READMEs and handoff reports, verbatim.

Every file in this repo is traceable to its original location by md5 hash in [`PROVENANCE.tsv`](PROVENANCE.tsv).

## Source packages

| date | package (zip) | what it is | status in v4 |
|---|---|---|---|
| ≤2026-07-15 | `3body_quantization_v3-master` | original repo: S-wave QC core (`QC_functions.h` etc.), Python fitter `QC_fitter_1.py`, K3iso/F3iso/pole makers, GPU experiments; `main_fitter/v33f_…_v3_package` = hot-window FCN checkpoint | old S-wave core and the 726-line `F2_functions_v2.h` are **superseded** (p-wave bugs); main_fitter superseded by L20 |
| 2026-06-15 | `KKpi-v32w-v32y-repro-master` | v32w trusted single-L fitter (χ²/ν = 1.44, 12 levels, L20, cutoff 0.335), v32y memory-bounded spectrum vs L, tools, 71 archived old versions | v32w kept in `legacy/v32w_single_L_fitter/`; tools in `tools/legacy_v32/` |
| 2026-06-29 | `full_cache_generator_package_20260629` | v33d GPU (v32zu) and CPU component-cache generators, `legacy_source/` with ~150 historical sources | generators in `src/cachegen_*`; selected history in `legacy/root_finding_lineage/` |
| 2026-07-05 | `KKpi_v33o_production_fitter` | v33f/v33n production fitter snapshot, cutoff scan 0.315–0.365 (L20+L24, full grid, `raw_sign_only`) | superseded; scan driver kept as `scripts/fit/run_v33n_cutoff_scan_fits.py` |
| 2026-07-09 | `kkpi-v33n-corrected-gpu-cache-reader` | **combined raw cache reader fix** (complex convention), CPU subset validation | v33g runtime-cache builder taken from here (plain read) |
| 2026-07-14 | `KKpi-v33p-cache-anchored-classifier-fitter` | cache-anchored classifier restart, **human-label workflow**, TP/FP scoring, GPU det backend audit | classifier scripts in `scripts/classifier/` |
| 2026-07-15 | `KKpi-L20-user-validated-hotwindow-fitter` | **newest**: projected-basis hot-window cache, L20 pass-only fits at 0.315–0.345, zero tests | **main source of `src/`, `scripts/`, `data/`, `results/`** |
| 2026-08-21 | `threebody_fullmatrix_hpc_cuda_bench_v3_codex` | infinite-volume ladder equation GPU benchmark (not finite-volume QC) | not copied; see `docs/RELATED_CODES.md` |

## Phases (condensed)

1. **Pre-v31 (v7–v30).**
   - GPU/CPU $F_3$ pipelines with sign-flip root searches (v7–v24).
   - Cached Minuit fitter (v25–v26).
   - Eigenbranch, SVD and eigen-min root finders (v29–v30).
   - **v30m/n: corrected Cartesian $\ell=1$ projector.** This is still used.
2. **v31.**
   - Determinant underflow fixed by slogdet (v31u).
   - Projection audits showed the projector and the basis order are not the cause of extra sign flips (v31q, v31w).
   - Blind pole classifier (v31y).
   - Projected $F_3^{-1}$ zero finder (v31z–zd).
   - Toy-mass and scattering-parameter mapping fixes (v31ze–zm).
   - digonto classifier v1 (v31zn–zw).
3. **v32.**
   - Iterative refinement (v32–v32d).
   - $\mathcal K_{\rm df,3}$ fitter with full $F_3^{-1}$ cache (v32e–q). **v32e–p read L24 data at L20: invalid.**
   - Classifier v2 / norm refine (v32r–v).
   - **v32w fixed coarse-grid pollution: refinement points were being merged into the sign-flip grid.**
   - v32x/y spectrum vs L (v32y limits memory).
4. **v33.**
   - Separate F2/K2 GPU kernels and the 20,000-point grid (v33d/e).
   - Multi-L fitter (v33f).
   - Runtime $\mathcal K_3$-basis cache (v33g).
   - Classifier sweeps (v33i–m).
   - Cutoff scan (v33n).
   - Reader fix (v33n).
   - Human labels (v33p).
   - **Hot-window projected-basis cache (L20): FCN 1.4 s → 24 ms.**
5. **After v33.**
   - The CUDA full-matrix benchmark, a separate infinite-volume project.
   - The consolidation into this repo (2026-10-21), which found the covariance issue (`OPEN_ISSUES.md` issue 1).

## Lessons that must not be relearned

- Use slogdet via LU, never a raw determinant.
- Project $F_3^{-1}$; never invert a projected $F_3$.
- Scale the matrix by $(L\xi)^{-6}$; never rescale the determinant by $N_{\rm proj}$.
- A sign change is not a level. Reject dimension jumps first. About 85 % of the remaining sign changes are poles or artefacts.
- Never let refinement points enter the coarse sign-flip grid (the v32w bug).
- The complex read convention of the raw cache is a known trap (`OPEN_ISSUES.md` issue 3).
- Coverage (all levels found) is not correctness. Different classifiers gave 32/32 levels with different assignments and different χ².
- Check the lattice file's $L$ against the computed $L$. Check the `.jack` sample convention (`OPEN_ISSUES.md` issue 1).
- Do not hold all irreps' caches in memory at once (v32x → v32y). Load caches once, outside the FCN.

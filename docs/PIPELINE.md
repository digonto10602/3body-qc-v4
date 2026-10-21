# PIPELINE — from physics parameters to a $\mathcal K_{\rm df,3}$ fit

This document gives the stage-by-stage data flow, the program that does each stage, and the file formats between stages. The equations are in [`METHODS.md`](METHODS.md); the problems are in [`OPEN_ISSUES.md`](OPEN_ISSUES.md).

```text
 (A) cache generation  ── GPU src/cachegen_gpu (production)  /  CPU src/cachegen_cpu (reference)
        │   writes, per (L/a_s, irrep): raw F3inv+Vsel cache  *_F3inv_Vsel_gpu.bin  + .meta.json   [GBs, external]
        ▼
 (B) determinant scan ── src/cache_io/v33h_patched_gpu_cache_oldscale_det_scan.cpp (trusted CPU/OpenMP)
        │   det grid CSV (row, Ecm, Nfull, Nproj, sign, log|det|) + candidate sign-change brackets CSV
        ▼
 (C) labelling        ── scripts/classifier/*  (plot det family, v3/v4 predictions, score vs labels)
        │   human user_label=true/false per bracket  → data/labels/*.csv
        ▼
 (D) accepted windows ── scripts/windows/build_l20_user_validated_windows.py → data/windows/<L>/<irrep>/accepted_windows.csv
        ▼
 (E) hot-window cache ── src/cache_build/v33p_build_projected_basis_hotwindow_cache.cpp
        │   per window row r: F~_r, B~(1..4)_r  (projected, scaled)  → projected_basis_hotwindows.bin + .json
        ▼
 (F) fit              ── src/fitter/qc_fitter_norm_refine_v2_multiL.cpp  (Minuit2; ~24 ms / FCN)
        │   driven by scripts/fit/run_l20_full_cutoff_fits_v2.py (prefit audit → fcn-once → fit → after-fit → reports)
        ▼
 (G) audits/plots     ── scripts/audit/*, tools/audit/*, scripts/plot/*
```

## (A) Cache generation

- **GPU, production:** `./build_gpu.sh cachegen` builds the all-object writer (`-DV32ZU_ALLOBJ_QC_WRITER=1`). Run it with `scripts/cachegen/gpu_original/run_gpu_component_cachegen.sh <Lbyas>` and `configs/cachegen_gpu/config_gpu_component_cachegen_sep_f2k2_coarse20000.sh`. The original scripts assume the package layout, so adapt the paths.
- **CPU, correctness reference:** `./build.sh cachegen_cpu`, config `configs/cachegen_cpu/config_cpu_component_cachegen.in`. It writes component caches only, not the combined file.
- **Hard-coded settings:** physics constants, parity and `eig_tol` are compiled into the GPU binary (`OPEN_ISSUES.md` issue 8).
- **Grid:** $E_i=0.26310+i\Delta E$, $N_c=20000$, up to 0.36, in units of $a_t^{-1}$. The sidecar `.meta.json` is the source of truth for the window.
- **Location of the existing caches (outside the repo):** `/media/digonto/Data/F3inv_cache/v33e_coarse20000_sep_f2k2_L<L>_<irrep>/v32zu_Lbyas<L>_gpu_cache/v32zu_Lbyas<L>_xi3p444_irrep<irrep>_coarse20000_F3inv_Vsel_gpu.bin`

### Combined raw cache format (`*_F3inv_Vsel_gpu.bin`)

| part | type | content |
|---|---|---|
| header | 39 B ASCII | `KKPI_F3INV_VSEL_GPU_V32ZT_EXACT_V32Y 1\n` |
| record +0 | uint64 | magic `0x5653325a4f524543` |
| +8 | int32×5 | `grid_i, N1, N2, N(=N1+N2), Nproj` |
| +28 | double×2 | `Ecm, Elab` |
| +44 | complex128 $N^2$ | $F_3^{-1}$, column-major, bytes $(\mathrm{Re},\mathrm{Im})$ |
| … | complex128 $N\,N_{\rm proj}$ | $V_{\rm sel}$, column-major |

- **Ordering and index:** records are **grouped by matrix shape, not sorted by energy**, and there is no index. Sort by `grid_i` after reading.
- **Complex reading:** the bytes are stored in natural $(\mathrm{Re},\mathrm{Im})$ order. The trusted path reads them through `read_complex` (one swap) plus `swap_real_imag` (a second swap), which cancel. See `OPEN_ISSUES.md` issue 3.
- **Size:** about 2.6 GB for L20/100_A2.
- **Component caches:** format `KKPI_COMPONENT_CACHE_V001`, one file per object (`plm_config, klm_config, K2inv, F2, G, F3, F3inv, Vsel, projected_F3inv`), stored row-major. Full byte layout: report 002 §3.5.

## (B) Determinant scan

- **Build:** `./build.sh scan` → `bin/det_scan_raw_cache_v33h`.
- **Arguments:** `--gpu-cache-root --Lbyas --irrep --Emin --Emax --coarseN --threads`, the fixed K3 parameters, `--old-scaling true`, and `--complex-read-convention variant_04_real_imag_swapped`.
- **What it does:**
  1. validates the sidecar metadata;
  2. sorts the rows by energy;
  3. builds $Q$ and computes the LU determinant per row (OpenMP);
  4. writes the determinant grid, a summary, the dimension-jump list, Hermiticity checks and the worst rows.
- **Experimental GPU alternative:** `src/gpu_det/` (`./build_gpu.sh detlu`). It is not validated against the CPU, so it must not drive fits.

## (C) Labelling (the current gate for "which roots exist")

1. Plot the determinant family: `scripts/classifier/plot_v33p_classifier_det_family_interactive.py --grid … --candidate-file … --n-scale 100 --show`.
2. Get classifier predictions: `scripts/classifier/run_v33p_classifier_sweep_L20_100_A2.py`, with modes `digonto_v3_window` and `digonto_v4_window`.
3. Digonto fills in `user_label` for every bracket. Columns: `bracket_id,E_left,E_right,user_label,user_zero_estimate,comment`.
4. Score the classifier against the labels with `scripts/classifier/score_v33p_classifier_against_user_labels.py` (TP/FP/FN/TN; tolerance $8\times10^{-6}$). A sector is accepted only if FP = FN = 0.

The labels for L20 (5 sectors) and L24 (000_A1m, 100_A2; provisional) are in `data/labels/`. The L20/100_A2 file uses the old grid start 0.26301; it is remapped in step (D).

## (D) Accepted windows

- **Script:** `scripts/windows/build_l20_user_validated_windows.py` keeps only `user_label=true` rows and writes `accepted_windows.csv` per sector (copies in `data/windows/L20/`).
- **Window size:** $[\min(c-h,\ell),\max(c+h,r)]$ with $h$=`window_half_rows`=50, expanded once to $\pm H$=`max_window_half_rows`=250. $\pm50$ rows correspond to about $\pm2.4\times10^{-4}$ in $a_tE$.

## (E) Hot-window projected-basis cache

- **Build:** `./build.sh hotwindow` → `bin/build_projected_basis_hotwindow_cache`. Parallel driver: `scripts/windows/run_v33p_projected_basis_hotwindow_parallel.py`. Halo rows are merged with `scripts/windows/merge_projected_basis_rows.py`.
- **What it reads and writes:** it reads the raw cache through the v33h reader and writes the five projected matrices per row, already divided by $(L\xi)^6$.
- **Self-check:** it compares against direct assembly at the reference $\vec c$; the reported agreement is $\le4.7\times10^{-10}$, with 100 % sign agreement.
- **Cost:** 1006–2511 rows per sector, 300–1570 s per sector at about 1.7 rows/s with 18 threads.
- **Binary format:** magic `V33PHWB1`, `uint32 version=1`, `uint32 count`. Then per row: `uint64 row, double Ecm, int32 Nfull, int32 Nproj, uint32 dim`, followed by 5 column-major `dim×dim` complex128 matrices (`Fproj, B0, B1, BB, BE`).
- **JSON sidecar** (copies in `data/hot_windows_meta/L20/`), required fields:
  - `schema_version v33p_projected_basis_hotwindow_v1`
  - `scaling_convention divided_by_pow_Lxi_6`
  - the SHA-256 of the accepted-window CSV
  - `row_indices`
- **Older alternative:** `src/cache_build/v33g_build_runtime_k3basis_cache.cpp` builds the same five matrices for all 20,000 rows (`*_runtime_K3basis.bin`), unscaled. It was used by the v33g–v33n full-grid fits.

## (F) Fit

- **Build:** `./build.sh fitter` → `bin/k3df_fitter_multiL` (needs Minuit2).
- **Usage:** `bin/k3df_fitter_multiL <config> <mode>`, where mode is one of `fit` (default), `fcn-once`, `spectrum-only`, `benchmark-fcn --repeat N --warmup M`, `sensitivity`, `dump-targets`, `pipeline-test` or `audit-l20-hot-window`.
- **Production config template:** `configs/fitter/L20_hotwindow_base_user_validated.in`, with resolved per-cutoff copies `configs/fitter/L20_hotwindow_Ecm_cut_0p3{15,25,35,45}.in`. These still contain the source machine's absolute paths. Key settings:
  ```
  root_search_mode = accepted_windows      projected_basis_cache_mode = hot_windows
  fallback_full_scan = 0                   classifier_mode = digonto_v3_window
  window_half_rows = 50                    max_window_half_rows = 250
  chi_square_mode = raw_cov_inv            use_lattice_covariance = 1
  lattice_jackknife_dir = …                (→ data/lattice/L20)
  accepted_zeros_file_L20_<irrep>, candidate_brackets_file_L20_<irrep>, det_grid_file_L20_<irrep>
  K3iso0_guess … K3E_guess, *_step, float_*
  ```
- **Hot-window inputs:** in hot-window mode the fitter reads `det_grid_file_*` only for row/Ecm/Nfull/Nproj metadata, plus the hot-window `.bin`. The 3 GB raw caches are never opened.
- **Driver:** `scripts/fit/run_l20_full_cutoff_fits_v2.py` (cutoffs 0.315–0.345, warm-started; prefit audit, fcn-once, fit, after-fit, sensitivity, reports). Pass-only wrapper: `scripts/fit/run_l20_user_validated_passonly.py`. Zero tests (fixing $K_B=K_E=0$): `scripts/fit/run_l20_k3b_k3e_zero_tests.py`.
- **Reference config:** `configs/fitter/v33m_multiL_L20_L24_fullgrid_reference.in`, the full-grid L20+L24 production config, for reference only.

## (G) Audits and plots

| script | role |
|---|---|
| `scripts/audit/audit_l20_user_validated_prefit.py` | level count = window count, and rows present, per sector and cutoff |
| `scripts/audit/audit_l20_chi2_reconstruction.py` | independent χ² recomputation from saved outputs |
| `tools/audit/chi2_covariance_convention_check.py` | recomputes χ² under both `.jack` readings (`OPEN_ISSUES.md` issue 1) |
| `scripts/audit/check_v33p_production_scatter_params.py` | checks that configs carry the production two-body parameters |
| `scripts/plot/L20_spectrum_plotter_model_true_zeros_yaxis_fixed.py` | spectrum vs $L$ with lattice, model and non-interacting levels (needs historical inputs, see `OPEN_ISSUES.md` issue 4) |
| `tools/nonint3body_general_doublecover/` | non-interacting three-particle levels per irrep |
| `tools/legacy_v32/make_qc_state_lines_from_scattered_spectrum.py` | turns scattered $E(L)$ points into continuous state lines |
| `tools/legacy_v32/inspect_v32_f3inv_cache_*.py` | inspectors for the older v32 merged caches |

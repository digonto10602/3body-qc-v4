# C. Fitter lineage and fast-fit caching design (K K-bar pi, RFT three-body QC)

> **Research report — consolidation investigation.** Date: 2026-10-21. Git commit: none (written before the repo existed; first commit of `3body-qc-v4`). Prompt in force: `prompts/000-consolidation-and-roadmap.md`. Related: `reports/research/001`–`006`.
> Absolute paths below refer to the unpacked original zips (`/home/digimonk/Projects/three-body-qc-codes/extracted/...`) and to a temporary scratchpad. For each file's new location in this repo, see `docs/PROVENANCE.tsv`.
> The system is two kaons + one pion treated as a 2+1 system (identical kaons). Where this report says "K K-bar pi", read "KKπ"; that wording came from the investigation brief and is not established. See `docs/OPEN_ISSUES.md` #7.


Date: 2026-10-21. Read-only investigation; nothing was built or run except
`md5sum`/`diff`/`awk`/`grep` and two tiny arithmetic one-liners (listed in Section 9).
No graphify graph exists in any of the packages examined, so raw files were read directly.

## 0. Path abbreviations used in this report

All packages live under `/home/digimonk/Projects/three-body-qc-codes/extracted/<pkg>/<pkg>/`.

| short name | full path (package root) | date (from brief) |
|---|---|---|
| `V32W` | `.../KKpi-v32w-v32y-repro-master/KKpi-v32w-v32y-repro-master/latest_working_versions/v32w_QC_fitter_norm_refine_v2_fixedcoarse` | 15 Jun |
| `V33O` | `.../KKpi_v33o_production_fitter-main/KKpi_v33o_production_fitter-main` | 5 Jul |
| `V33N` | `.../kkpi-v33n-corrected-gpu-cache-reader-main/kkpi-v33n-corrected-gpu-cache-reader-main` | 9 Jul |
| `V33P` | `.../KKpi-v33p-cache-anchored-classifier-fitter-main/KKpi-v33p-cache-anchored-classifier-fitter-main` | 14 Jul |
| `MF` | `.../3body_quantization_v3-master/3body_quantization_v3-master/main_fitter/v33f_k3df_v33e_cache_classifier_v3_package` | 15 Jul |
| `L20` | `.../KKpi-L20-user-validated-hotwindow-fitter-main/KKpi-L20-user-validated-hotwindow-fitter-main` | 15 Jul |

`...` = `/home/digimonk/Projects/three-body-qc-codes/extracted`.

Jargon defined once here:

- **FCN**: the function Minuit minimizes (here: the chi-square as a function of the four
  $\mathcal{K}_{\rm df,3}$ parameters). One "FCN call" = one evaluation of chi-square at one parameter point.
- **Minuit2 / Migrad / Hesse**: the ROOT minimization library; Migrad is its quasi-Newton
  minimizer, Hesse computes the second-derivative matrix to give parameter errors.
- **EDM** ("estimated distance to minimum"): Migrad's convergence measure.
- **Cache row / grid row**: one energy point $E_r$ of a fixed 20,000-point grid in
  $E_{\rm cm}$, for which $F_3^{-1}(E_r)$ and the projector $V_{\rm sel}(E_r)$ were precomputed on GPU.
- **Projector $V_{\rm sel}$**: a rectangular complex matrix (columns = basis vectors of the chosen
  finite-volume irrep, e.g. `100_A2`) that projects the full matrix onto one irrep.
- **Hot window** (the project's own term, not a physics term): a contiguous block of grid rows
  around one human-validated zero of the quantization-condition determinant; only these rows are
  cached in projected form for fitting.
- **True zero / pole / uncertain**: project labels for a sign change of $\det$ that is,
  respectively, a genuine root, a sign change caused by a divergence, or undecided.
- **TP/FP/FN/TN**: true positive / false positive / false negative / true negative counts when
  classifier output is compared with human labels.
- Energies are in temporal lattice units, $a_t E_{\rm cm}$ (the code uses `atmK` $= a_t m_K = 0.09698$, `xi` $=\xi=a_s/a_t=3.444$).

---

## 1. What was asked

Six questions about the fitter: (1) which multi-L fitter source is newest and what changed;
(2) exactly how the fast FCN works (precomputation, determinant assembly, root location, level
matching, chi-square, Minuit setup, masks, timings); (3) the classifier and the human-label
workflow; (4) fit results and blockers; (5) stated next steps, verbatim where possible;
(6) a minimal set of fitter files for a consolidated repo, with best-copy paths.

## 2. Concrete actions taken

1. `md5sum` over all copies of `qc_fitter_norm_refine_v2*.cpp`, `K3df_minuit_fit_v32f_*.hpp`,
   `v33p_build_projected_basis_hotwindow_cache.cpp`, `v33h_*`, `v33g_*`, `gpu_cache_reader_convention.hpp`
   in the six packages, plus the content-hash list `scratchpad/hashes.txt`.
2. `diff` between every distinct copy of the multi-L fitter (results in Section 3).
3. Read the newest fitter (`L20/source/qc_fitter_norm_refine_v2_multiL.cpp`, 2107 lines) almost in full,
   the library it includes (`qc_fitter_norm_refine_v2.cpp`, 452 lines), and the relevant parts of
   `K3df_minuit_fit_v32f_fullF3inv_QCfull_cached_classifier.hpp` (1595 lines),
   `v33p_build_projected_basis_hotwindow_cache.cpp`, `v33g_runtime_k3basis_cache.hpp`,
   `gpu_cache_reader_convention.hpp`, `lattice_data_covariance_cpp.hpp`.
4. Read the driver scripts (`run_l20_full_cutoff_fits_v2.py`, `run_l20_user_validated_passonly.py`,
   `materialize_v33p_hot_windows.py`, `build_l20_user_validated_windows.py`,
   `score_v33p_classifier_against_user_labels.py`, head of `run_v33n_cutoff_scan_fits.py`).
5. Read the handoff docs (`L20/README.md`, `FINAL_SHUTDOWN_HANDOFF_REPORT.md`,
   `FULL_L20_FIT_AND_VALIDATION_REPORT.md`, `L24_NEXT_STEPS.md`, `MF/HOT_WINDOW_FCN_CHECKPOINT.md`,
   `V33N/GPU_CACHE_READER_FIX.md`, `VALIDATION_STATUS.md`, `CLASSIFIER_RESTART_NOTES.md`,
   `V33P/README.md`, `AGENTS.md` key sections, v33n/v33m/v33o reports, v32w README/VERSION_SUMMARY).
6. Read fit outputs: `L20/output/l20_full_cutoff_fits_user_validated_0p36_passonly/*`
   (summary CSV, per-cutoff `fit_levels`, covariance, sensitivity), `L20/output/l20_k3b_k3e_zero_tests/*`,
   hot-window cache JSON metadata, user-label CSVs.
7. Recomputed the spread of one lattice input file with `awk` (Section 6.1) because the reported
   lattice errors looked implausible.

---

## 3. Q1 - Which fitter source is newest, and what changed

### 3.1 Hash/size table of the multi-L fitter `qc_fitter_norm_refine_v2_multiL.cpp`

| package | lines | md5 (first 10) | notes |
|---|---:|---|---|
| `V32W` | - | - | no multi-L file; single-L fitter only (`qc_fitter_norm_refine_v2.cpp`, 250 lines, md5 `7f677876c1`) |
| `V33O` | 1373 | `5bbee1f05c` | identical to V33N and V33P |
| `V33N` | 1373 | `5bbee1f05c` | identical |
| `V33P` | 1373 | `5bbee1f05c` | identical |
| `MF` | 2021 | `4b24bbe5e2` | adds the hot-window FCN |
| `L20` | 2107 | `cfbc500f4e` | `MF` + chi-square trace audit |
| `full_cache_generator_package_20260629-main/legacy_source` | 1856 | `64480ac491` | side branch (v32x with `scan_backend=cpu_kbasis|gpu_batched_det`, includes `v32x_fast_backends.hpp`); 1737 differing lines vs V33N; not on the main line |

The companion library `qc_fitter_norm_refine_v2.cpp` (included by the multi-L file with
`V32W_NO_MAIN`, `L20/source/qc_fitter_norm_refine_v2_multiL.cpp:15-17`) is 452 lines and
identical between `MF` and `L20` (md5 `674e2cb5f4`); `V33N/O/P` share a different 452-line copy
(md5 `ce7e3dacdd`); `V32W` has the original 250-line single-L version.

The shared header `K3df_minuit_fit_v32f_fullF3inv_QCfull_cached_classifier.hpp` is identical in
`V33O/V33N/V33P/L20` (md5 `e74d0696da`, 1595 lines); the `V32W` copy (md5 `9d6f534b68`, 1517 lines)
differs mainly by a Minuit2-availability shim (`V32F_HAS_MINUIT2` with a stub `FCNBase`).

**Answer: the newest and most complete fitter is
`L20/source/qc_fitter_norm_refine_v2_multiL.cpp` (2107 lines).**
Ordering is unambiguous from the diffs: `V33O = V33N = V33P` (1373) $\subset$ `MF` (2021) $\subset$ `L20` (2107).

### 3.2 What changed V33P (1373) $\to$ MF (2021): the hot-window FCN (+648 / -~60 lines)

`diff V33P MF` hunks (MF line numbers):

- New config keys `max_fcn_evals`, `root_search_mode`, `det_backend`, `projected_basis_cache_mode`,
  `projected_basis_cache_root`, `fallback_full_scan`, `window_half_rows`, `max_window_half_rows`,
  per-sector `accepted_zeros_file_L*_*`, `candidate_brackets_file_L*_*`, `det_grid_file_L*_*`,
  and the parameter mask `float_K3iso0/1/B/E` (MF 149-160, 224-250; L20 149-163, 228-256).
- GPU-cache reader rewritten to use `gpu_cache_reader_convention::require_gpu_cache_format_header`
  and the shared `read_gpu_cache_matrix_col_major`, **plus a second real/imag swap** to agree with
  the trusted `v33h` determinant scan (L20 664-666, 693-705; see Section 6.3).
- Basis precompute no longer done for all 20,000 rows at load time (L20 737-740).
- New structs `RuntimeBlock::BasisRow` (projected $F$ and four projected $\mathcal K$ basis matrices
  per row), `AcceptedZeroWindow`, `WindowGridRow` (L20 792-839).
- New loaders: `load_window_grid`, `load_projected_basis_hotwindows` (binary `V33PHWB1` format with
  schema/hash checks), `load_window_brackets`, `load_accepted_windows` (L20 873-1072).
- New root finder `evaluate_accepted_window` (L20 1074-1114) and the accepted-window branch of
  `MultiLFCN::model_for` (L20 1395-1408).
- Masked Minuit parameters, FCN-call guard, `sensitivity` run mode, `audit-l20-hot-window` run
  mode (`run_l20_hot_window_audit`, L20 1667-1797), timing counters.

### 3.3 What changed MF (2021) $\to$ L20 (2107): chi-square audit trace only (88 diff lines)

- Config keys `chi2_trace_path`, `chi2_level_trace_path`, `chi2_trace_test_name` (L20 154-156, 232-234).
- `AcceptedZeroWindow` gains `source_label_file`, `source_label_row` (L20 828-829, 1040-1049, 1068).
- Constructor opens two CSV traces and precomputes the audit weight matrix (L20 1318-1336).
- `operator()` calls `write_chi2_audit` every FCN (L20 1365); `write_chi2_audit` (L20 1516-1564)
  writes per-FCN and per-level chi-square contributions.
- No change to physics, classifier, root finding or chi-square itself.

### 3.4 Other L20-only refinements outside the fitter

- `v33p_build_projected_basis_hotwindow_cache.cpp`: L20 (477 lines) adds a one-row "halo"
  `kFitterRowHalo = 1` around each max window (L20 22-24, 190-193) so that the fitter's
  independent nearest-row rounding never requests an uncached row. MF copy (471 lines) lacks it.
- `v33h_patched_gpu_cache_oldscale_det_scan.cpp`: L20 adds a bridge API
  `assemble_scaled_projected_QC_for_row` (L20 510-532) used by the hot-window builder; V33P copy lacks it.
- `scripts/run_l20_full_cutoff_fits_v2.py`: L20 (430 lines) vs MF (422) adds `WINDOW_ROOT`/`CACHE_ROOT`
  module variables and copies the minimization-phase trace before the after-fit `fcn-once`
  truncates it (L20 lines 25-26, 320-324).

### 3.5 Earlier evolution (from docs)

| step | key change | evidence |
|---|---|---|
| v31zo | `digonto_classifier_v1` (sign flip + Gaussian-like peaks in $\min|\lambda|$ and min singular value), plotting-side only | `V32W/README_v31zo_110A2_DIGONTO_CLASSIFIER_V1.md:16-27` |
| v32s/v32t | `digonto_classifier_v2`: 10 new CPU points inside each sign-flip bracket, recursive split, local rescale to $[-1,1]$; v32t accepts one decreasing shoulder | `V32W/README_v32s_TRUE_REFINED_CLASSIFIER.md`, `README_v32t_oneShoulder.md` |
| v32w | trusted single-L fitter; fixed coarse grid kept separate from refinement cache; 4-parameter Minuit fit | `V32W/README_v32w_FIXED_COARSE.md:1-10` |
| v32x | multi-L (L=20,24), `.jack` input from `Ecm_data/data`, block-diagonal covariance, GPU coarse cache per (L, irrep) | `L20/source/qc_fitter_norm_refine_v2_multiL.cpp:1-11` |
| v33a/v33f | `digonto_v3_window` coarse-grid classifier (no refinement), `v4` merge | `qc_fitter_norm_refine_v2.cpp:31-50, 146-278` |
| v33g | runtime cache `*_runtime_K3basis.bin`: all 20,000 rows stored as projected $F$ + four projected $\mathcal K$ basis matrices | `L20/source/v33g_runtime_k3basis_cache.hpp:196-285` |
| v33m/v33n | classifier speed search, `raw_sign_only` chosen; cutoff scan 0.315-0.365 | `V33N/reports/v33m_production_classifier_benchmark_and_fit_full_report.md`, `v33n_production_sanity_cutoff_scan_full_report.md` |
| v33n fix | combined raw GPU cache reader convention fixed | `V33N/GPU_CACHE_READER_FIX.md` |
| v33p | human-label scoring, cache-anchored windows, all-sector validation plan | `V33P/README.md:9-24` |
| MF | accepted-window hot FCN, one-FCN checkpoint | `MF/HOT_WINDOW_FCN_CHECKPOINT.md` |
| L20 | projected-basis hot-window caches, pass-only cutoff fits, zero tests, chi-square traces | `L20/README.md` |

---

## 4. Q2 - Exactly how the fast FCN works (L20 code)

### 4.1 The quantization condition as implemented

At each grid energy $E_r$ and for one $(L,\Gamma)$ sector (box size $L/a_s$ = `Lbyas`, irrep $\Gamma$):

$$
Q(E_r;\vec c) \;=\; \frac{1}{(L\xi)^6}\, V_r^\dagger\!\left[F_3^{-1}(E_r) + \mathcal K_{\rm df,3}(E_r;\vec c)\right] V_r ,
\qquad D_r(\vec c) = \det Q(E_r;\vec c),
$$

with $\vec c = (K^{\rm iso}_0, K^{\rm iso}_1, K_B, K_E)$ = config names `K3iso0, K3iso1, K3B, K3E`.
Source of the formula: `L20/README.md:44-49`, `hpp:860-912` (`assemble_QC`), scale
`block_cnorm` = $(L\xi)^{6}$ (`multiL.cpp:271`, `const_norm_power=6`). The code states
"No determinant-level Nproj rescaling" (`L20/README.md:51`).

Because $\mathcal K_{\rm df,3}$ is linear in $\vec c$ (the code passes $\{K^{\rm iso}_0, K^{\rm iso}_1\}$ as a
coefficient vector and $K_B, K_E$ as scalars to `k3_2plus1::K3mat_2plus1`, `hpp:887-909`):

$$
\mathcal K_{\rm df,3}(E;\vec c) = \sum_{i=1}^{4} c_i\, \mathcal K^{(i)}(E), \qquad
\mathcal K^{(i)} \equiv \mathcal K_{\rm df,3}(E; \vec c = \hat e_i).
$$

Hence

$$
Q(E_r;\vec c) = \tilde F_r + \sum_{i=1}^{4} c_i\, \tilde B^{(i)}_r, \qquad
\tilde F_r = \frac{V_r^\dagger F_3^{-1}(E_r) V_r}{(L\xi)^6},\quad
\tilde B^{(i)}_r = \frac{V_r^\dagger \mathcal K^{(i)}(E_r) V_r}{(L\xi)^6}.
$$

**Yes - your guess is exactly what the code does.** Each basis matrix is projected once, offline.

- Builder (offline, expensive): `v33p_build_projected_basis_hotwindow_cache.cpp:245-258`
  computes `result.F = project(F3inv)` and `result.B0/B1/BB/BE = project(k3_matrix(entry,par,z0..zE))`
  with `project(M) = (Vsel^dagger M Vsel)/scale`, `scale = pow(L*xi,6)` (line 212). It also checks the
  basis sum against a direct full assembly at the reference $\vec c$ (lines 259-267):
  max matrix difference $\le 4.66\times10^{-10}$ (abs), determinant sign agreement 100 % in all five
  L20 sectors (`FULL_L20_FIT_AND_VALIDATION_REPORT.md:52-58`).
- An older equivalent existed in memory only: `precompute_projected_k3_basis` (`hpp:688-717`)
  fills `K3_proj_basis[0..3]` (unscaled), and `assemble_QC` fast path (`hpp:863-876`) sums them.
  The v33g runtime cache serialized the same five matrices for all 20,000 rows
  (`v33g_runtime_k3basis_cache.hpp:228-249`), unscaled (scale applied later in `eval_entry_QC`,
  `qc_fitter_norm_refine_v2.cpp:90`).

### 4.2 What is precomputed and cached (three cache levels)

| level | file | content per row | size/cost (from files) |
|---|---|---|---|
| raw GPU cache (external) | `/media/digonto/Data/F3inv_cache/v33e_coarse20000_sep_f2k2_L20_<irrep>/v32zu_Lbyas20_gpu_cache/v32zu_Lbyas20_xi3p444_irrep<irrep>_coarse20000_F3inv_Vsel_gpu.bin` | full $F_3^{-1}$ ($N\times N$), $V_{\rm sel}$ ($N\times N_{\rm proj}$) | e.g. 3,260,096,935 bytes for L20/200_A2 (`hot_windows.../L20_200_A2/projected_basis_hotwindows.json`) |
| v33g runtime cache | `cache/v33g_runtime_k3basis/v33g_Lbyas<L>_xi3p444_irrep<irrep>_coarse20000_runtime_K3basis.bin` (+`.meta.json`) | $V^\dagger F_3^{-1} V$ and four $V^\dagger\mathcal K^{(i)} V$, all rows, unscaled | load 6.23 s for 10 sectors (`V33N/reports/v33n_final_sanity_10fcn.md`) |
| hot-window projected-basis cache | `L20/output/hot_windows_L20_user_validated_0p36/L20_<irrep>/projected_basis_hotwindows.{bin,json}` | $\tilde F_r$, $\tilde B^{(1..4)}_r$ already divided by $(L\xi)^6$, only window rows | 1006-2511 rows per sector; 150,850,108 bytes for 200_A2; build 300-1572 s per sector at about 1.7 rows/s with 18 OpenMP threads (`MF/output/.../L20_USER_VALIDATED_0P36_ACCEPTED_WINDOW_FIX_REPORT.md`, JSON `openmp` field) |

Hot-window binary format (`multiL.cpp:975-986`): magic `V33PHWB1`, `uint32 version=1`, `uint32 count`,
then per row `uint64 row, double Ecm, int32 Nfull, int32 Nproj, uint32 dim`, then five column-major
`dim x dim` complex128 matrices (`Fproj, B0, B1, BB, BE`, read with `comp(re,im)`, line 915).
The JSON sidecar must contain `schema_version v33p_projected_basis_hotwindow_v1`,
`complex_convention variant_04_real_imag_swapped`, `scaling_convention divided_by_pow_Lxi_6`,
`K3_basis_terms` (line 956), the SHA-256 of the accepted-window CSV (checked against the file on
disk, lines 958-971) and `row_indices` that must match the binary exactly (lines 988-992).

### 4.3 Startup (outside the FCN)

`main` (`multiL.cpp:1810-2106`) $\to$ `multiconfig_from_config` (192-260) $\to$
`load_multil_jack_targets` (329-448) $\to$ `MultiLFCN` constructor (1263-1337):

1. For each $(L,\Gamma)$ block, `build_runtime_block` (1218-1259). With
   `projected_basis_cache_mode=hot_windows` (1225-1235): read the determinant-grid CSV only for its
   row/Ecm/Nfull/Nproj metadata (`load_basis_metadata_grid`, 938-947) and load the hot-window binary
   (949-994). **The 3 GB raw caches are never opened.**
2. With `root_search_mode=accepted_windows` (1295-1317): `fallback_full_scan` must be 0 (else
   throw, 1296); `det_backend=gpu_batched` is rejected (1297); `load_accepted_windows` builds the
   window list per block (1021-1072).
3. Timings: cache load 0.64-0.95 s per cutoff run (`load_time_sec` in
   `L20/output/l20_full_cutoff_fits_user_validated_0p36_passonly/Ecm_cut_0p3*/fit_metrics.csv`).
   For comparison, the MF raw-cache hot-window path needed 350.26 s and 7.65 GB RSS to load two
   20,000-row sectors (`MF/output/v33p_minimal_fitter_000_A1m/COLD_START_LOAD_TIME_NOTE.md`).

### 4.4 Accepted windows (how each lattice level gets its search window)

`load_accepted_windows` (`multiL.cpp:1021-1072`):

- Reads three per-sector files from the config: accepted zeros CSV (`accepted_zeros_file_L20_<irrep>`),
  candidate brackets CSV (`bracket_id,row_left,row_right`), and the determinant grid CSV.
- Keeps only rows with `user_label` true (line 1047) and inside the cutoff (1048-1050).
- **Hard gate**: number of accepted zeros must equal the number of lattice levels in that block (1052).
- Sorts zeros by energy (1053) and pairs the $k$-th zero with the $k$-th lattice level (1063).
- Window for each zero, with $c$ = grid row nearest the zero estimate, $[\ell,r]$ = bracket rows,
  $h$ = `window_half_rows` (default 50) and $H$ = `max_window_half_rows` (default 250) (1059-1066):

$$
\text{initial rows} = [\min(c-h,\ell),\ \max(c+h,r)],\qquad
\text{max rows} = [c-H,\ c+H].
$$

  With grid spacing $\Delta E = (0.36-0.26310)/19999 = 4.845\times10^{-6}$ (lattice units), the
  initial window is about $\pm 2.4\times10^{-4}$ and the maximal window about $\pm 1.2\times10^{-3}$
  around the human-validated zero.
- `previous_model_zero` is set to the human zero estimate (1067) and is **never updated** during
  the fit (the windows are `const` inside the FCN).

### 4.5 Inside one FCN call

`MultiLFCN::operator()` (`multiL.cpp:1339-1390`):

1. Expand masked parameters to 4 values (fixed ones from the config guesses) (1340-1349).
2. Hard FCN-call guard `max_fcn_evals` (1351-1353; throws).
3. `model_for` (1392-1428), accepted-window branch (1395-1408): for every block and every window
   call `evaluate_accepted_window` (1074-1114):
   - for each row $r$ in the window: $A_r = \tilde F_r + c_1\tilde B^{(1)}_r + c_2\tilde B^{(2)}_r + c_3\tilde B^{(3)}_r + c_4\tilde B^{(4)}_r$ (line 1090);
     `det_info` = Eigen `PartialPivLU` determinant (`qc_fitter_norm_refine_v2.cpp:85`); the
     classifier observable is $y_r=\operatorname{Re}D_r$ (1092).
   - run `find_QC_zeros_v3_from_coarse_grid` on the window rows (1099), optionally the v4 merge (1100);
   - keep only `true_zero` candidates; take the one closest to `previous_model_zero` (1101-1106);
   - if none, expand once to the max window (1108-1111) and retry; if still none, return empty.
   - root position inside a flip: linear interpolation between the two rows,
     $E_0 = E_i - y_i\,(E_{i+1}-E_i)/(y_{i+1}-y_i)$ (`linzero`, `qc_fitter_norm_refine_v2.cpp:92`).
4. `assign_block_model` (1116-1129): sort the block's zeros, assign $k$-th zero to the $k$-th lattice
   level; missing ones get model energy 0.
5. `found` = number of positive finite model energies; chi-square from `chi_square_v32f`; if
   `found < targets` the FCN returns `failure_penalty` = $10^{100}$ (1360-1362).
6. Exceptions inside the FCN are caught and converted to the penalty (1386-1389).

There is **no full-grid scan and no fallback** in this branch. Note that the logged
`fallback_full_scan=0 full_scan=0` is a hard-coded string (1371) and the counters
`full_scan_occurred`/`fallback_occurred` (775-776) are never set anywhere, so "fallback/full_scan
0/0" in the reports is true by construction, not a measured outcome.

The legacy branch (`root_search_mode=full_scan`, 1409-1425) is `find_QC_zeros_refined`
(`qc_fitter_norm_refine_v2.cpp:280-382`): determinant over all 20,000 rows (OpenMP), then
v3/v4 window classifier, or the v2 recursive refinement that builds new $F_3^{-1}$ rows on the CPU.

### 4.6 Chi-square and covariance

`chi_square_v32f` (`hpp:474-513`), with residual $r_i = E^{\rm lat}_i - E^{\rm model}_i$:

$$
\chi^2_{\rm raw\_cov\_inv} = r^{T} C^{+} r,\qquad
\chi^2_{\rm corr\_inv} = \tilde r^{T} R^{+}\tilde r,\ \tilde r_i = r_i/\sigma_i,\qquad
\chi^2_{\rm diagonal} = \sum_i (r_i/\sigma_i)^2 .
$$

Production uses `chi_square_mode = raw_cov_inv`. $C^{+}$ is `symmetric_pseudoinverse_v32f`
(`hpp:381-394`): symmetrize, `SelfAdjointEigenSolver`, invert eigenvalues with
$|\lambda| > \max(10^{-12}\max|\lambda|, 10^{-300})$. It is recomputed on every FCN call (cheap for
$\le 40$ levels, but could be done once). Any non-positive or non-finite model energy returns the
penalty (`hpp:166`).

Covariance construction in the multi-L fitter (`multiL.cpp:411-434`): files
`<L>_<irrep>_n<state>.jack` in `lattice_jackknife_dir`, one header line skipped, column 2 read
(`jack_energy_column=2`), each sample converted $E_{\rm cm}=\sqrt{E_{\rm lab}^2-P^2}$,
$P = 2\pi|\vec n|/(\xi L)$ (312-327), then

$$
C_{ij} = \frac{N-1}{N}\sum_{a=1}^{N}\bigl(x^{(i)}_a-\bar x^{(i)}\bigr)\bigl(x^{(j)}_a-\bar x^{(j)}\bigr),
$$

within one $L$; zero across different $L$ (block diagonal). **See Section 6.1 for a likely problem
with this formula.**

$\nu = n_{\rm levels} - n_{\rm floating}$ (`multiL.cpp:1604`).

### 4.7 Minuit2 setup and parameter mask

`multiL.cpp:2062-2097`: `MnUserParameters` gets only floating parameters (2064-2067) with steps from
config; optional limits (2068); `MnMigrad migrad(fcn,u)`, `migrad(max_fcn)` if `max_fcn>0` (2070-2072);
`MnHesse` only if the minimum is valid (2073); best values re-evaluated with `model_for` and
chi-square recomputed (2086-2090); `Up()=1` (1338). The cutoff driver sets
`max_fcn=400`, `max_fcn_evals=500`, steps $(10^3, 10^4, 10^3, 10^4)$ for
$(K^{\rm iso}_0,K^{\rm iso}_1,K_B,K_E)$ and warm-starts each cutoff from the previous best fit
(`scripts/run_l20_full_cutoff_fits_v2.py:178-183, 327`).

Latent bugs (my reading, not reported in the docs):

- `SetLimits(i,...)` loops `i=0..3` regardless of the mask (2068); with fewer than 4 floating
  parameters and `use_parameter_limits=1` this addresses non-existent parameters.
- `minuit_covariance_to_eigen_v32f(min,4)` (2091) fills the first $n_{\rm float}$ rows in Minuit
  order but labels them `K3iso0 K3iso1 K3B K3E` (1624-1625); correct only when the fixed parameters
  are the last ones (true for the zero tests, which fix K3B, K3E).

### 4.8 Run modes of the binary `bin/v33f_k3df_fitter_multiL_v33e <config> <mode>`

`fit` (default), `fcn-once`/`fcn`, `spectrum-only`, `benchmark-fcn --repeat --warmup`, `sensitivity`,
`dump-targets`, `pipeline-test`, `audit-l20-hot-window` (`multiL.cpp:1834-2060`). Build:
`L20/build.sh` (g++ -std=c++17 -O3 -fopenmp, Eigen, Minuit2, `Faddeeva.cc`).

### 4.9 Timings per FCN (all from files)

| version / mode | targets | time per FCN | source |
|---|---:|---:|---|
| v33m/v33n, v33g runtime cache, full 20,000-row scan, 10 sectors, `raw_sign_only` | 32 | $1.37$-$1.43$ s (10-call avg; std not given; range 1.32-1.51 s over modes) | `V33N/reports/v33m_production_classifier_benchmark_and_fit_full_report.md:31-39` |
| MF hot window, raw cache in memory, 000_A1m L20+L24 | 4 | $7.36\times10^{-4}$ s (404 rows) | `MF/HOT_WINDOW_FCN_CHECKPOINT.md:8-13` |
| MF reduced fit, same | 4 | avg $5.44\times10^{-4}$ s, max $7.39\times10^{-4}$ s | `MF/output/v33p_minimal_fitter_000_A1m/FITTER_000_A1M_REDUCED_MINIMIZATION_REPORT.md` |
| L20 projected-basis hot window, cutoff 0.315 | 6 | $3.22\times10^{-3}$ s | `L20_PASSONLY_FIT_SUMMARY.csv` |
| same, 0.325 / 0.335 / 0.345 | 9 / 12 / 14 | $6.50\times10^{-3}$ / $1.86\times10^{-2}$ / $2.38\times10^{-2}$ s | same |
| full L20 fit wall time per cutoff | 6-14 | 1.12 / 2.31 / 2.79 / 3.75 s (88/175/86/108 FCN calls) | same |

The per-FCN loop over windows is serial (no OpenMP pragma in `evaluate_accepted_window` or the
accepted-window part of `model_for`), and the v3 classifier prints one line per window to stdout
(`qc_fitter_norm_refine_v2.cpp:224`); both are easy speed-ups.

---

## 5. Q3 - The classifier and the human-label workflow

### 5.1 Classifier generations

| name | where | rule (paraphrased from code/docs) |
|---|---|---|
| `digonto_classifier_v1` | v31zo plotting scripts | accept a $\det$ sign flip only if $\min|\lambda|$ and the minimum singular value both show a Gaussian-like local peak (`V32W/README_v31zo_...md:16-27`) |
| `digonto_classifier_v2` (`v2_refine`) | `qc_fitter_norm_refine_v2.cpp:102-144, 280-382` | for each coarse flip, build `v2_refine_points=10` new CPU $F_3^{-1}$ points inside, split recursively if several flips (depth $\le 8$), then `classify`: true zero if a shoulder value drops to $\le$ `zero_ratio`$=0.8$ of the bracket-edge value |
| `digonto_v3_window` | `qc_fitter_norm_refine_v2.cpp:146-226`, doc `V33P/docs/DIGONTO_CLASSIFIER_V3_ALGORITHM.md` | coarse grid only; for flip between rows $i,i+1$ use shoulders $(i-2,i-1,i)$ and $(i+3,i+2,i+1)$; rescale by local max; a shoulder is zero-like if every step toward the flip drops $|y|$ by $\ge 15\,\%$, pole-like if it rises by $\ge 15\,\%$; 3-to-2 trimming if the middle point is a local extremum (tolerance 2 %); true zero if at least one shoulder zero-like and not both pole-like |
| `digonto_v4_window` | `qc_fitter_norm_refine_v2.cpp:228-278` | v3 then merge true zeros closer than `v4_merge_tol`$=10^{-3}$ (keep the one with the smallest scaled $|y|$ next to the flip; energy = running mean) |
| `raw_sign_only`, `raw_sign_clustered`, `algo_v6_*`, `algo_v7_*` | dispatch table `qc_fitter_norm_refine_v2.cpp:61-74` | v33m speed-search variants; in this code `raw_sign_only` and `algo_v7_*` dispatch to the v3 path, `raw_sign_clustered`/`algo_v6` to v4 (lines 41-50) |

The determinant observable for all of them is $\operatorname{Re}\det Q$ on the grid. Sign changes
across a change of $N_{\rm proj}$ (a "dimension jump") must not be classified as zeros
(`V33N/CLASSIFIER_RESTART_NOTES.md:32-34`).

### 5.2 Why human labels were needed

- On the corrected L20/100_A2 grid in $[0.34,0.36]$, 8 sign-change brackets were found
  (`V33N/CLASSIFIER_RESTART_NOTES.md:15-24`); on the full window $[0.2631,0.36]$ 25 brackets,
  of which the human marked only 2 as true zeros (`L20/input/v33p_user_truezero_labels_L20_100_A2_E026301_0360.csv`).
  Overall the L20 label files contain 19 true vs 117 false brackets (3/7, 2/23, 4/32, 5/25, 5/30 for
  000_A1m/100_A2/110_A2/111_A2/200_A2, `MF/output/v33p_projected_basis_cache/l20_full_cutoff_fits_user_validated_0p36/L20_USER_VALIDATED_0P36_ACCEPTED_WINDOW_FIX_REPORT.md`).
  So most sign changes are poles or spurious, and shape heuristics alone were not trusted.
- The v33m speed search showed different classifiers produce different root assignments with the same
  coverage (`v4`: $\chi^2=0.0848$ vs $0.2544$ for `raw_sign_only`, both 32/32), i.e. coverage is not
  correctness (`V33N/reports/v33m_...md:17-29`).
- Reader/convention bugs had produced spurious imaginary parts and broken sign logic
  (`V33N/GPU_CACHE_READER_FIX.md:5-12`), so earlier "accepted" zeros were not reliable.
- The stated rule became: "Raw-sign modes may generate candidate brackets for diagnostics but are not
  production classifiers" (`V33P/README.md:125-126`) and "Accept a sector only if user-labeled `FP=0` and
  `FN=0`" (`V33P/docs/NEXT_VALIDATION_PLAN_ALL_L_IRREPS_CACHE_ANCHORED.md:10`).

### 5.3 Workflow (from `L20/README.md:64-75` and scripts)

1. CPU/OpenMP determinant scan of the raw cache (`v33h_patched_gpu_cache_oldscale_det_scan.cpp`)
   or GPU LU scan (`v33p_cpuassembled_gpu_det_scan.cpp`, `v33p_gpu_det_lu.cu`) $\to$ `det_grid_*.csv`.
2. Candidate sign-change brackets $\to$ `candidate_sign_change_brackets*.csv`.
3. Interactive determinant-family plot (`scripts/plot_v33p_classifier_det_family_interactive.py`,
   curves $\operatorname{Re}D(E)\cdot 10^k$) and v3/v4 predictions.
4. Human fills `user_label` (true/false) per bracket: CSV columns
   `bracket_id,E_left,E_right,user_label,user_zero_estimate,comment`.
5. Scoring (`V33P/scripts/score_v33p_classifier_against_user_labels.py:44-61`): a predicted true zero
   is TP if the nearest labeled bracket midpoint is true and within `tolerance` $=8\times10^{-6}$; else FP;
   FN = true labels not matched; TN = false labels minus FP; score $=TP - FP - 0.25\,FN$; acceptable iff
   $FP=FN=0$. Result for the reference sector L20/100_A2: v3 and v4 both $TP=2, FP=0, FN=0, TN=23$,
   accepted $E=0.291556099879,\ 0.300441497547$ (`V33P/README.md:16-21`).
6. `scripts/build_l20_user_validated_windows.py` writes `accepted_windows.csv` from true rows only;
   `materialize_v33p_hot_windows.py` is the older equivalent with explicit window arithmetic (lines 130-145).
7. `v33p_build_projected_basis_hotwindow_cache` builds the hot-window cache from those windows;
   `merge_projected_basis_rows.py` merges halo rows.
8. Prefit audit: per cutoff and sector, lattice level count = user-true window count, and all requested
   rows present in the cache (`run_l20_full_cutoff_fits_v2.py:249-282`).

Important: in the L20 production path, the human labels decide which roots exist; the v3 classifier
inside the FCN only re-finds the root in a 101-row window as $\vec c$ changes.

---

## 6. Problems found during this review (read before trusting any L20/v33n chi-square)

### 6.1 The multi-L covariance is probably inflated by a factor $(N-1)^2 \approx 3.6\times10^{5}$

Observation (computed by me, command in Section 9):

| file | $N$ | mean | sample s.d. | $\sqrt{(N-1)\,\mathrm{var}}$ | s.d.$/\sqrt N$ |
|---|---:|---:|---:|---:|---:|
| `L20/input/lattice_spectrum/20_000_A1m_n0.jack` | 603 | 0.269227 | 0.012643 | 0.310209 | 0.000515 |
| `L20/input/lattice_spectrum/20_001_A2_n0.jack` | 603 | 0.305566 | 0.011828 | 0.290199 | 0.000482 |

The fitter's reported lattice error for 000_A1m level 0 is `0.31020860761372221`
(`Ecm_cut_0p335/l20_full_cut_0p335_fit_levels_allL.dat`, row 0), i.e. $\sigma/E \approx 115\,\%$.
The covariance eigenvalues in the MF audit are $0.0137$ to $0.264$
(`MF/output/v33p_minimal_fitter_000_A1m/COVARIANCE_INVERSE_AUDIT.md`).

Cause: the legacy v32w reader first forms leave-one-out averages and then applies the
jackknife formula (`L20/source/lattice_data_covariance_cpp.hpp:126-163, 293-298, 373-380`):

$$
\tilde x_a = \frac{1}{N-1}\sum_{b\ne a} x_b \;\Rightarrow\; \tilde x_a-\bar{\tilde x} = -\frac{x_a-\bar x}{N-1},
\qquad C^{\rm v32w}_{ij} = \frac{N-1}{N}\sum_a(\tilde x^{(i)}_a-\bar{\tilde x}^{(i)})(\tilde x^{(j)}_a-\bar{\tilde x}^{(j)}) = \frac{C^{\rm multiL}_{ij}}{(N-1)^2}.
$$

The multi-L reader (`multiL.cpp:290-311, 414-428`) applies the same $\frac{N-1}{N}\sum$ formula
directly to the file values, skipping the leave-one-out step. If the `.jack` files contain
per-configuration values (which their spread and the 115 % error strongly suggest; the header line is
`603 1 0 0 1`), then $C^{\rm multiL} = (N-1)^2\, C^{\rm correct}$ and

$$
\chi^2_{\rm multiL} = \frac{\chi^2_{\rm correct}}{(N-1)^2} = \frac{\chi^2_{\rm correct}}{362{,}404}\quad (N=603).
$$

Rescaling the reported L20 pass-only values by $362{,}404$ gives $\chi^2 = 11.8,\,34.5,\,46.5,\,60.6$ and
$\chi^2/\nu = 5.9,\,6.9,\,5.8,\,6.1$ at cutoffs 0.315/0.325/0.335/0.345, versus the trusted v32w
$\chi^2/\nu=1.44$ at cutoff 0.335 (12 levels). These rescaled numbers are only indicative: the
minimizer stopped early (next point), so they are not minima of the corrected objective.

Consequences, all visible in the outputs:

- Minuit2 stops when EDM $< 0.002\times\text{tol}\times\text{Up}$ (default tol $0.1$, so $2\times10^{-4}$;
  this threshold is from my knowledge of Minuit2, please verify). With a total $\chi^2\sim10^{-4}$ that
  test passes almost immediately. Indeed $K_E$ never moved from its start value
  $-1226756.7068845264$ in any L20 fit, $K_B$ did not move at 0.315-0.335, $K^{\rm iso}_1$ did not
  move at 0.315/0.325 (`L20_PASSONLY_FIT_SUMMARY.csv`), and in v33n only $K^{\rm iso}_0$ moved at
  0.315-0.335 (`V33N/reports/v33n_cutoff_scan/cutoff_scan_summary.md`).
- Hesse errors are defined by $\Delta\chi^2=\mathrm{Up}=1$, so they are inflated by $(N-1)=602$:
  e.g. $\sigma(K^{\rm iso}_0)=1.13\times10^{8}$ at 0.335 would become about $1.9\times10^{5}$.
  The value "2.00" for all four errors at 0.325 and in the zero tests is a degenerate Hesse output.
- The docs explicitly attribute the tiny values to "covariance scales around 0.3--0.4" and call them
  "not evidence of double division" (`FULL_L20_FIT_AND_VALIDATION_REPORT.md:107-112`); that check
  verified arithmetic consistency, not the input-format assumption.
- The `AGENTS.md` "starting point" `chi2=0.00027552010191639938 model_found=32/32`
  (`V33P/AGENTS.md:78-86`) shows the issue entered with v32x (the first multi-L reader).

What I am NOT sure of: the true content of the `.jack` files. Digonto should confirm whether column 2
holds per-configuration values or jackknife samples. A quick check: compare
$\sigma$ for one level between the v32w `masses` file
(`/home/digonto/Codes/Practical_Lattice_v2/3body_quantization/lattice_data/KKpi_interacting_spectrum/twoptvar_analysis/masses/<ensemble>_<mom>_state_<n>`) and the `.jack` file.

### 6.2 Fit-quality caveats independent of 6.1

- Level matching is "k-th accepted zero to k-th lattice level" per $(L,\Gamma)$ (`multiL.cpp:1053, 1124-1127`;
  rule restated in `V33P/AGENTS.md:88-97`). It cannot represent a model level that has no lattice
  counterpart; this is exactly the 0.355 blocker (Section 7.3).
- If the fit moves a root more than about $1.2\times10^{-3}$ (max window) or the v3 test rejects it,
  the FCN jumps to $10^{100}$; this caused the zero-test root loss at 0.335. A smooth penalty or window
  re-centring on the last good model zero would help.
- Root resolution relies on linear interpolation over $\Delta E = 4.8\times10^{-6}$; fine if the true
  $\sigma\sim5\times10^{-4}$, but should be checked once the covariance is fixed.
- `parameter_impact_summary.csv` at 0.335 reports `max_abs_dE_dparam = 0.3068451968442849` for all four
  parameters, which equals the model energy of level 3 (0.30684519684428491); with steps $10^3$-$10^4$
  a derivative of that size is not plausible. The plus-step of $K^{\rm iso}_1$ lost a root
  ($\Delta\chi^2_+ = 10^{100}$). Treat that file as unreliable.

### 6.3 Complex read convention: two swaps that cancel

`gpu_cache_reader_convention::read_complex` returns `comp(b,a)` for stored doubles `(a,b)`
(`L20/source/gpu_cache_reader_convention.hpp:179-186`); the fitter then swaps again
(`multiL.cpp:693-705`), and so does v33h (`v33h_patched_gpu_cache_oldscale_det_scan.cpp:383-392`).
Net effect: stored doubles are read as $(\mathrm{Re},\mathrm{Im})$, which is exactly what the V33N/O/P
fitter did with a plain loop (`V33N/source/qc_fitter_norm_refine_v2_multiL.cpp:653-660`).
However `L20/source/v33g_build_runtime_k3basis_cache.cpp:216-217` uses the single-swap reader without
the second swap, whereas the V33N copy (`V33N/source/v33g_build_runtime_k3basis_cache.cpp:215-222`)
reads plainly. **Any `*_runtime_K3basis.bin` built with the L20-era builder would be inconsistent.**
For the consolidated repo: read plainly, delete the swap helpers, keep the hermiticity and
$\max|\operatorname{Im}\det|$ checks as a guard.

### 6.4 Portability of the release tree

- The hot-window JSON stores an absolute `accepted_windows_path` under `/home/digonto/...MF...`
  and the fitter verifies its SHA-256 (`multiL.cpp:964-971`), so the release fails unless that path
  exists or is rewritten.
- `run_l20_user_validated_passonly.py:13-15` and `build_l20_user_validated_windows.py:11-12` point at
  `output/v33p_gpu_classifier_all_irrep/...` and `output/v33p_projected_basis_cache/...` (the MF
  layout), not the L20 release layout `output/accepted_windows_L20_user_validated_0p36`,
  `output/hot_windows_L20_user_validated_0p36`.
- `fitter_input_base_user_validated.in` contains absolute MF paths and relative paths to
  `output/v33p_fasttrack_cpu_existing_000_A1m/...` grids not shipped in L20.
- The five `.bin` hot-window payloads are not in Git (three exceed 100 MB)
  (`L20/LARGE_FILES_AND_CACHE_POLICY.md:1-16`).
- `K3df_minuit_fit_v32f...hpp:60-62` pulls in a whole debug program
  (`debug_v30r_svd_validated_cartesian_l1_100_A2_zero_scan.cpp`, 2203 lines) with `main` renamed;
  the include closure of the fitter is about 20 files (Section 8).

---

## 7. Q4 - Fit results achieved

### 7.1 v32w (trusted single-L baseline, 15 Jun)

L=20, cutoff 0.335, 12 levels, 4 parameters, `raw_cov_inv` with the v32w masses-file covariance,
`v2_refine` classifier, coarse grid 10,000 rows, start values all $0.1$:

$$
\chi^2 = 11.489748255611568,\quad \nu=8,\quad \chi^2/\nu = 1.436218531951446,\quad 12/12 \text{ found},
$$
$K^{\rm iso}_0=34470.228772793678$, $K^{\rm iso}_1=-1035173.9140440958$, $K_B=253921.14757953485$,
$K_E=-961680.91711685632$ (`KKpi-v32w-v32y-repro-master/.../README.md:60-74`,
`docs/VERSION_SUMMARY.md:13-25`). No parameter errors are quoted there.

### 7.2 v33n cutoff scan (L20+L24, 9 Jul; computed with `raw_sign_only`, v33g runtime cache, full scan)

From `V33N/reports/v33n_cutoff_scan/cutoff_scan_summary.md` (all CONVERGED, start = v32w values):

| cutoff | targets | final $\chi^2$ | $\chi^2/\nu$ | $K^{\rm iso}_0$ | $K^{\rm iso}_1$ | $K_B$ | $K_E$ | FCN calls | fit wall (s) |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 0.315 | 18 | 0.05942 | 0.00424 | 60236.20 | -1035173.91 | 253921.15 | -961680.92 | 179 | 247.27 |
| 0.325 | 24 | 0.06908 | 0.00345 | 60136.15 | -1035173.91 | 253921.15 | -961680.92 | 179 | 261.58 |
| 0.335 | 32 | 0.25175 | 0.00899 | 60243.85 | -1035173.91 | 253921.15 | -961680.92 | 179 | 259.98 |
| 0.345 | 36 | 0.61945 | 0.01936 | 29067.16 | -1007464.29 | 260656.20 | -917434.50 | 271 | 394.76 |
| 0.355 | 39 | 0.66314 | 0.01895 | 28317.16 | -1035173.91 | 253921.15 | -961643.05 | 194 | 287.05 |
| 0.365 | 40 | 1.60577 | 0.04460 | 28316.86 | -1035173.91 | 253921.15 | -961670.13 | 178 | 262.37 |

Wall times from `V33N/diagnostics/v33o/cutoff_timing_summary.csv` (reconstructed from logs).
Minuit warned "Initial matrix not pos.def." at 0.315-0.335. These roots come from `raw_sign_only`,
later declared non-production, so these fits should be considered superseded. Under the 6.1
hypothesis the 0.335 value would correspond to $\chi^2 \approx 9.1\times10^{4}$, i.e. a poor fit.

### 7.3 L20 pass-only projected-basis hot-window fits (15 Jul)

L=20 only, sectors 000_A1m, 100_A2, 110_A2, 111_A2, 200_A2; 19 user-true windows through
$E_{\rm cm}=0.36$; `digonto_v3_window`, `accepted_windows`, `hot_windows`; start values
$(73735.840894011912, -972421.14060757787, 347174.05548116949, -1226756.7068845264)$, warm-started
cutoff to cutoff. From `L20/output/l20_full_cutoff_fits_user_validated_0p36_passonly/L20_PASSONLY_FIT_SUMMARY.csv`:

| cutoff | levels | $\chi^2$ (as reported) | $\nu$ | $\chi^2/\nu$ | $K^{\rm iso}_0$ | $K^{\rm iso}_1$ | $K_B$ | $K_E$ | FCN | Hesse status |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---|
| 0.315 | 6 | $3.2621\times10^{-5}$ | 2 | $1.6311\times10^{-5}$ | 151422.46 | -972421.14 | 347174.06 | -1226756.71 | 88 | VALID, errors $1.3$-$88\times10^{8}$ |
| 0.325 | 9 | $9.5309\times10^{-5}$ | 5 | $1.9062\times10^{-5}$ | 125968.22 | -972421.14 | 347174.06 | -1226756.71 | 175 | all errors 2.00 (degenerate) |
| 0.335 | 12 | $1.2833\times10^{-4}$ | 8 | $1.6041\times10^{-5}$ | 32040.84 | -1074716.65 | 347174.06 | -1226756.71 | 86 | errors $0.8$-$6.1\times10^{8}$ |
| 0.345 | 14 | $1.6711\times10^{-4}$ | 10 | $1.6711\times10^{-5}$ | 40607.52 | -1067567.87 | 376187.14 | -1226756.71 | 108 | errors $0.3$-$11\times10^{8}$ |

All 4 found every root; "fallback/full scan 0/0" (see 4.5 caveat). An earlier non-user-validated run
(`MF/output/v33p_projected_basis_cache/l20_full_cutoff_fits_v2/L20_FULL_ECM_CUTOFF_FIT_SUMMARY.md`)
gave $K^{\rm iso}_0 = 229109, 203655, 69164$ at 0.315/0.325/0.335 and was blocked at 0.345/0.355.

### 7.4 K3B/K3E zero tests (`L20/output/l20_k3b_k3e_zero_tests/L20_K3B_K3E_ZERO_TEST_FINAL_REPORT.md:7-61`)

| test | cutoff | $\chi^2/\nu$ | found | $K^{\rm iso}_0$ | $K_B$, $K_E$ | status |
|---|---:|---:|---|---:|---|---|
| fixed $K_B=K_E=0$ | 0.315 | $1.0543\times10^{-5}$ | 6/6 | 87875.46 | 0, 0 | PASS |
| fixed | 0.325 | $1.2474\times10^{-5}$ | 9/9 | 102015.40 | 0, 0 | PASS |
| fixed | 0.335 | - | 11/12 | - | 0, 0 | BLOCKED (root lost) |
| float from 0 | 0.315 | $2.0395\times10^{-5}$ | 6/6 | 99189.40 | 0, 0 ($\sigma=2.0$) | PASS |
| float from 0 | 0.325 | $1.6929\times10^{-5}$ | 9/9 | 124643.52 | 0, 0 ($\sigma=2.0$) | PASS |
| float from 0 | 0.335 | - | 11/12 | - | - | BLOCKED |

$K^{\rm iso}_1$ stayed at $-972421.14$ in every zero test. Independent reconstruction of $\chi^2$
agreed with the fitter to about $10^{-14}$ relative (Table B), which validates the arithmetic only.
That $K_B, K_E$ "float from zero" never left zero is consistent with the early-stop explanation in 6.1.

### 7.5 Blockers

1. **0.355**: L20/200_A2 has 4 lattice levels below 0.355 but 5 user-true zeros
   (0.28624, 0.29918, 0.33098, 0.34798, 0.35346 from `L20/input/v33p_gpu_user_truezero_labels_L20_200_A2.csv`
   rows 5, 9, 18, 27, 32), so the count gate fails (`PREFIT_AUDIT_PASSONLY.md`, last table row).
2. **L24 inactive**: determinant scans for L24/110_A2 and 111_A2 paused; windows not validated
   (`L24_NEXT_STEPS.md:3-9`). The MF one-FCN test used L24/000_A1m (3 roots) successfully.
3. **L24/100_A2 bracket 36**: user label `false` at $E\in[0.340769, 0.340774]$
   (`L20/input/v33p_user_truezero_labels_L24_100_A2_E026310_0360.csv`, row 37), described as a
   classifier false positive to be resolved before that sector is used.
4. (Mine) covariance normalization (6.1) and the read-convention inconsistency (6.3).

---

## 8. Q5 - Stated next steps (verbatim)

`L20/FINAL_SHUTDOWN_HANDOFF_REPORT.md:96-103` ("Exact next action"):

> Resolve the L20/200_A2 one-to-one lattice/window mapping, then rerun only the
> 0.355 prefit audit and fit if the mapping is unambiguous. Separately finish
> L24 determinant scans, manually validate its windows, exclude L24/100_A2 until
> the bracket-36 false positive is resolved, and only then prepare combined L20+L24
> fits. Never promote false or blank user labels and never enable original-cache
> fallback inside the FCN.

`L20/README.md:124-128`:

> Resolve L20/200_A2's extra user-true window before attempting .355. Then
> resume L24 determinant validation, resolve L24/100_A2 bracket 36, build only
> user-validated L24 windows/caches, and run L24 prefit/one-FCN checks before a
> combined fit. Never silently fall back to original full-cache reads inside an
> FCN.

`L20/L24_NEXT_STEPS.md:5-15`:

> 1. Confirm/resume the paused L24/110_A2 and L24/111_A2 determinant scans.
> 2. Verify the cache metadata Ecm window and corrected raw-reader provenance.
> 3. Generate determinant-family plots and candidate brackets for each L24 irrep.
> 4. Manually label every candidate; only explicit `user_label=true` rows may be promoted.
> 5. Resolve L24/100_A2 bracket 36 before including that sector.
> 6. Build `accepted_windows_L24_user_validated_*` and expanded projected-basis hot-window caches.
> 7. Run row-coverage and metadata/hash audits.
> 8. Run one-FCN tests with `fallback_full_scan=0` before any minimization.
> 9. Run combined L20+L24 fits only after sector-level classifier and model-root counts pass.
>
> Do not run cachegen, use false/blank windows, or enable original-cache fallback.

`L20/FULL_L20_FIT_AND_VALIDATION_REPORT.md:109-112` (chi-square reporting):

> Report the new quantity as
> `fitter-reported raw covariance chi2/dof` until the legacy fit is recomputed
> with the same covariance and residual convention.

`MF/HOT_WINDOW_FCN_CHECKPOINT.md:33-35` (earlier checkpoint):

> Run only the bounded short `000_A1m` minimization. Do not start production
> minimization, all-sector validation, cachegen, or the GPU scanner from this
> checkpoint.

`V33P/README.md:135-141` (GPU determinant backend):

> The required next implementation is a chunked double-complex GPU determinant
> backend behind the corrected reader, with grouped dimensions, LU pivot parity,
> partial temporary CSVs, progress, CPU fallback, and CPU-reference parity.
> ... Until that gate passes, the CPU
> scanner remains the trusted reference and GPU results must not drive fitter
> validation.

`V32W/../docs/VERSION_SUMMARY.md` gives no next step; `V33N/VALIDATION_STATUS.md:55-56` lists
"The classifier still needs user-labeled true zeros before algorithm tuning can continue" and
"L25 non-interacting coverage is still partial".

My recommended additions (not in the docs): (a) settle the `.jack` covariance question (6.1) before
any further fit; (b) re-run the four L20 pass-only fits with a corrected covariance and a tighter
Migrad tolerance, starting from v32w values; (c) read the raw cache plainly and rebuild/verify any
runtime cache made by the L20-era v33g builder (6.3).

---

## 9. Q6 - Recommended minimal fitter file set for the consolidated repo

Canonical root for everything below unless stated: `L20 = /home/digimonk/Projects/three-body-qc-codes/extracted/KKpi-L20-user-validated-hotwindow-fitter-main/KKpi-L20-user-validated-hotwindow-fitter-main`.

### 9.1 C++ sources (fitter executable `v33f_k3df_fitter_multiL_v33e`)

Include closure computed by following `#include "..."` from the two entry files (commented-out
includes excluded):

| file (under `L20/source/`) | lines | role |
|---|---:|---|
| `qc_fitter_norm_refine_v2_multiL.cpp` | 2107 | **main fitter**: config, data, caches, windows, FCN, Minuit, run modes |
| `qc_fitter_norm_refine_v2.cpp` | 452 | classifier v2/v3/v4, `det_info`, `eval_entry_QC`, config parser (included as a library) |
| `K3df_minuit_fit_v32f_fullF3inv_QCfull_cached_classifier.hpp` | 1595 | cache entry struct, `assemble_QC`, basis precompute, `chi_square_v32f`, pseudoinverse |
| `v33g_runtime_k3basis_cache.hpp` | 287 | runtime cache I/O; also forwards the GPU complex reader |
| `gpu_cache_reader_convention.hpp` | 199 | raw GPU cache header/format checks (simplify per 6.3) |
| `lattice_data_covariance_cpp.hpp` | 511 | legacy jackknife reader (reference for 6.1) |
| `debug_v30r_svd_validated_cartesian_l1_100_A2_zero_scan.cpp` | 2203 | pulled in by the `.hpp` for physics helpers; should be split into a header |
| `F3_cpu_openmp_v25_K3QC_cached_core.hpp`, `fv_projector_cartesian_l1_v30m.hpp`, `K3_functions_2plus1.hpp`, `F2_functions_v2.h`, `G_functions_v2.h`, `K2_functions_v2.h`, `QC_functions_v2.h`, `functions.h`, `projections.hpp`, `projections_v1.hpp`, `projections_from_config.hpp`, `real_wigner_d.hpp`, `spherical_functions.h`, `eigenvalue_tracker.hpp`, `dig_tools.hpp`, `Faddeeva.cc`, `Faddeeva.hh` | - | physics core (shared with the cache-generator area; consolidate with that report) |

Offline cache builders and reference scanner:

| file | role |
|---|---|
| `L20/source/v33p_build_projected_basis_hotwindow_cache.cpp` | **hot-window projected-basis builder** (L20 copy has the row halo) |
| `L20/source/v33h_patched_gpu_cache_oldscale_det_scan.cpp` + `L20/source/compare_gpu_cache_vs_cpu_openmp_det_oldscale.cpp` | trusted CPU/OpenMP determinant scan of the raw cache; L20 copy has the bridge API needed by the builder |
| `V33N/source/v33g_build_runtime_k3basis_cache.cpp` (`/home/digimonk/Projects/three-body-qc-codes/extracted/kkpi-v33n-corrected-gpu-cache-reader-main/kkpi-v33n-corrected-gpu-cache-reader-main/source/v33g_build_runtime_k3basis_cache.cpp`) | full-grid runtime cache builder; prefer the V33N copy (plain read) over L20's single-swap copy, or port to plain reading |
| `L20/source/v33p_gpu_det_lu.cu`, `v33p_gpu_det_lu.hpp`, `v33p_cpuassembled_gpu_det_scan.cpp` | optional GPU determinant scan used to produce the 110/111/200_A2 grids (not validated as production per `V33P/README.md:128-141`) |

Build script: `L20/build.sh`.

### 9.2 Python scripts (under `L20/scripts/`)

| script | role | note |
|---|---|---|
| `run_l20_full_cutoff_fits_v2.py` | cutoff-suite driver: prefit audit, fcn-once, fit, after-fit, sensitivity, reports | fix paths (6.4) |
| `run_l20_user_validated_passonly.py` | wrapper restricting to 0.315-0.345 | fix paths |
| `build_l20_user_validated_windows.py` | user labels $\to$ accepted windows | fix paths |
| `materialize_v33p_hot_windows.py` | generic per-sector window materializer | simplest reference for window arithmetic |
| `merge_projected_basis_rows.py` | merges halo rows into hot-window caches | |
| `audit_l20_user_validated_prefit.py`, `audit_l20_chi2_reconstruction.py` | prefit and independent $\chi^2$ audits | |
| `plot_v33p_classifier_det_family_interactive.py` | determinant-family QA plot for labeling | |
| `/home/digimonk/Projects/three-body-qc-codes/extracted/KKpi-v33p-cache-anchored-classifier-fitter-main/KKpi-v33p-cache-anchored-classifier-fitter-main/scripts/score_v33p_classifier_against_user_labels.py` | TP/FP/FN/TN scoring | only in V33P |
| `.../KKpi-v33p-cache-anchored-classifier-fitter-main/.../scripts/run_v33p_classifier_sweep_L20_100_A2.py` | CSV-only classifier harness | only in V33P |
| `run_l20_k3b_k3e_zero_tests.py`, `finalize_l20_k3b_k3e_zero_tests.py` | zero tests | optional |

### 9.3 Configs and data

| item | path |
|---|---|
| base fitter config actually used | `L20/output/l20_full_cutoff_fits_user_validated_0p36_passonly/fitter_input_base_user_validated.in` (needs path rewrite) |
| per-cutoff resolved configs | `L20/output/l20_full_cutoff_fits_user_validated_0p36_passonly/Ecm_cut_0p3{15,25,35,45}/fitter_input.txt` |
| v33m/v33n full-scan config (reference) | `L20/configs/config_v33m_multiL_all_irreps_fastest_production.in`, `L20/configs/v33n_cutoff_scan/*.in` |
| v32w trusted config | `/home/digimonk/Projects/three-body-qc-codes/extracted/KKpi-v32w-v32y-repro-master/KKpi-v32w-v32y-repro-master/latest_working_versions/v32w_QC_fitter_norm_refine_v2_fixedcoarse/configs/config_v32w_QC_fitter_norm_refine_v2_fixedcoarse.in` |
| human labels (L20 x5, L24 x4) | `L20/input/v33p_*user_truezero_labels_*.csv` |
| accepted windows | `L20/output/accepted_windows_L20_user_validated_0p36/L20_<irrep>/accepted_windows.csv` |
| hot-window metadata/validation | `L20/output/hot_windows_L20_user_validated_0p36/L20_<irrep>/projected_basis_hotwindows{.json,_validation.csv}` (+ `.bin` from the full zip, not in Git) |
| lattice inputs used | `L20/input/lattice_spectrum/20_*.jack` (19 files, 603 samples each) |
| fit outputs (results to preserve) | `L20/output/l20_full_cutoff_fits_user_validated_0p36_passonly/`, `L20/output/l20_k3b_k3e_zero_tests/` |

### 9.4 Not needed in the consolidated fitter

All `debug_v29*`-`v31*`, `F3_device_resident_*`, `F3_gpu_omp_*`, `test_*`, `printer_*`,
`K3df_minuit_fit_v24*/v25/v26/v31l` headers, `scan_*` programs, and the `plot_v31*`/`v32*` plotting
scripts are historical; the cache generation programs belong to the cache-generator report.

### 9.5 Suggested consolidation design for a fast FCN (based on what works)

1. Keep the projected-basis decomposition $Q=\tilde F_r+\sum_i c_i\tilde B^{(i)}_r$; it is validated
   (matrix differences $\le 4.7\times10^{-10}$, 100 % sign agreement).
2. One cache format: per sector, rows $\to (\tilde F_r, \tilde B^{(1..4)}_r)$, scale applied, plain
   complex128, metadata with physics parameters and SHA-256; the hot-window file is a row subset of it.
3. Precompute $C^{+}$ once; parallelize windows with OpenMP; no stdout in the FCN.
4. Re-centre windows on the last accepted model zero and use a finite, smooth penalty when a root is
   lost, instead of $10^{100}$.
5. Unit tests: (a) basis sum vs direct assembly at a reference $\vec c$; (b) chi-square against an
   independent NumPy computation; (c) covariance from a synthetic per-configuration file reproduces
   $\sigma=\mathrm{s.d.}/\sqrt N$ (guards 6.1); (d) reader hermiticity and $\max|\operatorname{Im}\det|$.

---

## 10. What the results mean (summary)

- The fast-FCN design is sound and fast: an FCN on 14 levels costs about 24 ms and a full Minuit fit
  2-4 s, versus about 1.4 s per FCN and 250-400 s per fit with the earlier full-grid scan. The
  expensive part moved offline into the hot-window cache build (300-1570 s per sector).
- The roots are anchored by human labels; the in-FCN classifier only re-finds them locally.
- The L20 parameter values and chi-square values should not be quoted as physics results yet:
  the multi-L covariance appears to be inflated by $(N-1)^2$, which makes $\chi^2\sim10^{-4}$, stops Migrad
  almost at the start point, and inflates Hesse errors by about 600. The v32w fit
  ($\chi^2/\nu=1.44$) used a different (apparently correct) covariance construction.

## 11. Problems and assumptions

- I did not run any fit or build; timings and fit numbers are copied from output files.
- The covariance diagnosis (6.1) assumes the `.jack` files contain per-configuration values; this
  must be confirmed by Digonto. The Minuit2 EDM threshold value quoted is from my knowledge of Minuit2.
- The "two swaps cancel" conclusion (6.3) follows from reading the code; I did not read cache bytes.
- Dates of packages were taken from the brief, not from file metadata (all files carry the extraction
  mtime `Jul 15 10:53`).
- Commands I ran for numbers in this report:
  `awk 'NR>1{n++; s+=$2; ss+=$2*$2} END{...}' L20/input/lattice_spectrum/20_000_A1m_n0.jack` (and `20_001_A2_n0.jack`);
  `python3 -c` arithmetic for $(N-1)^2=362404$, rescaled $\chi^2$, $\Delta E=4.845\times10^{-6}$;
  `md5sum` / `diff` as in Section 2.

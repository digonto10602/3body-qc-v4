# KKpi L20 user-validated hot-window fitter

This release packages the L20 determinant/classifier validation, explicit
human true-zero windows, projected-basis hot-window caches, pass-only K3df
fits, true-zero spectra through Ecm=0.36, and chi-square/zero-test audits.

## Current status

`L20_PASSONLY_FITS_COMPLETE_ALL_TRUE_ZEROS_SAVED`

Completed fit cutoffs are 0.315, 0.325, 0.335, and 0.345. Cutoff 0.355 is
blocked by the L20/200_A2 lattice=4 versus user_true=5 mapping mismatch.
L24 is inactive.

## Included and excluded

Included: source, headers, Python/shell workflow scripts, configurations,
labels, the 19 user-true accepted windows, lattice jackknife inputs used by
the L20 fits, final hot-window cache metadata/validation and binaries, fit
outputs, zero-test outputs, reports, and representative plots.

Excluded from Git: the five large `.bin` cache payloads and all original full
F3inv/Vsel cache data. The `.bin` payloads are included in the full zip. The
original cache root must be restored separately at `/media/digonto/Data/F3inv_cache`
or configured through the fitter input.

## Layout

```text
source/                         C++/headers for reader, fitter, basis cache
scripts/                        validation, fitting, plotting, audit scripts
configs/                        historical/configured inputs
input/                          user labels and packaged L20 jackknife files
output/accepted_windows_*       only explicit user-true windows
output/hot_windows_*            projected basis matrices and metadata
output/l20_full_cutoff_*        four completed pass-only fits
output/l20_k3b_k3e_zero_tests   fixed/float zero tests and chi2 audit
plots/                          representative classifier and spectrum plots
tests/                           package smoke-test reports/logs
```

## Physics and production defaults

```text
QC = F3inv + K3df
projected_QC = Vsel.adjoint() * QC * Vsel
scaled_projected_QC = projected_QC / pow(Lbyas * xi, 6.0)
det = determinant(scaled_projected_QC)
```

No determinant-level Nproj rescaling is allowed. Defaults:
`xi=3.444`, `atmK=0.09698`, `atmpi=0.06906`, `eta_1=1`, `eta_2=0.5`,
`alpha=0.5`, `epsilon_h=0`, `max_shell_num=20`, `tolerance=1e-12`,
`parity=-1`, `eig_tol=0.05`, `norm_tol=1e-12`, `proj_tol=1e-10`,
`waves_vec_1=0 1`, `waves_vec_2=0`, `scatter1_00=4.04`,
`scatter1_10=-43.2`, `scatter2_00=4.12`.

Baseline K3df:

```text
73735.840894011912  -972421.14060757787  347174.05548116949  -1226756.7068845264
```

## Workflow

1. Scan the corrected raw cache and plot determinant families.
2. Generate candidate sign-change brackets.
3. Apply only explicit human `user_label=true` rows.
4. Build expanded hot-window projected-basis rows using
   `v33p_build_projected_basis_hotwindow_cache.cpp`.
5. Validate row coverage, hashes, matrix reconstruction, determinant signs,
   and metadata.
6. Fit with `root_search_mode=accepted_windows`,
   `projected_basis_cache_mode=hot_windows`, and `fallback_full_scan=0`.
7. Audit chi-square and preserve all fitted model true zeros.

Important functions/classes are indexed in `CODE_FUNCTION_INDEX.md`. The main
fitter entry is `v32x_multiL::MultiLFCN::operator()` in
`source/qc_fitter_norm_refine_v2_multiL.cpp`; chi-square is
`chi_square_v32f`; covariance inversion is
`symmetric_pseudoinverse_v32f`.

## Build

Dependencies: C++17 compiler, Eigen, OpenMP, Minuit2, Python 3, NumPy, and
Matplotlib. For the original cache reader, the external raw cache and its
metadata are also required. Run:

```bash
./build.sh
./RUN_SMOKE_TESTS.sh
```

The build does not run cachegen or fitting.

## Reproduction commands

```bash
./REPRODUCE_L20_PASSONLY.sh       # expensive; not part of smoke tests
./REPRODUCE_ZERO_TESTS.sh         # expensive; .355 remains excluded
python3 scripts/audit_l20_chi2_reconstruction.py
```

The scripts use only L20 accepted sectors and preserve the external cache
requirement. Do not run them merely to inspect this release.

## Plotting

Interactive determinant QA:

```bash
python3 scripts/plot_v33p_classifier_det_family_interactive.py \
  --grid <det_grid.csv> --candidate-file <candidate.csv> \
  --prediction-file <predictions.csv> --n-scale 100 --show
```

The spectrum plotter is `L20_spectrum_plotter_model_true_zeros_yaxis_fixed.py`.
It requires the historical old-QC and lattice/noninteracting inputs; the
release includes the new model-root CSV and a representative previously
generated PDF/PNG, but not absent historical old-isotropic files.

## Known blockers and next action

Resolve L20/200_A2's extra user-true window before attempting .355. Then
resume L24 determinant validation, resolve L24/100_A2 bracket 36, build only
user-validated L24 windows/caches, and run L24 prefit/one-FCN checks before a
combined fit. Never silently fall back to original full-cache reads inside an
FCN.

See `FINAL_SHUTDOWN_HANDOFF_REPORT.md`, `FULL_L20_FIT_AND_VALIDATION_REPORT.md`,
and `L24_NEXT_STEPS.md` for the full handoff.

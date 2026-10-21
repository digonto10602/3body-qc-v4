# Full L20 fit and validation report

## Status

`L20_PASSONLY_FITS_COMPLETE_ALL_TRUE_ZEROS_SAVED`

This report summarizes the reproducible L20 user-validated classifier,
projected-basis hot-window cache, fitting, zero-test, and chi-square audit
artifacts. L24 was not used.

## Classifier validation

The determinant-family plots show the scaled projected determinant as a family
of curves `det_real(Ecm)*10^k`, with sign-change candidates and classifier
markers overlaid. Human validation materializes one `user_label` per candidate.
Only explicit `user_label=true` rows become accepted windows; false and blank
rows are excluded. Both `digonto_v3_window` and `digonto_v4_window` were used
as the accepted classifier modes in the validated workflow.

The representative plot requested by the historical workflow,
`100_A2_L24_det.pdf`, is not present in this workspace. The packaged
representative is `plots/classifier_validation/L20_100_A2_classifier_det_family.pdf`.
Interactive generation is provided by `scripts/plot_v33p_classifier_det_family_interactive.py`.

## Accepted-window package

Source labels:

```text
input/v33p_user_truezero_labels_L20_000_A1m_E026310_0360.csv
input/v33p_user_truezero_labels_L20_100_A2_E026301_0360.csv
input/v33p_gpu_user_truezero_labels_L20_110_A2.csv
input/v33p_gpu_user_truezero_labels_L20_111_A2.csv
input/v33p_gpu_user_truezero_labels_L20_200_A2.csv
```

User-true counts are 3, 2, 4, 5, and 5 respectively: 19 accepted windows
through Ecm=0.36. The package preserves the original false rows but the
accepted-window package contains only true rows.

## Projected-basis hot-window caches

Each final cache stores `scaled_projected_F3inv` and the four projected K3df
basis matrices, together with row indices, Ecm/Nproj/dimension metadata,
accepted-window hashes, complex convention, and validation results. Row
coverage uses the expanded `max_row_left..max_row_right` ranges, including the
row-boundary fix that prevents the fitter from requesting an uncached endpoint.
Generation is OpenMP row-parallel with deterministic serial output ordering.

Actual final saved rows (including expansion/halo rows) are:

| sector | saved rows | matrix max diff | determinant relative diff | sign agreement |
|---|---:|---:|---:|---:|
| L20/000_A1m | 1507 | 4.66e-10 | 2.23e-9 | 1507/1507 |
| L20/100_A2 | 1006 | 4.55e-13 | 5.43e-10 | 1006/1006 |
| L20/110_A2 | 2012 | 4.66e-10 | 1.01e-5 | 2012/2012 |
| L20/111_A2 | 2511 | 4.66e-10 | 6.52e-5 | 2511/2511 |
| L20/200_A2 | 2509 | 2.91e-11 | 4.69e-4 | 2509/2509 |

The metadata reports accepted-window hash match and zero missing requested
rows. The older logical window counts were 1503, 1002, 2004, 2505, and 2505;
the larger saved counts are the required expanded fitter coverage.

## Physics invariant

```text
QC = F3inv + K3df
projected_QC = Vsel.adjoint() * QC * Vsel
scaled_projected_QC = projected_QC / pow(Lbyas * xi, 6.0)
det = determinant(scaled_projected_QC)
```

Raw determinant in the scaled projected basis is the classifier/fitter
quantity. No Nproj or determinant-level rescaling is applied.

## L20 pass-only fits

| cutoff | levels | model_found | chi2/dof | K3iso0 | K3iso1 | K3B | K3E | roots <=0.36 | status |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|---|
| .315 | 6 | 6/6 | 1.6311e-5 | 151422.46 | -972421.14 | 347174.06 | -1226756.71 | 19 | PASS |
| .325 | 9 | 9/9 | 1.9062e-5 | 125968.22 | -972421.14 | 347174.06 | -1226756.71 | 19 | PASS |
| .335 | 12 | 12/12 | 1.6041e-5 | 32040.84 | -1074716.65 | 347174.06 | -1226756.71 | 19 | PASS |
| .345 | 14 | 14/14 | 1.6711e-5 | 40607.52 | -1067567.87 | 376187.14 | -1226756.71 | 19 | PASS |

All four are actual minimization results using projected-basis hot-window mode,
with fallback/full-scan `0/0`. The .355 fit was not run because L20/200_A2 has
lattice=4 versus user_true=5.

## Old isotropic versus new s+p spectrum

The old isotropic QC files named in the historical workflow were not present
in this workspace, nor were the three old-vs-new CSV comparison files. The
new model-root source is
`output/l20_full_cutoff_fits_user_validated_0p36_passonly/L20_PASSONLY_TRUE_ZEROS_ALL_WINDOWS_ALL_CUTOFFS.csv`;
the packaged spectrum plot is `plots/spectrum/KKpi_L20_model_spectrum_fitcut_0p335.pdf`.
Consequently, an exact level-by-level old-isotropic shift table cannot be
reconstructed from this package alone. The prior analysis reported that levels
below .335 shifted mostly upward by O(1e-4), while the largest visible
above-cutoff shift was the high 000_A1m level moving downward by about 0.00108.

## Chi-square audit

The current objective is raw covariance chi-square,
`(E_lattice-E_model)^T C^+ (E_lattice-E_model)`, with
`chi2/dof=chi2/(nlevels-nfloating)`. Independent full-covariance reconstruction
matched the fitter output to about 1e-14 relative for completed zero tests.
The tiny ~1.6e-5 values are explained by residuals around 1e-4--1e-3 and
covariance scales around 0.3--0.4; they are not evidence of double division.
The old isotropic ~1.5 value has an unverified objective convention, so the
comparison is not apples-to-apples. Report the new quantity as
`fitter-reported raw covariance chi2/dof` until the legacy fit is recomputed
with the same covariance and residual convention.

Exact source points are documented in
`output/l20_k3b_k3e_zero_tests/L20_K3B_K3E_ZERO_TEST_FINAL_REPORT.md`.

## K3B/K3E zero test

The fixed-zero branch (`K3B=K3E=0`, only K3iso0/K3iso1 floating) passes .315
and .325 but loses one required root at .335. The float-from-zero branch also
keeps K3B and K3E at zero in the completed .315/.325 tests; its HESSE errors
near 2 are suspiciously degenerate and not meaningful significance estimates.
Neither branch reached .345 after the .335 root loss. The validated four-
parameter pass-only branch is the only tested branch that succeeds through
.345.

## Errors and correlations

Per-cutoff `final_parameters.csv`, `parameter_errors.csv`, covariance, and
correlation files are preserved under the pass-only and zero-test output roots.
Pass-only HESSE status is recorded as VALID, but unusually small errors in some
cutoffs, especially .325, require caution. Zero-test floating K3B/K3E errors
near 2 are clearly degenerate diagnostics, not production uncertainties.
K3E remains fixed/unmoved in the tested branches; K3iso0 absorbs most of the
fixed-zero compensation.

Pass-only HESSE error values recorded in `L20_PASSONLY_FIT_SUMMARY.csv` are:

| cutoff | sigma(K3iso0) | sigma(K3iso1) | sigma(K3B) | sigma(K3E) | covariance status |
|---:|---:|---:|---:|---:|---|
| .315 | 342579611.35 | 128986425.01 | 134948646.91 | 8786443671.24 | VALID, ill-conditioned scale |
| .325 | 2.00 | 2.00 | 2.00 | 2.00 | VALID flag, suspicious/degenerate |
| .335 | 113211204.38 | 84176025.16 | 219464496.69 | 608399252.10 | VALID, large uncertainties |
| .345 | 32989189.93 | 224279530.06 | 317856751.35 | 1116819088.83 | VALID, large uncertainties |

The native covariance and correlation files are under each cutoff directory,
for example `Ecm_cut_0p335/l20_full_0p335_parameter_covariance.dat` and
`l20_full_0p335_parameter_correlation.dat`. The zero-test covariance and
correlation tables are under `output/l20_k3b_k3e_zero_tests/`; their
float-from-zero K3B/K3E errors near 2 are diagnostic-only.

## Exact reproducibility inputs and outputs

The fitter reads L20 jackknife files from the original directory
`/home/digonto/Codes/KKpi_I2/spectrum/Ecm_data/data`; the 19 files used by the
four completed cutoffs are copied into `input/lattice_spectrum/` in this
package. The external raw cache is `/media/digonto/Data/F3inv_cache`.

Key files are:

```text
output/l20_full_cutoff_fits_user_validated_0p36_passonly/L20_PASSONLY_FIT_SUMMARY.csv
output/l20_full_cutoff_fits_user_validated_0p36_passonly/L20_PASSONLY_TRUE_ZEROS_ALL_WINDOWS_ALL_CUTOFFS.csv
output/l20_full_cutoff_fits_user_validated_0p36_passonly/L20_PASSONLY_LATTICE_VS_MODEL_ALL_CUTOFFS.csv
output/l20_k3b_k3e_zero_tests/L20_K3B_K3E_ZERO_TEST_FINAL_REPORT.md
output/l20_k3b_k3e_zero_tests/INDEPENDENT_CHI2_RECONSTRUCTION.csv
```

## L24 next steps

1. Finish or verify L24 determinant scans.
2. Validate L24 windows manually from determinant plots.
3. Promote only explicit user-true rows.
4. Keep L24/100_A2 excluded until bracket 36 is resolved.
5. Build the L24 accepted-window package and expanded hot-window caches.
6. Run L24 prefit coverage and one-FCN checks.
7. Only then run L24 or combined L20+L24 fits.
8. Never include false/blank windows or allow original-cache fallback.

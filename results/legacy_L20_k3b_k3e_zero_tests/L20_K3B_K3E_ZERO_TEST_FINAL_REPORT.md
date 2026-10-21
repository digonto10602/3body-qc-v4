# L20 K3B/K3E zero tests and chi-square audit

CHECKPOINT_STATUS = `L20_K3B_K3E_ZERO_TESTS_COMPLETE_WITH_ROOT_LOSS_AT_0P335`

L20 only; cutoffs .315, .325, .335, .345 requested; L24 inactive; .355 not run; projected-basis hot-window mode only; fallback/full_scan=0/0 in completed fits.

## Table A

| test | cutoff | initial chi2/dof | final chi2/dof | FCN evals | model_found | status |
|---|---:|---:|---:|---:|---|---|
| fixed_K3B0_K3E0 | 0.315 | 1.1013989783328644e-05 | 1.0542666800363224e-05 | 93 | 6/6 | PASS |
| fixed_K3B0_K3E0 | 0.325 | 1.2765039793309932e-05 | 1.247421508822199e-05 | 74 | 9/9 | PASS |
| fixed_K3B0_K3E0 | 0.335 | 1e+99 | nan | 1 | 11/12 | BLOCKED_MODEL_FOUND |
| fixed_K3B0_K3E0 | 0.345 | nan | nan | 0 | 0/14 | BLOCKED_MODEL_FOUND |
| float_from_K3B0_K3E0 | 0.315 | 2.202797956665729e-05 | 2.039453905725892e-05 | 225 | 6/6 | PASS |
| float_from_K3B0_K3E0 | 0.325 | 1.754065510253554e-05 | 1.692892130183815e-05 | 175 | 9/9 | PASS |
| float_from_K3B0_K3E0 | 0.335 | 1.25e+99 | nan | 1 | 11/12 | BLOCKED_MODEL_FOUND |
| float_from_K3B0_K3E0 | 0.345 | nan | nan | 0 | 0/14 | BLOCKED_MODEL_FOUND |

## Table B

| test | cutoff | reported chi2/dof | independent chi2/dof | diagonal-only chi2/dof | full-cov chi2/dof | conclusion |
|---|---:|---:|---:|---:|---:|---|
| fixed_K3B0_K3E0 | 0.315 | 1.0542666800363224e-05 | 1.0542666800363573e-05 | 1.7023018238863943e-06 | 1.0542666800363573e-05 | PASS |
| fixed_K3B0_K3E0 | 0.325 | 1.247421508822199e-05 | 1.2474215088222325e-05 | 2.7697293276560407e-06 | 1.2474215088222325e-05 | PASS |
| fixed_K3B0_K3E0 | 0.335 | nan | nan | nan | nan | BLOCKED_OR_NOT_RUN |
| fixed_K3B0_K3E0 | 0.345 | nan | nan | nan | nan | BLOCKED_OR_NOT_RUN |
| float_from_K3B0_K3E0 | 0.315 | 2.039453905725892e-05 | 2.0394539057259572e-05 | 3.3073140717690547e-06 | 2.0394539057259572e-05 | PASS |
| float_from_K3B0_K3E0 | 0.325 | 1.692892130183815e-05 | 1.692892130183859e-05 | 3.921556870125577e-06 | 1.692892130183859e-05 | PASS |
| float_from_K3B0_K3E0 | 0.335 | nan | nan | nan | nan | BLOCKED_OR_NOT_RUN |
| float_from_K3B0_K3E0 | 0.345 | nan | nan | nan | nan | BLOCKED_OR_NOT_RUN |

## Table C

| test | cutoff | K3B_final | sigma_K3B | pull_K3B | K3E_final | sigma_K3E | pull_K3E | moved_from_zero? |
|---|---:|---:|---:|---:|---:|---:|---:|---|
| fixed_K3B0_K3E0 | 0.315 | 0.0 | 0.0 | NA | 0.0 | 0.0 | NA | no |
| fixed_K3B0_K3E0 | 0.325 | 0.0 | 0.0 | NA | 0.0 | 0.0 | NA | no |
| fixed_K3B0_K3E0 | 0.335 | 0.0 | 0.0 | NA | 0.0 | 0.0 | NA | no |
| fixed_K3B0_K3E0 | 0.345 | 0.0 | 0.0 | NA | 0.0 | 0.0 | NA | no |
| float_from_K3B0_K3E0 | 0.315 | 0.0 | 1.99999988079072 | 0.0 | 0.0 | 1.99999988079072 | 0.0 | no |
| float_from_K3B0_K3E0 | 0.325 | 0.0 | 1.9999998807907167 | 0.0 | 0.0 | 1.9999998807907156 | 0.0 | no |
| float_from_K3B0_K3E0 | 0.335 | nan | nan | NA | nan | nan | NA | no |
| float_from_K3B0_K3E0 | 0.345 | nan | nan | NA | nan | nan | NA | no |

## Table D

| cutoff | old isotropic chi2/dof ~1.5 | previous 4-param chi2/dof | fixed-zero chi2/dof | float-from-zero chi2/dof | apples-to-apples? |
|---:|---:|---:|---:|---:|---|
| 0.315 | 1.5 | 1.6310700374147595e-05 | 1.0542666800363224e-05 | 2.039453905725892e-05 | NO versus legacy isotropic; YES versus previous pass-only definition |
| 0.325 | 1.5 | 1.9061793244739764e-05 | 1.247421508822199e-05 | 1.692892130183815e-05 | NO versus legacy isotropic; YES versus previous pass-only definition |
| 0.335 | 1.5 | 1.604077078944618e-05 | nan | nan | NO versus legacy isotropic; YES versus previous pass-only definition |
| 0.345 | 1.5 | 1.6711008405490147e-05 | nan | nan | NO versus legacy isotropic; YES versus previous pass-only definition |

## Chi-square normalization and FCN audit

The shared FCN calls `chi_square_v32f` with `chi_square_mode=raw_cov_inv`; that function forms the raw Ecm residual and contracts it with the Eigen self-adjoint pseudoinverse of the jackknife covariance. The fitter reports `chi2` directly, then computes `chi2_dof=chi2/(nlevels-nfloating)`. Fixed parameters are excluded from `nfloating`. No second division was found. The trace and independent reconstruction agree for completed fits.

The tiny values are numerically explained by the saved data: residuals are roughly 1e-4--1e-3 Ecm while covariance diagonal square roots are roughly 0.3--0.4. This is a real value for this covariance, not evidence that the old isotropic ~1.5 is comparable. The exact legacy isotropic objective/source was not found in the current tree, so future reports should call the current quantity `raw covariance chi2/dof` until that legacy path is audited.

Fixed-zero tests compensate mainly through K3iso0; float-from-zero .315/.325 did not move K3B or K3E away from zero, but their HESSE errors are suspiciously degenerate (~2) and should not be interpreted as meaningful significance. Both tests lose one required root at .335 from the prescribed sequential start; .345 is consequently not run.

## Exact source audit points

- `source/qc_fitter_norm_refine_v2_multiL.cpp:1359-1365`: `MultiLFCN::operator()` builds the model, calls `chi_square_v32f`, applies the missing-root penalty, and writes the FCN/level traces.
- `source/K3df_minuit_fit_v32f_fullF3inv_QCfull_cached_classifier.hpp:381-393`: `symmetric_pseudoinverse_v32f` uses a symmetrized Eigen self-adjoint eigendecomposition with `rcond=1e-12`.
- `source/K3df_minuit_fit_v32f_fullF3inv_QCfull_cached_classifier.hpp:474-512`: `chi_square_v32f`; raw mode contracts `(E_lattice-E_model)` with the covariance pseudoinverse.
- `source/qc_fitter_norm_refine_v2_multiL.cpp:1604-1615`: `ndof=ndata-npar` and `chi2_dof=chi2/ndof`.
- `source/lattice_data_covariance_cpp.hpp:421-425`: jackknife covariance/correlation prefactor is `(Njack-1)/Njack`.

## Required output files

- `ZERO_TEST_FCN_TRACE_ALL.csv` and `ZERO_TEST_FCN_LEVEL_TRACE_ALL.csv`
- `INDEPENDENT_CHI2_RECONSTRUCTION.csv/.md`
- `CHI2_COMPARISON_TO_OLD_ISOTROPIC.csv/.md`
- per-cutoff `fcn_trace.csv`, `fcn_level_trace.csv`, and `chi2_scale_audit.csv/.md`
- `ZERO_TEST_FIT_SUMMARY_COMPARISON.csv/.md` and `ZERO_TEST_ROOT_SHIFT_COMPARISON.csv/.md`

No commit/push was performed.

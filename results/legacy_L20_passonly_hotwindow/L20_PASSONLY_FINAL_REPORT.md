# L20 user-validated pass-only fit final report

CHECKPOINT_STATUS = L20_PASSONLY_FITS_COMPLETE_ALL_TRUE_ZEROS_SAVED

## Scope

L20 only; sectors `000_A1m`, `100_A2`, `110_A2`, `111_A2`, and `200_A2`. Projected-basis hot-window mode was used with no original-cache fallback and no full scan. L24 remained inactive.

## Results

See `L20_PASSONLY_FIT_SUMMARY.csv`. All four cutoffs 0.315, 0.325, 0.335, and 0.345 converged with complete model_found counts and 0/0 fallback/full scan.

## All fitted-model true zeros

For each fitted cutoff, `model_true_zeros_all_windows.csv` contains all 19 explicit user=true windows through Ecm=0.36. Roots are obtained by linear interpolation of the raw determinant between bracketing sampled rows; `PASS_LINEAR_DET_ROOT` is intentionally recorded. No logdet is used.

## Required artifacts

- `L20_PASSONLY_TRUE_ZEROS_ALL_WINDOWS_ALL_CUTOFFS.csv`
- `L20_PASSONLY_MODEL_ECM_LEVELS_ALL_CUTOFFS.csv`
- `L20_PASSONLY_LATTICE_VS_MODEL_ALL_CUTOFFS.csv`
- `L20_PASSONLY_PARAMETER_DRIFT.csv`

## Blocked next cutoff

0.355 was not run: L20/200_A2 has four lattice levels but five user-true windows, requiring explicit mapping resolution.

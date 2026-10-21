# L20 pass-only completed checkpoint

CHECKPOINT_STATUS = `L20_PASSONLY_FITS_COMPLETE_ALL_TRUE_ZEROS_SAVED`

## Scope

L20 only, using projected-basis hot-window mode and no original-cache fallback.
L24 remained inactive. No fitting or cache generation is part of this checkpoint push.

Fitted cutoffs: `0.315`, `0.325`, `0.335`, `0.345`.

Skipped cutoff: `0.355`.

Reason skipped: L20/200_A2 has lattice=4 versus user_true=5, so its one-to-one lattice/window mapping requires resolution before fitting.

## Fit results

| cutoff | model_found | chi2/dof | K3iso0 | K3iso1 | K3B | K3E | true roots through 0.36 |
|---:|---|---:|---:|---:|---:|---:|---:|
| 0.315 | 6/6 | 1.6310700374147595e-05 | 151422.46172441466 | -972421.1406075779 | 347174.0554811695 | -1226756.7068845264 | 19 |
| 0.325 | 9/9 | 1.9061793244739764e-05 | 125968.22184657381 | -972421.1406075779 | 347174.0554811695 | -1226756.7068845264 | 19 |
| 0.335 | 12/12 | 1.6040770789446181e-05 | 32040.84189923505 | -1074716.6512728098 | 347174.05548116943 | -1226756.7068845264 | 19 |
| 0.345 | 14/14 | 1.6711008405490147e-05 | 40607.515418119 | -1067567.866187696 | 376187.1404411546 | -1226756.7068845264 | 19 |

All fits have `fallback/full_scan = 0/0`. The 19 roots per fitted cutoff are restricted to explicit user-validated true windows and saved through Ecm=0.36.

Covariance files are finite positive-definite Minuit-HESSE outputs. The 0.325 errors are unusually small and should be treated cautiously.

## Key artifacts

- `L20_PASSONLY_FINAL_REPORT.md`
- `L20_PASSONLY_FIT_SUMMARY.csv`
- `L20_PASSONLY_TRUE_ZEROS_ALL_WINDOWS_ALL_CUTOFFS.csv`
- `L20_PASSONLY_LATTICE_VS_MODEL_ALL_CUTOFFS.csv`

Per-cutoff parameter, covariance, model-root, and lattice-vs-model CSV/JSON/MD artifacts are stored under the four `Ecm_cut_*` directories. Large projected-basis binaries, original external F3inv/Vsel caches, build artifacts, and logs are intentionally excluded from Git.

## Next action

Resolve the L20/200_A2 one-to-one lattice/window mapping before attempting the 0.355 fit.

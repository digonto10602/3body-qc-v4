# Final shutdown handoff: L20 user-validated hot-window fitter

## Checkpoint

`CHECKPOINT_STATUS = L20_PASSONLY_FITS_COMPLETE_ALL_TRUE_ZEROS_SAVED`

This checkpoint contains completed L20 projected-basis hot-window fits at
Ecm cutoffs 0.315, 0.325, 0.335, and 0.345. Every completed fit found all
selected model roots and used `fallback_full_scan=0` and `full_scan=0`.

The 0.355 fit was intentionally not run: L20/200_A2 has four lattice levels
but five explicit user-true windows. L24 was inactive and must remain inactive
until its determinant scans and user validation are complete.

## Restore

The clean release tree is `release/KKpi-L20-user-validated-hotwindow-fitter/`.
The full local package is `KKpi-L20-user-validated-hotwindow-fitter_FULL_WORKING_PACKAGE.zip`.
The original full raw caches remain external at `/media/digonto/Data/F3inv_cache/`.
Copy the zip contents to a new workspace, set the external cache root in the
fitter configuration, install Eigen, OpenMP, Minuit2, Python, NumPy, and
Matplotlib, then run `build.sh` and `RUN_SMOKE_TESTS.sh`.

## Physics and defaults

```text
QC = F3inv + K3df
projected_QC = Vsel.adjoint() * QC * Vsel
scaled_projected_QC = projected_QC / pow(Lbyas * xi, 6.0)
det = determinant(scaled_projected_QC)
```

No determinant-level Nproj rescaling is used. Production defaults are
`xi=3.444`, `atmK=0.09698`, `atmpi=0.06906`, `eta_1=1`, `eta_2=0.5`,
`alpha=0.5`, `epsilon_h=0`, `max_shell_num=20`, `tolerance=1e-12`,
`parity=-1`, `eig_tol=0.05`, `norm_tol=1e-12`, `proj_tol=1e-10`,
`waves_vec_1=0 1`, `waves_vec_2=0`, `scatter1_00=4.04`,
`scatter1_10=-43.2`, and `scatter2_00=4.12`.

Starting K3df values are:

```text
K3iso0 = 73735.840894011912
K3iso1 = -972421.14060757787
K3B    = 347174.05548116949
K3E    = -1226756.7068845264
```

## Completed artifacts

Accepted L20 sectors are `000_A1m`, `100_A2`, `110_A2`, `111_A2`, and
`200_A2`. The accepted-window root is
`output/accepted_windows_L20_user_validated_0p36/`; the hot-window projected
basis root is `output/hot_windows_L20_user_validated_0p36/`; pass-only fits
are under `output/l20_full_cutoff_fits_user_validated_0p36_passonly/`; and
zero tests are under `output/l20_k3b_k3e_zero_tests/`.

The five final cache binaries are included in the zip but excluded from the
normal Git commit because three exceed GitHub's 100 MB file limit. Their JSON
metadata, validation CSVs, and reports are tracked. The package does not
contain `/media/digonto/Data/F3inv_cache`.

The five user label files are the two `input/v33p_user_truezero_labels_*`
files for 000_A1m/100_A2 and the three `input/v33p_gpu_user_truezero_labels_*`
files for 110_A2/111_A2/200_A2. The L20 lattice jackknife source is
`/home/digonto/Codes/KKpi_I2/spectrum/Ecm_data/data`; the 19 used files are
also copied under `input/lattice_spectrum/` in the release package.

Pass-only K3df results are:

| cutoff | model_found | chi2/dof | K3iso0 | K3iso1 | K3B | K3E |
|---:|---:|---:|---:|---:|---:|---:|
| .315 | 6/6 | 1.6311e-5 | 151422.46 | -972421.14 | 347174.06 | -1226756.71 |
| .325 | 9/9 | 1.9062e-5 | 125968.22 | -972421.14 | 347174.06 | -1226756.71 |
| .335 | 12/12 | 1.6041e-5 | 32040.84 | -1074716.65 | 347174.06 | -1226756.71 |
| .345 | 14/14 | 1.6711e-5 | 40607.52 | -1067567.87 | 376187.14 | -1226756.71 |

Each cutoff has 19 fitted model true-zero rows through Ecm=0.36. The
zero-test report is diagnostic: fixed-zero passes .315/.325 but loses a root
at .335; float-from-zero leaves K3B/K3E at zero with suspicious HESSE errors.

## Restart commands

```bash
cd release/KKpi-L20-user-validated-hotwindow-fitter
./build.sh
./RUN_SMOKE_TESTS.sh
python3 scripts/audit_l20_user_validated_prefit.py
python3 scripts/audit_l20_chi2_reconstruction.py
```

Do not run the reproduction wrappers until the external raw cache and lattice
paths are restored. Do not resume L24 or attempt .355 without resolving the
L20/200_A2 mapping mismatch.

## Exact next action

Resolve the L20/200_A2 one-to-one lattice/window mapping, then rerun only the
0.355 prefit audit and fit if the mapping is unambiguous. Separately finish
L24 determinant scans, manually validate its windows, exclude L24/100_A2 until
the bracket-36 false positive is resolved, and only then prepare combined L20+L24
fits. Never promote false or blank user labels and never enable original-cache
fallback inside the FCN.

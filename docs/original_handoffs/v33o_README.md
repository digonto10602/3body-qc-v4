# KKpi v33o production fitter

This package is the packaged v33f/v33n production fitter snapshot for the K3df multi-`Lbyas` workflow.

## What this package does

- loads the v33g runtime caches and the v33e/v33f fitter binary
- validates target-loader parity, classifier dispatch, and parameter masks
- runs the v33n cutoff scan for `Ecm_cutoff = 0.315, 0.325, 0.335, 0.345, 0.355, 0.365`
- writes QC spectra, spectrum plots, determinant plots, and summary reports
- packages the source, scripts, configs, docs, and validation reports without bundling large generated caches

## Canonical workflow

1. Compile the fitter.
2. Validate target counts.
3. Run the benchmark-fcn smoke test.
4. Use `scripts/run_v33n_cutoff_scan_fits.py` to reproduce the cutoff scan.
5. Inspect `reports/v33n_production_sanity_cutoff_scan_full_report.md` and `reports/v33n_cutoff_scan/cutoff_scan_summary.md`.

## Build and smoke test

```bash
bash scripts/compile_v33g_all.sh
python3 scripts/validate_fitter_target_counts_v33k.py
stdbuf -oL -eL bin/v33f_k3df_fitter_multiL_v33e configs/config_v33m_multiL_all_irreps_fastest_production.in benchmark-fcn --repeat 10 --warmup 2
```

## Main scripts

- `scripts/run_v33n_cutoff_scan_fits.py`: orchestrates the v33n cutoff scan, fit attempts, QC spectrum generation, spectrum plots, determinant plots, and reports.
- `scripts/compile_v33g_all.sh`: builds the fitter binary.
- `scripts/validate_fitter_target_counts_v33k.py`: checks the target-loader parity against the production config.
- `scripts/verify_classifier_dispatch_v33m.py`: confirms classifier dispatch behavior.
- `scripts/run_v33m_classifier_speed_search.py`: benchmark helper for classifier search modes.
- `scripts/plot_classifier_algorithm_sweep_v33l.py`: classifier sweep plotting and diagnostics.

## Logic flow

`run_v33n_cutoff_scan_fits.py` follows this order:

- compile the binary
- run target parity validation
- run a warm 10-FCN benchmark
- dump targets for each cutoff
- try spectrum-only validation first
- run the full fit for each cutoff
- write QC spectrum outputs and the non-interacting overlay file
- generate spectrum and determinant plots
- write per-cutoff reports and the combined cutoff summary
- write the final production report and start-here handoff

## Important parameters

- canonical production config: `configs/config_v33m_multiL_all_irreps_fastest_production.in`
- v33n cutoff values: `0.315, 0.325, 0.335, 0.345, 0.355, 0.365`
- default irreps: `000_A1m`, `100_A2`, `110_A2`, `111_A2`, `200_A2`
- default `Lbyas`: `20`, `24`

## How to modify the workflow

- Edit `configs/config_v33m_multiL_all_irreps_fastest_production.in` to change the canonical production seed.
- Edit `scripts/run_v33n_cutoff_scan_fits.py` to change the cutoff list, target selection, plotting, or report layout.
- Edit `source/qc_fitter_norm_refine_v2_multiL.cpp` if you need fitter logging or parameter-mask behavior changes.
- Keep the physics formulas and classifier rules unchanged unless you are intentionally changing the scientific model.

## Generated outputs

- `output_v33n/cutoff_scan/Ecmcut_*/`: converged QC spectra, used lattice spectrum, and non-interacting line lists
- `plots/v33n_cutoff_scan/Ecmcut_*/`: spectrum and determinant plots
- `reports/v33n_cutoff_scan/`: per-cutoff and combined production reports
- `diagnostics/v33o/`: v33o verification CSVs

## Known limits

- The current recent non-interacting cache stops at L24.
- Minuit warnings can appear in converged fits.
- Historical v33n cutoff wall times were backfilled from saved logs; future runs should parse wall time directly.

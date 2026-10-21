# KKpi v33o Production Fitter

This repository packages the current v33n checkpoint of the KKpi / K3df multi-`Lbyas` workflow.

Current state:

- combined raw GPU cache reader bug fixed and validated
- corrected old-scale determinant scan passed for `Lbyas=20`, `irrep=100_A2`, `Ecm=[0.34, 0.36]`, `coarseN=20000`
- CPU OpenMP subset validation passed against the corrected cache path
- corrected `fix1` interactive plot regenerated from the right CSV
- classifier restart bundle prepared, waiting only on user-labeled true zeros

The packaging goal is reproducibility, not new physics or new fitting rules.

## What this repository contains

- source code for the fitter, cache builders, classifier utilities, and plotting helpers
- configs for the v33f, v33g, v33m, and related workflows
- validation reports and checkpoint notes
- packaging docs for the v33o handoff

Large generated caches, build trees, and binary outputs are intentionally excluded.

## Critical checkpoint

The key fix was the combined raw `F3inv+Vsel` GPU cache reader. The target cache is:

```text
/media/digonto/Data/F3inv_cache/v33e_coarse20000_sep_f2k2_L20_100_A2/v32zu_Lbyas20_gpu_cache/v32zu_Lbyas20_xi3p444_irrep100_A2_coarse20000_F3inv_Vsel_gpu.bin
```

The corrected path reconstructs the raw payload consistently and yields:

- `max |det_imag| = 0`
- `F3inv` hermiticity error at roundoff
- `projected_QC` hermiticity error at roundoff
- CPU/GPU subset agreement at sub-percent level

## Exact old scaling

Use matrix-level old scaling only:

```cpp
const double const_norm_scale = std::pow(Lbyas * xi, 6.0);
Eigen::MatrixXcd QC = F3inv + K3df;
Eigen::MatrixXcd projected_QC = Vsel.adjoint() * QC * Vsel;
Eigen::MatrixXcd scaled_projected_QC = projected_QC / const_norm_scale;
std::complex<double> det = scaled_projected_QC.determinant();
```

Do not apply determinant-level `Nproj` rescaling.

## Validation numbers

The current checkpoint is the one documented in:

- `reports/v33o_final_production_verification.md`
- `reports/v33n_final_sanity_target_parity.md`
- `reports/v33n_final_sanity_classifier_dispatch.md`
- `reports/v33n_final_sanity_parameter_mask.md`
- `reports/v33n_final_sanity_10fcn.md`

Representative values:

- `model_found = 32/32`
- `chi2 = 0.2543510436151185` for the 10-FCN smoke benchmark
- `max percent_diff_F3inv_frob ~ 3.2005e-05 %`
- `max percent_diff_QC_frob ~ 3.2005e-05 %`
- `max percent_diff_projected_QC_frob ~ 1.3901e-05 %`
- `max percent_diff_scaled_projected_QC ~ 1.3901e-05 %`
- `max stable percent_diff_det_complex ~ 6.0955e-07 %`
- `max projector_reldiff ~ 9.35e-16`

## Build and smoke test

```bash
bash scripts/compile_v33g_all.sh
python3 scripts/validate_fitter_target_counts_v33k.py
stdbuf -oL -eL bin/v33f_k3df_fitter_multiL_v33e configs/config_v33m_multiL_all_irreps_fastest_production.in benchmark-fcn --repeat 10 --warmup 2
```

## Main workflow commands

Run the production cutoff scan:

```bash
python3 scripts/run_v33n_cutoff_scan_fits.py
```

Generate the runtime cache validation plots:

```bash
python3 scripts/plot_v33g_runtime_qc_det_compare.py
```

Validate classifier dispatch:

```bash
python3 scripts/verify_classifier_dispatch_v33m.py
```

## Current classifier state

The classifier is not being tuned in this checkpoint. It is parked at:

- `output/v33n_classifier_restart_from_fix1_L20_100_A2/`

That bundle contains candidate sign-change brackets and segment summaries for the eight candidate roots in the corrected `fix1` CSV.

## Cache format notes

Two related cache formats exist in this repo:

1. Combined raw `F3inv+Vsel` GPU cache
   - file suffix: `_F3inv_Vsel_gpu.bin`
   - used by the corrected reader path
   - stores raw complex payloads for `F3inv` and `Vsel`

2. Runtime cache
   - file suffix: `_runtime_K3basis.bin`
   - stores projected matrices and K3 basis pieces for the hot fitter path
   - paired with `.meta.json` metadata

The reader must be format-aware. File name alone is not enough.

## Repository layout

- `source/`: fitter, cache, and classifier implementation
- `scripts/`: build, validation, plotting, and scan drivers
- `configs/`: runtime and production configs
- `docs/`: packaging notes and code-map documents
- `reports/`: validation reports and checkpoint notes
- `diagnostics/`: small verification CSVs

## Reproducibility

Large runtime caches live outside the repo and are mounted from:

```text
/media/digonto/Data/F3inv_cache
```

The repo itself carries the source, scripts, configs, docs, and validation evidence needed to recreate the checkpoint.

## Known limits

- L25 non-interacting data is not present in the current recent cache tree.
- This package does not include large binary cache artifacts.
- The classifier still needs user-labeled true zeros before any tuning step.

## See also

- `QUICKSTART.md`
- `CODE_STRUCTURE.md`
- `FUNCTION_MAP.md`
- `VALIDATION_STATUS.md`
- `GPU_CACHE_READER_FIX.md`
- `CLASSIFIER_RESTART_NOTES.md`
- `REPRODUCIBILITY.md`

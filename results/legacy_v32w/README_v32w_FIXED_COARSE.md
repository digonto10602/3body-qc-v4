# v32w QC fitter: fixed coarse grid + separate refinement cache

This package fixes the v32v bug where refinement points were merged back into the initial coarse sign-flip grid. In v32w:

- `coarse_cache` is fixed at exactly `coarseN` rows per irrep and is the only grid used for the initial sign-flip scan.
- `refined_cache` stores newly built refinement F3inv/Vsel points and is reused inside refinement only.
- The saved output cache is a merged cache (`coarse + refined`) for future reuse, but on the next run it is split again into fixed coarse rows and refined rows.
- The final recomputed chi-square is forced to `failure_penalty` if the final saved spectrum does not contain all model levels.

Run:

```bash
mamba activate minuit2-fit
export OMP_NUM_THREADS=18
export EIGEN_INC=/usr/include/eigen3
bash run.sh fit 10000
```

Important outputs:

```text
output_v32w_QC_fitter_norm_refine_v2_fixedcoarse/debug_v32w_QC_fitter_norm_refine_v2_fixedcoarse_fit_summary.dat
output_v32w_QC_fitter_norm_refine_v2_fixedcoarse/debug_v32w_QC_fitter_norm_refine_v2_fixedcoarse_fit_levels.dat
output_v32w_QC_fitter_norm_refine_v2_fixedcoarse/debug_v32w_QC_fitter_norm_refine_v2_fixedcoarse_bestfit_QC_spectrum.dat
output_v32w_QC_fitter_norm_refine_v2_fixedcoarse/debug_v32w_QC_fitter_norm_refine_v2_fixedcoarse_all_QC_candidates.dat
output_v32w_QC_fitter_norm_refine_v2_fixedcoarse/debug_v32w_QC_fitter_norm_refine_v2_fixedcoarse_parameter_covariance.dat
output_v32w_QC_fitter_norm_refine_v2_fixedcoarse/debug_v32w_QC_fitter_norm_refine_v2_fixedcoarse_parameter_correlation.dat
```

Cache output:

```text
cache/v32w_merged_F3inv_Vsel_cache.bin
```


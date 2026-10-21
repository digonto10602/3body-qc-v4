# KKpi v32w + v32y reproducibility package

This package consolidates the two latest usable KKpi workflow versions plus tools, example outputs, and an archive of older available packages.

## What is inside

```text
latest_working_versions/
  v32w_QC_fitter_norm_refine_v2_fixedcoarse/
    Trusted K3df fitter baseline. Uses fixed coarse cache separation and normalized determinant refinement.
  v32y_QC_spectrum_vs_L_streamcache_fixedK3df_compilefix/
    Memory-bounded spectrum-vs-L runner. Processes one Lbyas/irrep cache at a time and writes per-irrep cache files.

tools/
  cache inspectors, state-line interpolation tools, and plotting scripts.

outputs_examples/
  Small output zips from v32w/v32x/v32y and generated continuous QC state-line files.

old_versions_archive/
  All version/package zips available in this conversation from v31/v32 plus noninteracting/source packages.

docs/
  Source/reference PDF and this package’s version guides.
```

## Environment

These commands assume your usual local environment:

```bash
mamba activate minuit2-fit
export OMP_NUM_THREADS=18
export EIGEN_INC=/usr/include/eigen3
```

The scripts auto-detect the conda/mamba compiler and Minuit2 library from `$CONDA_PREFIX` when possible.

## Main trusted fitter: v32w

Use this version for fitting `K3df` to the L=20 ensemble. It is the latest workably good fitter baseline.

```bash
cd latest_working_versions/v32w_QC_fitter_norm_refine_v2_fixedcoarse
mamba activate minuit2-fit
export OMP_NUM_THREADS=18
export EIGEN_INC=/usr/include/eigen3
bash run.sh fit 10000
```

Equivalent explicit script:

```bash
bash run_v32w_QC_fitter_norm_refine_v2_fixedcoarse.sh \
  fit \
  10000 \
  configs/config_v32w_QC_fitter_norm_refine_v2_fixedcoarse.in
```

Known good v32w fit result:

```text
chi2 = 11.489748255611568
ndata = 12
npar = 4
ndof = 8
chi2/dof = 1.436218531951446
model_found = 12/12

K3iso0 = 34470.228772793678
K3iso1 = -1035173.9140440958
K3B    = 253921.14757953485
K3E    = -961680.91711685632
```

Main v32w output directory:

```text
output_v32w_QC_fitter_norm_refine_v2_fixedcoarse/
```

Best-fit spectrum file:

```text
output_v32w_QC_fitter_norm_refine_v2_fixedcoarse/debug_v32w_QC_fitter_norm_refine_v2_fixedcoarse_bestfit_QC_spectrum.dat
```

## Memory-bounded spectrum-vs-L: v32y

Use this version to generate QC spectra as a function of integer `Lbyas` using fixed K3df parameters. This version replaces memory-heavy v32x for long L scans.

Run only Lbyas=24 first:

```bash
cd latest_working_versions/v32y_QC_spectrum_vs_L_streamcache_fixedK3df_compilefix
mamba activate minuit2-fit
export OMP_NUM_THREADS=18
export EIGEN_INC=/usr/include/eigen3
bash run.sh spectrum_vs_L 10000 24 24 \
  34470.228772793678 \
  -1035173.9140440958 \
  253921.14757953485 \
  -961680.91711685632
```

Full integer Lbyas scan, Lbyas=16...24:

```bash
bash run.sh spectrum_vs_L 10000 16 24 \
  34470.228772793678 \
  -1035173.9140440958 \
  253921.14757953485 \
  -961680.91711685632
```

v32y uses the cache directory from its config. The default is:

```text
/media/digonto/Data/F3inv_cache
```

v32y cache naming is per L and per irrep:

```text
v32y_Lbyas24_xi3.444_irrep110_A2_coarse10000_F3inv_Vsel_merged.bin
```

This is intentionally different from v32x, which used one very large all-irrep merged cache.

Main v32y output directory:

```text
output_v32y_QC_spectrum_vs_L_streamcache_fixedK3df/
```

## Cache inspection

Check one cache file:

```bash
python3 tools/inspect_v32_f3inv_cache_v32f_fixed.py \
  --cache /media/digonto/Data/F3inv_cache/v32x_Lbyas20_xi3.444_coarse10000_F3inv_Vsel_merged.bin \
  --irrep 110_A2 \
  --coarse-n 10000 \
  --ecm-min 0.261 \
  --ecm-max 0.360 \
  --write-energies cache_110_A2_energies.dat
```

For a big cache, use the progress/mmap inspector:

```bash
python3 tools/inspect_v32_f3inv_cache_mmap_progress.py \
  --cache /media/digonto/Data/F3inv_cache/v32x_Lbyas24_xi3.444_coarse10000_F3inv_Vsel_merged.bin \
  --irrep 110_A2 \
  --coarse-n 10000 \
  --ecm-min 0.261 \
  --ecm-max 0.360 \
  --progress-every 250
```

Interpretation:

```text
coarse_present = 10000 and missing_coarse = 0  => coarse grid complete
refined_or_extra > 0                          => refinement points were saved
refined_or_extra = 0                          => coarse-only cache, or no sign flips/refinement for that irrep
```

## Continuous QC state-line interpolation

Convert scattered QC spectra into continuous line-like files over Lbyas:

```bash
python3 tools/make_qc_state_lines_from_scattered_spectrum.py \
  outputs_examples/output_v32x_QC_spectrum_vs_L_fixedK3df.zip \
  outputs_examples/output_v32y_QC_spectrum_vs_L_streamcache_fixedK3df.zip \
  --outdir qc_state_lines_from_v32x_v32y \
  --lmin 16 \
  --lmax 24 \
  --interp pchip \
  --dense-n 401 \
  --ymax 0.36 \
  --plot
```

Output columns:

```text
label track_id interp_method n_observed Lbyas spatial_L Ecm
```

Use `--interp linear` for strict piecewise-linear state lines, or `--interp pchip` for cubic-like smooth lines with fewer overshoot problems than ordinary cubic splines.

## Plotting helpers

The `tools/` folder includes the recent spectrum overlay scripts, including the fixed-left-y-axis version:

```bash
python3 tools/L20_old_lattice_recent_QC_spectrum_plotter_recent_nonint_L16to24_yaxis_fixed.py \
  --cutoff 0.335 \
  --recent-fit-file output_v32w_QC_fitter_norm_refine_v2_fixedcoarse/debug_v32w_QC_fitter_norm_refine_v2_fixedcoarse_bestfit_QC_spectrum.dat \
  --nonint-script /path/to/nonint3body_general_doublecover_package/source/nonint3body_general_doublecover.py \
  --output output_plots/KKpi_L16to24_oldQC_L20_recentFit_recentNonint_yaxis_fixed.png \
  --output-pdf output_plots/KKpi_L16to24_oldQC_L20_recentFit_recentNonint_yaxis_fixed.pdf \
  --show
```

## Version guidance

- Use `v32w` for actual K3df fitting.
- Use `v32y` for memory-bounded QC spectrum vs L generation with fixed K3df parameters.
- Keep `v32x` only as a reference. It can build huge all-irrep cache files and can fill RAM/swap.
- Use `tools/inspect_*` before rerunning large jobs to confirm whether cache exists and is complete.
- Use `tools/make_qc_state_lines_from_scattered_spectrum.py` to turn scattered Lbyas spectra into line-like branches.

See `docs/VERSION_SUMMARY.md`, `latest_working_versions/*/QUICKSTART_LOCAL.md`, and `old_versions_archive/README_OLD_VERSIONS.md` for more details.

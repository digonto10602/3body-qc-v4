# Tools quick guide

## Cache inspectors

- `inspect_v32_f3inv_cache_v32f_fixed.py`: main binary cache inspector; handles `success` as int32.
- `inspect_v32_f3inv_cache_mmap_progress.py`: progress-printing inspector for very large cache files.
- `inspect_v32_f3inv_cache_robust.py`: robust parser that reports incomplete/corrupt offsets.
- `inspect_v32_f3inv_cache.py`: earlier inspector kept for reference.

Example:

```bash
python3 inspect_v32_f3inv_cache_v32f_fixed.py \
  --cache /media/digonto/Data/F3inv_cache/v32x_Lbyas20_xi3.444_coarse10000_F3inv_Vsel_merged.bin \
  --irrep 110_A2 \
  --coarse-n 10000 \
  --ecm-min 0.261 \
  --ecm-max 0.360 \
  --write-energies cache_110_A2_energies.dat
```

## Continuous QC state-line maker

- `make_qc_state_lines_from_scattered_spectrum.py`: reads scattered QC spectrum outputs and tracks/interpolates branches over Lbyas.

Example:

```bash
python3 make_qc_state_lines_from_scattered_spectrum.py \
  ../outputs_examples/output_v32x_QC_spectrum_vs_L_fixedK3df.zip \
  ../outputs_examples/output_v32y_QC_spectrum_vs_L_streamcache_fixedK3df.zip \
  --outdir qc_state_lines_from_v32x_v32y \
  --lmin 16 --lmax 24 --interp pchip --dense-n 401 --plot
```

## Plotting scripts

The `L20_old_lattice_recent_*` scripts overlay old QC, L=20 lattice data, recent v32w fit, and noninteracting levels. The latest fixed-y-axis version is:

```text
L20_old_lattice_recent_QC_spectrum_plotter_recent_nonint_L16to24_yaxis_fixed.py
```

## v32x all-irrep cache readiness checker

- `check_v32x_L24_cache_ready.py`: checks a large v32x-style all-irrep merged cache for stable file size, parseability, complete coarse grid, and failed entries.

Example:

```bash
python3 check_v32x_L24_cache_ready.py \
  --cache /media/digonto/Data/F3inv_cache/v32x_Lbyas24_xi3.444_coarse10000_F3inv_Vsel_merged.bin \
  --coarse-n 10000 \
  --ecm-min 0.261 \
  --ecm-max 0.360
```

This does not generate the QC spectrum; it only checks cache readiness.

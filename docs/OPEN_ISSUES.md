# OPEN ISSUES — read before trusting any number

Issues found while consolidating the code (2026-10-21). Evidence and line numbers are in `reports/research/001`–`006`.
Severity:
- **BLOCKER**: invalidates results.
- **HIGH**: must be fixed before the next production fit.
- **MED/LOW**: should be fixed; does not invalidate results.

Nothing listed here has been fixed in the code yet. The code in `src/` is copied verbatim from the newest packages.

**Digonto's answers to decisions D1–D5 (2026-10-21)** are recorded in `notes/STATUS.md`; the issues below are annotated accordingly.

---

### 1. BLOCKER (CONFIRMED by D1) — lattice covariance too large by $(N-1)^2=602^2\approx3.6\times10^5$

- **Where:** `src/fitter/qc_fitter_norm_refine_v2_multiL.cpp`, `load_multil_jack_targets` (about lines 329–448).
- **What:** it applies the jackknife formula $C=\frac{N-1}{N}\sum_a(x_a-\bar x)(x_a-\bar x)$ directly to the 603 values in each `data/lattice/L20/*.jack` file. The older reader `src/fitter/lattice_data_covariance_cpp.hpp` (used by the trusted v32w fit) first builds leave-one-out averages, i.e. it treats the file values as **per-configuration** values.
- **Evidence:**
  1. The resulting error on the lowest level is $\sigma\approx0.31$, about 115 % of the energy itself ($a_tE\approx0.269$). The per-configuration reading gives $\sigma\approx5\times10^{-4}$.
  2. Neighbouring file entries differ by about 0.01, typical of single configurations, not of jackknife samples.
  3. The ratio of the two readings is exactly 602 for all 19 files.
  4. In every L20 fit, $K_E$ never moved from its start value, Hesse errors are $10^8$–$10^9$, and $\chi^2/\nu\sim10^{-5}$.
- **Re-evaluation:** `python3 tools/audit/chi2_covariance_convention_check.py .` re-evaluates the *legacy* best-fit model energies (`results/legacy_L20_passonly_hotwindow/`) with the per-configuration covariance:

  | cutoff $a_tE^*$ | levels | $\chi^2$ as fitted | $\chi^2$ per-config cov | $\chi^2/\nu$ |
  |---|---|---|---|---|
  | 0.315 | 6 | $3.3\times10^{-5}$ | 11.8 | 5.9 |
  | 0.325 | 9 | $9.5\times10^{-5}$ | 34.5 | 6.9 |
  | 0.335 | 12 | $1.3\times10^{-4}$ | 46.5 | 5.8 |
  | 0.345 | 14 | $1.7\times10^{-4}$ | 60.6 | 6.1 |

  These are not fit results. The fits never converged with the inflated covariance, so they must be redone.
- **Answer (D1, 2026-10-21):** the `.jack` files hold the energy ("mass") of one state on each of $N$ gauge configurations: $N=603$ for L20 and $N=553$ for L24. They are **per-configuration values**, so this issue is real.
- **Fix:** build leave-one-out jackknife averages first, as `lattice_data_covariance_cpp.hpp` does. Add a unit test: a synthetic per-configuration file must reproduce $\sigma=\text{s.d.}/\sqrt N$. Then redo all L20 fits.
- **Consequence:** every χ², parameter value and parameter error in `results/legacy_L20_*` and in the v33n cutoff scan is invalid. The v32w fit ($\chi^2/\nu=1.44$) used the older reader and is not affected by this issue.

### 2. HIGH — p-wave $\mathcal K_2^{-1}$ differs from the reference implementation

- **Where:** `src/core/K2_functions_v2.h:102,106`, and the GPU copy `src/cachegen_gpu/K2_functions_gpu_safe.cuh:414-415`.
- **What:** these compute $B_1=-1/(a_1q^{*2})+\dots$ plus $C_1=q^{*2}|q^*|(1-H)$. The two terms carry different powers of momentum.
- **Reference:** QC3_release (`base_code/F3/K2i_mat.py:26-27`) uses $-1/a_1+q^{*2}|q^*|(1-H)$, which is the commented-out line 103.
- **Consequence:** `scatter1_10=-43.2` does not mean the reference's $a_1$, and the near-threshold p-wave behaviour is different. The size of the effect is not yet quantified.
- **Answer (D2, 2026-10-21):** the two-body inputs come from Digonto's earlier work and from Hadron Spectrum Collaboration papers.
- **Still open (planner research):** find those papers and the parametrisation in which $a_1=-43.2$ is quoted, then use the matching $\mathcal K_2^{-1}$ form. Until then, the coded form is the production convention and must not be changed silently.
- **Planner finding (2026-10-21, `reports/research/007` §3.6):** the values are HadSpec's $I=3/2$ $\pi K$ global fit (arXiv:1411.2004 Sec. VI, ERE $k^{2\ell+1}\cot\delta_\ell=1/a_\ell+\dots$): $a_0=(-4.04\pm0.05\pm0.15)\,a_t$, $a_1=(43.2\pm3.7\pm15.4)\,a_t^3$. In the code's $Q_0$ convention (2111.12734 Eqs. 15, 17, 20–21) the consistent element is $-1/a_1+q^{*2}|q^*|(1-H)$ with the existing $a_1=-43.2$ (the commented-out line 103), i.e. QC3_release's form. The coded form overstates the p-wave $\mathcal K_2^{-1}$ by $1/q^{*2}$. Both forms are now named (`coded_v2`, `q0_consistent`, `docs/METHODS.md` §5); production stays `coded_v2` until **decision D6** (after `prompts/001` step 9 reports the level shifts).

### 3. HIGH — complex-number read convention (`variant_04_real_imag_swapped`)

- **What the writer does:** the GPU cache writer stores $(\mathrm{Re},\mathrm{Im})$ in natural order.
- **What the readers do:**
  - `gpu_cache_reader_convention::read_complex` (`src/cache_io/gpu_cache_reader_convention.hpp:179`) returns `comp(b,a)`, i.e. swaps once.
  - The trusted scanner `v33h_…det_scan.cpp` and the fitter then call `swap_real_imag` again, so the net read is natural and **correct**.
- **Origin of the confusion:** a debug reader `comp(read_f64(), read_f64())`. C++ does not fix the order in which the two reads happen.
- **Copies to watch:**
  - The L20 copy of `v33g_build_runtime_k3basis_cache.cpp` swaps only once. This repo therefore ships the **v33n copy**, which reads naturally.
  - `src/cache_io/compare_gpu_cache_vs_cpu_openmp_det_oldscale.cpp` swaps once. Do not use its read path; it is kept only because `v33h` `#include`s it for helpers.
- **Fix:** one plain reader, e.g. NumPy `np.frombuffer(buf,'<c16')`, plus tests:
  - $\|F_3^{-1}-F_3^{-1\dagger}\|/\|F_3^{-1}\|\approx10^{-15}$;
  - $\max|\mathrm{Im}\det|=0$;
  - $\|V^\dagger V-1\|\approx10^{-15}$.

### 4. HIGH — inputs needed to re-run fits are not in this repo

| missing | needed for | where it was |
|---|---|---|
| raw cache root (`*_F3inv_Vsel_gpu.bin` + `.meta.json`, about 2.6–3.3 GB per sector) | building det grids and hot-window caches | `/media/digonto/Data/F3inv_cache/v33e_coarse20000_sep_f2k2_L<L>_<irrep>/` |
| hot-window binaries `projected_basis_hotwindows.bin` (about 150 MB per sector) | running the fast fitter | "full zip" of the L20 release; only the JSON/CSV metadata is here (`data/hot_windows_meta/`) |
| 20,000-row det-grid CSVs and candidate-bracket CSVs (all sectors except L20/100_A2 brackets) | fitter metadata (`det_grid_file_*`, `candidate_brackets_file_*`) | `main_fitter/.../output/v33p_*` on the source machine |
| L24 lattice spectrum `24_*.jack` | combined L20+L24 fit | `/home/digonto/Codes/KKpi_I2/spectrum/Ecm_data/data` |
| old `V20_*`, `QC_states_jackknifed_*` | spectrum plotter overlays | — |

**Answer (D4, 2026-10-21):** these binaries will not be provided. The new direction is to make the QC evaluation fast enough that the cache files are not needed (`prompts/000`, "Update 2026-10-21"). The L24 `.jack` files are still needed for M3.

### 5. HIGH — level matching and label-count mismatches

- **Matching rule:** the fitter pairs the $k$-th accepted zero with the $k$-th lattice level and requires equal counts per sector. Any missed low zero shifts every match above it.
- **Rule from Digonto (D3, 2026-10-21):** sort both lists ascending. With $n$ QC levels and $m\le n$ lattice levels, use only the lowest $m$ QC levels; the extra QC levels are not an error.
  - The fitter still needs this change: it currently requires equal counts.
  - $n<m$ (a missing QC root) remains an error.
- **Known blockers:**
  - L20/200_A2 has 5 user-true windows but 4 lattice levels at cutoff 0.355.
  - L20/000_A1m has 3 vs 2 at 0.36.
  - L24/100_A2 bracket 36 is unresolved.
- **Effect of D3 on the blockers:** the first two blockers (5 vs 4 and 3 vs 2) are resolved by the rule. Bracket 36 still needs a root decision.

### 6. MED — FCN design issues (`src/fitter/qc_fitter_norm_refine_v2_multiL.cpp`)

- **Discontinuous χ²:** a lost root returns $\chi^2=10^{100}$. Migrad assumes a smooth function.
- **Fixed windows:** windows are never re-centred on the current model zero (`previous_model_zero` is fixed).
- **Logging:** "fallback/full_scan 0/0" is a hard-coded log string (about line 1371), not a measurement.
- **Speed:**
  - $C^+$ is recomputed every FCN call;
  - the window loop is serial;
  - the classifier prints one line per window.
- **Parameter mask:** `SetLimits(i,…)` loops over all 4 parameters even when some are fixed (about line 2068). The parameter-covariance labels assume the fixed parameters come last (about lines 1624, 2091).

### 7. MED — physics inputs without provenance or uncertainty

- **Isospin:** confirmed by Digonto (D5): KKπ with $I=2$, i.e. $K^+K^+\pi^+$. The 2+1 identical-pair treatment applies.
- **Two-body parameters:** 4.04, −43.2 and 4.12 come from Digonto's earlier work and from HadSpec papers (D2). The exact citations and the uncertainties are still to be written down. Older drivers carried $a_0^{KK}=4.07$ with "0.07 stat" and $a_0^{K\pi}=4.04$ with "0.05 stat, 0.15 syst".
- **Fixed constants:** $\xi$, $a_tm_K$, $a_tm_\pi$ are constants; no uncertainty is propagated.

### 8. MED — cache generator robustness (`src/cachegen_gpu/F3_gpu_cachegen_v32zu.cu`)

- **Row order:** rows are grouped by matrix size, not sorted by energy, and there is no index. Readers must sort by `grid_i`.
- **Missing metadata:** parity, `eig_tol` and the scattering parameters are compiled into the binary (not run-time settable); parity and `eig_tol` are also absent from the sidecar metadata. `projection_source` names `generate_swave_cpu_cache_qc.cpp` (kept for reference in `legacy/cachegen_reference/`).
- **Silent failures:** rows are written even when the LU reports `info>0` (singular).
- **CPU vs GPU:** the GPU $F$ uses a power-series erfi and plain summation; the CPU uses Faddeeva and Neumaier summation. Agreement on projected objects was reported as $\le3.2\times10^{-5}$ % (v33n reports), but only on a subset.
- **Memory guard:** in the all-object build it always forces single-task waves, so streams never overlap.

### 9. MED — GPU determinant backend not validated

`src/gpu_det/` (v33p) consumes CPU-assembled matrices. It has not passed CPU-reference parity, so its results must not drive fits (`docs/original_handoffs/v33p_GPU_*`).

### 10. LOW — minor physics-core bugs (not on the production path)

- `functions.h:1564`: `threebody_Gpoles` uses `pole_energy2*pole_energy1` where `pole_energy2` squared is intended.
- `config_maker_2` puts zeros into the wrong rows (`functions.h:473-475`).
- The `Q0norm=false` path in $G$ divides by $q^{*2\ell}$ instead of $q^{*\ell}$. Production uses `Q0norm=true`.
- The `F3iso` diagnostic vector puts weight on p-wave rows.

### 11. LOW — portability

- Scripts and configs contain absolute paths (`/home/digonto/Codes/...`, `/media/digonto/Data/F3inv_cache`).
- `build.sh`/`build_gpu.sh` have not been run in this layout; Eigen and Minuit2 were not available on the consolidation machine.

### 12. Historical warnings (do not reuse these results)

- v32e–v32p read the L=24 lattice files while computing at L=20, so their fits are invalid.
- The v33n cutoff-scan fits used the `raw_sign_only` classifier, later ruled "not production".
- Old (pre-`F2_functions_v2.h` 2376-line) p-wave $F$ had a wrong $\sqrt{\pi/\alpha^3}$ factor and a missing $\sqrt2$ on $G_{12}$. Any p-wave numbers from `3body_quantization_v3/source` (726-line version) are wrong.

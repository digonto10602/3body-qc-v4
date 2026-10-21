# D: Input data and parameters for the KKπ three-body QC fitter

> **Research report — consolidation investigation.** Date: 2026-10-21. Git commit: none (written before the repo existed; first commit of `3body-qc-v4`). Prompt in force: `prompts/000-consolidation-and-roadmap.md`. Related: `reports/research/001`–`006`.
> Absolute paths below refer to the unpacked original zips (`/home/digimonk/Projects/three-body-qc-codes/extracted/...`) and to a temporary scratchpad. For each file's new location in this repo, see `docs/PROVENANCE.tsv`.
> The system is two kaons + one pion treated as a 2+1 system (identical kaons). Where this report says "K K-bar pi", read "KKπ"; that wording came from the investigation brief and is not established. See `docs/OPEN_ISSUES.md` #7.


- **Date:** 2026-10-21
- **Git commit:** n/a. `extracted/` is not a git repository; everything below refers to the unpacked zips under `/home/digimonk/Projects/three-body-qc-codes/extracted/`.
- **Prompt in force:** coordinator brief "Area D: input data and parameters" (no `prompts/NNN` file exists yet for this consolidation).
- **Related reports:** the sibling area reports (A–C, E…) in the same `scratchpad/reports/` folder.
- **Mode:** read only. The only things I ran were cheap Python/NumPy reads of the 19 small `.jack` files, plus one χ² recomputation on 14 levels (script: `scratchpad/scripts/chi2_check.py`). No cache generation or fits.

Short names used below:
- `PKG` = `/home/digimonk/Projects/three-body-qc-codes/extracted/KKpi-L20-user-validated-hotwindow-fitter-main/KKpi-L20-user-validated-hotwindow-fitter-main`. This is the newest package (July 2026) and the one with the most complete data.
- `REPRO` = `/home/digimonk/Projects/three-body-qc-codes/extracted/KKpi-v32w-v32y-repro-master/KKpi-v32w-v32y-repro-master`
- `NEST` = `/home/digimonk/Projects/three-body-qc-codes/extracted/_nested_old_versions`

Terms used in this report:
- **Jackknife:** a resampling method. With $N$ gauge configurations you build $N$ "leave-one-out" averages. The spread of these averages gives the statistical error, and their joint spread gives the covariance between energy levels.
- **Ensemble file / "rescaled" ensemble:** the HadSpec/ADAT file format (header `Ncfg Nt type ...`). It stores one number per configuration. These numbers do **not** have to be jackknife samples. See §3.2.
- **Irrep:** an irreducible representation of the cubic little group, e.g. `100_A2` means total momentum $\vec P = \frac{2\pi}{L}(1,0,0)$, irrep $A_2$.
- **Ecm:** centre-of-mass energy in temporal lattice units, $a_t E_{\rm cm}$.
- **Lbyas:** spatial box size $L/a_s$ (20 or 24).
- **True-zero window:** an $E_{\rm cm}$ bracket in which a person (Digonto) confirmed that the sign change of the QC determinant is a real zero, not a pole. "Window", "bracket" and "hot window" are this codebase's own labels, not physics terms.

---

## 1. What was asked

1. Inventory the lattice spectrum data (jackknife files): format, ensembles, irreps, sample count, reader code, covariance construction, units, cutoff selection, duplicates. Say which files make up the complete fitting set.
2. Explain the format and meaning of the user true-zero label CSVs and the accepted windows.
3. Enumerate and explain the config keys of the latest production configs, the physics defaults, and where the scattering parameters come from.
4. Give the exact $K_{\rm df,3}$ parametrization ($K_{3,\rm iso}^{(0)}, K_{3,\rm iso}^{(1)}, K^B, K^E$) and the baseline values.
5. Explain the non-interacting-energy packages: which is newest and which the plotters use.
6. List the plotting and analysis scripts worth keeping.
7. Explain what `QC3_release` and `source` (from `source.zip`) in `_nested_old_versions` are.

## 2. Actions taken

- Listed every package. Searched all of `extracted/` for `*.jack`, `*state*`, `szscl21*`, `V20_*`, label CSVs, `*.in` and `*.bin`.
- Read the data readers:
  - `PKG/source/lattice_data_covariance_cpp.hpp` (legacy reader, used by the single-L v32w fitter).
  - `load_multil_jack_targets` in `PKG/source/qc_fitter_norm_refine_v2_multiL.cpp`, lines 273–447 (the reader in production use).
- Parsed all 19 `.jack` files with NumPy (`python3 -I`). Computed means and errors under both possible sample conventions.
- Recomputed the reported χ² of the four L20 "pass-only" fits from `Ecm_cut_*/lattice_vs_model_levels.csv` plus the raw `.jack` files.
- Read the K2 code (`K2_functions_v2.h`, `K2_functions_gpu_safe.cuh`), the K3 code (`K3_functions_2plus1.hpp`), the config parser (`settings_from_config` in `qc_fitter_norm_refine_v2.cpp`, `multiconfig_from_config` in `..._multiL.cpp`) and the QC3_release Python (`K2i_mat.py`, `K3main.py`).
- Checked duplicates with `md5sum`. The supplied `hashes.txt` does **not** cover `_nested_old_versions` (0 lines match), so I hashed those files directly.

---

## 3. Findings

### 3.1 Lattice spectrum data: what exists

**Only one ensemble's per-level data is present anywhere in `extracted/`:** 19 files, all at $L/a_s = 20$.

`PKG/input/lattice_spectrum/` (136 407 bytes total, about 7.2 kB per file):

- File name pattern: `20_<P><irrep>_n<k>.jack`. The frame label uses z-ordering (`001`, `002`, `011`, `111`). The fitter maps these with `internal_alias`: `001→100`, `002→200`, `011→110`.
- No other copy exists. Each `.jack` file name is unique across the whole tree.
- The original location was `/home/digonto/Codes/KKpi_I2/spectrum/Ecm_data/data/`, on Digonto's machine and not packaged. The `I2` in that path presumably means isospin 2 (see §5).

**File format.** One header line, `603 1 0 0 1`, followed by 603 lines of `0 <value>`:
- In the ADAT ensemble header convention: $N_{\rm cfg}=603$, time extent 1, real data.
- Column 1 is always 0 (the time-slice index). Column 2 is the energy.
- The energy is the **lab-frame** energy $a_t E_{\rm lab}$. The config key `lattice_jack_energy_type = En_lab` says so, and $E_{\rm lab}=E_{\rm cm}$ holds exactly for the `000` files.

**Conversion to the centre-of-mass frame**, done per sample in `convert_samples_to_ecm`:

$$a_tE_{\rm cm}=\sqrt{(a_tE_{\rm lab})^2-(a_tP)^2},\qquad a_tP=\frac{2\pi}{\xi\,(L/a_s)}\,|\vec n_P|,\qquad \xi=3.444 .$$

**Table 1: all 19 levels.**
- $\sigma_{\rm JK}$ is the correct jackknife error *if* the file holds per-configuration ("rescaled") values: $\sigma=\mathrm{std}/\sqrt N$.
- $\sigma_{\rm fitter}$ is what the production multi-L fitter actually uses (see §3.2).
- Command: the inline `python3 -I` loop over `*.jack` in this report's actions.

| file | $\lvert\vec n_P\rvert^2$ | $a_tP$ | mean $a_tE_{\rm lab}$ | mean $a_tE_{\rm cm}$ | $\sigma_{\rm JK}$ ($10^{-3}$) | $\sigma_{\rm fitter}$ |
|---|---|---|---|---|---|---|
| `20_000_A1m_n0.jack` | 0 | 0.000000 | 0.269227 | 0.26923 | 0.515 | 0.3102 |
| `20_000_A1m_n1.jack` | 0 | 0.000000 | 0.345606 | 0.34561 | 0.506 | 0.3047 |
| `20_001_A2_n0.jack` | 1 | 0.091219 | 0.305566 | 0.29161 | 0.506 | 0.3045 |
| `20_001_A2_n1.jack` | 1 | 0.091219 | 0.313625 | 0.30004 | 0.586 | 0.3525 |
| `20_001_A2_n2.jack` | 1 | 0.091219 | 0.374492 | 0.36320 | 0.519 | 0.3123 |
| `20_002_A2_n0.jack` | 4 | 0.182439 | 0.339872 | 0.28664 | 0.612 | 0.3684 |
| `20_002_A2_n1.jack` | 4 | 0.182439 | 0.350626 | 0.29933 | 0.587 | 0.3533 |
| `20_002_A2_n2.jack` | 4 | 0.182439 | 0.377970 | 0.33094 | 0.613 | 0.3691 |
| `20_002_A2_n3.jack` | 4 | 0.182439 | 0.392604 | 0.34738 | 1.093 | 0.6578 |
| `20_011_A2_n0.jack` | 2 | 0.129004 | 0.333348 | 0.30733 | 0.525 | 0.3161 |
| `20_011_A2_n1.jack` | 2 | 0.129004 | 0.341758 | 0.31644 | 0.516 | 0.3109 |
| `20_011_A2_n2.jack` | 2 | 0.129004 | 0.343375 | 0.31814 | 0.774 | 0.4661 |
| `20_011_A2_n3.jack` | 2 | 0.129004 | 0.353563 | 0.32915 | 0.519 | 0.3122 |
| `20_011_A2_n4.jack` | 2 | 0.129004 | 0.393186 | 0.37139 | 0.597 | 0.3594 |
| `20_111_A2_n0.jack` | 3 | 0.157996 | 0.357254 | 0.32036 | 0.563 | 0.3390 |
| `20_111_A2_n1.jack` | 3 | 0.157996 | 0.369695 | 0.33413 | 0.772 | 0.4649 |
| `20_111_A2_n2.jack` | 3 | 0.157996 | 0.373642 | 0.33854 | 0.572 | 0.3441 |
| `20_111_A2_n3.jack` | 3 | 0.157996 | 0.378740 | 0.34416 | 0.582 | 0.3503 |
| `20_111_A2_n4.jack` | 3 | 0.157996 | 0.384992 | 0.35103 | 0.576 | 0.3465 |

**Thresholds** (from the masses below):
- $a_t(2m_K+m_\pi)=0.26302$ (the $KK\pi$ threshold).
- $a_t(2m_K+2m_\pi)=0.33208$ (the $KK\pi\pi$ threshold, drawn in the plotter).

**Levels kept per cutoff.** The fitter keeps a level when the jackknife mean of its $E_{\rm cm}$ samples is $\le$ `Ecm_cutoff`. My Table 1 reproduces the fitter's counts exactly:

| `Ecm_cutoff` | 000_A1m | 100_A2 | 110_A2 | 111_A2 | 200_A2 | total | fitter log |
|---|---|---|---|---|---|---|---|
| 0.315 | 1 | 2 | 1 | 0 | 2 | 6 | 6 ✓ |
| 0.325 | 1 | 2 | 3 | 1 | 2 | 9 | 9 ✓ |
| 0.335 | 1 | 2 | 4 | 2 | 3 | 12 | 12 ✓ |
| 0.345 | 1 | 2 | 4 | 4 | 3 | 14 | 14 ✓ |
| 0.355 | 2 | 2 | 4 | 5 | 4 | 17 | blocked (200_A2 mismatch) |
| 0.360 | 2 | 2 | 4 | 5 | 4 | 17 | — |

Two levels are never used below 0.36: `20_001_A2_n2` (0.36320) and `20_011_A2_n4` (0.37139).

**Ensemble identity.**
- The legacy reader's default ensemble string is `szscl21_20_128_b1p50_t_x4p300_um0p0840_sm0p0743_n1p265`: $L/a_s=20$, $T/a_t=128$, $\beta=1.50$, bare anisotropy 4.3, light/strange bare masses −0.0840/−0.0743.
- The masses in use are $a_tm_\pi=0.06906$ and $a_tm_K=0.09698$. From memory, these match HadSpec's published values for its $m_\pi\approx391$ MeV anisotropic ensembles. This is my recollection and should be checked against the HadSpec papers.

**Other ensembles: L=24 data is missing.**
- The production config `config_v33m_multiL_all_irreps_fastest_production.in` uses `Lbyas_values = 20 24` and its header says "32 targets" at cutoff 0.335. That means 20 L=24 levels existed on the original machine.
- The only trace of them in the tree is three means in `3body_quantization_v3-master/.../output/v33p_minimal_fitter_000_A1m/LEVEL_ZERO_MATCH_REPORT.md`:
  - `24_000_A1m_n0` = 0.266883
  - `24_000_A1m_n1` = 0.320893
  - `24_000_A1m_n2` = 0.329533
- **No `24_*.jack` file is present.**

**Legacy-format data is also missing.** The legacy `szscl21_..._<irrep>_state_<k>` files and the `V20_<irrep>` summary files (needed by the plotter) are not present either. The legacy path was `.../3body_quantization{,_v2}/lattice_data/KKpi_interacting_spectrum/twoptvar_analysis/masses/`.

**Complete set needed for fitting.**
- **L=20:** the 19 files above are complete through $E_{\rm cm}=0.36$.
- **L=20+L=24 (the stated production target):** the L=24 `.jack` files must be retrieved from `/home/digonto/Codes/KKpi_I2/spectrum/Ecm_data/data/24_*.jack` on Digonto's machine. At 0.335 there should be 20 of them below cutoff; how many exist above that is unknown.

### 3.2 Critical: the production fitter inflates the data errors by a factor of 602

Two readers exist. They treat the file values differently.

**(a) Legacy reader** `covariance_between_states_szscl21_based_result` (`PKG/source/lattice_data_covariance_cpp.hpp`, used by the single-L v32w fitter). It treats the file values $x_c$ as per-configuration data. It first builds jackknife samples, then the jackknife covariance:

$$y_c=\frac{1}{N-1}\sum_{c'\neq c}x_{c'},\qquad C_{ij}=\frac{N-1}{N}\sum_{c=1}^{N}\big(y^{(i)}_c-\bar y^{(i)}\big)\big(y^{(j)}_c-\bar y^{(j)}\big).$$

The conversion $E_{\rm lab}\to E_{\rm cm}$ is applied to each jackknife sample.

**(b) Multi-L reader** `load_multil_jack_targets` (`..._multiL.cpp` lines 330–447). This is used by every v33f–v33p and L20-hot-window fit. It treats the file values directly as jackknife samples. There is no resampling step:

$$C_{ij}=\frac{N-1}{N}\sum_{c}\big(x^{(i)}_c-\bar x^{(i)}\big)\big(x^{(j)}_c-\bar x^{(j)}\big)\quad\text{(computed per }L\text{; cross-}L\text{ blocks set to }0\text{)}.$$

**The evidence says the `.jack` files are per-configuration ("rescaled") ensembles, so (a) is right and (b) is wrong.**

1. The sample spread inside one file is about 0.01 in $a_tE$, while neighbouring levels are about 0.01 apart. With reader (b), every level would have $\sigma\approx0.30$–$0.66$ (Table 1, last column), which is larger than the energies themselves.
2. Under (b) the covariance is exactly $(N-1)^2=602^2\approx3.6\times10^5$ times the covariance under (a). Table 1's ratio $\sigma_{\rm fitter}/\sigma_{\rm JK}=602.0$ holds for all 19 files.
3. The saved fitter covariance `Ecm_cut_0p345/l20_full_cut_0p345_covariance_allL.dat` has diagonal $\approx0.096$, i.e. $\sigma\approx0.31$. This is the (b) value.
4. **χ² recomputation.** Using the saved model energies, `chi2_check.py` reproduces the fitter-reported χ² to all printed digits under (b). Under (a) it gives values about $3.6\times10^5$ times larger:

| cutoff | levels | χ² reported (b) | χ² recomputed (b) | χ² with (a) | dof | χ²/dof with (a) |
|---|---|---|---|---|---|---|
| 0.315 | 6 | 3.262e-5 | 3.2621e-5 | 11.82 | 2 | 5.91 |
| 0.325 | 9 | 9.531e-5 | 9.5309e-5 | 34.54 | 5 | 6.91 |
| 0.335 | 12 | 1.283e-4 | 1.2833e-4 | 46.51 | 8 | 5.81 |
| 0.345 | 14 | 1.671e-4 | 1.6711e-4 | 60.56 | 10 | 6.06 |

   Command: `python3 -I scratchpad/scripts/chi2_check.py $PKG`. Under (a), the correlation-matrix condition number is $6.5\times10^2$ to $1.9\times10^3$, and the largest single-level pull is 1.1–1.7σ.
5. **The fits did not really minimise.** Because the χ² surface is about $3.6\times10^5$ times too flat, Minuit's stopping test (expected distance to minimum $<10^{-3}\times$ tolerance) is met almost at once. This is consistent with the fit table:
   - $K^E$ is unchanged from its start value (−1226756.7068845264) at every cutoff.
   - $K^B$ is unchanged at three of four cutoffs.
   - The fitted parameter errors are $10^8$–$10^9$, or exactly 2.00 at cutoff .325.
   
   The "χ²/dof ≈ 1.6e-5" in `FULL_L20_FIT_AND_VALIDATION_REPORT.md` is therefore an artefact of the error scale, not a statement about fit quality. The old isotropic fit's χ²/dof ≈ 1.5 used reader (a), which explains the "not apples-to-apples" remark in that report.
6. **Why the internal audit missed it.** `PKG/scripts/audit_l20_chi2_reconstruction.py` reads the fitter's own saved covariance, so it could not detect the scale error.

**Recommendation for the new repo.**
- Use one reader with an explicit config key, e.g. `jack_sample_convention = ensemble | jackknife`, defaulting to `ensemble`, i.e. resample first as in (a).
- Add a unit test (call it an acceptance check) that `20_000_A1m_n0.jack` gives $a_tE_{\rm cm}=0.269227$ with $\sigma=5.15\times10^{-4}$ (tolerance $10^{-6}$ relative).
- Every L20 fit result so far must be redone after this fix.
- The tester or planner should confirm the file convention with whoever produced the `.jack` files, because that is the one remaining assumption.

### 3.3 User true-zero label CSVs and accepted windows

**Label CSVs** live in `PKG/input/`. Columns:

```
bracket_id,E_left,E_right,user_label,user_zero_estimate,comment
```

- Each row is one sign-change bracket of the scaled projected determinant $\det\!\big[(\xi L/a_s)^{-6}\,V_{\rm sel}^\dagger(F_3^{-1}+K_{\rm df,3})V_{\rm sel}\big]$, evaluated at the baseline $K_{\rm df,3}$ (§3.6).
- `E_left`/`E_right` are $a_tE_{\rm cm}$ grid points. They are consecutive nodes of `linspace(E0, 0.36, 20000)`, so the bracket width is $(0.36-E_0)/19999\approx4.845\times10^{-6}$. I checked that the bracket edges sit at integer grid indices.
- The tag `E026310_0360` means $E_0=0.26310$. `E026301` means $E_0=0.26301$.
- `user_label` is `true` (a real QC zero, i.e. a model energy level), `false` (a pole or spurious sign flip), or blank.
- Only `true` rows are used. The other two columns are always empty.

| file (all in `PKG/input/`) | bytes | rows | true | role |
|---|---|---|---|---|
| `v33p_user_truezero_labels_L20_000_A1m_E026310_0360.csv` | 653 | 10 | 3 (0.26932, 0.34538, 0.35795) | **used** |
| `v33p_user_truezero_labels_L20_100_A2_E026301_0360.csv` | 1554 | 25 | 2 (0.29156, 0.30044) | **used** (grid $E_0=0.26301$, remapped to the 0.26310 grid when windows are built) |
| `v33p_gpu_user_truezero_labels_L20_110_A2.csv` | 2212 | 36 | 4 | **used**; identical (md5 `e80cd5b8…`) to `..._L20_110_A2_E026310_0360.csv` |
| `v33p_gpu_user_truezero_labels_L20_111_A2.csv` | 1851 | 30 | 5 | **used**; identical (md5 `15a08872…`) to `..._L20_111_A2_E026310_0360.csv` |
| `v33p_gpu_user_truezero_labels_L20_200_A2.csv` | 2151 | 35 | 5 (0.28624, 0.29918, 0.33098, 0.34797, 0.35346) | **used** |
| `v33p_user_truezero_labels_L20_100_A2.csv` | 384 | 8 | 0 | obsolete partial file (only $E>0.34$); drop |
| `v33p_user_truezero_labels_L24_000_A1m_E026310_0360.csv` | 1253 | 20 | 3 | L=24, not yet used |
| `v33p_gpu_user_truezero_labels_L24_100_A2.csv` | 3167 | 39 | 7 | L=24; identical (md5 `f79dafed…`) to `..._L24_100_A2_E026310_0360.csv`; "bracket 36" unresolved |

In total the five L20 files have 3+2+4+5+5 = **19 accepted windows** through 0.36.

**Mismatches between windows and lattice levels.** These are not errors in the files, but they block fits:
- **L20/200_A2.** The window at 0.35346 has no lattice partner (only four `002` files). At cutoff 0.355 there are 5 windows but 4 lattice levels. This is the documented blocker.
- **L20/000_A1m.** The window at 0.35795 has no lattice partner (two `000` files). This will block cutoff 0.36 the same way.

**Accepted windows** live in `PKG/output/accepted_windows_L20_user_validated_0p36/L20_<irrep>/accepted_windows.csv` (7 661 bytes in total). Columns:

```
Lbyas, irrep, lattice_level_index, lattice_Ecm, bracket_id, E_left_bracket, E_right_bracket,
zero_estimate_initial, center_row, row_left, row_right, max_row_left, max_row_right,
previous_model_zero, has_previous_model_zero, inside_Ecm_cutoff, source_label_file, source_label_row
```

- Each row pairs a true window with a lattice level in sorted order.
- `center_row` is the grid index; `row_left/right` = centre ± 50 (`window_half_rows`); `max_row_left/right` = centre ± 250 (`max_window_half_rows`).
- The projected-basis "hot window" caches in `PKG/output/hot_windows_L20_user_validated_0p36/` (1.97 MB of JSON and CSV) store the row lists, provenance hashes and validation results. Rows saved: 1507, 1006, 2012, 2511 and 2509.
- **The `.bin` payloads are not in `extracted/`.** `find -name '*.bin'` finds 0 files. They must come from the "full zip" or be rebuilt from the raw cache at `/media/digonto/Data/F3inv_cache`.

### 3.4 Config files and keys

The parsers are `settings_from_config` (`PKG/source/qc_fitter_norm_refine_v2.cpp` line 82) and `multiconfig_from_config` (`..._multiL.cpp` lines 182–257). The format is `key = value`, with `#` starting a comment. Lists are space- or comma-separated. Defaults are hard-coded in those two functions.

**Production configs (all byte-identical between the v33n, v33o, v33p and L20 packages; `diff -rq` on `configs/` is empty):**

| file | bytes | role |
|---|---|---|
| `PKG/configs/config_v33m_multiL_all_irreps_fastest_production.in` | 2355 | locked v33m production: L20+L24, cutoff 0.335, `cache_backend = v33g_runtime`, `classifier_mode = raw_sign_only` |
| `PKG/output/l20_full_cutoff_fits_user_validated_0p36_passonly/fitter_input_base_user_validated.in` | 5326 | **latest working config** (L20 hot-window, `root_search_mode = accepted_windows`) |
| `PKG/output/.../Ecm_cut_0p3xx/fitter_input.txt` | ~5 kB each | the exact per-cutoff inputs used |
| `PKG/configs/v33n_cutoff_scan/config_v33n_Ecmcut_0p3{15,25,35,45,55,65}.in` | 13 290 total | cutoff scan; same keys as v33m |
| `PKG/configs/config_v32w_QC_fitter_norm_refine_v2_fixedcoarse.in` | 2478 | single-L v32w; uses the legacy szscl21 reader |
| `full_cache_generator_package_20260629-main/.../cpu_component_cachegen/configs/config_cpu_component_cachegen.in` | 786 | cache generator; same physics keys plus `Q0norm = 1`, `sort_orbit_flag = 0` |

**Physics keys.** All are in temporal lattice units. Each value is identical in every production config, and `PKG/scripts/check_v33p_production_scatter_params.py` checks them all.

| key | value | meaning (from code) |
|---|---|---|
| `xival` (alias `xi`) | 3.444 | renormalised anisotropy $\xi=a_s/a_t$. Momenta are $a_t\vec p=\frac{2\pi}{\xi L/a_s}\vec n$ |
| `atmK` | 0.09698 | $a_tm_K$; "M1" (two identical particles) |
| `atmpi` | 0.06906 | $a_tm_\pi$; "M2" |
| `eta_1` | 1.0 | symmetry factor $\eta_1$ in $K_2^{-1}$ for flavour 1: kaon spectator, $K\pi$ pair (distinguishable) |
| `eta_2` | 0.5 | $\eta_2$ for flavour 2: pion spectator, $KK$ pair (identical particles) |
| `alpha` | 0.5 | $\alpha$ of the smooth exponential UV cutoff in the $F_2$ sum–integral difference ("alpha_KSS" in QC3_release `constants.py`) |
| `epsilon_h` | 0.0 | $\epsilon_H$ shift in the argument of the spectator cutoff function $H$ (`cutoff_function_1`) |
| `max_shell_num` | 20 | integer range $\lvert n_{x,y,z}\rvert\le20$ for the spectator-momentum search and the $F_2$ lattice sum |
| `tolerance` | 1e-12 | values of $H$ below this are set to 0 when building configurations |
| `parity` | −1 | parity factor used by the irrep projector (three pseudoscalars) |
| `eig_tol` | 0.05 | keep projector eigenvectors with eigenvalue in $[1-0.05,\,1+0.05]$ to form $V_{\rm sel}$ |
| `norm_tol` / `proj_tol` | 1e-12 / 1e-10 | normalisation and projection checks on $V_{\rm sel}$ |
| `waves_vec_1` | `0 1` | partial waves in the $K\pi$ pair: s and p |
| `waves_vec_2` | `0` | partial waves in the $KK$ pair: s only |
| `scatter1_00` | 4.04 | $K\pi$ s-wave $a_0$ (units of $a_t$) |
| `scatter1_10` | −43.2 | $K\pi$ p-wave $a_1$ (units of $a_t^3$; see the caveat in §3.5) |
| `scatter2_00` | 4.12 | $KK$ s-wave $a_0$ (units of $a_t$) |
| `scatterI_ab` (other entries) | 0 | key format is `scatter<channel>_<ell><k>`: $\ell=0..3$; $k=0$ is $a_\ell$, $k=1$ is $r_\ell$, $k=2$ is $P_\ell$ |

**Data, grid and fitting keys** in the multi-L fitter:

| key | production value | meaning |
|---|---|---|
| `Lbyas_values` | `20 24` (v33m) / `20` (L20) | box sizes to fit |
| `irreps_L20`, `irreps_L24`, `list_of_mom` | `000_A1m 100_A2 110_A2 111_A2 200_A2` | irreps per L |
| `lattice_jackknife_dir` | `/home/digonto/Codes/KKpi_I2/spectrum/Ecm_data/data` | location of the `.jack` files; **must be repointed** |
| `lattice_jack_energy_type` | `En_lab` | file holds $E_{\rm lab}$, converted per sample |
| `jack_skip_header_lines`, `jack_energy_column` | 1, 2 | header lines to skip; 1-based column to read |
| `Ecm_cutoff` (= `energy_cutoff`) | 0.335 / 0.315–0.345 | keep levels with mean $E_{\rm cm}\le$ cutoff |
| `max_total_lattice_levels`, `max_lattice_levels_per_block` | 0 | truncation for smoke tests (0 = off) |
| `scan_E0`, `scan_E1`, `coarseN` | 0.26301 or 0.26310, 0.36, 20000 | coarse $E_{\rm cm}$ grid of the precomputed $F_3^{-1}$/$V_{\rm sel}$ cache |
| `gpu_cachegen_*` | same as above | used only when `auto_build_missing_coarse=1` |
| `cache_backend` | `v33g_runtime` (v33m) / `raw_v33e` (L20) | cache-reader choice |
| `coarse_cache_root` | `/media/digonto/Data/F3inv_cache` | external raw cache (multi-GB; for example the 200_A2 cache is 3.26 GB) |
| `projected_basis_cache_mode`, `projected_basis_cache_root` | `hot_windows`, path | use the hot-window projected caches |
| `root_search_mode` | `accepted_windows` | search for model zeros only inside user-true windows |
| `accepted_zeros_file_L<L>_<irrep>`, `candidate_brackets_file_…`, `det_grid_file_…` | paths | per-sector window and diagnostic inputs |
| `window_half_rows`, `max_window_half_rows` | 50, 250 | grid rows searched around each window centre |
| `fallback_full_scan` | 0 | 1 would allow a slow full-cache scan inside the minimisation function (forbidden in production) |
| `classifier_mode` | `digonto_v3_window` (L20) / `raw_sign_only` (v33m) | zero-versus-pole classifier, Digonto's own heuristic |
| `v3_*`, `v2_*`, `v4_merge_tol` | see file | classifier tunables |
| `const_norm_power` | 6 | QC matrix divided by $(\xi L/a_s)^{6}$ before the determinant |
| `det_backend` | `auto` / `cpu_openmp` | determinant engine |
| `K3iso0_guess`…`K3E_guess`, `*_step` | see §3.6 | Minuit start values; steps $10^3$/$10^4$ |
| `float_K3iso0`…`float_K3E` | 1 | 0 fixes the parameter |
| `use_parameter_limits`, `param_lower/upper` | 0, ±1e8 | parameter bounds |
| `use_lattice_covariance`, `chi_square_mode` | 1, `raw_cov_inv` | $\chi^2=r^TC^{+}r$ with $C^+$ a pseudoinverse (eigenvalues below $10^{-12}\max$ dropped); alternatives are `diagonal` and `corr_inv` |
| `failure_penalty` | 1e100 | χ² returned when a level is missing |
| `max_fcn`, `max_fcn_evals` | 0 / 400 / 500 | call budgets |
| `refined_cache_*`, `save_refined_cache` | disabled | legacy refinement caches |
| `output_dir`, `output_tag`, `print_*`, `omp_threads`, `debug` | — | input/output and threading |

**v32w-only keys** (legacy szscl21 reader):
- `Lval`, `ensemble`, `masses_path`, `max_state = 8` (state indices 0..7 tried), `lattice_energy_type`, `print_found_files`.
- `load/save/require_existing_binary_f3inv_cache`, `binary_f3inv_cache_file`, `write_cache_grid`.

### 3.5 Two-body input ($K_2^{-1}$) and where the scattering parameters come from

The code is `K2_inv_ERE_ang_mom` (`K2_functions_v2.h`) and `K2inv_ell_gpu` (`K2_functions_gpu_safe.cuh`, the production cache generator). With $q^{*2}$ the pair relative momentum squared and $H$ the spectator cutoff function:

$$\big[K_2^{-1}\big]_{\ell=0}=\frac{\eta_i}{8\pi\sqrt{\sigma_i}}\Big[-\frac{1}{a_0}+\tfrac12 r_0q^{*2}+P_0q^{*4}+|q^*|\,(1-H)\Big],$$

$$\big[K_2^{-1}\big]_{\ell=1}=\frac{\eta_i}{8\pi\sqrt{\sigma_i}}\Big[-\frac{1}{a_1q^{*2}}+\tfrac12 r_1+P_1q^{*2}+q^{*2}|q^*|\,(1-H)\Big].$$

The full matrix is block-diagonal, $\mathrm{diag}\big(K_{2,1}^{-1},\,2K_{2,2}^{-1}\big)$, with $(\eta_1,\eta_2)=(1,\tfrac12)$.

**Sign convention.** $q\cot\delta_0=-1/a_0$, so $a_0>0$ means $\delta_0<0$ (repulsive). Both s-wave values are positive.

**Caveat (outside my area, flagged for the auditor).** The two terms inside the $\ell=1$ bracket do not have matching powers of $q^*$:
- The scattering-volume term is $q\cot\delta_1=-1/(a_1q^{*2})$.
- The threshold term is $q^{*2}|q^*|(1-H)$, i.e. the form that goes with $q^3\cot\delta_1$.
- QC3_release (`base_code/F3/K2i_mat.py`) uses $q^3\cot\delta_1+q^{*2}|q^*|(1-H)$ with `qcot_fit_p = -1/a1`.

So the C++ p-wave $a_1$ term differs from QC3_release by a factor $1/q^{*2}$ in the scattering-volume term. This may be intentional if the C++ harmonics carry no $|q^*|^\ell$ factors, but then the threshold term would be inconsistent. It needs a direct numerical check of C++ against QC3_release for $K_2^{-1}$ at one $(E,\vec k)$.

**Origin of 4.04, −43.2 and 4.12.** I found no fit, script or citation that produces them. They appear only as hard-coded defaults:
- `scatter_params_1[4][3] = {{4.04,0,0},{-43.2,0,0},…}` in `K3df_minuit_fit_v32f_…hpp` line 151, and in `F3_cpu_openmp_v25_K3QC_cached_core.hpp` lines 133–136.
- In the docs (`V33E_CACHE_AND_FITTER_LOGIC.md`) they are named `piK_swave_a0`, `piK_pwave_a1`, `KK_swave_a0`.

That doc also plans a future "two-body + K3df" mode with these three as floating parameters. No two-body fit data (πK or KK jackknife files) is present in the tree. **They should be treated as externally supplied numbers whose source and uncertainties are unknown, and Digonto should be asked.**

### 3.6 $K_{\rm df,3}$ parametrization (exact, from `PKG/source/K3_functions_2plus1.hpp`)

This file is a line-by-line C++ port of QC3_release `base_code/Kdf3/K3main.py`, `K3B.py` and `K3E.py` (Blanton, Romero-López, Sharpe, arXiv:2111.12734). It is identical (md5 `276f568e…`) in all seven post-v32w packages.

Define $M\equiv2m_K+m_\pi$ and

$$\Delta=\frac{E_{\rm cm}^2}{M^2}-1 .$$

For a matrix element between spectator flavour $i$ (momentum $\vec p$, pair wave $\ell'm'$) and flavour $j$ ($\vec k$, $\ell m$):

$$\big[K_{\rm df,3}\big]_{i\vec p\ell'm';\,j\vec k\ell m}=c_{ij}\Big\{\delta_{\ell'0}\delta_{\ell0}\big(K_{3,\rm iso}^{(0)}+K_{3,\rm iso}^{(1)}\Delta\big)+K^B\,\mathcal B_{ij}+K^E\,\mathcal E_{ij}\Big\},$$

with block factors $c_{11}=1$, $c_{12}=c_{21}=1/\sqrt2$, $c_{22}=1/2$.

The two linear-order kinematic structures are `K3B_element` ($\mathcal B_{ij}$) and `K3E_element` ($\mathcal E_{ij}$). They have s–s, s–p, p–s and (for $\mathcal E_{11}$) p–p pieces built from boosted momenta, and all are divided by $M^2$. Two examples:

$$\mathcal B_{22}^{ss}=\frac{2}{M^2}\Big[E_{\rm cm}^2-E(\omega_p+\omega_k)+\vec P\cdot(\vec p+\vec k)+m_\pi^2-4m_K^2\Big],\qquad \mathcal E_{22}^{ss}=\frac{2}{M^2}\Big[m_\pi^2-\omega_p\omega_k+\vec p\cdot\vec k\Big].$$

All other entries are in `K3_functions_2plus1.hpp` lines 140–433.

In the QC, $K_{\rm df,3}$ enters as

$$\det\!\Big[\tfrac{1}{(\xi L/a_s)^6}\,V_{\rm sel}^\dagger\big(F_3^{-1}+K_{\rm df,3}\big)V_{\rm sel}\Big]=0 .$$

There is no extra $1/L^3$ or $1/(2\omega)$ factor on $K_{\rm df,3}$ in the C++, so the parameter values are tied to the code's $F_3$ normalisation, in $a_t$ units. Because $K_{\rm df,3}$ is linear in its parameters, the hot-window caches store four projected basis matrices $B_0,B_1,B_B,B_E$ (`v33g_runtime_k3basis_cache.hpp`).

**Parameter values** (dimensionless code units; nothing in the code attaches physical units):

| set | $K_{3,\rm iso}^{(0)}$ | $K_{3,\rm iso}^{(1)}$ | $K^B$ | $K^E$ | source |
|---|---|---|---|---|---|
| "baseline" (label generation) | 73735.840894011912 | −972421.14060757787 | 347174.05548116949 | −1226756.7068845264 | `PKG/README.md` |
| v33m production start | 34470.228772793678 | −1035173.9140440958 | 253921.14757953485 | −961680.91711685632 | `config_v33m_…production.in` |
| L20 fit, cut .315 | 151422.46 | −972421.14 | 347174.06 | −1226756.71 | `L20_PASSONLY_FIT_SUMMARY.csv` |
| L20 fit, cut .325 | 125968.22 | −972421.14 | 347174.06 | −1226756.71 | same |
| L20 fit, cut .335 | 32040.84 | −1074716.65 | 347174.06 | −1226756.71 | same |
| L20 fit, cut .345 | 40607.52 | −1067567.87 | 376187.14 | −1226756.71 | same |

The fit rows are not trustworthy minima (§3.2). They are only start values for a refit.

The legacy isotropic approach (`generate_K3iso.cpp`, `K3iso_maker.py`) instead computed $K_{3,\rm iso}=-1/F_{3,\rm iso}$ at each lattice energy and fitted $K^{(0)}+K^{(1)}\Delta$. Those files are history only.

### 3.7 Non-interacting three-body energies

All of these use

$$a_tE_{\rm lab}=\sum_{a=1}^{3}\sqrt{m_a^2+q^2|\vec n_a|^2},\quad q=\frac{2\pi}{\xi L/a_s},\quad \vec n_1+\vec n_2+\vec n_3=\vec n_P,\quad E_{\rm cm}=\sqrt{E_{\rm lab}^2-q^2|\vec n_P|^2}.$$

| package (`NEST/…`) | language | what it adds | status |
|---|---|---|---|
| `nonint3body_energies_package` | C++ (12 kB) | raw energy list plus a raw count of momentum triplets; no irreps | superseded |
| `nonint3body_irrep_energies_package` | C++ (20 kB) | permutation-representation irrep multiplicities (treats particles as scalars, so wrong $A_1/A_2$ assignment for pseudoscalars) | superseded |
| `nonint3body_irrep_pseudoscalar_package` | C++ (22 kB) | character $\det(g)^{N_{ps}}\times$(number of fixed points), $N_{ps}=3$: the correct pseudoscalar irreps (e.g. moves levels from `100_A1` to `100_A2`) | correct; C++ alternative |
| `nonint3body_general_jpc_package` | C++ (17 kB) | general $J^{PC}$ per particle from config | intermediate |
| `nonint3body_general_doublecover_package` | **Python** (23 kB) | general $J^{PC}$, double cover, spinorial irreps; one file per (P, irrep) with columns `level Ecm group_degeneracy multiplicity irrep_dim raw_energy_multiplicity example_n1..3 …` | **newest; the one the plotters call** |

The double-cover script has a single version everywhere (md5 `8a84183b…`, about 60 copies). Use `PKG/tools/nonint3body_general_doublecover_package/` (26 925 bytes, including configs and README).

The L20 plotter writes its own config for each $L/a_s=16..24$ and calls this script, so `source/nonint3body_general_doublecover.py` plus NumPy is all it needs. I checked a sample output: at $L/a_s=20$, `P001_A2`, the lowest level is $E_{\rm cm}=0.284934$ (one kaon with $\vec n=(0,0,1)$).

### 3.8 Plotting and analysis scripts worth keeping

| path | bytes | what it does | external inputs it still needs |
|---|---|---|---|
| `PKG/L20_spectrum_plotter_model_true_zeros_yaxis_fixed.py` | 23 160 | 5-panel L20 spectrum: old QC (x=18), lattice (x=20), model true zeros (x=22), non-interacting curves for L=16–24, thresholds $KK\pi$=0.26302 and $KK\pi\pi$=0.33208 | **`V20_<irrep>` lattice files and old `QC_states_jackknifed_*` files (both missing)**; paths hard-coded to `/home/digonto/.../3body_quantization_v2/`. Should be rewritten to read the `.jack` files directly. |
| `PKG/scripts/plot_v33p_classifier_det_family_interactive.py` | 3 900 | interactive determinant-family plot (curves of $\det\times10^k$ with candidates and markers) for labelling true zeros | det-grid CSVs (not packaged) |
| `PKG/scripts/build_l20_user_validated_windows.py` | 10 844 | labels to accepted windows (true rows only, sha256 provenance); alias map `001→100` etc. | candidate-bracket CSVs (not packaged) |
| `PKG/scripts/score_v33p_classifier_against_user_labels.py` | 4 089 | scores an automatic classifier against the user labels | — |
| `PKG/scripts/check_v33p_production_scatter_params.py` | 6 472 | fails if a config deviates from the production physics constants | — (good as a CI check, i.e. an automatic test run on every change) |
| `PKG/scripts/audit_l20_chi2_reconstruction.py` | 4 618 | χ² re-contraction | reuses the fitter covariance, so blind to §3.2; fix before reuse |
| `REPRO/tools/make_qc_state_lines_from_scattered_spectrum.py` | 21 103 | tracks QC levels across $L/a_s$ (assignment plus PCHIP/spline interpolation) into continuous state lines | SciPy; v32x/v32y spectrum-vs-L outputs |
| `REPRO/tools/inspect_v32_f3inv_cache_v32f_fixed.py` | 11 920 | main inspector for `V32F_F3INV_FULL_CACHE_1` binary caches (entries per irrep, grid completeness) | raw caches |
| `REPRO/tools/check_v32x_L24_cache_ready.py` | 7 832 | cache readiness check | raw caches |
| `PKG/plots/` | 821 061 | reference PDF/PNG (`KKpi_L20_model_spectrum_fitcut_0p335.pdf`, classifier det family) | — |

The `REPRO/tools/` files `inspect_v32_f3inv_cache{,_robust,_mmap_progress}.py` and the older `L20_old_lattice_recent_*` plotters are earlier variants and can be dropped.

### 3.9 `QC3_release` and `source` in `_nested_old_versions`

**`NEST/QC3_release/QC3_release`** (7.6 MB working tree; `base_code/` 0.55 MB):
- This **is** the public Python code of T. Blanton, F. Romero-López and S. Sharpe ("Release of the implementation of the QC3 … Accompanying code to arXiv:2111.12734").
- It is a git clone of `https://github.com/ferolo2/QC3_release.git` at commit `9bc4b84` ("Update F_mat.py to fix error in product").
- It has local, uncommitted edits by Digonto (+248/−8 lines across `F3_mat.py`, `F_mat.py`, `G_mat.py`, `K3main.py`, `defns.py`, `group_theory_defns.py`, `projections.py`, `solution_2pt.py`, `test.py`) and debug files, including a 6.4 MB `F2_python_debug_log.txt`.

Relevance:
1. It is the reference ("oracle") for $F_2$, $G$, $K_2^{-1}$ and $K_{\rm df,3}$ of 2+1 systems. The C++ $K_{\rm df,3}$ is a direct port of it.
2. `NEST/source/source/compare_K3_cpp_vs_python.py` with `test_K3mat_2plus1_cpp_dump.cpp` and the dump `k3_cpp_matrix.dat`/`k3_cpp_configs.dat` is a ready-made regression test for $K_{\rm df,3}$ against QC3. It uses $K^{(0)}=200$, $K^{(1)}=400$, $K^B=400$, $K^E=300$, waves s+p / s.
3. It is the natural tool for the p-wave $K_2^{-1}$ check in §3.5.

There is no LICENSE file in the snapshot, so check upstream before redistributing. My recommendation: vendor pristine upstream `base_code/` at `9bc4b84` under `third_party/` for tests only, and keep Digonto's diffs as a patch file.

**`NEST/source/source`** (119 files, 4.5 MB): a C++ snapshot (Feb–May 2026) of the GPU/CPU pipeline v3–v19 era. Compared with `PKG/source`:
- 3 files differ (`lattice_data_covariance_cpp.hpp` is older, without the audit rows; `functions.h`; `F3_cpu_openmp_v25_K3QC_cached_core.hpp`).
- 5 files are only here: two compiled test binaries, `build.log`, and the K3 dump `.dat` files.

Only the K3 regression items above are worth keeping; everything else is superseded.

---

## 4. Recommended files to copy into the new repo

Sizes are in bytes (`du -sb`). "Copy" means verbatim. "Rewrite" means take it as the specification and re-implement.

**Data (all needed for L20 fits):**

| # | source path | bytes | target |
|---|---|---|---|
| 1 | `PKG/input/lattice_spectrum/` (19 × `20_*.jack`) | 136 407 | `data/lattice/L20/` |
| 2 | `PKG/input/v33p_user_truezero_labels_L20_000_A1m_E026310_0360.csv` | 653 | `data/labels/` |
| 3 | `PKG/input/v33p_user_truezero_labels_L20_100_A2_E026301_0360.csv` | 1 554 | `data/labels/` |
| 4 | `PKG/input/v33p_gpu_user_truezero_labels_L20_110_A2.csv` | 2 212 | `data/labels/` (rename without `gpu_`) |
| 5 | `PKG/input/v33p_gpu_user_truezero_labels_L20_111_A2.csv` | 1 851 | `data/labels/` |
| 6 | `PKG/input/v33p_gpu_user_truezero_labels_L20_200_A2.csv` | 2 151 | `data/labels/` |
| 7 | `PKG/input/v33p_user_truezero_labels_L24_000_A1m_E026310_0360.csv` | 1 253 | `data/labels/` (L24, provisional) |
| 8 | `PKG/input/v33p_user_truezero_labels_L24_100_A2_E026310_0360.csv` | 3 167 | `data/labels/` (L24, bracket 36 unresolved) |
| 9 | `PKG/output/accepted_windows_L20_user_validated_0p36/` | 7 661 | `data/windows/L20/` |
| 10 | `PKG/output/hot_windows_L20_user_validated_0p36/` (JSON+CSV only; `.bin` absent) | 1 973 236 | `data/hot_windows/L20/` (metadata; rebuild the binaries) |

Skip: `v33p_user_truezero_labels_L20_100_A2.csv` (obsolete), and the `_E026310_0360` duplicates of items 4, 5 and 8.

**Missing, must be requested from Digonto:**
- `24_*.jack` (L=24 spectrum).
- `V20_<irrep>` and old `QC_states_jackknifed_*` (plotter only).
- The hot-window `.bin` payloads (or the "full zip").
- The raw $F_3^{-1}/V_{\rm sel}$ cache root (`/media/digonto/Data/F3inv_cache`).
- The provenance of the two-body parameters.

**Configs:**

| # | source path | bytes | note |
|---|---|---|---|
| 11 | `PKG/output/l20_full_cutoff_fits_user_validated_0p36_passonly/fitter_input_base_user_validated.in` | 5 326 | template for the new config; strip absolute paths |
| 12 | `PKG/configs/config_v33m_multiL_all_irreps_fastest_production.in` | 2 355 | L20+L24 reference |
| 13 | `PKG/configs/config_v32w_QC_fitter_norm_refine_v2_fixedcoarse.in` | 2 478 | legacy single-L reference (documents the szscl21 reader keys) |
| 14 | `full_cache_generator_package_20260629-main/full_cache_generator_package_20260629-main/cpu_component_cachegen/configs/config_cpu_component_cachegen.in` | 786 | cache-generation physics keys (`Q0norm`, `sort_orbit_flag`) |

**Code (specification or tests):**

| # | source path | bytes | action |
|---|---|---|---|
| 15 | `PKG/source/lattice_data_covariance_cpp.hpp` | 18 449 | copy; its resample-then-covariance path is the correct convention |
| 16 | `PKG/source/K3_functions_2plus1.hpp` | 16 452 | copy (the $K_{\rm df,3}$ definition) |
| 17 | `PKG/source/K2_functions_v2.h` | 10 847 | copy, after resolving the §3.5 p-wave question |
| 18 | `NEST/QC3_release/QC3_release/base_code/` | 550 612 | `third_party/QC3_release` (pristine upstream `9bc4b84`); test oracle |
| 19 | `NEST/source/source/compare_K3_cpp_vs_python.py` | 5 837 | `tests/` |
| 20 | `NEST/source/source/test_K3mat_2plus1_cpp_dump.cpp` | 4 766 | `tests/` |
| 21 | `NEST/source/source/k3_cpp_matrix.dat`, `k3_cpp_configs.dat` | 94 911 + 2 555 | `tests/data/` (reference dump) |

**Tools and scripts:**

| # | source path | bytes |
|---|---|---|
| 22 | `PKG/tools/nonint3body_general_doublecover_package/` | 26 925 |
| 23 | `PKG/L20_spectrum_plotter_model_true_zeros_yaxis_fixed.py` (rewrite data readers) | 23 160 |
| 24 | `PKG/scripts/plot_v33p_classifier_det_family_interactive.py` | 3 900 |
| 25 | `PKG/scripts/build_l20_user_validated_windows.py` | 10 844 |
| 26 | `PKG/scripts/score_v33p_classifier_against_user_labels.py` | 4 089 |
| 27 | `PKG/scripts/check_v33p_production_scatter_params.py` | 6 472 |
| 28 | `REPRO/tools/make_qc_state_lines_from_scattered_spectrum.py` | 21 103 |
| 29 | `REPRO/tools/inspect_v32_f3inv_cache_v32f_fixed.py` | 11 920 |
| 30 | `NEST/nonint3body_irrep_pseudoscalar_package/` (optional C++ cross-check of #22) | 176 728 |

**Previous results, kept as reference values for regression tests only:**
- `PKG/output/l20_full_cutoff_fits_user_validated_0p36_passonly/` (1 524 965 bytes): especially `L20_PASSONLY_TRUE_ZEROS_ALL_WINDOWS_ALL_CUTOFFS.csv` (27 164) and `Ecm_cut_*/lattice_vs_model_levels.csv`.
- `PKG/output/l20_k3b_k3e_zero_tests/` (7 342 313 bytes) is optional.

Total for items 1–29 without the optional ones: about 3.0 MB, of which hot-window metadata is 1.97 MB.

---

## 5. Problems and assumptions

1. **Error-scale bug (§3.2).** I assume the `.jack` files are per-configuration ensembles. The evidence is strong: the factor of exactly 602, errors larger than the energies under the alternative, the legacy reader's choice, and the stuck parameters. It still needs confirmation from the file's producer. If confirmed, all v33f–L20 fit χ² values and parameter errors are invalid by factors $602^2$ and $602$ respectively, and the "fits" are near their start points.
2. **System naming.** The brief says "K K-bar π". The code treats the two kaons as **identical** ($\eta_2=\tfrac12$, the $KK$ pair has s-wave only, the 2+1 formalism of identical-pair plus one), and the data path is `KKpi_I2`. That reads as maximal isospin $I=2$, $K^+K^+\pi^+$. Please confirm; a $K\bar K\pi$ system would need a different formalism (coupled channels, no identical-pair factor).
3. **Two-body parameters** (4.04, −43.2, 4.12): their origin, uncertainties and correlations are not documented anywhere in the tree. The p-wave $K_2^{-1}$ convention may also be inconsistent (§3.5). Both need the auditor or Digonto.
4. **Missing L=24 spectrum** and missing hot-window binaries. The L20+L24 production target cannot be reproduced from `extracted/` alone.
5. **Window/level mismatches:** 200_A2 (5 vs 4 at 0.355) and 000_A1m (3 vs 2 at 0.36). The fitter pairs them by sorted order, so these sectors block higher cutoffs until a matching rule is decided (planner's call).
6. **Grid mismatch:** the 100_A2 labels use $E_0=0.26301$ while the hot-window grid uses $E_0=0.26310$. `build_l20_user_validated_windows.py` remaps through the candidate-bracket CSV (bracket 7: 0.2915556 → 0.2915513). The new repo should use one grid.
7. **Mass and anisotropy uncertainties** are not propagated anywhere: $a_tm_\pi$, $a_tm_K$ and $\xi$ are fixed constants. The uncertainty of $\xi=3.444$ is not in the tree.
8. **Ensemble identification** ($m_\pi\approx391$ MeV HadSpec) is from my memory of published values, not from a file in the tree.
9. `hashes.txt` does not include `_nested_old_versions`. Duplicate checks there used fresh `md5sum` runs.

# F — Development history v31 → v32y: lessons, bugs fixed, dead ends

> **Research report — consolidation investigation.** Date: 2026-10-21. Git commit: none (written before the repo existed; first commit of `3body-qc-v4`). Prompt in force: `prompts/000-consolidation-and-roadmap.md`. Related: `reports/research/001`–`006`.
> Absolute paths below refer to the unpacked original zips (`/home/digimonk/Projects/three-body-qc-codes/extracted/...`) and to a temporary scratchpad. For each file's new location in this repo, see `docs/PROVENANCE.tsv`.
> The system is two kaons + one pion treated as a 2+1 system (identical kaons). Where this report says "K K-bar pi", read "KKπ"; that wording came from the investigation brief and is not established. See `docs/OPEN_ISSUES.md` #7.


- **Date:** 2026-10-21
- **Git commit:** none (the `extracted/` tree is not a git repository; this is a read-only archive review)
- **Prompt in force:** coordinator brief "Area F: development history v31 → v32y"
- **Related reports:** the other consolidation-area reports (A–E) written in the same scratchpad `reports/` folder
- **Mode:** read only. Nothing was compiled or run. Every statement below comes from a README, a patch-notes file, a config file or a source line, and the file is named. Where I give my own interpretation, it is labelled **[my inference]**.

---

## 0. Glossary (read first)

| Term | Meaning in this code base |
|---|---|
| QC | Quantization condition: the equation whose solutions $E$ are the finite-volume energy levels. Here (relativistic field theory, "RFT", formalism) it is $\det\!\big[V^\dagger\,(F_3^{-1}+\mathcal K_{\mathrm{df},3})\,V\big]=0$. |
| $F_2$, $G$ | Finite-volume kinematic matrices. $F_2$ is the two-particle loop sum and $G$ is the one-particle-exchange ("switch") matrix. Both have poles at non-interacting energies. |
| $K_2$ | Two-particle K-matrix, built from effective-range-expansion (ERE) parameters `scatter{1,2}_{ell}{index}`. |
| $H$, $F_3$ | $H=K_2^{-1}+F_2+G$ and $F_3=\tfrac13F_2-F_2H^{-1}F_2$ (as written in `README_v31zl_SCATTER_PARAM_FIX.md`). |
| $\mathcal K_{\mathrm{df},3}$ (K3df) | Short-distance three-particle K-matrix; the fit parameters `K3iso0, K3iso1, K3B, K3E` live in it. |
| irrep | Irreducible representation of the finite-volume symmetry group, labelled `nP_irrep`, e.g. `110_A2` = total momentum $\vec n_P=(1,1,0)$, irrep $A_2$. |
| $P_I$, `Vsel` | $P_I$ is the group-theory projector onto an irrep; `Vsel` is an orthonormal basis of its eigenvalue-1 subspace, so the projected matrix is $V^\dagger M V$. |
| non-interacting energy | Three free particles in the box: $E_{\rm lab}=\sum_i\sqrt{m_i^2+\vec p_i^{\,2}}$, $\vec p_i=2\pi\vec n_i/(\xi L)$, $E_{\rm cm}=\sqrt{E_{\rm lab}^2-\vec P^2}$. |
| slogdet | "Signed log-determinant", $\operatorname{sgn}(\operatorname{Re}\,\mathrm{phase}\det M)\cdot\log|\det M|$, computed from the LU (lower-upper) factorization diagonal so it does not underflow. |
| sign flip | Two neighbouring energy grid points where a real quantity (usually $\operatorname{Re}\det$) changes sign. It can be a true zero or a pole. |
| pole / zero | A pole is where a quantity passes through $\pm\infty$; a zero is where it passes through 0. Both flip the sign. |
| coarse grid / refinement | The coarse grid is the uniform energy scan (`coarseN` points). Refinement adds extra energy points inside a sign-flip bracket. |
| cache | Binary file holding $F_3^{-1}$ (and `Vsel`) per energy and irrep, because building $F_3$ is the expensive step. |
| Minuit / Migrad / Hesse | ROOT Minuit2 minimizer (Migrad = minimization algorithm, Hesse = error-matrix step). |
| $\chi^2$ | $\chi^2=r^{\mathsf T}C^{+}r$, $r_i=E_i^{\rm lat}-E_i^{\rm model}$, $C$ = jackknife data covariance, $C^+$ its symmetric pseudo-inverse (`chi_square_v32f`, mode `raw_cov_inv`). |
| OpenMP | Shared-memory thread parallelism (`OMP_NUM_THREADS=18` everywhere). |
| `Lbyas`, $\xi$ | Spatial box size in spatial lattice units, and the anisotropy $\xi=a_s/a_t=3.444$. In temporal lattice units $L=\texttt{Lbyas}\cdot\xi$ (68.88 for `Lbyas`=20). |
| lattice units | All lattice-run energies are in units of $a_t^{-1}$ (e.g. `atmK`$=a_tm_K=0.09698$, `atmpi`$=a_tm_\pi=0.06906$). Toy runs use $M_1=1$. |
| "digonto_classifier_v1/v2", "isozero", "oldviso", "oneShoulder", "p010/p100" | Labels invented during development (not physics terms). Each is defined where it first appears below. |

---

## 1. What was asked

1. A chronological table: version → change → reason (bug or feature) → outcome, for v31 → v32y.
2. The key physics/numerics lessons (projection audit, pole vs zero, slogdet, "isozero", waves0 vs waves01, p010 vs p100, oldviso, classifier v1 → v2 → norm_refine v2 → oneShoulder, iterative refinement, coarse-grid pollution bug (v32w), covariance fix (v32l), cacheonly fix, namespace fix, memory v32x → v32y).
3. Whether a gradient/slope/derivative of the determinant or eigenvalues was used to find roots, and how well it worked.
4. The open questions recorded at the time.

## 2. What I did

- Read the context-transfer package `_nested_old_versions/KKpi_repro_context_transfer_v31u_package/KKpi_repro_context_transfer_v31u/` in full: `CONVERSATION_CONTEXT.md`, `VALIDATION_AND_DEBUGGING.md`, `NEXT_CHATGPT_PROMPT.md`, `MISSING_OR_PLACEHOLDER.md`, the 37 legacy READMEs/Codex prompts in `docs/legacy/`, and the head of `source/F3_GPU_pipeline_version_history.md`.
- Read the guides `KKpi-v32w-v32y-repro-master/.../old_versions_archive/{ARCHIVE_VERSION_GUIDE.md,README_OLD_VERSIONS.md}`, `docs/VERSION_SUMMARY.md`, top-level `README.md`, `tools/README_TOOLS.md`, and every `README_v*`, `PATCH_NOTES_*`, `FIX_NOTES_*`, `CHANGES_*`, `QUICKSTART_LOCAL.md` in `latest_working_versions/v32w_*` and `v32y_*`.
- For each of the 71 unpacked packages under `_nested_old_versions/`, listed its markdown/txt files and read the version-specific ones (`CONVERSATION_CONTEXT.md` changes only at v31v/v31x/v31y/v31z; after v31z the per-package README is stale and the real notes are the accumulated `README_v3*.md` files).
- Where the docs were silent, diffed main sources between consecutive versions: v31z → compilefix → v31za → v31zb (`projected_F3inv_zero_finder_v31z.cpp`), the four v31zm variants, the `K2inv_EREord2_2plus1_mat(...)` call site in v31ze–v31zl, v32j → v32k cache logic, and `ensemble` / `lattice_energy_type` keys in the configs of v32e → v32w. Read the v31zn classifier (`plot_v31zi_F3inv_6panel.py`, lines 112–205) and the v32w fitter (`source/qc_fitter_norm_refine_v2.cpp`, 250 lines).
- Grepped all docs and sources for `gradient|slope|derivative|newton|secant`.

---

## 3. Pre-v31 background (from the v31u context-transfer package)

This history explains why v31 starts where it does.

| Era | What was tried | Outcome recorded |
|---|---|---|
| GPU v7–v8 | Sign flips of normalized $\det F_3$ classified by thresholds `small_y_threshold`, `spike_ratio_threshold`, `slope_ratio_threshold` | Replaced in v9 ("Removed dependence on" those thresholds) |
| GPU v9 | "Inward-shape" rule: if $|f(\text{left})|+|f(\text{right})|$ decreases moving inward → `likely_zero`, if it increases → `likely_pole` | Became the conceptual ancestor of classifier v2 (Section 5.6) |
| CPU v21–v25 | "side_abs_inward_trend" (both sides must decrease, `trend_fraction_required = 0.65`); v24 QC with $\mathcal K_{\rm df,3}$; v25–v26 cached Minuit fitter | Fitter infrastructure later reused |
| v29p–v29r | $\det\big((V^\dagger F_3V)^{-1}\big)$ eigen-branches; Gram–Schmidt re-orthonormalization of `Vsel` | `Vsel` was not exactly orthonormal before GS |
| v30, v30c | Config maker re-ordered "momentum-first" vs "partial-wave-first" | Ordering only; later proven a pure permutation (v31q) |
| v30m/v30n | **Corrected Cartesian $\ell=1$ finite-volume projector** (real-harmonic basis $m=+1\to x$, $0\to z$, $-1\to y$; signed-permutation action). Legacy `P_irrep_projection_2plus1` must not be used for production `Vsel` | Projector health: $\|P_I^2-P_I\|\sim10^{-16}$, leakage $\lesssim10^{-13}$ |
| v30q–v30z | Blind `100_A2`/`110_A2` scans: SVD minimum, eigen-min bisection, residual-checked bisection, pole reject, local minimum ($\min|\lambda|<10^{-6}$) | "Simple sign flips could mix true zeros and pole/asymptote artifacts"; closest-eigenvalue envelope switching was a problem |
| v31a–v31b | Coarse-only rule on the closest-to-zero eigenvalue: $+\to-$ = ZERO, $-\to+$ = POLE, midpoint, no refinement ($N=1501$) | Orientation rule later abandoned (v31zn found zeros with both orientations) |
| v31e–v31l | Minuit fitter with lattice covariance, RAM/binary $F_3^{-1}$ cache, linear-interpolated zero, Hungarian (global) eigen-branch tracking | Infrastructure for v32e–v32w |
| v31m–v31q | Inverse/det conventions; old vs new config maker; full s/p mixing audit; **permutation test** | v31q: $F_2,G,K_2,H,F_3$ and the projector all transform correctly under the permutation; row order is **not** the cause of bad zeros |
| v31r–v31u | Non-interacting shells + group-theory degeneracies; $F_2^{\rm iso},G^{\rm iso}$; raw det → scaled det → **slogdet** | Raw $\det F_2$, $\det G$ underflowed to exactly 0 for `waves_vec_1={0,1}` (v31t); fixed by LU slogdet (v31u) |

---

## 4. Chronological table v31v → v32y (Question 1)

"Bug" = something was wrong; "Feature" = new capability or new test; "Infra" = build/run plumbing. Outcomes are as recorded in the docs; "not recorded" means no result file or note survives.

### Phase A — $F_2$/$G$/$FG$ diagnostics, projection audit, blind pole classifier (default `111_A2` → `110_A2`, $L=20$)

| Version | What changed | Why | Outcome |
|---|---|---|---|
| v31v (+compilefix) | Combined $FG=F_2+G$; 4 panels: $\mathrm{Re}\,FG^{\rm iso}$, slogdet$(FG)$ for waves $\{0\}$ and $\{0,1\}$, slogdet$(V^\dagger FG V)$ | Feature: user wanted the combined object | $FG^{\rm iso}$ poles match non-interacting energies, but slogdet$(FG)$ has **extra** sign changes. Compilefix: `nnP_config` had to be `std::vector<comp>` not `std::vector<int>` for `get_projector_cached_v30q()` |
| v31w | Projection audit file (`*_projection_matrix_Vsel_audit.dat`) naming the function that builds $P_I$, `Vsel`, `Pproj`; sign-flip audit vs $FG^{\rm iso}$ | Feature: find what distinguishes pole-matched flips from extra flips | $P_I$, `Vsel` healthy; plot column map was broken (fixed in v31x) |
| v31x | Event audit per flip: $FG^{\rm iso}$ jump ratio, $\|FG\|$ growth, min/max $|\lambda|$, reciprocal of max-eigenvalue branch, `combined_event_class`; `plot_show=0`; `merge_asof` summary | Bug (plot/summary parsing) + feature | Lessons recorded (Section 5.2): extra flips exist **before** projection; a determinant flip alone is not a pole; reject dimension jumps first |
| v31y | **Blind** pole classifier (does not use the non-interacting list), thresholds `blind_pole_maxeig_abs_min=2e-5`, growth $\ge1.05$, `blind_zero_minabs_max=1e-4`, merge tol $2\times10^{-4}$; oracle comparison afterwards | Feature | No numbers recorded; the v31z context calls it the "technology" to transfer |

### Phase B — projected $F_3^{-1}$ zero finder (`110_A2`, $L=24$, waves $\{0,1\}$)

| Version | What changed | Why | Outcome |
|---|---|---|---|
| v31z (+compilefix) | Target becomes zeros of $\det\big((V^\dagger F_3V)^{-1}\big)$; blind zero = slogdet flip + no dimension jump + finite eigen/singular-value zero signature; terminal prints non-interacting energies and degeneracies | Feature | Compilefix: helper names `trim`, `Options` collided with headers → renamed `trim_v31z`, `OptionsV31z` |
| v31za | Unbuffered stdout, per-stage prints, optional `smoke_test` timing of 3 energies | Bug: run looked frozen | Diagnostic only |
| v31zb | Fixed field names in smoke print (`dd.ok`, `min_abs_eig(...)`) | Bug (compile) | — |
| v31zc | `debug` build `-O0 -g` → `-O2 -g`; no endpoint re-evaluation per flip (`audit_eigenbranch_detail=0`); refine only accepted roots (`refine_iter=28`); `smoke_test=0` | Bug (speed): "-O0 made each matrix evaluation extremely slow" | Usable runtime |
| v31zd | 5-function plot: slogdet$(F_3^{-1})$, slogdet(projected), $\min|\lambda|$, tracked branches, $\sigma_{\min}$; raw eigenvalues saved so Python re-plots without re-running C++ | Feature | Plot decoupled from compute |

### Phase C — toy-mass sanity tests and scattering-parameter mapping

| Version | What changed | Why | Outcome |
|---|---|---|---|
| v31ze | Toy $M_1=1$, $M_2=0.9$, $L=4$, $\xi=1$, $E\in[2.9,3.8]$, 1501 pts, matching double-cover non-interacting tool | Feature: a small, controllable test system | — |
| v31zf | **Swapped** `scatter_params_1/2` at the `K2inv_EREord2_2plus1_mat` call | Suspected bug in channel mapping | Reverted next version |
| v31zg | Swap **undone**; convention fixed as `scatter_params_i_{ell}{index}`, $i$ = spectator: $i=1$ → spectator $M_1$, pair $(M_1,M_2)$; $i=2$ → spectator $M_2$, pair $(M_1,M_1)$; legacy $K_2^{-1}=\mathrm{diag}(K_{2,1}^{-1},\,2K_{2,2}^{-1})$ | Bug (v31zf's swap was wrong) | Mapping documented in `SCATTER_PARAMS_USAGE_v31zg.md` |
| v31zh | Two-case overlay, $M_2=0.5$, $a=0.01$ vs $0.1$ | Feature | — |
| v31zi | Adds $F_{3,\rm iso}^{-1}=1/(v_{\rm iso}^{\mathsf T}F_3v_{\rm iso})$ as 6th panel | Feature: a scalar reference curve | Used as visual "truth" for classifier v1 |
| v31zj, v31zk | Irrep `000_A1m`; then $L=6$, 3501 pts | Feature | — |
| v31zl | **Both runs had silently used the old KKπ default ERE values** (`scatter1_00=4.04`, `scatter1_10=-43.2`, `scatter2_00=4.12`); now zeroed then config applied; values printed | Bug | Terminal check line `[v31zl-scatter] ...` added |
| v31zm (4 variants) | `coarseN=10501`; 1000 extra points inside every coarse flip of four signals; merge + sort by $E$; compilefix ×2 (literal newlines inside strings); plotfix (mask values outside y-window before matplotlib, avoids overflow with $\times10^{300}$ scaling) | Feature + bugs | — |

### Phase D — classifier v1, s vs s+p waves, p-wave strength, $v_{\rm iso}$

| Version | What changed | Why | Outcome |
|---|---|---|---|
| v31zn | First blind classifier for "$F_{3,\rm iso}^{-1}$ zeros" from projected data only ("isozero" = my reading of the package name: zeros of $F_{3,\rm iso}^{-1}$); grey dashed lines | Feature | Rule in Section 5.5; orientation "any" because zeros appeared with both orientations |
| v31zo (+plotfix) | Renamed to **`digonto_classifier_v1`**; `110_A2` | Feature | Called "successful" in v31zq README |
| v31zp | `000_A1m`, waves $\{0\}$ vs $\{0,1\}$, p-wave `scatter1_10=0.01` | Feature: does s+p change the levels? | Not recorded |
| v31zq | Same for `110_A2` | Feature | Not recorded |
| v31zr, v31zs | `000_A1m`, `scatter1_10=0.10` ("p010"); v31zs adds 10 % progress prints | Feature | Not recorded |
| v31zt | **p-wave-safe $v_{\rm iso}$**: only $\ell=0,m=0$ rows non-zero (weights 1 and $1/\sqrt2$) | Bug: legacy $v_{\rm iso}$ put weight on p-wave rows | Affects only the diagnostic $F_{3,\rm iso}^{-1}$, not the classifier |
| v31zu | **"oldviso"**: legacy $v_{\rm iso}$ restored (all `plm` rows weight 1, all `klm` rows $1/\sqrt2$, p-wave rows not zeroed) | User choice (reason not recorded) | Kept in all later toy runs |
| v31zv | `000_A1m`, `scatter1_10=1.00` ("p100") | Feature | Not recorded |
| v31zw | `110_A2`, both runs waves $\{0,1\}$, p010 vs p100 | Feature | Not recorded |

### Phase E — iterative refinement and lattice masses

| Version | What changed | Why | Outcome |
|---|---|---|---|
| v32 | Iterative refinement of $\det(V^\dagger F_3V)^{-1}$ flips: `coarseN=100`, 10 interior points/round, stop when $|E^{(k)}-E^{(k-1)}|<10^{-6}$, max 40 rounds | Feature: cheaper than 10501-point grid | — |
| v32a | Plot file-name fix (`waves_0_1_p010` → `waves_0_1`) | Bug | — |
| v32b | `coarseN=200`; Python sorts all tables by $E$ before classifying | Bug: unsorted merged tables | — |
| v32c | `coarseN`, `refineN` from command line; refinement points OpenMP-parallel | Infra | — |
| v32d | Lattice masses `atmK=0.09698`, `atmpi=0.06906`, $\xi=3.444$, `Lbyas=20`, $E\in[0.26,0.36]$, fixed ERE `4.04, -43.2, 4.12`, `coarseN=1000`, `refineN=20` | Feature: move to the real system | — |

### Phase F — K3df fitter, caching, diagnostics

| Version | What changed | Why | Outcome |
|---|---|---|---|
| v32e | Minuit2 $\mathcal K_{\rm df,3}$ fitter; cached projected $F_3^{-1}$; C++ classifier v1; cutoff $E_{\rm cm}\le0.335$ | Feature | Config used ensemble prefix `szscl21_24_...` (see v32q) |
| v32f | **Convention change**: cache full $F_3^{-1}$ and `Vsel`; $QC_{\rm proj}=V^\dagger(F_3^{-1}+\mathcal K_{\rm df,3})V$; $F_{3,\rm proj}^{-1}=V^\dagger F_3^{-1}V$ (not $(V^\dagger F_3V)^{-1}$) | Feature: correct place to add $\mathcal K_{\rm df,3}$ | Standard from here on |
| v32g, v32h | Minuit2 include/lib discovery; prefer active `$CONDA_PREFIX` | Infra | — |
| v32i | **Namespace fix**: unqualified `K3mat_2plus1(...)` → `k3_2plus1::K3mat_2plus1(...)` with `Kiso={K3iso0,K3iso1}` vector | Bug (compile) | — |
| v32j | `coarseN=10000`, `refineN=0`; load-if-exists-else-build cache | Feature | — |
| v32k | **Cache-only fix**: `require_existing_binary_f3inv_cache=1`; search config path, `cache/cache.bin`, two legacy names, then newest `cache/*.bin`; abort instead of rebuilding | Bug: an hours-long $F_3^{-1}$ rebuild could start when the expected file name was not found | Fitter never silently rebuilds |
| v32l | **Covariance fix**: `MnUserCovariance(row,col)` asserts if row/col $\ge$ `Nrow`; now reads `Nrow()` and zero-fills with a warning | Bug (crash after Migrad) | — |
| v32m | Classifier v1 applied to slogdet of projected QC; guesses all 0.1; 3-panel best-fit plot | Test requested after v32l | Not recorded |
| v32n | No fit; fixed $\mathcal K_{\rm df,3}=0.1$; plot QC diagnostics + lattice + non-interacting | Diagnostic | — |
| v32o | Fixed $\mathcal K_{\rm df,3}=10^{-5}$; 4 panels slogdet$(F_{3,\rm proj}^{-1})$, slogdet$(QC_{\rm proj})$, $\min|\lambda|$, $\sigma_{\min}$; strict irrep filtering | Diagnostic | — |
| v32p | `lattice_energy_type=Ecm` default (no $\sqrt{E^2-P^2}$) | Suspected double conversion | Reversed in v32s |
| v32q | **Ensemble prefix fix**: `szscl21_24_...` → `szscl21_20_...` (v32e–v32p read the $L=24$ mass files while computing at `Lbyas`=20); lattice audit files | Bug | Expected `000_A1m` averages 0.267867, 0.343913, 0.355552, 0.405368, 0.416452 $a_t^{-1}$ (only state 0 below cutoff 0.335) |

### Phase G — classifier v2, normalization, true refinement, fitter baseline

| Version | What changed | Why | Outcome |
|---|---|---|---|
| v32r (+compilefix) | **`digonto_classifier_v2`** on $\det F_{3,\rm proj}^{-1}$: local rescale to $[-1,1]$, `true_zero` if values trend toward 0 near the flip, `pole` if toward $\pm1$; refinement uses 10 neighbouring **cached** points. Compilefix: `clean_label` → `k3df_fit_v32f::clean_label` | Feature | Neighbour points are not new points (fixed in v32s) |
| v32s (+compilefix) | **norm_refine v2**: $M_{\rm scaled}=F_{3,\rm proj}^{-1}/(\texttt{Lbyas}\,\xi)^6$; **true** refinement: 10 new energies per bracket, build new $F_3^{-1}$/`Vsel`, recursive split when several flips, `v2_max_split_depth`; refined cache saved; `lattice_energy_type=En_lab` | Feature + convention reversal | Both shoulders had to decrease |
| v32t | **oneShoulder**: `true_zero` if **either** valid shoulder decreases toward zero | Both-shoulder rule rejected real zeros (inferred from the change) | Rule used in the fitter |
| v32u (+compilefix) | v32t classifier inside the Minuit fitter; missing refinement points built on the fly. Compilefix: link `source/Faddeeva.cc` (`Faddeeva::erfi`) | Feature | — |
| v32v | Verbose progress; per-FCN true_zero/pole/uncertain counts; final recomputation marks fit valid only if all levels found; writes `minuit_fval` and `recomputed_final_chi2` | Debug | Exposed the coarse-grid pollution bug |
| **v32w** | **Fixed coarse grid**: the loaded cache is split into exactly `coarseN` coarse rows (match tolerance $\max(10^{-10},10^{-3}\Delta E)$) and a separate refined cache; initial flip scan uses coarse rows only; failure penalty $10^{100}$ if any level is missing | Bug (Section 5.9) | **Trusted baseline**: $\chi^2=11.4897$, $n_{\rm data}=12$, 4 parameters, $\chi^2/\mathrm{dof}=1.4362$, 12/12 levels found (`docs/VERSION_SUMMARY.md`) |

### Phase H — spectrum vs $L$

| Version | What changed | Why | Outcome |
|---|---|---|---|
| v32x | Integer `Lbyas` scan, fixed $\mathcal K_{\rm df,3}$ = v32w best fit; global cache folder `/media/digonto/Data/F3inv_cache`; one all-irrep cache per `Lbyas` | Feature | Fills RAM and swap at high `Lbyas` |
| v32y (+compilefix) | One `(Lbyas, irrep)` at a time; per-irrep cache files; split coarse/refined by **move** not copy; stream merged cache to a temporary file then atomic rename; free memory before next irrep. Compilefix: `clean_label` qualification again | Bug (memory) | Preferred spectrum-vs-$L$ tool |

---

## 5. Key physics and numerics lessons (Question 2)

### 5.1 Use slogdet, never the raw determinant

For `waves_vec_1={0,1}` the raw $\det F_2$ and $\det G$ came out as exactly 0 (v31t, `VALIDATION_AND_DEBUGGING.md`). This is floating-point underflow (the product of many small LU pivots), not a physical zero. The fix (v31u) is

$$
\log|\det M|=\sum_i\log|U_{ii}|,\qquad
\mathrm{phase}(\det M)=\mathrm{sgn}(\pi)\prod_i\frac{U_{ii}}{|U_{ii}|},\qquad
\mathrm{slogdet}(M)=\operatorname{sgn}\!\big(\operatorname{Re}\,\mathrm{phase}\big)\,\log|\det M| ,
$$

with $PM=LU$ and $\pi$ the row permutation. The v32w fitter keeps both: it decides sign flips from $\operatorname{Re}\det$ and keeps $\log|\det|$ from the LU diagonal (`det_info()` in `qc_fitter_norm_refine_v2.cpp`).

### 5.2 Projection audit: the projector was not the problem

- v31q: old and new config makers produce the same basis up to a permutation; $F_2,G,K_2,H,F_3$ and the projector all transform correctly. Row order does **not** cause bad zeros.
- v31w/v31x audit: $\|P_I^2-P_I\|$, $\|P_I^\dagger-P_I\|$, $\|V^\dagger V-I\|$ and $\|VV^\dagger-P_I\|$ are at machine precision (expected $\sim10^{-16}$ per the v30n/v30r prompts).
- v31x (`v31y/CONVERSATION_CONTEXT.md`): with `waves_vec_1={0}` the **extra** determinant sign flips are already present in the full (unprojected) $FG$. Projection does not create them.
- Mixed s/p is physical: in moving frames with unequal masses, odd and even $\ell$ share an irrep (Wilson, Dudek, Edwards, Thomas, arXiv:1411.2004, cited in `README_v31p_all_QC_mixed_sp_diagnostics.md`). Non-zero $F_2$/$G$ s–p blocks are expected; $K_2$ must stay block-diagonal in $\ell$.
- **Do not** use the legacy `P_irrep_projection_2plus1` for production; use the Cartesian $\ell=1$ projector `fv_projector_cartesian_l1_v30m.hpp` (v30n).

### 5.3 Pole vs zero: a sign flip is not a level

Recorded list of what a sign flip can be (v31u `VALIDATION_AND_DEBUGGING.md`): a true smooth zero, a pole/asymptote, a basis/shell **dimension jump** (the number of momentum states below the cut changes, so the matrix size changes), a closest-eigenvalue branch switch, or an underflow/overflow artifact. Recorded signatures (v31x/v31y):

- reject any interval where raw $N$ or projected dimension changes;
- pole-like: the largest $|\lambda|$ and $\|M\|$ grow, and $1/\lambda_{\max}$ passes through zero;
- zero-like: the smallest $|\lambda|$ (and $\sigma_{\min}$) goes to zero.

The $F_2/G$ degeneracy work (v31r–v31v): $F_2^{\rm iso}$, $G^{\rm iso}$ and $FG^{\rm iso}$ have poles at the non-interacting energies; the group-theory multiplicity of each shell in the irrep is $m_\Lambda=\frac1{|G|}\sum_g\chi_\Lambda(g)^*\chi_{\rm shell}(g)$. In v31t all six `111_A2` shells in the scan had multiplicity 0, yet v31u still plots every shell because shell positions correlate with spikes and dimension jumps. The context marks this character computation as "should be cross-checked" against the projector machinery.

### 5.4 Inverse convention: project the QC, not $F_3$

v31p compared $(V^\dagger F_3V)^{-1}+V^\dagger\mathcal K V$ with $V^\dagger(F_3^{-1}+\mathcal K)V$. From v32f on the code uses

$$
QC_{\rm proj}(E)=V^\dagger\big(F_3^{-1}(E)+\mathcal K_{\rm df,3}(E)\big)V,
$$

and caches the **full** $F_3^{-1}$. The two forms agree only when $V$ spans an exactly invariant subspace (zero leakage). The second form is the one where $\mathcal K_{\rm df,3}$ can be added before projection.

### 5.5 Classifier v1 and "isozero" (v31zn → v32m)

Rule (`classify_projected_iso_zeros_v31zn`, renamed `digonto_classifier_v1` in v31zo). Accept the grid interval $[E_i,E_{i+1}]$ when all of these hold:

1. both points succeeded, with the same raw $N$ and the same projected dimension;
2. $\operatorname{sgn}\operatorname{Re}\det$ of the projected $F_3^{-1}$ flips, **either orientation**;
3. $R_{\lambda}\ge50$ and $R_{\sigma}\ge50$, where for $y=\min|\lambda|$ or $y=\sigma_{\min}$
$$
R=\frac{\max_{\rm core}y}{\operatorname{median}_{\rm shoulders}y},\qquad \text{core}=i-2\ldots i+3,\quad \text{shoulders}=\text{up to 80 points each side outside the core}.
$$

The root is the midpoint. The rule was tuned visually in toy runs ($M_1=1$, $M_2=0.5$, $L=4$ or 6) so that it reproduces the zeros of the scalar $F_{3,\rm iso}^{-1}$ without using that curve as input. The v31zq README calls it "the same successful digonto_classifier_v1".

**Caution [my inference, not stated in the docs]:** v1 asks for a **peak** in $\min|\lambda|$ and $\sigma_{\min}$ at the flip. A zero of the projected $F_3^{-1}$ should give a **dip** toward 0. A peak at a sign flip looks more like the whole small-dimensional matrix becoming large, which is pole-like. Classifier v2 (Section 5.6) uses the opposite logic ("toward zero" = zero, "toward $\pm1$" = pole). The archive does not explain why v1 was dropped, but v32m–v32o were "old classifier" tests followed directly by v2. Before reusing anything tuned on v1 output, check whether v1 was selecting poles of $F_3^{-1}$ (zeros of $F_3$).

### 5.6 Classifier v2 → norm_refine v2 → oneShoulder

Final rule in v32w (`classify()` in `qc_fitter_norm_refine_v2.cpp`):

1. Coarse bracket $[B_L,B_R]$ = adjacent coarse points with opposite $\operatorname{sgn}\operatorname{Re}\det(QC_{\rm scaled})$.
2. Insert `v2_refine_points`=10 equally spaced interior energies. Reuse coarse or refined cache entries when $|E-E'|\le10^{-13}\max(1,|E|)$, otherwise build $F_3^{-1}$ and `Vsel` (OpenMP).
3. If the refined mesh has 0 flips → `uncertain` (`lost_flip_after_refine`). If it has more than 1 flip → split at midpoints between flips and queue the children, up to depth `v2_max_split_depth`=8. If exactly 1 flip at $[F_L,F_R]$ → classify.
4. Rescale by $m=\max_{\rm mesh}|y|$, $\tilde y=y/m$. With $r=$`v2_zero_ratio`$=0.80$:
$$
\text{left}=\big(B_L<F_L\big)\wedge\big(|\tilde y(F_L)|\le r\,|\tilde y(B_L)|\big),\qquad
\text{right}=\big(F_R<B_R\big)\wedge\big(|\tilde y(F_R)|\le r\,|\tilde y(B_R)|\big).
$$
`true_zero` if left **or** right (v32t "oneShoulder"; v32s needed both), else `pole`.
5. Root by one secant/linear step: $E^\ast=F_L-y_L\,(F_R-F_L)/(y_R-y_L)$.

Resolution: coarse spacing $\Delta E=(0.36-0.261)/9999=9.90\times10^{-6}\,a_t^{-1}$ (config `config_v32w_...fixedcoarse.in`). One refinement round gives $\Delta E/11\approx9.0\times10^{-7}\,a_t^{-1}$ before interpolation. There is no further convergence loop in v32w, unlike v32 (which iterated to $10^{-6}$).

**Normalization (norm_refine).** $QC_{\rm scaled}=QC_{\rm proj}/c$ with $c=(\texttt{Lbyas}\,\xi)^{6}=68.88^6\approx1.07\times10^{11}$ (`const_norm_power=6`). For a $d\times d$ block, $\det(M/c)=c^{-d}\det M$ with $c>0$, so **sign flips and the local max-rescaled shape are unchanged**. The normalization only keeps $\det$ inside double-precision range; it does not change physics or classification. The exponent 6 has no recorded derivation.

**Level assignment.** For each irrep the $k$-th lowest `true_zero` is matched to the $k$-th lowest lattice level (`model_for()`). Nothing tracks branches. If a candidate is missing, the model energy is 0 and $\chi^2$ is set to the penalty $10^{100}$.

### 5.7 waves0 vs waves01, p010 vs p100, oldviso

- "waves0" = `waves_vec_1={0}` (s-wave only in the $\pi K$ pair whose spectator is a K). "waves01" = `{0,1}` (s+p). `waves_vec_2={0}` (KK pair) throughout.
- "p010"/"p100" are labels for the p-wave ERE parameter `scatter1_10` = 0.10 / 1.00 (toy units). v31zp/zq used 0.01. Note that the file `README_v31zt.txt` is titled "v31zv … a=0.10" while its body says 1.00. That doc is inconsistent; the configs are correct (v31zt: 0.10, v31zv: 1.00).
- "oldviso" = the legacy $v_{\rm iso}$ (Section 4, v31zu). It only changes the diagnostic curve $F_{3,\rm iso}^{-1}=1/(v_{\rm iso}^{\mathsf T}F_3v_{\rm iso})$. With waves $\{0,1\}$ the legacy vector also weights p-wave rows, so it is not a pure s-wave projection. Classifier v1 does not use it.
- Outcomes (level shifts from adding p-wave) were **not written to any surviving file**.
- Production lattice settings: waves $\{0,1\}$, `scatter1_00=4.04`, `scatter1_10=-43.2`, `scatter2_00=4.12` (fixed, not fitted).

### 5.8 Iterative refinement (v32 → v32c) vs one-shot refinement (v32s → v32w)

v32: start from 100–200 coarse points; each round places 10 points inside the current flip window, re-finds the flip, linearly interpolates, and stops when $|E^{(k)}-E^{(k-1)}|<10^{-6}$ (max 40 rounds). The bracket shrinks by a factor 11 per round. The fitter (v32u/w) instead uses a dense fixed coarse grid (10000) plus one 10-point refinement. The classifier needs **shoulders**, and shoulders only exist if the bracket stays at coarse-grid width (Section 5.9).

### 5.9 The coarse-grid pollution bug (v32v → v32w)

From `PATCH_NOTES_v32w_fixed_coarse_cache.txt`: in v32u/v32v, refinement points built during one Minuit function call were inserted into the main irrep grid. The next call used this denser, mixed grid to find the initial flips. The brackets then collapsed to tiny intervals, the shoulders vanished, and the final recomputation classified **all** QC sign flips as poles. The fit therefore depended on its own history. Fix: an immutable coarse grid of exactly `coarseN` rows for the initial scan, plus a separate refined cache used only inside refinement. The saved file holds both, and they are split again on load.

**General lesson:** any shape-based classifier, meaning a classifier that looks at how values change on each side of a flip, has an energy scale built in. The bracket width must be fixed, or the result is not reproducible.

### 5.10 Covariance fix (v32l), cache-only fix (v32k), namespace fixes (v31z, v32i, v32r, v32y)

- v32l: after Migrad, `MnUserCovariance(row,col)` asserted because the covariance had fewer rows than 4 (for example when Hesse failed or a parameter was fixed). The fix reads `Nrow()` and zero-fills. **Risk:** a zero-filled parameter covariance can look like a valid error; check `valid` and the Minuit status.
- v32k: load before build; abort if no cache. **Risk:** the search falls back to the **newest `cache/*.bin`**, and the binary header holds only a magic string `V32F_F3INV_FULL_CACHE_1`. There is no $L$, $\xi$, masses, ERE parameters, waves or `coarseN` in it (`load_binary_f3inv_cache_v32f`). Since $F_3^{-1}$ depends on all of these, a stale cache is loaded silently.
- Name clashes recurred four times (`trim`/`Options` in v31z, `K3mat_2plus1` in v32i, `clean_label` in v32r and v32y). The reason is that source files are `#include`d into each other (for example `debug_v30r_...cpp` is included by the v32f header). The consolidated repo needs proper namespaces and one compilation unit per tool.
- v32u: `Faddeeva.cc` must be compiled and linked (`Faddeeva::erfi`).

### 5.11 Lattice input conventions (v32p → v32s)

- Ensemble prefix: v32e–v32p used `szscl21_24_...` with `Lval=20`. Any fit from those versions used the wrong volume's data. Fixed in v32q (`szscl21_20_...`).
- Energy frame: v32e–v32o always converted lab → CM ($E_{\rm cm}=\sqrt{E^2-P^2}$, `lattice_data_covariance_cpp.hpp:E_to_Ecm`). v32p–v32r treated the files as already CM. v32s–v32w went back to `En_lab` (convert). The final trusted fit uses **`En_lab`**. For `000_A1m` the two agree ($P=0$), so the v32q audit numbers cannot tell the conventions apart.
- $\chi^2$ uses a pseudo-inverse of the jackknife covariance (`symmetric_pseudoinverse_v32f`) and the cutoff $E_{\rm cm}\le0.335\,a_t^{-1}$, giving 12 levels over `000_A1m 100_A2 110_A2 111_A2 200_A2`.

### 5.12 Memory (v32x → v32y)

Each cache entry stores $F_3^{-1}$ ($N\times N$ complex), $F_{3,\rm proj}^{-1}$ ($d\times d$), `Vsel` ($N\times d$) and the basis lists. An approximate size is
$$
\text{bytes}\approx N_E\cdot16\,(N^2+Nd+d^2)\ \text{per irrep},
$$
with $N_E\ge10^4$ and $N$ growing with `Lbyas`, because more momentum shells fall below the cutoff. v32x held all 5 irreps and then **copied** the cache while splitting coarse from refined, so the peak was at least twice the cache. v32y processes one irrep at a time, moves instead of copying, streams the file out, and frees memory before the next irrep. The archive gives no measured RAM numbers. `tools/check_v32x_L24_cache_ready.py` and `inspect_v32_f3inv_cache_mmap_progress.py` exist because the v32x files were too large to parse naively.

### 5.13 Scattering-parameter plumbing

Three separate bugs: the swap in v31zf (wrong) and its revert in v31zg; the old KKπ defaults that leaked into toy runs (v31zl); and the identical-particle factor $2K_{2,2}^{-1}$ inside the legacy `K2inv_EREord2_2plus1_mat`. Lesson: zero every ERE entry, apply the config, then print the active values. v31zl does this.

---

## 6. Gradient / slope / derivative use (Question 3)

**Short answer.** No version from v31u to v32y uses an analytic or finite-difference derivative of $\det$ or of an eigenvalue to **locate** roots, and there is no Newton or gradient root finder. The root is always found by bracketing a sign flip and then using the midpoint (v31b, v1), bisection (v30s–v30y, to $10^{-9}$), or a single secant/linear step $E^\ast=E_L-y_L(E_R-E_L)/(y_R-y_L)$ (v31l onward, v32w). Slopes appear only in **classification**:

| Where | What "slope" means | Result |
|---|---|---|
| GPU v7/v8 `classify_sign_flips_zero_vs_pole` (in `F3_gpu_omp_cublas_pipeline_v8_*.cu`, still shipped in v32l `v32l_work/source/`) | `slope_ratio` = max/min of $|\Delta y/\Delta E|$ over a $\pm5$-point window, threshold $10^2$ (or $10^1$ at one call site), combined with `spike_ratio` and `small_y_threshold` | Dropped in GPU v9 ("Removed dependence on … slope_ratio_threshold"). The thresholds were hand-tuned and fragile. No numbers recorded |
| GPU v9 "inward shape"; CPU v21 "side_abs_inward_trend" | Sign of $d|y|/dE$ moving inward from both sides; v21 needs $\ge65\,\%$ of steps decreasing on **both** sides | Superseded; the same idea came back in v32s |
| v31x/v31y | Growth ratios of $\lambda_{\max}$ and $\|M\|$ between bracket endpoints ($\ge1.05$), and $1/\lambda_{\max}$ | A ratio test, not a derivative. Blind pole classifier, no outcome numbers |
| v31zn/v1 | Peak/shoulder ratio ($\ge50$) of $\min|\lambda|$ and $\sigma_{\min}$ | Works in toy runs; see the caution in 5.5 |
| v32s/v32t/v32w (v2) | A two-point finite-difference monotonicity test: $|\tilde y|$ must fall by at least 20 % ($r=0.8$) from bracket edge to final sub-bracket on at least one side | Best evidence: v32w finds 12/12 levels and $\chi^2/\mathrm{dof}=1.436$. It fails when the bracket is too narrow (Section 5.9) |
| Very old `pole_searching.h` | `det_F3mat_poles_secant_method` — secant iteration on $\det F_3$ | Legacy, not used in v31–v32 |

**How well it worked.** Slope-threshold classifiers were unstable and were abandoned twice: GPU v9, and again implicitly when v2 dropped peak ratios. The robust part is the monotone-toward-zero test on a **fixed-width** bracket. Whether it is truly robust is limited by three things: one round of refinement, a single interpolation step, and no check of the final residual $|\det|$ or $\sigma_{\min}$, which v30t did have (`REJECT_FINAL_EIGEN_RESIDUAL_TOO_LARGE`, `REJECT_SVD_CROSSCHECK_FAILED`). A natural "gradient" upgrade that was never tried is the sign of $d\lambda/dE$ for the tracked eigenvalue that crosses zero, using Hungarian tracking (v31l) or the eigenvector-overlap bisection from v30s.

---

## 7. Open questions (Question 4)

### 7.1 Recorded in the archive at the time

1. Which sign flips are physical levels: smooth zero, pole, shell jump, or closest-eigenvalue envelope switch? (v31u `NEXT_CHATGPT_PROMPT.md`)
2. Overlay QC zero/pole candidates, dimension jumps, tracked eigen-branches and slogdet$(F_2/G)$; reject candidates near poles or shell jumps unless smooth branch behaviour confirms them. (v31u)
3. Cross-check the little-group character degeneracy code against the exact projector for all irreps. (v31u, known issue 4)
4. Needed for a production solver: final zero-acceptance criteria, validated lattice paths, chosen $\mathcal K_{\rm df,3}$ parameterization, target irreps, **expected number of levels per irrep**. (`MISSING_OR_PLACEHOLDER.md`)
5. For `110_A2` at $L=24$: do the blind candidates match non-interacting energies, and are the extra flips finite zeros, numerical artifacts or projection/basis events? (v31z `NEXT_CHATGPT_PROMPT.md`)
6. The earlier v30d/v30e target: s+p must give the same **six** `111_A2` levels as s-wave only (reference $E_{\rm cm}$ = 0.319305, 0.334219, 0.337600, 0.343683, 0.349761, 0.356261 $a_t^{-1}$, `README_v30d_whole_region.md`), with any extras kept apart as diagnostics.
7. v32x/v32y: produce continuous state lines over `Lbyas`=16–24 (the tools exist; whether lines connect correctly across avoided crossings is not checked).

### 7.2 Noticed during this review (**[my observations]**, need confirmation by running)

1. **v32w `flip()` does not compare `proj_dim`** across the two points. Earlier classifiers (v31y, v31zn) rejected dimension jumps first. A shell opening inside the coarse grid could therefore become a candidate.
2. v1 (peak) and v2 (toward zero) are close to logically opposite (Section 5.5). Results from v31zn–v32q need to be re-examined before they are trusted.
3. The cache header has no physics metadata, and there is a "newest *.bin" fallback (Section 5.10).
4. The fit uses Migrad, a gradient-based minimizer, on a $\chi^2$ that is piecewise constant in places and jumps to $10^{100}$ when a level is lost. The best-fit parameters are $O(10^4$–$10^6)$ starting from guesses of 0.1, and no parameter errors are recorded in `VERSION_SUMMARY.md`. Parameter uncertainties and fit stability (start-point scan) are untested.
5. The lab vs CM energy convention changed three times; the last word (`En_lab`) has no recorded justification. For moving frames this changes the data by up to $E-\sqrt{E^2-P^2}$.
6. Level matching is by order (the $k$-th zero goes to the $k$-th lattice level). A missed or spurious low zero shifts every level above it.
7. The exponent 6 in the normalization only guards against overflow. Its value should not be read physically.
8. "p-wave-safe" vs legacy $v_{\rm iso}$: the legacy definition was kept without a recorded reason.

---

## 8. Recommendations for the consolidated repo (short)

- Keep exactly two production paths: the v32w fitter and the v32y spectrum vs $L$. Archive the rest. Port the fixed-coarse/refined split and the streaming cache.
- Add physics metadata and a hash to the cache header, and refuse to load on mismatch. Remove the "newest *.bin" fallback.
- Restore the dimension-jump rejection in `flip()`. Add a final residual check (relative $\sigma_{\min}$ or $|\lambda|_{\min}$) to `true_zero`, as v30t did.
- Write the conventions into `docs/METHODS.md`: QC projection form, ERE index mapping and the factor 2, lab→CM conversion, $v_{\rm iso}$ definition, normalization.
- Use namespaces and separate compilation; stop `#include`-ing `.cpp` files.
- Acceptance tests: the six `111_A2` s-wave levels (v30d list); the v32q lattice audit averages; the v32w fit ($\chi^2=11.4897$, 12/12 levels); a toy $M_1=1,M_2=0.5,L=6$ `000_A1m` case with known non-interacting levels (`000_A1m_L6_upto10_nonint_levels_summary.txt`).

---

## 9. Problems and assumptions

- No run outputs, fit logs or plots survive in the unpacked packages. `outputs_examples/` is referenced but absent. Every "outcome" is taken from README text, so no numbers could be re-derived. The only quantitative results are the v32w fit numbers, the v32q expected averages and the v30d reference levels, each cited to its file.
- The per-package `README.md` is stale from v31z on (it still describes v31z or v31zo). I used the accumulated `README_v*.md`, `PATCH_NOTES_*` and configs instead.
- The meaning of "isozero" comes from the v31zn package name and the docstring "Blind F3iso-zero candidate classifier"; I found no separate definition.
- The memory formula in 5.12 is a structural estimate from the cache fields, not a measurement.
- Interpretations marked **[my inference]** or **[my observations]** were not checked by running code.

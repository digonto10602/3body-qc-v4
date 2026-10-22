# 002 — M1 step 7 re-planned: cold projector from orbits, threshold basis builder, pair kernels for $\mathcal K_{\rm df,3}$ and $G$, C++ tail, and pole/zero discrimination in the root refiner

- **Date:** 2026-10-21
- **Git commit at planning:** `a25b1e1`
- **Trigger:** `notes/STATUS.md` "Next action" 1 and 4; `reports/tests/006-m1-step6-fast-path-v1.md` (stage timings, notes 1–5); coordinator addendum of 2026-10-21 on step 9 (`results/k2_pwave_shift_L20.json`, `scripts/m1/k2_pwave_shift.py`): `qc3.roots.refine_root` converges to poles of $\det Q$ as well as zeros.
- **Research:** `reports/research/008-step7-speed-replan-orbit-projector-and-kernels.md` (all prototype numbers quoted below come from it unless a command is given).
- **Supersedes:** `prompts/001` **step 7 only**. Steps 8 and 9 and the milestone gate of `prompts/001` §3 stay as written, with the one addition in §3.3 below (step 8's χ² proxy must use the classifying root refiner). Everything else in `prompts/001` (method §2, tolerances §4, interface §5, do-not-retry §6, benchmark protocol §7) remains in force; this prompt only adds to it.
- **Acceptance:** the existing 33 tests of `tests/acceptance/m1/` keep their checks, tolerances and data; the only edits are the `coded_v2` pins of the 2026-10-22 amendment (§4). This plan adds `tests/acceptance/m1/test_m1_fast_v2.py` (9 tests), `tests/acceptance/m1/test_m1_roots.py` (4 tests) and `tests/acceptance/m1/test_m1_production_default.py` (1 test) → 47 tests, to be locked as `acceptance-v4`.
- **Amended 2026-10-22** for Digonto's decisions **D7** (production p-wave form → `q0_consistent`; sub-step 7.8, §4 amendment) and **D8** (exact free-level solutions are artefacts, excluded from D3 matching; §9). Step numbering unchanged.

Jargon: **cold** call = the first fast-path call for a basis (nothing for it in memory yet); **warm** = everything energy-independent for that basis is in memory; **orbit** of a spectator momentum $\mathbf k$ = $\{R\mathbf k: R\in LG(\mathbf P)\}$, $LG$ the little group (the cubic rotations/reflections that leave $\mathbf P$ fixed); **signature** = the ordered list of basis rows $(\text{flavour},\mathbf n,\ell,m)$; **pole** = an energy where $\det Q\to\pm\infty$ because $F_3^{-1}$ has a pole, as opposed to a **zero** (a finite-volume energy level); **pair kernel** = a routine that computes everything shared by one pair of spectator momenta once and then fills the whole $(\ell'm';\ell m)$ block.

## 0. Diagnosis — why step 7 is re-planned

1. **The plan's rung 5 is no longer the bottleneck.** After rungs 1–4, $F$ is 2–16 ms of a 26–100 ms warm point; class tables could save ≤10 ms. Prototyped class counts (research 008 §6) give only a 4.2–5.4× reduction of terms, i.e. ≤1.3× on the total: below the 1.5× stop rule.
2. **The cold projector build dominates real use.** Every basis change (58 distinct bases over the production range: 5/12/16/13/12 per irrep) pays `build_cached_projector_v30q`: it forms dense $N\times N$ complex $S(R)W(R)$ products for every group element and then `diagnose_convention` multiplies all $|LG|^2$ pairs of them (2304 dense products for $O_h$), then eigen-decomposes. 5.4 s at `000_A1m` $N=103$; fast cold vs oracle only 1.9×.
3. **Energy-independent work is redone per point.** The basis builder scans $41^3$ candidates × 2 channels (21 ms, constant) although the thresholds $E^+_{\min}(\mathbf n)$ do not depend on $E$. `K3mat_2plus1` is called 4 times and recomputes the pair kinematics for each of up to 16 entries of a block (43 ms at $N=160$); `G_ij_lm` does the same (20 ms).
4. **Correctness debts.** (a) The legacy projector cache key omits `eig_tol`, `norm_tol`, `proj_tol` (report 006 note 2; stale hit demonstrated). (b) `refine_root` accepts any sign change of $\det Q$: step 9 found two poles returned as roots and three windows with a pole–zero pair and no end-point sign change. (c) The NumPy tail uses OpenBLAS threads unless pinned (note 1).

## 1. Goal of step 7 (revised)

Warm single-core time per point ≈5 ms at $N\approx56$ and ≤25 ms at $N=160$ (≥150× vs the oracle at the largest anchor, expected ≈400×), **cold** time ≤ warm + 0.1 s for every basis, no BLAS threads in the fast path, a complete cache key, and a root refiner that returns zeros only. All existing acceptance tests stay green at every sub-step.

## 2. Method

### 2.1 Projector from orbits (rung 7.2) — `src/qc3_core/fast_projector.hpp`

Conventions are those of the oracle (`fv_projector_cartesian_l1_v30m.hpp`, convention `c0`, arXiv:2111.12734 Eqs. 55–58):
$$
[\hat U(R)]_{\mathbf p\ell'm'i;\,\mathbf k\ell mj}=\delta_{ij}\,\delta_{\mathbf p,R\mathbf k}\,\delta_{\ell'\ell}\,\Pi(R)\,\mathcal D^{(\ell)}_{m'm}(R),\qquad
\hat P_\Gamma=\frac{d_\Gamma}{|LG|}\sum_{R\in LG}\chi_\Gamma(R)\,\hat U(R),
$$
with $R\mathbf k$ = `cubic_transf(k, R)`, $\Pi(R)=-1$ for improper $R$ (`!is_in_rotations_list(R)`, applied for **every** $\ell$ including $\ell=0$ when `parity == -1`), $\mathcal D^{(0)}=1$, $\mathcal D^{(1)}_{m'm}(R)$ the Cartesian signed permutation with axes $(m=1,0,-1)\leftrightarrow(x,z,y)$ (`cartesian_signed_perm_element_v30m`), $U=S\,W$, $\chi_\Gamma$ from `chi()`, $LG$ from `little_group()` — all taken from `src/core/projections_v1.hpp` / the v30m header, **not** re-typed. $\hat P_\Gamma$ is block diagonal in flavour, orbit and $\ell$ with eigenvalues exactly 0 or 1 (2111.12734 App. A.3 Eqs. 104–106).

For a one-dimensional irrep (all five production sectors: $A_1^-\equiv A_{1u}$ of $O_h$, $A_2$ of $C_{4v},C_{2v},C_{3v}$), the rearrangement lemma gives $\hat P_\Gamma\hat U(R')=\chi_\Gamma(R')\hat P_\Gamma$, hence the image of $\hat P_\Gamma$ on a (flavour, orbit, $\ell$) block is
$$
\mathrm{Im}\,\hat P_\Gamma\big|_{\rm block}=\mathrm{span}_{m=-\ell}^{\ell}\;\hat P_\Gamma\,e_{(\mathbf k_0,\ell,m)},\qquad \mathbf k_0=\text{the first momentum of the orbit in basis order},
$$
and $\hat P_\Gamma e_j=\frac1{|LG|}\sum_R\chi(R)\Pi(R)\mathcal D^{(\ell)}_{m'm}(R)\,e_{(R\mathbf k_0,\ell,m')}$ is a vector with $\le|LG|$ non-zero entries, each a multiple of $1/|LG|$. **Algorithm:** for each (flavour, orbit, $\ell$) block, form the $N\times(2\ell+1)$ candidate matrix $C$ by sparse accumulation (cost $O(|LG|)$ per column), take its SVD (`Eigen::JacobiSVD`, thin $U$), keep the left singular vectors with $s_i>10^{-8}s_{\max}$ (measured: smallest kept $0.2887$, largest dropped $1.0\times10^{-16}$), append them to $V_{\rm sel}$ in block order. Total cost $O(|LG|\,N)$; no dense $P$, no eigen-solve. **Self-check (mandatory, always on):** apply the sparse $\hat P_\Gamma$ to the finished $V_{\rm sel}$ and require $\|\hat P_\Gamma V-V\|_F\le10^{-10}\sqrt{N_{\rm proj}}$ and $\|V^TV-1\|_F\le10^{-12}$; otherwise throw (`fast: PROJECTOR_SELF_CHECK`). The check costs $O(|LG|\,N\,N_{\rm proj})$ ≈ 0.1 ms once per basis and makes a wrong candidate set impossible to miss. For $d_\Gamma>1$ (not in production) the candidate set is **all** $|O|(2\ell+1)$ columns of the block; implement the branch, and it is covered only by the self-check.

`eig_tol`, `norm_tol`, `proj_tol` are not used by this construction (the eigenvalues are exactly 0/1; any `eig_tol` in $(10^{-12},0.5)$ selects the same subspace). They stay in the cache key (§2.6) and in the metadata. The eigen method remains available as `options["projector"] = "eigen"` for cross-checks; it must use the same key. $V_{\rm sel}$ columns differ from the oracle's by a unitary mixing inside each block, so comparisons use $VV^\dagger$ (unique), $N_{\rm proj}$, the eigenvalues of $Q$, and $\ln|\det Q|$ (all invariant). Prototype evidence: $N_{\rm proj}$ equal at 20 points, $\|VV^T-P\|/\|P\|\le2.6\times10^{-16}$.

### 2.2 Threshold basis builder (rung 7.3) — `src/qc3_core/fast_basis.hpp`

The oracle keeps $\mathbf n$ iff `En.real() > En_min_plus_for_config(p, P, mi, mj, mk, 0.02, eps_h).real()`. Per (sector, channel), compute once the array $t_c[\mathbf n]$ of these 68 921 doubles with the **same function and the same arguments** (`comp` arithmetic, `0.02`), stored in loop order ($n_x$ outer … $n_z$ inner); per energy, keep $\mathbf n$ iff `E_lab > t_c[n]` — the identical double comparison, so the selection and order are bit-identical to `config_maker_4_momentum_first` by construction. Then expand the $(\ell,m)$ list exactly as the oracle. Memory 0.55 MB per (sector, channel); first use of a sector costs ≈21 ms once (the same scan the oracle does per call). The sorted distinct thresholds converted to $E^*=\sqrt{t^2-P^2}$ are the **dimension jumps**, exposed as `qc3.fast.dimension_jumps` (used by M2's classifier: a sign change of $\det Q$ across a jump is not a root).

### 2.3 Pair kernels for $\mathcal K_{\rm df,3}$ and $G$ (rungs 7.4, 7.5) — `src/qc3_core/fast_k3.hpp`, `fast_g.hpp`

Same formulas as `src/core/K3_functions_2plus1.hpp` (`K3B_element`, `K3E_element`, `K3_element_raw`, `two_plus_one_block_factor`) and `src/core/G_functions_v2.h` (`G_ij_lm`), restructured as: for each **unordered** pair of basis momenta $(\mathbf p,i),(\mathbf k,j)$ compute once the pair kinematics ($\omega$'s, $\sigma$'s, $q^{*2}$'s, the boosts $\mathbf p^*_k,\mathbf k^*_p$ (for $G$) and `pp_vec1`, `kk_vec1`, `psk`, `ksp`, $\gamma_p,\gamma_k$, $V,V'$, the tensor $t$ (for $\mathcal K_E$), the $K^{\rm iso}$ scalars $1$ and $\Delta$), then fill the $(2\ell'+1)\times(2\ell+1)$ block of **all four** basis matrices (`iso0`, `iso1`, `B`, `E`) resp. of $G$ (including the $\sqrt2(-1)^\ell$ flavour factors and the zero $(2,2)$ block), and mirror to the transposed block. Symmetry is a property of the formulas; measured on the oracle: `iso0`, `iso1` exactly symmetric, `B` $\le1.0\times10^{-17}$, `E` $\le8.2\times10^{-17}$, $G$ $\le5.4\times10^{-14}$ (relative Frobenius, worst at `200_A2` 0.34738, $\kappa(F_3)=4.8\times10^8$). The new kernels are **verified, not trusted**: `test_m1_fast_vs_oracle.py` ($10^{-12}$ vs the oracle) and `test_m1_reference.py` logic (the planner's independent `qc3ref` gives the same numbers) both exercise them. $\mathcal K_3(\vec c)=\sum_ic_iB_i$ is formed from the four matrices (no separate pass).

### 2.4 C++ tail, LU solve, lazy conversion (rung 7.6)

In `fast_point.hpp`: $H=\mathcal K_2^{-1}+F+G$ (real), $X_1=H^{-1}F$ by `PartialPivLU`, $F_3=F/3-F X_1$ with the product using $F$'s block-diagonal structure ($4\times4$ and $1\times1$ blocks, $O(16N^2)$); scale $s=\max|F_{3,ab}|$; LU of $F_3/s$; $X=(F_3/s)^{-1}V_{\rm sel}/s$ by `solve` with $N_{\rm proj}$ right-hand sides (the projection of the inverse, `prompts/001` §2.1); $Q=\big(V^TX+V^T\mathcal K_3V\big)/L^6$ with $L=\xi L/a_s$ (invariant 1: scale the matrix, no $N_{\rm proj}$ factor); $(\mathrm{sign},\ln|\det Q|)$ from `PartialPivLU` of $Q$: $\ln|\det|=\sum_k\ln|U_{kk}|$, sign $=\det(\mathsf P)\prod_k\mathrm{sgn}\,U_{kk}$ (invariant 2). The full $F_3^{-1}=$`solve(Identity)`$/s$ is computed **only if** `"F3inv"` is in `want`. Matrices are converted to NumPy only if requested. The fast Python wrapper no longer calls `qc3.qc` (so no BLAS call and no thread issue); `qc3.qc` stays for the oracle and the identity tests. Keys of `timings` unchanged (`"projection"`, `"lu"` now measured in C++).

### 2.5 Pole/zero discrimination (rung 7.1) — `src/python/qc3/roots.py`

`docs/ROOT_FINDING.md` §3: with $\ell(E)=\ln|\det Q(E)|$, zeros are local minima and poles local maxima ($g=\ell'$ jumps $-\infty\to+\infty$ at a zero, $+\infty\to-\infty$ at a pole). After Brent converges to $E^*$ (bracket width ≤ `xtol`), evaluate
$$
\Delta_h(E^*)=\ell(E^*)-\tfrac12\big[\ell(E^*-h)+\ell(E^*+h)\big],\qquad h=10^{-6}\ (0.2\ \text{grid rows}).
$$
For a simple zero at $E_0$, $\ell\approx\ln|E-E_0|+c$ near $E_0$, so $\Delta_h\approx\ln(|E^*-E_0|/h)\le\ln(\text{xtol}/h)=-11.5$; for a simple pole the sign flips; higher order only increases $|\Delta_h|$; the smooth background contributes $O(h\,|\ell'_{\rm bg}|)\ll1$. Decision: **zero** if $\Delta_h\le-D$, **pole** if $\Delta_h\ge+D$, **ambiguous** otherwise, $D=4$ (a factor $e^4\approx55$ in $|\det|$). Measured on the step-9 cases (research 008 §7): zeros $-11.6\ldots-21.6$, poles $+11.9,+12.7$ → margin ≥2.9× in $\ell$ units. $N$ must be equal at $E^*\pm h$ (else `DimensionJumpError`). `refine_root` keeps Brent on $f(E)=\mathrm{sign}\cdot e^{\ell-\ell_0}$ and, by default (`kind="zero"`), raises `PoleError` for a pole and `AmbiguousRootError` for an ambiguous result (both subclasses of `MissingRootError`); `kind="any"` restores the old behaviour. `find_sign_changes` scans $[E_L,E_R]$ on `n_scan` equally spaced points, Brent-refines every sign-change bracket and classifies each; `find_zeros` keeps the zeros. The scan is for the D7 table and M2, **not** for the fit function. Cost of classification: 3 extra evaluations per root.

**Step 8 consequence:** the χ² proxy of `prompts/001` step 8 refines the 12 labelled roots with `refine_root` (default `kind="zero"`) and must report `kind == "zero"` for all 12 (`coded_v2`, $\vec c_{\rm base}$; step 9 shows dips of 13–28 in $\ell$ there). Its evaluation count includes the 3 classification calls per root. It must **not** use `find_zeros` (a scan is M2's root finder and must first pass the FP = FN = 0 gate of `docs/ROOT_FINDING.md` §3.3).

### 2.6 Complete cache key (rung 7.0)

The fast path gets its **own** in-memory projector cache in `src/qc3_core/fast_projector.hpp` (a `std::unordered_map` behind a mutex), keyed by `irrep`, `parity`, `nnP`, `waves_1`, `waves_2`, `eig_tol`, `norm_tol`, `proj_tol`, `projector_method`, and the full basis signature (both channels' $(\mathbf n,\ell,m)$ lists). `qc3.fast.clear_caches()` clears it (and the threshold tables of §2.2); `qc3.fast.cache_info()` reports entry/hit/miss counts. The fast path **stops using** the legacy `g_projector_cache_v30q`. The legacy key is **not** edited (it lives in the oracle's include chain; the oracle clears that cache on every call so it is unaffected); its defect is recorded in `docs/OPEN_ISSUES.md` by the coordinator.

### 2.7 Optional $F$ work (rung 7.7, only under the stop rule)

If after 7.6 $F$ is ≥50 % of the warm time at $N=160$ and a prototype shows ≥1.5× on the **total**: (i) `double` accumulators with pairwise or Neumaier summation instead of x87 `long double` (pairwise error bound $O(\log_2n\,u)$, Higham, SIAM J. Sci. Comput. 14 (1993) 783–799, doi:10.1137/0914050), vectorised inner $z$ loop; (ii) the class tables of `prompts/001` rung 5 (measured term reduction 4.2–5.4×). The unit truncation test ($10^{-14}$ vs `trunc_eps=1e-40`) and the $10^{-12}$ oracle comparison must still pass; an $\exp$ recurrence along $z$ is allowed only with its error ($\le2k\,u\approx1.3\times10^{-14}$ for $k\le60$ steps) demonstrated to stay below $10^{-13}$ on the gamma-stress points of report 006.

## 3. Sub-steps (replacing `prompts/001` §3 row 7)

Order is by correctness first, then by gain on realistic use. Test ids are `tests/acceptance/m1/<file>::<test>`. Builder → tester (mode `step`) → commit, as before. Each sub-step records its stage timings (1 core, `OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 taskset -c 0`, the 4 points of report 006's "warm stage breakdown" + the cold `000_A1m` 0.34538) in the tester report; the formal benchmark stays at the gate (`prompts/001` §7, step 8).

| sub-step | what | done-check | acceptance ids to turn green |
|---|---|---|---|
| **7.0** stubs + cache key | Add the stubs of §5 (raise `NotImplementedError("qc3.<module>.<fn>: not implemented (prompts/002 step 7.k)")`); `options["projector"]`, `fast.cache_info`, `fast.dimension_jumps`, the new `qc3.roots` names, and `refine_root(..., kind=...)` passed explicitly all raise `NotImplementedError` until built; `refine_root` without `kind` keeps its current behaviour until 7.1. Then implement §2.6: fast-owned cache with the complete key, `projector="eigen"` using the **existing** `build_cached_projector_v30q` builder through the new cache, `cache_info`, `clear_caches`. Unit test in `tests/unit/test_fast_path.py`: `eig_tol=0.05` then `1e-30` in one process gives $N_{\rm proj}$ 17 then 5 at `100_A2` 0.2915561 (the report-006 case) and `cache_info()["projector_entries"] == 2`. | `pytest --collect-only -q tests/acceptance` collects 46 tests; new tests fail only with `NotImplementedError` except the one listed; existing 33 green; unit tests green | `test_m1_fast_v2.py::test_projector_cache_key_complete` |
| **7.1** pole/zero discrimination | §2.5 in `qc3.roots`: `RootInfo`, `PoleError`, `AmbiguousRootError`, `DimensionJumpError`, `classify_sign_change`, `refine_root(kind="zero")` default, `refine_root_info`, `find_sign_changes`, `find_zeros`. Unit tests: synthetic $\ell(E)=\ln|E-E_0|$ and $-\ln|E-E_p|$ callables classify correctly; `kind="any"` reproduces the old result. `scripts/m1/k2_pwave_shift.py` switches to `find_zeros`/`refine_root_info` (no own `_is_zero`). | `pytest -q tests/acceptance/m1/test_m1_roots.py -m "not slow"`; `test_fast_root_vs_oracle_100_A2_windows` (locked) still green with the new default | `test_m1_roots.py::test_classify_sign_change_on_step9_cases`, `::test_refine_root_rejects_pole_110_A2_n1_q0`, `::test_refine_root_classifies_12_coded_v2_zeros`; slow: `::test_find_zeros_recovers_hidden_zeros` |
| **7.2** orbit projector | §2.1 in `fast_projector.hpp`; `options["projector"]` default becomes `"orbit"`; metadata `"projector_method"`; builder updates `tests/unit/test_fast_path.py::test_cached_projector_identical_to_cold_build` to compare $VV^\dagger$ (not `array_equal` of $V$) and adds a unit test that `"orbit"` and `"eigen"` give equal $N_{\rm proj}$ and $VV^\dagger$ to $10^{-12}$ at the 3 unit points. | cold `000_A1m` 0.34538: `timings["projector"]` < 1 ms (was 5.4 s); all 33 locked tests green (incl. slow `test_fast_vs_oracle_*`) | `test_m1_fast_v2.py::test_projector_orbit_matches_eigen_and_oracle_subspace`, `::test_projector_cold_build_is_fast` (first assertion) |
| **7.3** threshold basis builder | §2.2 in `fast_basis.hpp`; `make_kinematics` for the fast path uses it (the oracle keeps `config_maker_4_momentum_first`); `qc3.fast.dimension_jumps`. Unit test: basis equal to the oracle's at 10 seeded energies per irrep. | `timings["basis"]` < 0.5 ms warm; cold `000_A1m` 0.34538 total ≤ warm + 0.1 s | `test_m1_fast_v2.py::test_dimension_jumps_match_rule_and_N`, `::test_projector_cold_build_is_fast` (all assertions); slow: `::test_fast_basis_all_signatures_vs_reference` |
| **7.4** $\mathcal K_3$ pair kernel | §2.3 `fast_k3.hpp`; symmetric fill; one pass for the four matrices. | `timings["K3"]` ≤ 5 ms at `111_A2` 0.35103; `test_fast_vs_oracle_*` green (fast and slow) | (locked ids of step 6, re-run) |
| **7.5** $G$ pair kernel | §2.3 `fast_g.hpp`; symmetric fill. | `timings["G"]` ≤ 3 ms at `111_A2` 0.35103; same tests | (locked ids of step 6, re-run) |
| **7.6** C++ tail | §2.4; `qc3.fast` no longer imports `qc3.qc`; metadata `"tail": "cpp"`; conversions on demand. | `test_fast_metadata_and_timings` green; process CPU/wall ≤ 1.1 without `OPENBLAS_NUM_THREADS` (tester repeats the report-006 check); 111_A2 0.35103 warm ≤ 25 ms | `test_m1_fast_v2.py::test_fast_tail_consistent_with_numpy_and_want_subset`; slow: `::test_fast_speedup_single_thread_largest_point_ge_150x` |
| **7.7** optional $F$ | §2.7, only if the stop rule allows. | stop rule | (locked ids, re-run) |
| **7.8** D7: production default → `q0_consistent` (**mandatory, before step 8**) | `qc3.constants.K2_PWAVE_FORM_DEFAULT = "q0_consistent"`; `Sector.k2_pwave_form` defaults to it (so `Sector.production` and plain `Sector(...)` both give it); `oracle.K2_PWAVE_FORMS` unchanged (both selectable); comments/docstrings in `oracle.py`, `fast.py`, `sector.py` updated; `docs/METHODS.md` §5 already says so. Builder updates the unit tests `tests/unit/test_sector_oracle.py::test_sector_properties` (line 15) and `::test_k2_pwave_forms_differ_only_in_channel1_pwave_diagonal` (line 72) to the new default. `scripts/m1/k2_pwave_shift.py` is unaffected (explicit forms). **Step 8 consequence:** `qc3.bench.run` gets a `k2_pwave_form` argument defaulting to the production default and records it in `results/bench_*.json`; the χ² proxy (12 labelled roots) must pin `k2_pwave_form="coded_v2"`, because the labels are `coded_v2` zeros (step 9: with `q0_consistent` 5 of 12 windows contain no end-point sign change or a pole); timing is form-independent. | `python3 -c "from qc3.sector import Sector; print(Sector('100_A2').k2_pwave_form)"` → `q0_consistent`; all 46 other acceptance tests green (they pin `coded_v2`); unit tests green | `test_m1_production_default.py::test_production_default_k2_pwave_form_is_q0_consistent` |

**Stop rule (unchanged):** a rung must gain ≥1.5× on the warm total at the point where it applies, or it is dropped and recorded. **Step 7 done** when 7.0–7.6 and 7.8 are green; then `prompts/001` steps 8 (with §2.5's and 7.8's consequences) and 9 (done, `reports/tests/007`; its script adopts `find_zeros` in 7.1) and the milestone gate.

## 4. Acceptance tests added (`acceptance-v4`; nothing existing changes)

| test | checks | tolerance | justification |
|---|---|---|---|
| `test_m1_fast_v2.py::test_projector_cache_key_complete` | `projector="eigen"`: `100_A2` 0.2915561 with `eig_tol` 0.05 then $10^{-30}$, and in reverse order after `clear_caches()`: $N_{\rm proj}$ order-independent, equal to the oracle's (17) at 0.05; `cache_info()["projector_entries"] == 2` | exact integers | a key omission is a discrete fault (report 006: 17 vs 5) |
| `::test_projector_orbit_matches_eigen_and_oracle_subspace` | 4 small anchor points, $\vec c_{\rm base}$: $N_{\rm proj}$ equal (orbit, eigen, oracle); $\|V_oV_o^\dagger-V_eV_e^\dagger\|_F/\|V_eV_e^\dagger\|_F$ and vs the oracle's $VV^\dagger$; sign equal, $\lvert\Delta\ln\lvert\det Q\rvert\rvert$ orbit vs eigen; `metadata["projector_method"]=="orbit"`; $\|V^\dagger V-1\|$ | $10^{-12}$; $10^{-6}$; $10^{-12}$ | $VV^\dagger$ is unique; entries of $\hat P$ are multiples of $1/|LG|$; measured $6.5\times10^{-16}$ (margin 1500). $\ln\det$: same LU round-off floor as `prompts/001` §4 (`test_fast_vs_oracle_small_points`), measured $1.9\times10^{-9}$ worst (margin 500). Orthonormality: measured $1.2\times10^{-15}$ |
| `::test_projector_cold_build_is_fast` | `clear_caches()`; `000_A1m` 0.34538 ($N=103$, $|LG|=48$, the 5.4 s case): first call is a miss, `timings["projector"]` ≤ 0.02 s; second call hits; cold total − warm total ≤ 0.1 s | 0.02 s; 0.1 s | orbit build is $O(|LG|N)\approx5\times10^3$ operations plus a $103\times3$ SVD → <1 ms even in debug builds (margin ≥20); cold extras after 7.3 = threshold tables (≈21 ms, the oracle's own scan) + build, margin 4. Machine-independent to within these margins; the old value is 270× above the bound |
| `::test_dimension_jumps_match_rule_and_N` | per irrep: `dimension_jumps(sector, 0.26310, 0.36)` equals the thresholds computed in-test from the §2.2 rule (NumPy, $E^*=\sqrt{t^2-P^2}$), sorted, strictly inside; counts 4/11/15/12/11 reproduced by the rule, not hard-coded; $N$ from `evaluate_point` at $E_j\pm2\times10^{-9}$ increases across every jump | $10^{-12}$ relative on energies | same float64 formula, different operation order → $\lesssim10\,u\approx10^{-15}$ (margin 1000); $2\times10^{-9}$ is $10^{7}$ above the round-trip error of $E^*\leftrightarrow E_{\rm lab}$ and 2000× below the grid step |
| `::test_fast_basis_all_signatures_vs_reference` (slow) | at the midpoint of every one of the 58 constant-$N$ segments, `evaluate_point` basis equals `qc3ref.basis.expected_basis` | exact | the generalisation requirement: every basis the production range can produce |
| `::test_fast_tail_consistent_with_numpy_and_want_subset` | 4 small points, $\vec c_{\rm base}$: `qc3.qc.slogdet(f.Q)` vs `(f.sign, f.lnabsdet)`; `want=("Q",)` call (LU solve path, no `F3inv`) vs full `want` (which also forms $F_3^{-1}$): same sign, $\lvert\Delta\ln\lvert\det\rvert\rvert$; `F3inv is None` for the subset; `metadata["tail"]=="cpp"` | $10^{-6}$ | Eigen vs LAPACK LU and solve-vs-inverse differ by round-off only: `prompts/001` §4 floor $1.4\times10^{-9}$, measured $1.9\times10^{-9}$ at the worst anchor (margin 500) |
| `::test_fast_speedup_single_thread_largest_point_ge_150x` (slow) | `111_A2` 0.35103 ($N=160$), $\vec c_{\rm base}$, 1 thread, warm: median of 3 fast vs 1 oracle call (the oracle's own projector build excluded via its `timings`) | ≥150 | expected ≈400× (8.4 s / ≤25 ms); margin 2.7× for machine noise; ratio measured in one process on one machine |
| `::test_fast_single_thread_cpu_time` | `000_A1m` 0.26932 warm, 20 calls: `time.process_time` / wall ≤ 1.3 with **no** `OPENBLAS_NUM_THREADS` pinning inside the test | 1.3 | report 006 measured 3.9 with the NumPy tail and 0.96–0.99 without; 1.3 leaves room for the interpreter |
| `::test_fast_options_projector_values` | `options={"projector":"bogus"}` raises `ValueError`; metadata has `projector_method`, `projector_cache_hit`, `tail` | exact | interface |
| `test_m1_roots.py::test_classify_sign_change_on_step9_cases` | `classify_sign_change` at the 9 step-9 points (7 zeros, 2 poles; `q0_consistent`, $\vec c_{\rm base}$, fast backend): kind as recorded, $\lvert\Delta_h\rvert\ge8$, $N$ equal at $E^*\pm h$ | kind exact; $\lvert\Delta_h\rvert\ge8$ | measured $11.6$–$21.6$ (research 008 §7); exact-arithmetic bound 11.5 for xtol $10^{-11}$, $h=10^{-6}$; 8 leaves 1.45× for a less-converged $E^*$ |
| `::test_refine_root_rejects_pole_110_A2_n1_q0` | window = accepted ±250 rows of `110_A2` level 1, `q0_consistent`: `refine_root` raises `PoleError`; `refine_root_info(kind="any")` returns kind `"pole"` at $E=0.3161291111\pm10^{-8}$ | $10^{-8}$ | same backend, Brent xtol $10^{-11}$; a $10^{-6}$ change of $\ln\lvert\det\rvert$ moves a simple root by $\lesssim10^{-13}$; $10^{-8}$ is 500× below the grid step |
| `::test_refine_root_classifies_12_coded_v2_zeros` | the 12 accepted windows inside cutoff 0.335 (`inside_Ecm_cutoff`), `coded_v2`, $\vec c_{\rm base}$: `refine_root_info` kind `"zero"`, $E$ within $5\times10^{-6}$ of `previous_model_zero`, $\Delta_h\le-8$ | as listed | the step-8 proxy's precondition; $5\times10^{-6}$ = label resolution (`prompts/001` §4) |
| `::test_find_zeros_recovers_hidden_zeros` (slow) | widened windows of `100_A2` n0, `110_A2` n0, `200_A2` n2 (`q0_consistent`): `refine_root` raises `MissingRootError` (no end-point sign change); `find_zeros(n_scan=241)` returns exactly one zero at the step-9 value; `find_sign_changes` contains ≥1 pole | $10^{-8}$ on $E$ | step-9 values (`results/k2_pwave_shift_L20.json`, commit recorded there) are golden values with provenance; same justification as above |

| `test_m1_production_default.py::test_production_default_k2_pwave_form_is_q0_consistent` | `Sector(...)`, `Sector.production(...)` and `qc3.constants.K2_PWAVE_FORM_DEFAULT` give `"q0_consistent"`; both forms stay selectable; a fast-path result records `metadata["k2_pwave_form"] == "q0_consistent"`; vs an explicit `coded_v2` Sector the channel-1 $s$-wave diagonal of $\mathcal K_2^{-1}$ is identical and the $p$-wave diagonal differs | exact | D7; the two forms differ only in the channel-1 $\ell=1$ diagonal (`prompts/001` §2.2) |

No existing tolerance changes. If any of the 9 step-9 points comes out **ambiguous**, stop and report (do not change $D$ or $h$).

### Amendment 2026-10-22 (acceptance-v4): D7 pins in the locked tests, new default test

Follows Digonto's decision **D7**: production `k2_pwave_form` → `q0_consistent`. The locked tests must keep testing what they test today, so `m1_helpers.sector()` pins `k2_pwave_form="coded_v2"` unless overridden (anchors JSON, labelled windows and all step-6 comparisons were produced with `coded_v2`); `test_m1_oracle.py` keeps its `coded_v2` metadata assertion (now the pinned form; comment added); `test_m1_reference.py` passes the form explicitly in both branches; `test_m1_fast_vs_oracle.py::test_fast_root_vs_oracle_100_A2_windows` pins `coded_v2` explicitly. The new tests of this plan: `test_m1_roots.py` was already explicit everywhere (`q0_consistent` for the step-9 cases, `coded_v2` for the 12 labelled zeros); `test_m1_fast_v2.py` goes through the pinned helper (form-agnostic checks). Full wording and the "does not weaken" argument: `prompts/001` §4, amendment of 2026-10-22. No tolerance or data change.

## 5. Interface additions (stubs in 7.0; existing names unchanged)

### 5.1 `qc3.fast`
```python
DEFAULT_OPTIONS = {"trunc_eps": 1e-16, "use_classes": False, "threads": 1, "projector": "orbit"}   # "orbit" | "eigen"
def evaluate_point(sector, Ecm, c=(0,0,0,0), want=ALL, options=None) -> PointResult   # as before; unknown option keys -> ValueError, bad "projector" value -> ValueError
def dimension_jumps(sector, E_min: float, E_max: float) -> np.ndarray   # sorted E_cm strictly inside (E_min, E_max) where N changes (either channel)
def cache_info() -> dict      # {"projector_entries": int, "projector_hits": int, "projector_misses": int, "threshold_entries": int}
def clear_caches() -> None    # projector cache and threshold tables
```
`metadata` gains `"projector_method"` (`"orbit"`/`"eigen"`), `"tail"` (`"cpp"`), keeps `"projector_cache_hit"`, `"use_classes"`, `"trunc_eps"`, `"threads"`, `"k2_pwave_form"`.

`qc3.constants` (sub-step 7.8): `K2_PWAVE_FORM_DEFAULT: str = "q0_consistent"`; `qc3.sector.Sector.k2_pwave_form` defaults to it (`PRODUCTION` dict unchanged; `oracle.K2_PWAVE_FORMS` unchanged). `timings` keys unchanged. `PointResult.F3inv` is `None` unless `"F3inv"` in `want`; `Q` is always computed and returned when `"Q"` in `want`.

### 5.2 `qc3.roots`
```python
class PoleError(MissingRootError): ...
class AmbiguousRootError(MissingRootError): ...
class DimensionJumpError(MissingRootError): ...
@dataclass
class RootInfo:
    E: float; kind: str                 # "zero" | "pole" | "ambiguous"
    lnabsdet: float; lnabsdet_minus: float; lnabsdet_plus: float; discriminant: float   # Delta_h
    N: int; bracket: tuple; n_evaluations: int
def classify_sign_change(evaluate, sector, c, E_star, h=1e-6, D=4.0) -> RootInfo
def refine_root(evaluate, sector, c, E_left, E_right, xtol=1e-11, max_iter=200, kind="zero", h=1e-6, D=4.0) -> float
    # kind="zero": PoleError / AmbiguousRootError if the sign change is not a zero; kind="any": old behaviour (returns E of any sign change)
def refine_root_info(evaluate, sector, c, E_left, E_right, xtol=1e-11, max_iter=200, kind="zero", h=1e-6, D=4.0) -> RootInfo
def find_sign_changes(evaluate, sector, c, E_left, E_right, n_scan=241, xtol=1e-11, h=1e-6, D=4.0) -> list[RootInfo]   # all, sorted by E
def find_zeros(evaluate, sector, c, E_left, E_right, n_scan=241, xtol=1e-11, h=1e-6, D=4.0) -> list[RootInfo]          # kind == "zero"
```
`evaluate` is called as `evaluate(sector, E, c, want=("Q",))`, as now.

### 5.3 C++ (`src/qc3_core/`)
`fast_projector.hpp` (`struct FastProjector {Eigen::MatrixXd V; int Nproj; std::string method;}`, `fast_projector_cached(const Kin&, const Input&, method) -> const FastProjector&`, cache map + counters, `clear_fast_caches()`), `fast_basis.hpp` (`threshold_table(sector-key, channel) -> const std::vector<double>&`, `fast_kinematics(const Input&) -> Kin` producing the same `Kin` as `make_kinematics`, `dimension_jumps(...)`), `fast_k3.hpp` (`k3_basis_pairs(const Kin&, std::array<MatR,4>&)`), `fast_g.hpp` (`G_pairs(const Kin&, MatR&)`). `module.cpp`: `fast_point(params, trunc_eps, projector, want_flags)`, `dimension_jumps`, `cache_info`, `clear_caches`. The oracle files `oracle_point.hpp`, `src/core`, `src/fitter`, `src/cachegen_cpu` stay byte-unchanged (tester checks `git diff --quiet`).

## 6. Do-not-retry (additions to `prompts/001` §6)

| approach | why not |
|---|---|
| Rebuilding the projector with `build_cached_projector_v30q` (dense $S W$ products + `diagnose_convention` closure test + eigen-solve) on every new basis | 0.03–5.4 s per basis; the closure diagnostics are $O(|LG|^2N^3)$ and the eigen-solve is unnecessary (eigenvalues are exactly 0/1) |
| Precomputing the 58 projectors to disk | D4 (no caches); unnecessary once a build is <1 ms |
| Re-scanning $41^3$ candidates per energy in the basis builder | thresholds are energy-independent (§2.2) |
| Four `K3mat_2plus1` calls with per-entry kinematics; per-entry boosts in `G_ij_lm` | 16× redundant kinematics per block, no symmetry |
| NumPy/OpenBLAS tail inside the fast path | spawns threads unless pinned (report 006 note 1) |
| Full $N\times N$ inverse when only $Q$ is wanted | $N_{\rm proj}$ right-hand sides suffice (`prompts/001` §2.1) |
| Accepting any sign change of $\det Q$ as a root; window-end "dip" tests with a width-dependent threshold | poles (step 9); the local $\Delta_h$ test is width-independent |
| Class tables as the next rung for $F$ | term reduction only 4.2–5.4×; ≤1.3× on the total |

## 7. Risks and stop conditions

- **Orbit projector disagrees with the eigen method at some production basis** ($N_{\rm proj}$ or $VV^\dagger$): that would mean $\hat U$ is not a representation for that basis (a convention or physics problem, not a speed problem) → stop, report the basis, ask Digonto. The always-on self-check turns this into an exception at the first such point.
- **$d_\Gamma>1$ irreps** are implemented but untested (none in production); say so in the metadata if ever used.
- **Symmetric fill of $G$** differs from the oracle by ≤$5.4\times10^{-14}$ at ill-conditioned points (margin 18 against $10^{-12}$). If a future point exceeds $10^{-12}$, compute both orderings instead of mirroring (no tolerance change).
- **Transcription errors in the pair kernels** are caught at $10^{-12}$ by the locked tests; the builder must not "fix" a mismatch by adjusting the kernel to the oracle's numbers but by finding the formula error (both the oracle and `qc3ref` are available).
- **Timing tests** (`cold_build`, `ge_150x`, `cpu_time`) can fail on a heavily loaded machine; the tester re-runs once pinned before declaring FAIL.
- **Ambiguous classification** ($|\Delta_h|<4$) at any labelled window: stop and report (pole and zero closer than $\sim10^{-6}$, or a double root); do not tune $D$ or $h$.
- **Stop and ask Digonto:** any of the above that is not a code bug; any need to change invariant-3 constants; if step 8's proxy finds a labelled window whose root is not a zero under `coded_v2`.

## 8. Expected cost

- Builder: 7.0 ≈ 1 h; 7.1 ≈ 2 h; 7.2 ≈ 3 h; 7.3 ≈ 2 h; 7.4 ≈ 4 h (the largest transcription); 7.5 ≈ 2 h; 7.6 ≈ 2 h.
- Run time (1 core): warm ≈5 ms at $N=56$, 18–25 ms at $N=160$ (from report 006: basis 21→0.05, $\mathcal K_3$ 43→≤5, $G$ 20→≤3, tail 0.6→0.2, rest unchanged); cold = warm + ≤22 ms (first use of a sector) or + <1 ms (new basis of a known sector). Memory: +5.5 MB threshold tables for all 5 sectors; projector cache ≤ 58 × $N\,N_{\rm proj}$ doubles ≈ 3 MB.
- New tests: fast part +≈25 s (12 coded_v2 refinements ≈ 130 fast calls, 9 classifications, cache tests); slow part +≈2 min (3 scans of 241 points, 58 reference bases ≈ 30 s, one oracle call at $N=160$ ≈ 8 s).

## 9. M2 outline additions (D7, D8; `prompts/001` §9 otherwise unchanged)

- **D7 (2026-10-22):** all M2 fits run with the production default `q0_consistent`. The human-labelled windows (invariant 7) are `coded_v2` zeros; step 9 (`reports/tests/007`, `results/k2_pwave_shift_L20.json`) shows the `q0_consistent` zeros shift by $-0.1$ to $-0.4\,\sigma_{\rm lat}$ and, in 5 of 12 windows, leave the ±250-row window or share it with a pole. M2's root finder therefore cannot start from the labelled brackets with `q0_consistent`; it must find and classify roots itself (`find_sign_changes` + the $g$-jump/`$\Delta_h$` test) and pass the FP = FN = 0 gate of `docs/ROOT_FINDING.md` §3.3 scored against the labels **with `coded_v2`** (the form the labels belong to), before being used with `q0_consistent`.
- **D8 (2026-10-22):** exact free-level solutions of the QC are artefacts and are **excluded from D3 level matching**. Example (`reports/tests/007` finding 2): in `110_A2` at $E^*=0.319468828$, $|\det Q|\propto(E-E^*)^4$ with no sign change, exactly the non-interacting energy of $K(\mathbf 0)K(0,1,0)\pi(1,0,0)$ with $\mathbf P=(1,1,0)$; the projected $\mathcal K_{\rm df,3}$ and $\mathcal K_2$ do not lift it, and it is the same for both forms. Rule for M2: the detector of even-order zeros ($g=\ell'$ jumps $-\to+$ with no sign change) must compare every candidate with the non-interacting spectrum of the sector (`tools/nonint3body_general_doublecover/`, `docs/METHODS.md` §14) and flag a candidate within $10^{-7}$ of a free level as `free_level`; `match_levels` (D3) receives only `zero` kinds. The labels stay ground truth (they contain no free-level solutions). No M1 test is added for D8 (it needs the free-level enumerator, an M2 interface); the M1 classifier already reports the kind, so M2 only adds the `free_level` kind.

# Tester report 002: M1 step 2, oracle wrapper

- **Date:** 2026-10-21
- **Git commit:** `59b07e1` (working tree clean except `graphify-out/`; `_core.git_commit()` = `59b07e1`)
- **Prompt in force:** `prompts/001-m1-cache-free-single-energy-qc.md` (§3 row 2, §4, §5.2-5.4)
- **Mode:** step
- **Related:** `reports/tests/001-m1-step1-build-stubs.md`, `reports/research/007-single-energy-qc-methods-and-reference.md`
- **Raw numbers:** `results/tester_002_m1_step2_oracle.json` (commit, parameters, machine; no RNG, so no seed)
- **Machine:** 4 vCPU Linux 6.18, g++ 13.3.0 `-O3 -fopenmp`, Python 3.13.16, NumPy with scipy-openblas 0.3.34 (Haswell kernels)

## Verdict: PASS (with two notes for the auditor, section "Notes")

## Commands run
```
rm -f src/python/qc3/_core*.so; ./build.sh core
python3 -m pytest -q tests/acceptance/m1/test_m1_oracle.py -p no:cacheprovider --durations=0
python3 -m pytest -q tests/acceptance/m1 -m "not slow" -p no:cacheprovider -rA --tb=line
python3 -m pytest -q tests/acceptance/m1/test_m1_reference.py::test_ingredients_vs_reference_large_points -p no:cacheprovider
python3 -m pytest -q tests --ignore=tests/acceptance -p no:cacheprovider -v
<scratchpad>/check6.py, thr.py, tim.py   (independent scripts; not added to tests/)
OPENBLAS_NUM_THREADS={1,4}, OMP_NUM_THREADS={1,4}, both=4  with test_oracle_anchors_k3_baseline_small_points
.claude/bin/acceptance-check
```
`graphify affected "Sector"` lists `oracle.py`, `fast.py`, `point.py`, `lattice.py`, `qc.py`, `m1_helpers.sector`, and the 5 unit tests. Everything listed is exercised by the runs above. (`graphify affected "evaluate_point"` / `"oracle_point"` give "No unique node match".)

## Results

| # | check | expected | got | tolerance | result |
|---|---|---|---|---|---|
| 1 | clean rebuild `./build.sh core` | compiles | OK, **90.9 s** (step 1 empty module: 7.4 s; oracle pulls in the whole fitter header chain) | n/a | PASS |
| 2 | `test_m1_oracle.py` incl. slow | 7 green | **7 passed in 72.3 s** | n/a | PASS |
| 2a | per-test time | | anchors_all_points (slow, 12 pts) 49.5 s; k3_baseline_small 9.9 s; anchors_small 9.8 s; basis_rule_and_order 3.1 s; other 3 tests <5 ms (cached oracle results) | n/a | info |
| 3 | `m1 -m "not slow"` tally | 9 pass, rest NotImplementedError, one ValueError | **9 passed, 18 failed, 6 deselected.** Failures: 17 `NotImplementedError` (lattice 8, qc_identities 4, fast 3, portability 2) + 1 `ValueError` (`test_K2inv_pwave_both_forms_vs_reference`: "q0_consistent not available yet (step 3)") | n/a | PASS (as expected; the 18 are PENDING, not failures) |
| 3a | informational: `test_ingredients_vs_reference_large_points` (slow, step 3 id) | pass | **passed** (20.3 s) | $10^{-12}$ | PASS (info) |
| 4 | unit tests | 5 pass | **5 passed** (4.3 s) | n/a | PASS |
| 5 | code review of `oracle_point.hpp` | same calls, same order, no copied formulas | see below | n/a | PASS |
| 6 | Q (NumPy) vs Q (Eigen) | round-off | $\|Q_{\rm np}-Q_{\rm Eig}\|/\|Q\| = 9\times10^{-20}\ldots2.7\times10^{-16}$ over all 16 (point, $c$) cases | $\sim10^{-15}$ | PASS |
| 6b | $\lvert\ln\lvert\det\rvert_{\rm np}-\ln\lvert\det\rvert_{\rm Eig}\rvert$ vs $\kappa(Q)\,\varepsilon$ | $\ll \kappa\varepsilon$ | worst: 100_A2 0.2915561 base $1.1\times10^{-9}$ vs $\kappa\varepsilon=4.1\times10^{-4}$; 200_A2 0.34738 $5.2\times10^{-10}$ vs $5.2\times10^{-8}$; all others $\le 8\times10^{-13}$ | $\kappa\varepsilon$ | PASS |
| 6c | 100_A2 0.2915561 at $\vec c_{\rm base}$: $\kappa(Q)$ and margin | builder: $3.8\times10^{12}$, ~7x margin | $\kappa(Q)=3.76\times10^{12}$ confirmed; $\Delta\ln\lvert\det\rvert$(np-anchor) $=-1.378\times10^{-10}$ vs $10^{-9}$: **margin 7.3x** | $10^{-9}$ | PASS |
| 6d | same test vs thread counts | identical | OPENBLAS=1, OPENBLAS=4, OMP=1, OMP=4, both=4: $\Delta=-1.378\times10^{-10}$ and $\ln\lvert\det\rvert=-118.960684995831542$ in all five; test passes in all five | $10^{-9}$ | PASS |
| 7 | single thread, restored | cpu/wall ~ 1; ICV restored | cpu/wall = 0.99-1.18 with default env and with OMP/OPENBLAS=4; `omp_get_max_threads()` = 4 before and after a full run (restored) | n/a | PASS |
| 7b | projector cache cleared per call | same projector time on a repeat | 100_A2 0.2915561: 35.3 ms then 34.8 ms (no cache hit) | n/a | PASS |
| 8 | `acceptance-check` | PASS | "PASS: tests/acceptance matches acceptance-v1" | n/a | PASS |

Counts: step-2 ids 7/7 pass; whole m1 fast set 9 pass / 18 pending (0 unexpected); unit 5/5; new failures vs last report: 0.

### Check 5: review of `oracle_point.hpp` against `build_component_entry` (cachegen L223-288)
Same functions, same order, same arguments (read line by line):
`Ecm_to_E` -> `config_maker_4_momentum_first` x2 (flavour 1 then 2) -> `F2_2plus1_mat` -> `K2inv_EREord2_2plus1_mat` -> `G_2plus1_mat` -> $H=K_2^{-1}+F+G$ -> `PartialPivLU(H).solve(F2)` -> $F_3=F/3-F\,X$ -> `get_projector_cached_v30q` -> $F_3/\max|F_3|$, `PartialPivLU`, `solve(1)`, `/= f3_scale`. `make_base_physics(FitSettings)` is used for the parameters. The four $K_3$ basis matrices call `k3_2plus1::K3mat_2plus1` exactly as `precompute_projected_k3_basis` does (unit `Kiso`/`KB`/`KE`, `'n'`). Only differences: the oracle projects $F_3^{-1}+K_3$ together instead of $V^\dagger F_3^{-1}V + V^\dagger K_3 V$ separately (algebraically identical; $\|Q_{\rm rebuild}-Q\|/\|Q\|\le 2\times10^{-16}$).
- No formulas copied; all physics comes from `#include` of the existing fitter header.
- Constants not exposed in `Input` (`tolerance`=1e-12, `norm_tol`=1e-12, `proj_tol`=1e-10, `Q0norm`=true, `sort_orbit_flag`=false, remaining scatter entries 0) stay at the `FitSettings` defaults. These equal the values in `configs/cachegen_cpu/config_cpu_component_cachegen.in`. They are not written into `metadata` (only `Q0norm`, `z_min` via PRODUCTION).
- Threads: RAII `SingleThread` sets `omp_set_num_threads(1)` and `Eigen::setNbThreads(1)`, restores both in the destructor (also on exception). The physics headers reached by the oracle contain no `#pragma omp`; the only parallel loops are in fitter functions not called here. Verified numerically (cpu/wall, ICV restored).
- Projector cache: `g_projector_cache_v30q.clear()` under its mutex at the start of every call; verified by the repeat-call timing.
- Invariants: (1) $Q=L^{-6}V^\dagger(F_3^{-1}+K_3)V$, inverse projected, matrix scaled by $L^{-6}$ (`np` line in `oracle.py`); (2) `slogdet` LU, $\det$ not rescaled, no $N_{\rm proj}$ factor; (3) all `PRODUCTION` keys in `metadata` (test passes). Anchors reproduced bit-for-bit with $\vec c=0$ (see below), so the $L^{-6}$ placement matches the anchor producer.

### Check 6 table (16 cases; full data in the JSON)
$\Delta_{\rm np}$ = NumPy minus anchor; $\Delta_{\rm Eig}$ = Eigen LU minus anchor.

| point | $c$ | $N$ / $N_{\rm proj}$ | $\kappa(Q)$ | $\Delta_{\rm np}$ | $\Delta_{\rm Eig}$ | $\kappa\varepsilon$ |
|---|---|---|---|---|---|---|
| 100_A2 0.2915561 | 0 | 56/17 | 9.0e6 | 0 | 0 | 1.0e-9 |
| 100_A2 0.2915561 | base | 56/17 | **3.8e12** | **-1.4e-10** | **-1.2e-9** | 4.1e-4 |
| 000_A1m 0.26932 | base | 35/5 | 3.6e6 | -1.1e-13 | -8.5e-14 | 3.9e-10 |
| 200_A2 0.28664 | base | 47/16 | 8.9e4 | 0 | 4.3e-14 | 9.8e-12 |
| 200_A2 0.34738 | 0 | 121/32 | 4.8e8 | 0 | **5.2e-10** | 5.2e-8 |
| 111_A2 0.32036 | 0 | 64/20 | 1.0e7 | 0 | -1.7e-13 | 1.1e-9 |
| 100_A2 0.36 | 0 | 108/30 | 7.6e4 | 0 | 4.9e-12 | 8.4e-12 |
| (7 further zero-$c$ cases) | 0 | | 3.6e4-5.6e5 | 0 | $\le 8\times10^{-13}$ | |

Signs agree (NumPy, Eigen, anchor, `eigvalsh` product) in all 16 cases.

### Check 7: single-thread timings (`PointResult.timings`, $\vec c=\vec c_{\rm base}$, seconds)
| point | $N$/$N_{\rm proj}$ | total | $F$ | projector | projector share | other |
|---|---|---|---|---|---|---|
| 100_A2 0.29156 | 56/17 | 3.14 | 3.075 | 0.034 | 1.1% | basis 0.023 |
| 100_A2 0.29180 | 56/17 | 2.68 | 2.610 | 0.035 | 1.3% | |
| 000_A1m 0.26932 | 35/5 | 2.03 | 1.759 | 0.244 | 12.0% | |
| 000_A1m 0.34538 | 103/9 | **10.42** | 5.077 | **5.292** | **50.8%** | |
| 100_A2 0.3599 | 108/30 | 6.05 | 5.727 | 0.249 | 4.1% | |
| 100_A2 0.36 | 108/30 | 6.01 | 5.724 | 0.224 | 3.7% | |
| 110_A2 0.30733 | 68/27 | 3.63 | 3.568 | 0.021 | 0.6% | |
| 110_A2 0.35 | 108/41 | 5.75 | 5.620 | 0.075 | 1.3% | |
| 111_A2 0.32036 | 64/20 | 2.56 | 2.498 | 0.031 | 1.2% | |
| 111_A2 0.35103 | 160/40 | 8.95 | 8.402 | 0.437 | 4.9% | |
| 200_A2 0.28664 | 47/16 | 2.12 | 2.071 | 0.022 | 1.0% | |
| 200_A2 0.34738 | 121/32 | 6.71 | 6.336 | 0.294 | 4.4% | |

$F$ is 87-98% of the time except at 000_A1m 0.34538, where the projector (5.3 s, 51%) matters as much as $F$. $G$, $K_2^{-1}$, $F_3$, $F_3^{-1}$, $K_3$, projection and LU together are $\le 0.12$ s everywhere. Range 2.0-10.4 s per point (research 007 said 2.8-6.5 s; the 000_A1m 0.34538 outlier is 60% above its upper end, driven by the projector). A second run with `OMP_NUM_THREADS=OPENBLAS_NUM_THREADS=4` gave the same totals within ~5%.

## Notes (not failures; for the auditor / planner)
1. **Anchor tolerances sit below the LU round-off floor of the pure-Eigen path.** With Eigen's `PartialPivLU` for the last step (the "unchanged" code), the anchors would fail: 200_A2 0.34738 gives $5.2\times10^{-10}>10^{-10}$ (slow test), and 100_A2 0.2915561 at $c_{\rm base}$ gives $1.2\times10^{-9}>10^{-9}$. Prompt §4 justified the tolerance by "Eigen vs LAPACK ~1e-13 measured"; that holds only for $\kappa(Q)\lesssim10^{7}$. The builder's choice (final $Q$ and $\ln\lvert\det\rvert$ in NumPy, the same call that made the anchors) is what makes the tests pass: $\Delta_{\rm np}=0$ exactly at all 12 zero-$c$ anchors. It is legitimate (invariants 1-2 respected, Eigen values kept in metadata, Q agrees to $\le 3\times10^{-16}$) but means these tests check bit-reproducibility of the same LAPACK, not the physical value: at 100_A2 0.2915561/base the independent `eigvalsh` product differs from the LU value by $4\times10^{-7}$ (true uncertainty $\lesssim\kappa\varepsilon=4\times10^{-4}$). The result is stable against BLAS/OMP thread count (check 6d) but may shift on a different CPU/BLAS kernel (this machine: OpenBLAS Haswell kernels). Also note the 7x margin on the baseline test. Tolerance changes are the planner's call; I changed nothing.
2. `test_oracle_reality_hermiticity_isometry` compares `slogdet(r.Q)` with `r.lnabsdet`, which is the same NumPy call, so that last assertion is tautological. Harmless.
3. `metadata["git_commit"]` is the commit at build time (`59b07e1`), not a dirty-tree marker. The `.so` must be rebuilt after each commit if provenance is wanted.
4. `K3_basis` is returned whenever it was computed, regardless of `want`. `want` excludes `K3` only when `c=0`.

## Failures
None. PENDING (expected `NotImplementedError`, later steps): 17 tests (lattice 8, qc_identities 4, fast 3, portability 2) plus `test_K2inv_pwave_both_forms_vs_reference` (`ValueError`, step 3). The 6 slow tests are not counted in the 18; of these, `test_oracle_anchors_all_points` and `test_ingredients_vs_reference_large_points` were run and pass.

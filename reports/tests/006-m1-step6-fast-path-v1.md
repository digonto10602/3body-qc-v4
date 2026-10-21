# Tester report 006: M1 step 6, fast path v1 (rungs 1-4)

- **Date:** 2026-10-21
- **Git commit:** 269aba6 (HEAD when `_core` was built; the only change since 812f9b9 "WIP (under test): M1 step 6" is `graphify-out/`, so all source is that of 812f9b9)
- **Prompt in force:** `prompts/001-m1-cache-free-single-energy-qc.md` (§3 row 6, §2.3 rungs 1-4, §4 tolerances; METHODS §3 truncation rule, §15). Acceptance tag `acceptance-v3`.
- **Mode:** step
- **Related:** `reports/tests/002-m1-step2-oracle-wrapper.md`, `003-m1-step3-reference-pwave.md`, `005-m1-step5-portability.md`, `reports/research/007-single-energy-qc-methods-and-reference.md`
- **Raw numbers (all in `results/`, each with commit, parameters, seed, machine where applicable):**
  `m1_step6_oracle_unchanged.json`, `m1_step6_truncation_check.json`, `m1_step6_gamma_stress.json`, `m1_step6_trunc_eps_scan.json`, `m1_step6_oracle_imag_parts.json`, `m1_step6_cachekey_check.json`, `m1_step6_speedtest_protocol.json`, `m1_step6_timings_1thread.json`
- **Scripts (scratchpad, not in repo):** `/tmp/claude-0/-home-user-3body-qc-v4/756827e6-de12-5b97-9946-026a972abb7d/scratchpad/{dump_oracle,cmp_oracle,trunc_test,gamma_stress,cachekey,threads,speedtest,timing,strace_run}.py`
- **Machine:** Intel Xeon @ 2.10 GHz, 4 cores, Linux 6.18, no GPU. Timing runs: `OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 taskset -c 0`.

## Verdict: PASS (5 non-blocking notes, none touches the physics)

Counts: acceptance m1 **33 passed / 0 failed / 0 pending** (incl. the 6 slow); unit tests **27 passed / 0 failed**; `acceptance-check` PASS (matches acceptance-v3). New failures vs. report 005: none.

## Commands run
```
graphify update . ; graphify affected fast_point ; graphify query "fast_point truncation projector cache"
rm src/python/qc3/_core*.so; ./build.sh core                       # 97 s, _core 0.3.0-step6, git_commit 269aba6
git diff 48d1319 -- src/qc3_core/oracle_point.hpp                  # reviewed (see check 2)
git diff --quiet 8b3e7d8 HEAD -- src/core src/fitter src/cachegen_cpu   # exit 0: byte-unchanged
python3 -m pytest -q tests/acceptance/m1/test_m1_oracle.py -p no:cacheprovider          # 7 passed (71 s)
python3 dump_oracle.py <old build 48d1319> old.pkl ; python3 dump_oracle.py <HEAD build> new.pkl ; python3 cmp_oracle.py
python3 -m pytest -q tests/acceptance/m1/test_m1_fast_vs_oracle.py -p no:cacheprovider --durations=0   # 7 passed
python3 -m pytest -q tests/acceptance/m1 -p no:cacheprovider                            # 33 passed (229 s)
python3 -m pytest -q tests --ignore=tests/acceptance -p no:cacheprovider                # 27 passed (19 s)
.claude/bin/acceptance-check                                                            # PASS
trunc_test.py, gamma_stress.py, cachekey.py, threads.py, strace_run.py, speedtest.py, timing.py   (see Scripts)
```
(The old oracle was a throw-away `git worktree` of 48d1319, built with `./build.sh core`, and removed afterwards.)

## Results

| # | check | expected | got | tolerance | result |
|---|---|---|---|---|---|
| 1 | clean rebuild | builds | `build.sh core` OK, 97 s, `_core` 0.3.0-step6 | - | PASS |
| 2a | `test_m1_oracle.py` incl. slow | 7/7 | 7 passed (71 s; `anchors_all_points` 48.7 s) | - | PASS |
| 2b | refactored oracle vs oracle built from 48d1319: $F,G,\mathcal K_2^{-1},F_3,F_3^{-1},V_{\rm sel},\mathcal K_3$ (4 basis + sum), $Q$, basis, sign, $\ln\lvert\det Q\rvert$ at 12 anchors x {$\vec c=0$, baseline} = 24 points | bit-identical | `np.array_equal` True for every array at all 24 points; $\lvert\Delta\ln\lvert\det Q\rvert\rvert=0.0$ | exact | PASS |
| 2c | refactored oracle vs `oracle_anchors_L20.json`, 12 anchors, $\vec c=0$ | $\Delta\ln\lvert\det Q\rvert=0$ | $\lvert\Delta\rvert=0.0$ at all 12, signs equal | exact | PASS |
| 2d | `git diff 48d1319 oracle_point.hpp` is a pure refactor | no formula edit | Reviewed: (i) `Kin`/`make_kinematics` = the removed block moved verbatim (same `FitSettings` fields, same `make_base_physics`, same `config_maker_4_momentum_first` argument lists); the cache clear still precedes the build; (ii) `compute_k3` = the removed K3 loop (same `K3mat_2plus1` args: `En_c`, `plm/klm_config`, `total_P`, `atmK`, `atmpi`, `Kiso`, `KB`, `KE`, `'n'`; same `K3 += c_i B_i` order); (iii) added `stats` map only. Confirmed by 2b | - | PASS |
| 2e | `src/core`, `src/fitter`, `src/cachegen_cpu` byte-unchanged since 8b3e7d8 | no diff | `git diff --quiet` exit 0 | exact | PASS |
| 3 | `test_m1_fast_vs_oracle.py`, all 7 | 7 green | 7 passed in 171 s (listed below) | - | PASS |
| 4a | full `tests/acceptance/m1` | 33 green | 33 passed (229 s), 0 pending | - | PASS |
| 4b | unit tests | green | 27 passed | - | PASS |
| 4c | `acceptance-check` | PASS | "PASS: tests/acceptance matches acceptance-v3" | - | PASS |
| 5a | truncation: $F$(fast, $\epsilon=10^{-16}$) vs $F$(fast, $10^{-40}$), 35 points (5 seeded random + $E=0.2631$ + $E=0.36$ per irrep; seed 20261021) | $\le10^{-14}$ | max $1.07\times10^{-16}$ (110_A2, 0.28258); typical $10^{-18}$ | $10^{-14}$ | PASS |
| 5a | same points, $F$(fast) vs $F$(oracle) rel. Frobenius | $\le10^{-12}$ | max $7.34\times10^{-15}$ (200_A2, 0.36); `fast40` vs oracle identical to 3 digits, so this is **not** truncation (guess: summation order/cancellation in $B-C$) | $10^{-12}$ (acceptance) | PASS (margin 136) |
| 5a | worst single spectator-momentum block, rel. | informational | $3.5\times10^{-14}$ (111_A2, 0.31854) | - | info |
| 5a | $\gamma$ stress: $E$ chosen to maximise $\gamma_{\max}$ over the basis, per irrep (000_A1m 5.8; 100_A2 7.4; 110_A2 8.1; 111_A2 8.1; 200_A2 9.98) | $\le10^{-14}$ / $\le10^{-12}$ | trunc: $\le4.8\times10^{-18}$; vs oracle $\le8.5\times10^{-16}$ (block $\le6.5\times10^{-15}$) | as above | PASS |
| 5a | $\epsilon$-scan, 111_A2 0.32036, vs $10^{-40}$ | error $\propto\epsilon$ | rel. err: $\epsilon=10^{-4}$: $7.4\times10^{-6}$; $10^{-8}$: $6.1\times10^{-10}$; $10^{-12}$: $5.1\times10^{-14}$; $10^{-16}$: $3.3\times10^{-18}$; also $\epsilon=0.99,0.5$ run without error | - | PASS (error $\approx0.01$-$0.07\,\epsilon$) |
| 5b | oracle $F,G,\mathcal K_2^{-1},F_3,F_3^{-1},V_{\rm sel},Q,\mathcal K_3$: $\max\lvert\mathrm{Im}\rvert/\max\lvert\cdot\rvert$ over 24 points | $\approx0$ | exactly $0.0$ for all of them, incl. the four $\mathcal K_3$ basis matrices and baseline $\mathcal K_3$; also $F$ at all 35 points of 5a | $\le10^{-15}$ | PASS |
| 5b | $\lVert F_3^{-1}-F_3^{-1\dagger}\rVert/\lVert F_3^{-1}\rVert$ (oracle) | $\le10^{-14}$ | max $2.5\times10^{-15}$ | $10^{-14}$ | PASS |
| 5c | no file I/O in the fast path | none | `strace -f -e trace=file,read,write,openat,open,pread64,mmap,newfstatat` between markers, cold call + 2 warm calls (110_A2, 111_A2): **0** syscalls. `grep` for `fopen/ifstream/ofstream/cout` in `fast_point.hpp`, `oracle_point.hpp`, `module.cpp`: none | 0 | PASS |
| 5c | projector cache keyed correctly: 85 calls shuffled across 5 irreps, 45 energies, 14 values of $N$ shared between irreps, no clearing (48 hits) vs a cold build per point | identical | $V_{\rm sel}$, $\ln\lvert\det Q\rvert$, sign, $N$, $N_{\rm proj}$ all bit-equal: 0 mismatches | exact | PASS |
| 5c | key completeness | all inputs in key | **fails for `eig_tol`**: 100_A2 0.2915561 after a default call, `eig_tol=1e-30` gives a cache hit with $N_{\rm proj}=17$ (cold build gives 5) | - | NOTE 2 (not production) |
| 5d | invariants 1-3 | | Code read: $F_3^{-1}$ is the explicit LU inverse of $F_3/s$ then $/s$ (same as the oracle); $Q=V^\dagger(F_3^{-1}+\mathcal K_3)V/L^6$ with $L=\xi L/a_t=68.88$ in `qc.project_Q` (matrix scaled, no $N_{\rm proj}$ factor, shared with the oracle); $\ln\lvert\det\rvert$ from LAPACK LU `slogdet`; metadata contains every invariant-3 constant (printed: $\xi,a_tm_K,a_tm_\pi,\eta,\alpha,\epsilon_H$, `max_shell_num`, waves, parity, `eig_tol`, scatter, `Q0norm`, `k2_pwave_form`, `git_commit`, `trunc_eps`) | - | PASS |
| 5e | single thread | 1 thread | C++ part: `SingleThread` guard forces `omp_set_num_threads(1)` and `Eigen::setNbThreads(1)`; with `OPENBLAS_NUM_THREADS=1` process CPU/wall $=0.96$-$0.99$. **Without it the NumPy tail (`project_Q`, `slogdet`) runs on OpenBLAS threads: CPU/wall $=3.9$** (wall unchanged, 4 cores kept busy) | - | NOTE 1 |
| 6 | timings, see below | informational | see below | - | info |

### Check 3: step-6 test ids (pinned to one core, `OMP_NUM_THREADS=1`)

| test | time |
|---|---|
| `test_fast_vs_oracle_small_points` | 18.84 s |
| `test_fast_basis_identical_to_oracle` | < 5 ms (cached in the process) |
| `test_fast_metadata_and_timings` | < 5 ms |
| `test_fast_vs_oracle_all_anchor_points` (slow) | 47.75 s |
| `test_fast_vs_oracle_random_energies_all_irreps` (slow) | 38.70 s |
| `test_fast_root_vs_oracle_100_A2_windows` (slow) | 49.88 s |
| `test_fast_speedup_single_thread_ge_20x` (slow) | 16.09 s |

**Speed-up measured with the acceptance protocol** (110_A2, $E=0.35$, $N=108$, median of 3, 1 thread, pinned; `speedtest.py`): fast $59.6$ ms (59.7, 58.7, 60.3), oracle $5.525$ s; ratio **$92.7\times$** (threshold 20). Excluding the oracle's own cold projector build (0.07 s here) the ratio is $91.4\times$.

### Check 6: 12 anchors, $\vec c$=baseline, 1 thread, pinned (`results/m1_step6_timings_1thread.json`)

Medians: fast warm = 5 repeats, oracle = 3 repeats ($N<100$) or 2.

| point | $N$ | $N_{\rm proj}$ | fast cold | fast warm | oracle (of which projector) | speed-up warm | excl. oracle projector | vs fast cold | $\lvert\Delta\ln\lvert\det\rvert\rvert$ |
|---|---|---|---|---|---|---|---|---|---|
| 100_A2 0.29156 | 56 | 17 | 0.089 s | 0.033 s | 2.68 s (0.06) | 82x | 81x | 30x | $1.9\times10^{-9}$ |
| 100_A2 0.29180 | 56 | 17 | 0.113 | 0.034 | 2.56 (0.03) | 75x | 74x | 23x | 0 |
| 000_A1m 0.26932 | 35 | 5 | 0.255 | 0.026 | 2.07 (0.23) | 80x | 71x | 8x | $7.8\times10^{-14}$ |
| 000_A1m 0.34538 | 103 | 9 | 5.481 | 0.054 | 10.16 (5.26) | 189x | 91x | 1.9x | $2.7\times10^{-12}$ |
| 100_A2 0.35990 | 108 | 30 | 0.268 | 0.059 | 5.57 (0.21) | 94x | 91x | 21x | $1.4\times10^{-13}$ |
| 100_A2 0.36000 | 108 | 30 | 0.262 | 0.057 | 5.40 (0.20) | 95x | 91x | 21x | $5.7\times10^{-12}$ |
| 110_A2 0.30733 | 68 | 27 | 0.057 | 0.038 | 3.11 (0.02) | 82x | 82x | 55x | $8.5\times10^{-14}$ |
| 110_A2 0.35000 | 108 | 41 | 0.134 | 0.059 | 5.35 (0.07) | 91x | 90x | 40x | $8.5\times10^{-14}$ |
| 111_A2 0.32036 | 64 | 20 | 0.064 | 0.036 | 2.44 (0.03) | 68x | 67x | 38x | 0 |
| 111_A2 0.35103 | 160 | 40 | 0.525 | 0.103 | 8.41 (0.42) | 82x | 78x | 16x | $2.0\times10^{-13}$ |
| 200_A2 0.28664 | 47 | 16 | 0.050 | 0.031 | 1.97 (0.04) | 63x | 62x | 40x | 0 |
| 200_A2 0.34738 | 121 | 32 | 0.354 | 0.068 | 6.35 (0.28) | 94x | 90x | 18x | $5.1\times10^{-10}$ |

Speed-up warm: min $63\times$, median $82\times$, max $189\times$ (the max is inflated by the oracle's 5.3 s projector build); excluding the oracle projector: min $62\times$, median $81\times$, max $91\times$. Fast cold vs oracle: median $22\times$, min $1.9\times$ (000_A1m 0.34538, see Note 4). All $\lvert\Delta\ln\lvert\det\rvert\rvert\le5.1\times10^{-10}$ (acceptance $10^{-6}$); the largest ones are at the ill-conditioned points already noted in the STATUS.

**Warm stage breakdown (ms, one call; `fast.timings`)**:

| point | basis | F | K2inv | G | F3 | proj. | F3inv | K3 | proj. + LU | total |
|---|---|---|---|---|---|---|---|---|---|---|
| 100_A2 0.29156 ($N$=56) | 21.0 | 3.4 | 0.04 | 2.3 | 0.14 | 0.03 | 0.08 | 4.9 | 0.13 | 32.0 |
| 000_A1m 0.34538 (103) | 21.1 | 5.6 | 0.12 | 7.8 | 0.67 | 0.08 | 0.44 | 16.7 | 0.25 | 52.7 |
| 110_A2 0.35000 (108) | 21.3 | 6.5 | 0.11 | 8.6 | 0.78 | 0.08 | 0.45 | 18.6 | 0.43 | 56.9 |
| 111_A2 0.35103 (160) | 21.3 | 10.6 | 0.28 | 19.6 | 2.30 | 0.12 | 1.48 | 42.7 | 0.63 | 99.5 |

Oracle at 111_A2 0.35103: $F$ 7.83 s of 8.36 s (94%); fast $F$ 10.6 ms, i.e. **$F$ itself is faster by about $740\times$**. In the fast path, $F$ is now only 11-17% of the time. Lattice terms kept are 5.0-8.7% of the oracle's cube.

## Notes (non-blocking; none changes a physics number)

1. **NumPy tail is not single-threaded by itself.** `fast.py`/`oracle.py` call `project_Q` and `slogdet` through NumPy/OpenBLAS. In the default environment the process burns about 4 cores (CPU/wall $=3.9$, measured idle) although wall time is the same; under load the wall time doubled (0.44 s vs 0.18 s at 110_A2). The "1 core" claim (D6) holds only with `OPENBLAS_NUM_THREADS=1` (or `OMP_NUM_THREADS=1`) set before NumPy is imported, or with `taskset`. The benchmark protocol should state the environment. (guess: a `threadpoolctl.threadpool_limits(1)` around the tail, or doing the projection in C++, fixes it.)
2. **Projector cache key omits `eig_tol` (and `norm_tol`, `proj_tol`).** Stale hit demonstrated (17 vs 5 columns at `eig_tol=1e-30`). Production uses the fixed 0.05 and the oracle clears the cache, so nothing in M1 is affected; it matters if Sector overrides are ever mixed in one process. The key does include irrep, parity, $\vec P$ and the full $(\vec n,\ell,m)$ list of both channels, so different-basis collisions did not occur in 85 shuffled calls.
3. **Rung 4 is only partly implemented**: $G$ and $\mathcal K_3$ use the oracle functions with no symmetric fill, and the projection is the explicit inverse plus NumPy products (not an LU solve with $N_{\rm proj}$ right-hand sides). This is allowed by invariant 1 and the acceptance tests pass, but the stage table shows that after rungs 1-3 the time is dominated by the **shared basis builder (21 ms, independent of $N$)**, $\mathcal K_3$ (5-43 ms) and $G$ (2-20 ms). The step-7 stop criterion (<= 20 ms at $N\approx160$) cannot be reached by class tables for $F$ alone: warm total is 99.5 ms at 111_A2 0.35103.
4. **A cold projector build is not "ms-10 ms"** (plan §2.3 rung 4): it is 5.4 s at 000_A1m 0.34538 ($N$=103) and 0.2 s at $N$=108 in 100_A2, so each new basis size in a root search pays it once (fast cold is only 1.9x the oracle there). The projector depends only on the basis and is a candidate for the later speed loop.
5. **Wording of the truncation bound.** The bound $3R^2e^{-\alpha(R^2-x^2)}/(R^2-x^2)\le\epsilon$ is valid per omitted term (re-derived: $\lvert Y_aY_b\rvert\le3r^2$ for $\ell\le1$ needs only $r^2\ge1/3$; the logarithmic derivative $1/r^2-1/(r^2-x^2)-\alpha<0$ for $x^2\ge0$ and $\le\lvert x^2\rvert/(r^2(r^2+\lvert x^2\rvert))-\alpha<0$ for $x^2<0$ and $r^2>R^2\ge x^2+1/\alpha$, so the tail is decreasing; the bounding box and $\pm1$ margins are correct; the ellipsoid is clipped by the same $[-20,20]^3$ cube as the oracle). The **total** dropped sum is not $\le\epsilon$: the ellipsoid volume gives $\approx2\pi\gamma R/\alpha\cdot\epsilon\approx10^{3}\epsilon\approx10^{-13}$ at $\gamma=10$ as a worst case using the crude $3r^2$ bound. Measured effect is 3-4 orders of magnitude smaller ($\le1.1\times10^{-16}$), because real $\lvert Y_aY_b\rvert$ is far below $3r^2$, so the stated $10^{-14}$ holds with margin 90 and the acceptance $10^{-12}$ with margin $\approx10^4$. Hence the unit test's $10^{-14}$ is justified empirically, not by the per-term bound alone.

## Failures

None.

## Verdict line
`VERDICT: PASS`. 33/33 acceptance m1 passed (0 failed, 0 pending), 27/27 unit passed, `acceptance-check` PASS, oracle bit-identical to the pre-refactor build at 24 points, fast path 63-95x (excluding the oracle projector) at the 12 anchors with warm cache, 92.7x in the acceptance-protocol test. The formal benchmark report (§7) is not part of this run.

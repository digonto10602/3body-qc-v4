# Test report 009: M1 step 7.2 (orbit/character projector) and 7.3 (threshold basis builder, `dimension_jumps`)

- **Date:** 2026-10-22
- **Git commit:** `e53627d` (HEAD). The code under test is `9b6cef5`; `git diff 9b6cef5 HEAD -- src tests` is empty (only `.gitignore` and graphify files differ). `_core` was rebuilt from scratch in this tree (1 m 49 s); its `metadata["git_commit"]` reads `46c6111`, the HEAD at build time, whose `src/` is the same.
- **Prompt in force:** `prompts/002-m1-step7-speed-replan-projector-basis-k3-g.md` (§2.1, §2.2, §3 rows 7.2/7.3, §4, §5, §7). Acceptance tag `acceptance-v4` (47 tests).
- **Mode:** step.
- **Related:** `reports/tests/008-m1-step7-0-1-8-cache-roots-default.md`, `reports/tests/006-m1-step6-fast-path-v1.md`, `reports/research/008-step7-speed-replan-orbit-projector-and-kernels.md`.
- **Verdict: PASS for 7.2 and 7.3.** 45 passed, 2 failed, 0 new failures. Both failures are the tests that wait for 7.4–7.6 (the cause was checked, see "Classification"). Compared with report 008: 40 -> 45 passed; the 7 earlier `NotImplementedError` tests are now 5 green and 2 pending-on-7.6.

## Commands run (repo root; scripts are copied to `results/tester_009_scripts/`)

```
graphify update . ; graphify affected "build_orbit_V" / "fast_projector_cached" / "fast_kinematics" / "dimension_jumps"
mv src/python/qc3/_core*.so <scratch> ; ./build.sh core                                   # clean rebuild
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 python3 -m pytest -q tests/acceptance/m1 -rA --durations=0      # 339 s
env -u OMP_NUM_THREADS -u OPENBLAS_NUM_THREADS python3 -m pytest -q tests/acceptance/m1 -rA              # 400 s (unpinned)
python3 -m pytest -q tests --ignore=tests/acceptance ; .claude/bin/acceptance-check ; pytest --collect-only -q tests/acceptance
git diff 8f8e9d3 -- src/qc3_core/oracle_point.hpp src/core src/fitter src/cachegen_cpu
# own checks (results/tester_009_scripts): projcheck.py 12 | lndet.py | scanN.py | ulp.py | stale.py | io_probe.py (strace -f) | tailcheck.py
# C++ harness (g++ -O1, includes the production headers; mutants shadow fast_projector.hpp by -I order): harness.cpp
# timings, taskset -c 0, OMP=OPENBLAS=1: timing.py | fresh.py | coldwarm.py ; cpuwall.py (pinned and unpinned)
```
Raw numbers: `results/tester_009_{summary,projector_check,lndet,basis_scan,basis_bitident,timing,coldwarm_111A2}.json`. Seed for all random draws: 20261022. All own scripts use `coded_v2`, the pinned form of the locked tests.

## Summary table

| check | expected | got | tolerance | result |
|---|---|---|---|---|
| clean rebuild `_core` | builds | builds, 2.14 MB | n/a | PASS |
| `tests/acceptance/m1` pinned (47) | 45 pass, 2 pending on 7.4-7.6 | 45 passed, 2 failed (below) | exact | PASS |
| `tests/acceptance/m1` unpinned | same set | 45 passed, same 2 failed | exact | PASS |
| unit tests | all green | 44 passed (23 s) | exact | PASS |
| `acceptance-check` | PASS | "PASS: tests/acceptance matches acceptance-v4"; 47 collected | exact | PASS |
| oracle untouched | empty diff | 0 lines (`oracle_point.hpp`, `src/core`, `src/fitter`, `src/cachegen_cpu` vs `8f8e9d3`) | exact | PASS |
| 5 tests of this step | green | `test_projector_orbit_matches_eigen_and_oracle_subspace`, `test_projector_cold_build_is_fast`, `test_dimension_jumps_match_rule_and_N`, `test_fast_basis_all_signatures_vs_reference` (slow), `test_fast_single_thread_cpu_time` all PASSED | locked | PASS |
| 7.2 $N_{\rm proj}$: my dense $P$ vs orbit vs eigen, **58 bases** (all) | equal | equal at all 58; oracle also equal at 14 bases (incl. 111_A2 $N$=160, 000_A1m $N$=135 and $N$=103, 200_A2 $N$=137) | exact | PASS |
| 7.2 $\lVert V_oV_o^\dagger-P\rVert/\lVert P\rVert$, 58 bases | $\ll10^{-12}$ | max orbit $4.0\times10^{-16}$; eigen $6.9\times10^{-16}$; oracle (14 bases) $6.9\times10^{-16}$ | $10^{-12}$ | PASS |
| 7.2 own $P$ is a projector | $P^2=P$, $P=P^T$ | $\lVert P^2-P\rVert/\lVert P\rVert\le1.4\times10^{-16}$, $P_{\rm herm}=0$ (58 bases) | $10^{-12}$ | PASS |
| 7.2 orthonormality $\lVert V^\dagger V-1\rVert_F$ | $\le10^{-12}$ | orbit max $2.2\times10^{-15}$; eigen max $1.6\times10^{-15}$; $V$ exactly real | $10^{-12}$ | PASS |
| 7.2 sign of $\det Q$, 10 points | equal | equal for orbit, eigen, oracle and my eigenvalue path at all 10 | exact | PASS |
| 7.2 $\lvert\Delta\ln\lvert\det Q\rvert\rvert$ orbit vs eigen | $\le10^{-6}$ | max $1.5\times10^{-9}$ (100_A2 0.2915561, $\kappa(Q)=2.6\times10^{12}$) | $10^{-6}$ | PASS |
| same, orbit vs oracle / eigen vs oracle | $\le10^{-6}$ | max $1.5\times10^{-9}$ / $5.1\times10^{-10}$ | $10^{-6}$ | PASS |
| same, orbit vs my own `eigvalsh` path with my own $V$ | $\le10^{-6}$ | max $3.6\times10^{-7}$ (the same ill-conditioned point); $\le1.2\times10^{-7}$ elsewhere | $10^{-6}$ | PASS (round-off floor) |
| 7.2 other irreps, orbit vs existing eigen builder (harness) | equal | 25 (group, irrep) cases incl. 2-dim ($E,E_g,E_u$) and 3-dim ($T_{1g},T_{1u},T_{2g},T_{2u}$): $N_{\rm proj}$ equal, max rel$(VV^\dagger)=2.4\times10^{-15}$ | $10^{-12}$ | PASS |
| 7.2 self-checks always on, throw | throw, no fallback | read: no `try/catch` in `src/qc3_core/fast_*.hpp`, `module.cpp`, `fast.py`. Harness: dropped row $\Rightarrow$ throw 4/4; mutant $\times2$ scale $\Rightarrow$ throw 5/5; mutant truncated candidates $\Rightarrow$ throw 4/5 | n/a | PASS |
| 7.3 number of bases (oracle's own basis by bisection on the 20000-point grid) | 5/12/16/13/12 | 5/12/16/13/12 distinct $N$ (N values listed in `tester_009_basis_scan.json`) | exact | PASS |
| 7.3 `dimension_jumps` count, no missed / invented jump | 4/11/15/12/11 | 4/11/15/12/11; each of the 53 grid cells with a change contains exactly one jump, and no jump lies outside such a cell | exact | PASS |
| 7.3 jump position, oracle bisection to double precision vs `dimension_jumps` | grid resolution | $\le1.1\times10^{-16}$ absolute ($4.8\times10^{-6}$ is the grid step) | $4.8\times10^{-6}$ | PASS |
| 7.3 jumps vs my NumPy threshold formula | $10^{-12}$ rel | max relative difference $0.0$ (bitwise equal) at all 53 | $10^{-12}$ | PASS |
| 7.3 basis bit-identical to `config_maker_4_momentum_first` | identical | 1881 energies (300 random per irrep, both ends, and $\pm1,2,3$ ulp around every jump), 0 mismatches of momenta (complex doubles), $n$, $\ell$, $m$, order, $N_1$, $N_2$; the flip of $N$ sits within 2 ulp of the `dimension_jumps` value | bitwise | PASS |
| 7.3 threshold-table key complete | no stale hit | 132 cases varying $L/a$, $\xi$, $a_tm_K$, $a_tm_\pi$, $\epsilon_H$, `max_shell_num`, $\alpha$, $\eta$ with a warm cache: 0 mismatches | bitwise | PASS |
| 7.3 no disk I/O | 0 syscalls | `strace -f` of cold `dimension_jumps` + 3 evaluations: 0 `openat/read/write/stat/...` between markers; `cache_info` = 4 threshold entries (2 sectors x 2 channels) | exact | PASS |
| timing: cold projector, all 58 bases (5 repeats each) | $<20$ ms (test); plan $<1$ ms | max 1.20 ms (111_A2 $N$=148); median 0.33 ms; 000_A1m 0.34538: 0.63-0.94 ms (was 5.4 s) | 20 ms | PASS |
| timing: cold $-$ warm total, 12 anchors, wall | $\le0.1$ s | 6.7-29.4 ms (one jitter case $-8.2$ ms) | 0.1 s | PASS |

### Failures (both pending on later sub-steps)

| test | got | classification |
|---|---|---|
| `test_fast_tail_consistent_with_numpy_and_want_subset` | `assert 'numpy' == 'cpp'` on `metadata["tail"]` at the first point (100_A2 0.291556) | **pending 7.6.** My copy of the remaining assertions (`tailcheck.py`, 4 points, orbit projector) all hold: sign equal to `qc3.qc.slogdet` of $Q$ (difference $0$ because the tail still is that NumPy code), `F3inv is None` for `want=("Q",)`, subset equals full ($0$). |
| `test_fast_speedup_single_thread_largest_point_ge_150x` | pinned: 96.8x (oracle 7.51 s, 7.00 s without its projector; fast median 72.3 ms); unpinned: 44.2x (fast 168 ms) | **pending 7.4-7.6.** Warm stages at 111_A2 0.35103: $\mathcal K_3$ 36.5 ms, $G$ 17.4, $F$ 9.1 (63 of 68 ms); basis 0.40 and projector 0.11 are negligible. The 150x bound needs a fast time $\le47$ ms. |

## Classification of the three tests the caller flagged

- `test_fast_single_thread_cpu_time`: **passes both pinned and unpinned** (the expectation "passes only with BLAS pinned" is **not** confirmed). The test uses 000_A1m 0.26932 ($N$=35, $N_{\rm proj}$=5), where OpenBLAS never threads (CPU/wall = 0.996 unpinned, 0.983 pinned). So the locked test cannot detect the tail's thread use. The thread use is real and is the NumPy tail (`metadata["tail"] == "numpy"`): CPU/wall over 20 warm calls, unpinned: 111_A2 0.35103: **3.06**; 000_A1m 0.34538: **2.29**; pinned 0.99. Setting only `OPENBLAS_NUM_THREADS=1` gives 0.97-1.00, so no other thread source exists. This belongs to 7.6; its done-check ("CPU/wall $\le1.1$ without pinning") must be measured at a large point, not at the test's point.
- `test_fast_tail_consistent...` and `test_fast_speedup...`: as in the table above.

## 7.2 physics review

Code: `build_orbit_V` in `src/qc3_core/fast_projector.hpp`. It uses the oracle's own `little_group`, `chi`, `cubic_transf`, `parity_factor_for_l(parity, R, 0, 1)` ($\Pi(R)=-1$ for improper $R$, for every $\ell$) and `D_real_element_v30m`, and builds $U(R)e_{(k,\ell,m)}=\Pi(R)\sum_{m'}D^{(\ell)}_{m'm}(R)\,e_{(Rk,\ell,m')}$, which is $U=S\,W$ of the oracle's convention $c_0$ (`order_WS=false`, no inverse, no transpose, `parity_mode=1`). $P_\Gamma=\frac{d_\Gamma}{|LG|}\sum\chi U$ is the formula of the plan. For $d_\Gamma=1$ it uses the image spanned by $P e_{(k_0,\ell,m)}$ per (flavour, orbit, $\ell$) block (rearrangement lemma); for $d_\Gamma>1$ all block rows are candidates.

My independent derivation (own NumPy, from scratch, no `chi`, no `little_group`): the 48 signed permutation matrices $M$; $LG(\mathbf P)=\{M:M\mathbf P=\mathbf P\}$; for all five production sectors the irrep is one-dimensional with $\chi_\Gamma(M)=\det M$ ($A_{1u}$ of $O_h$; $A_2$ of $C_{4v},C_{2v},C_{3v}$: $+1$ on rotations, $-1$ on reflections); three pseudoscalars give $\Pi(M)=\det M$; so $\chi\Pi=1$ and $P=\frac1{|G|}\sum_M U_{\rm orb}(M)$ with $D^{(1)}=M$ in the axes $(m=1,0,-1)\leftrightarrow(x,z,y)$ and $D^{(0)}=1$. That is the $J^P=0^-$ sector (orbital part totally symmetric under $LG$). This $P$ is a Hermitian idempotent ($\le1.4\times10^{-16}$) and equals $V_oV_o^\dagger$, $V_eV_e^\dagger$ and $V_{\rm oracle}V_{\rm oracle}^\dagger$ to $\le7\times10^{-16}$ at all 58 bases (oracle at 14), so the oracle's character table and conventions agree with the textbook group theory for these five irreps, and the orbit code reproduces them. Not checked: whether $J^P=0^-$ is the physical target of "A2" at $\mathbf P\ne0$ beyond this derivation (it is the standard pseudoscalar assignment in $C_{4v}$).

Limit of the self-check (note, not a defect): it proves that $V$ lies in, and spans, the image of the *sparse* $P$ built from the same $U$, so it catches an incomplete candidate set and a non-idempotent $P$ but cannot catch a *convention* error (for example, a wrong $\Pi$ would still give a consistent projector). Convention errors are caught only by the comparisons with the oracle/eigen builder above, which are in the locked test at 4 points and in my 58-basis check.

### $\ln\lvert\det Q\rvert$, $\vec c=\vec c_{\rm base}$, `coded_v2`, 2 energies per irrep

| irrep | $E_{\rm cm}$ | $N$ | $N_{\rm proj}$ | orbit-eigen | orbit-oracle | eigen-oracle | orbit-mine | $\kappa(Q)$ |
|---|---|---|---|---|---|---|---|---|
| 000_A1m | 0.30000 | 47 | 6 | 0 | $1.1\times10^{-13}$ | $1.1\times10^{-13}$ | $1.2\times10^{-13}$ | $1.0\times10^{3}$ |
| 000_A1m | 0.34538 | 103 | 9 | $3.8\times10^{-12}$ | $6.5\times10^{-12}$ | $2.7\times10^{-12}$ | $8.7\times10^{-10}$ | $4.5\times10^{7}$ |
| 100_A2 | 0.2915561 | 56 | 17 | $1.5\times10^{-9}$ | $1.5\times10^{-9}$ | $9.5\times10^{-11}$ | $3.6\times10^{-7}$ | $2.6\times10^{12}$ |
| 100_A2 | 0.33500 | 100 | 28 | $1.1\times10^{-13}$ | $1.5\times10^{-12}$ | $1.6\times10^{-12}$ | $1.7\times10^{-11}$ | $6.1\times10^{5}$ |
| 110_A2 | 0.30733 | 68 | 27 | 0 | $8.5\times10^{-14}$ | $8.5\times10^{-14}$ | $1.1\times10^{-12}$ | $4.2\times10^{4}$ |
| 110_A2 | 0.34000 | 100 | 39 | $2.8\times10^{-14}$ | $2.8\times10^{-14}$ | 0 | $4.1\times10^{-12}$ | $3.4\times10^{4}$ |
| 111_A2 | 0.31000 | 52 | 17 | $1.4\times10^{-14}$ | $1.4\times10^{-14}$ | $2.8\times10^{-14}$ | $2.7\times10^{-13}$ | $9.5\times10^{5}$ |
| 111_A2 | 0.35103 | 160 | 40 | $8.5\times10^{-14}$ | $1.1\times10^{-13}$ | $2.0\times10^{-13}$ | $9.3\times10^{-12}$ | $7.4\times10^{5}$ |
| 200_A2 | 0.30500 | 79 | 22 | $9.9\times10^{-14}$ | $8.5\times10^{-14}$ | $1.8\times10^{-13}$ | $5.1\times10^{-13}$ | $8.3\times10^{3}$ |
| 200_A2 | 0.34738 | 121 | 32 | $8.5\times10^{-10}$ | $3.5\times10^{-10}$ | $5.1\times10^{-10}$ | $1.2\times10^{-7}$ | $5.2\times10^{8}$ |

Maxima: orbit-eigen $1.5\times10^{-9}$, orbit-oracle $1.5\times10^{-9}$, eigen-oracle $5.1\times10^{-10}$ (tolerance $10^{-6}$). The two largest values sit at the two ill-conditioned points known from STATUS (round-off floor of the LU, the same order as the $1.9\times10^{-9}$ in the plan). "mine" = $\sum\ln\lvert\lambda_i\rvert$ of $V_m^\dagger(F_3^{-1}+\mathcal K_3)V_m/s$ from the oracle's $F_3^{-1}$, $\mathcal K_3$ and my own $V_m$ (the scale $s$ from the oracle's $Q$ is $L^6$ to $4\times10^{-16}$, invariant 1).

## 7.3 review

- `make_setup` is a copy of the first half of `make_kinematics` (the oracle file is byte-unchanged). The comparison hook `compare_kinematics` confirmed the Kin is bit-identical (above). The comparison `E_lab > t[n]` is the oracle's own double comparison with the same function and arguments, so identity holds by construction and by test.
- Jump enumeration (own scan): for each irrep I bisected $N(E)$ of the **oracle's** `config_maker_4_momentum_first` on the production grid $E_i=0.26310+i\cdot0.0969/19999$, $i<20000$ (total oracle calls 201/533/715/568/529), then refined every jump in $E_{\rm cm}$ to the last double. Distinct bases 5/12/16/13/12 = 58; $\Delta N$ per jump is 1-48 (e.g. 000_A1m 12, 48, 8, 32), always positive. `dimension_jumps` dedupes thresholds equal within $10^{-12}$ relative (orbit partners); no real pair of distinct jumps is that close at L20.
- Memory and I/O: threshold tables only in memory (0 syscalls); 4 tables for the 2 sectors in the probe.

## Item 6: timings (taskset -c 0, OMP=OPENBLAS=1, $\vec c_{\rm base}$, `want=("Q",)`)

Stage times in ms (median of 7 warm calls; "cold" = after `clear_caches()`, which also drops the threshold tables):

| point | $N$ | $N_{\rm proj}$ | cold basis | cold proj | warm basis | warm proj | $F$ | $G$ | $\mathcal K_3$ | warm total | cold total $-$ warm |
|---|---|---|---|---|---|---|---|---|---|---|---|
| 100_A2 0.29156 | 56 | 17 | 8.74 | 0.22 | 0.44 | 0.036 | 2.9 | 2.1 | 4.5 | 10.7 | 8.8 |
| 100_A2 0.29180 | 56 | 17 | 8.00 | 0.19 | 0.41 | 0.038 | 3.1 | 2.0 | 5.6 | 12.0 | 6.6 |
| 000_A1m 0.26932 | 35 | 5 | 10.57 | 0.33 | 0.34 | 0.024 | 1.4 | 0.8 | 1.7 | 4.5 | 10.8 |
| 000_A1m 0.34538 | 103 | 9 | 7.60 | 0.94 | 0.57 | 0.104 | 7.4 | 11.1 | 22.4 | 44.1 | 9.9 |
| 100_A2 0.35990 | 108 | 30 | 13.33 | 0.54 | 0.52 | 0.083 | 5.5 | 7.8 | 17.1 | 32.7 | 28.6 |
| 100_A2 0.36000 | 108 | 30 | 8.08 | 0.39 | 0.55 | 0.114 | 7.5 | 12.1 | 24.7 | 47.1 | $-8.0$ (jitter) |
| 110_A2 0.30733 | 68 | 27 | 13.20 | 0.27 | 0.54 | 0.071 | 5.7 | 4.8 | 9.5 | 21.6 | 13.2 |
| 110_A2 0.35000 | 108 | 41 | 13.59 | 0.40 | 0.43 | 0.083 | 5.8 | 7.4 | 16.0 | 31.6 | 17.4 |
| 111_A2 0.32036 | 64 | 20 | 8.00 | 0.21 | 0.41 | 0.043 | 4.8 | 2.1 | 5.3 | 13.3 | 7.4 |
| 111_A2 0.35103 | 160 | 40 | 8.05 | 0.53 | 0.40 | 0.110 | 9.1 | 17.4 | 36.5 | 67.8 | 9.8 |
| 200_A2 0.28664 | 47 | 16 | 8.43 | 0.20 | 0.36 | 0.034 | 2.7 | 1.3 | 3.0 | 7.8 | 8.9 |
| 200_A2 0.34738 | 121 | 32 | 7.65 | 0.48 | 0.43 | 0.094 | 6.2 | 9.7 | 21.7 | 40.6 | 7.2 |

- Projector, all 58 bases, 5 cold repeats each: **max 1.20 ms**, median 0.33 ms (eigen builder on the same bases: up to 11.8 s). 000_A1m 0.34538: 0.63-0.94 ms (plan: $<1$ ms; was 5.4 s).
- Basis builder warm: 0.34-0.57 ms at the anchors, max 0.68 ms over the 58 bases. The plan's done-check says $<0.5$ ms; it is missed by 4 of 12 anchors (0.52, 0.55, 0.54, 0.57 ms) and not an acceptance criterion; it is below 1 % of the warm total. Cold basis (both threshold tables of a sector + basis) 7.4-16.6 ms (median per base), i.e. about 4 ms per table, not the 21 ms of the plan.
- Per-point warm totals now: 4.5-12 ms at $N\le64$, 33-68 ms at $N\ge100$ (dominated by $\mathcal K_3$ and $G$, targets of 7.4/7.5).

### The builder's "unexplained 52.6 ms cold $-$ warm" at 111_A2 0.35103 ($N$=160): not reproducible, it is timing jitter in $\mathcal K_3$/$G$, not the builder

`coldwarm.py`: 30 paired repetitions (`clear_caches()` then one cold call, then median of 3 warm calls) in one process. Median cold $-$ warm per stage, [quartiles], max:

| stage | cold $-$ warm median (ms) | quartiles | max |
|---|---|---|---|
| basis (two threshold tables) | **7.75** | [7.44, 8.53] | 18.6 |
| projector | **0.48** | [0.43, 0.54] | 0.76 |
| $K_2^{-1}$ | 0.33 | [0.30, 0.39] | 0.53 |
| $F$ | $-0.26$ | [$-0.40$, 0.12] | 5.8 |
| $G$ | 0.15 | [$-0.46$, 0.55] | 10.5 |
| $\mathcal K_3$ | $-0.19$ | [$-0.94$, 3.79] | **32.1** |
| total | **8.41** | [4.70, 13.25] | **57.5** |

So the deterministic cold cost is about 8.6 ms (7.75 + 0.48 + 0.33): the threshold tables. The 52.6 ms is a tail event: the maximum total (57.5 ms) is dominated by a $\mathcal K_3$ spike (+32 ms) and a $G$ spike (+10 ms) in the same call; those two stages do not depend on the cache state (median difference $\approx0$). Six fresh-process runs (first call ever vs warm): first $-$ warm wall 7.7-15.6 ms, so no one-off initialisation cost either. Guess (not verified): scheduling jitter on a shared machine (a second builder runs in another worktree).

## Notes (non-blocking)

1. **`test_fast_single_thread_cpu_time` is not discriminating** at its point ($N$=35); see above. If the planner agrees, a point with $N\ge100$ would make it detect the NumPy tail (measured 2.3-3.1 unpinned). I did not change it (locked).
2. The self-check cannot detect a wrong convention (above); the oracle comparison does, at 58 bases.
3. `dimension_jumps(sector, E_min, E_max)` evaluates `_params` at the interval midpoint only to fill the input structure; the thresholds do not depend on it.
4. 2-dim and 3-dim irreps work (25 cases checked vs the eigen builder) although none is in production. The only irrep that failed was `000` as `A1g` for parity $-1$: the *eigen* builder throws `projection_produced_zero_columns` (the image is empty); the orbit builder returns 0 columns without complaint. Not a production case.

## Return

`VERDICT: PASS` (7.2 and 7.3). 45 passed / 2 failed / 0 other pending / 0 new failures (the 2 are pending on 7.6 and 7.4-7.6). Unit 44 passed, `acceptance-check` PASS, oracle files untouched.

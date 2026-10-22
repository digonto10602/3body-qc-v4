# Test report 010: M1 step 7.4 ($\mathcal K_{\rm df,3}$ pair kernel), 7.5 ($G$ pair kernel), 7.6 (C++ tail)

- **Date:** 2026-10-22
- **Git commit:** `338a0e8` (HEAD). The code under test is `3e26d29` (merge of `1e5269a`, `9bc2751`, `be34c57`); `git diff 3e26d29 338a0e8 -- src tests` is empty (only graphify files differ). `_core` was rebuilt from scratch in this tree (1 m 50 s); `metadata["git_commit"]` reads `338a0e8`.
- **Prompt in force:** `prompts/002-m1-step7-speed-replan-projector-basis-k3-g.md` (§2.3, §2.4, §3 rows 7.4–7.6, §4, §5; `docs/METHODS.md` §15). Acceptance tag `acceptance-v4` (47 tests).
- **Mode:** step.
- **Related:** `reports/tests/009-m1-step7-2-3-orbit-projector-basis.md`, `reports/tests/006-m1-step6-fast-path-v1.md`.
- **Verdict: PASS** (with 3 notes the coordinator should read: N1 the locked CPU-time test fails when run alone, N2 the oracle is less accurate than the fast path at $\kappa(Q)\gtrsim10^{13}$, N3 the G fill-mode deviation from the plan). 47/47 acceptance tests pass (pinned and unpinned), 63/63 unit tests pass, `acceptance-check` PASS. No pending, no new failures. Compared with report 009: 45 -> 47 passed; the 2 tests that were pending on 7.4–7.6 are green.

## Commands run (repo root; own scripts are in `results/tester_010_scripts/`, logs in `results/tester_010_logs/`)

```
graphify update . ; graphify affected "src/qc3_core/fast_point.hpp" / "fast.py"   # module.cpp, qc3.fast, roots.py (indirect), unit + acceptance tests: all run below
mv src/python/qc3/_core*.so <scratch> ; ./build.sh core                             # clean rebuild, 1m50s
env -u OMP_NUM_THREADS -u OPENBLAS_NUM_THREADS python3 -m pytest -q tests/acceptance/m1 -rA --durations=0     # 373 s, unpinned
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 python3 -m pytest -q tests/acceptance/m1 -rA                          # 300 s, pinned
python3 -m pytest -q tests --ignore=tests/acceptance ; .claude/bin/acceptance-check ; pytest --collect-only -q tests/acceptance
git diff b6ad563 -- src/qc3_core/oracle_point.hpp src/core src/fitter src/cachegen_cpu
python3 -m pytest -q tests/acceptance/m1/test_m1_fast_v2.py::test_fast_single_thread_cpu_time   # x3 unpinned, x3 with OPENBLAS_NUM_THREADS=1
blasprobe.py ; cpu_alone.py <delay> ; cpuwall.py ; cpuwall2.py <irrep> <E> [delay]              # CPU/wall probes
compare.py (47 energies, seed 20261022) | refcmp.py (vs qc3ref) | tailcheck2.py | offprod.py | illcond.py | gmode_hp.py | gmodes.py
g++ -O1 harness_tail.cpp ; g++ -O1 harness_complexV.cpp -I<mutant dir>                          # C++ harnesses (not production)
taskset -c 0 python3 timing.py ; taskset -c 0 python3 timing2.py ; taskset -c 0 python3 coldwarm.py ; taskset -c 0 python3 gmodes.py
```
Raw numbers: `results/tester_010_{summary,compare,refcmp,tail,timing,speedup,coldwarm_111A2,gmodes,gmode_hp,illcond,offprod}.json` (commit, parameters, seed 20261022 inside; all own scripts use `coded_v2`, the form the locked tests pin, and $\vec c=\vec c_{\rm base}$).

## Summary table

| check | expected | got | tolerance | result |
|---|---|---|---|---|
| clean rebuild `_core` | builds | builds, 2.24 MB, 1 m 50 s | n/a | PASS |
| `tests/acceptance/m1` unpinned (47, slow incl.) | 47 pass | **47 passed**, 373 s | exact | PASS |
| `tests/acceptance/m1` pinned | 47 pass | **47 passed**, 300 s | exact | PASS |
| unit tests (`tests --ignore=tests/acceptance`) | all green | 63 passed, 33.5 s (incl. 19 new in `test_fast_k3_g_tail.py`; 44 before) | exact | PASS |
| `acceptance-check` | PASS | "PASS: tests/acceptance matches acceptance-v4"; 47 collected | exact | PASS |
| the 2 tests that waited for 7.6 | green | `test_fast_tail_consistent_with_numpy_and_want_subset` and `test_fast_speedup_single_thread_largest_point_ge_150x` PASSED (pinned and unpinned); measured ratio below | locked | PASS |
| `test_fast_single_thread_cpu_time` **alone**, unpinned | pass | **FAIL 3/3**: CPU/wall = 4.00, 3.85, 3.68 (limit 1.3) | 1.3 | FAIL (note N1; cause verified: not the fast path) |
| same test alone with `OPENBLAS_NUM_THREADS=1` | pass | passes 3/3 (0.12–0.14 s) | 1.3 | PASS |
| oracle untouched | empty diff | 0 lines (`oracle_point.hpp`, `src/core`, `src/fitter`, `src/cachegen_cpu` vs `b6ad563`) | exact | PASS |
| $\mathcal K_3$ basis (iso0, iso1, B, E) vs oracle, 47 energies | $\le10^{-12}$ rel. Frobenius | max $8.9\times10^{-17}$ (anchors $5.8\times10^{-17}$, random $7.1\times10^{-17}$, next to jumps $8.9\times10^{-17}$); $\mathcal K_3(\vec c)$ max $9.2\times10^{-17}$ | $10^{-12}$ | PASS |
| $\mathcal K_3$ basis vs independent `qc3ref`, 8 points ($N$ = 35…160) | $\le10^{-12}$ | max $3.5\times10^{-15}$ (iso0/iso1 at 110_A2 0.35) | $10^{-12}$ | PASS |
| $\mathcal K_3$ symmetry | symmetric | fast matrices exactly symmetric (0.0); oracle's own asymmetry: iso0, iso1 0.0, B $1.2\times10^{-17}$, E $1.3\times10^{-16}$ | n/a | PASS |
| $G$ (default, every ordered pair) vs oracle, 47 energies | $\le10^{-12}$ (builder: exactly 0) | **exactly 0.0** (max abs 0.0) at all 47, and at 4 non-production parameter sets | $10^{-12}$ | PASS |
| $G$ (`g_symmetric=True`) vs oracle | $\le10^{-12}$ | max $1.03\times10^{-13}$ (111_A2 0.337697); anchors max $3.8\times10^{-14}$; $\le$ oracle's own asymmetry ($1.46\times10^{-13}$) at every point; exactly symmetric | $10^{-12}$ | PASS |
| $G$ vs independent `qc3ref`, 8 points | $\le10^{-12}$ | max $2.4\times10^{-16}$ | $10^{-12}$ | PASS |
| $(2,2)$ block of $G$ | exactly 0 | 0 (unit test; `G` at 47 points identical to oracle) | exact | PASS |
| $N$, $N_{\rm proj}$ fast = oracle, 47 energies incl. 10 pairs straddling a jump (offsets $1.6\times10^{-9}$ … $1.7\times10^{-7}$; $N$ differs by 1–24 across each pair) | equal | equal at all 47 | exact | PASS |
| $\lvert\Delta\ln\lvert\det Q\rvert\rvert$ fast (ordered $G$) vs oracle, 12 anchors | $\le10^{-6}$ | max $3.63\times10^{-9}$ (100_A2 0.2915561, $\kappa(Q)=3.8\times10^{12}$); others $\le3.5\times10^{-10}$ | $10^{-6}$ | PASS |
| same, 35 random/near-jump energies | $\le10^{-6}$ | 34 of 35 $\le1.5\times10^{-9}$ (random $\le2.6\times10^{-11}$); **one point $3.06\times10^{-6}$**, see N2: the oracle is the inaccurate one | $10^{-6}$ | PASS with note N2 |
| sign of $\det Q$ vs oracle, 47 + 47 (symmetric $G$) energies | equal | equal at all | exact | PASS |
| tail: $Q$ (LU solve, $N_{\rm proj}$ RHS) vs own NumPy $V^\dagger(F_3^{-1}+\mathcal K_3)V/L^6$ with **explicit `np.linalg.inv`** and $F_3$ rebuilt from fast's $F,G,\mathcal K_2^{-1}$, 12 anchors | $\le10^{-8}$ | $Q$ rel. diff max $1.7\times10^{-15}$; $F_3$ max $6.0\times10^{-15}$; $F_3^{-1}$ max $2.4\times10^{-14}$ | $10^{-8}$ | PASS |
| tail: $\ln\lvert\det Q\rvert$ vs NumPy `slogdet` of own $Q$ / vs `eigvalsh` | consistent | `slogdet`: max $2.2\times10^{-9}$ (at $\kappa(Q)=3.8\times10^{12}$), all others $\le1.8\times10^{-12}$; `eigvalsh`: max $1.6\times10^{-7}$ (same point), others $\le2\times10^{-10}$; sign equal in all three at 12/12 | $10^{-6}$ | PASS |
| tail: $Q$ Hermitian, real | Hermitian | $\lVert Q-Q^\dagger\rVert/\lVert Q\rVert\le2.7\times10^{-16}$; $\max\lvert{\rm Im}\,Q\rvert=0$; `sign_imag` 0 | $10^{-10}$ | PASS |
| tail: `slogdet_lu` on random real symmetric matrices ($n$ = 1…103, 180 cases, 94 negative-det) vs eigenvalues | exact sign, accurate ln | 0 sign errors; worst rel. diff $2.1\times10^{-13}$; complex Hermitian (60 cases): 0 sign errors, $1.0\times10^{-14}$, ${\rm Im}$ phase $\le7.8\times10^{-14}$; swap matrix $\det=-1$ gives sign $-1$; $10^{-200}I_3$ gives $-1381.55$ (exact) | n/a | PASS |
| tail: singular / non-finite | raises | zero-pivot matrix: `domain_error "Q is singular"`; NaN, Inf matrix: `"Q contains non-finite entries"`; via Python `c=(nan,0,0,0)` and `(inf,0,0,0)`: `ValueError: fast: Q contains non-finite entries` | n/a | PASS |
| tail: complex-$V$ branch (`projector="eigen"`) | exercised and correct | with `projector="eigen"` $V$ is exactly real at all 5 irreps (the branch is unreachable in production; the unit test `test_complex_projector_branch...` does not enter it). Forced through a mutant projector ($V\to VU$, $U$ a seeded random unitary; `harness_complexV.cpp`): $Q'=U^\dagger QU$ to $3\times10^{-16}$…$5.5\times10^{-16}$, $Q'$ Hermitian $\le5.7\times10^{-16}$, sign equal, $\Delta\ln\lvert\det\rvert$ $\le6.4\times10^{-10}$ at 4 of 5 points; at 100_A2 0.2915561 ($\kappa=3.8\times10^{12}$) $1.1\times10^{-5}$ (note N4) | $10^{-6}$ | PASS (production), note N4 |
| non-production parameters (eps_H 0.1–0.5, $\xi$ 3.0–3.2, other masses, $\eta_2=0.7$, $\alpha=0.4$, `max_shell_num`=12, `q0_consistent`), 4 points | K3, G, $F$, $\ln\lvert\det\rvert$ agree | $G=0.0$; $\mathcal K_3\le6.9\times10^{-17}$; $F\le8.3\times10^{-15}$; $\Delta\ln\lvert\det\rvert\le2.2\times10^{-12}$; signs equal | $10^{-12}$ | PASS |
| timing: $\mathcal K_3$ stage at 111_A2 0.35103 | $\le5$ ms | 0.92 ms (warm, median of 7; 0.75 ms in another run) | 5 ms | PASS |
| timing: $G$ stage at 111_A2 0.35103 | $\le3$ ms | 0.25 ms (0.44 ms in the interleaved A/B run, ordered) | 3 ms | PASS |
| timing: 111_A2 0.35103 warm total (wall) | $\le25$ ms | 13.8 ms (C++ total 13.4 ms; 15 reps: 13.75–15.9 ms) | 25 ms | PASS |
| speed-up at 111_A2 0.35103, 1 core, BLAS unpinned | $\ge150\times$ | **581×** with the oracle's projector build (oracle 8.13 s median of 7, fast 13.98 ms median of 15), **554×** without it (oracle projector 0.386 s) | 150× | PASS |
| CPU/wall unpinned, steady state, $N=160$ and $N=103$ | $\le1.1$ | $N=160$: 0.97–1.00 (5 runs with 0.5 s delay after import, 3 runs of `cpuwall.py`); $N=103$: 0.98–1.00 (5 + 3 runs). With `OPENBLAS_NUM_THREADS=1`: 0.94–1.00 | 1.1 | PASS |
| CPU/wall in the first ~0.3 s of a fresh process (no delay) | $\le1.1$ | $N=160$ (0.3–0.4 s window): 0.99, 1.00, 1.34, 1.05, 1.17; $N=103$ (0.15 s window): 1.80, 1.45, 1.51, 1.59, 1.98; with `OPENBLAS_NUM_THREADS=1`: 0.98–1.00 for both | 1.1 | note N1 (numpy import, not the fast path) |

## Notes

### N1. `test_fast_single_thread_cpu_time` alone: cause is OpenBLAS start-up, not the fast path (verified)

Run alone and unpinned the locked test fails (CPU/wall 4.00, 3.85, 3.68). It passes in the full suite (earlier tests burn the spin-up first), with `OPENBLAS_NUM_THREADS=1` (3/3), and when 0.3 s of idle time is inserted after the warm-up (ratio 0.978; 0.5 s: 0.97–1.00 at $N=160$ and $N=103$).
- **Bare probe, no qc3, no BLAS call:** `import numpy` followed by six 0.1 s sleeps gives CPU/wall of the first sleep $=1.95, 1.60, 2.22$ (3 runs; the import itself costs 0.21–0.30 s of CPU in 0.06–0.08 s of wall) and $\approx7\times10^{-4}$ for the next five. With `OPENBLAS_NUM_THREADS=1` the first window is $8\times10^{-4}$ and the import costs 0.060 s CPU in 0.061 s wall. So the extra CPU is the OpenBLAS worker threads (4 here) spinning for roughly 0.1–0.3 s after the library is loaded; the process then has 4 threads ("threads_in_process": 4; 1 with the env var).
- The fast path adds no BLAS call: `qc3.fast` no longer imports `qc3.qc` (AST check in the unit test, and I read `fast.py`/`fast_point.hpp`: Eigen only, no `-lblas`). Steady-state CPU/wall is $0.97$–$1.00$ at the large points (report 009 measured 3.06 at $N=160$ before 7.6).
- The test window is only $\sim0.04$ s ($20\times2$ ms at $N=35$), so it is inside the spin-down. Result: the locked test is order-dependent. It is the planner's decision whether to change it (for example by idling 0.5 s after the warm-up, or running in a subprocess with the env var); I did not touch it. The 7.6 done-check "CPU/wall $\le1.1$ without `OPENBLAS_NUM_THREADS`" holds in steady state and does not hold in the first 0.3 s after `import numpy`.

### N2. At $\kappa(Q)\gtrsim10^{13}$ the oracle is less accurate than the fast path (verified with 60 digits)

At 111_A2, $E_{\rm cm}=0.3382011860823352$ (a random draw $1.1\times10^{-7}$ above a dimension jump, $N=112$, $N_{\rm proj}=30$, $\kappa(Q)=1.5\times10^{13}$) fast and oracle differ by $\Delta\ln|\det Q|=3.06\times10^{-6}$, above the $10^{-6}$ tolerance of the locked comparisons. $\mathcal K_3$ and $G$ are identical to the oracle there ($6\times10^{-17}$, 0.0), so the kernels are not the cause. An mpmath evaluation at 60 digits from the same matrices (`illcond.py`) gives $\ln|\det Q|_{\rm HP}=-119.84557658436$:

| value | $\ln\lvert\det Q\rvert$ | minus HP |
|---|---|---|
| fast, orbit $V$ | $-119.84557658436087$ | $+1.2\times10^{-12}$ |
| fast, eigen $V$ | $-119.8455759450704$ | $+6.4\times10^{-7}$ |
| oracle | $-119.8455796476132$ | $-3.06\times10^{-6}$ |

So the fast path is the more accurate one; the $3\times10^{-6}$ discrepancy is the oracle's round-off (explicit NumPy inverse + `slogdet`). This point is not in the locked tests (they pass). It supports the earlier note for the auditor (tester 002): the $10^{-6}$ oracle tolerance is below the oracle's own accuracy for $\kappa(Q)\gtrsim10^{13}$.

### N3. $G$ fill mode: recommendation to keep the default "every ordered pair"

| quantity | ordered pairs (default) | symmetric (plan §2.3) |
|---|---|---|
| $G$ vs oracle, 47 energies | 0.0 (bit-identical) | max $1.03\times10^{-13}$ (worst: 111_A2 0.337697, where the oracle's own asymmetry is $1.46\times10^{-13}$) |
| $\Delta\ln\lvert\det Q\rvert$ vs oracle, 12 anchors | max $3.63\times10^{-9}$ | max $1.31\times10^{-8}$ (100_A2 0.2915561) |
| margin to the $10^{-6}$ tolerance | 275× | 76× |
| exact (60-digit) $\ln\lvert\det Q\rvert$ of its own matrices at 100_A2 0.2915561 | $-118.96068499405698$ | $-118.96068498195734$ (differs from the ordered one by $+1.21\times10^{-8}$) |
| double-precision error vs own exact value there | $1.9\times10^{-9}$ | $-7.8\times10^{-10}$ |
| $G$ stage, $N=160$ / 121 / 103 / 56 (interleaved, 41 reps) | 0.442 / 0.143 / 0.103 / 0.037 ms | 0.298 / 0.110 / 0.084 / 0.031 ms |

The difference between the two modes is real physics-level at this precision: for an ill-conditioned point ($\kappa(Q)=3.8\times10^{12}$) the 1e-14 asymmetry of the oracle's $G$ moves the exact $\ln|\det Q|$ by $1.2\times10^{-8}$. Both modes are inside the tolerance. **Recommendation: keep the default (ordered pairs).** It reproduces the oracle's $G$ exactly (so nothing to explain at any point), costs $+0.14$ ms at $N=160$ ($\approx1\%$ of the 13.8 ms warm total), and the plan's own §7 fall-back says to do exactly this if a deviation grows. The coordinator should record in `prompts/002` §2.3 / `docs/METHODS.md` §15 that the symmetric fill is optional (`g_symmetric=True`, reachable only through `qc3._core.fast_point`; `qc3.fast.evaluate_point` does not expose it). The $\mathcal K_3$ kernel does use the symmetric mirror (verified: $\le9\times10^{-17}$ vs oracle).

### N4. Minor

- The complex-$V$ tail branch is unreachable in production (both projector methods give real $V$ for all five irreps; no other irreps parse). It is correct when forced (table above), but for a complex $V$ at $\kappa(Q)\sim4\times10^{12}$ the complex LU leaves `sign_imag` $=1.1\times10^{-5}$ and $\Delta\ln|\det|=1.1\times10^{-5}$, and `qc3.fast` raises "determinant phase has imaginary part ... Q is not Hermitian" for `sign_imag > 1e-10` although $Q$ is Hermitian to $5\times10^{-16}$. Only relevant if a complex-$V$ sector is ever added; the unit test `test_complex_projector_branch_agrees_with_real_branch` passes without entering that branch.
- Builder's warm range "1.9–12.9 ms" is not reproduced on this machine: I measure 3.3–13.8 ms wall (12 anchors, median of 7, `taskset -c 0`); the $N=160$ value (13.8 ms) matches the claim (13.7 ms). Run-to-run noise here is large (the $F$ stage at $N=160$ ranged 8.8–12.0 ms between runs). All values are below the 25 ms bar.
- The builder's "545×" is not the number I get (581×/554× with/without the oracle projector); the order of magnitude agrees. The locked test only requires 150×.
- $K_2^{-1}$ is not changed by this step and is still computed by the oracle's function in complex arithmetic (0.3 ms at $N=160$); `fast_point.hpp` raises if its imaginary part exceeds $10^{-14}$ relative.

## Physics review (the formulas, not just the numbers)

`fast_k3.hpp` against `K3B_element`, `K3E_element`, `K3_element_raw`, `two_plus_one_block_factor` and `fast_g.hpp` against `G_ij_lm` / `G_2plus1_mat` (read line by line):
- **$\mathcal K_3$ B-part.** $(2,2)$: $2(E_{\rm cm}^2-E(\omega_p+\omega_k)+\mathbf P\cdot(\mathbf p+\mathbf k)+M_\pi^2-4M_K^2)/M_{\rm tot}^2$. $(1,1)$: $E(\omega_p+\omega_k)-\mathbf P\cdot(\mathbf p+\mathbf k)-6M_K^2+\omega_{p'_1}\,{\rm pms}_0+\omega_{k'_1}\,{\rm kms}_0$, with $-\tfrac23 Y_{1m}(\mathbf p'_1)$ in the $(p,s)$ and $(s,p)$ blocks. $(1,2)$: $E\omega_p-\mathbf p\cdot\mathbf P+E_{\rm cm}^2-2E\omega_k+2\mathbf P\cdot\mathbf k+\omega_{p'_1}{\rm pms}_0+M_\pi^2-7M_K^2$. All identical, same operation order; $(2,1)$ block taken as the transpose of $(1,2)$, which equals the oracle's $(2,1)$ formula with roles swapped (the numerical comparison confirms: $\le9\times10^{-17}$).
- **E-part.** The boosted vectors $\mathbf p^*_k,\mathbf k^*_p,\mathbf P2_{psk},\mathbf P2_{ksp}$ use the stored boost parameters of $(E-\omega,\mathbf P-\mathbf p)$ and reproduce `::boost()` (checked against `functions.h` lines 1281–1326: same $\hat{\mathbf P}$, $\beta$, $\gamma$, exact-zero shortcut, `common` expression). $\gamma_p,\gamma_k,\beta,\hat\beta$, $V$, $V'$, $t_{ij}$ including the $(I+2)\%3$ index map and the $-1$ diagonal are identical; $P2_{psk,0}=\sqrt{\sigma_p+|\mathbf P2_{psk}|^2}$ uses $\sigma$ of the row momentum, $P2_{ksp,0}$ that of the column momentum, as in the oracle. Spin-0/1 blocks: $-\tfrac13 Y_1(\mathbf P2_{ksp}+V')$ and $-\tfrac13 Y_1(\mathbf P2_{psk}+V)$, $-\tfrac23 t$ for $(p,p)$.
- **Flavour blocks and $K_{\rm iso}$.** `two_plus_one_block_factor` is called directly (1, $1/\sqrt2$, $1/\sqrt2$, 1/2). $K_{\rm iso}=\{1,0\}\to1$, $\{0,1\}\to\Delta=E_{\rm cm}^2/M_{\rm tot}^2-1$ in the $(s,s)$ entries only. The K3 has no cutoff $H$ (neither does the oracle) and no $\eta$, so those cannot differ; the Q0 normalisation $|p|^\ell Y_{\ell m}$ is `spherical_harmonics`: $Y_{1,-1}=\sqrt3\,y$, $Y_{10}=\sqrt3\,z$, $Y_{11}=\sqrt3\,x$, matched by `Yq0`.
- **$G$.** $\mathbf k^*=$ boost$(\omega_k,\mathbf k,E-\omega_p,\mathbf P-\mathbf p)$, $\mathbf p^*=$ boost$(\omega_p,\mathbf p,E-\omega_k,\mathbf P-\mathbf k)$; $Y_{\ell_f m_f}(\mathbf k^*)\,Y_{\ell_i m_i}(\mathbf p^*)$; $1/(2\omega L^3)$ per side; $H(\sigma_p;m_j,m_k)H(\sigma_k;m_i,m_k)$ with the oracle's $(m_j,m_k)$: $G_{11}$: $(M_K,M_\pi)$ both; $G_{12}$: $(M_\pi,M_K)$ and $(M_K,M_K)$; $G_{21}$: $(M_K,M_K)$ and $(M_\pi,M_K)$; the table `cut[v]` in `g_momentum` reproduces exactly these. Denominator $(E-\omega_p-\omega_k)^2-|\mathbf P-\mathbf p-\mathbf k|^2-m_k^2$ with $m_k=M_\pi$ for $G_{11}$ and $M_K$ otherwise, including the square-root-then-square of the momentum, as in the oracle. Flavour factors: $G_{12}$: $\sqrt2(-1)^{\ell_f}$, $G_{21}$: $\sqrt2(-1)^{\ell_i}$, with the oracle's multiplication order; $(2,2)=0$. This explains why the default is bit-identical (0.0).
- **Real-arithmetic re-expression.** Valid in the production range: all complex quantities have zero imaginary part for real $E$, real momenta and $\sigma>0$. For a kaon spectator $q^{*2}+M_{1,2}^2=(\sigma\pm(M_1^2-M_2^2))^2/(4\sigma)\ge0$, so ${\rm pms}_0$ is real, and $\gamma^{-2}=\sigma/(E-\omega)^2>0$. `check_real` throws on any non-zero imaginary part or NaN at every conversion, and the result matrices are checked with `allFinite`. **The `pms_0` omission for a pion spectator is correct:** in `K3B_element` ${\rm pms}_0$ (kaon-row) appears only for $i=1$, ${\rm kms}_0$ only for $j=1$ ($(2,2)$ uses neither, $(1,2)$ uses ${\rm pms}_0$ of the kaon row, $(2,1)$ uses ${\rm kms}_0$ of the kaon column); in `K3E_element` $(1,2)$ uses ${\rm pms}_0$ ($i=1$), $(2,1)$ uses ${\rm kms}_0$ ($j=1$), $(1,1)$ both kaon, $(2,2)$ neither. No branch reads a pion-spectator ${\rm pms}_0$/${\rm kms}_0$. The numbers confirm it at 47 energies.
- **Tail (invariants 1, 2).** $Q=L^{-6}V^T(X/s+\mathcal K_3V)$ with $X=(F_3/s)^{-1}V$ is exactly $L^{-6}V^T(F_3^{-1}+\mathcal K_3)V$ (project the inverse; matrix scaled by $s=\max|F_3|$, not the determinant; no $N_{\rm proj}$ factor); numerically confirmed against the explicit inverse ($1.7\times10^{-15}$). $(\mathrm{sign},\ln|\det|)$ come from the partial-pivot LU, sign $=\det\mathsf P\prod{\rm sgn}\,U_{kk}$, singular or non-finite raises (no fallback). `F3inv` is formed only when requested (tested: `want=("Q",)` leaves all other fields `None`; subset results equal the full ones).

## Per-stage timings (`taskset -c 0`, BLAS not pinned, medians of 7 warm calls; `results/tester_010_timing.json`)

| irrep | $E_{\rm cm}$ | $N$ | $N_{\rm proj}$ | warm wall ms | cold wall ms | cold $-$ warm ms | basis | F | $\mathcal K_2^{-1}$ | G | F3 | proj. | F3inv | $\mathcal K_3$ | projection | LU |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 100_A2 | 0.291556 | 56 | 17 | 6.6 | 22.0 | 15.4 | 0.51 | 5.28 | 0.06 | 0.07 | 0.18 | 0.05 | 0.05 | 0.12 | 0.06 | 0.01 |
| 100_A2 | 0.291798 | 56 | 17 | 6.6 | 22.0 | 15.4 | 0.49 | 5.27 | 0.06 | 0.07 | 0.17 | 0.05 | 0.05 | 0.12 | 0.07 | 0.01 |
| 000_A1m | 0.26932 | 35 | 5 | 3.3 | 18.9 | 15.6 | 0.47 | 2.38 | 0.04 | 0.03 | 0.07 | 0.03 | 0.02 | 0.05 | 0.02 | 0.00 |
| 000_A1m | 0.34538 | 103 | 9 | 6.4 | 22.9 | 16.5 | 0.36 | 4.81 | 0.09 | 0.11 | 0.43 | 0.04 | 0.12 | 0.22 | 0.05 | 0.00 |
| 100_A2 | 0.3599 | 108 | 30 | 6.6 | 18.6 | 12.0 | 0.37 | 4.86 | 0.09 | 0.11 | 0.45 | 0.05 | 0.13 | 0.22 | 0.16 | 0.01 |
| 100_A2 | 0.36 | 108 | 30 | 6.8 | 14.1 | 7.4 | 0.36 | 4.98 | 0.08 | 0.11 | 0.43 | 0.04 | 0.12 | 0.22 | 0.17 | 0.01 |
| 110_A2 | 0.30733 | 68 | 27 | 4.6 | 12.7 | 8.1 | 0.33 | 3.66 | 0.05 | 0.05 | 0.14 | 0.03 | 0.04 | 0.09 | 0.07 | 0.01 |
| 110_A2 | 0.35 | 108 | 41 | 7.8 | 15.0 | 7.2 | 0.34 | 5.89 | 0.10 | 0.12 | 0.45 | 0.05 | 0.14 | 0.24 | 0.25 | 0.02 |
| 111_A2 | 0.32036 | 64 | 20 | 9.9 | 18.2 | 8.3 | 0.43 | 8.51 | 0.07 | 0.08 | 0.23 | 0.05 | 0.06 | 0.16 | 0.09 | 0.01 |
| 111_A2 | 0.35103 | 160 | 40 | 13.8 | 38.5 | 24.7 | 0.47 | 8.76 | 0.38 | 0.25 | 1.55 | 0.14 | 0.51 | 0.92 | 0.48 | 0.02 |
| 200_A2 | 0.28664 | 47 | 16 | 3.5 | 11.5 | 8.1 | 0.36 | 2.67 | 0.03 | 0.03 | 0.07 | 0.03 | 0.02 | 0.05 | 0.03 | 0.00 |
| 200_A2 | 0.34738 | 121 | 32 | 7.9 | 16.2 | 8.3 | 0.34 | 5.57 | 0.11 | 0.15 | 0.61 | 0.05 | 0.17 | 0.31 | 0.21 | 0.01 |

Cold = `clear_caches()` first (projector and threshold tables rebuilt); the cold excess is dominated by the one-time threshold-table build (7.7 ms for both channels) plus projector (0.2–0.8 ms), all $\le24.7$ ms $\ll0.1$ s. At 111_A2 0.35103, 30 paired cold/warm repetitions: total $21.9$ vs $13.0$ ms, difference median $9.25$ ms (IQR 8.5–11.4, max 14.2). Cold projector build at all 58 production bases: max 0.95 ms. "F3inv" is the LU of $F_3/s$ (the explicit inverse is not formed when only $Q$ is wanted). At 111_A2 0.35103 the $F$ stage is now 64 % of the total; $\mathcal K_3$ + $G$ together are 1.2 ms (was 54 ms in report 009).

## Failing output

Only the standalone run of the locked CPU-time test:
```
E   AssertionError: (0.17123654, 0.042766096999912406, 4.0040254316485955)
E   assert (0.0427... > 0 and (0.1712... / 0.0427...) <= 1.3)       # cpu, wall, cpu/wall
```
Best guess (verified by the bare-`numpy` probe above, so I call it the cause rather than a guess): OpenBLAS worker threads spinning for ~0.1–0.3 s after `import numpy`; the test's 0.04 s window falls inside it.

## Items not done / limits

- Pinned and unpinned full suite were run once each; the speed-up test passed both times (its measured ratio is in the timing section, not printed by pytest).
- $d_\Gamma>1$ irreps and $\ell\ge2$ waves are rejected by the kernels (`only l <= 1 is implemented`), as in the plan; untested by design.
- `mpmath` was not installed; I ran `pip install mpmath` in the container for the 60-digit checks (`illcond.py`, `gmode_hp.py`); it is not a project dependency.
- No benchmark was run (mode `step`); no `reports/benchmarks` comparison.

# Tester report 012: M1 sub-step 7.7, class moment tables + shared (canonical $\mathbf n_P$, $|\mathbf p|^2$) sums for the $F$ lattice sum

- **Date:** 2026-10-22
- **Git commit:** `5792545` (HEAD while testing; the tested code is identical to the merge commit `d444181`, the only difference is `graphify-out/`, 5 files, `git diff --stat d444181 5792545`). `_core` reports version `0.5.1-step7.7`.
- **Prompt in force:** `prompts/002-m1-step7-speed-replan-projector-basis-k3-g.md` (row 7.7, stop rule, §2.7, §4); `prompts/001` rung 5; `docs/METHODS.md` §3, §15. Acceptance tag `acceptance-v4`.
- **Mode:** step.
- **Related:** `reports/tests/010-m1-step7-4-5-6-kernels-tail.md`, `reports/tests/011-m1-step8-bench-harness.md`, `reports/tests/006-m1-step6-fast-path-v1.md` (gamma-stress energies).
- **Raw numbers:** `results/tester_012_*.json|txt`, `results/bench_tester_012_fastonly_step7_7.json`; scripts in `results/tester_012_scripts/`.

**VERDICT: PASS** (with notes N1 to N4 below; none blocks the sub-step).

## Commands run

```
rm src/python/qc3/_core*.so && ./build.sh core                          # clean rebuild, 98 s
python3 -m pytest -q tests/acceptance/m1                                # 47 passed (slow included), 4 min 57 s
python3 -m pytest -q tests --ignore=tests/acceptance                    # 74 passed, 33.9 s
.claude/bin/acceptance-check                                            # PASS: matches acceptance-v4
git diff 73c4e2d -- src/qc3_core/oracle_point.hpp src/core src/fitter src/cachegen_cpu   # 0 bytes
scripts/acc.py        (F, ln|det Q| vs oracle and vs use_classes=False; 12 anchors + 15 random (seed 20261022) + 5 gamma-stress)
scripts/sweep.py 300  (classes vs direct only: 1606 energies, 5 irreps, seed 20261022, 300 random + 11..16 dimension jumps x (+-3e-9) per irrep)
scripts/outl.py       (oracle on the 4 largest sweep outliers)
scripts/momcheck.cpp  (C++ harness, long-double brute force vs class moments + rotation, 10 canonical c x 8 kinematic cases, every element of each orbit)
scripts/tabsize.cpp + scripts/tables.py + scripts/stracerun.py  (memory, file I/O)
scripts/facc.py + scratch build with -DQC3_FACC="long double"    (accumulator type)
taskset -c 0 env OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 python3 scripts/timing.py     (21 warm, 9 cold repeats per anchor, classes and direct in one session)
PYTHONPATH=src/python python3 -m qc3.bench --out results/bench_tester_012_fastonly_step7_7.json --no-oracle --repeats 7 --proxy-repeats 5
```

## 1. Gate checks

| check | expected | got | tolerance | result |
|---|---|---|---|---|
| clean rebuild of `_core` | builds | built, version `0.5.1-step7.7` | n/a | PASS |
| `tests/acceptance/m1` (incl. slow) | 47/47 | 47 passed | exact | PASS |
| unit tests (`tests --ignore=tests/acceptance`) | all pass | 74 passed (incl. the 2 new 7.7 tests, 6 parametrised class-vs-direct cases) | exact | PASS |
| `acceptance-check` | PASS | PASS (matches acceptance-v4) | n/a | PASS |
| oracle untouched | empty diff | 0 bytes | exact | PASS |
| `graphify affected F_channel` | covered | only `evaluate()` (fast path), exercised by all fast-vs-oracle tests | n/a | PASS |

## 2. Math review

**(a) Sharing.** The summand is $Y_a(\mathbf r)Y_b(\mathbf r)\,e^{\alpha(x^2-r^2)}/(x^2-r^2)$ with $\mathbf r=\mathbf n+\mathbf n_P\,\mathrm{prod1}(\mathbf n\!\cdot\!\mathbf n_P)$, $\mathrm{prod1}(d)=\tfrac{d}{n_P^2}(\tfrac1\gamma-1)-\tfrac{\xi}{\gamma}$. The inputs $x^2,\gamma,\xi$ come from $\sigma_p=(E-\omega_p)^2-|\mathbf P-\mathbf p|^2$, i.e. from $(E,|\mathbf p|^2,n_P^2)$ and the channel masses; $\mathbf p$ enters nowhere else in the sum (the cutoff, the prefactor $A$ and the integral term $C_{ab}$ are outside and still per block). For a signed permutation $R$ with $R\mathbf c=\mathbf n_P$ the cube $[-n_{\max},n_{\max}]^3$ maps to itself and $\mathbf n=R\mathbf m$ gives $\mathbf n\!\cdot\!\mathbf n_P=\mathbf m\!\cdot\!\mathbf c$, $\mathbf r=R\mathbf r_{\rm can}$, hence $T_0$ is invariant, $T_1=R\,T_1^{\rm can}$, $T_2=R\,T_2^{\rm can}R^T$. The code does `T1[idx[j]] = sgn[j] T1c[j]`, `T2[idx[j]][idx[k]] = sgn[j] sgn[k] T2c[j][k]`, which is exactly $R$ with $R e_j=\mathrm{sgn}_j e_{{\rm idx}_j}$ and $R\mathbf c=\mathbf n_P$ (checked: ties in $|n_P|$ components only permute equal entries). Real harmonics: $Y_{1,-1}=\sqrt3\,r_y$, $Y_{10}=\sqrt3\,r_z$, $Y_{11}=\sqrt3\,r_x$; the code's $v_a$ uses index 1, 2, 0 for $m=-1,0,+1$, the same as `Yreal` and the oracle's `spherical_harmonics`. Numerically (`momcheck`): for all 10 canonical $c$ (000, 100, 110, 111, 200, 210, 211, 321, 222, 300), every element of each orbit (1144 pairs with 8 kinematic cases incl. $x^2<0$, $x^2>0$, $\gamma$ 1 to 9.98) the rotated class moments agree with a long-double brute force over the actual $\mathbf n_P$; all term counts are equal; worst error $3.0\times10^{-13}$ of $\sum|{\rm terms}|$ (see N1: the double-precision direct sum has $4.1\times10^{-13}$ in the same test, so this is conditioning, not a rotation error). Whole-$F$ agreement with the oracle at 32 points (all spectator orbits of the bases occur) confirms the $\ell=1$ rotation end to end.
*Share key.* The cache `fcache` is a local variable of `F_channel`, so it is per flavour channel (masses cannot mix). The key is the exact integers $(c_0,c_1,c_2,|\mathbf p|^2)$; the $10^{-10}$ relative tolerance on $x^2,\gamma,\xi$ is only a consistency guard: for equal $(c,|\mathbf p|^2)$ these are mathematically identical and differ by $\sim10^{-16}$, while different $c$ or $|\mathbf p|^2$ can never merge because the integers differ. Measured: at all 12 anchors `f_sums` equals the number of distinct (flavour, canonical $c$, $|\mathbf p|^2$) in the basis, i.e. no group was split and none merged (`results/tester_012_share_counts.json`; e.g. 111_A2 0.35103: 64 blocks, 20 sums).

**(b) Class moments.** $r^2=n^2+2\,\mathrm{prod1}\,(\mathbf n\!\cdot\!\mathbf n_P)+\mathrm{prod1}^2n_P^2$ depends on $\mathbf n$ only through $(n^2,\mathbf n\!\cdot\!\mathbf n_P)$, so $w$ is constant on a class. With $\ell_a,\ell_b\le1$, $Y_aY_b$ is at most quadratic in $\mathbf r$, so count, $\sum\mathbf n$, $\sum\mathbf n\mathbf n^T$ per class are sufficient and exact; the code stores integer sums and forms the centred $C=(\mathrm{cnt}\sum nn^T-\sum n\sum n)/\mathrm{cnt}$ from exact `long long` (maximum about $2\times10^{12}$, no overflow); $T_2=\sum w(C+\mathrm{cnt}\,\bar{\mathbf r}\bar{\mathbf r}^T)$ is algebraically exact with $\bar{\mathbf r}=\bar{\mathbf n}+\mathbf n_P\mathrm{prod1}$. $\mathrm{perp2}=(n^2n_P^2-d^2)/n_P^2$ is an exact rational (integer numerator); $r_\parallel=(d+n_P^2\mathrm{prod1})/|\mathbf n_P|$; $r^2=\mathrm{perp2}+r_\parallel^2$. $e^{\alpha(x^2-r^2)}=e^{\alpha(x^2-r_\parallel^2)}e^{-\alpha\,\mathrm{perp2}}$: one $\exp$ per group. Underflow guard: used only if $\alpha\,\mathrm{avail}<575$ (so the table factor $\ge e^{-575}\approx10^{-250}$), `eg` finite and $<10^{300}$, else one $\exp$ per class. In production $\alpha\,\mathrm{avail}\lesssim40$, so the fallback is not reached; I forced it with `trunc_eps` $=10^{-300}$ ($\alpha\,\mathrm{avail}\approx690$): finite, classes vs direct $\le1.9\times10^{-15}$, identical term counts (3 points, `trunc_eps` 1e-40, 1e-100, 1e-300). Sanity: every table has 68921 $=41^3$ points (`tabsize`).

**(c) Truncation.** The class path selects $\mathrm{perp2}\le R^2-r_\parallel^2$ by binary search, the same point set as the direct test $r^2\le R^2$ (up to rounding of a term at the boundary, which is $\le\epsilon$ by the plan). Result: `lattice_terms` (classes) $=$ `lattice_terms` (direct) at **all 12 anchors, all 15 random, all 5 gamma-stress points and all 1606 sweep energies** (0 mismatches; sweep includes $\pm3\times10^{-9}$ around every dimension jump), and in the 9 high-precision cases above.

## 3. Independent numerics (coded_v2, $\vec c_{\rm base}$, seed 20261022)

| check | expected | got (max) | tolerance | result |
|---|---|---|---|---|
| $\lVert F_{\rm cls}-F_{\rm oracle}\rVert/\lVert F_{\rm oracle}\rVert$, 12 anchors | small | $5.1\times10^{-15}$ | $10^{-12}$ (acceptance) | PASS |
| same, 15 random energies | small | $5.4\times10^{-15}$ | $10^{-12}$ | PASS |
| same, 5 gamma-stress points (report 006) | small | $1.9\times10^{-15}$ | $10^{-12}$ | PASS |
| per spectator block, anchors / random / gamma-stress | small | $4.4\times10^{-14}$ / $8.2\times10^{-15}$ / $1.5\times10^{-14}$ | $10^{-12}$ | PASS (N2) |
| $F_{\rm cls}$ vs $F_{\rm direct}$ (`use_classes=False`), anchors / random / gamma | small | $3.9\times10^{-15}$ / $3.8\times10^{-15}$ / $1.8\times10^{-15}$ | $10^{-13}$ (unit test) | PASS |
| per block cls vs direct, anchors / random / gamma | small | $3.8\times10^{-14}$ / $8.2\times10^{-15}$ / $2.0\times10^{-14}$ | n/a | info |
| direct vs oracle (reference for the above) | | $5.2\times10^{-15}$ / $3.4\times10^{-15}$ / $8.5\times10^{-16}$ (whole $F$) | | info |
| `lattice_terms` equal (cls vs direct) | equal | 32/32 points and 1606/1606 sweep energies | exact | PASS |
| $F$ exactly symmetric, off-block entries exactly 0 | yes | yes (32/32) | exact | PASS |
| sign of $\det Q$ vs oracle | equal | equal at 32/32 | exact | PASS |
| $\lvert\Delta\ln\lvert\det Q\rvert\rvert$ cls vs oracle: anchors / random / gamma | $\le10^{-6}$ | $7.1\times10^{-10}$ / $4.5\times10^{-10}$ / $5.6\times10^{-9}$ (the last at 200_A2, $\kappa(Q)=5.4\times10^8$, identical for the direct path) | $10^{-6}$ | PASS |
| $\lvert\Delta\ln\lvert\det Q\rvert\rvert$ cls vs direct | | $2.9\times10^{-9}$ (anchors; 100_A2 0.2915561, $\kappa=3.8\times10^{12}$), $\le2.8\times10^{-13}$ otherwise | info | PASS |
| sweep, classes vs direct, 1606 energies: $\max\lVert\Delta F\rVert/\lVert F\rVert$ per irrep 000 / 100 / 110 / 111 / 200 | | $1.7\times10^{-15}$ / $4.2\times10^{-15}$ / $1.7\times10^{-14}$ / $1.0\times10^{-13}$ / $1.1\times10^{-13}$; 99th percentile $1.3\times10^{-14}$ | $10^{-12}$ | PASS (N1) |
| long double vs double accumulators (classes path), 21 points | | $\lVert\Delta F\rVert/\lVert F\rVert\le1.9\times10^{-15}$; $\lvert\Delta\ln\lvert\det\rvert\rvert\le8.7\times10^{-11}$ | info | see N3 |

## 4. Memory and file I/O

| check | expected | got | result |
|---|---|---|---|
| strace (`-f -e trace=file,read,write,mmap,munmap,brk`) around a cold call after `clear_caches()` and 5 warm calls | no file I/O | 0 `open/read/write/stat` calls, 0 `mmap` (only `brk`/3 `munmap`); tables are built in process memory | PASS |
| tables for the full scan, 5 irreps x 58 bases (5+12+16+13+12) | few | 8 canonical $c$: 000, 100, 110, 111, 200, 210, 211, 300 (basis enumeration = `cache_info()["f_table_entries"]`, 8) | PASS |
| size per table (`sizeof(FClass)=96`, 24 per group) | | 0.07, 0.78, 1.28, 0.95, 0.78, 2.13, 2.21, 0.78 MB (classes 727 to 23004) | |
| total of F tables | | **8.98 MB** (largest 2.21 MB); process RSS growth over the whole scan 21 MB | PASS |

## 5. Timings (`taskset -c 0`, OMP and OPENBLAS = 1, median of 21 warm / 9 cold, `want=("Q",)`, coded_v2)

Classes vs direct in the same session (`results/tester_012_timing.json`):

| anchor | N | warm wall med / min (ms) | C++ total (ms) | F stage cls vs direct (ms) | total speed-up med | cold $-$ warm (ms) |
|---|---|---|---|---|---|---|
| 100_A2 0.2915561 | 56 | 0.93 / 0.84 | 0.78 | 0.158 vs 2.83 | 3.8 | 15.2 |
| 100_A2 0.2917984 | 56 | 0.78 / 0.69 | 0.66 | 0.131 vs 2.91 | 4.7 | 14.6 |
| 000_A1m 0.2693200 | 35 | 0.49 / 0.39 | 0.39 | 0.031 vs 1.42 | 3.9 | 10.2 |
| 000_A1m 0.3453800 | 103 | 1.57 / 1.46 | 1.41 | 0.133 vs 6.78 | 5.8 | 16.0 |
| 100_A2 0.3599000 | 108 | 1.90 / 1.82 | 1.74 | 0.264 vs 4.74 | 3.4 | 21.8 |
| 100_A2 0.3600000 | 108 | 1.84 / 1.77 | 1.70 | 0.266 vs 6.56 | 4.9 | 20.6 |
| 110_A2 0.3073300 | 68 | 1.08 / 1.00 | 0.95 | 0.232 vs 3.85 | 4.4 | 20.2 |
| 110_A2 0.3500000 | 108 | 2.07 / 1.94 | 1.91 | 0.370 vs 5.33 | 3.4 | 28.5 |
| 111_A2 0.3203600 | 64 | 1.03 / 0.99 | 0.91 | 0.278 vs 4.55 | 5.3 | 25.2 |
| 111_A2 0.3510300 | 160 | 3.89 / 3.71 | 3.74 | 0.440 vs 8.49 | 3.1 | 27.9 |
| 200_A2 0.2866400 | 47 | 0.82 / 0.71 | 0.70 | 0.177 vs 3.62 | 5.6 | 20.1 |
| 200_A2 0.3473800 | 121 | 2.33 / 2.20 | 2.20 | 0.370 vs 8.06 | 4.9 | 28.2 |

- Warm wall: $0.49$ to $3.89$ ms (claim 0.45 to 4.1: reproduced). Cold minus warm: $10.2$ to $28.5$ ms (claim 10 to 29: reproduced). The cold excess is probably the F-table build plus the first threshold/projector tables (a guess; not separated).
- **Total warm speed-up, same-session A/B: $3.1\times$ to $5.8\times$** (min-based $3.2$ to $5.6\times$); against the recorded report-010 timings (`results/tester_010_timing.json`) $3.4\times$ to $9.6\times$. The quoted upper end $8.2\times$ is only reached against old reference times for some points (N4).
- **F stage speed-up: $14\times$ to $51\times$** (claim $27$ to $68\times$, not reproduced at either end, N4). $F$ is now 8 to 30 % of the C++ total over the 12 anchors (12 % at $N=160$, was 70 %).
- Stop rule (prompts/002 §2.7): $F$ was 70 % of the warm total at $N=160$ and the gain on the total there is $3.1\times\ge1.5$: PASS.
- **$\chi^2$ proxy** (`python -m qc3.bench --no-oracle`, 12 labelled roots, coded_v2, 1 core, no caches): warm median **0.140 s** (min 0.112, max 0.171, n=5), cold 0.173 s; report 011 had 0.524 s (0.487 to 0.598) warm, 0.633 s cold: **3.7x faster**. Still 102 evaluations, 12/12 `kind="zero"`, $\chi^2=53.890297$ (report 011: 53.890297; difference $9\times10^{-12}$), roots identical to report 011 to $1.1\times10^{-16}$. Per evaluation 1.37 ms. Against the old direct FCN (1.4 s): 10x.

## Notes (none blocks)

- **N1 (precision near an $F$ pole).** Class vs direct differs by up to $1.06\times10^{-13}$ of $\lVert F\rVert$ in 4 of 1606 energies (200_A2 0.346802592, 111_A2 0.337572574, 0.330752399, 0.337780731). Cause (measured with a scratch build of `fast_point.hpp` that adds a long-double brute-force reference, `results/tester_012_scripts/fast_point_debug_scratch.hpp`): at 200_A2 0.346802592 the block (mass $m_K$, $|\mathbf p|^2=2$, $c=(1,0,1)$, $x^2=0.2619$) has $\min|x^2-r^2|=2.3\times10^{-4}$, i.e. the energy lies within $2\times10^{-4}$ (in $r^2$ units) of a free-level pole of $F$. The error of a double-precision $r^2$ is amplified by $x^2/\min|\ldots|\approx10^3$. Against the exact sum (same double inputs) the classes give $6.2\times10^{-13}$ and the direct sum $1.4\times10^{-13}$ of $\sum|{\rm terms}|$; the direct path agrees with the oracle to $10^{-17}$ only because it reproduces the oracle's rounding. These energies have $\kappa(Q)\sim3\times10^{12}$ (200_A2: fast vs oracle $5.8\times10^{-6}$ in $\ln|\det Q|$ for both paths, the known N2-type ill-conditioning). The unit test's $10^{-13}$ class-vs-direct tolerance is only exercised at 6 benign points; a point like these four would exceed it, while the acceptance limit $10^{-12}$ has margin $\ge9$.
- **N2.** Per-block class-vs-oracle error is $4.4\times10^{-14}$ at the anchor 110_A2 0.35 (direct: $1.1\times10^{-14}$) and $1.53\times10^{-14}$ at gamma-stress (the builder quoted $\le1.5\times10^{-14}$); both far below $10^{-12}$.
- **N3 (accumulators).** `QC3_FACC=long double` (scratch build) changes $F$ by $\le1.9\times10^{-15}$ and costs 0 to $1.8\times$ on the $F$ stage (0.03 to 0.3 ms absolute); the near-pole discrepancies of N1 are the same with long double. So `double` accumulators are adequate; the limiting error is the rounding of $r^2$.
- **N4 (builder's quoted speed-ups).** F stage "27 to 68x" and total "3.0 to 8.2x" are not reproduced by the same-session A/B (14 to 51x and 3.1 to 5.8x). The sub-step passes its own stop rule with large margin either way; please quote the A/B numbers.
- `fast.cache_info()` has the new key `f_table_entries`; `clear_caches()` empties the F tables (unit test passes).
- HEAD moved from `d444181` to `5792545` (graph refresh only) during testing; the `_core` I tested was built from `5792545` and the source files are the same.

## Return summary

PASS. Passed: gate checks 6/6; acceptance m1 47/47; unit 74/74; independent numerical checks all within tolerance. Failed: 0. Pending: 0. New failures vs report 011: none.

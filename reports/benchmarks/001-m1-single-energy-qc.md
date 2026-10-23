# Benchmark report 001: M1 single-energy QC, fast path vs oracle (milestone gate)

- **Date:** 2026-10-23
- **Git commit:** `c0334de` (`qc3._core.git_commit()` == `c0334de`; the tree was clean apart from untracked report/result files; `git_tree_dirty: False` in the JSON). `_core` version string `0.5.1-step7.7`.
- **Prompt in force:** `prompts/002` (7.0 to 7.8), protocol `prompts/001` §7; acceptance tag `acceptance-v5`. Auditor verdict PASS WITH NOTES (`reports/audits/001`), tester gate report `reports/tests/013-m1-gate-full.md` PASS.
- **Mode:** benchmark. This is the first file in `reports/benchmarks/`, so there is no previous benchmark report. The comparison is with the history in `reports/tests/006`, `010`, `012`, `011`, and the legacy numbers (1.4 s direct FCN, 24 ms hot-window FCN).
- **Raw numbers:** `results/bench_m1_gate.json` (the official run), `results/bench_m1_gate_fastonly_rerun1.json`, `..._rerun2.json` (21 repeats, fast only), `results/bench_m1_gate_threads4_fastonly.json`, `results/bench_m1_gate_minbased_speedups.json`.

**VERDICT: PASS.** Speed-up over the oracle on one core: median 2551x (min 1841x, max 8372x) with the oracle's projector build counted; median 2408x (1761x to 4176x) without it. No point is slower than the previous measurement (report 012) by more than 10 % on minima (worst +5.3 %). The physics output is unchanged: $\ln|\det Q|$ bit-identical to the report-012 run at all 27 points, $\chi^2$ proxy value identical to all printed digits.

## Machine and build
| item | value |
|---|---|
| CPU | Intel Xeon @ 2.10 GHz (`lscpu`), 4 virtual cores online, 15 GiB RAM, virtual machine; no GPU |
| OS / Python | Linux 6.18.44-fc-v153, glibc 2.39, Python 3.13.16, NumPy 2.5.3, SciPy 1.18.1 |
| compiler | `g++ (Ubuntu 13.3.0-6ubuntu2~24.04.1) 13.3.0` |
| flags (`./build.sh core`) | `-std=c++17 -O3 -fopenmp -shared -fPIC -DQC3_GIT_COMMIT="c0334de" -DV31L_DISABLE_MINUIT=1 -DV31R_DISABLE_MINUIT=1 -DV31S_DISABLE_MINUIT=1 -DV31T_DISABLE_MINUIT=1 -DV32F_HAS_MINUIT2=0 -DV32F_DISABLE_MINUIT` + Eigen 3 + pybind11 headers + `src/core/Faddeeva.cc`. No `-march=native`, no `-ffast-math`. |
| physics constants | $\xi=3.444$, $a_tm_K=0.09698$, $a_tm_\pi=0.06906$, $\eta=(1,\tfrac12)$, $\alpha=0.5$, $\epsilon_H=0$, `max_shell_num`=20, waves $\{0,1\}/\{0\}$, parity $-1$, `eig_tol`=0.05, scatter $(4.04,-43.2,4.12)$, $L/a_s=20$, $\vec c=\vec c_{\rm base}$; `k2_pwave_form` = `q0_consistent` (grid and anchor points), `coded_v2` for the $\chi^2$ proxy only (the labels are `coded_v2` zeros). Recorded in the JSON `constants`. |
| seed | none: fixed energies and windows, no random numbers (JSON `seeds`) |

## Commands
```
# official run (core 0, one thread everywhere; peak RSS from a 6-line getrusage wrapper because /usr/bin/time is absent)
PYTHONPATH=src/python taskset -c 0 env OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 \
    python3 <scratch>/runbench.py --out results/bench_m1_gate.json --core 0      # runbench.py = runpy.run_module("qc3.bench") + resource.getrusage; 117 s total
# noise check, fast path only, 21 repeats, same pinning (twice)
... python3 -m qc3.bench --no-oracle --no-proxy --repeats 21 --core 0 --out results/bench_m1_gate_fastonly_rerun{1,2}.json
# not the metric (protocol §7): 4 threads requested, 4 cores allowed
python3 -m qc3.bench --threads 4 --no-oracle --out results/bench_m1_gate_threads4_fastonly.json
```
Protocol: 15 grid points (5 irreps x $E^*\in\{0.2700,0.3100,0.3599\}$) + 12 anchors = 27 points. Fast: 1 warm-up + 5 timed repeats, plus one cold call after `clear_caches()`. Oracle: one call per point (it takes 1 to 16 s; its noise is therefore not measured, see caveats). Speed-up = oracle total stage / fast warm total stage (medians). "Total" columns are the C++ stage `total`; "wall" includes the Python/pybind call overhead (about 0.1 ms).

## Per-point table (official run; core 0, 1 thread)
| set | irrep | $E^*$ | $N$ | $N_{\rm proj}$ | oracle total [ms] (of which projector) | fast warm total, median [ms] | fast warm total, min [ms] | fast warm wall, median [ms] | fast warm wall, min of 3 runs [ms] | fast cold wall [ms] | speed-up total/total, with oracle projector | same, without oracle projector | oracle wall / fast wall min-of-3 | $|\Delta\ln|\det Q||$ fast vs oracle |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| grid | 000_A1m | 0.270000 | 35 | 5 | 1908 (229) | 0.402 | 0.352 | 0.506 | 0.391 | 11.1 | 4746 | 4176 | 4884 | 0.0e+00 |
| grid | 000_A1m | 0.310000 | 95 | 8 | 8222 (3960) | 1.105 | 1.063 | 1.223 | 1.184 | 14.3 | 7439 | 3856 | 6943 | 7.1e-15 |
| grid | 000_A1m | 0.359900 | 135 | 11 | 16203 (10134) | 1.935 | 1.813 | 2.143 | 1.996 | 16.2 | 8372 | 3136 | 8118 | 1.4e-14 |
| grid | 100_A2 | 0.270000 | 50 | 14 | 2236 (17) | 0.549 | 0.462 | 0.650 | 0.535 | 11.3 | 4072 | 4041 | 4182 | 0.0e+00 |
| grid | 100_A2 | 0.310000 | 64 | 20 | 2467 (40) | 0.780 | 0.733 | 0.911 | 0.842 | 14.6 | 3162 | 3110 | 2933 | 2.8e-14 |
| grid | 100_A2 | 0.359900 | 108 | 30 | 4722 (166) | 1.589 | 1.509 | 1.726 | 1.643 | 20.4 | 2971 | 2866 | 2874 | 2.3e-13 |
| grid | 110_A2 | 0.270000 | 28 | 13 | 986 (3) | 0.387 | 0.333 | 0.472 | 0.421 | 12.3 | 2549 | 2541 | 2342 | 2.8e-14 |
| grid | 110_A2 | 0.310000 | 68 | 27 | 2814 (20) | 0.877 | 0.828 | 0.999 | 0.944 | 17.3 | 3208 | 3185 | 2983 | 4.5e-13 |
| grid | 110_A2 | 0.359900 | 116 | 43 | 4790 (113) | 2.135 | 2.079 | 2.298 | 2.213 | 26.4 | 2243 | 2191 | 2166 | 3.7e-13 |
| grid | 111_A2 | 0.270000 | 40 | 14 | 1797 (8) | 0.436 | 0.399 | 0.532 | 0.487 | 12.5 | 4120 | 4100 | 3694 | 2.8e-14 |
| grid | 111_A2 | 0.310000 | 52 | 17 | 1900 (14) | 0.734 | 0.656 | 0.870 | 0.754 | 18.2 | 2587 | 2568 | 2521 | 1.4e-14 |
| grid | 111_A2 | 0.359900 | 160 | 40 | 7098 (411) | 3.357 | 3.253 | 3.552 | 3.353 | 26.3 | 2114 | 1992 | 2118 | 3.1e-13 |
| grid | 200_A2 | 0.270000 | 43 | 15 | 1577 (14) | 0.567 | 0.533 | 0.693 | 0.582 | 16.8 | 2781 | 2756 | 2710 | 1.4e-14 |
| grid | 200_A2 | 0.310000 | 95 | 25 | 3862 (110) | 1.437 | 1.255 | 1.573 | 1.373 | 18.1 | 2688 | 2612 | 2813 | 4.3e-14 |
| grid | 200_A2 | 0.359900 | 137 | 35 | 6334 (424) | 3.222 | 2.720 | 3.428 | 2.583 | 26.8 | 1966 | 1834 | 2452 | 4.3e-12 |
| anchor | 100_A2 | 0.291556 | 56 | 17 | 2232 (28) | 1.021 | 0.983 | 1.218 | 0.694 | 14.6 | 2186 | 2159 | 3218 | 1.4e-14 |
| anchor | 100_A2 | 0.291798 | 56 | 17 | 2260 (34) | 0.953 | 0.946 | 1.190 | 0.683 | 19.4 | 2372 | 2336 | 3309 | 1.4e-14 |
| anchor | 000_A1m | 0.269320 | 35 | 5 | 1603 (205) | 0.580 | 0.526 | 0.756 | 0.395 | 15.3 | 2761 | 2408 | 4058 | 1.1e-13 |
| anchor | 000_A1m | 0.345380 | 103 | 9 | 9133 (5004) | 2.069 | 1.949 | 2.249 | 1.459 | 23.0 | 4414 | 1996 | 6260 | 3.3e-13 |
| anchor | 100_A2 | 0.359900 | 108 | 30 | 4710 (204) | 2.558 | 2.427 | 2.747 | 1.785 | 29.9 | 1841 | 1761 | 2640 | 2.3e-13 |
| anchor | 100_A2 | 0.360000 | 108 | 30 | 4791 (198) | 2.243 | 1.762 | 2.438 | 1.715 | 24.2 | 2136 | 2048 | 2794 | 3.9e-12 |
| anchor | 110_A2 | 0.307330 | 68 | 27 | 2769 (17) | 1.382 | 1.171 | 1.552 | 0.927 | 29.4 | 2004 | 1992 | 2987 | 8.5e-14 |
| anchor | 110_A2 | 0.350000 | 108 | 41 | 4467 (89) | 1.832 | 1.659 | 1.953 | 1.772 | 33.3 | 2438 | 2390 | 2522 | 8.5e-14 |
| anchor | 111_A2 | 0.320360 | 64 | 20 | 1977 (27) | 0.915 | 0.876 | 1.042 | 0.955 | 22.2 | 2161 | 2132 | 2072 | 1.3e-13 |
| anchor | 111_A2 | 0.351030 | 160 | 40 | 7495 (386) | 3.852 | 3.429 | 4.034 | 3.551 | 27.2 | 1946 | 1845 | 2112 | 0.0e+00 |
| anchor | 200_A2 | 0.286640 | 47 | 16 | 1624 (16) | 0.676 | 0.640 | 0.794 | 0.720 | 20.3 | 2403 | 2378 | 2256 | 1.4e-14 |
| anchor | 200_A2 | 0.347380 | 121 | 32 | 5335 (211) | 2.091 | 1.922 | 2.240 | 2.054 | 28.5 | 2551 | 2450 | 2599 | 4.1e-10 |

Notes on the columns. "min of 3 runs" is the minimum of the per-run minima over the official run (5 repeats) and the two 21-repeat reruns; it is used where noise matters. Wall-clock noise on this VM is large: the first 5 anchors in the official run were up to 55 % slower in the median than in the two reruns (e.g. 100_A2 0.291556: 1.22 ms vs 0.83 and 0.79 ms), a transient of the host, not of the code. Minima are therefore the comparison metric. All 27 points: sign of $\det Q$ agrees with the oracle; max $|\Delta\ln|\det Q||=4.1\times10^{-10}$ (200_A2 0.34738; tolerance $10^{-6}$).

Summary of the speed-up: total/total with oracle projector min 1841, median 2551, max 8372; without oracle projector min 1761, median 2408, max 4176; oracle wall / fast wall (min of 3 runs) min 2072, median 2813, max 8118. The plan threshold was 20x (step 6), 150x at the largest point (step 7); the largest point $N=160$ (111_A2 0.351030) gives 1946x (total/total) and 1845x without the oracle projector. Warm fast total per point: 0.39 to 3.85 ms; cold minus warm at most 31.4 ms (the cold call builds the F tables, thresholds and projector; cold wall 11 to 34 ms, never seconds; the old 5.4 s cold projector of report 006 is gone).

## Stage breakdown (ms; official run, warm median over 5 repeats, cold = first call after `clear_caches()`)
**N=56 (100_A2, $E^*=0.291556$, $N_{\rm proj}=17$)** (anchor)

| stage | fast median | fast min | fast cold | oracle (1 call) |
|---|---|---|---|---|
| basis | 0.340 | 0.320 | 7.555 | 16.764 |
| F | 0.229 | 0.190 | 6.117 | 2180.813 |
| G | 0.054 | 0.051 | 0.053 | 1.785 |
| K2inv | 0.026 | 0.026 | 0.037 | 0.033 |
| F3 | 0.126 | 0.118 | 0.161 | 0.427 |
| F3inv (LU of $F_3/s$ in fast) | 0.033 | 0.031 | 0.037 | 0.293 |
| projector | 0.036 | 0.034 | 0.244 | 28.020 |
| K3 | 0.105 | 0.084 | 0.095 | 4.006 |
| projection | 0.044 | 0.039 | 0.043 | 0.126 |
| LU / slogdet | 0.008 | 0.007 | 0.007 | 0.090 |
| **total** | **1.021** | **0.983** | **14.350** | **2232.359** |

**N=160 (111_A2, $E^*=0.351030$, $N_{\rm proj}=40$)** (anchor)

| stage | fast median | fast min | fast cold | oracle (1 call) |
|---|---|---|---|---|
| basis | 0.231 | 0.222 | 7.341 | 18.874 |
| F | 0.507 | 0.357 | 16.185 | 6998.016 |
| G | 0.237 | 0.231 | 0.321 | 27.432 |
| K2inv | 0.079 | 0.077 | 0.123 | 0.150 |
| F3 | 1.293 | 1.262 | 1.234 | 12.205 |
| F3inv (LU of $F_3/s$ in fast) | 0.320 | 0.294 | 0.304 | 6.994 |
| projector | 0.058 | 0.058 | 0.431 | 385.860 |
| K3 | 0.413 | 0.405 | 0.653 | 44.385 |
| projection | 0.483 | 0.410 | 0.420 | 0.598 |
| LU / slogdet | 0.019 | 0.016 | 0.018 | 0.184 |
| **total** | **3.852** | **3.429** | **27.031** | **7494.699** |

At $N=160$ the largest fast stages are now $F_3=F/3-F H^{-1}F$ (1.29 ms, 34 %), the $Q$ projection (0.48 ms, 13 %), $F$ (0.51 ms, 13 %, the median is 42 % above the minimum of 0.36 ms, i.e. noisy) and $\mathcal K_{\rm df,3}$ (0.41 ms). The oracle is 93 % $F$ (7.0 s of 7.5 s). Source: `results/bench_m1_gate.json`, `points[*].fast.stages`, `.oracle.stages`.

## $\chi^2$ proxy (12 labelled L20 roots, cutoff 0.335, `coded_v2`, fast backend, `refine_root_info(kind="zero")`, 1 core, no caches on disk)
| quantity | value |
|---|---|
| roots refined / classified as zeros | 12 / 12 |
| function evaluations | 102 total (median 8.5 per root) |
| warm wall (median of 3 / min / max) | **0.122 s** / 0.117 s / 0.159 s (CPU/wall 0.992) |
| cold wall (first call, caches cleared) | 0.176 s |
| time per evaluation (warm median) | 1.19 ms |
| 4-thread request run (not the metric) | warm 0.115 s, cold 0.159 s, CPU/wall 0.994 (the fast path is single-threaded) |
| $\chi^2$ at the fixed baseline $\vec c$ (12 levels, D1 covariance, $N_{\rm cfg}=603$) | 53.8903 (identical to report 011/012; **not a fit**) |
| vs old direct FCN 1.4 s | 11.5x faster (using the warm median) |
| vs hot-window FCN 24 ms (needed GB-scale caches) | 5.1x slower (0.197 ratio) |

## Memory
Peak RSS of the whole benchmark process (fast + oracle + proxy), `resource.getrusage` in the wrapper: **93.4 MiB** (JSON `peak_rss_mb`: 47.3 MiB after the fast stage, 60.8 MiB after the oracle, 93.4 MiB after the proxy; these are cumulative maxima). No cache binaries are read or written. The report-012 figure of 9 MB for the $F$ tables of a full scan is consistent with this.

## Comparison with the history (minima where noise matters)
| point | oracle (006 / 010 / gate) | fast warm, step 6 (006) | step 7.4-7.6 (010) | step 7.7 (012) median / min | gate: median (official) / min of 3 runs |
|---|---|---|---|---|---|
| 100_A2 0.291556, $N=56$ | 2.68 s / n.a. / 2.23 s | 33 ms | 6.6 ms | 0.93 / 0.84 ms (012 report table); 0.79 / 0.75 ms (012 JSON) | 1.22 / 0.69 ms |
| 111_A2 0.351030, $N=160$ | 8.41 s / 8.13 s / 7.50 s | 103 ms | 13.8 ms | 3.89 / 3.71 ms (012 report table); 5.41 / 4.39 ms (012 JSON) | 4.03 / 3.55 ms |
| $\chi^2$ proxy warm (min) | | | | 0.140 (0.112) s; step 8: 0.49 to 0.53 s | 0.122 (0.117) s |
| $\chi^2$ proxy cold | | | | 0.173 s | 0.176 s |

Fast-path regression check against `results/bench_tester_012_fastonly_step7_7.json` (same harness, same points), per-point ratio of minima (min over my three gate runs / report-012 min): worst 1.053 (000_A1m 0.34538 and 200_A2 0.28664 at +5 %, 000_A1m 0.26932 at +4 %; within noise), best 0.68; 27 of 27 points are within the 10 % threshold. Physics output vs the same file: $\max|\Delta\ln|\det Q||=0$ (bit-identical), signs identical, $\chi^2$ proxy value identical. The audit fixes (finiteness guards, checked `.real()`) therefore cost nothing measurable.

Oracle vs history: the oracle takes 7.50 s at $N=160$ here (8.41 s in 006, 8.13 s in 010); machine-to-machine and run-to-run spread of this size (about 10 %) is expected. The oracle's projector build is 0.39 s at that point and up to 10.1 s on the 000_A1m rows (a single call each, cold projector cache), which is why the "with oracle projector" speed-ups reach 8372x on 000_A1m $N=135$ while the "without" speed-ups are 3136x.

## Caveats
1. The oracle was run once per point (protocol allows this); its noise is not measured, so speed-ups carry roughly the 10 % oracle spread plus the fast-path noise (±30 %, up to 2x transients on this VM at 1 ms scale). Order of magnitude is secure: all 27 points exceed 1700x.
2. The first official run and the two reruns differ by up to 55 % on individual median entries; minima agree within about 10 %. I report both and base comparisons on minima.
3. The 4-thread run is not the metric; the fast path has no threads, so only OpenBLAS/OpenMP availability changes, and the timings are unchanged within noise (111_A2 0.351030: 3.87 ms median; 100_A2 0.291556: 0.80 ms).
4. `build_flags` in the JSON says "see build.sh"; the full flag list is recorded above (the module cannot report its own flags).

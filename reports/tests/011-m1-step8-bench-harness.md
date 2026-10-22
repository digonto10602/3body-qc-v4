# Tester report 011: M1 step 8, benchmark harness `qc3.bench` and the $\chi^2$ proxy

- **Date:** 2026-10-22
- **Git commit:** `d718530` (the code of step 8 is commit `53e7c31`; `d718530` differs from it only in `graphify-out/`). Working tree clean apart from graphify files and my new `results/` files.
- **Prompt in force:** `prompts/002-m1-step7-speed-replan-projector-basis-k3-g.md` (§2.5, step-8 note); step 8 itself is in `prompts/001` §3 row 8, §5.7, §7.
- **Mode:** step (not the formal benchmark; no report in `reports/benchmarks/`).
- **Related:** `reports/tests/010-m1-step7-4-5-6-kernels-tail.md`, `reports/tests/007-m1-step9-pwave-shift.md`.
- **Machine:** Intel Xeon @ 2.10 GHz, 4 cores online, Python 3.13.16, NumPy 2.5.3, SciPy 1.18.1, g++ 13.3.0.

## Verdict: PASS (with minor notes)

## Commands run

```
./build.sh core                                                   # _core now baked at d718530 (was 338a0e8)
python3 -m pytest -q tests --ignore=tests/acceptance              # 67 passed in 37 s
python3 -m pytest -q tests/acceptance -m "not slow"               # 38 passed, 9 deselected, 39 s
.claude/bin/acceptance-check                                      # PASS (matches acceptance-v4)
PYTHONPATH=src/python python3 -m qc3.bench --out results/bench_tester_smoke_step8.json --no-anchors --oracle-repeats 1 --repeats 5 --proxy-repeats 3
PYTHONPATH=src/python python3 -m qc3.bench --out results/bench_tester_fastonly_step8.json --no-oracle --repeats 7 --proxy-repeats 5
PYTHONPATH=src/python python3 -m qc3.bench --out results/bench_tester_threads4_smoke_step8.json --no-oracle --no-anchors --threads 4 --repeats 3 --proxy-repeats 2
python3 <scratchpad>/indep_chi2.py   -> results/tester_011_chi2_independent.json   # own NumPy chi^2 from the 12 roots
python3 <scratchpad>/count.py                                      # counts real evaluate_point calls inside the proxy
```

All timing runs are deterministic (fixed energies and windows, no random numbers; seed: none).

## Results

| # | check | expected | got | tolerance | result |
|---|---|---|---|---|---|
| 1a | unit tests | 67 pass | 67 passed | all | PASS |
| 1b | acceptance, fast part | all green | 38 passed, 9 slow deselected | all | PASS |
| 1c | `acceptance-check` | PASS | PASS (matches `acceptance-v4`) | exact | PASS |
| 1d | `_core` rebuilt | baked commit = HEAD | `core_built_at_commit` = `d718530`, JSON records it | exact | PASS |
| 2a | thread env set before numpy import, recorded | `numpy_imported_before_env=false`, all 5 vars = 1 | false; OMP/OPENBLAS/MKL/NUMEXPR/VECLIB = 1; proxy CPU/wall 0.985–0.99 | CPU/wall $\le1.1$ | PASS |
| 2b | core pinned | affinity before $\to$ after | `[0,1,2,3]` $\to$ `[0]` (`os.sched_setaffinity`, the `taskset -c 0` equivalent); `--threads 4` run: `[0..3]` and env = 4 | exact | PASS |
| 2c | points | 15 grid + 12 anchors = 27 | 27 points; the 12 anchors equal the locked anchor file in (irrep, $E$, $N$, $N_{\rm proj}$) | exact | PASS |
| 2d | warm-up / repeats / medians / min / max / n | 1 / 5; median, min, max per stage | fast: 1 warm-up, 5 (or 7) repeats; oracle 1 repeat, 0 warm-up (protocol allows "once"); all 11 stages carry median/min/max/n | — | PASS |
| 2e | cold vs warm | one cold call after `clear_caches` | present for fast; cold $-$ warm $\le$ 11.9 ms (smoke), 14.1 ms (fast-only run); report 010: $\le$ 24.7 ms | — | PASS |
| 2f | oracle optional | `--no-oracle` works | works; speed-up fields simply absent | — | PASS |
| 2g | JSON provenance | commit, dirty flag, `_core` commit, versions, machine, env, constants, $K_2^{-1}$ form | all present: `git_commit=d718530`, `git_tree_dirty=false`, `core_built_at_commit`, `core_version=0.5.0-step7.6`, numpy, scipy, cpu model, cores online/allowed, env + affinity, production constants + $\vec c$ + grid, `k2_pwave_form=q0_consistent` (timing sectors) and `k2_pwave_form_proxy=coded_v2`, seeds, peak RSS | present | PASS (note 1, 2) |
| 2h | numbers plausible vs report 010 | fast 3.3–13.8 ms, oracle 2–10 s, $N=160$ 581× | fast warm 1.9–13.6 ms (fast-only run, 27 points); anchor 111_A2 0.35103: 13.97 ms (010: 13.8 ms); oracle 1.0–18.1 s (the first oracle call per sector includes the projector build, up to 0.4 s of it listed in `projector`); speed-up total 367–2019 (median 716), without the oracle's projector 348–854 (median 701); $\vert\Delta\ln\vert\det Q\vert\vert$ fast vs oracle $\le4.7\times10^{-12}$ at the 15 grid points | order of magnitude | PASS |
| 3a | proxy pins `coded_v2`, uses `refine_root_info(kind="zero")`, no scan | as stated | code: `PROXY_K2_FORM="coded_v2"`, `refine_root_info(..., kind="zero")`; no `find_zeros`/`find_sign_changes` anywhere in `bench.py`; JSON `k2_pwave_form="coded_v2"`; my own run confirmed the sector form is `coded_v2` | — | PASS |
| 3b | 12 of 12 `kind=="zero"` | 12 | 12; $\Delta_h=-13.7\ldots-26.1$ (step 9: $-13.1\ldots-26.6$) | — | PASS |
| 3c | roots vs step-9 `coded_v2` roots (`results/k2_pwave_shift_L20.json`) | $\le10^{-11}$ | max $\vert\Delta E\vert=8.4\times10^{-15}$ | $10^{-11}$ | PASS |
| 3d | roots vs labelled `previous_model_zero` | $<5\times10^{-6}$ | max $4.35\times10^{-8}$ | $5\times10^{-6}$ | PASS |
| 3e | $\chi^2$ at the fixed baseline $\vec c$, own NumPy | 53.89 (builder) | independent: $r^{\!\top}C^{-1}r=53.890297$ (solve) $=53.890297$ (explicit inverse); bench: 53.890297; difference $1.3\times10^{-12}$ | $10^{-8}$ relative | PASS |
| 3f | residual$/\sigma$ range | $-1.77\ldots+1.21$ | $-1.774\ldots+1.213$ (12 values in `results/tester_011_chi2_independent.json`) | — | PASS |
| 3g | JSON states it is not a fit | `is_a_fit=false` + note | `is_a_fit=false`; note "evaluated at the label-generation vector C_BASELINE; no parameter is varied or minimised" | — | PASS |
| 3h | evaluation count | 12 roots $\times\approx$6–9 | 102 total (8–9 per root; includes the 3 classification calls). I counted real `fast.evaluate_point` calls with a wrapper: 102, no mismatch with `n_evaluations` on any root | exact | PASS |
| 3i | time per evaluation / proxy time | 0.5–3 s (plan) | warm 0.49–0.53 s (medians 0.503, 0.524, 0.526 s in three runs), cold 0.57–0.63 s, i.e. 4.8–5.1 ms per evaluation; 2.7× faster than the old direct FCN (1.4 s), 21× slower than the hot-window FCN (24 ms, which needed GBs of cache); CPU/wall 0.99 | — | PASS |
| 4a | silent exception catching | none in numerical code | the only `except` clauses in `bench.py` are `OSError` in `_cpu_model` and `_git` (metadata only) and `ImportError` for SciPy's version string; the proxy docstring says a failing root propagates, and `chi2_proxy` raises `RuntimeError` on a label mismatch | — | PASS (note 2) |
| 4b | no hard-coded numbers beyond protocol | — | protocol constants only (energies 0.2700/0.3100/0.3599, cutoff 0.335, reference times 1.4 s and 24 ms, grid $E_0$, $dE$ in the metadata block) | — | PASS (note 3) |

Raw numbers: `results/bench_tester_smoke_step8.json` (15 grid points, oracle + fast, proxy), `results/bench_tester_fastonly_step8.json` (27 points, fast only, proxy with 5 warm passes), `results/bench_tester_threads4_smoke_step8.json`, `results/tester_011_chi2_independent.json`. All carry commit `d718530`, parameters and the "no seed" statement.

## How the independent $\chi^2$ was computed

From `data/lattice/L20/*.jack` directly (no `qc3.lattice`): $y_c=(\sum x-x_c)/(N-1)$ with $N=603$, $E^*_c=\sqrt{y_c^2-P^2}$ with $P=2\pi|\mathbf d|/(\xi\,20)$, mean over $c$, $C_{ij}=\tfrac{N-1}{N}\sum_c\delta_{i,c}\delta_{j,c}$, keep levels with mean $E^*\le0.335$ (the same 12 labels as the windows), residual $r=E_{\rm root}-\bar E$ (D3: the sorted pairing is the identity for all 5 irreps, checked), $\chi^2=r^{\top}C^{-1}r$. Condition number of $C$: $1.81\times10^{3}$ (matches report 004). The lattice means and $\sigma_i$ stored in the proxy JSON agree with mine to $10^{-14}$ for all 12 levels.

## Notes (none blocking)

1. **Compiler and flags are not recorded**, only the string "see build.sh COMMON (-std=c++17 -O3 -fopenmp)". Protocol §7 asks for "compiler + flags". Guess: not introspectable from the module; easy to add (for example write the compiler version into the build step, or call `g++ --version` in `provenance()`).
2. `_git()` returns an empty string if `git` is missing, which would make `git_commit=""` and `git_tree_dirty=false` without any warning. It is metadata only and harmless in this repo, but a missing commit in a benchmark JSON should be loud.
3. `fast_threads_option: 1` and the grid numbers ($E_0=0.2631$, $dE=0.0969/19999$, $n=20000$) in the JSON are literals, not read from the code that uses them. Today they are correct.
4. The oracle is run with 0 warm-up calls, so its numbers contain the one-time projector build (0.4 s at $N=160$, more at some points). The harness records `projector` per stage and a second speed-up without it, so the information is there; the formal benchmark may want `--oracle-warmup 1` (not exposed in the CLI; only `run()` has it).
5. Timing noise on this machine: in one run the $N=160$ grid point (111_A2, 0.3599) read 21.9 ms warm and in the next run 12.8 ms; report 010 saw the same kind of noise (F stage 8.8–12 ms). For the formal benchmark use several repeats and compare minima.
6. A `--threads 4` run changes only the BLAS threads and the affinity; the fast path itself still runs one thread (`fast_threads_option: 1`), so that run is only a thread-safety smoke test.

## Counts

Passed: 67 unit + 38 fast acceptance + 1 acceptance-check + all rows above. Failed: 0. Pending: 0 (the 9 slow acceptance tests were not selected in this mode and are untouched by step 8). New failures vs report 010: none.

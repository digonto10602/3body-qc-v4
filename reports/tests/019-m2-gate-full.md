# Tester report 019: M2 gate, mode full

- Date: 2026-10-26. Commit `1065788` (tree clean; `_core` rebuilt at this commit, `qc3._core.git_commit()` = `1065788`). Plan in force: `prompts/003-m2-label-free-roots-tracking-fit.md`. Acceptance tag `acceptance-v7`.
- Related: `reports/tests/018-m2-step6-fits.md`. Raw numbers: `results/tester_019_hesse_repro.json`; scripts: `results/tester_019_scripts/`.
- Machine: Xeon @ 2.80 GHz, 4 cores (report 013 ran on a 2.1 GHz Xeon), OpenBLAS 0.3.34 (Haswell kernel), numpy 2.5.3, iminuit 2.33.0. All runs `OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1` unless stated.

## VERDICT: FAIL (78 of 79 acceptance tests pass; 1 fails, an environment-dependent anchor tolerance, not a code change)

## Commands and results
| check | expected | got | tolerance | result |
|---|---|---|---|---|
| `rm _core*.so; ./build.sh core`; `_core.git_commit()` vs HEAD | equal | 130 s build; `1065788` = HEAD | exact | PASS |
| `python3 -m pytest -q tests/acceptance -p no:cacheprovider -rA` (incl. slow) | 79/79 | **78 passed, 1 failed** in 1044.6 s | | **FAIL** |
| `python3 -m pytest -q tests --ignore=tests/acceptance` | all green | 154 passed (103 s) | | PASS |
| `.claude/bin/acceptance-check` | PASS | PASS (matches acceptance-v7) | | PASS |
| `scripts/m2/score_labels.py` (needs `PYTHONPATH=src/python`) | exit 0, FP=FN=0 | exit 0; TP/FP/FN/TN = 31/0/0/176 over 7 sectors; 53046 evals, 120.6 s | FP=FN=0 | PASS |
| `git diff 191530d -- src/qc3_core/oracle_point.hpp src/core src/fitter src/cachegen_cpu` | empty | empty | | PASS |
| graphify `affected`/`query` on `fit_cutoff` | covered | fit, error-gap, tracking, roots tests all in the acceptance and unit runs | | PASS |

Note: running `score_labels.py` rewrites the tracked `results/m2_label_gate.json` (provenance fields only: commit `da1c6fb`->`1065788`, built-at commit, timings). I reverted it with `git checkout`; the tree is clean.

### The one failure
`tests/acceptance/m1/test_m1_oracle.py::test_oracle_anchors_all_points` (slow). Point 200_A2, $E^*=0.34738$, $\vec c=0$:
oracle $\ln|\det Q|=-162.4167393752109$, anchor $-162.41673937565963$, $|\Delta|=4.49\times10^{-10}$ vs tolerance $10^{-10}$. The sign and the dimensions agree.

- Deterministic: identical with 1 thread, 4 threads, and with the thread variables unset.
- Cause established, not just guessed: it depends on the OpenBLAS CPU kernel. With `OPENBLAS_CORETYPE` = Haswell / SkylakeX / Zen the gap is $4.49\times10^{-10}$; with Sandybridge or Nehalem it is $2.27\times10^{-10}$ (still above $10^{-10}$). The anchors were produced on another CPU, and the oracle's `np.linalg.slogdet` is the LAPACK call. Report 006 had $|\Delta|=0$ exactly at all 12 anchors, on the 2.1 GHz machine.
- This is the risk STATUS recorded from tester 002 ("anchor tolerances are below the round-off floor at ill-conditioned points; may fail on a different CPU or BLAS"). $\kappa(F_3)=4.8\times10^{8}$ at this point. For scale, the oracle's independent Eigen-LU value is $-162.41673937466325$ (gap $10^{-9}$ from the anchor).
- The oracle source is byte-identical to `191530d`, so this is not a code regression. It is nevertheless a failing locked test on the current machine. Changing the tolerance is the planner's decision (e.g. a kappa-scaled bound). I did not touch it.

## Hesse reproducibility (Tables A and B)

Commands: `python3 -m qc3.fit --cutoffs <c> --form q0_consistent [--fix KB KE] --no-verify --out ...` for $c\in\{0.335,0.345\}$ and A (4 parameters) and B (KB=KE=0), twice each, each in a fresh process (`results/tester_019_scripts/run.sh`).

**Result: deterministic on a fixed machine.** All 8 runs are bit-identical to each other and to the committed JSONs at `9911237` (chi2, edm, $\vec c$, Hesse errors, FD-Hessian errors, error gap, flags, nfcn). Example, B/0.345: $\chi^2=55.2232755983$, EDM $2.33\times10^{-10}$, Hesse $(1.3586, 2.1207)\times10^5$, FD $(1.6287, 2.5356)\times10^5$, ratio 1.199/1.196 (under the 1.2 threshold, so no disagreement flag), `has_accurate_covar` True, nfcn 56. A/0.345: $\chi^2=35.5421543066$, EDM $9.05\times10^{-10}$, ratios (0.82, 1.29, 2.62, 2.95), all four flagged, `error_exceeds_value:Kiso1`. A/0.335 (`hesse_not_accurate`) and B/0.335 (no flags) are also unchanged.

**Why the flags changed against the previous JSONs.** The previous JSONs came from `e9d0215` (WIP, 2026-10-24 23:08); `fit.py` has since been developed. Between `95ebb88` and HEAD only the error-gap diagnostics were added (no change to the fitting path). A/0.335 and B/0.335 are bit-identical to the old JSONs. At the two 0.345 fits the minimum moved by about $10^{-10}$ in $\chi^2$ (A: 35.5421543069 -> ...066; B: 55.2232755979 -> ...983) and EDM $4.6\times10^{-9}\to9.0\times10^{-10}$ (A), $8.0\times10^{-9}\to2.3\times10^{-10}$ (B). Hesse then changed: B/0.345 Hesse $(5.77, 9.07)\times10^5\to(1.36, 2.12)\times10^5$; A/0.345 $K_0$ error $1.16\times10^5\to2.69\times10^5$.

**What drives it: the Hesse error is a function of where Migrad stops, because $\chi^2(\vec c)$ is piecewise smooth (poles and dimension jumps within $1\sigma$; reports 018 N2/N3).** Test: same fit, start vector $\vec c_{\rm start}$ perturbed by a relative amount $\epsilon$ (seeded `default_rng(seed)`, fresh trackers each time, `perturb.py`). Hesse error of $K_0$ in units of $10^3$, and the FD/Hesse ratio:

| case | $\epsilon$, seed | $\chi^2$ | EDM | Hesse $K_0$ / $K_1$ ($10^3$) | FD/Hesse | flags |
|---|---|---|---|---|---|---|
| B/0.345 | 0 | 55.22327560 | 2.3e-10 | 135.9 / 212.1 | 1.20 | none |
| B/0.345 | 1e-10, 1 | 55.22327560 | 1.7e-9 | 353.8 / 552.6 | 0.46 | corr>=0.99 |
| B/0.345 | 1e-10, 2 | 55.22327560 | 7.3e-10 | 219.0 / 343.1 | 0.74 | none |
| B/0.345 | 1e-8, 1 | 55.22327560 | 1.2e-11 | 230.8 / 359.4 | 0.71 | none |
| B/0.345 | 1e-8, 2 | 55.22327560 | 2.5e-9 | 204.4 / 316.2 | 0.80 | none |
| B/0.345 | 1e-6, 1 | 55.22327560 | 2.1e-9 | 108.0 / 167.6 | 1.51 | none |
| B/0.335 | all 6 | 36.58654145 | 7.5e-8 | 157.6 to 168.1 / 239.1 to 255.3 | 0.94 to 1.00 | none |
| A/0.345 | 0 | 35.54215431 | 9.0e-10 | 269 / 386 | 0.82, 1.29, 2.62, 2.95 | 4 disagree |
| A/0.345 | 1e-10, 1 | 35.54215431 | 8.6e-10 | 188 / 354 | 1.17, 1.40, 2.28, 2.26 | 3 disagree |
| A/0.345 | 1e-10, 2 | **3876.14** | 1.1e-8 | 204 / 925 | 0.29, nan, 0.07, nan | different local minimum |
| A/0.345 | 1e-8, 1 | **3876.14** | 1.2e-9 | 272 / 283 | 0.21, nan, 0.09, nan | different local minimum, `hesse_not_accurate` |
| A/0.345 | 1e-8, 2 | 35.54215431 | 1.1e-9 | 206 / 234 | 1.07, 2.12, 1.93, 1.95 | 3 disagree |
| A/0.345 | 1e-6, 1 | 35.54215431 | 1.1e-8 | 97 / 369 | 2.28, 1.35, 1.74, 1.65 | 4 disagree |

Findings:
- B/0.335 is stable (Hesse spread 7 % at $\epsilon\le10^{-8}$; FD/Hesse 0.94 to 1.00).
- B/0.345: $\chi^2$ is the same to $6\times10^{-10}$ in all six runs, yet the Hesse error varies by a factor 3.3 ($1.08\times10^5$ to $3.54\times10^5$) and FD/Hesse from 0.46 to 1.51. Hesse at 0.345 is not reproducible under $10^{-10}$ changes of the start; the committed value (ratio 1.20, no flag) is one draw from that spread. The FD-Hessian is the stable one: its $K_0$ is $1.63\times10^5$ (committed), and the six runs give $(1.63, 1.63, 1.63, 1.64, 1.63, 1.63)\times10^5$ from $\mathrm{fd}=\mathrm{gap}\times\mathrm{Hesse}$ (e.g. $0.46\times353.8=162.7$).
- A/0.345: FD $K_0$ stays at about $2.2\times10^5$, Hesse spans $0.97$ to $2.72\times10^5$. Two of five perturbed starts end in a different local minimum at $\chi^2=3876.14$ (a wall: penalised or missing level); the unperturbed start is in the good minimum, so the committed table is fine, but the minimum is not robust.
- Mechanism (a guess supported by the above): Minuit2's Hesse takes finite-difference steps from the Migrad state (current step sizes and EDM stopping point). Walls within $1\sigma$ make the second derivative depend on the step, so a $10^{-10}$ change in the stopping point changes the Hesse matrix. Nothing is random: no RNG or threading is involved (bit-identical repeats, thread count 1 vs 4 not needed). I did not test a different machine or BLAS for the fit.

## Return notes
- The `error_gap` flag is therefore useful, but its pass/fail at 0.345 is not stable under tiny start changes; B/0.345 "no longer flagged" is within the noise of this effect. The FD-Hessian errors should be quoted beside Hesse for the 0.345 fits.
- Recommend: planner decides the tolerance for `test_oracle_anchors_all_points` (kappa-scaled, or regenerate anchors); until then the gate on this machine is FAIL by one locked test.

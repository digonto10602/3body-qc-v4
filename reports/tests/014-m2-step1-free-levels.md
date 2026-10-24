# 014 M2 step 1: interface stubs + `qc3.free_levels` (tester, mode step)

Date 2026-10-24. Commit `e1016fc` (tree clean). Plan in force: `prompts/003` (§2.1, §3 step 1, §4, §5), tag `acceptance-v7`.
Raw numbers: `results/m2_step1_free_levels_bruteforce.json`. Independent script: scratchpad `bf.py` (own brute force, $|n_i|\le 8$ per component, all ordered kaon pairs filtered to $n_1\le n_2$, $E^*=\sqrt{E_{\rm lab}^2-P^2}$, constants from `qc3.constants.PRODUCTION`, merge $10^{-11}$).

## Commands
`pytest -q tests/acceptance/m2/test_m2_free_levels.py`; `pytest -q tests/acceptance/m2 --ignore=.../test_m2_free_levels.py`; `pytest -q tests --ignore=tests/acceptance`; `pytest -q tests/acceptance/m1 -m "not slow"`; `.claude/bin/acceptance-check`; `git diff 191530d -- src/qc3_core/oracle_point.hpp src/core src/fitter src/cachegen_cpu`; `pytest --collect-only -q tests/acceptance`.

## Results
| check | expected | got | tol | |
|---|---|---|---|---|
| m2 free-level tests | 11 green | 11 passed | - | PASS |
| rest of m2 | only `NotImplementedError` | 21 failed, all `NotImplementedError` (rootscan step 2 x10, tracking step 4 x6, fit step 5 x3, fit step 6 x2) | - | PASS (pending, not failures) |
| acceptance collect | 79 | 79 | - | PASS |
| m1 non-slow | green | 38 passed, 9 deselected | - | PASS |
| unit tests (all, incl. new `test_m2_step1.py`) | green | 91 passed (75 s) | - | PASS |
| `acceptance-check` | PASS | PASS (matches acceptance-v7) | - | PASS |
| oracle/core/fitter/cachegen_cpu diff vs `191530d` | empty | empty | - | PASS |
| independent free levels, counts L20 (000/100/110/111/200) | 2/3/4/6/5 | 2/3/4/6/5 | exact | PASS |
| same, L24 | 3/6/11/6/9 | 3/6/11/6/9 | exact | PASS |
| $\max\lvert E_{\rm prod}-E_{\rm mine}\rvert$, all 10 (L, irrep) | $\le10^{-12}$ | $5.6\times10^{-17}$ (L24 000_A1m; 0 elsewhere) | $10^{-12}$ | PASS |
| multiplicities and state sets, all levels | equal | equal for all 10 sectors | exact | PASS |
| `110_A2` L20 level | 0.319468828, mult 2 | 0.31946882819441086, mult 2 | $5\times10^{-10}$ | PASS |
| stubs name the right step | step 2/4/5/6 per §3 | rootscan->2, tracking->4, `QCChi2`/`make_trackers`->5, `fit_cutoff`/`verify_fit`/`cutoff_scan`/`save_results`/`main`->6 | - | PASS |
| silent exception catching | none | none; only `fast.evaluate_point` catches `ValueError`, re-raises `SingularQError` iff message has "Q is singular", else bare `raise` | - | PASS |

Multiplicity lists (mine = production) for the record: L20 000: [3,6]; 100: [1,1,4]; 110: [1,1,1,2]; 111: [1,3,1,3,3,3]; 200: [1,1,1,1,2]. L24 110: [1,1,1,2,2,2,2,2,4,16,9].

## Problems
None. Notes: (a) `free_levels` reports the mean of the merged energies, the test helper the first member; they agree to 1 ulp here. (b) the plan row says "collects 76"; 79 are now collected (v7 added 3 tests), consistent with STATUS.

VERDICT: PASS (11/11 free-level, 38 m1 non-slow, 91 unit; 21 m2 pending on `NotImplementedError`; 0 failures).

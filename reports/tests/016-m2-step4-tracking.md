# Tester report 016: M2 step 4 (qc3.tracking)

- Date: 2026-10-24. Commit: a93e308 (merge 6204201 + graph refresh). Plan in force: prompts/003 (§2.3, §2.5, §2.6 D9, §3 step 4, §4, §5). Tag acceptance-v7. Mode: step.
- `_core` was stale (c0334de); rebuilt with `./build.sh core`, now reports a93e308.
- Related: reports/tests/014, 015. Raw numbers: `results/m2_step4_tracking_jacobian_check.json`, `results/m2_step4_tracking_cost_check.json` (seed 20261024).
- Commands: `pytest -q tests/acceptance/m2/test_m2_tracking.py` (incl. slow); `pytest -q tests --ignore=tests/acceptance`; `pytest -q tests/acceptance/m2 -m "not slow"`; `pytest -q tests/acceptance/m2/test_m2_rootscan.py tests/acceptance/m2/test_m2_free_levels.py` (incl. slow); `pytest -q tests/acceptance/m2/test_m2_fit.py`; `pytest -q tests/acceptance/m1 -m "not slow"`; `.claude/bin/acceptance-check`; scratch scripts (chk.py, poles.py, cost.py in the session scratchpad).

## Results

| check | expected | got | tolerance | verdict |
|---|---|---|---|---|
| m2 test_m2_tracking.py (6, incl. slow, 54 s) | 6 pass | 6 pass | - | PASS |
| m2 rootscan + free_levels (incl. slow, 3m17s) | all pass | 21 pass | - | PASS |
| m2 test_m2_fit.py | NotImplementedError only | 5 fail, all `NotImplementedError` (make_trackers step 5 x3, cutoff_scan step 6 x2) | - | PENDING (5) |
| unit tests (tests/ minus acceptance) | green, ~122+ | 131 pass | - | PASS |
| m1 acceptance, non-slow | green | 38 pass | - | PASS |
| acceptance-check | PASS | PASS (matches acceptance-v7) | - | PASS |
| Jacobian eq.(1) vs own central FD of Brent roots, 000_A1m, 2 random c (3 % scatter on baseline, coded_v2) | rel <= 1e-5 | production 6.0e-8, 5.6e-8; own exact-root-null-vector version 5.5e-8, 5.1e-8 | 1e-5 | PASS |
| same, 100_A2, 2 random c | rel <= 1e-5 | production 3.5e-6, 3.4e-6; own version 3.9e-8, 3.8e-8 | 1e-5 | PASS |
| slope sign $u^TQ'u$ | > 0 | +0.0074 (000_A1m), +0.0082 (100_A2); production slope equals mine to 4 digits | - | PASS |
| null vector: $\lambda_{\min}(Q(E_0))$, next eigenvalue | ~0 vs gap | $\le 10^{-11}$ vs $2.7\times10^{-5}$ (000_A1m), $3.1\times10^{-5}$ (100_A2) | - | PASS |
| $Q$ linear in $c$: $\tilde B_i=\partial Q/\partial c_i$ | equal | max rel diff $7\times10^{-12}$ to $6\times10^{-11}$ (000_A1m), $\sim10^{-7}$ (100_A2, round-off at $\kappa$ large) | - | PASS |
| poles independent of c (100_A2, 7 poles, 20 % random c vs baseline) | same | max shift $1.5\times10^{-11}$ (scan resolution); zeros move 0.2915561 -> 0.2915579 | - | PASS |
| tracked vs fresh scan, 5 sectors x 3 updates (q0_consistent, E_max 0.335, 5/15/30 % c changes) | <= 1e-9 (builder claim <= 1.1e-12) | max $2.4\times10^{-12}$ (111_A2 / 200_A2 order 1e-12) | 1e-9 | PASS (builder's 1.1e-12 claim is not reproduced: 2.4e-12 seen; far inside 1e-9) |
| evaluations per zero per substep | 5-7 | 5.0 - 6.7 (incl. the 2 Jacobian points) | <= 12 | PASS |
| fallbacks (local / full rescans) in those 15 updates | none expected | 0 / 0 | - | PASS |

Note on the 100_A2 production-vs-FD value (3.5e-6, vs 6e-8 for 000_A1m): it comes from taking $u$ from $\bar Q=\tfrac12[Q(E_0+h)+Q(E_0-h)]$ rather than from $Q(E_0)$ (my version, 4e-8). Within the 1e-5 requirement and the plan's quoted 1.8e-6 order; for a nearby pole the $O(h^2)$ shift grows. It only affects the predictor, not the root.

## Code reading (robustness and loudness)
- `_check_loud` raises `AmbiguousRootError` for any `ambiguous` or `even_zero` candidate in full scans, rescans and local rescans; the `even_zero` at a free level is classed `free_level` by `scan_sector` and kept apart in `free_level_solutions` (D8; counts 0/2/4/3/3 for 000_A1m/100_A2/110_A2/111_A2/200_A2 at the baseline with q0_consistent, never tracked).
- `rescan()` raises `RootSetChangedError` on count mismatch or > 1e-9. Merged roots raise `RootSetChangedError`. Fallbacks (secant not converged, pole crossed, dimension jump, predictor safeguard, non-distinct) are counted, not silent; no path drops a root without a rescan. Unit tests cover them (12 in test_tracking.py).
- D9: `wrong_sign_roots()` flags `slope < 0` and never excludes.
- Minor note (not a failure): the constructor scans only $[E_{\min},E_{\max}]$, so a zero that starts above $E_{\max}$ and moves down into the range is only noticed by `rescan()`/the fit's verification (`RootSetChangedError`), not by `update()`. A full rescan inside `update()` likewise drops parked roots. By design in the plan (parking applies to roots leaving upward); the fit step should call `rescan()` as planned.

## Verdict
PASS. 0 failures; 5 PENDING (m2 fit, steps 5-6). No new failures vs report 015.

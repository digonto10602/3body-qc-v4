# Test report 015: M2 step 2, `qc3.rootscan.scan_sector`

Date 2026-10-24. Commit `35cc1a8`. Plan in force: `prompts/003-m2-label-free-roots-tracking-fit.md` (sections 2.1, 2.2 incl. 4b, 2.6 D9, 3 step 2, 4 v7 amendment, 5), tag `acceptance-v7`. Mode: step. Previous: `reports/tests/014`.
Commands (1 core, `OMP_NUM_THREADS=OPENBLAS_NUM_THREADS=1`): `pytest -q tests/acceptance/m2/test_m2_rootscan.py`; `pytest -q tests/acceptance/m2 --deselect ...test_m2_rootscan.py`; `pytest -q tests --ignore=tests/acceptance`; `pytest -q tests/acceptance/m1 -m "not slow"`; `.claude/bin/acceptance-check`; `git diff 191530d -- src/qc3_core/oracle_point.hpp src/core src/fitter src/cachegen_cpu`.
Own scripts (scratchpad, not in repo): independent label scoring, pre-snap free-level offsets, 110_A2 pair at fine resolution. Raw numbers: `results/m2_step2_tester_gate.json`.

## Verdict: PASS

| check | expected | got | tolerance | result |
|---|---|---|---|---|
| `test_m2_rootscan.py` (10 tests, slow incl.) | all pass | 10 passed, 178 s | - | PASS |
| other m2 tests | only NotImplementedError (steps 4-6) | 11 failed, all NotImplementedError (tracking 6, make_trackers 3, cutoff_scan 2); 11 passed (free_levels) | - | PENDING 11, 0 FAIL |
| unit tests (`tests --ignore=tests/acceptance`) | green | 110 passed | - | PASS |
| m1 acceptance, non-slow | green | 38 passed | - | PASS |
| `acceptance-check` | PASS | PASS (acceptance-v7) | - | PASS |
| oracle/core/fitter/cachegen_cpu diff vs `191530d` | empty | empty | - | PASS |
| own gate, coded_v2, C_BASELINE, slack 1/4 row | FP = FN = 0, no unlabelled zero | see table below, 0 / 0 / 0 in 7/7 | - | PASS |
| 110_A2 pair at C_LEGACY_0335 | zero 0.31595530708897, pole 0.31594715696 | zero 0.3159553070889677 (diff 7e-15), pole 0.31594715696866 (diff 1e-11) | 1e-9 | PASS |
| zero slope sign (independent, h = 1e-7, 1e-8, 1e-9) | > 0 | +8.349e-3 (all h); scan (h = 1e-6): +8.41e-3 | sign | PASS |

Own gate (labels read from the CSVs, `rootscan.LABEL_FILES` paths, own scoring):

| sector | TP | TN | FP | FN | unlabelled zeros | free_level | poles | evals | s (1 core) |
|---|---|---|---|---|---|---|---|---|---|
| L20 000_A1m | 3 | 7 | 0 | 0 | 0 | 2 | 7 | 5917 | 6.5 |
| L20 100_A2 | 2 | 23 | 0 | 0 | 0 | 3 | 23 | 7205 | 9.8 |
| L20 110_A2 | 4 | 32 | 0 | 0 | 0 | 4 | 32 | 8019 | 11.2 |
| L20 111_A2 | 5 | 25 | 0 | 0 | 0 | 6 | 25 | 7776 | 11.3 |
| L20 200_A2 | 5 | 30 | 0 | 0 | 0 | 5 | 30 | 8033 | 13.9 |
| L24 000_A1m | 3 | 17 | 0 | 0 | 0 | 3 | 17 | 6740 | 16.1 |
| L24 100_A2 | 9 | 42 | 0 | 0 | 0 | 5 | 43 | 9337 | 28.6 |

207 scored + 1 excluded = 208 brackets. No `even_zero` or `ambiguous` candidate remains in any sector. Free-level counts equal `qc3.free_levels` counts (L20 2/3/4/6/5).
Bracket 36 of L24/100_A2: excluded only through `LABEL_EXCLUDE` in the locked `m2_helpers.py` and the `exclude_ids` argument of `score_labels` (plan section 3 step 3 / acceptance table), not by code in the scan. The scan contains no special case. For the record: its human label is `false`, and the scan finds zero zeros there and one pole (0.3407707614, $\Delta_h=+34.7$), so it would have scored TN.

## Judgement of the builder's choices
(a) Snap of E to the free level: justified. I recorded E before the snap for all 20 L20 free-level candidates (and the same code path for L24): every one is an even-order (order 2 or 4), no-sign-change candidate located within $|E-E_f|\le9.8\times10^{-10}$ of the exact free level ($3\times10^{-10}$ for 110_A2 quartic), i.e. 100x inside the D8 tolerance $10^{-7}$ of section 2.1/2.2 item 7. The snap needs `kind in {zero, even_zero, ambiguous}` and $|E-E_f|\le10^{-7}$, sets kind = `free_level` and `free_level` = $E_f$ in the same step. It is the exact solution of arXiv:1901.07095 section 4.4.3; order, slope and $\Delta_h$ are kept. The plan (item 7) specifies the kind and `free_level` field but does not say that `E` is overwritten; this is a harmless extension, not a special case for the test. One candidate (200_A2, order 2) is `ambiguous` before the snap (pole-polluted order estimate) and is relabelled by the same rule. Note: a genuine simple zero within $10^{-7}$ of a free level would also be relabelled; this is the plan's rule (margin $\ge10^3$ to the nearest labelled zero).
(b) Probe centre at $E_f+2.5\times10^{-9}$: sound. $F$ has a pole exactly at $E_f$ (fast path raises NONFINITE_COMPONENT). An order-$m$ zero gives $\Delta_h\approx m\ln(2.5\times10^{-9}/10^{-6})=-6m\le-12$ against $D=4$; all free-level solutions were found and flagged ($\Delta_h$ from $-12.8$ to $-34$).
(c) Rescaled variable $t\in[0,1]$ for `minimize_scalar(bounded)`: correct. scipy's tolerance $\sqrt\epsilon|x|+x_{\rm atol}/3\approx4\times10^{-9}$ at $x=0.28$ would otherwise defeat `xatol`=1e-12; in $t$ it is $\sim6\times10^{-16}$. Non-convergence raises.
(d) Zero $3.8\times10^{-6}$ below a pole is `ambiguous` (order estimate 1.41, polluted by the pole): consistent with section 2.2 items 5, 7 (stated in the plan, h2 = 3e-6 chosen for exactly this). It is never a silent zero, but it is dropped from `zeros()`. Real consequence for steps 4-6: the tracker/fit must treat any `ambiguous` candidate as a loud event (section 2.6), otherwise a genuine weakly coupled level at $s\lesssim4\times10^{-6}$ from its pole would vanish from D3. At the baseline and the 0.335 legacy vector none occurs. The resolved case ($s=8.15\times10^{-6}$) is correct.

## Cost
Evaluations 5917-9337 per sector, 6.5-28.6 s (L20 000_A1m 6.5 s, plan bound 15 s; L24 100_A2 28.6 s has no time bound). The 10 000 cap exists only as the test assertion `n_evaluations <= MAX_EVALUATIONS` in the locked test; `scan_sector` has no cap and no loud branch for it (it simply runs on). L24 100_A2 uses 9337 (93 %). At other parameter vectors (fit steps, more candidates) the test cap can be exceeded; cost is about 3 evaluations per candidate plus fine rescans. Risk is only for a test at a different vector, not for correctness. Suggest the planner note that the cap is a test bound only (a guess: raise it only via the planner for L24 or a tracker rescan benchmark).

## Problems
None blocking. Notes: item (d) above; E-overwrite in the free-level snap is not in section 2.2 (suggest a one-line doc addition); no timing was measured for `test_fast_single_thread_cpu_time`. Graphify `affected scan_sector`: no dependents (new module); covered by the tests above.

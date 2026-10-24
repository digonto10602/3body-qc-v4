# M2 step 5 tester report: qc3.fit.QCChi2, make_trackers, lattice.cov_eig_pinv, tracker fallback

- Date: 2026-10-24. Commit: `15d351d` (tree clean apart from my two `results/m2_step5_tester_*.json` files). Mode: step.
- Prompt in force: `prompts/003-m2-label-free-roots-tracking-fit.md` §2.4-§2.6, §3 step 5, §4, §5. Acceptance tag `acceptance-v7`.
- Raw numbers: `results/m2_step5_tester_chi2.json`, `..._grad.json`, `..._cost.json` (seed 20261024). Scripts (scratchpad, not in repo): `indep_chi2.py`, `grad.py`, `pen.py`, `cost.py`.

## Commands
```
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 python3 -m pytest -q tests/acceptance/m2 (minus test_m2_fit.py)   # 27 passed
python3 -m pytest -q tests/acceptance/m2/test_m2_fit.py -rA                                            # 3 passed, 2 failed (step 6)
python3 -m pytest -q tests --ignore=tests/acceptance                                                   # 143 passed
python3 -m pytest -q tests/acceptance/m1 -m "not slow"                                                 # 38 passed, 9 deselected
.claude/bin/acceptance-check                                                                          # PASS (matches acceptance-v7)
git diff 191530d -- src/qc3_core/oracle_point.hpp src/core src/fitter src/cachegen_cpu | wc -l          # 0
```

## Results
| check | expected | got | tolerance | verdict |
|---|---|---|---|---|
| acceptance m2 (30 of 32 tests) | all green except 2 cutoff_scan | 30 passed incl. `test_qcchi2_at_baseline_...`, both legacy/v32w code-equivalence tests; 2 failed with `NotImplementedError: qc3.fit.cutoff_scan` (step 6) | - | PASS (2 PENDING) |
| unit tests | ~143 | 143 passed | - | PASS |
| m1 non-slow acceptance | green | 38 passed | - | PASS |
| acceptance-check | PASS | PASS (v7) | - | PASS |
| oracle / core / fitter / cachegen_cpu diff vs `191530d` | empty | 0 lines | - | PASS |
| independent chi2, coded_v2, C_BASELINE, 0.335, 12 levels | 53.890297 | 53.890296754032 (QCChi2 identical) | rel 1e-8 | PASS (rel 0) |
| independent chi2, q0_consistent, C_BASELINE | 47.9426 | 47.942625940431 vs QCChi2 47.942625940432 | 1e-8 | PASS (8e-15) |
| coded_v2, 2 random c (±3 %) | agree | rel 7.8e-10, 2.9e-9 | 1e-8 | PASS |
| q0_consistent, 2 random c (±3 %) | agree | rel 4.3e-10, 1.6e-9 | 1e-8 | PASS |
| gradient vs central FD of QCChi2 (4 points: 2 forms x 2 c, steps 1e-3 and 1e-4 relative) | <= 1e-4 | norm-rel 4.6e-6 to 1.8e-5; worst component 3.5e-5 | 1e-4 | PASS |
| penalty path, forced (q0, E_max=0.30, 8 of 12 levels lack a zero) | finite, counted, zero gradient rows | chi2 = 5.14e5 finite; `n_penalised_calls`=1; residual = 100.000 sigma for each of the 8; J rows 0; grad finite | - | PASS |
| D9/D8 (§2.6) | wrong-sign zeros flagged, stay in D3 matching; free levels excluded | code: `wrong_sign_roots()` = zeros with slope<0, kept in `energies()`; `QCChi2` records them (`matched`, `level_index`) and counts calls; unit tests `test_d9_*`, `test_real_trackers_d8_*` pass. Real trackers at C_BASELINE and C_LEGACY_0335, both forms, E_max 0.34: no wrong-sign zero occurs (list empty), so the real-data flag is exercised only by the fake-tracker unit test | - | PASS (see note 2) |
| tracker fix: regression tests present | exist and fail without the fix | I ran a copy of `qc3` in the scratchpad with the old `except (PoleError, DimensionJumpError)`: both `test_corrector_*` tests FAIL (111_A2 repro: `AmbiguousRootError: E0=0.33819589844776443 is not a zero (Delta_h = -1.731)`); with the fix both pass | - | PASS |
| `_check_loud` still raises for a genuinely ambiguous local-rescan candidate | raises | 2nd test injects an `ambiguous` candidate into the local scan: `pytest.raises(AmbiguousRootError)` passes; code path `_local_rescan -> _check_loud` read | - | PASS |
| cost per chi2 call, 1 core, q0_consistent, 40 calls each, random steps | builder: 86-112 ms (0.34), 187-241 ms (0.36) | E_max=0.34: 1e-3: 82/85/95 ms (min/med/max), 65 evals; 2e-2: 88/109/160 ms, 76 evals. E_max=0.36: 1e-3: 154/183/243 ms, 95 evals; 2e-2: 201/255/400 ms, 133 evals. 0 local, 0 full rescans, 0 penalised | informational | PASS (agrees with builder) |

## Notes
1. The chi2 differences between tracker and fresh scans (1e-9 relative at random c) come from the roots, which agree to about 1e-11; they are well below the 1e-8 bar.
2. D9 on real data: no wrong-sign zero appeared at the tested points; the reporting path is covered only by the synthetic test. Not a defect.
3. The FD-gradient mismatch (about 1.7e-5 for q0, 5e-6 for coded_v2) does not shrink with the step (1e-3 vs 1e-4), so it is a small systematic error in the Hellmann-Feynman $J$ (the $h=10^{-6}$ stencil), not FD noise. It is 6-20 times below the 1e-4 bar. Guess: the central-difference $Q'$ truncation. Worth remembering when Migrad/Hesse need the gradient to better than 1e-5.
4. PENDING: `test_cutoff_scan_q0_consistent_0335_two_stage_fit_verifies_and_records` and `..._fixed_KB_KE_two_parameter_fit_records_flags` (NotImplementedError, step 6).

New failures vs last report (016): none.

VERDICT: PASS (30 passed in the M2 acceptance run + 143 unit + 38 M1 non-slow; 0 failed; 2 pending, step 6).

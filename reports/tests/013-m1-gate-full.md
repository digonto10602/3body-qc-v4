# Tester report 013: M1 milestone gate, mode `full`

- **Date:** 2026-10-23
- **Git commit:** `c0334de` (HEAD). `qc3._core.git_commit()` == `c0334de` == `git rev-parse --short HEAD` after the final clean rebuild (see N1).
- **Prompt in force:** `prompts/002` (step 7 replan; steps 8, 9 and gate as in `prompts/001`). Acceptance tag `acceptance-v5` (47 tests).
- **Mode:** full. **Related:** `reports/audits/001-m1-milestone-audit.md`, `reports/tests/012-m1-step7-7-f-classes.md`; benchmark part: `reports/benchmarks/001-m1-single-energy-qc.md`.
- **Raw output:** `results/tester_013_m1_acceptance_pytest.txt`, `results/tester_013_unit_pytest.txt`, `results/tester_013_gate_full.json`.

**VERDICT: PASS** (47/47 acceptance, 79/79 unit, acceptance-check PASS, 0 pending, 0 failed).

## Commands run
```
graphify update .
rm src/python/qc3/_core*.so && ./build.sh core            # clean rebuild (98 s), twice, see N1
python3 -c "import qc3._core as c; print(c.git_commit())"  # PYTHONPATH=src/python
python3 -m pytest -q tests/acceptance/m1 -p no:cacheprovider -rA     # 4 min 32 s
python3 -m pytest -q tests --ignore=tests/acceptance -p no:cacheprovider -rA   # 71 s
.claude/bin/acceptance-check                                         # before and after the test runs
git diff b6ad563 -- src/qc3_core/oracle_point.hpp src/core src/fitter src/cachegen_cpu
python3 -m pytest -q <test id> -p no:cacheprovider   # 2x each, the two v5-corrected tests, run alone, no thread env vars set
```

## Results
| check | expected | got | tolerance | result |
|---|---|---|---|---|
| clean rebuild of `_core` | builds | built | n/a | PASS |
| `_core.git_commit()` vs `git rev-parse --short HEAD` | equal | `c0334de` == `c0334de` | exact | PASS (N1) |
| `tests/acceptance/m1`, all 47 incl. slow | 47 passed | 47 passed, 0 failed, 0 pending, 0 skipped (271.3 s) | exact | PASS |
| unit tests (`tests --ignore=tests/acceptance`) | all pass | 79 passed (71.1 s) | exact | PASS |
| `acceptance-check` | PASS | "PASS: tests/acceptance matches acceptance-v5" (before and after the runs) | n/a | PASS |
| `test_fast_single_thread_cpu_time` alone, unpinned, run 1 / 2 | pass | passed (1.62 s) / passed (1.60 s) | locked test | PASS |
| `test_oracle_reality_hermiticity_isometry` alone, unpinned, run 1 / 2 | pass | passed (7.85 s) / passed (8.64 s) | locked test | PASS |
| audit fix A: `test_fast_raises_on_nonfinite_K2inv_like_oracle` | pass | passed | n/a | PASS |
| audit fix B: `test_fast_raises_on_complex_prefactor_instead_of_dropping_imag[kw0,kw1]` | pass | both passed | n/a | PASS |
| audit fix N6: `test_xtol_guard_rejects_ambiguous_discriminant_settings` | pass | passed | n/a | PASS |
| audit fix N8: `test_pwave_form_is_always_passed_explicitly_to_the_core` | pass | passed | n/a | PASS |
| oracle sources unchanged since `b6ad563` (`oracle_point.hpp`, `src/core`, `src/fitter`, `src/cachegen_cpu`) | empty diff | empty diff (exit 0, 0 bytes) | exact | PASS |

I read the four audit-fix tests: each asserts a `RuntimeError`/`ValueError`/`KeyError` with a message match against an independent condition (for A, the oracle itself raises `BAD_F3_SCALE` on the same input), so none is tautological. The two timing tests passed with no `OMP_/OPENBLAS_/MKL_` variables set in the environment (checked with `env`).

## Notes
- **N1 (build/commit drift).** The first clean rebuild reported `6784003`; a WIP auto-commit (`c0334de`, only `graphify-out/manifest.json`, `git diff 6784003 c0334de --stat`: 1 file) moved HEAD during the build. I rebuilt again at `c0334de` and all tests above ran on that build. No `src/`, `tests/` or `results/` content differs between the two commits.
- Working tree clean after the runs (`git status --short` empty), so the tests wrote nothing into tracked files.
- `graphify affected refine_root` lists `test_refine_root_finds_simple_zero_and_raises_without_sign_change` and `test_fast_root_vs_oracle_100_A2_windows`; both ran and passed. `graphify affected evaluate_point` has no unique node; coverage taken from the full-suite run.
- Machine: 4-core Intel Xeon @ 2.1 GHz, 15 GiB, no GPU.

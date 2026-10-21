# Tester report 005: M1 step 5, portability + three no-silent-failure fixes

- **Date:** 2026-10-21
- **Git commit:** 5b396fc (`src/python/qc3/lattice.py` and `tests/unit/test_qc_lattice_roots.py` as committed there; I edited no source)
- **Prompt in force:** `prompts/001-m1-cache-free-single-energy-qc.md` (§3 row 5, §4 amendment; acceptance tag `acceptance-v3`)
- **Mode:** step
- **Related:** `reports/tests/004-m1-step4-qc-lattice.md` (notes 1-3 re-checked here)
- **Raw numbers:** `results/tester_005_m1_step5_portability.json` (commit, cwd, all outputs)
- **Script (scratchpad, not in repo):** `/tmp/claude-0/-home-user-3body-qc-v4/756827e6-de12-5b97-9946-026a972abb7d/scratchpad/chk_step5.py`; independent NumPy for the covariance, run from cwd `/` with absolute `PYTHONPATH`.
- **Scope exclusion (as briefed):** step 6 is being built concurrently. I did not run `test_m1_fast_vs_oracle.py` or the oracle tests, and did not rebuild `_core`.

## Verdict: PASS (2 non-blocking notes)

## Commands run
```
graphify update . ; graphify affected default_data_dir   # 3 tests, all run below
python3 -m pytest -q tests/acceptance/m1/test_m1_portability.py tests/acceptance/m1/test_m1_lattice.py tests/acceptance/m1/test_m1_qc_identities.py -p no:cacheprovider
python3 -m pytest -q tests/unit/test_qc_lattice_roots.py -p no:cacheprovider
cd / && PYTHONPATH=/home/user/3body-qc-v4/src/python python3 -W ignore chk_step5.py <scratch>/w results/tester_005_m1_step5_portability.json
grep -rIn -E "/home/digonto|/media/digonto|/home/digimonk" src/python build.sh src/qc3_core tests
.claude/bin/acceptance-check
```

## Results

| check | expected | got | tolerance | result |
|---|---|---|---|---|
| 1. acceptance portability + lattice + qc_identities | 14 green | 14 passed (0.15 s) | exact | PASS |
| 2. `tests/unit/test_qc_lattice_roots.py` | all green | 14 passed (0.49 s) | exact | PASS |
| 3a. `default_data_dir()`, `QC_DATA_ROOT` unset, cwd `/` | `<repo>/data`, exists | `/home/user/3body-qc-v4/data`, is_dir True | exact | PASS |
| 3b. `QC_DATA_ROOT=/some/where/else` | that path (no existence check) | `/some/where/else` | exact | PASS |
| 3c. `QC_DATA_ROOT=""` (empty) | falls back to `<repo>/data` | falls back | exact | PASS |
| 3d. `QC_DATA_ROOT=/nonexistent`, `load_levels` | loud, clear error | `FileNotFoundError: no lattice directory L20 under /nonexistent (tried lattice/L20, L20)` | | PASS |
| 4/fix 1. all 50 resamples with $E_{\rm lab}<P$ (synthetic 100_A2 n0, $E_{\rm lab}\approx0.05$, $P=0.0912$) | raise (was: level silently dropped) | `ValueError: <path>/20_001_A2_n0.jack: nonfinite converted Ecm in 50 of 50 jackknife resamples (E_lab^2 < P^2 = 0.00832096)` | | PASS |
| 4/fix 1. only some resamples nonfinite (mean near $P$, 55 of 60) | raise | same message, "55 of 60" | | PASS |
| 4/fix 2. $C=\mathrm{diag}(1,1,-0.5)$ (was: returned 2.0) | raise | `ValueError: ... not positive semi-definite (lambda_min/lambda_max = -0.5)` | | PASS |
| 4/fix 2. $C=0$ (was: returned 0.0) | raise | `ValueError: ... no positive eigenvalue (lambda_max <= 0)` | | PASS |
| 4/fix 2. $-\mathbb 1$; NaN in $C$; NaN in residual; asymmetric $C$; shape (3,) vs 4x4; non-square 3x2 | raise | all `ValueError` with specific messages (non-finite / not symmetric / shapes inconsistent) | | PASS |
| 4/fix 2. legit inputs not rejected: $\mathrm{diag}(1,1,0)$ (PSD singular, $r=(1,2,3)$); $\mathrm{diag}(1,1,-10^{-14})$ (round-off negative) | $\chi^2=1+4=5$ | 5.0 and 5.0 | exact | PASS |
| 4/fix 3. names `20_001_A2_n3x.jack`, `_n.jack`, `_n-1.jack`, `20_001_A2.jack`, `_n0_old.jack` next to a good `_n0` | raise | all 5 `ValueError: <path>: file name starts with '20_001_' but does not match '<Lbyas>_<ddd>_<irrep>_n<k>.jack' ...` | | PASS |
| 4/fix 3. unrelated files in the directory (`20_000_A1m_n0`, `24_001_A2_n0`, `readme.jack`) while loading 100_A2 | no error, 1 level | OK, 1 level | exact | PASS |
| 4/real. level counts at cutoff 0.315/0.325/0.335/0.345/0.355 | 6/9/12/14/17 | 6/9/12/14/17 | exact | PASS |
| 4/real. $\sigma$(20_000_A1m_n0) | $5.15296690388\times10^{-4}$ | $5.15296690388256\times10^{-4}$; rel. diff to report-004 value 0.0; to indep. sd/$\sqrt N$ $2.5\times10^{-14}$ | $10^{-10}$ | PASS |
| 4/real. mean $E^*$ of 20_000_A1m_n0 | 0.2692266019900497 | identical (diff 0.0) | $10^{-12}$ | PASS |
| 4/real. 12-level cov at 0.335 vs report 004 | same | $E^*$ diff 0.0, $\sigma_i$ rel. diff 0.0, $\lambda_{\min}=1.8231735849\times10^{-9}$, $\lambda_{\max}=3.2956162798\times10^{-6}$ (rel. diff 0.0), $\kappa=1.8076\times10^3$ | $10^{-12}$ | PASS |
| 4/real. 12-level cov vs my independent NumPy, labels | equal | labels equal; max rel. diff $1.2\times10^{-15}$ | $10^{-10}$ | PASS |
| 4/real. `chi2` on real cov, residual = 1e-3 * N(0,1), seed 7 | unchanged | 252.0499275024062; old formula (eigh pseudo-inverse, cut $10^{-12}$, no checks) gives the identical number (rel. diff 0.0); vs `solve`: $-1.0\times10^{-14}$ | $10^{-10}$ | PASS |
| 5. forbidden literals in `src/python`, `build.sh`, `src/qc3_core`, `tests/unit` | 0 | 0 files each | exact | PASS |
| 6. `acceptance-check` | PASS acceptance-v3 | "PASS: tests/acceptance matches acceptance-v3" | exact | PASS |

The only forbidden-literal hit under `tests/` is `tests/acceptance/m1/test_m1_portability.py:7`, the `FORBIDDEN` constant itself (the v3 amendment makes the test skip itself, as intended).

## Informational (not a failure): legacy files outside the M1 path still containing the literals
147 files in the repo (excluding `.git`, `graphify-out`); by directory:
`results/legacy_L20_passonly_hotwindow` 85, `configs/fitter` 12, `docs/original_handoffs` 7, `reports/research` 6, `data/windows` 5, `data/hot_windows_meta` 5, `results/legacy_L20_k3b_k3e_zero_tests` 3, `docs` 3, `configs/cachegen_gpu` 3, `src/fitter` 2, `src/cache_io` 2, `src/cache_build` 2, `scripts/windows` 2, `prompts` 2, `tests/acceptance` 1 (the FORBIDDEN constant), `tools/legacy_v32` 1, `scripts/plot` 1, `scripts/fit` 1, `scripts/classifier` 1, `legacy/root_finding_lineage` 1, `configs/legacy` 1, `configs/cache_build` 1.

## Notes (non-blocking)
1. `chi2` with an empty `(0,0)` covariance raises `ValueError: zero-size array to reduction operation maximum which has no identity` (from NumPy). It does raise, but the message is not informative. Guess: `np.abs(C).max()` is evaluated before a size check. Cosmetic; no real path produces an empty $C$ in `chi2` (zero levels are handled by `load_levels`).
2. `load_levels` raises on a nonfinite $E^*$ even for a level that is above the cutoff. This is conservative and matches the legacy C++ reader (`Nonfinite converted Ecm`); I record it as a design choice, not a defect.

## Failures
None. New failures vs. report 004: none. Pending: the step-6 tests (`test_m1_fast_vs_oracle.py` and the oracle-dependent tests) were deliberately not run, per the brief.

# Tester report 004: M1 step 4, QC helpers + lattice reader (D1, D3)

- **Date:** 2026-10-21
- **Git commit:** 9ad62c5 (`src/` and `tests/` are byte-identical to dc7147b, the commit named in the brief)
- **Prompt in force:** `prompts/001-m1-cache-free-single-energy-qc.md` (§3 row 4, §5.5–5.6; acceptance tag `acceptance-v2`)
- **Mode:** step
- **Related:** `reports/tests/003-m1-step3-reference-pwave.md`, `reports/research/004-lattice-data-parameters-and-inputs.md`, `reports/research/007-single-energy-qc-methods-and-reference.md`
- **Raw numbers:** `results/tester_004_m1_step4_qc_lattice.json` (commit, parameters, seeds, machine inside)
- **Scripts (not in repo, scratchpad):** `.../scratchpad/chk_lattice.py`, `.../scratchpad/chk_qc.py`; independent NumPy, they do not reuse builder code except as the object under test.

## Verdict: PASS (with 3 non-blocking notes)

## Commands run
```
python3 -m pytest -q tests/acceptance/m1/test_m1_qc_identities.py tests/acceptance/m1/test_m1_lattice.py -p no:cacheprovider
python3 -m pytest -q tests --ignore=tests/acceptance -p no:cacheprovider
python3 -m pytest -q tests/acceptance/m1 -m "not slow" -p no:cacheprovider -rA --tb=line
.claude/bin/acceptance-check
graphify update . ; graphify affected {project_Q,refine_root,load_levels,match_levels,chi2,slogdet}
python3 chk_lattice.py ; python3 chk_qc.py     # independent checks 4(a)-(g)
```

## Test suites

| check | expected | got | tolerance | result |
|---|---|---|---|---|
| step-4 ids (qc_identities 4 + lattice 8) | 12 green | 12 passed | exact | PASS |
| unit tests (`tests --ignore=tests/acceptance`) | all green | 16 passed | exact | PASS |
| M1 regression, not slow | 22 green, 5 NotImplementedError | 22 passed, 5 failed (all `NotImplementedError`) | exact | PASS |
| `acceptance-check` | PASS acceptance-v2 | "PASS: tests/acceptance matches acceptance-v2" | exact | PASS |

PENDING (not failures, `NotImplementedError`): `test_m1_fast_vs_oracle.py::{test_fast_vs_oracle_small_points, test_fast_basis_identical_to_oracle, test_fast_metadata_and_timings}` (step 6, `qc3.fast`); `test_m1_portability.py::{test_no_absolute_home_paths_in_m1_code, test_data_root_env_var}` (step 5, `default_data_dir`). 6 `slow` M1 tests deselected in step mode. New failures vs. report 003: none.

Graphify `affected` shows every changed symbol is covered by tests I ran (project_Q: 5 tests; load_levels: 6; match_levels: 2; refine_root: 1 unit test, plus `test_fast_root_vs_oracle_100_A2_windows`, which is pending on step 6).

## Independent physics / statistics checks

All "indep" numbers are my own NumPy (`np.loadtxt`, explicit formulas), not `qc3`.

| check | expected | got | tolerance | result |
|---|---|---|---|---|
| (a) mean $E$ of 20_000_A1m_n0 | $\bar x$ = 0.2692266019900497 | code = indep, identical to all digits | $10^{-12}$ rel | PASS |
| (a) $\sigma$ of 20_000_A1m_n0 | s.d./$\sqrt{603}$ = 5.15296690388243e-4 | 5.15296690388256e-4 (rel. diff $2.5\times10^{-14}$) | $10^{-10}$ rel | PASS |
| (a) $C_{00}$ vs $(N-1)/N\sum(y-\bar y)^2$ built by me | equal | rel. diff $2.2\times10^{-16}$ | $10^{-10}$ | PASS |
| (a) legacy factor (raw values treated as jackknife samples, `load_multil_jack_targets`: $(N-1)/N\sum(x-\bar x)^2$) | $(N-1)^2=602^2=362404$ in variance | variance ratio 362403.99999998; legacy $\sigma$ = 0.3102 vs corrected 5.153e-4 (ratio 602.0000) | $10^{-9}$ rel | PASS (confirms no inflation in the new code; quantifies the legacy bug) |
| (b) conversion applied per resample | $\sigma_{E^*}/\sigma_{E_{\rm lab}}\approx E_{\rm lab}/E^*$ (delta method) | 100_A2 n0: 1.04777 vs 1.04778; n1: 1.04518 vs 1.04519; 111_A2 n0: 1.11496 vs 1.11496; n1: 1.10610 vs 1.10610 (not 1.0, so $C$ is of $E^*$) | $10^{-4}$ rel | PASS |
| (b) code vs independent per-resample $E^*$ (8 moving-frame levels: 100_A2 n0–2, 111_A2 n0–4) | equal | mean rel. diff 0.0; $\sigma$ rel. diff $\le 4.4\times10^{-16}$ | $10^{-12}$ | PASS |
| (b) mean of resamples vs $\sqrt{\bar x^2-P^2}$ | differ slightly (nonlinear map) | differ by $\approx 6\times10^{-11}$ (100_A2 n0); code uses the resample mean, as the legacy reader does | info | PASS |
| (b) cutoff rule "keep if mean $E^*\le$ cutoff" | inclusive | cutoff = mean keeps the level (1); cutoff = nextafter(mean, 0) drops it (0) | exact | PASS |
| (c) counts at cutoff 0.315/0.325/0.335/0.345/0.355 | 6/9/12/14/17 (research 004; builder) | code 6/9/12/14/17; my recount from the 19 files 6/9/12/14/17 | exact | PASS |
| (c) 12-level cov at 0.335: symmetric | $C=C^T$ | max$|C-C^T|=0$ | $10^{-18}$ | PASS |
| (c) PSD | $\lambda_{\min}\ge0$ | $\lambda_{\min}=1.82\times10^{-9}$, $\lambda_{\max}=3.30\times10^{-6}$, condition number $\kappa=1.81\times10^{3}$ | $\lambda_{\min}\ge-10^{-12}\lambda_{\max}$ | PASS |
| (c) correlations | $0<\max\lvert\rho\rvert<1$ | min $\rho=0.491$, max $\rho=0.951$ (off-diagonal) | | PASS |
| (c) cov vs my independent $C$ (same level order) | equal | max rel. diff $1.2\times10^{-15}$; labels identical | $10^{-10}$ | PASS |
| (c) $\sigma$ of the 12 levels | $\sim5\times10^{-4}$ | 5.05e-4 to 7.70e-4; $N_{\rm cfg}=603$ | info | PASS |
| (d) D3 `match_levels`, 2000 random cases (seed 20261021, $m\in[1,7]$, $n\in[m,m+4]$, shuffled inputs) | lowest $m$ QC levels paired in sorted order with sorted lattice levels, indices into the original arrays, no duplicates | 0 mismatches | exact | PASS |
| (d) $n<m$ | `MissingRootError` | raised for (m,n) = (2,1), (5,4), (3,0): 3/3 | exact | PASS |
| (e) `chi2` on full-rank cov (3, 12; seed 7) | $r^TC^{-1}r$ | rel. diff 0.0 and $2.5\times10^{-15}$ ($\kappa$ = 22, 329); on the real 12x12 cov: $1.8\times10^{-14}$ | $10^{-10}$ | PASS |
| (e) relative cut semantics | drop eigenvalues $\le$ `eig_rel_cut`$\cdot\lambda_{\max}$ | with $\lambda=(1,10^{-11},10^{-13},0.5)$ and cut $10^{-9}$: matches the 2-direction sum to $1.7\times10^{-15}$; default $10^{-12}$ drops only $10^{-13}$ (match $3\times10^{-5}$, limited by `eigh` accuracy $\epsilon\lambda_{\max}/10^{-11}$ of my test matrix, not a code error); scale invariant ($C\to10^{-6}C$) | | PASS |
| (f) `slogdet` singular / non-finite | raise, never $-\infty$ | zeros, rank-1 (ones), NaN, Inf, empty, non-square, non-real phase: all `ValueError` | exact | PASS |
| (f) `slogdet` sign type and value | `int` $\pm1$; equals $\prod\mathrm{sign}(\lambda_i)$ | `int`; 30x30 symmetric (14 negative eigenvalues): sign 1 = eigenvalue sign; $\ln\lvert\det\rvert$ differs from $\sum\ln\lvert\lambda\rvert$ by $<10^{-13}$ | $10^{-10}$ | PASS |
| (f) `project_Q` vs `oracle.py` assembly, 100_A2 E=0.291556099879, baseline $\vec c$, $N=56$, $N_{\rm proj}=17$ | arrays identical | `np.array_equal(project_Q(F3inv,Vsel,K3,L), r.Q)` = True (complex128); `slogdet` = $(-1,\,-118.96068499583154)$ identical to the oracle (exact `==`); Hermiticity $5\times10^{-19}$ | exact | PASS |
| (g) `refine_root`, oracle, 100_A2 level 0 bracket $[0.29155126256,\,0.29155610781]$ (smallest $N$=56), baseline $\vec c$ | converged within xtol = $10^{-11}$ | root = 0.29155609960818824; 5 oracle calls (2 end points + 3 Brent steps), 13.8 s; $\det Q$ sign flips between root$\mp2\times10^{-11}$ ($\ln\lvert\det\rvert$ -121.57 at both sides vs -109 and -115 at the bracket ends); at the root $\ln\lvert\det\rvert$ = -131.7 | $10^{-11}$ | PASS |
| (g) informational | | root differs from the window file's `previous_model_zero` (0.29155609987916675) by $-2.7\times10^{-10}$ (that value comes from the older grid pipeline; the oracle root is the finer one) | info | n/a |
| 5. hard-coded expected numbers in `qc.py`, `roots.py`, `lattice.py` | none | grep for decimals with $\ge3$ digits: none (the only constants are tolerances/defaults `1e-12`, `1e-10`, `3.444` default `xi`) | | PASS |
| 5. silent exception catching in numerical code | none | one `try/except ValueError: continue` in `load_levels`, in file-name parsing only (see note 3). No other `except`. | | PASS |

## Notes (non-blocking; for the builder or auditor, no fixes done by me)

1. **NaN energies are dropped silently in `load_levels`.** If some resample has $E_{\rm lab}^2<P^2$ then $E^*=\mathrm{NaN}$; the mean is NaN and `nan <= cutoff` is False, so the level disappears with only NumPy's `invalid value encountered in sqrt` warning. Reproduced with a synthetic 100_A2 file at $E_{\rm lab}\approx0.05<P=0.0912$: the level was gone and no exception was raised. The legacy C++ reader throws `Nonfinite converted Ecm`. Not triggered by the real L20 files. Guess: missing `np.isfinite` check after `elab_to_ecm`.
2. **`chi2` drops strongly negative eigenvalues and returns 0.0 for an all-zero covariance.** `keep = w > eig_rel_cut*w.max()`; for `C=diag(1,1,-0.5)` it returns 2.0 (negative direction ignored), for `C=0` it returns 0.0. A covariance with a large negative eigenvalue or zero rank is a data error that should raise. Not triggered by the real data ($\lambda_{\min}=1.8\times10^{-9}>0$).
3. **File discovery skips unparsable names silently** (`except ValueError: continue` on the `_n<k>` suffix), e.g. `20_001_A2_n3x.jack`. Harmless for the 19 L20 files (test asserts 19), but a typo'd name would silently reduce the level count.
4. Minor test quality (unit test file): `test_jackknife_samples_input_gives_sqrt_N_minus_1_over_N_sum` compares the formula with itself (tautology). The acceptance tests cover the real behaviour, so no gap.
5. `refine_root` unit test uses a fake linear $\det$; the oracle check above (g) covers the real path. `test_fast_root_vs_oracle_100_A2_windows` stays PENDING on step 6.

## Failures
None.

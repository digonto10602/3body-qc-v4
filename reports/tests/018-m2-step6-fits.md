# Tester report 018: M2 step 6 (fit_cutoff / verify_fit / cutoff_scan / save_results / CLI / code equivalence)

- Date: 2026-10-25 (runs started 2026-10-24). Commit `95ebb88` (tree clean except `graphify-out/`). Mode: step.
- Plan in force: `prompts/003-m2-label-free-roots-tracking-fit.md` §2.7, §3 step 6, §4, §5, §7. Acceptance tag `acceptance-v7`.
- Related: `reports/tests/017-m2-step5-qcchi2.md`. Raw numbers: `results/tester_018_step6_independent_checks.json`; scripts: `results/tester_018_scripts/`.

## VERDICT: PASS (with notes N1 to N4)

## Commands
```
python3 -m pytest -q tests --ignore=tests/acceptance                 # 153 passed, 75 s
python3 -m pytest -q tests/acceptance/m1 -m "not slow"               # 38 passed, 9 deselected, 35 s
python3 -m pytest -q tests/acceptance/m2 -m "not slow"               # 20 passed, 12 deselected, 57 s
python3 -m pytest -q tests/acceptance/m2/test_m2_fit.py              # 5 passed (incl. the 3 slow step-6 tests), 344 s
.claude/bin/acceptance-check                                         # PASS: matches acceptance-v7
python3 results/tester_018_scripts/{run_min,hscan,hscan2,probe,ce,audit_cmp}.py   # independent checks below
```

## Test suites
| check | expected | got | PASS/FAIL |
|---|---|---|---|
| unit tests | all | 153/153 | PASS |
| m1 non-slow | all | 38/38 | PASS |
| m2 non-slow | all | 20/20 | PASS |
| `test_m2_fit.py` (slow) | 5/5 | 5/5 | PASS |
| acceptance-check | PASS | PASS (v7) | PASS |
| pending / new failures vs 017 | 0 | 0 / 0 | PASS |

## Independent chi^2 (own fresh scans, own NumPy D1 covariance from the `.jack` files, own D3; no tracker)
Covariance/E_lat agree with `qc3.lattice` to $1.2\times10^{-15}$ (relative). Tolerance $10^{-8}$ relative.

| point | reported $\chi^2$ | mine | rel. diff | PASS |
|---|---|---|---|---|
| A (4 free) 0.335 | 29.614722295 | 29.614722337 | $1.4\times10^{-9}$ | PASS |
| A 0.345 | 35.542154307 | 35.542154228 | $2.2\times10^{-9}$ | PASS |
| B ($K_B=K_E=0$) 0.335 | 36.586541445 | 36.586541398 | $1.3\times10^{-9}$ | PASS |
| B 0.345 (extra) | 55.223275598 | 55.223275725 | $2.3\times10^{-9}$ | PASS |
| code equiv. $\chi^2(\vec c_{\rm v32w})$ | 55.0420345549 | 55.0420345549 | $<10^{-12}$ | PASS |
| code equiv. $\chi^2(\vec c_{\rm legacy})$ | 50.0102930423 | 50.0102930423 | $<10^{-12}$ | PASS |
| code equiv. refit | 31.6326867062 | 31.6326866466 | $2\times10^{-9}$ | PASS |

## Is it a minimum? Gradient and $\Delta\chi^2$
Gradient: central differences of my $\chi^2$ with steps $0.003, 0.01, 0.03\,\sigma$. In every case $g_i\sigma_i$ scales as step$^2$ (A/0.335 comp. 0: $-4.0\times10^{-4}, -4.6\times10^{-3}, -4.2\times10^{-2}$ at 0.003/0.01/0.03), i.e. it is the cubic truncation term and the true gradient is $\lesssim10^{-4}$ in units of $\chi^2/\sigma$. Largest at the smallest step: $|g\sigma|=1.6\times10^{-3}$ (A/0.345 comp. 3, step $0.003\sigma_{\rm FD}$, again $\propto$ step$^2$). Minima are real. PASS.
A caution: the cubic term is large (A/0.335: $f'''\sigma^3\approx2\,f''\sigma^2$), so $\chi^2$ is not quadratic over 1$\sigma$.

$\Delta\chi^2$ test. Step $\delta=\pm s\,\mathrm{cov}[:,i]/\sigma_i$ (the 1-sigma conditional shift of the reported Hesse covariance; exactly $\Delta\chi^2=1$ at $s=1$ if the covariance is right). Column "implied $\sigma/\sigma_H$" $=1/\sqrt{\Delta\chi^2(0.1\sigma_H)/0.01}$ is the true error scale relative to Hesse.

| point, param | $\Delta\chi^2(\pm1\sigma_H)$ | mean | implied $\sigma/\sigma_H$ | FD-Hessian $\sigma_{\rm FD}/\sigma_H$ |
|---|---|---|---|---|
| A 0.335, $K_0$ | 0.619 / 0.862 | 0.74 | 1.187 | 1.191 |
| A 0.335, $K_1$ | 0.724 / 0.737 | 0.73 | 1.171 | 1.174 |
| A 0.335, $K_B$ | 0.725 / 0.846 | 0.79 | 1.135 | 1.136 |
| A 0.335, $K_E$ | 0.799 / 0.772 | 0.79 | 1.133 | 1.133 |
| A 0.345, $K_0$ | 0.370 / 0.398 | 0.38 | 1.64 | 1.89 |
| A 0.345, $K_1$ | 3.70 / 2.59 | 3.15 | 0.61 | 0.62 |
| A 0.345, $K_B$ | 2.48 / 3.83 | 3.15 | 0.62 | 0.67 |
| A 0.345, $K_E$ | 3.69 / 2.56 | 3.13 | 0.62 | 0.67 |
| B 0.335, $K_0$ | 0.883 / 1.193 | 1.04 | 0.991 | 0.992 |
| B 0.335, $K_1$ | 1.209 / 0.869 | 1.04 | 0.992 | 0.992 |
| B 0.345, $K_0$ | 8.36 / 23.3 | 15.9 | 0.279 | 0.283 |
| B 0.345, $K_1$ | 23.4 / 8.36 | 15.9 | 0.279 | 0.280 |

Only B/0.335 has reliable Hesse errors ($\Delta\chi^2\approx1.04$). Hesse errors are wrong at A/0.335 (too small by 14 to 19 %, `has_accurate_covar` False, flagged), at A/0.345 (K_0 too small ×1.6 to 1.9; the others too large ×1.5 to 1.6; the flags `error_exceeds_value`/`correlation_ge_0.99` are set there), and at B/0.345 (too large by 3.5, `hesse_not_accurate` is flagged). So the Hesse errors are not safe to quote at A/0.335, A/0.345, B/0.345.

## Builder deviations
**(a) LM-only `reject_ambiguous_trials`: safe for the reported minimum, with one caveat.** Code check: `rho()` returns the all-penalised residual only inside LM; `jac`, Migrad, the start, the final model and `verify_fit` still raise. The reported minimum is produced by Migrad from `c_lm` and then checked (my $\chi^2$ above, $|g\sigma|\lesssim10^{-3}$, fresh scan counts and roots), so rejection cannot make a non-minimum look like a minimum. Events: B 0.315 (1), 0.325 (2), 0.345 (1); all `ambiguous` candidates at 100_A2, $E\approx0.29130$, i.e. a zero within about $10^{-6}$ of a pole at 0.29130. The fit minima sit $6\times10^{-5}$ to $2\times10^{-4}$ below that pole, so the wall is not binding at the minima. Caveat (N2): the rejection acts as a wall and can steer LM to the side where the start is, i.e. it chooses a local minimum. At B/0.335 a probe along $K_1$ to $-3\sigma_H$ found no lower $\chi^2$ (monotone rise, +2.7 at $-0.3\sigma_H$ to +311 at $-2\sigma_H$). PASS, recorded not silent (`ambiguous_trials_lm`).

**(b) Hesse vs FD-Hessian: the FD-Hessian (finite difference of the analytic gradient) is the trustworthy one.** In the three cases with a real discrepancy the $\Delta\chi^2$-implied error scale equals $\sigma_{\rm FD}/\sigma_H$ to 0.3 % (A/0.335, B/0.345) or 5 to 10 % (A/0.345, where $\chi^2$ is strongly non-quadratic). My own second-difference Hessian of $\chi^2$ is ill-conditioned (H eigenvalues span $10^4$) and over-estimates by 20 to 50 %, so do not use it as a referee; the profile $\Delta\chi^2$ is the referee. Correction to the brief: the level 7e-6 above the dimension jump 0.3304785 is in **A/0.345** (200_A2 $n=2$ at 0.3304855), not B/0.345 (that level is at 0.33069). B/0.345's 3.5x is a plain Hesse failure.

**(c) `chi2_covariance_convention_check.py`: the only difference is converting each configuration instead of each jackknife resample, but the gap is larger than 0.75 for Table B.** I reproduced our $\chi^2$ with the D1 recipe to $10^{-8}$ (all 8 rows) and ran the tool's recipe on the same model energies. Tool minus ours: A $+0.352, -0.755, -0.569, -0.011$ (0.315 to 0.345); B $-1.588, -2.287, -1.633, -3.182$. Hybrids show both pieces matter: tool covariance with D1 means moves $\chi^2$ by $0.1$ to $0.2$ at A/0.335 (29.436 vs 29.615), the D1 covariance with the tool means moves it by $-0.4$ (29.207). The covariance differs by 1 % in the max element, the means by up to $1.2\times10^{-4}$ (0.19 $\sigma$) because $\sqrt{\bar x^2-P^2}$ is not linear in the noisy per-configuration value. Judged: yes the only difference (substituting resamples for configurations in the tool gives ours exactly), but the "up to 0.75" holds for Table A only; Table B goes to 3.2, so a plan line "same $\chi^2$ to $10^{-6}$" cannot hold with the tool's recipe as written.

## Fresh-scan verification (all 9 fits in the 3 JSONs)
Re-derived from `per_sector` tracked vs fresh arrays up to (highest matched root + 0.002): counts equal in all 9, `max_root_diff` $\le2.0\times10^{-12}$ (bar $10^{-9}$), $|\chi^2_{\rm fresh}-\chi^2|\le1.4\times10^{-7}$ (bound $1.7\times10^{-4}$ to $1.4\times10^{-3}$), matched roots equal the lowest $m$ tracked, `n_penalised_calls_final_stage` 0, `ok` True. PASS.

## Provenance (JSONs)
Present: commit `e9d0215`, `git_tree_dirty: true`, core built at `a93e308` (0.5.1-step7.7), numpy 2.5.3, scipy 1.18.1, iminuit 2.33.0, machine, command, threads env, production constants, form, c_start, fixed, D1/D3/D8/D9 strings, grid, tracking, seeds (deterministic, none). Notes: **(N1)** the recorded commit is the parent of the WIP commit, tree dirty, so the exact code is `95ebb88`'s `fit.py` (it differs from `e9d0215`); results should be regenerated or the dirty state noted at the milestone gate. Code-equivalence file: $\vec c_{\rm v32w}$ $\chi^2$ 55.0420, legacy $\vec c$ 50.0103, level diffs $\le6.1\times10^{-7}$ (all) and $4.0\times10^{-8}$ (no pole within 10 rows), as the plan states. PASS.

## Findings for the coordinator and Digonto
- **N2/N3 (physics, not a code bug):** at the A/0.335 minimum the 110_A2 $n=1$ model level (0.3154136) is $2.4\times10^{-5}$ (5 rows) above a pole (0.3153899). A zero crosses this pole at $+0.115\sigma_H$ in $K_0$ (level jumps from 0.315411 to 0.315388, $\chi^2$ jump about 0.2). Similar zero-pole or dimension-jump walls are within 1 $\sigma_H$ in A/0.345 (200_A2 level 7e-6 above the jump 0.3304785; ambiguous 111_A2 at 0.33765), B/0.335 and B/0.345 (100_A2 pole 0.29130). So $\chi^2(\vec c)$ is piecewise-smooth and the quoted Hesse $\pm1\sigma$ intervals cross walls. Recommend quoting the FD-Hessian errors beside Hesse for A/0.335, A/0.345, B/0.345, and mention the walls to Digonto.
- **N4:** `fd_hessian_rel` and `fd_hessian_note` are recorded (A/0.345 needed 0.002 after 0.02 hit a dimension jump); good, not silent.

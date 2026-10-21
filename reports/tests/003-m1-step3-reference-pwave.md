# Tester report 003: M1 step 3, reference checks and the `q0_consistent` p-wave form

- **Date:** 2026-10-21
- **Git commit:** `3763718` (HEAD; the step-3 sources are identical to those of `c4fcb81`; the working tree also holds the step-4 builder's uncommitted edits in `src/python/qc3/{lattice,qc,roots}.py` and `tests/unit/test_qc_lattice_roots.py`, none of which touch the step-3 code). `_core.git_commit()` = `3763718`.
- **Prompt in force:** `prompts/001-m1-cache-free-single-energy-qc.md` (§3 row 3, §2.2, §4)
- **Mode:** step
- **Related:** `reports/tests/002-m1-step2-oracle-wrapper.md`, `reports/research/007-single-energy-qc-methods-and-reference.md`
- **Raw numbers:** `results/tester_003_m1_step3_reference_pwave.json` (commit, parameters, machine; no RNG, so no seed)
- **Machine:** 4 vCPU Linux 6.18, g++ 13.3.0 `-O3 -fopenmp`, Python 3.13.16
- **Scope note (coordinator instruction):** a builder is working on M1 step 4 in parallel. `test_m1_qc_identities.py`, `test_m1_lattice.py` and the unit file `tests/unit/test_qc_lattice_roots.py` are step-4 territory and are **excluded from the step-3 verdict**. They are listed below for information only.

## Verdict: PASS

## Commands run
```
graphify update . ; graphify affected "K2inv_EREord2_2plus1_mat" ; graphify affected "Sector"
git diff --stat {f7613c6,8b3e7d8} HEAD -- src/core src/fitter src/cachegen_cpu src/cache_build
rm -f src/python/qc3/_core*.so; ./build.sh core                       # clean rebuild, 85.1 s
python3 -m pytest -q tests/acceptance/m1/test_m1_reference.py tests/acceptance/m1/test_m1_oracle.py -p no:cacheprovider --durations=0
python3 -m pytest -q tests --ignore=tests/acceptance -p no:cacheprovider -v
python3 -m pytest -q tests/acceptance/m1 -m "not slow" -p no:cacheprovider -rA --tb=line
.claude/bin/acceptance-check
python3 <scratchpad>/check_s3.py     # independent NumPy K2inv + oracle-vs-qc3ref sweep, both forms (not added to tests/)
```
`graphify affected "K2inv_EREord2_2plus1_mat"` lists the legacy and fitter callers (`build_component_entry`, `build_cache_entry`, ...) plus the new `K2inv_2plus1_by_form`. All the former are byte-unchanged (check 4a), so they need no new test. `graphify affected "Sector"` lists `oracle.py`, `fast.py`, `point.py`, `m1_helpers.sector` and the unit tests; all are run below.

## Results

| # | check | expected | got | tolerance | result |
|---|---|---|---|---|---|
| 1 | clean rebuild `./build.sh core` | compiles | OK, 85.1 s; `_core.git_commit()` = `3763718` | n/a | PASS |
| 2 | `test_m1_reference.py` (4 fast + 1 slow) | 5 green | **5 passed** (slow large-point test 19.9 s) | n/a | PASS |
| 2b | `test_m1_oracle.py` regression | 7 green | **7 passed** (12 passed together in 81.9 s) | n/a | PASS |
| 3 | unit tests, step-3 file `tests/unit/test_sector_oracle.py` | 6 pass (5 old + new `test_k2_pwave_forms_differ_only_in_channel1_pwave_diagonal`) | **6 passed** | n/a | PASS |
| 3b | unit tests, whole `tests/unit` (info) | | 16 passed (10 in step-4 file `test_qc_lattice_roots.py`, ignored for the verdict) | n/a | info |
| 4a | production code byte-unchanged | empty diff | `git diff --stat f7613c6 HEAD` and `8b3e7d8 HEAD` on `src/core src/fitter src/cachegen_cpu src/cache_build`: **empty** both times. Step-3 diff touches only `oracle.py`, `module.cpp`, `oracle_point.hpp`, new `k2inv_forms.hpp`, and the unit test. | exact | PASS |
| 4b | `q0_consistent` = §2.2 element, term by term | see below | matches | n/a | PASS |
| 4c | $\vec c=0$ anchors, `coded_v2`, $\Delta\ln\lvert\det Q\rvert$ | 0 | **0.0 exactly** at all 5 points tested here (100_A2 0.2915561, 100_A2 0.2917984, 000_A1m 0.26932, 111_A2 0.32036, 200_A2 0.28664); signs equal the anchor signs. All 12 anchors (incl. large) pass in `test_oracle_anchors_all_points`. | $10^{-10}$ | PASS |
| 4d | independent NumPy element from §2.2 vs oracle | $\lesssim10^{-14}$ | full-matrix rel. Frobenius $\le 2.6\times10^{-15}$ (`coded_v2`), $\le 2.3\times10^{-16}$ (`q0_consistent`); p-wave elementwise max relative error $\le 4.9\times10^{-15}$ / $\le 4.0\times10^{-16}$ | $10^{-12}$ | PASS |
| 5 | oracle vs `qc3ref` (table below) | $\le10^{-12}$ | all $\le 3.5\times10^{-15}$ | $10^{-12}$ | PASS |
| 6 | fast tally (see below) | 10 green (6 oracle + 4 reference) | **10 green** for the step-3 files | n/a | PASS |
| 7 | `acceptance-check` | PASS | "PASS: tests/acceptance matches acceptance-v1" | n/a | PASS |

Counts for step 3: step-3 ids 5/5 pass; oracle regression 7/7; step-3 unit 6/6; new failures vs report 002: 0.

### 4b. Term-by-term check of `K2_inv_ERE_ang_mom_q0_consistent` against prompts/001 §2.2
Target: $\frac{1}{2\omega_kL^3}\frac{\eta_1}{8\pi\sqrt\sigma}\big[-\frac1{a_1}+q^{*2}|q^*|(1-H)\big]$ (channel 2 multiplied by 2).

| term | in `k2inv_forms.hpp` | OK |
|---|---|---|
| $\eta/(8\pi\sqrt\sigma)$ | `A = eta_i/(8π sqrt(sigma_i))` | yes (same as production) |
| $-1/a_1 + \tfrac12 r_1q^2+P_1q^4$ | `B = -1/a1 + 0.5*r1*q2psq + P1*q2psq^2`; $r_1=P_1=0$ in production, so $B=-1/a_1$ | yes |
| $q^{*2}\lvert q^*\rvert(1-H)$ | `q2psq*C*D`, `C = abs(sqrt(q2psq_star))`, `D = 1 - cutoff_function_1(...)` (same helpers as production, $\epsilon_H=0$) | yes |
| $1/(2\omega_kL^3)$ | `1.0/(2 ω_k L L L)` with $\omega_k$ = spectator energy, in `..._i_mat_q0_consistent` | yes |
| channel masses | ch.1: `(m1,m1,m2)`; ch.2: `(m2,m1,m1)`, same as production | yes |
| channel-2 factor | `bottomRightCorner = 2.0*K2` | yes |
| other $\ell$ / non-diagonal | delegates to unchanged `K2_inv_ERE_ang_mom`; matrix is diagonal | yes |
| dispatcher | `"coded_v2"` calls the unchanged `K2inv_EREord2_2plus1_mat`; unknown string throws; Python also raises `ValueError` first | yes |

The coded $\ell=1$ branch, for comparison, is `-1/(a1 q2) + 0.5 r1 + P1 q2` then `A*(B + q2*C*D)`. So only the bracket $B$ differs; the threshold term is shared. Confirmed numerically: $\mathcal K_2^{-1}$ differs between forms **only** on channel-1 $\ell=1$ diagonal entries (zero difference everywhere else; $F$ and $G$ bit-identical between forms; `form_diff_rel_K2inv` $\approx1.00$ because the coded form is $\sim1/q^{*2}\gg1$ times larger).

Exercise of the threshold term: in the p-wave channel-1 rows, the number with $q^{*2}<0$ and $0<H<1$ (so that $(1-H)|q|$ matters) is 24/30 (100_A2 0.2915561), 18/21 (000_A1m 0.26932), 18/21 (200_A2 0.28664), 30/60 (110_A2 0.35), 72/96 (111_A2 0.35103). The $q^{*2}|q^*|(1-H)$ term and $H$ are therefore tested, not just the $-1/a_1$ term.

### 4d. Independent implementation
My own NumPy (scratchpad `check_s3.py`; not `qc3ref`, not C++) builds the diagonal from §2.2 directly: smooth step $J(z)=\exp[-e^{-1/(1-z)}/z]$ with $z=(\sigma-|m_j^2-m_k^2|)/((m_j+m_k)^2-|m_j^2-m_k^2|)$, $q^{*2}=(\sigma-(m_j+m_k)^2)(\sigma-(m_j-m_k)^2)/(4\sigma)$, $\sigma=(E-\omega_k)^2-(\mathbf P-\mathbf k)^2$, $\eta=(1,\tfrac12)$, channel-2 ×2. Agrees with the oracle at the levels in row 4d for both forms at all 6 points.

### 5. Max relative Frobenius deviation oracle vs `qc3ref` ($\lVert a-b\rVert_F/\lVert b\rVert_F$, reference in the denominator, $\vec c=0$)

Small set = the four anchor points with $N\le60$: 100_A2 0.2915561, 100_A2 0.2917984, 000_A1m 0.26932, 200_A2 0.28664 (note: 111_A2 0.32036 has $N=64$ and is outside the helper's "small" cut). Large = 110_A2 0.35000 ($N=108$) and 111_A2 0.35103 ($N=160$).

| quantity | small points (max of 4) | large points (max of 2) | tolerance |
|---|---|---|---|
| $F$ | $9.7\times10^{-16}$ | $7.7\times10^{-16}$ | $10^{-12}$ |
| $G$ | $2.4\times10^{-16}$ | $2.3\times10^{-16}$ | $10^{-12}$ |
| $\mathcal K_2^{-1}$, `coded_v2` | $1.2\times10^{-16}$ | $2.0\times10^{-16}$ | $10^{-12}$ |
| $\mathcal K_2^{-1}$, `q0_consistent` | $1.1\times10^{-16}$ | $1.1\times10^{-16}$ | $10^{-12}$ |
| $K_3$ iso0 | $1.7\times10^{-15}$ | $3.5\times10^{-15}$ | $10^{-12}$ |
| $K_3$ iso1 | $2.3\times10^{-15}$ | $3.5\times10^{-15}$ | $10^{-12}$ |
| $K_3$ B | $1.2\times10^{-15}$ | $2.5\times10^{-15}$ | $10^{-12}$ |
| $K_3$ E | $3.7\times10^{-16}$ | $4.8\times10^{-16}$ | $10^{-12}$ |

Margin to tolerance is $\ge 280$. Note that the locked acceptance test checks `q0_consistent` at one small point only (100_A2 0.2915561) and the large points only for `coded_v2`; the table above extends both to all six points, and `q0_consistent` passes everywhere. Per-point numbers are in the JSON.

### 6. Fast tally (`pytest tests/acceptance/m1 -m "not slow" -rA --tb=line`)
27 selected: **21 passed, 6 failed, 6 deselected (slow)**.
- Step-3-relevant, green: `test_m1_oracle.py` 6/6 fast + `test_m1_reference.py` 4/4 fast = **10 green, as expected**.
- Failing, all PENDING, not step-3: `test_m1_fast_vs_oracle.py` x3 (`NotImplementedError: qc3.fast.evaluate_point: step 6`), `test_m1_portability.py` x2 (`NotImplementedError: qc3.lattice.default_data_dir: step 5`), `test_m1_qc_identities.py::test_scale_applies_to_matrix_not_determinant` (step-4 file, in flux).
- Passing but step-4 territory and ignored for this verdict: `test_m1_lattice.py` 8/8 and 3 of the 4 `test_m1_qc_identities.py` fast tests (the builder's in-progress work).

## Failures
None for step 3.

## Notes (not failures)
1. The locked tests check `q0_consistent` only at one small point; I extended the check to six points (above). If the coordinator wants that permanently in the repo, it is a one-line extension of `tests/unit` (I did not add it).
2. `metadata["git_commit"]` is `3763718`, the HEAD at build time. Rebuild `_core` after the next commit before producing any result that leaves the project (as in report 002).
3. Which form matches HadSpec's $a_1$ remains the pending Digonto decision (D7 in `notes/STATUS.md`); this step only supplies the switch, default `coded_v2` unchanged.

# Tester report 020: M2 gate, mode full (rerun after audit-002 fixes)

- Date: 2026-10-27. Commit `9a37450` (HEAD; only `graphify-out/` differs from `c95a5c4`). `qc3._core.git_commit()` = `9a37450` after `rm _core*.so; ./build.sh core` (107 s). Acceptance tag `acceptance-v8`. Plan in force `prompts/003`.
- Related: `reports/tests/019-m2-gate-full.md` (previous gate, FAIL on one anchor tolerance), `reports/audits/002-m2-milestone-audit.md`, benchmark `reports/benchmarks/002-m2-fit.md`.
- Machine: Xeon @ 2.10 GHz, 4 cores, OpenBLAS, `OMP_NUM_THREADS=OPENBLAS_NUM_THREADS=1`; acceptance pinned `taskset -c 0`.

## VERDICT: PASS

| check | expected | got | tolerance | result |
|---|---|---|---|---|
| clean `_core` rebuild | `git_commit()` == HEAD | `9a37450` == `9a37450` | exact | PASS |
| `pytest -q tests/acceptance` (incl. slow) | 79/79 | 79 passed, 0 failed, 0 pending, 794 s | | PASS |
| `pytest -q tests --ignore=tests/acceptance` | all green | 172 passed (84 s) | | PASS |
| `.claude/bin/acceptance-check` | PASS | PASS (matches acceptance-v8) | | PASS |
| `scripts/m2/score_labels.py` | exit 0, FP=FN=0 | exit 0; TP/FP/FN/TN = 31/0/0/176 (7 sectors), 53191 evals, 96 s | FP=FN=0 | PASS |
| `git diff 191530d -- src/qc3_core/oracle_point.hpp src/core src/fitter src/cachegen_cpu` | empty | empty | | PASS |
| graphify `affected fit_cutoff` | covered | cutoff_scan, code_equivalence, main: covered by test_m2_step6_fit, test_audit002_notes, m2 acceptance fits | | PASS |

- The 019 failure (`test_oracle_anchors_all_points`, 200_A2) now passes under the kappa-scaled tolerance of acceptance-v8.
- One warning in the acceptance run (not a failure): `test_cutoff_scan_q0_consistent_0335_fixed_KB_KE_...` logs a full rescan in 111_A2 that changed the tracked set 3 -> 2 roots (the third, 0.33821, left the tracked range); this is the audit-N3 loud-event path working.
- `score_labels.py` was run with `--out results/tester_m2gate_label_gate.json` (not the tracked `results/m2_label_gate.json`, so nothing to revert; the tracked file is untouched).
- Raw: `results/tester_m2gate_acceptance.log`, `results/tester_m2gate_label_gate.json`.

## Extra physics check: is the 110_A2 zero at the free level E_f a genuine artefact?

$E_f = 0.31946882819441$ (multiplicity 2 non-interacting level, `free_levels`). Scripts `results/tester_m2gate_scripts/freecheck{,2,3}.py`; raw `results/tester_m2gate_freelevel_*.json`; seed 20261027 (numpy `default_rng`; two draws of $\vec c=\vec c_{\rm base}\,(1+U(-0.1,0.1))$ per component), 5 vectors $\vec c$ x 2 p-wave forms (`coded_v2`, `q0_consistent`). `scan_sector` on [0.317, 0.322], then dense probes around $E_f$.

First, an important correction to the premise. The "simple sign-change zero, order 0.99996, slope +0.319, $3.6\times10^{-13}$ from $E_f$" is the **synthetic** `Q` in `tests/unit/test_tracking.py::test_free_level_solution_is_excluded_from_roots_and_reported` (the fake `F` has a linear factor $(E-E_f)$ by construction; the warning text with slope 0.31882819 comes from that test). It is not a property of the real core. I did not find such a zero in the real 110_A2 at any of the 10 (form, $\vec c$) combinations.

What the real core shows at $E_f$ (all 10 combinations): one candidate, kind `even_zero`, estimated order $3.9999$ to $4.0001$ (relabelled `free_level`, no sign change), no sign-change zero within $2\times10^{-3}$ of $E_f$ other than ordinary levels.

| form | $\vec c$ | order at $E_f$ | exponent of $|\det Q|$ vs $|\delta|$, $\delta=\pm(10^{-7}..10^{-6})$ | zero position from quartic fit $E_0-E_f$ (+side / -side) | nearest ordinary zero $E-E_f$ |
|---|---|---|---|---|---|
| q0_consistent | baseline | 3.99999 | 3.999 / 4.001 | $+2.3\times10^{-12}$ / $-1.8\times10^{-11}$ | $-5.3\times10^{-4}$ |
| q0_consistent | A/0.335 min | 4.00000 | 3.999 / 4.001 | $-9.7\times10^{-11}$ / $+9.7\times10^{-11}$ | $-8.7\times10^{-4}$ |
| q0_consistent | rand1 | 3.99999 | 3.999 / 4.001 | $+2.8\times10^{-12}$ / $-1.8\times10^{-11}$ | $-5.3\times10^{-4}$ |
| q0_consistent | rand2 | 3.99999 | 3.999 / 4.001 | $-3.2\times10^{-13}$ / $-1.5\times10^{-11}$ | $-5.3\times10^{-4}$ |
| q0_consistent | $c=0$ | 4.00000 | 3.999 / 4.001 | $-4.4\times10^{-11}$ / $+3.6\times10^{-11}$ | $-6.4\times10^{-4}$ |
| coded_v2 | baseline | 3.99997 | 4.002 / 4.001 | $+2.5\times10^{-10}$ / $-3.1\times10^{-10}$ | $-3.1\times10^{-4}$ |
| coded_v2 | A/0.335 min | 4.00010 | 3.997 / 4.000 | $-5.5\times10^{-11}$ / $+2.4\times10^{-11}$ | $-6.7\times10^{-4}$ |
| coded_v2 | rand1 | 3.99997 | 4.001 / 4.001 | $+2.5\times10^{-10}$ / $-3.1\times10^{-10}$ | $-3.1\times10^{-4}$ |
| coded_v2 | rand2 | 4.00011 | 4.001 / 4.001 | $+2.4\times10^{-10}$ / $-3.0\times10^{-10}$ | $-3.2\times10^{-4}$ |
| coded_v2 | $c=0$ | 3.99999 | 3.999 / 4.001 | $+7.8\times10^{-11}$ / $-1.1\times10^{-10}$ | $-4.2\times10^{-4}$ |

- Method for the position: $|\det Q|^{1/4}$ is fitted by a quadratic in $\delta=E-E_f$ on 25 log-spaced points in $[2\times10^{-8},3\times10^{-5}]$ on each side; its root is $E_0$. This avoids the round-off region.
- Inside $|\delta|\lesssim2\times10^{-9}$, $\ln|\det Q|\approx-230$ is the round-off floor and the sign flips at 10 to 51 % of the 2001 points in $[E_f-2\times10^{-9},E_f+2\times10^{-9}]$ (208 to 1020 flips): any "sign-change zero" at $|\delta|\lesssim10^{-9}$ is noise on a quartic zero. This is why the scan's raw $E$ for the candidate moves by up to $9\times10^{-9}$ (`coded_v2`) between $\vec c$ values; the fit-based $E_0$ does not.
- The ordinary 110_A2 zeros in the same window move by $3\times10^{-4}$ to $9\times10^{-4}$ over the same $\vec c$ set, i.e. $\gtrsim10^{6}$ times more than the zero at $E_f$ ($\lesssim3\times10^{-10}$).

Verdict on D8: **artefact confirmed.** The zero at $E_f$ is a fourth-order zero ($|\det Q|\propto\delta^{4}$, exponent $4.00\pm0.003$ on both sides) fixed at $E_f$ to $\lesssim3\times10^{-10}$ (`coded_v2`) and $\lesssim10^{-10}$ (`q0_consistent`) for all five $\vec c$ including $\vec c=0$. It does not move with $\vec c$; D8's exclusion is correct here. Nothing to report loudly. Caveat: the check covers this one level (110_A2, $E_f=0.3194688$); the other 27 `free_level` candidates in `results/tester_m2gate_label_gate.json` are all within the `FREE_TOL`=$10^{-7}$ gate and have order 2.0 or 4.0 (one, 200_A2 0.34678, has estimated order 2.42, slope $-3\times10^{-9}$, not examined further).

## Problems / assumptions
- Order of the 200_A2 free-level candidate at 0.346780 is 2.42 (not an integer); it passes the gate (within $10^{-7}$ of a free level, no sign change) but its order estimate is unexplained (guess: a nearby pole contaminates the order estimate).

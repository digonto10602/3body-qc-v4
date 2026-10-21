# Report 007: M1 step 9, p-wave $\mathcal K_2^{-1}$ level-shift table (D7)

- **Date:** 2026-10-21
- **Git commit (HEAD when run):** `9e579e5`. The JSON under test records `git_commit = a25b1e1`, `git_tree_dirty = true`. Between `a25b1e1` and `9e579e5` only `scripts/m1/k2_pwave_shift.py` was added under `src/` and `scripts/` (`git diff --stat`), so the physics code is the same.
- **Prompt in force:** `prompts/001-m1-cache-free-single-energy-qc.md` (§3 row 9, §2.2).
- **Mode:** step.
- **Related:** `reports/tests/006-m1-step6-fast-path-v1.md`, `reports/tests/004-m1-step4-qc-lattice.md`.
- **Raw numbers:** `results/tester_k2_pwave_shift_check.json` (my oracle points, all 24 scans). No RNG is used.
- **Under test:** `scripts/m1/k2_pwave_shift.py`, `results/k2_pwave_shift_L20.json`.

## Commands

```
graphify update . ; graphify query "p-wave k2 q0_consistent coded_v2 refine_root scan"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 PYTHONPATH=src/python
python3 <inline>                      # check 1: JSON vs qc3.lattice.load_levels(default_data_dir(),20,IRREPS,0.335)
python3 <scratchpad>/oracle_zero.py   # check 2: ORACLE, all 12 rows x 2 forms (4 processes)
python3 <scratchpad>/pole.py          # check 3: refine_root on the +-250-row window, fast path
python3 <scratchpad>/scan481.py       # check 4: 481-point scan of the +-250-row window, 12 levels x 2 forms
.claude/bin/acceptance-check
```
The scratchpad is `/tmp/claude-0/-home-user-3body-qc-v4/756827e6-de12-5b97-9946-026a972abb7d/scratchpad`. These drivers are mine, not the builder's script. Nothing under `src/` or `tests/acceptance/` was touched. `_core` was not rebuilt.

## Results

| # | check | expected | got | tolerance | verdict |
|---|---|---|---|---|---|
| 1a | 12 rows, labels equal `load_levels(...).labels` in order | 12 | 12, identical (irrep, index) | exact | PASS |
| 1b | `lattice_Ecm_mean` vs `ls.Ecm`; `lattice_sigma` vs $\sqrt{\mathrm{diag}\,\mathrm{cov}}$ | equal | max diff 0 and 0 ($N_{\rm cfg}=603$) | exact | PASS |
| 1c | constants, `c_baseline`, commit, command, machine, python recorded | present | all present. `production_constants == PRODUCTION` after a JSON round trip. `c_baseline == C_BASELINE` | exact | PASS |
| 1d | internal arithmetic: `dE`, `dE/σ`, `coded−previous` | consistent | $\le10^{-15}$ | $10^{-12}$ | PASS |
| 2a | oracle: sign change of $\det Q$ across $E_{\rm table}\pm10^{-9}$, 12 rows × 2 forms (24 roots) | yes | 24 of 24 | — | PASS |
| 2b | oracle root (secant of the two $\pm10^{-9}$ points) vs table | $\le10^{-9}$ | max $\lvert\Delta\rvert=2.0\times10^{-12}$ (110_A2 n=2 q0). coded_v2 max $1.1\times10^{-12}$ (111_A2 n=0) | $10^{-9}$ | PASS |
| 2c | oracle dip: $\ln\lvert\det Q\rvert(E_0\pm10^{-5})-\ln\lvert\det Q\rvert(E_0\pm10^{-9})$ | $+\ln10^4=9.21$ for a simple zero; $<0$ for a pole | 9.03–9.21 for all 24, every one positive | $>3$ | PASS |
| 3 | `refine_root` in ±250 rows, q0, 110_A2 n=1 and 111_A2 n=0, converges to a pole | 0.316129 / 0.319554 | 0.31612911112995534 / 0.3195543911709816. $\ln\lvert\det\rvert$ rises by 2.30 per decade toward the root (110_A2, below the root at distances $10^{-5},\ldots,10^{-9}$: $-172.106,-169.803,-167.501,-165.198,-162.899$; above at $10^{-9}$: $-162.892$). Sign flips. This is $\propto1/\lvert E-E_p\rvert$, a pole | — | PASS (finding) |
| 4a | q0, 12 levels: validated zeros in the ±250-row window (481 points, spacing $5.06\times10^{-6}$ = 1.04 rows, 2.1× the builder's density) | 1 | 1 in all 12; equals the table (or scan-recovered value) to $\le2.0\times10^{-12}$ | $10^{-9}$ | PASS |
| 4b | same, coded_v2 (extra) | 1 | 1 in all 12; $\le1.1\times10^{-12}$ | $10^{-9}$ | PASS |
| 4c | sign-invisible zero or zero–pole pair closer than the spacing | none | one: 110_A2 n=2, both forms, at $E=0.3194688$ (details below) | — | FINDING, see below |
| 4d | nearest zero–pole distance | $>$ spacing | min $9.2\times10^{-5}$ (coded 111_A2 n=0), q0 min $1.2\times10^{-4}$ (100_A2 n=0); 18 spacings or more | $5.06\times10^{-6}$ | PASS |
| 5a | $\Delta E<0$ for all 12 levels | all negative | 7 of 7 plan-rule rows and 12 of 12 scan-assisted rows are negative. Range $-2.1\times10^{-8}$ (000_A1m) to $-2.97\times10^{-4}$ (111_A2 n=1). $\max\lvert\Delta E\rvert/\sigma_{\rm lat}=0.386$ | — | PASS (numbers only) |
| 5b | $\lvert E_{\rm coded}-E_{\rm prev.model}\rvert$ | $\le5\times10^{-8}$ | max $4.35\times10^{-8}$ (111_A2 n=1) | $5\times10^{-8}$ | PASS |
| 5c | oracle and fast agree at the two worst levels (111_A2 n=1, 110_A2 n=0) | exact | oracle root minus fast root $=1.3\times10^{-14}$ and $7.9\times10^{-15}$ | $10^{-9}$ | PASS |
| — | `.claude/bin/acceptance-check` | PASS | `PASS: tests/acceptance matches acceptance-v3` | — | PASS |

Check 2 covered all 24 (row, form) pairs, not only the 4 requested. That includes 110_A2 n=1, 200_A2 n=2, 100_A2 n=0, 111_A2 n=0 (the scan-recovered rows) and 000_A1m n=0, 110_A2 n=3 (plan rule worked). The oracle took about 1 min per row. For the 5 scan-recovered rows I used the builder's scan zero as $E_{\rm table}$.

Table of shifts. $\Delta E=E_{\rm q0}-E_{\rm coded}$ from `dE_scan_assisted`, which I reproduced. "plan" is yes when `refine_root` found the zero in the accepted window or the ±250-row window without the scan.

| irrep | n | plan | $E_{\rm coded}$ | $\Delta E$ | $\sigma_{\rm lat}$ | $\Delta E/\sigma$ |
|---|---|---|---|---|---|---|
| 000_A1m | 0 | yes | 0.2693227 | $-2.1\times10^{-8}$ | $5.15\times10^{-4}$ | $-4.1\times10^{-5}$ |
| 100_A2 | 0 | **no** | 0.2915561 | $-1.31\times10^{-4}$ | $5.05\times10^{-4}$ | $-0.259$ |
| 100_A2 | 1 | yes | 0.3004415 | $-1.25\times10^{-4}$ | $5.84\times10^{-4}$ | $-0.215$ |
| 110_A2 | 0 | **no** | 0.3067911 | $-1.56\times10^{-4}$ | $5.24\times10^{-4}$ | $-0.298$ |
| 110_A2 | 1 | **no** | 0.3158304 | $-1.85\times10^{-4}$ | $5.15\times10^{-4}$ | $-0.360$ |
| 110_A2 | 2 | yes | 0.3191549 | $-2.15\times10^{-4}$ | $7.70\times10^{-4}$ | $-0.280$ |
| 110_A2 | 3 | yes | 0.3287180 | $-7.21\times10^{-5}$ | $5.18\times10^{-4}$ | $-0.139$ |
| 111_A2 | 0 | **no** | 0.3194219 | $-1.78\times10^{-4}$ | $5.61\times10^{-4}$ | $-0.318$ |
| 111_A2 | 1 | yes | 0.3343321 | $-2.97\times10^{-4}$ | $7.70\times10^{-4}$ | $-0.386$ |
| 200_A2 | 0 | yes | 0.2862444 | $-1.89\times10^{-4}$ | $6.09\times10^{-4}$ | $-0.310$ |
| 200_A2 | 1 | yes | 0.2991783 | $-6.26\times10^{-5}$ | $5.83\times10^{-4}$ | $-0.107$ |
| 200_A2 | 2 | **no** | 0.3309818 | $-2.06\times10^{-4}$ | $6.11\times10^{-4}$ | $-0.337$ |

The 5 "no" rows are the ones where the plan rule (accepted window, then ±250 rows) gives `MissingRootError`. In 2 of them (110_A2 n=1, 111_A2 n=0) `refine_root` finds a sign change that is a pole (check 3). The builder's script catches this with its dip test. In the other 3 the ±250-row window has no sign change, because the shifted zero is outside it or the signs agree.

## Findings (not step-9 failures)

1. **`qc3.roots.refine_root` accepts poles (needs a fix).**
   - In the 110_A2 n=1 and 111_A2 n=0 q0 windows it returns a pole as if it were a root. Numbers are in check 3.
   - Its docstring says the caller must supply brackets that contain a zero. The step-9 script works around this with a dip test, but `qc3.roots` does not validate by itself.
   - A fitter or root finder that calls `refine_root` in widened windows can silently lock onto a pole.
   - Guess: add a post-check $\ln\lvert\det Q\rvert(E_r)<\min(\text{ends})-\delta$ (the builder's `DIP_MIN=3` works: all real zeros have dips 16–27) and raise `MissingRootError` otherwise.

2. **A hidden even-order zero at a free level in 110_A2, found by my 481-point scan.**
   - Where: $E^*=0.319468828$, inside the ±250-row window of 110_A2 n=2. It is the same for both forms.
   - What I saw: $\ln\lvert\det Q\rvert$ has a local minimum with the same sign on both sides ($+1$). With 25 points at $5\times10^{-7}$ spacing it is symmetric about the centre, $\lvert\det Q\rvert\propto(E-E^*)^4$ (fit slope $4.0$ in $\ln\lvert\det\rvert$ vs $\ln\lvert E-E^*\rvert$). The floor is $-211.16$ (q0) and $-211.83$ (coded).
   - Identification: this is exactly the non-interacting energy of $K(0)\,K(0,1,0)\,\pi(1,0,0)$ with $\mathbf P=(1,1,0)$ (my enumeration gives 0.319468828). It is a degenerate free level, and the projected $\mathcal K_{\rm df,3}$ and $\mathcal K_2$ do not lift it.
   - The sign-change scans cannot see it. `docs/ROOT_FINDING.md` §1, item 3 already lists even-order zeros as a known blind spot.
   - Is it a "second zero" that breaks the level assignment? Judged against the check-4 wording, no: it is not a zero–pole pair closer than the spacing, and it is 3.1e-4 above the coded n=2 zero and 5.3e-4 above the q0 n=2 zero, so it does not change the pairing of the sign-changing zeros. But it does change the QC level count (see the question below). I do not call FAIL.
   - Question for the coordinator and Digonto: are such free-level solutions to be counted as QC levels under D3 (sorted, lowest $m$)? The labelled windows (the ground truth, invariant 7) do not contain it. For 110_A2 it would be a 5th QC level between n=2 and n=3 below the cutoff, and D3 would then pair n=3 with it. The same blind spot exists in the legacy labels, so the table is consistent with them. Only that one hidden extremum occurred in the 24 scans at 1-row spacing; a quartic zero at a free level is easy to see at that spacing.

3. **Reader note (not a physics issue).** `lattice_Ecm_window_csv` (from `accepted_windows.csv`) differs from the D1 jackknife mean in 11 of 12 rows by $2.4\times10^{-5}$ to $1.2\times10^{-4}$ (up to $0.19\sigma$). Only the D1 mean is used in $\Delta E/\sigma$, so the table is unaffected.

4. **JSON records `git_tree_dirty = true` at `a25b1e1`** (the run happened before the commit). The source is unchanged from that commit (see header).

## Sign argument (check 5)

I verified only the numbers: all 12 $\Delta E<0$.

Reasoning I did not verify numerically: `q0_consistent` has $\mathcal K_2^{-1}\propto-1/a_1$ against `coded_v2` $\propto-1/(a_1q^{*2})$, so $\mathcal K_2^{-1}$ is smaller by $q^{*2}$ in the p-wave. A smaller $\mathcal K_2^{-1}$ means a larger p-wave $\mathcal K_2$ and a stronger p-wave interaction. HadSpec $a_1=+43.2\,a_t^3$ with $k^3\cot\delta_1=1/a_1>0$ is attractive (positive $\delta_1$), so the levels should shift down. That agrees with all signs but is not a test.

## Verdict

**PASS.** All checks 1, 2, 3, 4a, 4b, 4d and 5 pass. Checks 1 and 4 hold with the numbers above.
- Step 9 done-check: the JSON exists, has 12 rows, both forms and the lattice $\sigma$. Five rows have `dE = null` under the strict plan rule. The scan-assisted $\Delta E$ covers them, as the plan allows ("widen by ±250 rows and report").
- Two findings to pass on: `refine_root` accepts poles (finding 1), and the hidden quartic zero in 110_A2 (finding 2).

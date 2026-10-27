# Audit 002: M2 milestone (label-free roots, tracking, L20 fit)

- Date: 2026-10-27. Commit `294cc37` (diff `milestone-1..HEAD`, src/scripts/docs/prompts/tests m2+unit). Plan in force: `prompts/003` (v7 + 2026-10-26 as-built). Related: `reports/audits/001`, `reports/tests/014-019`.
- Auditor extra check: independent dense 1-row-grid sign-change scan (`qc3.roots.find_sign_changes`, split at jumps) vs `scan_sector` at the **fit minima** (q0_consistent, A/0.345: 000, 100, 200; A/0.335: 110, 111). Zeros and poles agree in all 5 sectors (zeros to <=4e-12, poles by count). Script: `/tmp/claude-0/-home-user-3body-qc-v4/756827e6-de12-5b97-9946-026a972abb7d/scratchpad/dense.py`. This covers the production form and non-baseline c, which the label gate (coded_v2, baseline c only) does not.

## VERDICT: PASS WITH NOTES (no blocking issue)

## Checklist
| item | status | where / note |
|---|---|---|
| Rootscan: sign changes, even-order via minima, 4b pairs | OK | `rootscan.py:238-266`; dense cross-check above |
| Free-level D8 rule and E snap | OK (note N2) | `rootscan.py:340-367`; snap justified (<=9.8e-10) |
| D9 slope sign recorded, never used to drop | OK | `rootscan.py:123-135`, `tracking.py:251`, `fit.py:149` |
| Jacobian dE0/dc_i = -u^T B~_i u / u^T Q' u | OK (note N5) | `tracking.py:102-117`; `L6=(xi*Lbyas)^6` = (L xi)^6 (`sector.py:38`); Q' at fixed c; B~ averaged at E0+-h |
| Tracking fallbacks loud/counted | OK (note N3) | `tracking.py:322-364`, counters in FitResult (all 0 full rescans in committed JSONs) |
| Zero entering from above E_max | OK | caught by `rescan`/`verify_fit` (`fit.py:489-538`); margin logic correct |
| chi2: D1 covariance (resample then E*, (N-1)/N) | OK | `lattice.py:50-82,122-132` |
| D3 matching incl. n>m, D8 exclusion, penalty, gradient | OK | `fit.py:100-190`; penalised rows J=0 is the exact gradient of the constant penalty |
| LM -> Migrad/Hesse; LM-only ambiguous rejection | OK (note N6) | `fit.py:405-412`; Migrad, start, minimum, verification still raise |
| Verification up to highest matched root + 0.002 | OK (note N4) | `fit.py:509-516` |
| Hesse vs FD-Hessian, error-gap flag | OK with caveat (N1) | `fit.py:293-333` |
| chi2/nu ~ 3.5 caveats | OK (N1) | `OPEN_ISSUES.md` #14 |
| float64, seeds, provenance | OK (N7) | fit JSONs: commit `9911237`, clean, core = commit, iminuit/numpy/scipy, machine, command, threads; deterministic (no RNG) |
| Tests exercise production path | OK | tester 019 (78/79; the 1 failure is the M1 anchor tolerance, out of scope) |

## Notes (all NON-BLOCKING)

**N1. Error reporting at 0.345 needs a stronger caveat (stats).** Tester 019: Hesse errors depend on where Migrad stops (chi2 piecewise smooth). B/0.345 Hesse K0 spans 1.1-3.5e5 under 1e-10 start perturbations at identical chi2; FD-Hessian stays 1.63e5. The committed B/0.345 has FD/Hesse = 1.199, just under the 1.2 flag threshold (`fit.py:317`), so `errors_disagree=False` is a lucky draw, not evidence of reliability. A/0.345: 2 of 5 perturbed starts reach a different local minimum (chi2=3876); no multi-start/global check exists for any cutoff, and parameters move a lot between cutoffs (A: K1 5.8e5 at 0.335 vs 9e4 at 0.345). Also: errors are nominal $\Delta\chi^2=1$ with $\chi^2/\nu=3.5$-4.6 (p ~ 2e-4 at A/0.335) and are not rescaled; several parameters have error ~ value. Fix: add the 019 findings to `OPEN_ISSUES.md` #14 and a `caveats` list in the result JSON `extra` (quote FD errors at 0.345; Hesse not reproducible there; no global-minimum guarantee; errors unscaled). Optional: print a WARNING when FD/Hesse is within 10% of the threshold at a cutoff known to be wall-limited, or add a 5-start spread to the flag.

**N2. Silent relabelling and dead band in `scan_sector`.** (a) Any zero-like candidate (including a sign-change `zero`) within `free_tol`=1e-7 of a free level becomes `free_level` and leaves `zeros()` (`rootscan.py:360-367`); it is kept in `free_level_solutions` but never printed. A physical zero crossing a free level during tracking would be dropped only in a scan (the verification count then mismatches, so it is loud there). (b) `edge`=2e-6 (0.4 row) next to every dimension jump and range end is not scanned (`rootscan.py:325`), and free levels within 3e-6 of a segment end are skipped (`:344-346`). A zero in the dead band is missed silently by the scan (a tracked one makes `verify_fit` fail the count). Tester 018 saw 200_A2 7e-6 above a jump at A/0.345. Suggested minimal fix: record `dead_band_hits` / `skipped_free_levels` in `ScanResult` and print them in `format_fit`; list `free_level_solutions` with slope sign in the fit JSON.

**N3. Full rescan replaces the tracked set without comparing** (`tracking.py:231-239`, called at `:308`). It is counted (`n_full_rescans`, 0 in all committed fits) but a root lost or gained by it is not reported. Suggested: log count change when `_load_full_scan` is entered from `_substep`.

**N4. `verify_fit` returns `ok=False` instead of raising** (`fit.py:536`). All committed fits have `ok=True` (max root diff <=2e-12), but a failing verification would only show as `ok=False` in the CLI line. Suggested: `format_fit` prints a WARNING, and the cutoff loop or the gate script asserts `ok`.

**N5. Jacobian accuracy near poles.** $Q'$ is a two-point central difference with $h=10^{-6}$; relative error $\sim(h/d_{\rm pole})^2$, i.e. $\sim2\times10^{-3}$ for a zero 5 rows (2.4e-5) from a pole (110_A2 at A/0.335). It affects the gradient and Hesse at the 0.1% level, not the minimum position or the verified chi2. FD checks in 016/017 were at benign points. Suggested: one FD-vs-analytic gradient check at a near-pole point (A/0.335, 110_A2).

**N6. LM-only rejection (Digonto-approved) is safe for the reported minimum.** The rejected trials only change which basin LM ends in; the reported chi2 is the true function at Migrad's minimum, verified by a fresh scan, with Migrad raising on any ambiguous candidate. Rejections are counted and printed (Table B: 1/2/0/1). It can bias the basin selection only (documented in #14). No action.

**N7. Provenance.** `results/m2_label_gate.json` records commit `da1c6fb` and core built at `c0334de` (older than HEAD); `src/` is unchanged since `9911237` and the core sources since `c0334de` (empty `git diff --stat`), so the numbers stand. The lattice input files' md5 are not recorded in the fit JSONs. Suggested: add `md5` of the `.jack` files used to `_document`.

**N8. Gate coverage.** The label gate (FP=FN=0, 207 brackets) is at baseline c with coded_v2 only; the q0_consistent form is gated by 12 oracle-verified roots plus my dense cross-check at the 5 fit minima. The 4b zero-pole-pair rule rests on a thin margin (dipole +0.65/-1.85 vs smooth second differences <=0.6, METHODS §13) and is blind for $s\lesssim2\times10^{-6}$; stated in the docs. Wide-c exploration by LM is only covered by the verification at the minimum.

## Answers to the focus questions
1. Rootscan classes missed: only the dead band (N2b), pairs with $s\lesssim2\times10^{-6}$ (stated), and two simple zeros closer than the fine step 1.2e-6 (would be `ambiguous`/loud). Dead band is silent in the scan, loud at verification.
2. Tracking: Jacobian formula, scaling and fixed-c $Q'$ are correct; fallbacks are counted; `ambiguous`/`even_zero` raise; zeros from above $E_{\max}$ caught by verification.
3. chi2: D1 resample-then-transform covariance, D3 incl. n>m, D8, penalty and gradient are consistent.
4. Fit: LM-only rejection cannot bias beyond basin selection; Hesse not reproducible at 0.345 (N1); chi2/nu ~ 3.5 is flagged in #14 ("do not quote as physics"), but errors are unscaled and no global search exists.
5. Provenance good (N7).

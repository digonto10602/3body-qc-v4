# Test report 008: M1 step 7.0 (projector cache key), 7.1 (pole/zero discrimination), 7.8 (D7 default q0_consistent)

- **Date:** 2026-10-22
- **Git commit:** `e5d72a5` (HEAD). The code under test is `00dbb32`; `git diff 00dbb32 e5d72a5 -- src tests` is empty (only graphify files differ).
- **Prompt in force:** `prompts/002-m1-step7-speed-replan-projector-basis-k3-g.md` (§2.5, §2.6, rows 7.0, 7.1, 7.8). Acceptance tag `acceptance-v4` (47 tests).
- **Mode:** step.
- **Related:** `reports/tests/006-m1-step6-fast-path-v1.md`, `reports/tests/007-m1-step9-pwave-shift.md`, `reports/research/008-step7-speed-replan-orbit-projector-and-kernels.md`.
- **Verdict: PASS.** 40 passed, 0 failed, 7 pending (all `NotImplementedError` naming 7.2 or 7.3), 0 new failures.

## Commands run (repo root, `OMP_NUM_THREADS=1`)

```
graphify update . ; graphify affected "refine_root" ; graphify affected "build_cached_projector_v30q"
mv src/python/qc3/_core*.so <scratch>/_core_old.so ; ./build.sh core            # clean rebuild, 1m50s
python3 -m pytest -q tests/acceptance/m1                                        # incl. slow, 358 s
python3 -m pytest -q tests --ignore=tests/acceptance                            # unit
.claude/bin/acceptance-check
python3 -m pytest --collect-only -q tests/acceptance
git diff 1e90c43 -- src/qc3_core/oracle_point.hpp src/core src/fitter src/cachegen_cpu | wc -l
# independent probes (scripts in the session scratchpad; outputs in results/tester_008_*.json)
stale.py  sameN.py  keyharness.cpp (g++ -O1, includes fast_projector.hpp)  strace -e trace=openat,write,read ...
disc.py   synth.py  real.py  xtol.py  k2_rerun.py (copy of scripts/m1/k2_pwave_shift.py writing to scratch)
```

## Summary table

| check | expected | got | tolerance | result |
|---|---|---|---|---|
| clean rebuild of `_core` | builds | builds, `.so` 2.0 MB | n/a | PASS |
| `tests/acceptance/m1` (47 collected) | 40 pass, 7 fail only with `NotImplementedError` naming 7.2/7.3 | 40 pass, 7 fail, all `NotImplementedError` (list below) | exact | PASS |
| unit tests | all green | 33 passed (18 s) | exact | PASS |
| `acceptance-check` | PASS (v4) | "PASS: tests/acceptance matches acceptance-v4" | exact | PASS |
| `pytest --collect-only tests/acceptance` | 47 | 47 | exact | PASS |
| oracle untouched | empty diff | 0 lines (`oracle_point.hpp`, `src/core`, `src/fitter`, `src/cachegen_cpu` vs `1e90c43`); oracle anchor tests in `test_m1_oracle.py`, `test_m1_fast_vs_oracle.py` green | exact | PASS |
| 7.0 key completeness: differential test, 130 cases | 0 stale hits | 0 mismatches (details below) | bitwise `V_sel`, $N_{\rm proj}$, sign, $\ln\lvert\det Q\rvert$ | PASS |
| 7.0 key completeness: C++ key harness, 17 variants | each changed input gives a new key | 15 new keys; the 2 "same" are the unchanged control and a no-op swap (both `n_x=-1`) | exact | PASS |
| 7.0 no disk I/O | none in cold and hit calls | 0 `openat/read/write/unlink/rename` syscalls between markers (strace) | exact | PASS |
| 7.1 discriminant, 26 points reproduced independently | kind as recorded, formula $\ln(\delta/h)$ | 24 zeros $\Delta_h\in[-26.61,-13.13]$, 2 poles $+11.91,+12.65$; 0 kind mismatches; recorded $\Delta_h$ equal to $<10^{-9}$ | kind exact; $\lvert\Delta_h\rvert\ge8$ | PASS |
| 7.1 consistency with $\delta=\lvert E^*-E_{\rm root}\rvert$ | $\lvert E^*-E_{\rm fine}\rvert=h\,e^{-\lvert\Delta_h\rvert}$ | ratio measured/implied $=0.9995\ldots1.0001$ (10 points with $\delta>10^{-13}$); smaller $\delta$ at the double-precision floor | 1e-3 | PASS |
| 7.1 12 `coded_v2` labelled zeros | all `zero` | 12/12 zero, $\Delta_h\in[-26.61,-13.73]$ | $\Delta_h\le-8$ | PASS |
| 7.1 `kind="any"` equals old `refine_root` | bitwise | 36 real windows (12 `coded_v2`, 12 `q0` accepted, 12 `q0` widened): identical $E$ (`repr`) and identical evaluation count; 11 of 12 synthetic cases identical | exact | PASS (one intended deviation, below) |
| 7.1 `DimensionJumpError` | raised iff $N$ differs within $\pm h$ | synthetic and real (100_A2 jump at $0.2653295153$, 110_A2 at $0.2721879708$): raised at $\pm0.5\,h$, not at $\pm1.5\,h$; synthetic edge $0.99h$ raised, $1.01h$ not | exact | PASS |
| step-9 re-run vs `47425c9` | root values unchanged | 50 numeric values (E_root, scan E, dE, dE_scan): max difference $0.0$ (bitwise); kinds added; same for my own re-run from the HEAD script | $\le10^{-11}$ | PASS |
| 7.8 defaults | `q0_consistent` | `Sector('100_A2')`, `Sector.production(...)`, `K2_PWAVE_FORM_DEFAULT` all `q0_consistent`; `PRODUCTION` dict lines unchanged in the diff | exact | PASS |
| 7.8 metadata | records the form | `metadata["k2_pwave_form"] == "q0_consistent"` (acceptance test green) | exact | PASS |
| 7.8 locked tests pin `coded_v2` | no silent cross-form comparison | `m1_helpers.sector` pins `coded_v2`; `test_m1_reference` explicit; `test_m1_fast_vs_oracle` explicit; all unit tests compare same-form quantities | n/a | PASS |

### Pending (7 failures, all `NotImplementedError`)

| test | message |
|---|---|
| `test_m1_fast_v2.py::test_projector_orbit_matches_eigen_and_oracle_subspace` | `projector='orbit' not implemented (prompts/002 step 7.2)` |
| `::test_projector_cold_build_is_fast` | `projector='orbit' ... (step 7.2)` |
| `::test_dimension_jumps_match_rule_and_N` | `qc3.fast.dimension_jumps: not implemented (step 7.3)` |
| `::test_fast_basis_all_signatures_vs_reference` | `qc3.fast.dimension_jumps ... (step 7.3)` |
| `::test_fast_tail_consistent_with_numpy_and_want_subset` | `projector='orbit' ... (step 7.2)` |
| `::test_fast_single_thread_cpu_time` | `projector='orbit' ... (step 7.2)` |
| `::test_fast_speedup_single_thread_largest_point_ge_150x` | `projector='orbit' ... (step 7.2)` |

The two tests assigned to this step pass: `test_projector_cache_key_complete`, `test_fast_options_projector_values`; so do `test_m1_roots.py` (4, including the slow `find_zeros` one) and `test_m1_production_default.py` (1).

## 7.0 details (cache key)

Code read: `fast_projector.hpp` keys on `basis_signature_v30q` (irrep, parity, $\mathbf P$, full $(\mathbf n,\ell,m)$ lists of both channels in order) plus `waves_1`, `waves_2`, `eig_tol`, `norm_tol`, `proj_tol` (17 digits) and the method. The builder `build_cached_projector_v30q` reads only these inputs plus two file-scope constants (`g_vsel_reorthonormalize`, `g_vsel_orth_tol`), which are not per-call state.

1. **Differential test** (`stale.py`, seed 20261022, `results/tester_008_cache_key_differential.json`). 130 inputs run in one process with a shared warm cache, each compared with a cold run (`clear_caches()` first). For every irrep: `eig_tol` in $\{10^{-30},10^{-3},0.3,0.49,0.05\}$, five waves variants, parity $\pm1$, `Lbyas` in $\{16,20,24\}$, `max_shell_num=10`, `epsilon_H=0.1`; plus 40 random mixtures. Result: 0 mismatches. 32 genuine hits, 86 entries, 12 inputs raise the same exception in both orders (`projection_produced_zero_columns` 7, "more than 4 (l,m)" 5). The test is sensitive: at `100_A2` 0.2915561, `eig_tol` $10^{-30}$ gives $N_{\rm proj}=5$ while 0.05 gives 17; parity $+1$ gives 2; waves variants give 11, 27, 13.
2. **Same basis lists, different irrep:** none occurs in the production range (0 shared basis arrays over 5 irreps $\times$ 40 energies), so a stale hit by irrep alone cannot arise from the lists; the C++ harness shows the key still changes with irrep and with $\mathbf P$ at fixed lists.
3. **`norm_tol`, `proj_tol` are not exposed in Python** (`Sector` and `module.cpp` do not pass them). They were varied in a small C++ harness (`keyharness.cpp`, compiled against `fast_projector.hpp`): $\texttt{norm\_tol}\times(1+10^{-7})$, $10^{-3}$, and the same for `proj_tol` each produce a distinct key and a separate cache entry (4 entries after 4 different builds; base again is a hit). `eig_tol` $0.05\to0.05+10^{-12}$ also gives a new key.
4. **No disk I/O:** `strace -f -e trace=openat,open,creat,write,read,unlink,rename` around a cold build plus a hit shows zero calls.
5. Not an issue but noted: `fast.clear_caches()` no longer clears the legacy `g_projector_cache_v30q`; the oracle clears it per call, so this is harmless. `cache_info` counters are cumulative (the unit test accounts for this).

## 7.1 physics review

**Derivation.** Near a simple zero $E_0$, $\ell(E)=\ln\lvert E-E_0\rvert+b(E)$, $b$ smooth. With $\delta=E^*-E_0$ and $h=10^{-6}$:
$$
\Delta_h=\ln\lvert\delta\rvert-\tfrac12\ln\lvert\delta^2-h^2\rvert+O(h^2b'')=\tfrac12\ln\frac{\delta^2}{h^2-\delta^2}\ \xrightarrow{\ \delta\ll h\ }\ \ln\frac{\delta}{h}\ <0 .
$$
For a simple pole the sign flips: $\Delta_h\approx+\ln(h/\delta)$. The linear part of $b$ cancels exactly in the symmetric stencil; the remainder is $O(h^2b'')\sim10^{-12}$. With Brent's $\lvert\delta\rvert\le x_{\rm tol}=10^{-11}$ the weakest value is $\lvert\Delta_h\rvert\ge\ln(h/x_{\rm tol})=11.51$. A pure-$\ln$ synthetic check agrees with the exact expression to $10^{-8}$ relative for $\delta=10^{-12}\ldots9\times10^{-7}$ (e.g. $\delta=10^{-11}$: $-11.512925$ vs $-11.512925$).

**Measured values** (`results/tester_008_discriminants.json`, `disc.py` re-evaluates $\ell$ at $E^*,E^*\pm h$ with `fast.evaluate_point`, not through `qc3.roots`): the 12 `coded_v2` zeros give $-13.73\ldots-26.61$; the 14 `q0_consistent` entries (12 zeros, 2 poles) give zeros $-13.13\ldots-25.03$ and poles $+11.91$ (111_A2 0.3195544) and $+12.65$ (110_A2 0.3161291). Overall zeros $-13.13\ldots-26.61$ (the builder's quoted lower end $-13.2$ is rounded; the actual smallest magnitude is $13.13$), poles $+11.91,+12.65$. All $\ge11.5$, as required.

**Consistency with the convergence of Brent.** I located each singularity independently by bisection on the sign of $\det Q$ down to the double spacing ($5.5\times10^{-17}$). The measured $\lvert E^*-E_{\rm fine}\rvert$ equals $h\,e^{-\lvert\Delta_h\rvert}$ to $0.05\%$ for every point where $\delta>10^{-13}$ (e.g. 100_A2 q0, $\delta=1.871\times10^{-12}$ measured vs $1.87\times10^{-12}$ implied; the pole 111_A2: $6.743\times10^{-12}$ vs $6.74\times10^{-12}$), and all $\delta\le6.7\times10^{-12}<x_{\rm tol}$. So the discriminant is a direct measurement of $\ln(\delta/h)$, and the observed spread (13 to 27) is just the spread of Brent's actual accuracy between $10^{-12}$ and $10^{-17}$.

**Limit when a zero and a pole are close (synthetic, `synth.py`).** A second singularity at distance $d$ of the opposite type adds $\mp f(d)$ with $f(d)=-\tfrac12\ln\lvert1-h^2/d^2\rvert$. For $d\ll h$ the result is $\Delta_h\approx\ln(\delta/d)$ for a zero and $\ln(d/\delta)$ for a pole. Hence misclassification or ambiguity needs $d\lesssim e^{4}\delta\approx55\,\delta\le5.5\times10^{-10}$, independent of $h$ (the stencil sees the other singularity only through $f$). Measured with the zero at $E_0$, the pole at $E_0+d$ and the bracket isolating the target (Brent $\delta\approx10^{-12}$): $d\ge3\times10^{-10}$ classifies correctly for both kinds ($\Delta_h=-6.87$ zero, $+4.14$ pole); $d=10^{-10}$ the zero is still correct ($-4.17$) and the pole is `ambiguous` ($+2.66$); for $d\le3\times10^{-11}$ both are `ambiguous`. No case classified wrongly (always `ambiguous`, never the opposite kind, down to $d=10^{-13}$). For every $d\ge3\times10^{-10}$, including $d=h(1\pm10^{-3})$ and $d=h$ exactly (stencil on the other singularity), all cases are classified correctly. Stated limit: the rule is reliable for opposite-type singularities separated by more than about $55\,x_{\rm tol}\approx5.5\times10^{-10}\ (\approx10^{-4}$ of a grid row, grid step $4.845\times10^{-6})$; closer pairs are reported ambiguous. Two same-type roots within $\approx10^{-13}$ of exactly $h$ could in principle be misread (the stencil hits a singularity), a measure-zero case; with a real evaluator that point raises `ValueError` ("Q is singular", see below), so it is loud.

**Other synthetic findings.**
- Required `xtol`: the rule needs $x_{\rm tol}\lesssim h\,e^{-4}=1.8\times10^{-8}$. With $x_{\rm tol}=10^{-7}$ a synthetic pole came out `ambiguous` ($\Delta_h=3.10$) and with $10^{-6}$ also ($0.43$). Not a bug (loud `AmbiguousRootError`), but nothing guards or documents it for callers who loosen `xtol`.
- The branch for $\ell(E^*)=-\infty$ in `_classify` is unreachable with the real evaluators: `qc3.qc.slogdet` raises `ValueError("Q is singular (zero pivot in LU)")` on an exactly singular $Q$. It is loud, not silent. My own bisection hit it once (at a representable root), so it is a real, if rare, path.
- `DimensionJumpError`: synthetic $N$ jump at $E^*+x$: raised for $x\in\{-0.5,-0.1,0,0.1,0.5,0.99\}\times10^{-6}$ and not for $-1.5\times10^{-6}$, $+1.01\times10^{-6}$, $+1.5\times10^{-6}$. Real: jump at 0.26532951530559 (100_A2, $N$ 34$\to$50) and 0.27218797076179 (110_A2, 28$\to$36), located by bisection; `classify_sign_change` raises at $E^*=$ jump $\pm0.5h$ and does not at $\pm1.5h$ and $\pm3h$.
- **`kind="any"` vs old `refine_root`:** 36 real windows, the returned $E$ is identical to the last digit and the call count is identical. One intended deviation: if a bracket end point is *exactly* a pole ($\ell=+\infty$), the old code raised `ValueError` ("function value is NaN") and the new one handles it (maps $l_0$ to 0 and $\pm\infty$ to $\pm10^{300}$ / 0). It occurred in 1 of 12 synthetic cases (2 vs 5 evaluations); it cannot arise with the real evaluator because `slogdet` raises first.

## Step-9 script re-run (check 5)

`results/k2_pwave_shift_L20.json` (HEAD) vs `git show 47425c9:results/k2_pwave_shift_L20.json`: 50 numeric root and shift values (`E_root`, scan zeros `E`, `dE`, `dE_scan_assisted`, 12 levels, both forms) differ by at most $0.0$ (limit $10^{-11}$). Added keys are `kind`, `discriminant` (and the scan discriminants). The recorded kinds: `coded_v2` 12 $\times$ zero; `q0_consistent` 7 $\times$ zero, 5 $\times$ no root by the plan rule (unchanged). The pole notes now read "PoleError ... $\Delta_h$ = 12.65 / 11.91" instead of the old ad-hoc dip message. I also re-ran a copy of `scripts/m1/k2_pwave_shift.py` (output to scratch, 154 s): all 50 values and all discriminants are bit-identical to the committed HEAD file (`results/tester_008_step9_rerun_vs_47425c9.json`). `n_evaluations` per root went from 8 to 10 (the 3 classification calls minus 1 re-used end point), as the plan expects.

## 7.8 details (D7)

- `Sector('100_A2').k2_pwave_form`, `Sector.production('100_A2').k2_pwave_form`, `constants.K2_PWAVE_FORM_DEFAULT`: all `q0_consistent`. `oracle.K2_PWAVE_FORMS == ("coded_v2","q0_consistent")` unchanged. The `PRODUCTION` dict has no changed lines in `git diff 47425c9 HEAD`.
- Fast-path result metadata contains `k2_pwave_form` (asserted by the locked test); the Python-level `PRODUCTION` dict does not contain the form, so the record lives in `PointResult.metadata`, as in the plan.
- The differences matter: with the default (q0) vs pinned `coded_v2`, $\ln\lvert\det Q\rvert$ at `c=0` differs by $0.011$ to $1.14$ at four small anchor points (e.g. 100_A2 0.2915561: $-105.991$ vs $-104.849$), so any silent cross-form comparison with the anchors (`coded_v2`) would fail well outside the $10^{-6}$ tolerances. Because the 40 locked tests pass, none does this. The pinning chain: `m1_helpers.sector()` uses `setdefault("k2_pwave_form","coded_v2")`; `oracle_point`/`fast_point` helpers and every locked comparison go through it; `test_m1_reference` passes the form explicitly in both branches; `test_m1_fast_vs_oracle::test_fast_root_vs_oracle_100_A2_windows` pins it; `test_m1_oracle` asserts the pinned `coded_v2` in metadata; the new roots test uses `q0_consistent` for the step-9 cases and `coded_v2` for the 12 labelled zeros (both explicit). Unit tests that use the new default compare backend against backend with the same sector, so they are form-consistent. No acceptance file was edited by me (`acceptance-check` PASS).

## Failures

None.

## Observations for the coordinator (not failures; guesses labelled)

1. No guard on `xtol` versus $h\,e^{-D}$ in `refine_root_info` (a caller with `xtol` $\gtrsim10^{-7}$ gets `AmbiguousRootError` on valid roots). Loud, so acceptable; a docstring note would help (guess: the builder assumed the default 1e-11).
2. The `l0 == -inf` branch in `roots._classify` is unreachable with the real evaluators (they raise `ValueError` on a singular pivot). Harmless dead code.
3. `norm_tol` and `proj_tol` cannot be set from Python, so the key coverage for them rests on code reading plus my C++ harness; the acceptance test cannot exercise them.
4. `PRODUCTION` does not contain the p-wave form (it is in `K2_PWAVE_FORM_DEFAULT` and in `metadata`). If invariant 3's "record the constants in every result's metadata" should include it for stored JSON, step 8's `bench.run` must write it (the plan already says so).

## Raw numbers

`results/tester_008_summary.json`, `tester_008_cache_key_differential.json`, `tester_008_discriminants.json`, `tester_008_synthetic_roots.json`, `tester_008_real_roots_checks.json`, `tester_008_xtol_scan.json`, `tester_008_step9_rerun_vs_47425c9.json`, `tester_008_anchor_form_check.json` (commit `00dbb32` code, seed 20261022 where random, machine Linux 6.18.44 x86_64, Python 3.13.16, 1 thread).

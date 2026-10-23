# 001 — M1 milestone audit: cache-free single-energy QC, D1 covariance, D3 matching, D7 p-wave form

- **Date:** 2026-10-22
- **Git commit audited:** `2fa090f` (HEAD of `claude/gifted-newton-6ys9qm`); the compiled `qc3._core` reports `git_commit() = 5792545`, and `git diff 5792545..HEAD -- src` is empty, so the binary audited is the source audited.
- **Commit range:** `cd16663..2fa090f`, restricted to `src scripts tests docs prompts` (49 files, +6548/−11).
- **Prompts in force:** `prompts/001-m1-cache-free-single-energy-qc.md` (steps 1–6, 8, 9, gate), `prompts/002-m1-step7-speed-replan-projector-basis-k3-g.md` (step 7, incl. the §4 amendment for acceptance-v4; v2/v3 amendments in `prompts/001` §4).
- **Acceptance tag:** `acceptance-v4`; `.claude/bin/acceptance-check` → `PASS: tests/acceptance matches acceptance-v4`.
- **Related reports:** tests `reports/tests/001`–`012` (esp. 006, 007, 008, 010 N1/N2, 011, 012), research `reports/research/007`, `008`.
- **Auditor scratch scripts (not in git):** `/tmp/claude-0/-home-user-3body-qc-v4/756827e6-de12-5b97-9946-026a972abb7d/scratchpad/aud/offprod.py`, `.../aud/k2check.py`, output `.../aud/offprod.log`.

## Verdict: PASS WITH NOTES

I found no blocking physics or numerics error. The fast path computes the same equations as the oracle, which is the existing C++ code. Invariants 1 and 2 hold. The D7 p-wave form matches HadSpec's published convention, which I checked against the paper itself. The D1 covariance and the D3 matching are as decided. There are 3 guard gaps (non-blocking, none reachable with production inputs), 2 defective acceptance tests that need a planner fix approved by Digonto, and some stale docs and provenance items.

## 1. What I did (commands)

1. Read STATUS, prompts 001/002, METHODS, ROOT_FINDING, OPEN_ISSUES, and tester reports 010 and 012. Read every new C++ header (`src/qc3_core/fast_{point,fclass,pairs,g,k3,projector,basis}.hpp`, `k2inv_forms.hpp`, `oracle_point.hpp`, `module.cpp`) and every Python module (`qc3/{fast,oracle,sector,constants,point,qc,roots,lattice,bench}.py`). Compared them line by line with the oracle formulas: `src/core/K2_functions_v2.h`, `K3_functions_2plus1.hpp` (`K3B_element`, `K3E_element`) and `fv_projector_cartesian_l1_v30m.hpp`.
2. `git diff --quiet cd16663..HEAD -- src/core src/fitter src/cachegen_cpu` returned true: the oracle's physics sources are byte-unchanged.
3. **Independent $\mathcal K_2^{-1}$ check** (`aud/k2check.py`). I re-derived every diagonal entry from the formula in NumPy, for both forms, both channels, at 4 points (100_A2 0.2915561, 111_A2 0.30, 000_A1m 0.27, 200_A2 0.345): $\eta/(8\pi\sqrt{\sigma})\,[B_\ell+C_\ell]/(2\omega L^3)$, with a factor 2 for channel 2. Maximum relative deviation from `qc3.fast`: $2.3\times10^{-16}$ (`q0_consistent`), $4.6\times10^{-15}$ (`coded_v2`, at a point with $q^{*2}\to0$).
4. **HadSpec check.** I downloaded arXiv:1411.2004 (PDF, 29 pages) into the scratchpad and read it with `pdftotext`. Eq. (9) is $k^{2\ell+1}\cot\delta_\ell=\frac{1}{a_\ell}+\frac12r_\ell k^2+O(k^4)$. The Sec. VI $I=3/2$ fit gives $a_0=(-4.04\pm0.05\pm0.15)\,a_t$ and $a_1=(43.2\pm3.7\pm15.4)\,a_t^3$.
5. **Off-production fast vs oracle** (`aud/offprod.py`, 1 thread). Cases: parity $+1$; `waves_1=(0,)`; $L/a_s=24$; $\alpha=0.8$; $E^*=0.26320$ (just above threshold); $E^*=0.3599$ (top of range, $N=160$); `q0_consistent`; `waves_2=(0,1)`. Each case compares $F,G,\mathcal K_2^{-1},\mathcal K_3,F_3^{-1}$ (relative Frobenius norm), $VV^\dagger$, the sign, and $\Delta\ln|\det Q|$ for default options and for `projector="eigen", use_classes=False`.
6. **Independent check of the D1 covariance.** I linearised it from the raw per-configuration values: $C\approx\frac{1}{N(N-1)}J\,\sum_c(x_c-\bar x)(x_c-\bar x)^TJ$ with $J=\bar x/E^*$. This agrees with `qc3.lattice.load_levels(...,0.335).cov` to $1.06\times10^{-5}$ relative. The difference is the expected nonlinear term. Also $\sigma(\text{000\_A1m n0})=5.153\times10^{-4}$, $\kappa(C)=1.81\times10^3$, $N_{\rm cfg}=603$.
7. Ran `test_fast_single_thread_cpu_time` alone, with and without `OPENBLAS_NUM_THREADS=1`.

## 2. Checklist

| item | status | where |
|---|---|---|
| Equations match METHODS term by term ($F$, $G$, $\mathcal K_2^{-1}$, $\mathcal K_{\rm df,3}$, $F_3$, $Q$) | OK | `fast_point.hpp:241-250` ($F$: $(2\pi/L)^{\ell'+\ell}S-\delta I_\ell$, prefactor $H/(16\pi^2L^4\omega(E-\omega))$); `fast_g.hpp:185-205` (flavour factors $\sqrt2(-1)^{\ell'}$, $\sqrt2(-1)^{\ell}$, zero (2,2) block); `fast_k3.hpp:116-243` vs `K3_functions_2plus1.hpp:140-430` (every branch transcribed, coefficients $-2/3,-1/3$, $c_{ij}$ from `two_plus_one_block_factor`); `fast_point.hpp:323-328` ($F_3=F/3-FH^{-1}F$) |
| $\mathcal K_2^{-1}$ p-wave (D7) vs HadSpec | OK | `k2inv_forms.hpp:236-247`; see §3.2 |
| Invariant 1 (project the inverse, $L^{-6}$ on the matrix, no $N_{\rm proj}$ factor) | OK | `fast_point.hpp:371-393`: $X=(F_3/s)^{-1}V$, $Q=V^T(X/s+\mathcal K_3V)/L^6$; `kin.L = par.L()` $=\xi L/a_s$ (`fast_basis.hpp:46`) |
| Invariant 2 (LU sign, $\ln|\det|$) | OK | `fast_point.hpp:263-278` (`PartialPivLU`, permutation sign × pivot signs; throws on a zero pivot) |
| Basis / index conventions ($m\leftrightarrow$ axis, ordering, flavour blocks) | OK | $m=(-1,0,1)\leftrightarrow(y,z,x)$ in `Yreal` (`fast_point.hpp:73-78`), `Yq0` (`fast_pairs.hpp:100-104`) and the class path (`fast_point.hpp:182`); $\mathcal K_E$ tensor index map `(I+2)%3` identical to the oracle (`fast_k3.hpp:210-219`); basis bit-identical to `config_maker_4_momentum_first` (`fast_basis.hpp:120-163`; tester 009: 1881 energies) |
| Boundary conditions / wrap-around / off-by-one | OK | ellipsoid box with $\pm1$ slack plus an exact $r^2>R^2$ test (`fast_point.hpp:200-224`); class tables cover the same cube $[-N,N]^3$ (`fast_fclass.hpp:74-90`); threshold table in the oracle's loop order |
| Moving frames (boost of the lattice sum) | OK | $\mathbf r=\mathbf n+\mathbf n_P\,[(\mathbf n\cdot\mathbf n_P/n_P^2)(\gamma^{-1}-1)-\xi/\gamma]$ (`fast_point.hpp:219-221`), $r^2=r_\perp^2+r_\parallel^2$ class split (`fast_fclass.hpp:172-186`), signed-permutation back-rotation algebra checked (`fast_point.hpp:149-177`); off-production $L/a_s=24$, $\alpha=0.8$ agree to $\le6\times10^{-16}$ in $F$ |
| Real-arithmetic re-expression | OK in production, guards incomplete | ISSUE A, ISSUE B |
| Projector (orbit/character) | OK | `fast_projector.hpp:72-179`: same `little_group`/`chi`/`cubic_transf`/`D_real_element_v30m`/`parity_factor_for_l(…,1)` as convention c0. Always-on checks $\|PV-V\|$, $\|V^TV-1\|$ and **$\mathrm{tr}P=N_{\rm proj}$** (this one detects a missing candidate). Parity $+1$ 100_A2: $VV^\dagger$ vs oracle $4.5\times10^{-16}$ |
| Truncation ($j_{\max}$-type cut-offs) | OK (note) | `max_shell_num=20`, `trunc_eps=1e-16`; per-term bound only, see Note N3 |
| Normalisation / reality | OK | $V^TV=1$ is checked; the determinant phase is checked in `fast.py:66-67`; $\mathcal K_2^{-1}$ imaginary part is checked (but see ISSUE A) |
| Precision (float64 reference, accumulators) | OK | class path uses `double` accumulators (`fast_fclass.hpp:35-38`), direct path `long double`; $F$ vs oracle $\le5.4\times10^{-15}$ (tester 012) |
| Randomness / seeds / error bars | OK | test seed `20261021` (`test_m1_fast_vs_oracle.py:87`); jackknife errors in `qc3.lattice` |
| Tests exercise the production path | OK | `m1_helpers.fast_point` uses default options (`use_classes=True`, `projector="orbit"`, ordered-pair $G$) |
| Tests really test the claim / tolerances | 2 ISSUES (test defects) | ISSUE C, ISSUE D; notes N4, N5 |
| No hard-coded answers in source | OK | golden values live only in tests (`test_m1_roots.py:22-32`) and `results/` |
| Statistics: D1 jackknife, covariance, pseudo-inverse $\chi^2$ | OK | `lattice.py:57-76,93-138,155-177`; independent linearised check in §1 item 6 |
| D3 matching | OK | `lattice.py:141-152`; χ² proxy pairing is trivially the identity (`bench.py:204-214`) |
| Roots: $\Delta_h$ discriminant, `DimensionJumpError`, `kind="zero"` default | OK (note) | `roots.py:58-73,130-150`; Note N6 |
| D8 recorded | OK | `docs/METHODS.md` §11, `docs/ROOT_FINDING.md` §4, `prompts/002` §9; no M1 code (correct, it is M2 work) |
| Provenance (commit, constants, form) | OK with notes | Note N7 |

## 3. Findings

### 3.1 Fast path vs oracle (priority 1)

**Agreement outside the production settings** (`aud/offprod.log`, relative Frobenius norm unless stated):

| case | $N$ / $N_{\rm proj}$ | $F$ | $\mathcal K_3$ | $F_3^{-1}$ | $VV^\dagger$ | $\Delta\ln|\det Q|$ | sign |
|---|---|---|---|---|---|---|---|
| 100_A2 0.30, parity $+1$ | 60 / 2 | $2.1\times10^{-16}$ | $7\times10^{-18}$ | $1.6\times10^{-15}$ | $4.5\times10^{-16}$ | $1.8\times10^{-15}$ | equal |
| 110_A2 0.31, waves$_1=\{0\}$ | 32 / 16 | $1.2\times10^{-14}$ | $7\times10^{-17}$ | $6\times10^{-16}$ | $5.8\times10^{-16}$ | $1.4\times10^{-13}$ | equal |
| 111_A2 0.30, $L/a_s=24$ | 112 / 30 | $6\times10^{-16}$ | $7\times10^{-17}$ | $4\times10^{-15}$ | $5\times10^{-16}$ | $1.1\times10^{-13}$ | equal |
| 200_A2 0.32, $\alpha=0.8$ | 99 / 26 | $5\times10^{-16}$ | $4\times10^{-17}$ | $9\times10^{-15}$ | $6\times10^{-16}$ | $2.0\times10^{-13}$ | equal |
| 100_A2 0.26320 (threshold) | 34 / 11 | $3\times10^{-16}$ | $2\times10^{-17}$ | $5\times10^{-16}$ | $3\times10^{-16}$ | 0 | equal |
| 111_A2 0.3599 (top, $N$ max) | 160 / 40 | $6\times10^{-16}$ | $9\times10^{-17}$ | $3.5\times10^{-13}$ | $6\times10^{-16}$ | $2.8\times10^{-13}$ | equal |
| 110_A2 0.33, `q0_consistent` | 92 / 36 | $6\times10^{-16}$ | $3\times10^{-17}$ | $6\times10^{-16}$ | $6\times10^{-16}$ | $1.9\times10^{-10}$ | equal |
| 000_A1m 0.33, parity $+1$ | — | oracle `projection_produced_zero_columns`; fast `ZERO_PROJECTED_DIM` | | | | both loud | |
| 000_A1m 0.34, waves$_2=\{0,1\}$ | — | **oracle `BAD_F3_SCALE`; fast returns a finite $\ln|\det Q|=-75.89$** | | | | ISSUE A | |

$G$ and $\mathcal K_2^{-1}$ were bit-identical (0.0) in every case that ran. The near-threshold point has imaginary $q^*$ for many spectators; the top point has the largest boosts. Other irreps cannot be reached: `parse_label` accepts only the 5 production labels (`src/fitter/debug_v30r_…zero_scan.cpp:66-74`). The $d_\Gamma>1$ projector branch is therefore unreachable and untested, apart from its self-checks.

**Is anything correct only in the production range?** The re-expression of the pair kernels in real arithmetic rests on $\sigma_k>\sigma_{\min}>0$ for every basis momentum. The threshold basis guarantees this, so $\beta<1$, $\gamma$ is real, and $q^{*2}+m^2\ge0$ for kaon spectators. The pion-spectator `pms_0`, which is imaginary below the $KK$ threshold, is correctly never used: I checked every branch of `K3B_element`/`K3E_element`, and the (2,1) branches use only the kaon's `kms_0`. The guards are almost complete, apart from the gaps below:

#### ISSUE A — $\mathcal K_2^{-1}$ guard is NaN-blind and $\mathcal K_2^{-1}$ is not checked for finiteness (NON-BLOCKING)
- **Where:** `src/qc3_core/fast_point.hpp:311-316` and `:351`.
- **What:** `if (sc > 1e-14 * scr) throw …` with `sc = max|Im K2inv|`. If an entry is $-\infty$, complex multiplication gives $\mathrm{Im}=\mathrm{NaN}$. Then `NaN > x` is false and the check passes. `K2inv.allFinite()` is never tested (the oracle tests it, `oracle_point.hpp:247`).
- **What happens:** $H$ gets a $-\infty$ diagonal entry, LU effectively decouples that row, and $F_3$ comes out finite. The fast path then returns a number where the oracle stops.
- **Reproduced:** `Sector.production("000_A1m", waves_2=(0,1))` at $E^*=0.34$. The $KK$ p-wave has $a_1=0$, so $\mathcal K_2^{-1}=-\infty+\mathrm{NaN}\,i$ on 81 rows. The oracle gives `BAD_F3_SCALE`; fast returns $\ln|\det Q|=-75.89$ silently.
- **Why not blocking:** production cannot reach it. Waves$_2=\{0\}$ is invariant 3, and odd $\ell$ for identical $K^+K^+$ is forbidden by Bose symmetry. The only other trigger is `coded_v2` at exactly $q^{*2}=0$, where $-1/(a_1q^{*2})$ is infinite; that is a measure-zero energy.
- **Fix (one line each):** `if (!K2inv.allFinite()) throw std::runtime_error("fast: NONFINITE K2inv");` and write the imaginary-part test as `if (!(sc <= 1e-14 * scr))`. Also add `H`/`X1` to the finiteness check at `:351`. **Test:** the `waves_2=(0,1)` call above must raise.

#### ISSUE B — two places drop imaginary parts with `.real()` and no check (NON-BLOCKING)
- **Where:** `fast_pairs.hpp:73-85` (`make_boost_par`: `beta2`, `gam2 = 1/sqrt(1-beta2^2)` in `comp`, then `.real()`) and `fast_point.hpp:113-125, 246-248` ($A$, $\gamma$, $x^2$, $\xi/\gamma$, $I_\ell$).
- **Problem:** the header comment `fast_pairs.hpp:14-15` claims a negative square root "would yield NaN … cannot fail silently". That holds for `double` sqrt, but these are complex square roots, so a non-real value would be truncated silently.
- **Production impact:** none, because $\sigma>\sigma_{\min}>0$ makes all these quantities real.
- **Fix:** `check_real(gam2,"boost gamma")` in `make_boost_par`. In `F_channel`, a per-momentum relative check $|\mathrm{Im}\,A|,|\mathrm{Im}\,I_\ell|\le10^{-14}|\cdot|$ (cost: one comparison per momentum).

Everything else checked out:
- **$F$ class path.** Exact algebra; same point set as the direct sum. Agrees with the oracle to $\le5.4\times10^{-15}$ and per block to $\le4.4\times10^{-14}$ (tester 012).
- **Projector.** The trace check makes an incomplete orbit candidate set impossible to miss.
- **Tail.** $FX_1$ uses $F$'s block-diagonal structure, which is exact because the groups partition each channel by momentum.
- **Real/complex branch.** Real Hermitian $V$ uses $V^T$; complex $V$ (`eigen`) uses $V^\dagger$.

### 3.2 p-wave $\mathcal K_2^{-1}$ (D7) (priority 2)

In the $Q_0$ convention, $F$ and $G$ carry $|q^*|^\ell$ in their harmonics (METHODS §3–4), so the matching two-body element is $q^{*2\ell}\mathcal K_2^{-1}$:
$$
[\mathcal K_2^{(i)}]^{-1}_{\ell}=\frac{\eta_i}{8\pi\sqrt{\sigma}}\Big[q^{*2\ell+1}\cot\delta_\ell+q^{*2\ell}|q^*|(1-H)\Big]
\quad(\text{arXiv:2111.12734, Eqs. 15, 17, 20--21}).
$$

**Convention.** HadSpec arXiv:1411.2004 Eq. (9) is $k^{2\ell+1}\cot\delta_\ell=1/a_\ell+\dots$. The code uses $-1/a$, so the code's $a_\ell$ is $-a_\ell^{\rm HadSpec}$:
- s wave: $q^*\cot\delta_0=-1/(+4.04)=1/(-4.04)$ ✓
- p wave (`q0_consistent`): $q^{*3}\cot\delta_1=-1/(-43.2)=+1/43.2=1/a_1^{\rm HadSpec}$ ✓
- units: $a_t^{-3}$ on both sides ✓

`k2inv_forms.hpp:245` has $B=-1/a_1+\tfrac12r_1q^{*2}+P_1q^{*4}$ and returns $A(B+q^{*2}|q^*|(1-H))$. This is identical to the commented-out line `K2_functions_v2.h:103`.

**`coded_v2`** (`K2_functions_v2.h:102-105`) has $B=-1/(a_1q^{*2})+\tfrac12r_1+P_1q^{*2}$. That is exactly $1/q^{*2}$ times the `q0_consistent` scattering-volume term, so the factor-$1/q^{*2}$ claim is **confirmed**. The form is also inconsistent in units ($a_t^{-1}$ against $a_t^{-3}$ for the threshold term) and is singular at the $K\pi$ threshold $q^{*2}=0$. That gives two more reasons, besides HadSpec, for D7.

**s wave and channel 2 are unchanged.** `K2_inv_ERE_ang_mom_q0_consistent` hands every non-$(\ell=1,\ell'=1,m=m')$ entry to the unchanged `K2_inv_ERE_ang_mom`. The mass arguments $(m_1,m_1,m_2)$ / $(m_2,m_1,m_1)$ and the factor $2$ on channel 2 match `K2inv_EREord2_2plus1_mat`. My independent formula check (§1 item 3) reproduces both forms in both channels to $\le4.6\times10^{-15}$. The locked `test_K2inv_pwave_both_forms_vs_reference` asserts that all non-p-wave rows (channel 2 included) are identical.

**`coded_v2` is still available and pinned correctly:**
- `m1_helpers.sector` uses `setdefault("k2_pwave_form","coded_v2")`.
- `test_m1_reference`, `test_m1_roots` and `test_fast_root_vs_oracle_100_A2_windows` pin the form explicitly.
- `test_m1_production_default.py` bypasses the helper.
- The χ² proxy pins `coded_v2` (`bench.py:42,192`).

Minor: the C++ `Input` default is still `"coded_v2"` (`oracle_point.hpp:56`). Python always passes the form explicitly, so this matters only to a C++-only caller (Note N8).

### 3.3 Statistics (priority 3)

`qc3.lattice` implements D1 exactly:
- leave-one-out $y_c=(\sum x-x_c)/(N-1)$;
- $E^*_c=\sqrt{y_c^2-P^2}$ for each resample;
- $\bar E^*$ = mean of $E^*_c$;
- $C=\frac{N-1}{N}\sum_c(E^*_c-\bar E^*)(E^*_c-\bar E^*)^T$, symmetrised.

The χ² code is robust:
- eigenvalue pseudo-inverse that drops $\lambda\le10^{-12}\lambda_{\max}$;
- raises if $C$ is not positive semi-definite, not symmetric, or non-finite.

My linearised recomputation agrees to $10^{-5}$ relative, which is the expected size of the nonlinear $\sqrt{\cdot}$ correction. `match_levels` is D3 (it raises if $n<m$).

**The χ² proxy value 53.890** (12 levels, cutoff 0.335) is labelled "not a fit" in the code and the JSON. Two more points matter:
- It was computed with `coded_v2`, not the production form.
- It is computed only at the fixed label-generation vector $\vec c_{\rm base}$.

So it must not be quoted as a physics χ². It is a timing proxy with a consistency check (independent recomputation $1.3\times10^{-12}$, tester 011).

### 3.4 Roots (priority 4)

- **Discriminant.** $\Delta_h=\ell(E^*)-\tfrac12[\ell(E^*-h)+\ell(E^*+h)]$ with $h=10^{-6}$ and $D=4$. For a simple zero it is bounded by $\ln(x_{\rm tol}/h)=-11.5$; for a pole the sign flips.
- **A zero next to a pole.** I checked analytically that a zero within $h$ of a pole still gives $\Delta_h\approx+10.8$ (pole). Opposite-type pairs come out "pole" or "ambiguous", never the wrong kind.
- **Dimension jumps.** `_classify` compares $N$ at $E^*\pm h$. Brent converging onto a jump therefore always raises `DimensionJumpError`, because $h\gg x_{\rm tol}$.
- **Default.** `kind="zero"` is the default and the χ² proxy uses it.
- **D8.** Recorded correctly as M2 work.

### 3.5 Tests and tolerances (priority 5)

#### ISSUE C — `test_fast_single_thread_cpu_time` is an acceptance-test defect (NON-BLOCKING for physics; needs a planner fix approved by Digonto)
- **Where:** `tests/acceptance/m1/test_m1_fast_v2.py:224-233`.
- **Reproduced at HEAD:** `python3 -m pytest -q tests/acceptance/m1/test_m1_fast_v2.py::test_fast_single_thread_cpu_time` **FAILS** with CPU $0.0342$ s / wall $0.0105$ s $=3.26$ on 4 cores. With `OPENBLAS_NUM_THREADS=1` it passes, and it also passes inside the full suite.
- **Cause:** the measured window is about 10 ms, and OpenBLAS worker threads spin for 0.1–0.3 s after `import numpy`. Those threads, not the fast path, are what the test counts.
- **Why the test does not discriminate:** at $N=35$ even a NumPy tail would not thread. The test therefore neither proves its claim (no BLAS threads in the fast path) nor passes deterministically.
- **Suggested replacement (planner):**
  - run in a fresh subprocess with no thread variables;
  - `time.sleep(1)` after imports;
  - measure at 111_A2 0.35103 ($N=160$) for at least 1 s of wall time;
  - same bound of 1.3.
  
  Alternatively assert with `threadpoolctl` that no BLAS call happens inside `evaluate_point`.

#### ISSUE D — tautological assertion in `test_oracle_reality_hermiticity_isometry` (NON-BLOCKING; planner fix)
- **Where:** `tests/acceptance/m1/test_m1_oracle.py:173-175`.
- **Problem:** `r.lnabsdet` is itself `np.linalg.slogdet(Q)` (`src/python/qc3/oracle.py:65,90`), so `abs(l - r.lnabsdet) <= 1e-10` compares NumPy with itself. The shape and phase checks on the lines before it are meaningful.
- **Fix:** compare with the independent Eigen LU in `r.metadata["lnabsdet_eigen_lu"]` and `["sign_eigen_lu"]`, at a tolerance justified by $\kappa(Q)$ at that point.

#### Notes on tolerances (acceptable for M1; record for the planner)
- **N4 (tester 002).** The oracle-vs-anchor tolerances $10^{-10}$/$10^{-9}$ (`test_m1_oracle.py:120-121`) are below the round-off floor at $\kappa(Q)=3.8\times10^{12}$. They test reproducibility on the same LAPACK, not accuracy, and may fail on another CPU or BLAS. This is acceptable while the machine is fixed. Before any other machine, the planner should make them $\kappa$-scaled or compare roots instead (needs Digonto's approval: it is a tolerance change).
- **N5 (tester 010 N2).** At $\kappa(Q)=1.5\times10^{13}$ (111_A2 $E^*=0.3382011860823352$) the oracle is off by $3\times10^{-6}$ from a 60-digit value; the fast path is off by $1.2\times10^{-12}$. The $10^{-6}$ fast-vs-oracle $\ln|\det|$ tolerance is therefore unattainable by the *oracle* near dimension jumps or near roots. The locked random test passes only because its seed avoids such points (it also skips $\pm10^{-5}$ around jumps). This is not a defect today. For M2, root-position comparisons (tolerance $10^{-9}$) are the physically meaningful check.
- The unit test for classes vs direct sum at $10^{-13}$ uses 6 benign points (STATUS N1). The locked $10^{-12}$ test against the oracle covers the class path at every comparison point, so this is acceptable.

### 3.6 Provenance (priority 6) and documentation

- **N6 — xtol coupling not enforced** (`roots.py:130-150`). With `xtol` $>he^{-D}\approx1.8\times10^{-8}$ every simple zero is "ambiguous". The failure is loud, never wrong, but a `ValueError` guard plus a docstring note is cheap. Also, `find_sign_changes` aborts the whole scan on a sign change at a jump (`DimensionJumpError`); M2 should skip brackets that contain a `dimension_jumps` entry.
- **N7 — provenance:**
  - (a) `_core` was built at `5792545`; HEAD is `2fa090f`; no `src` change in between, so results are valid. Rebuild `_core` at the `milestone-1` commit before the benchmark run so that `core_built_at_commit` equals the tag.
  - (b) `results/k2_pwave_shift_L20.json` has `git_commit 97ead1d` and **`git_tree_dirty: true`**, and no `core_built_at_commit`. Its energies are golden values in the locked `test_m1_roots.py:22-32`. Tester 008 re-ran it against `47425c9`, which supports them, but it should be regenerated from a clean tree with the core commit recorded.
  - (c) The bench JSON does not record compiler flags (known, tester 011).
  - (d) Invariant 3 in `CLAUDE.md` does not list `k2_pwave_form`. The metadata does record it, but since D7 it is a production constant that changes every level. The coordinator should ask Digonto whether to add it to invariant 3.
- **N8 — C++ default form.** `oracle_point.hpp:56` defaults `k2_pwave_form="coded_v2"`. Make it `q0_consistent` or required (no default).
- **N3 — truncation claim.** METHODS §3 and `fast_point.hpp:23-25` bound each omitted term by $\epsilon_{\rm tr}$, not the omitted *sum*. In the boosted frame the number of $\mathbf n$ per $r$-shell grows by $\gamma$, so the omitted sum is roughly $\epsilon_{\rm tr}\cdot2\pi\gamma R/\alpha$ ($\sim10^{-13}$ at $\gamma\approx10$, absolute). The measured $\le7.3\times10^{-15}$ relative at $\gamma\approx10$ (tester 006) is fine, but "<$10^{-14}$ for any $\gamma$" should read "per term; total measured …".
- **Stale docs (planner/coordinator at the gate):**
  - METHODS §15 (line 217) still says "symmetric fills of $G$" (the default is ordered pairs, a coordinator decision) and does not describe the 7.7 class path as the default with `double` accumulators.
  - `fast_point.hpp:14-15` says the projector cache is "the existing get_projector_cached_v30q map" (it is now the fast-owned cache).
  - `docs/OPEN_ISSUES.md:47` ends with "production stays `coded_v2` until decision D6". This contradicts D7 two lines above, and D6 is a different decision. Issue #2 should be marked closed.

## 4. Summary for the caller

- **Blocking issues:** none.
- **Non-blocking, code (builder, after the gate or at the start of M2):**
  - ISSUE A, `fast_point.hpp:311-316,351`: add a finiteness check and a NaN-safe test for $\mathcal K_2^{-1}$;
  - ISSUE B, `fast_pairs.hpp:73-85`, `fast_point.hpp:113-125,246-248`: `check_real` the $\gamma$ of the boost and the per-momentum $F$ scalars;
  - N6: `xtol` guard;
  - N8: C++ default form.
- **Acceptance-test defects (planner fix; Digonto must approve):**
  - ISSUE C, `test_m1_fast_v2.py:224-233`: order-dependent and does not discriminate (fails alone at CPU/wall 3.26);
  - ISSUE D, `test_m1_oracle.py:173-175`: tautological.
  - N4 and N5 are tolerance notes, not defects.
- **Docs and provenance:** N3, N7, stale lines in METHODS §15, `fast_point.hpp` header, OPEN_ISSUES #2.

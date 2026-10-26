# ROOT FINDING for the three-body QC

Finding the energies $E_n$ where $\det Q(E;\vec c)=0$ (`docs/METHODS.md` §6) is the hardest and slowest part of the fit. This file covers:
- what has been tried;
- what works now;
- why the problem is hard;
- the derivative-based ("gradient change") direction Digonto identified as the likely best method.

The full lineage with file:line references is in
[`reports/research/002-cache-generation-and-root-finding-lineage.md`](../reports/research/002-cache-generation-and-root-finding-lineage.md) §5–6. Code for the historical methods is in `legacy/root_finding_lineage/`.

Notation:
- $y(E)=\mathrm{Re}\det Q(E)$ on the 20,000-point grid.
- A **bracket** is a pair of neighbouring rows with $y_iy_{i+1}<0$.
- **true zero / pole / dimension jump** are the project's labels for what caused a bracket.

## 1. Why it is hard (measured)

1. **Most sign changes are not energy levels.** The human labels cover 208 brackets (L20 and L24, dimension jumps already removed), and only 31 (15 %) are true zeros. The rest are poles of $F_3^{-1}$, which sit near the free three-particle energies, or other artefacts.

| sector | brackets | true |
|---|---:|---:|
| L20/000_A1m | 10 | 3 |
| L20/100_A2 | 25 | 2 |
| L20/110_A2 | 36 | 4 |
| L20/111_A2 | 30 | 5 |
| L20/200_A2 | 35 | 5 |
| L24/000_A1m | 20 | 3 |
| L24/100_A2 | 52 | 9 |

2. **Dimension jumps.** $N$ and $N_{\rm proj}$ jump whenever a new momentum shell enters; L20/100_A2 has 11 jumps in $[0.2631,0.36]$. Determinants of different sizes cannot be compared, so brackets across a jump are discarded. A true root sitting at a jump edge would therefore be missed.
3. **Even-order or very close root pairs** produce no sign change, so every sign-flip method misses them.
4. **Following "the eigenvalue closest to zero"** jumps between eigenvalue branches and creates fake crossings.
5. **Raw determinants underflow** to 0. Always use the LU-based $(\mathrm{sign},\log|\det|)$.

## 2. What was tried (oldest → newest)

| method | idea | fate |
|---|---|---|
| secant on $\det F_3$ (`pole_searching.h`) | Newton-like step from two guesses | needs good guesses; old S-wave code only |
| sign flip + spike/slope-ratio thresholds (`file_processor.hpp`, GPU v7–v8) | large slope ratio near a flip = pole | dropped in v9: depends on thresholds |
| inward-shape / zoom refine (v9, v16–v18) | $|y|$ trend toward the flip, with mesh refinement | mixed the left and right sides |
| side-$|\det|$ trend (v21–v24) | each side separately: decreasing $|y|$ toward the flip = zero, increasing = pole | basis of later classifiers |
| eigenbranch / SVD / eigen-min bisection (v29–v30t) | track eigenvalues; accept at a minimum of $\sigma_{\min}$ | branch switching; small $\lambda$ without singularity |
| crossing direction of $\lambda_{\rm closest}$ (v31b) | $+\to-$ = zero, $-\to+$ = pole | never validated; no outputs survive |
| digonto classifier v1 (v31zo) | flip + peak of $\min|\lambda|$ and $\sigma_{\min}$ | logic opposite to v2; abandoned |
| v2 "norm refine" (v32r–w) | add 10 CPU-rebuilt points inside each bracket; shoulder test | worked for v32w; expensive (rebuilds $F_3^{-1}$ inside FCN) |
| **`digonto_v3_window`** (v33) | 6 coarse points; a shoulder is zero-like if each step toward the flip drops $|y|$ by ≥15 % | **only automated method with a clean labelled score: L20/100_A2 TP 2, FP 0, FN 0, TN 23** |
| `digonto_v4_window` | v3 + merge zeros within $10^{-3}$ | same score on L20/100_A2 |
| `raw_sign_only` and sweep modes (v33i–l) | assign flips to lattice levels at minimum cost | not blind; ruled non-production |
| **human labels** (v33p, current gate) | Digonto labels every bracket true/false | **19 accepted L20 windows; this is what the fits use** |

**Current production rule:**
1. Human labels decide which roots exist.
2. Inside each FCN call, `digonto_v3_window` only re-finds the labelled root in a ±50-row window (expanded once to ±250).
3. The root position is linear interpolation inside the bracket.

## 3. The "gradient change" idea, stated precisely

Every working classifier above is a discrete test of the **sign of $\frac{d}{dE}\ln|\det Q|$** on the two sides of a bracket. Making that exact gives a test with no thresholds. By Jacobi's formula,

$$
g(E)\equiv\frac{d}{dE}\ln\det Q(E)=\mathrm{Tr}\big[Q^{-1}(E)\,Q'(E)\big]=\sum_k\frac{\lambda_k'(E)}{\lambda_k(E)} .
$$

- **Zero of order $m$ at $E_0$:** $g\approx+\dfrac{m}{E-E_0}$, so $g$ jumps from $-\infty$ to $+\infty$ as $E$ crosses $E_0$.
- **Pole of order $m$ at $E_p$:** $g\approx-\dfrac{m}{E-E_p}$, so $g$ jumps from $+\infty$ to $-\infty$.
- **Consequences:**
  - The direction of the jump separates zeros from poles, whatever the sign of $\det$.
  - It also detects even-order zeros, where $\det$ does not change sign but $g$ still jumps from $-$ to $+$.
  - Equivalently, zeros are local minima of $\ln|\det Q|$ and poles are local maxima.

### 3.1 How to compute it cheaply with the existing caches

In a hot window $Q_r(\vec c)=\tilde F_r+\sum_ic_i\tilde B^{(i)}_r$ (`docs/METHODS.md` §9). Two options:

1. **From cached rows (no new physics):** take central differences of $\ln|\det Q_r|$ between neighbouring rows of the same $N_{\rm proj}$ segment, $g_r\approx(\ell_{r+1}-\ell_{r-1})/(2\Delta E)$ with $\ell=\ln|\det|$, $\Delta E=4.845\times10^{-6}$. These numbers already exist in every determinant scan, so the method can be scored against the 208 labelled brackets **today, offline**.
2. **Jacobi form:** $g_r=\mathrm{Tr}[Q_r^{-1}Q_r']$, with $Q'_r$ from differences of neighbouring cached matrices, or cached derivative matrices $\tilde F'_r,\tilde B'^{(i)}_r$ added to the hot-window format.
   - **Caveat:** $V_{\rm sel}$ eigenvectors may rotate inside degenerate subspaces from row to row, so matrix differences need a fixed $V_{\rm sel}$ per segment. $\ln|\det|$ differences are immune to this, because $\det(U^\dagger QU)=\det Q$ for unitary $U$.

### 3.2 Smooth root position for Minuit (proposal)

The current root position is linear interpolation of $y$ between two rows. It is continuous in $\vec c$ only while the bracket stays the same, and a lost root gives $\chi^2=10^{100}$. Both hurt Migrad, which assumes a smooth $\chi^2$. Two options:

- **Newton step on $\ln\det$:** $E\leftarrow E-1/g(E)$, started from the labelled zero. This converges quadratically for a simple zero and stays on the same root branch as $\vec c$ changes.
- **Matrix pencil inside a bracket:** with $Q(E_r+t\Delta E)\approx Q_r+t(Q_{r+1}-Q_r)$, the roots are the real eigenvalues $t\in[0,1]$ of the generalized problem $Q_rv=-t(Q_{r+1}-Q_r)v$. These are smooth functions of $\vec c$.
  - **Caveat:** near a pole the linear model is poor and gives a spurious root, so it must be combined with the $g$-jump test.

### 3.3 Acceptance gate for any new method

A new root finder may replace human labels in a sector only if, against the labelled brackets:
- FP = 0 and FN = 0 in **every** labelled sector;
- it runs inside the fit's time budget.

**Status (M2, `prompts/003`):** the gate is locked as `tests/acceptance/m2/test_m2_rootscan.py` (scored with `coded_v2` at $\vec c_{\rm base}$, the form and parameters the labels were made with; a labelled TRUE bracket widened by ¼ row must contain exactly one `zero`, a FALSE bracket none, and no `zero` may lie outside all brackets). The planner's prototype of the §5 method passed it on all 208 brackets of the 7 labelled sectors (L20 ×5, L24 ×2; `reports/research/009` §3.1). Inside the fit the roots are not re-found by scanning but tracked (§5.3), with a fresh scan at the minimum as the check.

## 4. Rules that apply to every method

- Work on $Q=(L\xi)^{-6}V^\dagger(F_3^{-1}+\mathcal K_{\rm df,3})V$, the projected inverse. Never use $(V^\dagger F_3V)^{-1}$.
- Use LU-based $(\mathrm{sign},\ln|\det|)$, never a raw $\det$.
- Never compare across an $N_{\rm proj}$ change. Handle segment edges explicitly.
- Check Hermiticity first: $\max|\mathrm{Im}\det|$ must be 0 to rounding, otherwise the reader convention is wrong.
- **Level matching (Digonto, D3, 2026-10-21):** sort the QC levels and the lattice levels ascending in each sector and pair them in order. If there are $n$ QC levels and $m\le n$ lattice levels, use only the lowest $m$ QC levels. This is only safe when no true low-lying root is missed, because a missed root shifts every pair above it. Complete root detection is therefore part of matching (`docs/OPEN_ISSUES.md` issue 5).
- **Poles are not roots (prompts/002 §2.5, M1 step 7.1):** every sign change found by Brent is classified with $\Delta_h=\ell(E^*)-\tfrac12[\ell(E^*-h)+\ell(E^*+h)]$, $\ell=\ln|\det Q|$, $h=10^{-6}$: zero if $\Delta_h\le-4$, pole if $\ge+4$, otherwise ambiguous (raised). `qc3.roots.refine_root` returns zeros only by default. Measured on step 9's windows: zeros $-11.6\ldots-21.6$, poles $+11.9,+12.7$ (`reports/research/008` §7).
- **Exact free-level solutions are artefacts (Digonto, D8, 2026-10-22):** the QC can have a solution exactly at a non-interacting three-particle energy that the interaction does not lift (seen as an even-order zero, $|\det Q|\propto(E-E^*)^4$, no sign change: `110_A2` at $E^*=0.319468828=E_{\rm free}[K(\mathbf 0)K(0,1,0)\pi(1,0,0)]$, `reports/tests/007` finding 2, identical for both p-wave forms). Such solutions are **excluded from D3 matching**. Rule for any root finder (M2): the even-order-zero detector (the $g$-jump $-\to+$ without a sign change, §3) must compare each candidate with the sector's non-interacting spectrum (`docs/METHODS.md` §14, `tools/nonint3body_general_doublecover/`) and label candidates within $10^{-7}$ of a free level as `free_level`; only `zero` kinds go to `match_levels`. The labels (invariant 7) contain no free-level solutions and remain the ground truth.

## 5. The M2 method (plan `prompts/003`; research `reports/research/009`)

### 5.1 Label-free scan — `qc3.rootscan.scan_sector`
1. Split $[E_{\min},E_{\max}]$ at the basis-dimension jumps (`qc3.fast.dimension_jumps`), keep $2\times10^{-6}$ away from each jump.
2. Coarse grid, 4 rows ($1.94\times10^{-5}$): $(\mathrm{sign},\ell=\ln|\det Q|)$ at ≈5000 points per sector.
3. Sign change → Brent → $\Delta_h$ (§4): zero / pole / ambiguous.
4. Grid-local minimum of $\ell$ with equal signs → bounded minimisation of $\ell$ → $\Delta_h\le-4$ = even-order zero, else a smooth minimum between two poles (discarded).
5. Order of every zero-like candidate from two stencil widths, $m=(\Delta_{h_1}-\Delta_{h_2})/\ln(h_2/h_1)$ ($h_1=10^{-6},h_2=3\times10^{-6}$); kind `zero` needs odd $m$ and a sign change, `even_zero` even $m$ and none, otherwise `ambiguous`.
6. Slope $u^TQ'u$ of the vanishing eigenvalue recorded (physical zeros: $>0$, arXiv:1901.07095 Eq. (64); all 31 labelled zeros are).
7. Every free level (`qc3.free_levels`) probed with $\Delta_h$; every zero-like candidate within $10^{-7}$ of a free level → kind `free_level` (D8, §4). The prototype found 27 such solutions in the 7 sectors (orders 2 and 4 = one or two degenerate free states), none elsewhere.
8. Fine step (¼ row) within one coarse step of every candidate: resolves opposite-type pairs down to $1.2\times10^{-6}$ (4× finer than the labels). Remaining blind spot: a zero–pole pair closer than that with neither partner seen.
Cost: ≈5100–5500 evaluations and 7–12 s per L20 sector (15–23 s at L24) in the prototype, cap 10 000 in the tests (legacy: 20 000 rows). Why "gradient change" became $\Delta_h$ plus a minimum search: both are the sign of $\ell'$ on the two sides of a candidate (§3), made threshold-free by Brent's $10^{-11}$ bracket.

### 5.2 Why the hard cases of §1 are handled
- Poles (§1 item 1): $\Delta_h\ge+4$; measured $\ge+11.6$.
- Dimension jumps (§1 item 2): segments; candidates within $2\times10^{-6}$ of a jump are not produced (labelled zeros are $\ge10^{-4}$ away).
- Even-order zeros (§1 item 3): the minimum search; all found ones are free-level solutions (D8).
- Branch switching (§1 item 4): no eigenvalue tracking is used.
- Underflow (§1 item 5): everything is in $\ell$; the secant/Brent variable is $\mathrm{sign}\,e^{\ell-\ell_{\rm ref}}$.

### 5.3 Root position inside the fit — `qc3.tracking.SectorRoots`
The §3.2 proposal is realised as a predictor–corrector continuation: predictor $E_0+J\cdot\delta\vec c$ with the analytic Jacobian $J_i=-(u^T\tilde B_iu)/(u^TQ'u)$ (`docs/METHODS.md` §17), corrector = secant to $10^{-11}$, sub-steps of at most 20 rows of predicted move, re-classification with $\Delta_h$ at every accepted root, safeguards → windowed rescan (4 → 1 → ¼ rows) → full rescan, roots above $E_{\max}+0.005$ parked. Measured: ≈6 evaluations per tracked zero per χ² call, tracked vs independent roots $\le2.4\times10^{-12}$. The fit function is smooth and has an exact gradient; a sector with fewer QC zeros than lattice levels returns a finite penalty ($100\sigma$ per missing level, counted) instead of $10^{100}$. At every minimum a fresh §5.1 scan must reproduce the tracked set (`verify_fit`).

### 5.4 Rules added in M2 (to §4)
- Kinds `zero`, `pole`, `even_zero`, `free_level`, `ambiguous`; only `zero` enters `match_levels`. `even_zero` or `ambiguous` away from a free level stops the fit (loud; Digonto decides).
- Never evaluate $\det Q$ exactly at a converged root (exact zero pivot → `SingularQError`, meaning an exact root).
- Record the slope $u^TQ'u$ of every zero; a negative slope is unphysical by arXiv:1901.07095 Eq. (64) — recorded and reported loudly, never used to drop roots (D9, decided by Digonto 2026-10-24: flag only).

### 5.5 As built (M2 steps 2–6; amendment 2026-10-26, `prompts/003` §2/§4)
- **Free-level snap.** A zero-like candidate (`zero`, `even_zero`, `ambiguous`) within $10^{-7}$ of a free level becomes `free_level` **and its $E$ is overwritten by the exact free-level energy** $E_f$ (order, slope, $\Delta_h$ kept). Justified: before the snap all 20 L20 candidates lie within $9.8\times10^{-10}$ of $E_f$ ($3\times10^{-10}$ for the 110_A2 quartic), 100× inside the tolerance, and the exact solution of the truncated QC sits at $E_f$ (arXiv:1901.07095 Sec. 4.4.3; `reports/tests/015` 3a). A genuine simple zero within $10^{-7}$ of a free level would be relabelled too (rule of the plan; the nearest labelled zero is $\ge10^3$ times farther).
- **Free-level probe.** $F$ has a pole exactly at $E_f$ (the fast path raises `NONFINITE_COMPONENT`), so $\Delta_h(E_f)$ is evaluated with the centre point at $E_f+2.5\times10^{-9}$ (`_PROBE_OFFSET`) and the side points at $E_f\pm h$.
- **Rescaled variable.** The bounded Brent minimisation (item 4) runs in $t\in[0,1]$ over the bracket; in $E$, scipy's tolerance $\sqrt\epsilon|x|+x_{\rm atol}/3\approx4\times10^{-9}$ at $x=0.28$ would defeat the requested `xatol` (`reports/tests/015` 3c).
- **`ambiguous` near a pole.** A simple zero closer than $\approx4\times10^{-6}$ to a pole gets a polluted order estimate (e.g. 1.41 for a zero $3.8\times10^{-6}$ below a pole) and is classed `ambiguous`, dropped from `zeros()`, never silent. At $\vec c_{\rm base}$ and $\vec c_{\rm legacy,0.335}$ none occurs; the resolved pair ($s=8.15\times10^{-6}$) is a normal `zero`. Tracking and the fit treat any `ambiguous` or non-free-level `even_zero` as a loud event (§5.4; §5.6).
- **`qc3.roots.find_sign_changes` / `find_zeros`** take `skipped=None`. Brackets containing a dimension jump are skipped (audit N6); the list is filled if given, otherwise a `RuntimeWarning` is issued, so a skip is never silent.
- **Label gate** (`scripts/m2/score_labels.py`, `results/m2_label_gate.json`): besides FP = FN = 0 a sector fails on (i) any zero outside every scored bracket, (ii) any remaining `ambiguous` or `even_zero` candidate (every even-order zero must be a free-level solution), (iii) a `free_level` candidate farther than $10^{-7}$ from an entry of `qc3.free_levels`. Result: TP 31, FP 0, FN 0, TN 176 (207 scored + bracket 36 of L24/100_A2 excluded), 7/7 sectors, 5917–9337 evaluations and 6.5–28.6 s per sector (`reports/tests/015`).

### 5.6 Tracking and fit behaviour as built (steps 4–6)
- **Poles do not move with $\vec c$.** $Q=F_3^{-1}(E)+\sum_ic_i\tilde B_i(E)$ with $\tilde B_i$ finite, so poles are $\vec c$-independent (measured $1.5\times10^{-11}$, `reports/tests/016`). Only zeros are tracked; the pole list from the scan is reused.
- **Corrector failures fall back to a local rescan.** `AmbiguousRootError`, `PoleError` or `DimensionJumpError` raised in the corrector (reproduced with `111_A2`) triggers the windowed rescan of `prompts/003` §2.5 instead of aborting; tracking cost 5–6.7 evaluations per zero per sub-step with 0 fallbacks in the acceptance runs.
- **Zeros entering from above $E_{\max}$** are not seen by the corrector (it follows existing roots). They are caught only by `SectorRoots.rescan()` and by `verify_fit`, so every fit minimum is checked with a fresh scan.
- **Ambiguous trial points.** Inside the Levenberg–Marquardt stage only, a trial point whose tracking raises `AmbiguousRootError` is rejected (all-penalised residual, counted in `ambiguous_trials_lm`; Digonto 2026-10-26). The start, the Jacobian, Migrad, the final model and `verify_fit` still raise. The rejection acts as a wall and can steer LM to the local minimum on the side of the start (`reports/tests/018` N2).
- **$\chi^2$ is piecewise smooth.** A zero crossing a pole, or a level at a dimension jump, makes $\chi^2(\vec c)$ jump or kink; walls lie within $1\sigma$ of several minima (`reports/tests/018` N2/N3; list in `docs/OPEN_ISSUES.md` #14). Parabolic (Hesse) errors are then approximate.

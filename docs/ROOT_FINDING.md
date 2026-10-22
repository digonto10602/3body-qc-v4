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

### 3.3 Acceptance gate for any new method (proposal)

A new root finder may replace human labels in a sector only if, against the labelled brackets:
- FP = 0 and FN = 0 in **every** labelled sector;
- it runs on the hot-window or full grid inside the FCN time budget (now about 24 ms per FCN for 14 levels).

The planner turns this into locked acceptance tests (`tests/acceptance/`).

## 4. Rules that apply to every method

- Work on $Q=(L\xi)^{-6}V^\dagger(F_3^{-1}+\mathcal K_{\rm df,3})V$, the projected inverse. Never use $(V^\dagger F_3V)^{-1}$.
- Use LU-based $(\mathrm{sign},\ln|\det|)$, never a raw $\det$.
- Never compare across an $N_{\rm proj}$ change. Handle segment edges explicitly.
- Check Hermiticity first: $\max|\mathrm{Im}\det|$ must be 0 to rounding, otherwise the reader convention is wrong.
- **Level matching (Digonto, D3, 2026-10-21):** sort the QC levels and the lattice levels ascending in each sector and pair them in order. If there are $n$ QC levels and $m\le n$ lattice levels, use only the lowest $m$ QC levels. This is only safe when no true low-lying root is missed, because a missed root shifts every pair above it. Complete root detection is therefore part of matching (`docs/OPEN_ISSUES.md` issue 5).
- **Poles are not roots (prompts/002 §2.5, M1 step 7.1):** every sign change found by Brent is classified with $\Delta_h=\ell(E^*)-\tfrac12[\ell(E^*-h)+\ell(E^*+h)]$, $\ell=\ln|\det Q|$, $h=10^{-6}$: zero if $\Delta_h\le-4$, pole if $\ge+4$, otherwise ambiguous (raised). `qc3.roots.refine_root` returns zeros only by default. Measured on step 9's windows: zeros $-11.6\ldots-21.6$, poles $+11.9,+12.7$ (`reports/research/008` §7).
- **Exact free-level solutions are artefacts (Digonto, D8, 2026-10-22):** the QC can have a solution exactly at a non-interacting three-particle energy that the interaction does not lift (seen as an even-order zero, $|\det Q|\propto(E-E^*)^4$, no sign change: `110_A2` at $E^*=0.319468828=E_{\rm free}[K(\mathbf 0)K(0,1,0)\pi(1,0,0)]$, `reports/tests/007` finding 2, identical for both p-wave forms). Such solutions are **excluded from D3 matching**. Rule for any root finder (M2): the even-order-zero detector (the $g$-jump $-\to+$ without a sign change, §3) must compare each candidate with the sector's non-interacting spectrum (`docs/METHODS.md` §14, `tools/nonint3body_general_doublecover/`) and label candidates within $10^{-7}$ of a free level as `free_level`; only `zero` kinds go to `match_levels`. The labels (invariant 7) contain no free-level solutions and remain the ground truth.

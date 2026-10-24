# 009 — M2 research: label-free root finding, free-level (D8) detection, root tracking for the fit, and the minimiser

- **Date:** 2026-10-24
- **Git commit:** `c67854b` (tag `milestone-1` + graph refresh); `qc3._core` built at `c0334de`, version `0.5.1-step7.7`
- **Prompt in force:** `prompts/002` (M1, done); this report prepares `prompts/003` (M2).
- **Related:** `reports/tests/007` (hidden even-order zero, `refine_root` accepts poles), `008` (Δ_h discriminant), `011` (χ² proxy), `013` (gate), `reports/audits/001` (N4–N6), `reports/benchmarks/001`, `reports/research/002` §5 (208 labelled brackets), `docs/ROOT_FINDING.md`.
- **Prototypes (git-ignored, lost with the container):** `scratch/planner/m2/{free_levels_proto,scan_proto,tracker_proto,fit_proto}.py`, logs `scan_L20_coded_v2_4.log`, `scan_L24_coded_v2_4.log`, `tracker_coded_v2.log`, `fit_q0_lm*.log`, `fit_v32w.log`, JSON next to them. Commands: `OMP_NUM_THREADS=1 python3 scratch/planner/m2/scan_proto.py <Lbyas> <irrep> 4 coded_v2` etc. 1 core, Xeon 2.1 GHz, no GPU.

Jargon: **sector** = one $(L/a_s,\mathbf d,\Gamma)$; **candidate** = an energy where $\det Q$ may vanish or diverge; **pole** = divergence of $\det Q$ (a zero of $F_3$); **even-order zero** = $|\det Q|\propto(E-E^*)^{m}$ with $m$ even, no sign change; **free level** = non-interacting three-particle energy; **tracking/continuation** = following a root while the fit parameters $\vec c$ change; **predictor–corrector** = guess the new root from a derivative, then correct it with a root solver; **Jacobian** $J$ = $\partial E_0/\partial c_i$; **LM** = Levenberg–Marquardt least-squares minimiser; **Migrad/Hesse** = MINUIT's minimiser and its error (second-derivative) estimator.

## 1. What was asked

Plan M2: (1) a root finder that needs no human labels, scored FP = FN = 0 against all labelled brackets (invariant 7); (2) fewer QC evaluations per root inside a fit; (3) the L20 fit of $\vec c=(K^{\rm iso}_0,K^{\rm iso}_1,K_B,K_E)$ with iminuit, D1 covariance, D3 matching, D8 exclusion, production form `q0_consistent`, cutoffs 0.315–0.345, and a v32w cross-check; (4) open items from STATUS.

## 2. Candidate methods compared

### 2.1 Finding every zero of $\det Q(E)$ on $[E_{\min},E_{\max}]$ (label-free)

| method | guarantees | cost | verdict |
|---|---|---|---|
| A. uniform sign scan at the legacy grid (20 000 rows) + human labels | complete for pairs farther apart than one row; needs a human for pole/zero | 20 000 evaluations/sector (24–80 s) | baseline (what the labels are) |
| B. **segment scan + Brent + $\Delta_h$ + local-minimum search + free-level probe** (this plan) | pole/zero decided by $\Delta_h$ (`prompts/002` §2.5, margin 2.9×); even-order zeros found as local minima of $\ell=\ln\lvert\det Q\rvert$ (depth $m\ln(h/x_{\rm tol})\ge23$ vs $\lesssim10^{-4}$ for a smooth minimum); every free level probed explicitly (D8 complete by construction); blind spot: an opposite-type pair closer than the fine step | 5100–5500 evaluations/sector at a 4-row coarse step (7–12 s at L20, 15–23 s at L24); +≈30 % for the fine step around singularities | **chosen**; prototype: 0 mismatches on 208/208 brackets |
| C. eigenvalue tracking ($\lambda_k(E)$ branches, zero crossings vs. $\pm\infty$) | branch switching at avoided crossings (`docs/ROOT_FINDING.md` §1 item 4; v29–v31 history) | similar | rejected (history) |
| D. argument principle / contour integration of $\partial_E\ln\det Q$ in the complex $E$ plane (counts zeros minus poles inside a contour; e.g. Delves–Lyness, *Math. Comp.* 21 (1967) 543, DOI 10.2307/2004990) | counts $N_z-N_p$ only, needs complex-$E$ evaluations the fast path does not provide | high | rejected for M2; noted for M4 as a completeness check |

Why B's local-minimum criterion is complete at grid resolution: for a zero of even order $m$ at $E^*$ and a grid point at distance $\delta\le\Delta/2$, the two neighbours at $\delta\pm\Delta$ are higher by $m\ln(|\delta\pm\Delta|/\delta)\ge0$, so the grid point is a local minimum of $\ell$ whenever the smooth background slope over $\Delta$ is below that contrast; the background slope is $\lesssim\Delta/d_{\rm pole}\approx0.2$ for the nearest pole at $d_{\rm pole}\ge10^{-4}$. Smooth minima of $\ell$ between two poles are common (24–37 per sector); they are rejected by $\Delta_h$ after a bounded minimisation (`scipy.optimize.minimize_scalar(method="bounded")`, Brent's `fmin`, Forsythe–Malcolm–Moler 1977 / Brent 1973 ch. 5) with `xatol`$=10^{-8}$: a smooth minimum gives $\Delta_h\approx-\tfrac12\ell''h^2\approx-10^{-4}$, an even zero $\Delta_h\le-m\ln(h/x_{\rm atol})=-9.2$ for $m=2$; $D=4$ separates them with margin 2.3 in $\ell$ units.

Order estimate (new): with $\Delta_{h}=\ell(E_m)-\tfrac12[\ell(E_m-h)+\ell(E_m+h)]$ at two stencil widths $h_1<h_2$,
$$m=\frac{\Delta_{h_1}-\Delta_{h_2}}{\ln(h_2/h_1)}+O\!\big(\ell''_{\rm bg}(h_2^2-h_1^2)/\ln(h_2/h_1)\big),$$
because the singular part contributes $m\ln(|E_m-E^*|/h)$ and the smooth part cancels to second order. Prototype ($h_1=10^{-6},h_2=10^{-5}$): all 27 even zeros gave $m=2.0000\pm0.0004$ or $4.0000\pm0.0008$; the one case with a pole $3.8\times10^{-6}$ away gave 1.59 (stencil straddles the pole), so the plan uses $h_2=3\times10^{-6}$ and declares the estimate consistent only if $|m-\mathrm{round}(m)|\le0.3$.

Physical-slope check (new, from the literature): at a simple zero the vanishing eigenvalue of $F_3^{-1}+\mathcal K_{\rm df,3}$ must increase through zero, $\lambda'(E_0)>0$ (arXiv:1901.07095 Sec. 4.4.1 Eq. (64), from the positivity of the residues of the finite-volume correlator, their Eqs. (60)–(63)); zeros with $\lambda'<0$ are unphysical ("wrong residue sign"). Prototype: all 12 labelled L20 zeros have $u^{T}Q'u\in[+1.3\times10^{-4},+3.0\times10^{-2}]$; the even zeros have $|u^TQ'u|\le3\times10^{-9}$ (a double zero has zero slope). Recorded as `slope` in every candidate; not yet used to drop roots (a physics rule for Digonto, §6).

### 2.2 Exact free-level solutions (D8)

arXiv:1901.07095 Sec. 4.4.3: with $\mathcal K_2$ truncated in partial waves, the QC has "free solutions" exactly at non-interacting energies (their Eq. (49)) for finite-volume states to which the truncated interactions do not couple; in every case they found these are **double zeros** of $F_3^{-1}$ ("only touches the axis") or single zeros with the wrong residue sign, they are unphysical, and the authors' prescription is to **ignore all solutions at free energies** (also arXiv:2111.12734 Sec. 2.1: "simply ignoring all solutions at free energies"). This is exactly Digonto's D8. The prototype found 27 even-order zeros in the 7 labelled sectors, every one within $6\times10^{-9}$ of a free level (orders 2 and 4; order 4 = two degenerate free states, e.g. $K(\mathbf 0)K(0,1,0)\pi(1,0,0)$ and $K(\mathbf 0)K(1,0,0)\pi(0,1,0)$ at $\mathbf P=(1,1,0)$, multiplicity 2, $E^*=0.319468828194$), and no even zero away from a free level. The nearest labelled genuine zero to a free level is $1.09\times10^{-4}$ (L24/100_A2), $3.1\times10^{-4}$ at L20, so a tolerance of $10^{-7}$ (`docs/ROOT_FINDING.md` §4) is safe by $\ge10^{3}$. One free-level solution (200_A2, $E^*=0.346779823$) sits $3.8\times10^{-6}$ below a pole; a 4-row sign scan then sees one net sign change in the cell and Brent converged into the even zero's round-off basin (det values $\sim10^{-22}$ relative there), classifying it "zero". Hence the plan (i) probes every free level explicitly with the $\Delta_h$ stencil, (ii) requires a consistent order estimate for the kind "zero", and (iii) adds a fine step ($\tfrac14$ row) within one coarse step of every singularity.

### 2.3 Fewer evaluations per root inside a fit

| method | evaluations per root per $\vec c$ change | accuracy | verdict |
|---|---|---|---|
| a. Brent on the labelled bracket (M1 χ² proxy) | 8–9 (incl. 3 for classification) | $10^{-11}$ | baseline; needs a bracket that still contains the root |
| b. **predictor–corrector**: predict $E_0+J\cdot\delta\vec c$ with the analytic $J$ (eq. below), correct with a secant on $g=\mathrm{sign}\,e^{\ell-\ell_{\rm ref}}$, re-classify and update $J$ from $Q(E_0\pm h)$ | 4–7 (secant) + 2 (classification + $J$) | $\le2.4\times10^{-12}$ vs Brent (36 cases, 5 % and 30 % parameter changes) | **chosen** |
| c. in-memory 3-point matrix stencil $\{\tilde F(E_j),\tilde B_i(E_j)\}$ with quadratic interpolation in $E$ and μs-cost LU for any $\vec c$ ("hot window rebuilt on demand") | 0 while the root stays inside the stencil | $10^{-15}$–$4\times10^{-11}$ for half-width $h\le3\times10^{-6}$, up to $10^{-8}$ at $h=10^{-5}$ next to a pole at $9\times10^{-5}$ | optional later rung: the fit's steps move roots by $10^{-5}$–$10^{-3}$, far outside a stencil that is accurate to $10^{-9}$ |
| d. Newton on the smallest eigenvalue with eigenvector continuity | 3–4 | — | not prototyped; eigenvalue identification risk near $\lambda_{\rm next}\approx2\times10^{-7}$ (110_A2 n=2) |

The analytic Jacobian: $Q(E;\vec c)=\tilde F(E)+\sum_ic_i\tilde B_i(E)$ exactly (prototype: $1.7\times10^{-18}$ relative; $V_{\rm sel}$ is identical at all energies of one basis segment because the orbit projector depends only on the basis). With $u$ the unit null vector of the real symmetric $Q(E_0)$, the Hellmann–Feynman theorem $\partial\lambda/\partial x=u^T(\partial Q/\partial x)u$ (Feynman, Phys. Rev. 56 (1939) 340, DOI 10.1103/PhysRev.56.340; valid for a simple eigenvalue) and the implicit function theorem on $\lambda(E_0(\vec c),\vec c)=0$ give
$$\frac{\partial E_0}{\partial c_i}=-\frac{u^T\tilde B_iu}{u^TQ'(E_0)u},\qquad Q'(E_0)\approx\frac{Q(E_0+h)-Q(E_0-h)}{2h},\ h=10^{-6}.$$
Prototype at the 12 labelled roots: $J$ agrees with central finite differences of the Brent root (relative step $10^{-3}$) to $4\times10^{-8}$–$1.8\times10^{-6}$; 3 evaluations. $u$ must be taken from $\bar Q=\tfrac12[Q(E_0+h)+Q(E_0-h)]$ and $\det Q$ must never be evaluated exactly at a converged root: Eigen's LU throws on an exactly zero pivot (`fast: Q is singular`), which happened once in the prototype at a root converged to $10^{-11}$ (tester 008 saw the same in a bisection). Sensitivities: $J\in[10^{-10},4\times10^{-9}]$ per unit of $c$.

Continuation: the fit's legitimate steps are of the order of the parameter uncertainties ($\sigma_{\rm GN}\approx(1.8\times10^5,8\times10^5,8\times10^5,2.4\times10^6)$ at the baseline), which move roots by up to $J\sigma\approx10^{-3}$ (200 grid rows), more than the distance to the nearest pole ($\ge9\times10^{-5}$). The plan therefore splits a parameter change into sub-steps with predicted moves $\le20$ rows ($10^{-4}$; prediction error $\le15\%$ of the move $=3$ rows) — a tangent predictor–corrector path follower in the sense of Allgower & Georg, *Numerical Continuation Methods* (Springer 1990, DOI 10.1007/978-3-642-61257-2, ch. 2–3). Prototype fit at cutoff 0.335: median 114 evaluations per χ² call (19 tracked zeros × 6), 223 ms; 6 local rescans and 0 full rescans in 108 calls.

### 2.4 Minimiser

| method | observed | verdict |
|---|---|---|
| iminuit Migrad with analytic gradient, initial `errors` = legacy steps $(10^3,10^4,10^3,10^4)$ | first line search jumps to $30\times$ the start values; a tracked root is lost | rejected as the first-stage minimiser |
| Migrad with `errors` = Gauss–Newton σ | still visits $\vec c$ where 110_A2 has 3 zeros for 4 lattice levels → `MissingRootError` | rejected as first stage |
| **`scipy.optimize.least_squares(method="lm")`** (MINPACK `lmder`, Moré 1977, Lecture Notes in Math. 630, DOI 10.1007/BFb0067700) on the whitened residual $\rho=L^T r$ with $C^{-1}=LL^T$, analytic Jacobian $-L^TJ_{\rm model}$, `x_scale` = Gauss–Newton σ | converges in 19 residual + 12 Jacobian evaluations at 0.335 (χ² 47.9 → 29.6); trust region rejects steps that lose a level when the FCN returns a finite penalty | **chosen as stage 1** |
| **iminuit Migrad + Hesse** started at the LM minimum with `errors` = Gauss–Newton σ, `strategy=1`, `errordef=1` | EDM $1.3\times10^{-11}$, valid; Hesse errors within 15–30 % of the Gauss–Newton ones (model non-linearity) | **stage 2 and the reported errors/correlations** (Digonto's requirement: iminuit) |

iminuit facts used (https://scikit-hep.org/iminuit/reference.html, 2.33): `Minuit(fcn, x0, grad=callable)`; `errordef = Minuit.LEAST_SQUARES = 1.0`; Migrad stops at $\mathrm{EDM}<0.002\cdot\mathrm{tol}\cdot\mathrm{errordef}=2\times10^{-4}$ (default `tol=0.1`); `migrad()`/`hesse()` return the object; `fmin.edm`, `fmin.is_valid`, `fmin.has_accurate_covar`, `nfcn`, `ngrad`; `covariance`, `errors`. MINUIT algorithms: James & Roos, Comput. Phys. Commun. 10 (1975) 343, DOI 10.1016/0010-4655(75)90039-9; iminuit DOI 10.5281/zenodo.3949207.

## 3. Prototype results (numbers)

### 3.1 Label gate (coded_v2, $\vec c_{\rm base}$, $[0.26310,0.36]$, coarse step 4 rows $=1.938\times10^{-5}$)

| sector | brackets | true | mismatches | zeros | poles | even zeros (all at free levels, orders) | evaluations | s |
|---|---:|---:|---:|---:|---:|---|---:|---:|
| L20/000_A1m | 10 | 3 | 0 | 3 | 7 | 2 (2, 4) | 5126 | 7.2 |
| L20/100_A2 | 25 | 2 | 0 | 2 | 23 | 3 (4, 2, 4) | 5333 | 9.2 |
| L20/110_A2 | 36 | 4 | 0 | 4 | 32 | 4 (4, 2, 2, 4) | 5438 | 10.7 |
| L20/111_A2 | 30 | 5 | 0 | 5 | 25 | 6 (4, 4, 2, 4, 4, 2) | 5429 | 10.7 |
| L20/200_A2 | 35 | 5 | 0 | 5 (+1 at a free level, see §2.2) | 29 | 4 (2, 4, 4, 2) | 5394 | 11.9 |
| L24/000_A1m | 20 | 3 | 0 | 3 | 17 | 3 (2, 4, 2) | 5242 | 15.1 |
| L24/100_A2 | 52 | 9 | 0 | 9 | 41 | 5 (4, 2, 4, 4, 4) | 5533 | 22.5 |

Total 208 brackets, 31 true, **FP = FN = 0**, no zero outside a labelled bracket once free-level solutions are excluded. Smallest $|\Delta_h|$ over all zeros and poles: 11.60 (bound $\ln(h/x_{\rm tol})=11.51$). Labelled zeros lie inside their bracket by $1.5\times10^{-8}$ (200_A2 bracket 18) to $2.4\times10^{-6}$; fast-vs-legacy root differences are $\le4.4\times10^{-8}$ (report 011), so a scoring slack of $\tfrac14$ row ($1.2\times10^{-6}$) is safe and unambiguous (neighbouring brackets $\ge7$ rows apart). Distance of labelled zeros to the nearest dimension jump $\ge1.0\times10^{-4}$. L24/100_A2 bracket 36 (documented unresolved) is labelled `false`; the scan finds no zero in it (consistent), it is excluded from the locked score.

### 3.2 Tracking (coded_v2, 12 labelled roots, seed 20261024)

$J$ vs finite differences: relative $4.1\times10^{-8}$–$1.8\times10^{-6}$. Secant tracking after random 5 % / 30 % changes of all four parameters: root error vs an independent Brent root $\le2.4\times10^{-12}$; 4–7 evaluations (2 start points + iterations + 1 confirmation); predictor error $\le5\%$ of the move at 5 %, $\le15\%$ at 30 %. Stencil (optional method c): see §2.3.

### 3.3 Fits (prototype, 1 core; the production numbers are M2 step 5's job)

- `q0_consistent`, cutoff 0.335, 12 levels, from $\vec c_{\rm base}$: LM 19 residual evaluations → χ² = 29.615, χ²/ν = 3.70; Migrad from there: EDM $1.3\times10^{-11}$, valid; $\vec c=(6.2\times10^5,5.8\times10^5,-1.60\times10^6,4.96\times10^6)$, Hesse errors $(2.4\times10^5,4.9\times10^5,5.0\times10^5,1.6\times10^6)$, Gauss–Newton errors $(2.8,6.8,5.9,19)\times10^5$; correlation $K_B$–$K_E$ $=-0.955$; $\max|r/\sigma|=2.43$. 19 655 evaluations, 40 s for the fit; 110 s for the five initial scans; fresh-scan verification at $\vec c_{\rm fit}$: all 12 tracked roots agree with the independent scan to $1.6\times10^{-12}$ (47 s).
- v32w code-equivalence check (`coded_v2`, D1 covariance from the repo's `.jack` files, 12 levels at 0.335): χ² at $\vec c_{\rm v32w}$ = **55.04**, not 11.49; the residuals are small ($|r/\sigma|\le1.70$, $\sum(r/\sigma)^2=8.0$) so the difference is in the correlation structure of the covariance (v32w read a different data product, the `masses` files of the two-point analysis, not in the repo). The v32w number can therefore **not** be reproduced from the repo's inputs; see §6.
- `coded_v2` refit from $\vec c_{\rm v32w}$ (12 levels, 0.335): χ² 55.04 → 31.63 (χ²/ν = 3.95), $\vec c=(4.4\times10^5,1.34\times10^6,-1.87\times10^6,6.19\times10^6)$, Hesse errors $(2.2,5.9,5.0,16)\times10^5$; 17 local and 2 full rescans in 120 calls; at the minimum the tracker had lost one root in 111_A2 (0.33818) and one in 200_A2 (0.35358) that the fresh scan found — the prototype's local rescan used a single 4-row step and "nearest zero" selection, which is why the plan's tracker goes to finer local steps (4 → 1 → 1/4 rows), parks roots above $E_{\max}$, and makes `rescan()` raise on any change.
- `q0_consistent`, cutoff 0.315 (6 levels, ν = 2), from $\vec c_{\rm base}$: LM 20 evaluations → χ² 14.97 → 2.64 (χ²/ν = 1.32), but $\vec c=(3.1\times10^6,2.2\times10^7,-6.0\times10^6,-6.9\times10^6)$ with errors of 60–100 % and correlations 0.95–0.98: unconstrained. Warm-starting 0.325 from this point failed (`MissingRootError`, many full rescans at $|c|\sim10^7$). Hence the plan starts every cutoff from $\vec c_{\rm base}$ and asks Digonto about fixing $K_B,K_E$ at low cutoffs (prompts/003 §7.2).
- `q0_consistent`, cutoff 0.345 (14 levels, ν = 10), from $\vec c_{\rm base}$: LM 28 evaluations → χ² 35.54 (χ²/ν = 3.55; Migrad stage not finished at the time of writing). The 0.325 prototype run stopped on a prototype-only bug (a local-rescan window placed beyond $E_{\max}$ for a root that had moved above 0.36; the plan parks such roots).

### 3.4 χ² tolerance bookkeeping

For the 12-level covariance at 0.335: $\lambda_{\min}(C)=1.82\times10^{-9}$, $\lambda_{\max}=3.30\times10^{-6}$, $\kappa=1.8\times10^3$. If every model energy is off by at most $\delta E$, Cauchy–Schwarz in the $C^{-1}$ metric gives $|\delta\chi^2|\le2\sqrt{\chi^2}\sqrt{q}+q$ with $q=m\,\delta E^2/\lambda_{\min}$: $\le1.3\times10^{-3}$ for $\delta E=10^{-9}$ (worst cutoff 0.345, $\chi^2\le50$) and $\le1.3\times10^{-5}$ for $\delta E=10^{-11}$.

## 4. What the results mean

- The $\Delta_h$ classifier plus a bounded local-minimum search and an explicit free-level probe reproduces every human label in both volumes with no tunable threshold beyond $D=4$ (margin 2.9×). The labels can be retired from the fit path once the locked tests pass (invariant 7).
- D8 is not a heuristic: it is the published prescription for truncated QCs, and the detector sees exactly the predicted structure (double zeros, multiplicity = degeneracy).
- One χ² evaluation costs ≈100–250 ms with exact roots (≈6 evaluations per tracked zero at ≈1.2–2 ms); the old hot-window FCN's 24 ms is not reached per call, but the old pipeline needed 300–1570 s of cache building per sector first. End to end: ≈2 min of scans + ≈1 min of fit per cutoff on one core, no caches.
- The fit itself is poorly constrained ($\sigma_{\rm GN}$ comparable to the parameter values; $K_E$ consistent with zero), and the first production fit gives χ²/ν ≈ 3.7 at 0.335 — a physics result for Digonto to judge, not a code question.

## 5. Problems and assumptions

- The v32w cross-check cannot be a tight locked test (inputs differ); it is kept as a reported code-equivalence number plus a refit.
- The blind spot of any sign scan (opposite-type pairs closer than the fine step, $1.2\times10^{-6}$) is stated, not removed; method D (argument principle) would close it later.
- The L24 labels were made with the same baseline $\vec c$ and `coded_v2` (research 004 §3.5); the prototype confirms this empirically (0 mismatches).
- Timings are single-run, 1 core; noise up to 2× (report 011).

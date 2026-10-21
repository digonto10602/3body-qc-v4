# Research 008 — M1 step 7 re-plan: cold projector, basis builder, $\mathcal K_3$/$G$ kernels, tail, and pole/zero discrimination

- **Date:** 2026-10-21
- **Git commit:** `a25b1e1`
- **Prompt in force:** `prompts/001-m1-cache-free-single-energy-qc.md` (step 7 is superseded by `prompts/002`)
- **Related:** `reports/tests/006-m1-step6-fast-path-v1.md` (stage timings, notes 1–5), `reports/research/007-single-energy-qc-methods-and-reference.md`, `results/m1_step6_timings_1thread.json`, `results/k2_pwave_shift_L20.json` (step 9, poles), `docs/ROOT_FINDING.md` §3.
- **Prototype scripts (git-ignored, lost with the container):** `scratch/planner/basis_signatures.py`, `proj_orbit_proto.py`, `symmetry_check.py`, and three inline scripts (singular values, class count, pole/zero discriminant) whose outputs are quoted below.

Jargon: **cold** = first call for a basis (projector not yet in memory); **warm** = projector in memory; **orbit** = the set $\{R\mathbf k: R\in LG\}$ of a spectator momentum under the little group $LG$ of $\mathbf P$; **signature** = the ordered list of $(\text{flavour},\mathbf n,\ell,m)$ rows of a basis; **pole** = an energy where $\det Q$ diverges because $F_3^{-1}$ has a pole (near free three-particle energies), as opposed to a **zero** (an energy level).

## 1. Where the time goes now (report 006, 1 core)

| stage | $N=56$ | $N=160$ | energy-dependent? |
|---|---|---|---|
| basis builder (`config_maker_4_momentum_first`, $41^3$ candidates × 2 channels) | 21.0 ms | 21.3 ms | thresholds are not |
| $F$ (rungs 1–3) | 3.4 | 10.6 | yes |
| $\mathcal K_{\rm df,3}$ (4× `K3mat_2plus1`) | 4.9 | 42.7 | yes |
| $G$ (`G_2plus1_mat`) | 2.3 | 19.6 | yes |
| $F_3$, $F_3^{-1}$ | 0.2 | 3.8 | yes |
| NumPy tail (`project_Q`, `slogdet`) | 0.2 | 0.6 | yes |
| **warm total** | **32** | **99.5** | |
| cold projector build (once per basis) | 33–230 ms | 0.07–5.4 s | no |

The cold build is `build_cached_projector_v30q` (`src/fitter/debug_v30r_svd_validated_cartesian_l1_100_A2_zero_scan.cpp` L534): it forms dense $N\times N$ complex $S(R)W(R)$ products for every $R\in LG$, then calls `diagnose_convention`, which multiplies all pairs $U(R_a)U(R_b)$ ($|LG|^2=2304$ dense products for $O_h$) to measure closure, then an eigen-decomposition. For `000_A1m` at $N=103$ this is 5.4 s, i.e. the fast path is only $1.9\times$ faster than the oracle on a cold call.

## 2. Projector from orbits (candidate (a)) — method and evidence

**Definition in use** (arXiv:2111.12734 Eqs. 55–58; code convention `c0` of `fv_projector_cartesian_l1_v30m.hpp`): $[\hat U(R)]_{\mathbf p\ell'm'i;\mathbf k\ell mj}=\delta_{ij}\delta_{\mathbf p,R\mathbf k}\delta_{\ell'\ell}\Pi(R)\mathcal D^{(\ell)}_{m'm}(R)$ with $\Pi(R)=-1$ for improper $R$ (applied for every $\ell$, including $\ell=0$), $\mathcal D^{(0)}=1$, $\mathcal D^{(1)}$ the Cartesian signed permutation in the axis order $(m=1,0,-1)\leftrightarrow(x,z,y)$; $\hat P_\Gamma=\frac{d_\Gamma}{|LG|}\sum_R\chi_\Gamma(R)\hat U(R)$. App. A.3 Eqs. 104–106 of the same paper state the block structure: $\hat P_\Gamma$ is block diagonal in flavour, orbit and $\ell$, its eigenvalues are exactly 0 or 1, and $\mathrm{tr}\hat P_\Gamma$ is the irrep multiplicity.

**Identity used** (rearrangement lemma; valid because the $U(R)$ form a representation, which the oracle's closure diagnostic confirms): for a one-dimensional irrep, $\hat P_\Gamma\hat U(R')=\frac1{|LG|}\sum_R\chi(R)\hat U(RR')=\chi(R'^{-1})\hat P_\Gamma=\chi(R')\hat P_\Gamma$ (real characters $\pm1$). Hence for any $\mathbf k=R\mathbf k_0$ in the orbit of $\mathbf k_0$, $\hat P_\Gamma|\mathbf k,\ell,m\rangle\in\mathrm{span}_{m''}\hat P_\Gamma|\mathbf k_0,\ell,m''\rangle$, so the image of $\hat P_\Gamma$ on the (flavour, orbit, $\ell$) block is spanned by the $2\ell+1$ vectors $\hat P_\Gamma e_{(\mathbf k_0,\ell,m)}$ of **one representative** $\mathbf k_0$. Each such vector has $\le|LG|$ non-zero entries, all multiples of $1/|LG|$; a rank-revealing orthonormalisation of these $\le3$ columns gives the block of $V_{\rm sel}$. Cost $O(|LG|\,N)$ (no dense matrix, no eigen-solve). For $d_\Gamma>1$ irreps (none in production) the same construction needs all $|O|(2\ell+1)$ columns of the block; this is the general fallback.

**Evidence** (`proj_orbit_proto.py`, 20 points: 12 anchors + 2 seeded random energies per irrep, seed 20261021): $N_{\rm proj}$ identical to the fast path's eigen-decomposition at all 20 points; $\|P^2-P\|/\|P\|\le1.4\times10^{-16}$; $\|V_fV_f^\dagger-P\|/\|P\|\le6.5\times10^{-16}$ (eigen method); $\|VV^T-P\|/\|P\|\le2.6\times10^{-16}$ (orbit method); $\|V^TV-1\|\le1.2\times10^{-15}$. Singular values of the candidate columns: smallest kept $0.2887$ ($=1/\sqrt{12}$, `000_A1m` orbit of size 12), largest dropped $1.0\times10^{-16}$: a threshold of $10^{-8}$ sits $10^{7}$ above the dropped and $3\times10^{7}$ below the kept values. Orbit sizes seen: 1, 2, 3, 4, 6, 8, 12. Python prototype time 1–27 ms (pure-Python loops); in C++ this is microseconds.

**How many bases exist?** (`basis_signatures.py`) The kept set is $\{\mathbf n: E_{\rm lab}>E^+_{\min}(\mathbf n)\}$ and $E^+_{\min}$ does not depend on $E$, so the number of distinct bases over $E^*\in[0.26310,0.36000]$ is the number of distinct thresholds inside the range plus one: `000_A1m` 5, `100_A2` 12, `110_A2` 16, `111_A2` 13, `200_A2` 12 — **58 in total**. $N$ runs from 28–40 at $E^*=0.2631$ to 108–160 at $0.36$. With the orbit construction a cold build costs about as much as one $G$ row, so precomputing all 58 is unnecessary; the per-signature in-memory cache stays.

## 3. Basis builder (candidate (b))

The same observation gives the builder: precompute $t_c(\mathbf n)=E^+_{\min}(\mathbf n)$ once per (sector, channel) with the oracle's own `En_min_plus_for_config` (68 921 values × 2 channels × 5 sectors ≈ 5.5 MB), then per energy keep $\mathbf n$ in loop order iff `E_lab > t_c(n)` — the identical double comparison the oracle performs, so the basis is bit-identical by construction. NumPy timing of the per-energy selection: 0.036–0.047 ms (vs 21 ms). The sorted thresholds also give the dimension jumps directly (needed by M2's root classifier).

## 4. $\mathcal K_{\rm df,3}$ and $G$ kernels (candidates (c), (d))

`K3mat_2plus1` evaluates `K3_element_raw` for every (row, col) and the fast path calls it four times ($\vec c=e_i$). Inside, `K3B_element` and `K3E_element` recompute per element the pair kinematics (two $\sigma$, two $q^{*2}$, four square roots, four boosts, two $\gamma$ factors, the vectors $V,V'$ and the $3\times3$ tensor $t$) that depend only on the momentum pair, not on $(\ell'm';\ell m)$ — up to 16 entries per channel-1 pair. The $K^{\rm iso}$ pieces are scalars per block. A **pair kernel** computes the kinematics once per unordered momentum pair and fills the $(2\ell'+1)\times(2\ell+1)$ block of all four basis matrices, then mirrors. Measured symmetry of the oracle's basis matrices (`symmetry_check.py`, 5 points, $N$ up to 160): `iso0`, `iso1` exactly symmetric; `B` $\le1.0\times10^{-17}$; `E` $\le8.2\times10^{-17}$ relative Frobenius. Expected: 42.7 ms → 2–4 ms at $N=160$ (≈16× fewer kinematics evaluations, 2× symmetry, 4 passes → 1).

`G_ij_lm` likewise recomputes the two boosts $\mathbf k^*_p,\mathbf p^*_k$, both cutoffs and the denominator for each of the 16 entries of a channel-1 pair. Measured asymmetry of the oracle's $G$: $8.1\times10^{-16}$ to $5.4\times10^{-14}$ (`200_A2` 0.34738, the ill-conditioned point, $\kappa(F_3)=4.8\times10^8$). A symmetric fill therefore differs from the oracle by at most $5.4\times10^{-14}$ relative Frobenius (margin 18 against the $10^{-12}$ acceptance). Expected: 19.6 ms → 1–2 ms.

## 5. Tail (candidate (e))

`qc3.fast` does $Q=V^\dagger(F_3^{-1}+\mathcal K_3)V/L^6$ and `slogdet` in NumPy; OpenBLAS spawns threads unless pinned (report 006 note 1). Doing it in C++ removes the only BLAS call from the fast path. Projection of the inverse via the LU solve $F_3X=V$ ($N_{\rm proj}$ right-hand sides; Eigen `PartialPivLU::solve` accepts a matrix right-hand side, https://eigen.tuxfamily.org/dox/classEigen_1_1PartialPivLU.html) replaces the $N\times N$ inverse when $F_3^{-1}$ itself is not requested; `prompts/001` §2.1 already allows this form. $F\,X_1$ in $F_3=F/3-FX_1$ uses the block-diagonal structure of $F$ ($4\times4$ and $1\times1$ blocks): $O(16N^2)$ instead of $O(N^3)$. The matrices are converted to NumPy only when requested (`want`), which also removes ~1 ms of element-wise copies at $N=160$. The $\ln|\det|$ from an Eigen LU differs from LAPACK's by round-off only (the oracle tests at $10^{-10}$ are unaffected because the oracle keeps its NumPy tail).

## 6. $F$ (optional)

Class tables (`prompts/001` rung 5) were estimated at 3–10× on $F$. Counting classes $(n^2,\mathbf n\cdot\mathbf n_P)$ inside the truncation ellipsoid for every spectator momentum (inline script): `111_A2` 0.35103: 342 454 terms → 80 737 classes (4.2×); `100_A2` 0.3599: 4.9×; `000_A1m` 0.34538: 5.4×. Since $\mathbf n_P=\mathbf d-\mathbf n_p$ is a general integer vector, the reduction is modest; the per-class work (a $4\times4$ quadratic form in the class moments) is larger than a per-term exponential, so the realistic gain on $F$ is ≈2–3×, i.e. ≤1.3× on the total after the other rungs. It stays optional under the 1.5× stop rule. A cheaper alternative with the same accuracy, if $F$ ends up dominant: `double` accumulators with pairwise/compensated summation (error bound $O(\log_2 n\,u)$ for pairwise, Higham, SIAM J. Sci. Comput. 14 (1993) 783–799, doi:10.1137/0914050; used here only for that standard bound) instead of the x87 `long double` accumulators, and vectorisation of the inner $z$ loop.

## 7. Pole/zero discrimination (coordinator addendum, step 9 findings)

`qc3.roots.refine_root` accepts any sign change of $\det Q$. Step 9 (`results/k2_pwave_shift_L20.json`, `q0_consistent`, $\vec c_{\rm base}$) found two poles taken as roots (`110_A2` n=1 at 0.31612911112995534; `111_A2` n=0 at 0.3195543911709816) and three windows with a pole+zero pair and no end-point sign change (`100_A2` n=0, `110_A2` n=0, `200_A2` n=2). `docs/ROOT_FINDING.md` §3: zeros are local minima of $\ell(E)=\ln|\det Q(E)|$, poles local maxima ($g=\ell'$ jumps $-\infty\to+\infty$ at a zero, $+\infty\to-\infty$ at a pole).

**Discriminant tested** (inline script, fast path): $\Delta_h(E^*)=\ell(E^*)-\tfrac12[\ell(E^*-h)+\ell(E^*+h)]$. For a simple zero at $E_0$ with $|E^*-E_0|=\delta\le$ `xtol`: $\ell\approx\ln|E-E_0|+c$, so $\Delta_h\approx\ln(\delta/h)\le\ln(\text{xtol}/h)$; for a simple pole the sign flips. With xtol $=10^{-11}$, $h=10^{-6}$: $|\Delta_h|\ge11.5$ in exact arithmetic for simple zeros/poles. Measured at the 7 step-9 zeros: $\Delta_h=-11.6,-13.2,-13.9,-14.4,-20.4,-21.1,-21.6$; at the 2 poles: $+11.9,+12.7$; $\ell(E^*\pm h)$ agree to $<0.02$ at every point (the local picture is symmetric). $N$ is constant across $E^*\pm h$ at all 9 points. A decision threshold $D=4$ (factor $e^4\approx55$ in $|\det|$) has margin $\ge2.9\times$ in $\ell$ units; $|\Delta_h|<D$ is reported as ambiguous (pole and zero closer than $\sim h$; none occurs here). The builder's window-end "dip" test (DIP_MIN = 3 against the window ends) worked for these cases but depends on the window width; the local test does not.

## 8. Ranking (gain on realistic use = a scan or a fit evaluation over a window)

| rank | rung | current | expected | risk | verification |
|---|---|---|---|---|---|
| 0 | fast-owned projector cache with complete key | stale hit possible with mixed `eig_tol` | correct | none | unit + acceptance test |
| 1 | pole/zero discrimination in `qc3.roots` | poles accepted as roots | zeros only | low (ambiguous case raises) | 5 real cases + 12 coded_v2 zeros |
| 2 | orbit projector | 0.03–5.4 s per new basis | < 1 ms | low ($N_{\rm proj}$ equality + $PV=V$ self-check) | oracle subspace $VV^\dagger$ to $10^{-12}$ |
| 3 | threshold basis builder | 21 ms | < 0.1 ms | none (bit-identical rule) | basis equality, 58 bases vs reference |
| 4 | $\mathcal K_3$ pair kernel | 5–43 ms | 0.5–4 ms | transcription (caught at $10^{-12}$) | oracle + reference |
| 5 | $G$ pair kernel | 2–20 ms | 0.3–2 ms | symmetric-fill margin 18 | oracle |
| 6 | C++ tail, LU solve, lazy conversion | 0.5–2 ms + threads | 0.2–1 ms, no BLAS | none | oracle, thread check |
| 7 | $F$ (optional) | 3–11 ms | 1.5–5 ms | accuracy if the exp recurrence is used | unit truncation test, oracle |

Expected warm time after rungs 2–6: ≈5 ms at $N=56$, ≈18–22 ms at $N=160$ (speed-up vs oracle ≈500× and ≈400×); cold ≈ warm + 1 ms.

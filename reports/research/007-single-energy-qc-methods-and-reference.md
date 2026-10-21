# 007 — Single-energy QC: independent reference, speed-up candidates, p-wave $\mathcal K_2$ convention

- **Date:** 2026-10-21
- **Git commit:** `1eaf1fa`
- **Prompt in force:** `prompts/000-consolidation-and-roadmap.md` (this report feeds `prompts/001-m1-cache-free-single-energy-qc.md`)
- **Related:** `reports/research/001` (physics core), `002` (cache lineage), `003` (fitter), `004` (lattice data and inputs); `docs/OPEN_ISSUES.md` #1, #2, #4, #5.
- **Author:** planner (kickoff). Everything numerical below was run in this cloud container (4 cores, g++ 13.3, Eigen 3, NumPy 2.5.3, SciPy 1.18.1), 1 core unless stated. Scripts: `scratch/planner/*.py`, binary `scratch/planner/build/generate_cpu_component_cache`.

## 1. What was asked

Check the physics and feasibility of the starter plan, find the best numerical methods for a **cache-free, CPU-only, single-energy** evaluation of the three-body quantization condition (QC), choose an **independent reference** (QC3_release may only be cited), and settle the inputs needed for the first numbered plan: tolerances, the p-wave $\mathcal K_2^{-1}$ convention (issue #2), the covariance recipe (D1), and the minimiser library.

Jargon: **QC** = quantization condition $\det Q(E)=0$; **oracle** = the existing C++ code run on 1 core at 1 energy (Digonto's definition of the baseline); **ingredient** = one of $F$, $G$, $\mathcal K_2^{-1}$, $\mathcal K_{\rm df,3}$; **ERE** = effective-range expansion of a two-body phase shift.

## 2. Actions taken

1. Read STATUS, prompts/000, docs, research 001–006, and the production code paths with graphify (`build_component_entry`, `F2_2plus1_mat`, `G_2plus1_mat`, `K2inv_EREord2_2plus1_mat`, `K3mat_2plus1`, the v30m projector).
2. Compiled the existing CPU generator unchanged (`src/cachegen_cpu/generate_cpu_component_cache.cpp`, flags of `build.sh`) and ran it on 1 core at 12 energies over all 5 L20 sectors. Timings in §3.1; matrices dumped for the comparisons below.
3. Wrote an independent NumPy reference for $F$, $G$, $\mathcal K_2^{-1}$ from arXiv:2111.12734 (Eqs. 3–15, 20–21, App. A) and for $\mathcal K_{\rm df,3}$ by **numerical angular projection** of the Lorentz-invariant threshold expansion (Eqs. 28–31). Compared element-wise with the C++ dumps (§3.2–3.3).
4. Measured the sensitivity of $\ln|\det Q|$ to relative perturbations of the ingredients to fix tolerances (§3.4).
5. Measured the cost of the $F$ lattice sum in NumPy with and without truncation (§3.5) and derived the truncation bound.
6. Read (via ar5iv) arXiv:2111.12734 and arXiv:1411.2004 (HadSpec $\pi K$ at $m_\pi\approx391$ MeV) to resolve the p-wave convention (§3.6).
7. Installed `iminuit` 2.33.0 with pip to confirm the minimiser option (§3.7).

## 3. Results

### 3.1 The oracle: existing CPU path, 1 core, 1 energy

`generate_cpu_component_cache` (same physics path as the fitter's `build_cache_entry`: basis → $F$, $\mathcal K_2^{-1}$, $G$ → $H^{-1}F$ by LU → $F_3$ → $F_3^{-1}$ by scaled LU → $V_{\rm sel}$ → $V^\dagger F_3^{-1}V$), `OMP_NUM_THREADS=1`, production constants:

| sector | $E^*$ ($a_t^{-1}$) | $N$ | $N_{\rm proj}$ | $\mathrm{sign}$, $\ln\lvert\det Q\rvert$ at $\mathcal K_{\rm df,3}=0$ | wall per point |
|---|---|---|---|---|---|
| 100_A2 | 0.2915561 | 56 | 17 | $-1$, $-105.9913785$ | 2.8 s |
| 100_A2 | 0.3599 / 0.3600 | 108 | 30 | $-1$, $-198.58$ / $-198.81$ | 6.0 s |
| 000_A1m | 0.26932 / 0.34538 | 35 / 103 | 5 / 9 | $+1$ | 6.5 s (avg) |
| 110_A2 | 0.30733 / 0.35000 | 68 / 108 | 27 / 41 | $+1$ / $-1$ | 4.5 s |
| 111_A2 | 0.32036 / 0.35103 | 64 / 160 | 20 / 40 | $+1$ / $-1$ | 5.9 s |
| 200_A2 | 0.28664 / 0.34738 | 47 / 121 | 16 / 32 | $-1$ / $-1$ | 4.5 s |

($Q=(L\xi)^{-6}V^\dagger F_3^{-1}V$, LU `slogdet`; full values with 15 digits in `tests/acceptance/m1/data/oracle_anchors_L20.json`.) The cost is dominated by the $F$ lattice sum: for every spectator momentum and **every** $(\ell'm';\ell m)$ pair (16 per channel-1 momentum) the code loops over the full $41^3=68\,921$ cube in complex arithmetic and stores all terms for a sorted Neumaier sum. Hermiticity of $F_3^{-1}$: $\|F_3^{-1}-F_3^{-1\dagger}\|/\|F_3^{-1}\|\le 3\times10^{-19}$; all ingredients are exactly real (max $|\mathrm{Im}|=0$), as expected from the real-harmonic convention (2111.12734 App. A.1).

Condition numbers (NumPy, 2-norm): $\kappa(F_3)$ from $3.6\times10^4$ to $4.8\times10^8$ (200_A2 at 0.34738); $\kappa(H)\sim10^3$–$10^4$.

### 3.2 Independent reference for $F$, $G$, $\mathcal K_2^{-1}$ (NumPy, from the paper)

Implemented from arXiv:2111.12734 (JHEP 02 (2022) 098): Eq. (8) in its UV-regulated form of App. A.2,
$$
\big[\tilde F^{(i)}\big]_{p\ell'm';p\ell m}=\frac{H^{(i)}(\mathbf p)}{16\pi^2L^4\omega_p(E-\omega_p)}\Big[\Big(\tfrac{2\pi}{L}\Big)^{\ell'+\ell}\sum_{\mathbf n\in[-N,N]^3}\frac{\mathcal Y_{\ell'm'}(\mathbf r)\mathcal Y_{\ell m}(\mathbf r)\,e^{\alpha(x^2-r^2)}}{x^2-r^2}-\delta_{\ell'\ell}\delta_{m'm}I_\ell\Big],
$$
$\mathbf r=\mathbf n_\perp+\gamma^{-1}(\mathbf n_\parallel-\xi\,\mathbf n_P)$, $\mathbf n_P=(\mathbf P-\mathbf p)L/2\pi$, $\gamma=(E-\omega_p)/\sqrt{\sigma}$, $\xi=\tfrac12(1+(m_j^2-m_k^2)/\sigma)$, $x=q^*L/2\pi$ (complex below the pair threshold), $I_0,I_1$ as in `docs/METHODS.md` §3 with `scipy.special.erfi`; Eq. (12) for $G$ with the boosts; Eq. (15) for $\mathcal K_2^{-1}$ (both p-wave forms, §3.6); flavour factors of Eqs. (4)–(7). Basis = the oracle's momentum list, $N=20$ cube.

| object | relative Frobenius difference NumPy vs C++ (100_A2, two energies) |
|---|---|
| $F$ | $4.0\times10^{-16}$, $9.7\times10^{-16}$ (max element rel. $6\times10^{-15}$) |
| $G$ | $8\times10^{-17}$, $1.5\times10^{-16}$ |
| $\mathcal K_2^{-1}$ (coded p-wave form) | $1.2\times10^{-16}$ |
| $F_3$ / $F_3^{-1}$ (from reference ingredients) | $9\times10^{-15}$ / $7\times10^{-16}$, $1.8\times10^{-15}$ |
| $\ln\lvert\det Q\rvert$ | $-105.99137851854572$ (both), $\lvert\Delta\rvert\le1\times10^{-14}$ |

Command: `python3 -I scratch/planner/ref_proto.py`. The NumPy reference (vectorised over the cube, one pass per momentum for all $(\ell m)$ pairs) takes 0.2 s for the full $F$ at $N=56$ versus ≈2.5 s in the C++ — the first sign that the oracle's $F$ loop is the bottleneck.

### 3.3 Independent reference for $\mathcal K_{\rm df,3}$ (numerical angular projection)

Method: write $\mathcal K_{\rm df,3}$ as the Lorentz-invariant function of 2111.12734 Eqs. (28)–(29),
$$
\mathcal K_{\rm df,3}=\mathcal K^{\rm iso,0}+\mathcal K^{\rm iso,1}\Delta+\mathcal K^{B}\,(\Delta_2+\Delta_2')+\mathcal K^{E}\,\tilde t_{22},\quad
\Delta=\frac{s-M^2}{M^2},\ \Delta_2=\frac{(p_1+p_{1'})^2-4m_1^2}{M^2},\ \tilde t_{22}=\frac{(p_2-p_2')^2}{M^2},
$$
($M=2m_K+m_\pi$, $p_1,p_{1'}$ the two kaons, $p_2$ the pion, primes = final state), build the three four-momenta for a spectator $\mathbf k$ of flavour $j$ and a pair direction $\hat a^*$ (pair rest frame, $|\mathbf a^*|=q^*_k$, boosted to the lab with $\boldsymbol\beta=(\mathbf P-\mathbf k)/(E-\omega_k)$), and project numerically,
$$
[\mathcal K]^{Q_0}_{p\ell'm'i;k\ell mj}=\frac{c_{ij}}{4\pi\,(q^*_p)^{\ell'}(q^*_k)^{\ell}}\int d\Omega'\,d\Omega\;Y_{\ell'm'}(\hat a'^*)\,\mathcal K_{\rm df,3}\,Y_{\ell m}(\hat a^*),
$$
with real $Y_{\ell m}$ (App. A.1), $c_{11}=1$, $c_{12}=c_{21}=1/\sqrt2$, $c_{22}=\tfrac12$ (Eq. 30), and the $(q^*)^{\ell}$ division implementing the $\hat Q$ convention of Eqs. (20)–(21) ("delete every occurrence of $|\mathbf a^*|$"). The integrand is a polynomial of degree ≤ 2 in each direction, so a Gauss–Legendre (6) × uniform-$\phi$ (12) product rule is exact. $q^*$ is kept complex below the pair threshold (the $(q^*)^\ell$ factors cancel analytically). Convention found: $\mathbf a^*$ is the momentum of the **kaon** in the pair for channel 1 (spectator kaon); with that choice the s–p element equals $-\tfrac23\,\mathcal Y_{1m}(\mathbf p^*)/M^2$, which is what the C++ port computes.

Result on `tests/reference/data` ($E^*=0.28$, $\mathbf d=(0,0,1)$, $L/a_s=20$, $\mathcal K^{\rm iso}=(200,400)$, $\mathcal K^B=400$, $\mathcal K^E=300$, $N=51$): relative Frobenius difference $8.7\times10^{-16}$, max element difference $9\times10^{-12}$ on elements of size $8\times10^3$; per block s–s $2\times10^{-13}$, s–p $2\times10^{-12}$, p–p $9\times10^{-12}$; each of the iso, $B$, $E$ pieces is symmetric to $3\times10^{-13}$. Command: `python3 -I scratch/planner/k3_ref_proto.py` (46 s, pure Python loops). **This closes the "no stored pass result" gap of `docs/METHODS.md` §7 with a reference that is independent of QC3_release.**

### 3.4 Sensitivity of $\ln|\det Q|$ → tolerances

Random relative perturbations of all three ingredients, 10–20 draws, seed fixed (`scratch/planner/perturb_and_time.py` and the inline script in the session log):

| point | $\kappa(F_3)$ | $\varepsilon=10^{-14}$: max $\lvert\delta\ln\lvert\det Q\rvert\rvert$ | $\varepsilon=10^{-12}$ | $\varepsilon=10^{-12}$: rel. Frob. $\delta F_3^{-1}$ |
|---|---|---|---|---|
| 100_A2, 0.29156, $N=56$ | $9.0\times10^6$ | $1.1\times10^{-13}$ ($\varepsilon=10^{-15}$) | $1.0\times10^{-9}$ ($\varepsilon=10^{-11}$) | — |
| 111_A2, 0.32036, $N=64$ | $1.0\times10^7$ | $9\times10^{-13}$ | $2.0\times10^{-11}$ | $1.7\times10^{-12}$ |
| 110_A2, 0.35000, $N=108$ | $8.6\times10^4$ | $1.7\times10^{-13}$ | $1.7\times10^{-11}$ | $2.0\times10^{-11}$ |
| 111_A2, 0.35103, $N=160$ | $1.9\times10^6$ | $2\times10^{-13}$ | $1.0\times10^{-11}$ | $1.3\times10^{-12}$ |
| 200_A2, 0.34738, $N=121$ | $4.8\times10^8$ | $1.4\times10^{-9}$ | $1.1\times10^{-9}$ | $1.4\times10^{-12}$ |

Reading: $\delta\ln|\det Q|\approx 10$–$100\,\varepsilon$ at well-conditioned points, but at the worst-conditioned point there is a **floor of $\approx10^{-9}$ set by the LU round-off of $F_3$ itself** (independent of $\varepsilon$). Tolerances chosen in the plan: ingredients $10^{-12}$ relative Frobenius (margin $\ge10^3$ over the measured $10^{-15}$–$10^{-16}$), $F_3^{-1}$ $10^{-9}$ (margin 50 over $2\times10^{-11}$), $\lvert\Delta\ln|\det Q|\rvert\le10^{-6}$ with equal sign (margin $\approx700$ over the $1.4\times10^{-9}$ floor; covers $\kappa(F_3)$ up to $\sim10^{11}$), root position $\lvert\Delta E_0\rvert\le10^{-9}\,a_t^{-1}$ (a relative $10^{-10}$ change of $Q$ moves a simple root by $\lesssim10^{-13}$; lattice errors are $5\times10^{-4}$).

### 3.5 Cost of the $F$ sum and the truncation bound

NumPy, $N=56$ point, one pass per momentum (26 momenta): full cube 0.199 s; sphere $|\mathbf n|^2\le80$ (2969 points) 0.005 s; $\le120$ (5497 points) 0.009 s. But a sphere in $\mathbf n$ is **not** a safe truncation: the summand decays as $e^{-\alpha r^2}$ with $\mathbf r$ the *boosted* vector, and $r_\parallel=(n_\parallel-\xi n_P)/\gamma$ with $\gamma$ up to $\approx8$ for spectator momenta near the $H$ cutoff ($z\to0.02$). Measured sphere-vs-cube block errors were $3\times10^{-4}$ at $|\mathbf n|^2\le120$ for exactly those blocks. The correct rule is a bound on the **omitted terms of the oracle's own cube**: keep $\mathbf n$ in the cube with $r^2\le R^2$; every omitted term satisfies $|e^{\alpha(x^2-r^2)}/(x^2-r^2)|\le e^{-\alpha(R^2-\mathrm{Re}\,x^2)}/(R^2-\mathrm{Re}\,x^2)$, so with $R^2=\mathrm{Re}\,x^2+\ln(1/\epsilon)/\alpha$ and $\epsilon=10^{-16}$ ($R^2\approx x^2+74$) the truncation error relative to the kept sum is below $10^{-14}$ for any $\gamma$ (the ellipsoid in $\mathbf n$-space is simply clipped by the cube, exactly as in the oracle). $x^2\le1.3$ over the production range ($q^*\lesssim0.1$, $L=68.88$). This is the same exponential-regulator truncation argument as in 2111.12734 App. A.2 / 1901.07095 App. B.2 (sum cutoff $|\mathbf n|<n_{\max}$), applied in $\mathbf r$ rather than $\mathbf n$.

### 3.6 The p-wave $\mathcal K_2^{-1}$ convention (OPEN_ISSUES #2) — resolved on paper, decision for Digonto

- 2111.12734 Eq. (15): $[\mathcal K_2^{(i)-1}]_{\ell m}=\frac{\eta_i}{8\pi\sqrt{\sigma_i}}\{q^*\cot\delta_\ell+|q^*|[1-H^{(i)}]\}$ for every $\ell$, in the convention where $\tilde F,\tilde G$ carry $\mathcal Y_{\ell m}(\mathbf a^*)/(q^*)^\ell$ (Eq. 8). Eqs. (20)–(21): removing the $q^*$ denominators ($\hat F\to\hat Q\hat F\hat Q$, i.e. the C++ `Q0norm=true`) requires $\hat{\mathcal K}_2^{-1}\to\hat Q\hat{\mathcal K}_2^{-1}\hat Q$, so for $\ell=1$ the consistent element is $\frac{\eta}{8\pi\sqrt\sigma}\{q^{*3}\cot\delta_1+q^{*2}|q^*|(1-H)\}$. With the ERE $q^{*3}\cot\delta_1=-1/a_1+\dots$ (Eq. 17) this is $-1/a_1+q^{*2}|q^*|(1-H)$: the commented-out line 103 of `K2_functions_v2.h`.
- HadSpec, arXiv:1411.2004 (PRD 91 (2015) 054008), Sec. VI ($I=3/2$ $\pi K$, same ensemble: $a_tm_\pi=0.06906(13)$, $a_tm_K=0.09698(9)$, $\xi=3.444(6)$), ERE Eq. (9) $k^{2\ell+1}\cot\delta_\ell=1/a_\ell+\tfrac12r_\ell k^2$, global fit to the $\pi\pi\pi K$ threshold (75 points): $a_0=(-4.04\pm0.05\pm0.15)\,a_t$, $a_1=(43.2\pm3.7\pm15.4)\,a_t^3$, $a_2=(-1.13\pm0.14\pm0.58)\times10^3\,a_t^5$, $\chi^2/N_{\rm dof}=0.96$. The elastic-region fit gives $a_1=(-2.1\pm29.3\pm25.8)\,a_t^3$ (P-wave consistent with zero). This is the source of `scatter1_00 = 4.04` and `scatter1_10 = -43.2`: the code's $-1/a$ convention flips the HadSpec sign.
- **Consequence.** With $a_1^{\rm code}=-43.2$, the coded element $-1/(a_1q^{*2})+q^{*2}|q^*|(1-H)$ equals $q^*\cot\delta_1+q^{*2}|q^*|(1-H)$ with $q^{*3}\cot\delta_1=+1/43.2$ — the scattering-volume term is the HadSpec value but in the *$q$-normalised* convention, while the $(1-H)$ term and $F,G$ are in the $Q_0$ convention. The consistent $Q_0$ element is $-1/a_1^{\rm code}+q^{*2}|q^*|(1-H)=+1/43.2+q^{*2}|q^*|(1-H)$ — i.e. **the existing value $-43.2$ with the commented-out line reproduces HadSpec exactly**. The coded form overstates $\mathcal K_2^{-1}$ in the p wave by a factor $1/q^{*2}$ (≈400 at $q^*=0.05$), i.e. it suppresses the (already weak) p-wave interaction. Changing it changes production physics (invariant 3), so the plan quantifies the level shifts first and asks Digonto (D6). No uncertainty on $a_0^{KK}=4.12$ was found in HadSpec papers opened here; D2 attributes it to Digonto's own work.

### 3.7 Libraries

- **Minimiser:** `iminuit` 2.33.0 installs cleanly here (`pip install iminuit`) and wraps the same Minuit2 (MIGRAD, HESSE) as the C++ fitter — no Minuit2 build needed. Source: https://scikit-hep.org/iminuit/ . Used from M2 on.
- **C++ ↔ Python:** `pybind11` (https://pybind11.readthedocs.io, pip-installable) for the single-point evaluator; returns NumPy arrays without copies through the buffer protocol; avoids per-call process start-up of a subprocess design. Eigen for LU (partial pivoting, `PartialPivLU`), as now.

## 4. What this means (recommendations carried into prompts/001)

1. **Reference:** the planner-written `qc3ref` package (NumPy, from the papers) is the independent reference for all four ingredients; agreement with the oracle is at the $10^{-15}$ level. QC3_release is not needed.
2. **Oracle:** the existing path wrapped unchanged (`qc3.oracle`); 2.8–6.5 s per point on 1 core, anchors stored with provenance.
3. **Speed-up ladder (ranked by expected gain / certainty):** (i) one lattice-sum pass per spectator momentum for all $(\ell m)$ pairs — removes a 16× redundancy; (ii) ellipsoid-within-cube truncation with the $10^{-16}$ bound — 10–20× on $F$; (iii) real arithmetic throughout — 2–4× on everything; (iv) $V^\dagger F_3^{-1}V$ by LU solve with $N_{\rm proj}$ right-hand sides instead of the explicit inverse (still the projection of the inverse); (v) energy-independent **class moment tensors** per (sector, momentum): group $\mathbf n$ by $(n^2,\mathbf n\cdot\mathbf n_P)$, store $\sum 1$, $\sum\mathbf n$, $\sum\mathbf n\otimes\mathbf n$ once; per energy the sum is $O(\#\text{classes})$ — my own construction, in the spirit of the shell organisation of Lüscher zeta functions (Rummukainen–Gottlieb 1995; Göckeler et al. 2012); (vi) symmetric fill of $G$ and $\mathcal K_{\rm df,3}$, vectorised boosts. Expected end state: 5–20 ms per point versus 3–6 s (≥ 200×), single core, no caches, all irreps and energies. Each rung is verified against the oracle at the plan's tolerances before the next.
4. **Covariance:** per-configuration (rescaled-ensemble) files; resample leave-one-out, convert each resample $E_{\rm lab}\to E^*$, then $C=\frac{N-1}{N}\sum(y-\bar y)(y-\bar y)^T$ — identical to the trusted v32w reader; D3 matching as a separate pure function.
5. **Decision for Digonto (D6):** adopt the $Q_0$-consistent p-wave element ($-1/a_1+q^{*2}|q^*|(1-H)$, value $-43.2$ unchanged), after the plan's step reports the level shifts.

## 5. Problems and assumptions

- Sources opened and read in this session: arXiv:2111.12734 (ar5iv, both halves; Eq. (31) was only paraphrased by the page extract, the $4\pi$/$\hat Q$ normalisation was fixed by the numerical agreement of the iso term and of all p-wave blocks), arXiv:1411.2004 (ar5iv, Sec. VI numbers quoted verbatim), arXiv:1206.4141 (abstract only). Cited as background without opening: arXiv:2105.12094, arXiv:1901.07095 App. B, hep-lat/9503028, hep-lat/0507006. The method recommended does not depend on their equation numbers.
- The NumPy reference uses the oracle's basis list (momenta, order); the basis rule itself ($z>0.02$ via $E>E^+_{\min}$) is checked separately in the acceptance tests from the formula in `docs/METHODS.md` §2.
- The $10^{-9}$ floor at $\kappa(F_3)=4.8\times10^8$ means two correct implementations can differ at that level there; the tolerance $10^{-6}$ is set from this measurement, not from a desire to pass.
- Timings are from this container (unknown CPU model, shared); the benchmark protocol in the plan records `lscpu` and repeats.

# METHODS — numerical methods, equations, conventions, sources

This file lists every numerical method in use in `3body-qc-v4`. Each section gives the equations, the
code location, and the source. It is a condensed version of
[`reports/research/001-physics-core-library-audit.md`](../reports/research/001-physics-core-library-audit.md),
which has file:line references for every statement. Labels marked *(project term)* are names
used inside this code base, not standard physics vocabulary.

**Jargon used below**
- **QC**: quantization condition, the determinant equation whose zeros are the finite-volume energies.
- **RFT formalism**: the relativistic-field-theory three-particle finite-volume formalism (Hansen–Sharpe line of work).
- **2+1 system**: two identical particles plus one different particle.
- **irrep**: an irreducible representation of the little group of the finite-volume box.
- **spectator**: the particle not in the interacting pair.

## 0. System and formalism

- **System.** Two kaons + one pion, treated as a 2+1 system. The kaons are the identical pair: channel 1 has a kaon spectator and a $K\pi$ pair, with symmetry factor $\eta_1=1$; channel 2 has a pion spectator and a $KK$ pair, with $\eta_2=\tfrac12$. The data path in the configs is `KKpi_I2`, which suggests maximal isospin ($K^+K^+\pi^+$). **Open item:** confirm the isospin channel; see `docs/OPEN_ISSUES.md` issue 7.
- **References cited in the code:**
  - arXiv:2111.12734 (Blanton, Romero-López, Sharpe: implementation; $F$, $\mathcal K_2$, harmonics);
  - arXiv:2105.12094 (2+1 flavour-matrix structure of $\mathbf F_3$);
  - arXiv:1803.04169 (identical-particle test target);
  - the public Python code `QC3_release` (github `ferolo2/QC3_release`, commit `9bc4b84`). It is the reference for the $\mathcal K_{\rm df,3}$ port. It is cited only and not included in this repo (`docs/RELATED_CODES.md`).
- **Units.** Everything is in units of the temporal lattice spacing $a_t$.

| quantity | value | where set |
|---|---|---|
| anisotropy $\xi=a_s/a_t$ | 3.444 (no uncertainty in tree) | `xival` in configs |
| $a_t m_K$ | 0.09698 | `atmK` |
| $a_t m_\pi$ | 0.06906 | `atmpi` |
| box | $L = \xi\,(L/a_s)$, with $L/a_s$ = `Lbyas` = 20 (24 planned) | `PhysicsParams::L()` |
| momenta | $a_t\mathbf k = \frac{2\pi}{\xi L/a_s}\mathbf n$ | `src/core/functions.h` |
| 3-particle threshold | $a_t(2m_K+m_\pi)=0.26302$ | arithmetic |
| ensemble | `szscl21_24_128_b1p50_t_x4p300_um0p0840_sm0p0743_n1p265` (string in code) | `src/fitter/K3df_minuit_fit_v32f_*.hpp` |

## 1. Kinematics and the cutoff function $H$ — `src/core/functions.h`

$\sigma_k=(E-\omega_k)^2-(\mathbf P-\mathbf k)^2$, $\;\lambda(x,y,z)=x^2+y^2+z^2-2(xy+yz+zx)$, $\;q_k^{*2}=\lambda(\sigma_k,m_j^2,m_k^2)/(4\sigma_k)$.

$$
z=(1+\epsilon_H)\frac{\sigma-|m_j^2-m_k^2|}{(m_j+m_k)^2-|m_j^2-m_k^2|},\qquad
H=J(z)=\begin{cases}0&z\le0\\ \exp\!\left(-\frac1z e^{-1/(1-z)}\right)&0<z<1\\ 1&z\ge1\end{cases}
$$

Production uses $\epsilon_H=0$.

## 2. Basis — `config_maker_4_momentum_first` in `src/core/functions.h`

- **Rows:** spectator momentum $\mathbf k=\frac{2\pi}{L}\mathbf n$ with $\mathbf n\in[-N,N]^3$, $N$=`max_shell_num`=20. A momentum is kept iff $E>E^+_{\min}(\mathbf k)=\omega_k+\sqrt{\sigma_{\min}+(\mathbf P-\mathbf k)^2}$ with $\sigma_{\min}=|m_j^2-m_k^2|+2z_{\min}(m_j+m_k)\min(m_j,m_k)/(1+\epsilon_H)$, $z_{\min}=0.02$ (equivalent to $z>0.02$).
- **Order (locked by `tests/acceptance/m1`):** loop $n_x,n_y,n_z$ each from $-N$ to $N$ (outer → inner); for each kept momentum the $(\ell,m)$ list $(0,0),(1,-1),(1,0),(1,1)$ (channel 1) or $(0,0)$ (channel 2); all channel-1 rows first, then channel 2. The reference implementation is `tests/acceptance/m1/qc3ref/basis.py`.
- **Waves:** `waves_vec_1 = 0 1` (channel 1 has $s$ and $p$ waves); `waves_vec_2 = 0` (channel 2 has the $s$ wave only).
- **Threshold form of the rule (fast path, prompts/002 §2.2):** $E^+_{\min}(\mathbf n)$ does not depend on $E$, so the fast path computes once per (sector, channel) the table $t_c[\mathbf n]=E^+_{\min}(\mathbf n)$ for all $41^3$ candidates with the same function (`En_min_plus_for_config`, `comp` arithmetic, $z_{\min}=0.02$) and per energy keeps $\mathbf n$ iff `E_lab > t_c[n]` — the identical double comparison, so the basis is bit-identical to `config_maker_4_momentum_first`. The sorted distinct thresholds converted to $E^*=\sqrt{t^2-P^2}$ (merged within $10^{-12}$ relative) are the **dimension jumps** (`qc3.fast.dimension_jumps`); over $E^*\in[0.26310,0.36]$ there are 4/11/15/12/11 of them for `000_A1m`/`100_A2`/`110_A2`/`111_A2`/`200_A2`, i.e. 58 distinct bases in total (`reports/research/008` §2).
- **Harmonics:** real, unnormalised $\mathcal Y_{\ell m}(\mathbf p)=\sqrt{4\pi}|\mathbf p|^\ell Y_{\ell m}(\hat{\mathbf p})$, so $\mathcal Y_{1,-1}=\sqrt3p_y$, $\mathcal Y_{10}=\sqrt3p_z$, $\mathcal Y_{11}=\sqrt3p_x$ (`src/core/spherical_functions.h`).

## 3. $F$ — sum minus integral — `src/core/F2_functions_v2.h`

$$
F^{(i)}_{\mathbf p\ell'm';\mathbf k\ell m}=\delta_{\mathbf p\mathbf k}\frac{H_i(\mathbf p)}{16\pi^2L^4\omega_p(E-\omega_p)}\Big[S_{\ell'm';\ell m}-\delta_{\ell'\ell}\delta_{m'm}I_\ell\Big],\qquad
S=\Big(\tfrac{2\pi}{L}\Big)^{\ell'+\ell}\sum_{\mathbf n_a}\frac{\mathcal Y_{\ell'm'}(\mathbf r)\mathcal Y_{\ell m}(\mathbf r)e^{\alpha(x^2-r^2)}}{x^2-r^2}
$$

- **Inside the sum:** $x=q^*_pL/2\pi$ (complex below the pair threshold; every formula is evaluated with complex $x$ and the result is real), and $\mathbf r$ is the boosted lattice vector
  $$\mathbf r=\mathbf n_\perp+\gamma^{-1}\big(\mathbf n_\parallel-\xi_{jk}\,\mathbf n_P\big),\qquad \mathbf n_P=\tfrac{L}{2\pi}(\mathbf P-\mathbf p),\quad \gamma=\tfrac{E-\omega_p}{\sqrt{\sigma_p}},\quad \xi_{jk}=\tfrac12\Big(1+\tfrac{m_j^2-m_k^2}{\sigma_p}\Big),$$
  with $\mathbf n_{\parallel,\perp}$ the components along/perpendicular to $\mathbf n_P$ ($\mathbf r=\mathbf n$ when $\mathbf n_P=0$). Source: arXiv:2111.12734 App. A.2 (unequal-mass boosted sum as in Leskovec–Prelovsek / Göckeler et al. 2012).
- **Regulator:** $\alpha=0.5$ (the Gaussian damping in $S$).
- **Accuracy (oracle):** the sum is evaluated with a sorted Neumaier (compensated) sum in `long double` over the full cube.
- **Truncation rule for the fast path (prompts/001 §2.3; `src/qc3_core/fast_point.hpp` header):** keep only $\mathbf n$ in the cube with $r^2\le R^2$, where $R^2$ is the smallest value with $3R^2e^{-\alpha(R^2-\mathrm{Re}x^2)}/(R^2-\mathrm{Re}x^2)\le\epsilon_{\rm tr}$, $\epsilon_{\rm tr}=10^{-16}$ (the ellipsoid is clipped by the same $[-20,20]^3$ cube as the oracle). **The bound is per omitted term, not on the omitted total:** for $\ell\le1$, $|\mathcal Y_a\mathcal Y_b|\le3r^2$, and the right-hand side decreases for $r^2>\mathrm{Re}\,x^2+1/\alpha$, so each dropped term is $\le\epsilon_{\rm tr}$ in absolute value (re-derived in `reports/tests/006` item 5). What is known about the **total**: the number of dropped terms near the boundary grows with the boost, so the crude per-term bound times the ellipsoid shell count gives $\lesssim2\pi\gamma R\,\epsilon_{\rm tr}/\alpha\approx10^{-13}$ absolute at $\gamma\approx10$ (the production worst case, `200_A2`). The measured effect is 3–4 orders of magnitude smaller, because the real $|\mathcal Y_a\mathcal Y_b|$ is far below $3r^2$ on the dropped terms: $F(\epsilon_{\rm tr}=10^{-16})$ vs $F(10^{-40})$ differs by $\le1.1\times10^{-16}$ relative over 35 points and $\le4.8\times10^{-18}$ at the five $\gamma$-stress points ($\gamma$ up to 9.98), `reports/tests/006` rows 5a. So the unit test's $10^{-14}$ bar on the truncation is justified empirically (margin ≈90), not by the per-term bound alone; the acceptance bar $10^{-12}$ has margin $\approx10^4$. A sphere in $\mathbf n$ is **not** safe (spectator momenta near the $H$ cutoff have $\gamma$ up to ≈10; measured $3\times10^{-4}$ errors, `reports/research/007` §3.5).
- **Class moment tables (fast path default since prompts/002 sub-step 7.7; `src/qc3_core/fast_fclass.hpp`; `options["use_classes"]=False` selects the direct term-by-term sum as the cross-check):** write $\mathbf r=\mathbf n+\mathbf n_P\,\rho(\mathbf n\cdot\mathbf n_P)$ with $\rho(d)=\frac{d}{n_P^2}(\gamma^{-1}-1)-\xi_{jk}\gamma^{-1}$ (the same vector as above). Then
  $$r^2=n^2+2\rho\,(\mathbf n\cdot\mathbf n_P)+\rho^2n_P^2$$
  depends on $\mathbf n$ only through the integers $(n^2,\mathbf n\cdot\mathbf n_P)$, so the weight $w=e^{\alpha(x^2-r^2)}/(x^2-r^2)$ is constant on each **class** of cube points with equal $(n^2,\mathbf n\cdot\mathbf n_P)$. For $\ell\le1$ the harmonics are affine in $\mathbf r$ ($\mathcal Y_a=u_a+\mathbf v_a\cdot\mathbf r$), so $S_{ab}=u_au_bT_0+u_a\mathbf v_b\cdot\mathbf T_1+u_b\mathbf v_a\cdot\mathbf T_1+\mathbf v_a^{T}T_2\mathbf v_b$ with the three moments $T_0=\sum w$, $\mathbf T_1=\sum w\,\mathbf r$, $T_2=\sum w\,\mathbf r\mathbf r^T$. Per class the energy-independent data (count, mean $\bar{\mathbf n}$, centred second moments from exact integer sums, and $e^{-\alpha\,\mathrm{perp}^2}$ with $\mathrm{perp}^2=n^2-(\mathbf n\cdot\mathbf n_P)^2/n_P^2$) are tabulated once; per energy the sum over classes gives the moments with `double` accumulators (a `long double` build changes $F$ by $\le1.9\times10^{-15}$, `reports/tests/012` N3). The truncation $r^2\le R^2$ becomes $\mathrm{perp}^2\le R^2-r_\parallel^2$ per group of equal $\mathbf n\cdot\mathbf n_P$ (binary search), i.e. **the same point set** as the direct sum (term counts equal at all 12 anchors, 15 random, 5 $\gamma$-stress points and 1606 sweep energies). Tables are built only for the canonical representative $\mathbf c=(c_0\ge c_1\ge c_2\ge0)$ of the cubic orbit of $\mathbf n_P$; a signed permutation $R$ with $R\mathbf c=\mathbf n_P$ maps the cube to itself, so $T_0$ is invariant, $\mathbf T_1=R\,\mathbf T_1^{\rm can}$, $T_2=R\,T_2^{\rm can}R^T$ (back-rotation of the $\ell\le1$ moment tensors; brute-force checked over 1144 orbit elements × 8 kinematic cases, `reports/tests/012` (a)). Since $x^2,\gamma,\xi_{jk}$ depend on the spectator only through $(E,|\mathbf p|^2,n_P^2)$ and the channel masses, all spectator momenta with the same (flavour channel, $\mathbf c$, $|\mathbf p|^2$) **share one sum**, which is then rotated back (e.g. 64 blocks → 20 sums at `111_A2` 0.35103). Accuracy: $F$ vs the oracle $\le5.4\times10^{-15}$ relative, per block $\le4.4\times10^{-14}$; classes vs direct reaches $1.1\times10^{-13}$ at 4 of 1606 energies that lie within $2\times10^{-4}$ (in $r^2$ units) of a free-level pole of $F$, where the rounding of $r^2$ is amplified by $x^2/\min|x^2-r^2|\approx10^3$ (`reports/tests/012` N1). Cost: 14–51× on the $F$ stage, 3.1–5.8× on the warm total; 8 tables, 9 MB, for a full production scan.
- **Subtracted integral:** $I_\ell$ is analytic, using `Faddeeva::erfi` (MIT-licensed, `src/core/Faddeeva.cc`).
- **Normalisation:** `Q0norm=true`, so the harmonics carry $|q^*|^\ell$.
- **Flavour matrix:** $\mathbf F=\mathrm{diag}(F^{(1)}[M_1,M_1,M_2],\,F^{(2)}[M_2,M_1,M_1])$.

## 4. $G$ — one-particle exchange ("switch-state") — `src/core/G_functions_v2.h`

$$
G^{(ij)}_{\mathbf p\ell'm';\mathbf k\ell m}=\frac{\mathcal Y_{\ell'm'}(\mathbf k^*_p)H_i(\mathbf p)H_j(\mathbf k)\mathcal Y_{\ell m}(\mathbf p^*_k)}{2\omega_p^{(i)}L^3\big[(E-\omega_p-\omega_k)^2-(\mathbf P-\mathbf p-\mathbf k)^2-m_k^2\big]2\omega_k^{(j)}L^3},
\qquad
\mathbf G=\begin{pmatrix}G^{(11)}&\sqrt2(-1)^{\ell'}G^{(12)}\\ \sqrt2 G^{(21)}(-1)^{\ell}&0\end{pmatrix}.
$$

## 5. $\mathcal K_2^{-1}$ — effective-range expansion — `src/core/K2_functions_v2.h`

$$
[\mathcal K_2^{(i)}]^{-1}_{\mathbf k\ell m}=\frac{1}{2\omega_kL^3}\frac{\eta_i}{8\pi\sqrt{\sigma_k}}\big[B_\ell+C_\ell\big],\qquad
B_0=-\tfrac1{a_0}+\tfrac12r_0q^{*2},\quad C_0=|q^*|(1-H),
$$
$$
B_1=-\tfrac1{a_1q^{*2}}+\tfrac12r_1\ \ (\text{as coded}),\quad C_1=q^{*2}|q^*|(1-H).
$$

- **Flavour matrix:** $\mathbf K_2^{-1}=\mathrm{diag}(\mathcal K^{(1)-1},\,2\mathcal K^{(2)-1})$.
- **Production values** (lattice units; no uncertainties and no documented source):

| key | meaning | value |
|---|---|---|
| `scatter1_00` | $a_0$ of the $K\pi$ pair | 4.04 |
| `scatter1_10` | $a_1$ of the $K\pi$ pair | −43.2 |
| `scatter2_00` | $a_0$ of the $KK$ pair | 4.12 |

  All $r_\ell$ and $P_\ell$ are 0.
- **Two named p-wave forms (prompts/001 §2.2; `Sector.k2_pwave_form`):**
  - `coded_v2` (production default until D7; still selectable): $B_1+C_1=-\frac{1}{a_1q^{*2}}+q^{*2}|q^*|(1-H)$ — the scattering-volume term is in the $q$-normalised convention of 2111.12734 Eq. (15) while $F$, $G$ and the $(1-H)$ term are in the $Q_0$ convention of Eqs. (20)–(21): mixed.
  - `q0_consistent` (**production default from D7, Digonto 2026-10-22**; `qc3.constants.K2_PWAVE_FORM_DEFAULT`, prompts/002 sub-step 7.8): $-\frac1{a_1}+q^{*2}|q^*|(1-H)$, i.e. $\hat Q\,\mathcal K_2^{-1}\hat Q$ with $q^{*3}\cot\delta_1=-1/a_1$ (Eq. 17). Same $a_1=-43.2$. Level shifts vs `coded_v2` on the 12 L20 levels: $-0.1$ to $-0.4\,\sigma_{\rm lat}$ (`reports/tests/007`, `results/k2_pwave_shift_L20.json`). The locked M1 acceptance tests pin `coded_v2` because their anchors and labelled windows were produced with it (prompts/001 §4, amendment 2026-10-22).
- **Source of the values (research 007 §3.6):** HadSpec arXiv:1411.2004 (PRD 91, 054008) Sec. VI, $I=3/2$ $\pi K$ on this ensemble, ERE $k^{2\ell+1}\cot\delta_\ell=1/a_\ell+\tfrac12r_\ell k^2$ (their Eq. 9), global fit: $a_0=(-4.04\pm0.05\pm0.15)\,a_t$, $a_1=(43.2\pm3.7\pm15.4)\,a_t^3$; the code's $-1/a$ convention flips the sign. `q0_consistent` with $a_1=-43.2$ reproduces HadSpec exactly. $a_0^{KK}=4.12$: Digonto's earlier work (D2), no uncertainty recorded.
- **Below the pair threshold:** $|q^*|=|\sqrt{q^{*2}}|$ (real magnitude of an imaginary $q^*$), as in the reference code.

## 6. $F_3$, $F_3^{-1}$ and the QC — `src/fitter/F3_cpu_openmp_v25_K3QC_cached_core.hpp`, `src/fitter/K3df_minuit_fit_v32f_*.hpp`

$$
\mathbf H=\mathbf K_2^{-1}+\mathbf F+\mathbf G,\qquad \mathbf X=\mathbf H^{-1}\mathbf F\ (\text{partial-pivot LU}),\qquad \mathbf F_3=\tfrac13\mathbf F-\mathbf F\mathbf X .
$$

- **Inverse:** $F_3^{-1}$ is computed by scaling with $s=\max|F_{3,ab}|$ and LU-inverting $F_3/s$.
- **Quantization condition:** at each energy, per $(L,\Gamma)$ sector,

$$
Q(E;\vec c)=\frac{1}{(L\xi)^6}\,V_{\rm sel}^\dagger\big[F_3^{-1}(E)+\mathcal K_{\rm df,3}(E;\vec c)\big]V_{\rm sel},\qquad \det Q(E_n;\vec c)=0 .
$$

- **Hard rules:**
  1. Project $F_3^{-1}$, never invert a projected $F_3$.
  2. Apply the $(L\xi)^6$ scale to the matrix, not the determinant. No determinant-level $N_{\rm proj}$ rescaling.
  3. Use the corrected complex read convention (`docs/PIPELINE.md` §3).

## 7. $\mathcal K_{\rm df,3}$ — four-parameter threshold expansion — `src/core/K3_functions_2plus1.hpp`

$$
\mathcal K^{(ij)}_{\rm df,3}=c_{ij}\Big[\delta_{\ell'0}\delta_{\ell0}(K^{\rm iso}_0+K^{\rm iso}_1\Delta)+K_B\mathcal K_B^{(ij)}+K_E\mathcal K_E^{(ij)}\Big],\quad
\Delta=\frac{E^{*2}}{(2m_K+m_\pi)^2}-1,\quad c_{11}=1,\ c_{12}=c_{21}=\tfrac1{\sqrt2},\ c_{22}=\tfrac12 .
$$

- **Fit parameters:** `K3iso0, K3iso1, K3B, K3E`.
- **Pieces:** $\mathcal K_B$ and $\mathcal K_E$ are written directly in momenta; they are a port of `QC3_release/base_code/Kdf3/K3main.py`, `K3E.py`.
- **Linearity:** $\mathcal K_{\rm df,3}$ is linear in $\vec c$. The production interface returns the four basis matrices $B_{\rm iso0},B_{\rm iso1},B_B,B_E$ with $\mathcal K_{\rm df,3}(\vec c)=\sum_ic_iB_i$.
- **Lorentz-invariant definition (arXiv:2111.12734 Eqs. 28–30):** $\mathcal K_{\rm df,3}=\mathcal K^{\rm iso,0}+\mathcal K^{\rm iso,1}\Delta+\mathcal K^B(\Delta_2+\Delta_2')+\mathcal K^E\tilde t_{22}$ with $\Delta_2=[(p_1+p_{1'})^2-4m_K^2]/M^2$ (the two kaons of the initial state; primed = final) and $\tilde t_{22}=(p_2-p_2')^2/M^2$ (pion momentum transfer). Partial waves in the $Q_0$ convention: $[\mathcal K]_{p\ell'm'i;k\ell mj}=\frac{c_{ij}}{4\pi(q^*_p)^{\ell'}(q^*_k)^{\ell}}\int d\Omega'd\Omega\,Y_{\ell'm'}(\hat a'^*)\mathcal K\,Y_{\ell m}(\hat a^*)$, where $\mathbf a^*$ is the momentum of the pair **kaon** in the pair rest frame.
- **Cross-check (done, 2026-10-21):** the planner's independent reference `tests/acceptance/m1/qc3ref/k3.py` evaluates this projection numerically (Gauss–Legendre 6 × 12 azimuthal points, exact for the degree-≤2 integrand) and agrees with the C++ port to $8.7\times10^{-16}$ relative Frobenius on `tests/reference/data` and to $\le10^{-15}$ on the L20 oracle points (`reports/research/007` §3.3). The old `tests/reference/compare_K3_cpp_vs_python.py` (needs QC3_release) is superseded.

## 8. Irrep projector $V_{\rm sel}$ — `src/core/fv_projector_cartesian_l1_v30m.hpp`, `src/core/projections_v1.hpp`

$$
P_\Gamma=\frac{d_\Gamma}{|LG|}\sum_{R\in LG}\chi_\Gamma(R)\,S(R)W(R),\qquad P_\Gamma\leftarrow\tfrac12(P_\Gamma+P_\Gamma^\dagger).
$$

- **The two factors:** $S(R)$ permutes momenta. $W(R)$ applies $\Pi(R)D^\ell(R)$, with parity $-1$ for improper $R$; $D^1$ is the Cartesian signed-permutation matrix (the "v30m correction").
- **Building $V_{\rm sel}$ (oracle):** take the eigenvectors of $P_\Gamma$ with $|\lambda-1|\le$`eig_tol`$=0.05$, then re-orthonormalise them (modified Gram–Schmidt). The legacy builder also runs convention diagnostics (all $|LG|^2$ products $U(R_a)U(R_b)$ on dense $N\times N$ matrices), which is why a cold build costs 0.03–5.4 s.
- **Group action in use** (arXiv:2111.12734 Eqs. 55–58; code convention `c0` of `fv_projector_cartesian_l1_v30m.hpp`): $[\hat U(R)]_{\mathbf p\ell'm'i;\mathbf k\ell mj}=\delta_{ij}\delta_{\mathbf p,R\mathbf k}\delta_{\ell'\ell}\Pi(R)\mathcal D^{(\ell)}_{m'm}(R)$, $R\mathbf k$ = `cubic_transf(k,R)`, $\Pi(R)=-1$ for improper $R$ (applied for every $\ell$ when `parity = -1`), $\mathcal D^{(1)}$ = Cartesian signed permutation with axes $(m=1,0,-1)\leftrightarrow(x,z,y)$, $U=SW$; characters from `chi()`, little groups from `little_group()` (`src/core/projections_v1.hpp`).
- **Building $V_{\rm sel}$ from orbits (fast path, prompts/002 §2.1; `src/qc3_core/fast_projector.hpp`):** $\hat P_\Gamma$ is block diagonal in flavour, orbit $\{R\mathbf k\}$ and $\ell$, with eigenvalues exactly 0 or 1 (2111.12734 App. A.3 Eqs. 104–106). For a one-dimensional irrep the rearrangement lemma gives $\hat P_\Gamma\hat U(R)=\chi_\Gamma(R)\hat P_\Gamma$, so the image of $\hat P_\Gamma$ on a block is spanned by the $2\ell+1$ vectors $\hat P_\Gamma e_{(\mathbf k_0,\ell,m)}$ of one representative momentum $\mathbf k_0$; each has $\le|LG|$ non-zeros (multiples of $1/|LG|$) and is accumulated sparsely. A thin SVD of the $N\times(2\ell+1)$ candidate matrix keeps the left singular vectors with $s_i>10^{-8}s_{\max}$ (measured: smallest kept $0.2887$, largest dropped $10^{-16}$). Cost $O(|LG|N)$. Always-on self-check: $\|\hat P_\Gamma V-V\|_F\le10^{-10}\sqrt{N_{\rm proj}}$, $\|V^TV-1\|_F\le10^{-12}$. For $d_\Gamma>1$ the candidates are all $|O|(2\ell+1)$ columns of the block. Verified against the eigen method at 20 points: equal $N_{\rm proj}$, $\|VV^T-P\|/\|P\|\le2.6\times10^{-16}$ (`reports/research/008` §2). $V_{\rm sel}$ differs from the oracle's by a unitary mixing inside each block; $VV^\dagger$, $N_{\rm proj}$, the eigenvalues of $Q$ and $\ln|\det Q|$ are invariant and are what the tests compare. `eig_tol` is not used by this construction (kept in the metadata and the cache key).
- **Sectors:**

| label | $\mathbf d$ | little group | irrep |
|---|---|---|---|
| `000_A1m` | (0,0,0) | $O_h$ | $A_1^-$ |
| `100_A2` | (0,0,1) | $C_{4v}$ | $A_2$ |
| `110_A2` | (1,1,0) | $C_{2v}$ | $A_2$ |
| `111_A2` | (1,1,1) | $C_{3v}$ | $A_2$ |
| `200_A2` | (0,0,2) | $C_{4v}$ | $A_2$ |

- **Dimension jumps:** the number of kept columns $N_{\rm proj}$ jumps whenever a new momentum shell enters the basis. A determinant sign change across such a jump is **not** a root.

## 9. Fast fit: projected $\mathcal K_3$ basis ("hot windows") — `src/cache_build/v33p_build_projected_basis_hotwindow_cache.cpp`, `src/fitter/qc_fitter_norm_refine_v2_multiL.cpp` — **historical (D4): superseded by the cache-free path of §15**

Because $\mathcal K_{\rm df,3}=\sum_{i=1}^4c_i\mathcal K^{(i)}$:

$$
Q_r(\vec c)=\tilde F_r+\sum_{i=1}^4c_i\tilde B^{(i)}_r,\qquad
\tilde F_r=\frac{V_r^\dagger F_3^{-1}(E_r)V_r}{(L\xi)^6},\quad \tilde B^{(i)}_r=\frac{V_r^\dagger\mathcal K^{(i)}(E_r)V_r}{(L\xi)^6}.
$$

- **Offline:** the five small matrices are stored for every grid row $r$ inside a *hot window* *(project term)*, a block of rows around each human-validated zero.
- **Validation:** the basis sum agrees with direct assembly to $\le4.7\times10^{-10}$ (absolute), with 100 % determinant-sign agreement on all 5 L20 sectors.
- **Inside each fit-function call (FCN):**
  1. $\det Q_r$ by partial-pivot LU on window rows only;
  2. the classifier finds the sign change;
  3. linear interpolation gives $E_0$.
- **Cost:** about 24 ms per FCN for 14 levels, versus about 1.4 s with full-grid scans (report 003 §4.9).

## 10. Energy grid

- **Grid:** $E_i=0.26310+i\,\Delta E$, $i=0..19999$, $\Delta E=(0.36-0.26310)/19999=4.845\times10^{-6}$ in $a_t^{-1}$.
- **Exception:** the old L20/100_A2 labels used $E_0=0.26301$; `scripts/windows/build_l20_user_validated_windows.py` remaps them.

## 11. Data, covariance, χ² — `src/fitter/qc_fitter_norm_refine_v2_multiL.cpp` (`load_multil_jack_targets`), `src/fitter/lattice_data_covariance_cpp.hpp`, `chi_square_v32f`

- **Input files:** `data/lattice/L20/20_<irrep>_n<k>.jack`. One header line (`603 1 0 0 1`), then 603 lines `0 <a_t E_lab>`. They are converted per sample to $E^*=\sqrt{E_{\rm lab}^2-P^2}$ with $P=2\pi|\mathbf d|/(\xi L/a_s)$.
- **Irrep renaming in file names:** `001→100`, `002→200`, `011→110`.
- **χ²:** $\chi^2=r^TC^+r$, with $r_i=E^{\rm lat}_i-E^{\rm model}_i$. $C^+$ is an eigenvalue pseudoinverse that drops eigenvalues with $|\lambda|\le10^{-12}\lambda_{\max}$. The covariance is block-diagonal in $L$.
- **Level matching:** the $k$-th model zero in a sector is matched to the $k$-th lattice level.
- **Covariance convention (decided, D1; implemented in `qc3.lattice`, M1 step 4):** the file values $x_c$ are per-configuration (rescaled-ensemble) values. Recipe = the trusted v32w reader `lattice_data_covariance_cpp.hpp`:
  $$y_c=\tfrac1{N-1}\sum_{c'\ne c}x_{c'},\qquad y_c\to E^*_c=\sqrt{y_c^2-P^2},\qquad \bar E^*=\tfrac1N\sum_cE^*_c,\qquad C_{ij}=\tfrac{N-1}{N}\sum_c(E^{*(i)}_c-\bar E^{*(i)})(E^{*(j)}_c-\bar E^{*(j)}),$$
  which for a linear map equals $C=\frac1{N(N-1)}\sum_c(x_c-\bar x)(x_c-\bar x)^T$, i.e. $\sigma=\mathrm{s.d.}/\sqrt N$ ($\approx5\times10^{-4}$ for the L20 levels). The multi-L reader's formula (no resampling) is wrong by $(N-1)^2$ and must not be used. See `docs/OPEN_ISSUES.md` issue 1.
- **Level matching (D3; `qc3.lattice.match_levels`):** sort both lists ascending; with $m$ lattice and $n\ge m$ QC levels pair the lowest $m$; $n<m$ raises `MissingRootError`.
- **Excluded from matching (D8, Digonto 2026-10-22):** exact free-level solutions of the QC are artefacts, not levels. Example: `110_A2`, $E^*=0.319468828$, $|\det Q|\propto(E-E^*)^4$ without a sign change, exactly the non-interacting energy of $K(\mathbf 0)K(0,1,0)\pi(1,0,0)$ at $\mathbf P=(1,1,0)$ (`reports/tests/007` finding 2). Any root detector feeding `match_levels` must compare its candidates with the non-interacting spectrum (§14) and drop those within $10^{-7}$ of a free level (`docs/ROOT_FINDING.md` §4; M2, prompts/002 §9). The human labels contain no such solutions and stay the ground truth (invariant 7).

## 12. Minimiser

- **Library and calls (legacy C++ fitter):** Minuit2 `MnMigrad`, then `MnHesse` if the minimum is valid; $\mathrm{Up}=1$.
- **From M2 on:** the same Minuit2 through `iminuit` (https://scikit-hep.org/iminuit/, version 2.33 installs with pip in the cloud container; MIGRAD + HESSE, `errordef=1`), driving the cache-free fast path of §15.
- **Settings:** steps $(10^3,10^4,10^3,10^4)$, `max_fcn=400`.
- **Warm start:** each cutoff starts from the previous cutoff's best fit.
- **Penalty:** a missing root returns `failure_penalty`$=10^{100}$, which makes χ² discontinuous (see `docs/ROOT_FINDING.md`).

## 13. Root finding

Root finding has its own document: [`docs/ROOT_FINDING.md`](ROOT_FINDING.md). It covers the classifier generations, the human-label gate, and the proposed derivative-based ("gradient change") method.

- **Pole/zero discrimination (prompts/002 §2.5; `qc3.roots`):** a sign change of $\det Q$ is either a zero (energy level) or a pole of $F_3^{-1}$. With $\ell(E)=\ln|\det Q(E)|$, zeros are local minima and poles local maxima (`ROOT_FINDING.md` §3). After Brent converges to $E^*$ (bracket $\le$ `xtol` $=10^{-11}$), the discriminant
  $$\Delta_h(E^*)=\ell(E^*)-\tfrac12\big[\ell(E^*-h)+\ell(E^*+h)\big],\qquad h=10^{-6},$$
  is $\approx\ln(|E^*-E_0|/h)\le\ln(\text{xtol}/h)=-11.5$ for a simple zero and $\ge+11.5$ for a simple pole (higher orders only increase $|\Delta_h|$). Decision: zero if $\Delta_h\le-4$, pole if $\Delta_h\ge+4$, otherwise ambiguous (raised, never guessed). Measured on the step-9 cases (`results/k2_pwave_shift_L20.json`): zeros $-11.6\ldots-21.6$, poles $+11.9,+12.7$ (`reports/research/008` §7). `refine_root` returns zeros only by default (`PoleError`, `AmbiguousRootError`, `DimensionJumpError` are `MissingRootError`s); `find_sign_changes`/`find_zeros` scan a window on `n_scan` points, Brent-refine each sign-change bracket and classify (used for the D7 table and M2, not inside the fit function).

## 14. Non-interacting levels — `tools/nonint3body_general_doublecover/`

Free three-particle energies $E=\sum_a\sqrt{m_a^2+\mathbf k_a^2}$ with $\sum\mathbf k_a=\mathbf P$, projected onto irreps (double-cover groups). This is used for plot overlays and for identifying where $F_3^{-1}$ has poles.

## 15. Cache-free single-energy evaluation: oracle, fast path, tolerances — `src/python/qc3/`, `src/qc3_core/` (plan `prompts/001`)

- **Oracle** (`qc3.oracle.evaluate_point`): the existing C++ path of §§1–8 bound unchanged (basis → $F,\mathcal K_2^{-1},G$ → $X=H^{-1}F$ by LU → $F_3=F/3-FX$ → $F_3^{-1}$ by LU of $F_3/\max|F_{3,ab}|$ → $V_{\rm sel}$ → $Q$), 1 thread. Measured 2.8 s ($N=56$) to 6.5 s ($N=100$–160) per energy on 1 core. Anchor values with provenance: `tests/acceptance/m1/data/oracle_anchors_L20.json`.
- **Fast path** (`qc3.fast.evaluate_point`): same equations, same basis, with (i) one lattice-sum pass per spectator momentum for all $(\ell'm';\ell m)$, (ii) the $r^2\le R^2$ truncation of §3, (iii) real arithmetic, (iv) in-memory energy-independent pieces ($V_{\rm sel}$ per basis signature, momentum lists), $V^\dagger F_3^{-1}V$ as $V^\dagger X$ with $F_3X=V$ solved by LU (still the projection of the inverse), a symmetric (mirrored) fill of the four $\mathcal K_{\rm df,3}$ basis matrices, and (v) since sub-step 7.7 the class moment tables of §3 by default. No disk cache anywhere.
- **Fast path v2 (prompts/002, step 7):** (vi) projector from orbits (§8) in a fast-owned in-memory cache keyed by irrep, parity, $\mathbf d$, waves, `eig_tol`/`norm_tol`/`proj_tol`, method and the full basis signature (the legacy key omitted the tolerances, report 006 note 2); (vii) threshold basis builder and dimension jumps (§2); (viii) pair kernels for $\mathcal K_{\rm df,3}$ and $G$: the kinematics shared by one pair of spectator momenta (boosts, $\sigma$, $q^{*2}$, $\gamma$ factors, the $\mathcal K_E$ tensor) are computed once and the whole $(\ell'm';\ell m)$ block is filled. For $\mathcal K_{\rm df,3}$ the four basis matrices are filled for each **unordered** pair and mirrored (they are symmetric to $\le1.3\times10^{-16}$ in the oracle; the fast matrices are exactly symmetric and agree with the oracle to $\le8.9\times10^{-17}$). For $G$ the **production default fills every ordered pair** in the oracle's operation order, so $G$ is **bit-identical** to `G_2plus1_mat` (exactly 0.0 difference at 47 energies and 4 off-production parameter sets); the oracle's $G$ is itself symmetric only to $\le1.5\times10^{-13}$ relative (the $1/(\ldots)$ denominators are evaluated in a different order for the two index orders), and at $\kappa(Q)\approx4\times10^{12}$ that asymmetry moves $\ln|\det Q|$ by $1.2\times10^{-8}$, so bit-faithfulness was preferred to the plan's symmetric fill (coordinator decision, `reports/tests/010` N3; cost $+0.14$ ms at $N=160$). The symmetric fill (pairs $A\le B$ plus the transposed block, exactly symmetric, $\le1.03\times10^{-13}$ from the oracle) remains available as `g_symmetric=True` on `qc3._core.fast_point` only; `qc3.fast.evaluate_point` does not expose it; (ix) C++ tail: $F_3=F/3-FX_1$ with the block-diagonal structure of $F$, $Q=(V^TX+V^T\mathcal K_3V)/L^6$ with $X$ the LU solution of $F_3X=V$ ($N_{\rm proj}$ right-hand sides; still the projection of the inverse), $(\mathrm{sign},\ln|\det|)$ from Eigen's partial-pivot LU, the full $F_3^{-1}$ only on request; no BLAS call in the fast path, so no thread spawning (process CPU/wall $0.97$–$1.00$ unpinned in steady state at $N=160$ and $N=103$, `reports/tests/010`; the acceptance test `test_fast_single_thread_cpu_time` measures this at those two points over $\ge0.6$ s of calls after an idle settle, amendment acceptance-v5); (x) sub-step 7.7: the $F$ lattice sum from the class moment tables of §3 (per canonical $\mathbf n_P$ orbit, sums shared across spectator momenta with equal (channel, canonical $\mathbf n_P$, $|\mathbf p|^2$), back-rotation of the $\ell\le1$ moment tensors, `double` accumulators, identical truncation set; `reports/tests/012`). The orbit/character projector is described in §8 and the threshold basis builder in §2. Measured warm cost after 7.7: 0.49–3.9 ms per point on one core (3.4 ms at $N=160$); cold = warm + ≤25 ms (threshold tables + projector); χ² proxy (12 roots, 102 evaluations) 0.14 s.
- **Agreement required (prompts/001 §4, from the sensitivity study in `reports/research/007` §3.4):** ingredients $\le10^{-12}$ relative Frobenius; $F_3^{-1}$ $\le10^{-9}$; $\lvert\Delta\ln|\det Q|\rvert\le10^{-6}$ with equal sign (the LU round-off floor is $1.4\times10^{-9}$ at $\kappa(F_3)=4.8\times10^8$); root positions $\le10^{-9}\,a_t^{-1}$. Measured sensitivity: $\delta\ln|\det Q|\approx10$–$100\,\varepsilon$ for a relative ingredient perturbation $\varepsilon$ at well-conditioned points.
- **Root refinement for comparisons** (`qc3.roots.refine_root`): Brent's method on $f(E)=\mathrm{sign}\cdot e^{\ln|\det Q|-\ell_0}$ inside a bracket with a sign change ($\ell_0$ = the larger endpoint value, to avoid under/overflow). The classifier that decides *which* brackets are roots is M2 (`docs/ROOT_FINDING.md`).
- **Benchmark metric:** median single-core wall time per point (`OMP_NUM_THREADS=1`, pinned), 5 irreps × 3 energies + the 12 anchors, 5 repeats after a warm-up; the χ² proxy is the time to refine the 12 labelled L20 roots at cutoff 0.335 (`prompts/001` §7).

## 16. Independent reference implementation — `tests/acceptance/m1/qc3ref/` (planner-owned)

NumPy/SciPy, written from arXiv:2111.12734 (Eqs. 3–15, 20–21, 28–31, App. A), not from the C++: `basis.py` (basis rule and order), `fgk.py` ($F$ with the full cube and `scipy.special.erfi`, $G$, both $\mathcal K_2^{-1}$ forms), `k3.py` ($\mathcal K_{\rm df,3}$ by numerical angular projection of the Lorentz invariants). Agreement with the oracle at the 12 L20 anchor points: basis exact; $F\le2.5\times10^{-14}$, $G\le2.4\times10^{-16}$, $\mathcal K_2^{-1}\le7\times10^{-15}$, $\mathcal K_{\rm df,3}\le10^{-15}$ (relative Frobenius). Cost: ≈0.2–1 s per point in Python. QC3_release is not used.

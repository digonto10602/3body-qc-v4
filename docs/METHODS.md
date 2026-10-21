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
  - the public Python code `QC3_release` (github `ferolo2/QC3_release`, commit `9bc4b84`). It is the reference for the $\mathcal K_{\rm df,3}$ port; see `third_party/QC3_release/`.
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

- **Rows:** spectator momentum $\mathbf k=\frac{2\pi}{L}\mathbf n$ with $\mathbf n\in[-N,N]^3$, $N$=`max_shell_num`=20. A momentum is kept iff $z>0.02$.
- **Order:** momentum outermost, then $(\ell,m)$.
- **Waves:** `waves_vec_1 = 0 1` (channel 1 has $s$ and $p$ waves); `waves_vec_2 = 0` (channel 2 has the $s$ wave only).
- **Harmonics:** real, unnormalised $\mathcal Y_{\ell m}(\mathbf p)=\sqrt{4\pi}|\mathbf p|^\ell Y_{\ell m}(\hat{\mathbf p})$, so $\mathcal Y_{1,-1}=\sqrt3p_y$, $\mathcal Y_{10}=\sqrt3p_z$, $\mathcal Y_{11}=\sqrt3p_x$ (`src/core/spherical_functions.h`).

## 3. $F$ — sum minus integral — `src/core/F2_functions_v2.h`

$$
F^{(i)}_{\mathbf p\ell'm';\mathbf k\ell m}=\delta_{\mathbf p\mathbf k}\frac{H_i(\mathbf p)}{16\pi^2L^4\omega_p(E-\omega_p)}\Big[S_{\ell'm';\ell m}-\delta_{\ell'\ell}\delta_{m'm}I_\ell\Big],\qquad
S=\Big(\tfrac{2\pi}{L}\Big)^{\ell'+\ell}\sum_{\mathbf n_a}\frac{\mathcal Y_{\ell'm'}(\mathbf r)\mathcal Y_{\ell m}(\mathbf r)e^{\alpha(x^2-r^2)}}{x^2-r^2}
$$

- **Inside the sum:** $x=q^*_pL/2\pi$, and $\mathbf r$ is the boosted lattice vector.
- **Regulator:** $\alpha=0.5$ (the Gaussian damping in $S$).
- **Accuracy:** the sum is evaluated with a sorted Neumaier (compensated) sum in `long double`.
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
- **⚠ Open issue:** the coded $B_1$ differs from the reference code, which uses $-1/a_1+q^{*2}|q^*|(1-H)$. See `docs/OPEN_ISSUES.md` issue 2.

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
- **Linearity:** $\mathcal K_{\rm df,3}$ is linear in $\vec c$. This is what makes the fast fit possible (§9).
- **Cross-check:** C++ vs Python is tested by `tests/reference/compare_K3_cpp_vs_python.py`. No stored pass result was found.

## 8. Irrep projector $V_{\rm sel}$ — `src/core/fv_projector_cartesian_l1_v30m.hpp`, `src/core/projections_v1.hpp`

$$
P_\Gamma=\frac{d_\Gamma}{|LG|}\sum_{R\in LG}\chi_\Gamma(R)\,S(R)W(R),\qquad P_\Gamma\leftarrow\tfrac12(P_\Gamma+P_\Gamma^\dagger).
$$

- **The two factors:** $S(R)$ permutes momenta. $W(R)$ applies $\Pi(R)D^\ell(R)$, with parity $-1$ for improper $R$; $D^1$ is the Cartesian signed-permutation matrix (the "v30m correction").
- **Building $V_{\rm sel}$:** take the eigenvectors of $P_\Gamma$ with $|\lambda-1|\le$`eig_tol`$=0.05$, then re-orthonormalise them (modified Gram–Schmidt).
- **Sectors:**

| label | $\mathbf d$ | little group | irrep |
|---|---|---|---|
| `000_A1m` | (0,0,0) | $O_h$ | $A_1^-$ |
| `100_A2` | (0,0,1) | $C_{4v}$ | $A_2$ |
| `110_A2` | (1,1,0) | $C_{2v}$ | $A_2$ |
| `111_A2` | (1,1,1) | $C_{3v}$ | $A_2$ |
| `200_A2` | (0,0,2) | $C_{4v}$ | $A_2$ |

- **Dimension jumps:** the number of kept columns $N_{\rm proj}$ jumps whenever a new momentum shell enters the basis. A determinant sign change across such a jump is **not** a root.

## 9. Fast fit: projected $\mathcal K_3$ basis ("hot windows") — `src/cache_build/v33p_build_projected_basis_hotwindow_cache.cpp`, `src/fitter/qc_fitter_norm_refine_v2_multiL.cpp`

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
- **⚠ Covariance convention:** the multi-L reader treats the 603 values as jackknife samples, $C=\frac{N-1}{N}\sum(\cdot)(\cdot)$. The legacy reader `lattice_data_covariance_cpp.hpp` treats them as per-configuration values and resamples first. The two differ by $(N-1)^2$. See `docs/OPEN_ISSUES.md` issue 1.

## 12. Minimiser

- **Library and calls:** Minuit2 `MnMigrad`, then `MnHesse` if the minimum is valid; $\mathrm{Up}=1$.
- **Settings:** steps $(10^3,10^4,10^3,10^4)$, `max_fcn=400`.
- **Warm start:** each cutoff starts from the previous cutoff's best fit.
- **Penalty:** a missing root returns `failure_penalty`$=10^{100}$, which makes χ² discontinuous (see `docs/ROOT_FINDING.md`).

## 13. Root finding

Root finding has its own document: [`docs/ROOT_FINDING.md`](ROOT_FINDING.md). It covers the classifier generations, the human-label gate, and the proposed derivative-based ("gradient change") method.

## 14. Non-interacting levels — `tools/nonint3body_general_doublecover/`

Free three-particle energies $E=\sum_a\sqrt{m_a^2+\mathbf k_a^2}$ with $\sum\mathbf k_a=\mathbf P$, projected onto irreps (double-cover groups). This is used for plot overlays and for identifying where $F_3^{-1}$ has poles.

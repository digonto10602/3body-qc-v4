# A — Physics core of the KKπ three-body quantization-condition code

> **Research report — consolidation investigation.** Date: 2026-10-21. Git commit: none (written before the repo existed; first commit of `3body-qc-v4`). Prompt in force: `prompts/000-consolidation-and-roadmap.md`. Related: `reports/research/001`–`006`.
> Absolute paths below refer to the unpacked original zips (`/home/digimonk/Projects/three-body-qc-codes/extracted/...`) and to a temporary scratchpad. For each file's new location in this repo, see `docs/PROVENANCE.tsv`.
> The system is two kaons + one pion treated as a 2+1 system (identical kaons). Where this report says "K K-bar pi", read "KKπ"; that wording came from the investigation brief and is not established. See `docs/OPEN_ISSUES.md` #7.


- **Date:** 2026-10-21
- **Git commit:** none — the extracted folders are not a git repository (read-only inspection of `/home/digimonk/Projects/three-body-qc-codes/extracted/`).
- **Prompt in force:** consolidation task, area "A: core physics library" (no `prompts/NNN` file exists for this task).
- **Related reports:** sibling reports B, C, … written by other agents in the same `reports/` folder (GPU path, fitters, cache formats). This report does not repeat their content.
- **Hash list used:** `scratchpad/hashes.txt` (md5 of every source file). Note: that list does **not** contain `.hh` or `.cuh` files; I hashed those myself with `md5sum` (read-only).

Path shorthands used everywhere below (all absolute):

| Shorthand | Absolute path |
|---|---|
| `L20/` | `/home/digimonk/Projects/three-body-qc-codes/extracted/KKpi-L20-user-validated-hotwindow-fitter-main/KKpi-L20-user-validated-hotwindow-fitter-main/source/` |
| `CG/` | `/home/digimonk/Projects/three-body-qc-codes/extracted/full_cache_generator_package_20260629-main/full_cache_generator_package_20260629-main/` |
| `V3/` | `/home/digimonk/Projects/three-body-qc-codes/extracted/3body_quantization_v3-master/3body_quantization_v3-master/` |
| `V33O/` | `/home/digimonk/Projects/three-body-qc-codes/extracted/KKpi_v33o_production_fitter-main/KKpi_v33o_production_fitter-main/source/` |
| `REF/` | `/home/digimonk/Projects/three-body-qc-codes/extracted/_nested_old_versions/QC3_release/QC3_release/` (the published Python reference code, README line 3: "Accompanying code to arXiv:2111.12734", by T. Blanton, F. Romero-López, S. Sharpe) |

---

## 1. What was asked

1. Which header files form the canonical, newest physics core used by the L20 hot-window fitter, the v33 fitters and the CPU/GPU cache generator; best copy of each; hash-identical duplicates.
2. The equations actually implemented for $F_2$, $G$, $\mathcal K_2^{-1}$, $F_3$, $F_3^{-1}$, $\mathcal K_{\rm df,3}$, units, basis, shell cutoff, irreps and the projector `Vsel`.
3. Which files are older experiments / dead ends and what they tried.
4. Bugs and fixes noted in comments (plus any I found by reading).
5. A recommended minimal file set for a consolidated repository.

## 2. Actions taken

- Grouped every core header by md5 hash across all packages (table in §3.2).
- Traced the `#include` chain from the three production entry points:
  - L20 fitter: `L20/qc_fitter_norm_refine_v2_multiL.cpp` (built by `.../KKpi-L20-.../build.sh` lines 15-17),
  - CPU cache generator: `CG/cpu_component_cachegen/source/generate_cpu_component_cache.cpp` (line 24),
  - GPU cache generator: `CG/gpu_component_cachegen/source/F3_gpu_cachegen_v32zu.cu` (lines 27-43, 473).
- Read line-by-line: `functions.h`, `spherical_functions.h`, `F2_functions_v2.h` (core part, lines 621-1424), `K2_functions_v2.h`, `G_functions_v2.h`, `QC_functions_v2.h`, `K3_functions_2plus1.hpp`, `fv_projector_cartesian_l1_v30m.hpp`, the projector parts of `projections_v1.hpp`, the F3 assembly in `F3_cpu_openmp_v25_K3QC_cached_core.hpp` and `K3df_minuit_fit_v32f_fullF3inv_QCfull_cached_classifier.hpp`, the projector cache in `debug_v30r_svd_validated_cartesian_l1_100_A2_zero_scan.cpp`.
- Diffed every old `V3/source/*_v2.h` file against the newest copy.
- Cross-checked conventions against the bundled reference Python code `REF/base_code/` (F, G, K2, K3, cutoff function, harmonics).
- Skimmed the old V3 drivers, Python scripts and the version-history Markdown files in `CG/legacy_source/`.
- No builds, no runs, no edits to the extracted folders.

---

## 3. Results

### 3.1 Formalism and references (as stated in the code)

Jargon first:
- **RFT** = "relativistic field theory" three-particle finite-volume formalism (the Hansen–Sharpe line of work). The code never uses the word "RFT"; it cites the papers below.
- **QC** = quantization condition: the determinant equation whose zeros are the finite-volume energy levels.
- **2+1 system** = two identical-type particles plus one different particle. Here: two kaons + one pion.

References explicitly cited in the code:

| Citation in code | Where | Used for |
|---|---|---|
| arXiv:2111.12734 | `L20/F2_functions_v2.h:622`, `L20/K2_functions_v2.h:10` ("tilde_K2 … equations 2.15 and 2.14"), `L20/spherical_functions.h:10`, `V3/source/QC_functions.h:81` | $F$, $\mathcal K_2$, harmonics $\mathcal Y_{\ell m}$; same paper as the Python reference `REF/README.md:3` |
| arXiv:2105.12094 | `V3/source/QC_functions.h:159` ("F3 for 2+1 system, how it is defined in the paper") | the 2+1 flavour-matrix structure of $\mathbf F_3$ |
| arXiv:1803.04169 (Fig. 1) | `V3/source/QC_functions.h:82` | test target for identical-particle F3 ("Raul's paper", lines 83, 111) |
| "Blanton–Romero-López–Sharpe Sec. 2.4 / App. A.3, Eq. (2.54), (2.56), App. A.1" | `L20/fv_projector_cartesian_l1_v30m.hpp:5,140,172,276` | irrep projector and real $\ell=1$ harmonics |
| "QC3_release/base_code/Kdf3/K3main.py", "K3E.py" | `L20/F3_cpu_openmp_v25_K3QC_cached_core.hpp:115-117`, `L20/K3_functions_2plus1.hpp:135` | $\mathcal K_{\rm df,3}$ parametrisation; checked by `CG/legacy_source/compare_K3_cpp_vs_python.py` |

The arXiv numbers above are copied from the code. My own (unverified here) identification: 2105.12094 is the Blanton–Sharpe paper on the 2+1 ($\pi^+\pi^+K^+$-type) formalism and 2111.12734 is the Blanton–Romero-López–Sharpe implementation paper. I did not open the papers, so equation numbers are only the ones written in the code comments.

**Channel assignment** (`L20/K2_functions_v2.h:15`, `V3/source/QC_functions.h:159-162`):
- flavour/channel $i=1$: spectator = kaon, interacting pair = $K\pi$, symmetry factor $\eta_1 = 1$;
- channel $i=2$: spectator = pion, interacting pair = $KK$, symmetry factor $\eta_2 = 1/2$.

So the code treats $K\bar K\pi$ exactly like a $\pi^+\pi^+K^+$-type 2+1 system, with the two kaons playing the role of the identical pair (masses $M_1 = m_K$, $M_2 = m_\pi$). Whether this identical-pair treatment is appropriate for $K\bar K$ (a $G$-parity / isospin argument) is not discussed anywhere in the code; that is a physics assumption to state in the new repo's `docs/METHODS.md`.

### 3.2 Question 1 — the canonical, newest physics core

**Include chain of the newest CPU code** (L20 fitter, v33o/v33p fitter, CPU cache generator all share it):

```
qc_fitter_norm_refine_v2_multiL.cpp (L20)         generate_cpu_component_cache.cpp (CG, line 24)
  └─ qc_fitter_norm_refine_v2.cpp:1                  │
      └─ K3df_minuit_fit_v32f_fullF3inv_QCfull_cached_classifier.hpp  ◄┘
          └─ (line 61, with main renamed) debug_v30r_svd_validated_cartesian_l1_100_A2_zero_scan.cpp
              ├─ F3_cpu_openmp_v25_K3QC_cached_core.hpp  (lines 58-64)
              │    ├─ functions.h ── spherical_functions.h
              │    ├─ F2_functions_v2.h ── Faddeeva.hh, spherical_functions.h
              │    ├─ K2_functions_v2.h
              │    ├─ G_functions_v2.h
              │    ├─ QC_functions_v2.h
              │    ├─ projections_v1.hpp ── dig_tools.hpp, real_wigner_d.hpp
              │    └─ K3_functions_2plus1.hpp
              └─ fv_projector_cartesian_l1_v30m.hpp
link: Faddeeva.cc
```

The GPU cache generator (`CG/gpu_component_cachegen/source/F3_gpu_cachegen_v32zu.cu`) uses GPU ports of the same formulas (`F2_gpu_safe_builder.cuh`, `G_gpu_safe_builder.cuh`, `K2_functions_gpu_safe.cuh`, `projections_gpu_v32y_exact_device.cuh`), plus the CPU `functions.h` (via v16 file line 82), `K3_functions_2plus1.hpp` and `QC_functions_v2.h` (lines 42-43), and the same `fv_projector_cartesian_l1_v30m.hpp` (line 473). The GPU details belong to another report; §3.9 lists only the physics-relevant differences.

**Best copy and duplicates** (md5 prefix, line count). "All 7 newest" means the identical file sits in: `CG/legacy_source/`, `L20/`, `V33O/`, `kkpi-v33n-.../source/`, `KKpi-v33p-.../source/`, and both `KKpi-v32w-v32y-repro-master/.../v32w_*/source/` and `.../v32y_*/source/`.

| File | Best copy | md5 (first 8) | Lines | Identical copies | Older different copy |
|---|---|---|---|---|---|
| `functions.h` | `L20/functions.h` | 381dadb5 | 1588 | all 7 newest + `CG/gpu_component_cachegen/source/` | `V3/source/functions.h` (e7608dd0, 1127 lines; lacks `config_maker_4_momentum_first`, `print_config_5col`, rewritten `boost`) |
| `spherical_functions.h` | `L20/spherical_functions.h` | 45f1c181 | 118 | all 9 copies identical, incl. V3 and CG GPU | — |
| `F2_functions_v2.h` | `L20/F2_functions_v2.h` | d79a1b32 | 2376 | all 7 newest | `V3/source/F2_functions_v2.h` (c5d6846f, 726 lines) — see §3.3.1 |
| `K2_functions_v2.h` | `L20/K2_functions_v2.h` | c08fe1b9 | 326 | all 7 newest | `V3/source/K2_functions_v2.h` (048dd016, 324) |
| `G_functions_v2.h` | `L20/G_functions_v2.h` | 0c9532c6 | 478 | all 7 newest | `V3/source/G_functions_v2.h` (e21d11bf, 316) |
| `QC_functions_v2.h` | `L20/QC_functions_v2.h` | 0c131940 | 322 | all 7 newest | `V3/source/QC_functions_v2.h` (5bc56fd8, 101) |
| `K3_functions_2plus1.hpp` | `L20/K3_functions_2plus1.hpp` | 276f568e | 559 | all 7 newest | not in V3 |
| `projections_v1.hpp` | `L20/projections_v1.hpp` | 34ba293b | 2413 | all 7 newest | not in V3 |
| `real_wigner_d.hpp` | `L20/real_wigner_d.hpp` | 12f3da71 | — | all 7 newest + CG GPU | — |
| `dig_tools.hpp` | `L20/dig_tools.hpp` | fc67ff86 | 88 | all 7 newest | `V3/source/dig_tools.hpp` (149f2ca5, 56; lacks `print_function_output_to_file`) |
| `fv_projector_cartesian_l1_v30m.hpp` | `L20/fv_projector_cartesian_l1_v30m.hpp` | 6910d692 | 415 | all 7 newest + `CG/k3iso01_swave_cpu_cache_qc_omp_package/source/` | — |
| `Faddeeva.cc` | `L20/Faddeeva.cc` | d89c365b | 2516 | all 8 copies identical | — |
| `Faddeeva.hh` (C++ API, used) | `L20/Faddeeva.hh` | 1cde27c1 | — | all 8 identical | — |
| `Faddeeva.h` (C99 API, **not used**) | — | 6e77829e | 67 | all 8 identical | — |
| `F3_cpu_openmp_v25_K3QC_cached_core.hpp` | `L20/…` | e8abc13a | 1405 | L20, V33O, `CG/legacy_source/` | — |
| `debug_v30r_svd_validated_cartesian_l1_100_A2_zero_scan.cpp` | `L20/…` | bb8df5d2 | 2203 | L20, V33O, `CG/legacy_source/` | — |
| `K3df_minuit_fit_v32f_…_classifier.hpp` | `L20/…` (= V33O) | e74d0696 | 1595 | L20, V33O | `CG/legacy_source/` copy (9d6f534b) differs only in Minuit2 include guards and the precomputed $\mathcal K_3$ basis (non-physics) |

Old S-wave-only family, carried unchanged into every package but **not** in the newest include chain: `F2_functions.h` (a11e8322, 718), `G_functions.h` (4a1b54bb, 152), `K2_functions.h` (44386df3, 180), `QC_functions.h` (566c4dc2, 2817), `pole_searching.h` (2ce552df), `splines.h` (a978c475), `projections.hpp` (1c649130), `projections_from_config.hpp` (dd3f7ddd), `eigenvalue_tracker.hpp` (467d4baf). All identical everywhere, including V3.

**Answer:** the newest physics core is the `_v2` family plus `functions.h`, `spherical_functions.h`, `K3_functions_2plus1.hpp`, `projections_v1.hpp` (+`real_wigner_d.hpp`, `dig_tools.hpp`), `fv_projector_cartesian_l1_v30m.hpp`, and Faddeeva. Every copy in the newest packages is byte-identical, so `L20/` can be taken as the single source.

### 3.3 Question 2 — equations implemented

Notation (mine, chosen to match the code variable names): $E$ = lab-frame total energy (`En`), $\mathbf P$ = total momentum (`total_P`), $E^* = \sqrt{E^2-\mathbf P^2}$ (`Ecm`), spectator flavour $i$ with mass $m_i$ and momentum $\mathbf k$, pair masses $m_j, m_k$, $\omega_k = \sqrt{k^2+m_i^2}$.

Kinematics (`L20/functions.h`):
- $\sigma_k = (E-\omega_k)^2 - (\mathbf P - \mathbf k)^2$ — `sigma_pvec_based`, lines 37-60.
- $\lambda(x,y,z) = x^2+y^2+z^2-2(xy+yz+zx)$ — `kallentriangle`, lines 62-69.
- $q^{*2}_k = \lambda(\sigma_k, m_j^2, m_k^2)/(4\sigma_k)$ — `q2psq_star`, lines 71-78.
- Boost of a 4-vector $(p^0,\mathbf p)$ into the rest frame of $(E_2,\mathbf P_2)$, with $\boldsymbol\beta = \mathbf P_2/E_2$: $\mathbf p^* = \mathbf p + \left[(\gamma-1)\,\mathbf p\cdot\hat{\boldsymbol\beta} - \gamma\beta\,p^0\right]\hat{\boldsymbol\beta}$ — `boost`, lines 1281-1328 (returns $\mathbf p$ unchanged when $\mathbf P_2=0$).

#### 3.3.1 Cutoff function $H$ (`L20/functions.h:97-141`)

$$
z = (1+\epsilon_H)\,\frac{\sigma - |m_j^2-m_k^2|}{(m_j+m_k)^2 - |m_j^2-m_k^2|},\qquad
J(z)=\begin{cases}0 & z\le 0\\ \exp\!\left(-\tfrac{1}{z}e^{-1/(1-z)}\right) & 0<z<1\\ 1 & z\ge 1\end{cases},\qquad H = J(z).
$$

For $m_j=m_k$ and $\epsilon_H=0$ it reduces to $z=\sigma/(4m^2)$ (line 126-130). This matches the "new" cutoff in the reference `REF/base_code/defns.py:72-77` (the denominator equals $2\min(m_j,m_k)(m_j+m_k)$). One difference: the reference `jj` (`defns.py:61-69`) clips $J=0$ for $z<0.02$ and $J=1$ for $z\ge0.97$; the C++ `Jfunc` does not clip. Because the C++ basis already drops every spectator with $z\le 0.02$ (next item), the remaining difference is $|1-J(0.97)|\approx 10^{-15}$, i.e. negligible. Production value: $\epsilon_H = 0$.

#### 3.3.2 Basis: spectator momentum × $(\ell,m)$, shell cutoff (`L20/functions.h:976-1113`)

`config_maker_4_momentum_first` (production) loops over integer vectors $\mathbf n\in[-N,N]^3$ with $N$ = `max_shell_num` = 20, sets $\mathbf k = \frac{2\pi}{L}\mathbf n$, and **keeps** $\mathbf n$ when

$$E > E_{\min}^{+}(\mathbf k) = \omega_k + \sqrt{\sigma_{\min} + (\mathbf P-\mathbf k)^2},\qquad \sigma_{\min} = |m_j^2-m_k^2| + \frac{2\,z_{\min}(m_j+m_k)\min(m_j,m_k)}{1+\epsilon_H},\quad z_{\min}=0.02$$

(lines 1053-1056; formula in `En_min_plus_for_config`, lines 732-758). This is exactly "$z>0.02$", i.e. $H > J(0.02)\approx 1.5\times10^{-8}$. The value of $H$ and the `tolerance` argument are computed but **not used** for the selection (lines 1049-1051). For each kept momentum all $(\ell,m)$ from `waves_vec` are appended, momentum outermost (comment lines 991-1000). Each basis row stores $(k_x,k_y,k_z,\ell,m)$ plus the integer $(n_x,n_y,n_z,\ell,m)$.

Production waves: `waves_vec_1 = {0,1}` (channel 1 has $s$ and $p$ waves, $m=-1,0,+1$), `waves_vec_2 = {0}` (channel 2 $s$-wave only) — `L20/F3_cpu_openmp_v25_K3QC_cached_core.hpp:106-107`, `L20/K3df_minuit_fit_v32f_…hpp:149-150`. This matches the reference default `waves='sp'` for channel 1 and `'s'` for channel 2 (`REF/base_code/F3/K2i_mat.py:86-87`).

Harmonics (`L20/spherical_functions.h:15-89`): real, unnormalised "calligraphic" harmonics $\mathcal Y_{\ell m}(\mathbf p) = \sqrt{4\pi}\,|\mathbf p|^\ell Y_{\ell m}(\hat{\mathbf p})$ with
$\mathcal Y_{00}=1$, $\mathcal Y_{1,-1}=\sqrt3\,p_y$, $\mathcal Y_{10}=\sqrt3\,p_z$, $\mathcal Y_{11}=\sqrt3\,p_x$; $\ell=2$ also coded (lines 55-80). Identical to `REF/base_code/defns.py:424-432` (`y1real`).

#### 3.3.3 $F_2$ — sum minus integral (`L20/F2_functions_v2.h`)

For spectator flavour $i$ at momentum $\mathbf p$ (`F2_ang_mom`, lines 1130-1243), the code builds (all from the code):

$$
\gamma = \frac{E-\omega_p}{\sqrt{\sigma_p}},\quad x = \frac{q^*_p L}{2\pi},\quad \xi_{jk} = \tfrac12\left(1+\frac{m_j^2-m_k^2}{\sigma_p}\right),\quad \mathbf n_{P-p} = \frac{L}{2\pi}(\mathbf P-\mathbf p)
$$

($\xi_{jk}$ is my label for the code variable `xii`, line 947; it is unrelated to the anisotropy $\xi$.)

$$
\mathbf r(\mathbf n_a) = \mathbf n_a + \mathbf n_{P-p}\left[\frac{\mathbf n_a\cdot\mathbf n_{P-p}}{n_{P-p}^2}\left(\frac1\gamma - 1\right) - \frac{\xi_{jk}}{\gamma}\right]\quad(\text{lines }956\text{-}1029;\ \mathbf r=\mathbf n_a\text{ if }\mathbf n_{P-p}=0)
$$

$$
S_{\ell'm';\ell m} = \left(\frac{2\pi}{L}\right)^{\ell'+\ell}\sum_{\mathbf n_a\in[-N,N]^3}\frac{\mathcal Y_{\ell'm'}(\mathbf r)\,\mathcal Y_{\ell m}(\mathbf r)\,e^{\alpha(x^2-r^2)}}{x^2-r^2}\qquad(\text{`I\_sum\_ang\_mom`, lines 908-1126})
$$

with $N$ = `max_shell_num` = 20, i.e. $41^3 = 68\,921$ terms per element; the terms are collected and added with a magnitude-sorted Neumaier (compensated) sum in `long double` (`sorted_neumaier_sum`, lines 64-113, used at line 1104). $\alpha$ = 0.5 is the Gaussian ultraviolet damping parameter of the sum (and of the subtracted integral).

The subtracted principal-value integral is analytic (`I0F`, `I1F`, `I2F`, lines 629-748), using $\operatorname{erfi}$ from the MIT Faddeeva package (`Faddeeva::erfi`, lines 652, 695, 734):

$$
I_0 = 4\pi\gamma\left[-\tfrac12\sqrt{\tfrac{\pi}{\alpha}}\,e^{\alpha x^2} + \tfrac{\pi x}{2}\operatorname{erfi}\!\left(\sqrt{\alpha x^2}\right)\right]
$$
$$
I_1 = 4\pi\gamma\left[-\sqrt{\tfrac{\pi}{\alpha^3}}\,\tfrac{1+2\alpha x^2}{4}\,e^{\alpha x^2} + \tfrac{\pi x^3}{2}\operatorname{erfi}\!\left(\sqrt{\alpha x^2}\right)\right]\left(\tfrac{2\pi}{L}\right)^{2}
$$
$$
I_2 = 4\pi\gamma\left[-\sqrt{\tfrac{\pi}{\alpha^5}}\,\tfrac{3+2\alpha x^2+4\alpha^2x^4}{8}\,e^{\alpha x^2} + \tfrac{\pi x^5}{2}\operatorname{erfi}\!\left(\sqrt{\alpha x^2}\right)\right]\left(\tfrac{2\pi}{L}\right)^{4}
$$

(the $(2\pi/L)^{2\ell}$ factors are the `Q0norm=true` branch, lines 787-809). I re-derived $I_1$ and $I_2$ from $I_0$ by writing $r^2 = x^2-(x^2-r^2)$ and the Gaussian moments; they agree with the code and with `REF/base_code/F3/F_mat.py:149-157`. $I_\ell$ is diagonal: nonzero only for $\ell'=\ell$, $m'=m$ (lines 775-824).

$$
F^{(i)}_{\mathbf p\ell'm';\,\mathbf k\ell m} = \delta_{\mathbf p\mathbf k}\,\frac{H_i(\mathbf p)}{16\pi^2L^4\,\omega_p(E-\omega_p)}\Big[S_{\ell'm';\ell m} - \delta_{\ell'\ell}\delta_{m'm}\,I_\ell\Big]\qquad(\text{line }1217,\ 1240)
$$

Note that $S$ is **not** diagonal in $(\ell,m)$: at fixed $\mathbf p$ it mixes $s$ and $p$ waves when $\mathbf n_{P-p}\neq0$ or $m_j\ne m_k$. The prefactor matches `REF/base_code/F3/F_mat.py:192`.

Flavour matrix (`F2_2plus1_mat`, lines 1364-1424): $\mathbf F = \operatorname{diag}\big(F^{(1)}[m_i,m_j,m_k = M_1,M_1,M_2],\ F^{(2)}[M_2,M_1,M_1]\big)$, with no extra factor $1/2$ on channel 2 (same as the reference `F_full_2plus1_scratch`).

`Q0norm` meaning (lines 764-765 comment): `true` = harmonics carry physical momentum factors $|\mathbf a^*|^\ell$ (no division by $q^{*\ell}$); `false` = divide by $x^{\ell}$ to get the "$q^*$-normalised" convention. **Production uses `Q0norm = true`** (`F3_cpu_openmp_v25_…hpp:101`, `K3df_minuit_fit_v32f_…hpp:147`).

#### 3.3.4 $G$ — switch-state (exchange) propagator (`L20/G_functions_v2.h`)

"Switch-state" here means the one-particle-exchange term in which the spectator changes from particle $j$ (in) to particle $i$ (out). Element (`G_ij_lm`, lines 74-203):

$$
G^{(ij)}_{\mathbf p\ell'm';\,\mathbf k\ell m} =
\frac{\mathcal Y_{\ell'm'}(\mathbf k^*_{p})\;H_i(\mathbf p)\,H_j(\mathbf k)\;\mathcal Y_{\ell m}(\mathbf p^*_{k})}
{2\omega^{(i)}_pL^3\;\Big[(E-\omega^{(i)}_p-\omega^{(j)}_k)^2-(\mathbf P-\mathbf p-\mathbf k)^2-m_k^2\Big]\;2\omega^{(j)}_kL^3}
$$

where $\mathbf k^*_p$ = $\mathbf k$ boosted to the rest frame of $(E-\omega_p,\mathbf P-\mathbf p)$ (line 130) and $\mathbf p^*_k$ = $\mathbf p$ boosted to the rest frame of $(E-\omega_k,\mathbf P-\mathbf k)$ (line 131); $H_i(\mathbf p)$ uses pair masses $(m_j,m_k)$ and $H_j(\mathbf k)$ uses $(m_i,m_k)$ (lines 142-143). The exchanged particle has mass $m_k$. With `Q0norm=true` this is returned as is (line 188-191).

Flavour matrix (`G_2plus1_mat`, lines 245-473):

$$
\mathbf G = \begin{pmatrix} G^{(11)}[M_1,M_1,M_2] & \sqrt2\,(-1)^{\ell'}\,G^{(12)}[M_1,M_2,M_1]\\[2pt] \sqrt2\,G^{(21)}[M_2,M_1,M_1]\,(-1)^{\ell} & 0\end{pmatrix}
$$

(lines 315, 376-378, 441-443; the $(-1)^\ell$ is the parity matrix "$P_L$"). The (2,2) block is zero because there is only one pion. Same structure as `REF/base_code/F3/G_mat.py:128` (`G12 = sqrt(2) * PL_1 @ Gmat_ij(...)`).

#### 3.3.5 $\mathcal K_2^{-1}$ (`L20/K2_functions_v2.h`)

Per-row element (`K2_inv_ERE_ang_mom`, lines 37-169; `K2inv_EREord2_i_mat`, lines 176-256), diagonal in momentum and $(\ell,m)$:

$$
\big[\mathcal K_2^{(i)}\big]^{-1}_{\mathbf k\ell m} = \frac{1}{2\omega_kL^3}\cdot\frac{\eta_i}{8\pi\sqrt{\sigma_k}}\Big[B_\ell(q^{*2}) + C_\ell\Big],\qquad C\text{ built from } |q^*|(1-H_i(\mathbf k))
$$

- $\ell=0$ (lines 77-95): $B_0 = -\frac1{a_0} + \frac12 r_0 q^{*2} + P_0 q^{*4}$, $\;C_0 = |q^*|\,(1-H)$.
- $\ell=1$ (lines 96-114): $B_1 = -\frac{1}{a_1 q^{*2}} + \frac12 r_1 + P_1 q^{*2}$, $\;C_1 = q^{*2}\,|q^*|\,(1-H)$.
- $\ell=2,3$ coded (lines 115-148) but unused in production.

Here $|q^*|$ is `abs(sqrt(q2))`, i.e. a real magnitude even below the pair threshold (line 57). Parameter array layout: `scatter_params_i[ell][0..2]` = $(a_\ell, r_\ell, P_\ell)$. Config keys `scatter{i}_{ell}{n}` are read in `L20/debug_v30r_…cpp:315`. Hence:

| Key | Meaning (from code) | Production value (lattice units of $a_t$) |
|---|---|---|
| `scatter1_00` | channel 1 ($K\pi$ pair), $\ell=0$, $a_0$ | 4.04 |
| `scatter1_10` | channel 1 ($K\pi$ pair), $\ell=1$, $a_1$ | $-43.2$ |
| `scatter2_00` | channel 2 ($KK$ pair), $\ell=0$, $a_0$ | 4.12 |
| all $r_\ell$, $P_\ell$ | effective ranges / shape terms | 0 |

(`L20/F3_cpu_openmp_v25_…hpp:133-136`, `L20/K3df_minuit_fit_v32f_…hpp:151-152`, `.../KKpi-L20-.../README.md:55-56`.) Sign convention: $q^*\cot\delta_0 = -1/a_0$; with this sign $a_0>0$ means a repulsive $s$-wave (my reading of the convention, not stated in code). The old V3 drivers used $a_0^{KK}=4.07$ with the comment "total uncertainty 0.07 stat" and $a_0^{K\pi}=4.04$ with "total uncertainty 0.05 stat 0.15 systematic" (`V3/source/generate_eigen_based_F3inv.cpp:33-34`); no uncertainty is carried in the newest code.

Flavour matrix (`K2inv_EREord2_2plus1_mat`, lines 264-324):
$\mathbf K_2^{-1} = \operatorname{diag}\big(\mathcal K^{(1)-1},\ 2\,\mathcal K^{(2)-1}\big)$ (line 317-318), matching `REF/base_code/F3/K2i_mat.py:88`.

**Important consistency issue for $\ell=1$** — see §3.6 item B1.

#### 3.3.6 $F_3$ and $F_3^{-1}$

All production paths build (`L20/F3_cpu_openmp_v25_…hpp:427, 476-479`; `L20/K3df_minuit_fit_v32f_…hpp:632-635`):

$$
\mathbf H = \mathbf K_2^{-1}+\mathbf F+\mathbf G,\qquad \mathbf X = \mathbf H^{-1}\mathbf F\ (\text{partial-pivot LU solve}),\qquad \mathbf F_3 = \frac{\mathbf F}{3} - \mathbf F\,\mathbf X .
$$

This is the nondegenerate flavour-matrix version (block matrices above), identical to `REF/base_code/F3/F3_mat.py:27`.

$F_3^{-1}$ in the v32f/v33 path (`K3df_minuit_fit_v32f_…hpp:649-672`): scale $s=\max_{ab}|F_{3,ab}|$, invert $\mathbf F_3/s$ by partial-pivot LU, then $\mathbf F_3^{-1} = (\mathbf F_3/s)^{-1}/s$; cache both the full inverse and $V^\dagger\mathbf F_3^{-1}V$. The comments stress that the projected quantity is $V^\dagger F_3^{-1}V$, **never** $(V^\dagger F_3 V)^{-1}$ (lines 223-226, 649-651).

Quantization condition (`F3_cpu_openmp_v25_…hpp:12-15, 496-503`; L20 README lines 44-49):

$$
\det\Big[V_{\rm sel}^\dagger\big(\mathbf F_3^{-1}+\mathbf K_{\rm df,3}\big)V_{\rm sel}\Big]\Big/\;(\text{scale}) = 0,\qquad \text{scale}=(L_{a_s}\,\xi)^{6}\ \text{applied to the matrix before the determinant (README line 47)}.
$$

Multiplying the whole matrix by a constant does not move its zeros; the scale only avoids overflow. The v32f fitter also uses "eigenvalue of the projected QC matrix closest to zero" sign changes (`K3df_minuit_fit_v32f_…hpp:20-23, 552-568`) — the root-finding/classifier logic belongs to another report.

Old diagnostic quantity still present: $F_{3,\rm iso} = \mathbf v^T\mathbf F_3\mathbf v$ with $\mathbf v = (1,\dots,1,\ 1/\sqrt2,\dots,1/\sqrt2)$ (`L20/QC_functions_v2.h:226-250`).

#### 3.3.7 $\mathcal K_{\rm df,3}$ parametrisation (`L20/K3_functions_2plus1.hpp`)

$M = 2M_1+M_2 = 2m_K+m_\pi$; $\Delta = E^{*2}/M^2 - 1$ (line 455; same as `REF/base_code/Kdf3/K3main.py:15-17`). For row flavour $i$ (spectator momentum $\mathbf p$, $(\ell',m')$) and column flavour $j$ ($\mathbf k$, $(\ell,m)$):

$$
\mathcal K^{(ij)}_{\rm df,3} = c_{ij}\Big[\delta_{\ell'0}\delta_{\ell0}\big(K^{\rm iso}_0 + K^{\rm iso}_1\Delta\big) + K_B\,\mathcal K_B^{(ij)} + K_E\,\mathcal K_E^{(ij)}\Big],\qquad c_{11}=1,\ c_{12}=c_{21}=\tfrac1{\sqrt2},\ c_{22}=\tfrac12
$$

(lines 436-482, 484-491, 543; block factors identical to `K3main.py:73`). The four fit parameters are `K3iso0`, `K3iso1`, `K3B`, `K3E` (`K3df_minuit_fit_v32f_…hpp:97-102`). There is no $1/L^3$ or $1/(2\omega)$ factor on $\mathcal K_{\rm df,3}$.

The code does not use Mandelstam-$t$ names; it writes the $\mathcal K_B$, $\mathcal K_E$ pieces directly in momenta. Selected pieces (each divided by $M^2$; code lines given):

- $\mathcal K_B^{(22)}$, $s$–$s$ (lines 196-205): $\frac{2}{M^2}\big[E^{*2} - E(\omega_p+\omega_k) + \mathbf P\cdot(\mathbf p+\mathbf k) + M_2^2 - 4M_1^2\big]$.
- $\mathcal K_B^{(11)}$, $s$–$s$ (lines 207-214): $\frac{1}{M^2}\big[E(\omega_p+\omega_k) - \mathbf P\cdot(\mathbf p+\mathbf k) - 6M_1^2 + \omega_{p'}\,\Delta\omega_p + \omega_{k'}\,\Delta\omega_k\big]$, with $\Delta\omega_p = \sqrt{q^{*2}_p+M_1^2}-\sqrt{q^{*2}_p+M_2^2}$ (code `pms_0`, line 179) and $\omega_{p'}$ the energy of $\mathbf p$ boosted into the $(E-\omega_p,\mathbf P-\mathbf p)$ frame (lines 180-181).
- $\mathcal K_B^{(11)}$, $p$–$s$ / $s$–$p$ (lines 215-222): $-\frac{2}{3M^2}\mathcal Y_{1m'}(\mathbf p')$, $-\frac{2}{3M^2}\mathcal Y_{1m}(\mathbf k')$.
- $\mathcal K_B^{(12)}$, $s$–$s$ (lines 228-235): $\frac1{M^2}\big[E\omega_p - \mathbf p\cdot\mathbf P + E^{*2} - 2E\omega_k + 2\mathbf P\cdot\mathbf k + \omega_{p'}\Delta\omega_p + M_2^2 - 7M_1^2\big]$; $\mathcal K_B^{(21)}$ is the mirror (lines 245-252).
- $\mathcal K_E^{(22)}$, $s$–$s$ (lines 364-370): $\frac{2}{M^2}\big[M_2^2 - \omega_p\omega_k + \mathbf p\cdot\mathbf k\big]$.
- $\mathcal K_E^{(12)}$, $s$–$s$ (lines 401-409): $\frac1{M^2}\big[2M_2^2 - (E-\omega_p)\omega_k + (\mathbf P-\mathbf p)\cdot\mathbf k + \Delta\omega_p\,\omega_{k^*}\big]$.
- $\mathcal K_E^{(11)}$ has $s$–$s$, $s$–$p$, $p$–$s$ and $p$–$p$ pieces built from double boosts, $\gamma$ factors and a $3\times3$ tensor $t$ (lines 291-397); too long to transcribe here — the code says it follows `K3E.py` (line 135, 350).

Defaults: the fitter starts from guess $0.1$ for all four (`K3df_minuit_fit_v32f_…hpp:160`); the CPU core default is $K^{\rm iso}=(200,400)$, $K_B=K_E=0$ "simple isotropic example used in QC3_release" (`F3_cpu_openmp_v25_…hpp:127-131`); the L20 README baseline (lines 58-62) is $(K^{\rm iso}_0, K^{\rm iso}_1, K_B, K_E) = (73735.84,\ -972421.14,\ 347174.06,\ -1226756.71)$ in lattice units (no uncertainties given there).

C++ vs Python $\mathcal K_3$ agreement was checked with `CG/legacy_source/compare_K3_cpp_vs_python.py` + `test_K3mat_2plus1_cpp_dump.cpp`; I found no stored result of that comparison in this area.

#### 3.3.8 Units and constants

| Quantity | Value | Source |
|---|---|---|
| anisotropy $\xi = a_s/a_t$ | 3.444 ("found from lattice", no uncertainty) | `V3/source/generate_eigen_based_F3inv.cpp:29`; `F3_cpu_openmp_v25_…hpp:99` |
| $a_t m_K$ | 0.09698 | `F3_cpu_openmp_v25_…hpp:92` |
| $a_t m_\pi$ | 0.06906 | `F3_cpu_openmp_v25_…hpp:91` |
| $L/a_s$ | 20 (also 24 in some drivers) | `F3_cpu_openmp_v25_…hpp:100` |
| $L$ used in all formulas | $L = \xi\cdot(L/a_s)$, i.e. in units of $a_t$ | `PhysicsParams::L()`, line 139 |
| momenta | $a_t\,\mathbf k = \frac{2\pi}{\xi\,(L/a_s)}\mathbf n$ | `functions.h:1032-1034` with that $L$ |
| $\alpha$, $\epsilon_H$, $N$ | 0.5, 0, 20 | lines 94-96 |
| three-particle threshold | $a_t(2m_K+m_\pi) = 0.26302$ (my arithmetic from the masses) | — |
| fitted energy window | $a_tE^*\in[0.261,\ 0.36]$ | `K3df_minuit_fit_v32f_…hpp:124-125` |

Everything is in units of the temporal lattice spacing $a_t$; the anisotropy only converts $L/a_s$ into $a_t$ units. Lattice ensemble string: `szscl21_24_128_b1p50_t_x4p300_um0p0840_sm0p0743_n1p265` (`K3df_minuit_fit_v32f_…hpp:169`).

#### 3.3.9 Irreps and the projector `Vsel`

Irrep labels (`L20/debug_v30r_…cpp:66-74`):

| Label | $\mathbf d$ (units $2\pi/L$) | Little group (code list, `projections_v1.hpp:248-305`) | Irrep string passed |
|---|---|---|---|
| `000_A1m` | (0,0,0) | $O_h$, 48 signed permutations | `A1u` ($A_1^-$) |
| `100_A2` | (0,0,1) | $C_{4v}$, 8 elements | `A2` |
| `110_A2` | (1,1,0) | $C_{2v}$, 4 elements | `A2` |
| `111_A2` | (1,1,1) | $C_{3v}$, 6 elements | `A2` |
| `200_A2` | (0,0,2) | $C_{4v}$, 8 elements | `A2` |

Group elements are signed permutations $R$ of $(x,y,z)$ acting as `cubic_transf` (`projections_v1.hpp:174`); characters $\chi_I(R)$ come from `chi` (`projections_v1.hpp:611`).

Production projector (`L20/fv_projector_cartesian_l1_v30m.hpp`, called with a **fixed** convention by `L20/debug_v30r_…cpp:534-584`):

$$
P_I = \frac{d_I}{|LG|}\sum_{R\in LG}\chi_I(R)\,U(R),\qquad U(R)=S(R)\,W(R),\qquad P_I\leftarrow\tfrac12(P_I+P_I^\dagger)
$$

(lines 257-280), block-diagonal in flavour (`U_2plus1`, lines 240-255), where
- $S(R)$: $|\mathbf n,\ell,m\rangle\to|R\mathbf n,\ell,m\rangle$ (lines 145-168);
- $W(R)$: $|\mathbf n,\ell,m\rangle\to\Pi(R)\sum_{m'}D^\ell_{m'm}(R)|\mathbf n,\ell,m'\rangle$ (lines 198-226), with $\Pi(R)=-1$ for improper $R$ when `parity = -1` (line 136-143, "Eq. (2.54)"). `parity=-1` is the intrinsic parity of three pseudoscalars $(-1)^3$.
- $D^0=1$; $D^1$ is the **Cartesian signed-permutation matrix** in the real basis $m=+1\to x$, $m=0\to z$, $m=-1\to y$ (lines 171-196). This is the "v30m correction"; $\ell\ge2$ falls back to `real_wigner_d::D_real_element`.
- Convention used in production: index 0 = momentum $R$ (not $R^{-1}$), $D(R)$ not transposed, order $S\cdot W$, parity mode 1 (`debug_v30r_…cpp:542-548`; GPU port documents the same at `CG/gpu_component_cachegen/source/projections_gpu_v32y_exact_device.cuh:6-11`). The function `best_projector` (lines 334-358) that searches 48 conventions is diagnostic only.

`Vsel` construction (`projections_v1.hpp:1554-1640`, `build_projector_from_eigenvectors_near_one`, then `debug_v30r_…cpp:575-580`):
1. Hermitian eigen-decomposition of $P_I$.
2. Keep eigenvectors whose eigenvalue $\lambda$ satisfies $|\lambda-1|\le$ `eig_tol` $=0.05$.
3. Normalise each kept column if its norm exceeds `norm_tol` $=10^{-12}$.
4. Re-orthonormalise with modified Gram–Schmidt (`g_vsel_reorthonormalize = 1`, line 203).
5. `proj_tol` $=10^{-10}$ is only used to print whether $V V^\dagger$ is Hermitian and idempotent (debug output); it does not change `Vsel`. The GPU host path even marks it unused (`F3_gpu_cachegen_v32zu.cu:460`).

Diagnostics stored with each cache row: $\|P_I^2-P_I\|/\|P_I\|$, representation closure, unitarity, eigenvalue gap (`debug_v30r_…cpp:468-488, 549-573`). Projector results are cached by a string key of the basis (`basis_signature_v30q`, lines 517-532), since `Vsel` changes only when a new momentum shell enters.

### 3.4 Question 3 — variants and dead ends

| File(s) | What it tried | Status |
|---|---|---|
| `F2_functions.h`, `G_functions.h`, `K2_functions.h`, `QC_functions.h` (all packages, identical) | First $s$-wave-only implementation. `F3_ID` / `F3_ID_mat` (`QC_functions.h:34-159`) = identical-particle test vs arXiv:1803.04169, with $F$ multiplied by $0.5$ and $L^3$ rescalings. `F3_ND_2plus1_mat` (line 196) = first 2+1 version, already with $\operatorname{diag}(K^{-1}_1, 2K^{-1}_2)$ and $\sqrt2$ on $G_{12},G_{21}$. Two variants of the sum (`I00_sum_F` vs `I00_sum_F_1`, "_1 includes the lattice anisotropy xi", `F2_functions.h:194`). Many pole/"additional poles" test functions (lines 1460-2245). | Superseded by `_v2`. Cannot be included together with `_v2` (both define `I0F`, `LinearSolver_3`, …). |
| `V3/source/generate_eigen_based_F3inv.cpp` + `QC_fitter_1.py` (`QC3_bissection_eigen_based`, line 2163) | Python fitter calls a C++ binary per energy that prints $1/F_{3,\rm iso}$; bisection on $F_{3,\rm iso}^{-1} + K^{\rm iso}$. | Dead end (one process per energy; $s$-wave iso only). |
| `splines.h`, `splines_test.cpp`, `generate_spline_based_F3inv.cpp`, `printer_F3_omp.cpp:test_F3inv_with_splines` | Cubic-spline interpolation of $F_3^{-1}(E)$ on a grid ($F(x)=\sum_i S_{ij}(x)F(x_i)$, `splines.h:1-7`) to avoid recomputation during fits. | Replaced by linear `interp1d` in Python, later by binary caches. |
| `generate_K3iso.cpp` + `K3iso_maker.py` | "Inverse" approach: for each lattice level and jackknife sample compute $K^{\rm iso} = -F_{3,\rm iso}^{-1}(E_{\rm lat})$ (`generate_K3iso.cpp:76`), then fit $K^{\rm iso}(E)$. | Superseded by direct QC fits. |
| `generate_F3iso.cpp` + `F3iso_maker.py` | Tabulate $F_{3,\rm iso}$ at lattice energies. | Diagnostic. |
| `generate_pole.cpp`, `pole_maker.py` (1543 lines), `pole_searching.h` | Secant search for poles of $F_3^{-1}$ / $\det F_3$ (`pole_searching.h:11,116`); compare with free spectrum vs $L$ (line 258-260). | Dead end for the $s$+$p$ code. |
| `QC_fitter_1.py` (3114 lines) | iminuit fits of 1 or 2 parameters, $K^{\rm iso}=K_0 + K_1\,(s-s_{\rm th})/s_{\rm th}$ (line 97), multi-$L$, covariance from `lattice_data_covariance`. | Ancestor of the C++ Minuit2 fitter; same $\Delta$ definition. |
| `projections.hpp`, `projections_from_config.hpp`, `eigenvalue_tracker.hpp` | C++ translation of `REF` `projections.py` + `group_theory_defns.py` (`projections.hpp:1-3`), reorder `config_maker_4` output into momentum-major order (`projections_from_config.hpp:57-60`), track eigenvalues vs energy. Used by `V3/source/test_gpu.cpp`. | Superseded by `projections_v1.hpp` + `fv_projector_…v30m.hpp`. Conflicts with `projections_v1.hpp` (both define `PI`, `chi`, `irrep_dim`). |
| `projections_v1.hpp::P_irrep_projection_2plus1` (lines 1141-1468) | First production projector: real Wigner $D$ (`real_wigner_d.hpp`) plus an "orbit" equality check. | Replaced in v30n/v30m: "replaces legacy P_irrep_projection_2plus1 path by fv_projector_cartesian_l1_v30m.hpp — uses the fixed U(R) with exact closure/P_I idempotency for s+p mixing" (`debug_v30r_…cpp:35-37`). Still compiled (header supplies `little_group`, `chi`, `Vsel` builder). |
| `functions.h`: `dawson_func`, `ERFI_func` (lines 180-200) | $\operatorname{erfi}$ by a $10^6$-step Riemann sum. | Replaced by Faddeeva. |
| `functions.h`: `config_maker`, `_1`, `_2`, `_3`, `_4`, `_5` | Successive basis builders: $k<k_{\max}$; $H>$tol with anisotropy; with integer list; with $(\ell,m)$; with $E_{\min}$ rule ($\ell m$ outermost); `_5` "not finished yet!!! DO NOT USE" (line 1123). | Production uses `config_maker_4_momentum_first` only. |
| `printer_F3_acc.cpp` (OpenACC), `gpu_solvers_batched_streams*.cpp` (cuSOLVER batched QR), `gpu_varsize_batched_inverse.cu`, `test_gpu.cpp` | Early GPU attempts. | Superseded by the `CG` GPU pipeline (other report). |
| `non_int_spectrum_maker_3body_trivial_irrep.py`, `functions.h:threebody_non_int_spectrum*` | Free three-particle energies (comment: "the formulation … that REDSTAR uses", line 1373-1374). | Still useful for plots; not physics-core. |

### 3.5 Question 4a — fixes visible by diffing old (V3) vs newest

1. **$I_1$ integral fixed.** Old `V3/source/F2_functions_v2.h:81-82`: `term = pi/alpha^3; B = -sqrt(pi/term)*…` gives $\sqrt{\alpha^3}$ instead of $\sqrt{\pi/\alpha^3}$. New `L20/F2_functions_v2.h:685-686` uses `sqrt(term)` — correct (I re-derived it, §3.3.3). Any old $p$-wave $F$ numbers from V3 are wrong.
2. **Missing $\sqrt2$ on $G_{12}$ fixed.** Old `V3/source/G_functions_v2.h:250` had `G12 = projector*G…` (no $\sqrt2$) while $G_{21}$ had it (line 298) — a non-symmetric $\mathbf G$. New `L20/G_functions_v2.h:378` has $\sqrt2$.
3. **Index bug in $K_2^{-1}$ fixed.** Old `V3/source/K2_functions_v2.h:202-203` read `klm_config[3][i]` for the column; new lines 204-205 read `[j]`. Harmless in practice (only $i=j$ is used) but correct now.
4. **$\ell=1$ cutoff term changed** from $|q^*|(1-H)$ (old line 105) to $q^{*2}|q^*|(1-H)$ (new line 106, comment "stable equation"). See B1 below — this change made the $\ell=1$ formula only half-converted.
5. **Sign convention of the $\mathbf r$ shift.** Old used $\mathbf n_{p-P}$ with $+\xi_{jk}/\gamma$ (old lines 272, 315); new uses $\mathbf n_{P-p}$ with $-\xi_{jk}/\gamma$ (new lines 956-958, 1001-1002, comment "changed the sign … according to the python code"). The two are algebraically identical (flipping $\mathbf n_{P-p}$ flips only the $\xi$ term), so this is **not** a physics change. It matches `REF/base_code/F3/F_mat.py:81`, where `alpha_ij` $=\omega^*_j/\sqrt{\sigma}$ equals $\xi_{jk}$.
6. **Compensated summation** (Neumaier, sorted by magnitude, `long double`) added to the $F$ sum (new lines 64-113, 1104).
7. **Projector replaced** by the Cartesian $\ell=1$ version (§3.3.9, §3.4).
8. **Solver:** old `QC_functions_v2.h` formed $\mathbf H^{-1}$ explicitly with QR (`V3/source/QC_functions_v2.h`, `LinearSolver_4`); new code solves $\mathbf H\mathbf X=\mathbf F$ by LU and inverts a max-scaled $F_3$.

### 3.6 Question 4b — problems I found by reading (not noted in the code)

Severity labels are mine: **B** = can change physics output in production; **M** = minor / unused path; **C** = cosmetic.

- **B1 — $\ell=1$ $\mathcal K_2^{-1}$ is dimensionally mixed under `Q0norm=true`.** `L20/K2_functions_v2.h:102,106` (and the GPU copy `L20/K2_functions_gpu_safe.cuh:414-415`) compute
  $B_1 + C_1 = \left(-\frac{1}{a_1q^{*2}} + \frac{r_1}{2} + P_1q^{*2}\right) + q^{*2}|q^*|(1-H)$.
  The first bracket is the $q^*$-normalised form ($q^{*3}\cot\delta_1/q^{*2}$, units of momentum), the second is the `Q0norm` form (units of momentum$^3$). The reference code, written for the same `Q0` harmonics ("Q matrices already applied"), uses $q^{*3}\cot\delta_1 + q^{*2}|q^*|(1-H)$ with $q^{*3}\cot\delta_1 = -1/a_1$ (`REF/base_code/F3/K2i_mat.py:26-27`, `REF/base_code/F3/qcot_fits.py:25-26`). The commented-out line 103 (`-1/a1 + 0.5*r1*q2 + P1*q2^2`) is the reference-consistent form. With the default $r_1=P_1=0$ the C++ code effectively uses $q^{*3}\cot\delta_1 = -1/(a_1q^{*2})$ instead of $-1/a_1$, i.e. a different $p$-wave threshold behaviour, so `scatter1_10 = -43.2` does not mean the same as the reference's $a_1$. This affects every production $F_3$ through channel 1. **This needs Digonto's decision** (it may have been intentional for numerical stability — the comment "stable equation" hints at that — but it does not match the cited reference). I did not run anything to quantify the effect.
- **M1 — `Q0norm=false` branch of $G$ is inconsistent with $F$.** `L20/G_functions_v2.h:194-195` divides by $(q^{*2})^{\ell}$ = $q^{*2\ell}$ per side, while $F$ divides by $x^{\ell}$ per side ($q^{*\ell}$ up to $(2\pi/L)^\ell$). Unused in production (`Q0norm=true`); do not expose this flag in the new repo without fixing it.
- **M2 — $\ell=2,3$ $\mathcal K_2^{-1}$** (`K2_functions_v2.h:115-148`) are in the $q^*$-normalised convention only; inconsistent with `Q0norm=true`. Unused.
- **M3 — `config_maker_2` fallback bug**: `L20/functions.h:473-475` pushes three zeros into `n_config[0]` instead of `[0],[1],[2]`. Unused in production.
- **M4 — `threebody_Gpoles`**: `L20/functions.h:1564` computes `sqrt(pole_energy2*pole_energy1 - P^2)`; should be `pole_energy2*pole_energy2`. Diagnostic only.
- **M5 — old `tilde_K2_00`** (`V3/source/K2_functions.h:96`, identical in all packages) tests $|\mathbf p|=|\mathbf k|$ instead of $\mathbf p=\mathbf k$ for the Kronecker delta. Old $s$-wave family only.
- **M6 — `F3iso` diagnostic vector** (`L20/QC_functions_v2.h:226-250`) puts weight 1 on the $p$-wave rows too; with $s$+$p$ bases this is not an isotropic ($\ell=0$) projection. Only used in diagnostics/plots ("F3inv_isotropic" scripts).
- **M7 — CPU/GPU numerical differences.** GPU $F$ uses a power series for $\operatorname{erfi}$ (80 terms, tol $10^{-15}$; `CG/gpu_component_cachegen/source/F2_gpu_safe_builder.cuh:21-24, 313-354`), which loses accuracy for large $\alpha x^2$, and plain (uncompensated) summation (lines 555-596) versus Neumaier on CPU; the GPU header itself asks for element-wise validation. The GPU path also casts `prod1` to its real part (line 577). Any new repo needs a CPU–GPU agreement test with a stated tolerance.
- **M8 — K2 below pair threshold** uses $|q^*|$ = `abs(sqrt(q2))` (`K2_functions_v2.h:57`), i.e. the real magnitude of an imaginary $q^*$. This matches the reference (`q_abs = sqrt(abs(q2))`, `K2i_mat.py:18`), so it is a convention, not a bug — record it in METHODS.
- **C1** — debug prints show `G11(i,j)` while labelled G12/G21 (`L20/G_functions_v2.h:394, 459`).
- **C2** — `spherical_functions.h:20` defines `norm = sqrt(4π)` but never uses it (the $\sqrt{4\pi}$ is absorbed in the $\sqrt3$ coefficients — correct).
- **C3** — `projections_v1.hpp:2407-2413` has a commented-out `main()`; harmless.

### 3.7 Question 5 — recommended minimal physics-core file set

Take every file from `L20/` (byte-identical to the v33o/v33p/v33n/v32w/v32y and `CG/legacy_source` copies):

**Must have (header-only core + one compiled unit):**

| # | Full path | Role |
|---|---|---|
| 1 | `/home/digimonk/Projects/three-body-qc-codes/extracted/KKpi-L20-user-validated-hotwindow-fitter-main/KKpi-L20-user-validated-hotwindow-fitter-main/source/functions.h` | kinematics, $H$, boost, basis builder |
| 2 | `.../KKpi-L20-user-validated-hotwindow-fitter-main/source/spherical_functions.h` | real $\mathcal Y_{\ell m}$ |
| 3 | `.../KKpi-L20-user-validated-hotwindow-fitter-main/source/F2_functions_v2.h` | $F$ (keep lines 1-1424; lines 1426-2376 are debug helpers) |
| 4 | `.../KKpi-L20-user-validated-hotwindow-fitter-main/source/K2_functions_v2.h` | $\mathcal K_2^{-1}$ (fix or document B1) |
| 5 | `.../KKpi-L20-user-validated-hotwindow-fitter-main/source/G_functions_v2.h` | $G$ |
| 6 | `.../KKpi-L20-user-validated-hotwindow-fitter-main/source/QC_functions_v2.h` | solvers, invertibility checks, $F_{3,\rm iso}$ |
| 7 | `.../KKpi-L20-user-validated-hotwindow-fitter-main/source/K3_functions_2plus1.hpp` | $\mathcal K_{\rm df,3}$ |
| 8 | `.../KKpi-L20-user-validated-hotwindow-fitter-main/source/fv_projector_cartesian_l1_v30m.hpp` | production $P_I$ |
| 9 | `.../KKpi-L20-user-validated-hotwindow-fitter-main/source/projections_v1.hpp` | `little_group`, `chi`, `cubic_transf`, `is_in_rotations_list`, `build_projector_from_eigenvectors_near_one` (the old `P_irrep_projection_2plus1` inside can be deleted later) |
| 10 | `.../KKpi-L20-user-validated-hotwindow-fitter-main/source/real_wigner_d.hpp` | included by #9; $\ell\ge2$ fallback in #8 |
| 11 | `.../KKpi-L20-user-validated-hotwindow-fitter-main/source/dig_tools.hpp` | included by #9 (printing helpers only) |
| 12 | `.../KKpi-L20-user-validated-hotwindow-fitter-main/source/Faddeeva.cc` and `Faddeeva.hh` | $\operatorname{erfi}$; MIT licence (header of `Faddeeva.cc`: "Copyright (c) 2012 Massachusetts Institute of Technology", http://ab-initio.mit.edu/Faddeeva). `Faddeeva.h` (C99 API) is not needed. |

**Should accompany (to be refactored, not copied verbatim):**

| # | Full path | What to extract |
|---|---|---|
| 13 | `.../KKpi-L20-user-validated-hotwindow-fitter-main/source/F3_cpu_openmp_v25_K3QC_cached_core.hpp` | `PhysicsParams` (lines 89-140) and the $F_3$ assembly (lines 377-503). The rest is scan/refinement logic for another module. |
| 14 | `.../KKpi-L20-user-validated-hotwindow-fitter-main/source/debug_v30r_svd_validated_cartesian_l1_100_A2_zero_scan.cpp` | Only `parse_label` (66-74), the config reader (≈250-320) and the projector cache `CachedProjectorV30q` (468-607). Today these live in a 2203-line debug program that is `#include`d with `main` renamed (`K3df_minuit_fit_v32f_…hpp:60-62`); move them into a small header such as `projector_cache.hpp` (name is my suggestion). |
| 15 | `.../KKpi-L20-user-validated-hotwindow-fitter-main/source/K3df_minuit_fit_v32f_fullF3inv_QCfull_cached_classifier.hpp` | `build_cache_entry` (589-686) is the cleanest end-to-end "energy → $F_3^{-1}$, `Vsel`" function. |

**Reference material to keep (read-only, for tests):** `/home/digimonk/Projects/three-body-qc-codes/extracted/_nested_old_versions/QC3_release/QC3_release/` (Python reference by Blanton–Romero-López–Sharpe for arXiv:2111.12734) and `/home/digimonk/Projects/three-body-qc-codes/extracted/full_cache_generator_package_20260629-main/full_cache_generator_package_20260629-main/legacy_source/compare_K3_cpp_vs_python.py` — these allow element-wise acceptance tests of $F$, $G$, $\mathcal K_2^{-1}$, $\mathcal K_3$.

**Leave out of the core:** `F2_functions.h`, `G_functions.h`, `K2_functions.h`, `QC_functions.h`, `pole_searching.h`, `splines.h`, `projections.hpp`, `projections_from_config.hpp`, `eigenvalue_tracker.hpp`, `Faddeeva.h`, all `V3/source/generate_*.cpp`, `printer_*.cpp`, `test*.cpp`, Python scripts in `V3/source/` (archive them under `legacy/`).

---

## 4. What this means

- There is effectively **one** physics core in the newest packages; all copies are byte-identical, so consolidation is about removing old copies, not merging diverged code.
- The newest core reproduces the structure of the published reference (BRS `QC3_release`): same $H$, same harmonics, same $F$ prefactor and integrals, same $\sqrt2$/parity structure of $\mathbf G$, same $\operatorname{diag}(1,2)$ for $\mathbf K_2^{-1}$, same $\mathcal K_3$ blocks and $\Delta$.
- The one place where it departs from the reference in a way that changes production numbers is the $\ell=1$ $\mathcal K_2^{-1}$ (B1). Any consolidated repo should either fix it (with a re-run of the fits) or document it as a deliberate modelling choice.

## 5. Problems and assumptions

- I did not open arXiv:2105.12094 or arXiv:2111.12734; equation-number claims are copied from code comments. The paper identifications (authors/topics) in §3.1 are from my memory and flagged as such.
- I did not compile or run anything, so no numerical comparison (e.g. the size of the B1 effect, CPU vs GPU agreement) is provided; every number in this report is a constant read from the files cited next to it.
- The reference's $p$-wave convention was read from `REF/base_code/F3/K2i_mat.py` and `qcot_fits.py`; I assumed the reference's `Q0` harmonics are the same as the C++ `Q0norm=true` ones (both use $\sqrt3\,p_{x,y,z}$ without $q^*$ normalisation; `REF/base_code/defns.py:424-432` vs `L20/spherical_functions.h:36-47`).
- The physical justification for treating $K\bar K$ as an identical pair ($\eta_2=1/2$, $\sqrt2$ factors, $\mathcal K_3$ factors $1/\sqrt2$, $1/2$) is not documented in the code.
- `hashes.txt` omits `.hh`/`.cuh`; their hashes in §3.2 come from my own `md5sum` runs.
- Labels I introduced myself: $\xi_{jk}$ (code `xii`), "switch-state" in §3.3.4 only as a gloss of the user's term, severity tags B/M/C, and the suggested header name `projector_cache.hpp`.

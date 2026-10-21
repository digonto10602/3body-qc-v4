# 001 — M1: a cache-free, CPU-only, single-energy three-body QC with an oracle and an independent reference

- **Date:** 2026-10-21
- **Git commit at planning:** `1eaf1fa`
- **Trigger:** kickoff on `prompts/000-consolidation-and-roadmap.md` ("Update 2026-10-21", decisions D1–D5 in `notes/STATUS.md`) and Digonto's clarification relayed by the coordinator during the kickoff: *"graphical manipulation" means **symbolic/analytic manipulation**, not GPU; no GPU at all; optimise for 1 CPU core; compare every variant against an **oracle**; the solution must work for every energy in the production range and every irrep.*
- **Research:** `reports/research/007-single-energy-qc-methods-and-reference.md` (all numbers quoted below come from it unless a command is given).
- **Supersedes:** M1 of `prompts/000` (steps 1–8 there are re-cut below; step 6 "restore inputs" is gone per D4; no GPU step).

Jargon used here: **QC** = quantization condition $\det Q(E)=0$; **oracle** = the existing C++ path run unchanged on 1 core at 1 energy (Digonto's baseline); **ingredient** = one of $F$, $G$, $\mathcal K_2^{-1}$, $\mathcal K_{\rm df,3}$; **sector** = one $(L/a_s,\ \mathbf d,\ \Gamma)$, e.g. `L20/100_A2`; **point** = one $(\text{sector},E^*)$; **anchor** = a stored oracle value with provenance; **reference** = the planner's NumPy implementation written from the published formulas (`tests/acceptance/m1/qc3ref/`).

## 0. Diagnosis (why the old pipeline is slow and what is wrong)

1. **Cost.** On 1 core the oracle needs 2.8 s ($N=56$) to 6.5 s ($N=100$–160) per energy (`reports/research/007` §3.1). Almost all of it is the $F$ lattice sum: for every spectator momentum and every $(\ell'm';\ell m)$ pair (16 per channel-1 momentum) it loops over the full $41^3$ cube in complex arithmetic. A NumPy version doing one pass per momentum takes 0.2 s for the same $F$, and the cube is >95 % wasted because the summand decays as $e^{-\alpha r^2}$ (§3.5). The old answer was to precompute $F_3^{-1}$ on a 20 000-point grid per sector on a GPU and cache it (GBs); D4 removes that option.
2. **Correctness debts.** (a) The lattice covariance is wrong by $(N-1)^2$ in the multi-L reader (D1 confirmed). (b) The p-wave $\mathcal K_2^{-1}$ mixes two normalisation conventions (§3.6 of research 007: the scattering-volume term is $q$-normalised, the $(1-H)$ term and $F,G$ are $Q_0$-normalised). (c) No stored evidence that $\mathcal K_{\rm df,3}$, $F$, $G$ match an independent implementation. (d) Absolute `/home/digonto` paths.
3. **What is right and must be kept.** The physics core (`src/core`) reproduces the published formulas: the planner's independent NumPy reference agrees with it to $10^{-15}$–$10^{-16}$ (relative Frobenius) for $F$, $G$, $\mathcal K_2^{-1}$ and $\mathcal K_{\rm df,3}$ (research 007 §3.2–3.3). All matrices are real (max $|\mathrm{Im}|=0$), $F_3^{-1}$ is Hermitian to $3\times10^{-19}$.

## 1. Goal of M1

A production path `qc3.fast.evaluate_point(sector, Ecm, c)` that, **with no disk cache**, returns $(\mathrm{sign},\ln|\det Q|)$ and the ingredients at any energy in the production grid range $[0.26310, 0.36]$ for any of the five L20 sectors, reproduces the oracle within the tolerances of §4 and is $\ge 20\times$ faster on 1 core (target after all steps: $\ge 100\times$, i.e. 5–50 ms per point), plus the corrected lattice-data reader (D1, D3), the independent reference checks, and a benchmark protocol. Root finding without labels, the fitter and the L20 refits are **M2** (outline in §9).

## 2. Method

### 2.1 Equations and conventions (unchanged physics; see `docs/METHODS.md` §§1–8 for the full set)

All quantities in units of $a_t$; $L=\xi L/a_s=68.88$; $\mathbf P=\frac{2\pi}{L}\mathbf d$; $E=\sqrt{E^{*2}+\mathbf P^2}$. Channel 1: spectator $K$ ($m_i=m_K$), pair $(K,\pi)$, $\eta_1=1$, waves $\{0,1\}$. Channel 2: spectator $\pi$, pair $(K,K)$, $\eta_2=\tfrac12$, waves $\{0\}$. Source for everything in this subsection: Blanton, Romero-López, Sharpe, arXiv:2111.12734 (JHEP 02 (2022) 098), equation numbers in parentheses.

- Basis (`docs/METHODS.md` §2): spectator momenta $\mathbf n\in[-20,20]^3$ kept iff $E>E^+_{\min}(\mathbf k)=\omega_k+\sqrt{\sigma_{\min}+(\mathbf P-\mathbf k)^2}$, $\sigma_{\min}=|m_j^2-m_k^2|+2z_{\min}(m_j+m_k)\min(m_j,m_k)$, $z_{\min}=0.02$; **order:** loop $n_x,n_y,n_z$ each from $-20$ to $20$ (outer→inner), for each kept momentum the $(\ell,m)$ list $(0,0),(1,-1),(1,0),(1,1)$ for channel 1 and $(0,0)$ for channel 2; channel 1 rows first, then channel 2. The fast path must return **exactly this basis in this order** (then matrices are comparable element-wise).
- $F$ (Eq. 8 in the UV-regulated form of App. A.2):
$$
F^{(i)}_{\mathbf p\ell'm';\mathbf p\ell m}=\frac{H^{(i)}(\mathbf p)}{16\pi^2L^4\omega_p(E-\omega_p)}\Big[S_{\ell'm';\ell m}-\delta_{\ell'\ell}\delta_{m'm}I_\ell\Big],\qquad
S=\Big(\tfrac{2\pi}{L}\Big)^{\ell'+\ell}\sum_{\mathbf n\in[-20,20]^3}\frac{\mathcal Y_{\ell'm'}(\mathbf r)\mathcal Y_{\ell m}(\mathbf r)e^{\alpha(x^2-r^2)}}{x^2-r^2},
$$
$$
\mathbf r=\mathbf n_\perp+\gamma^{-1}(\mathbf n_\parallel-\xi_{jk}\mathbf n_P),\quad \mathbf n_P=\tfrac{L}{2\pi}(\mathbf P-\mathbf p),\quad \gamma=\tfrac{E-\omega_p}{\sqrt\sigma},\quad \xi_{jk}=\tfrac12\Big(1+\tfrac{m_j^2-m_k^2}{\sigma}\Big),\quad x=\tfrac{q^*L}{2\pi},
$$
with $\mathbf n_\parallel,\mathbf n_\perp$ the components of $\mathbf n$ along/perpendicular to $\mathbf n_P$ ($\mathbf r=\mathbf n$ if $\mathbf n_P=0$); $x$ is complex below the pair threshold ($x^2<0$) and every formula is evaluated with complex $x$ (the result is real). $I_0=4\pi\gamma[-\tfrac12\sqrt{\pi/\alpha}\,e^{\alpha x^2}+\tfrac{\pi x}{2}\mathrm{erfi}(\sqrt{\alpha x^2})]$, $I_1=4\pi\gamma[-\sqrt{\pi/\alpha^3}\tfrac{1+2\alpha x^2}{4}e^{\alpha x^2}+\tfrac{\pi x^3}{2}\mathrm{erfi}(\sqrt{\alpha x^2})](2\pi/L)^2$ (App. A.2; $I_0$ from arXiv:1901.07095 App. B). $\mathcal Y_{00}=1$, $\mathcal Y_{1,-1}=\sqrt3 r_y$, $\mathcal Y_{10}=\sqrt3 r_z$, $\mathcal Y_{11}=\sqrt3 r_x$ (App. A.1, $Q_0$ convention of Eqs. 20–21). $\alpha=0.5$, $\epsilon_H=0$, $H=J(z)$ with $z=(\sigma-|m_j^2-m_k^2|)/((m_j+m_k)^2-|m_j^2-m_k^2|)$ (Eqs. 22, 27).
- $G$ (Eq. 12), $\mathbf G$ flavour matrix with $\sqrt2$ and $(-1)^{\ell}$ (Eqs. 5, 7); $\mathcal K_2^{-1}$ (Eqs. 14–15) with the $1/(2\omega_kL^3)$ factor and $\mathrm{diag}(1,2)$ (Eq. 6); $F_3=F/3-FH^{-1}F$, $H=\mathcal K_2^{-1}+F+G$ (Eq. 3); $\mathcal K_{\rm df,3}$ (Eqs. 28–31) exactly as `docs/METHODS.md` §7. Projector $P_\Gamma$ and $V_{\rm sel}$ as §8 (Eqs. 54–59).
- **QC (invariants 1–2):** $Q=(L\xi)^{-6}V_{\rm sel}^\dagger\big(F_3^{-1}+\mathcal K_{\rm df,3}\big)V_{\rm sel}$; $(\mathrm{sign},\ln|\det Q|)$ by partial-pivot LU. $V^\dagger F_3^{-1}V$ may be computed as $V^\dagger X$ with $X$ the LU **solution** of $F_3X=V$ (still the projection of the inverse; never $(V^\dagger F_3V)^{-1}$). $F_3$ is scaled by $s=\max|F_{3,ab}|$ before the LU, as now.

### 2.2 p-wave $\mathcal K_2^{-1}$: two named forms, production default unchanged until D6

$$
[\mathcal K_2^{(1)-1}]_{\ell=1}=\frac{1}{2\omega_kL^3}\frac{\eta_1}{8\pi\sqrt\sigma}\times\begin{cases}
-\dfrac{1}{a_1q^{*2}}+q^{*2}|q^*|(1-H) & \texttt{coded\_v2}\ \text{(production now; } a_1=-43.2)\\[6pt]
-\dfrac{1}{a_1}+q^{*2}|q^*|(1-H) & \texttt{q0\_consistent}\ \text{(2111.12734 Eqs. 15, 17, 20–21 with } Q_0\text{; same } a_1=-43.2)
\end{cases}
$$
HadSpec (arXiv:1411.2004 Sec. VI, Eq. 9: $k^{2\ell+1}\cot\delta_\ell=1/a_\ell+\dots$) gives $a_0=(-4.04\pm0.05\pm0.15)\,a_t$, $a_1=(43.2\pm3.7\pm15.4)\,a_t^3$ for $I=3/2$ $\pi K$ on this ensemble; the code's $-1/a$ sign convention maps these to $4.04$ and $-43.2$. `q0_consistent` reproduces HadSpec exactly; `coded_v2` overstates the p-wave $\mathcal K_2^{-1}$ by $1/q^{*2}$. The interface exposes `Sector.k2_pwave_form`; **default stays `coded_v2`** (invariant 3) until Digonto decides D6 on the basis of step 9's level shifts. Both forms are acceptance-tested against the reference.

### 2.3 Speed-up candidates, ranked (expected gain, how verified)

Every rung must pass `tests/acceptance/m1/test_m1_fast_vs_oracle.py` before the next one is attempted; the benchmark (§7) records the gain of each rung separately. Rungs 1–4 are required (step 6), rung 5 is step 7, rungs 6–7 are optional within step 7.

| # | technique | what it exploits | expected gain | verification |
|---|---|---|---|---|
| 1 | **One lattice-sum pass per spectator momentum** producing all $(\ell'm';\ell m)$ entries at once (accumulate $\mathcal Y_a(\mathbf r)\mathcal Y_b(\mathbf r)w(\mathbf r)$ for the 4 harmonics in one loop) | the oracle recomputes the same sum 16× | ≈10–16× on $F$ | $F$ vs oracle $\le10^{-12}$ |
| 2 | **Ellipsoid-within-cube truncation:** keep $\mathbf n\in[-20,20]^3$ with $r^2\le R^2$, $R^2=\mathrm{Re}\,x^2+\ln(1/\epsilon_{\rm tr})/\alpha$, $\epsilon_{\rm tr}=10^{-16}$ | each omitted term is $\le e^{-\alpha(R^2-\mathrm{Re}x^2)}/(R^2-\mathrm{Re}x^2)$ in magnitude; the kept set is a subset of the oracle's cube, so the oracle is reproduced to $<10^{-14}$ for every $\gamma$ (research 007 §3.5). Enumerate $\mathbf n$ by bounding box of the ellipsoid, test $r^2\le R^2$ | ≈10–20× on $F$ ($\sim3000$–$6000$ terms instead of 68 921) | same; plus unit test: truncation error vs full cube $\le10^{-13}$ at 20 random $(E,\mathbf p)$ |
| 3 | **Real arithmetic** (double, not complex) for all ingredients, $H^{-1}F$, $F_3$, $F_3^{-1}$, projector, $Q$; LU `dgetrf`-type | every matrix is real symmetric (App. A.1; measured max $\lvert\mathrm{Im}\rvert=0$); complex $x$ only inside the scalar $I_\ell$ and $w$ (whose results are real) | 2–4× on everything | identities test (reality, symmetry) + oracle comparison |
| 4 | **In-memory energy-independent pieces:** momentum/basis lists per $E$-segment, $V_{\rm sel}$ per basis signature (as the oracle's `CachedProjectorV30q`), harmonics tables; $V^\dagger F_3^{-1}V$ via LU solve with $N_{\rm proj}$ RHS; symmetric fill of $G$, $\mathcal K_{\rm df,3}$ (compute $a\le b$, mirror) | energy independence; symmetry | removes the $P_\Gamma$ eigen-solve per point (ms–10 ms); ≈2× on $G,\mathcal K_3$, LU | oracle comparison |
| 5 | **Class moment tensors (symbolic reduction of the sum):** for fixed (sector, $\mathbf p$) group the cube points by the pair of integers $(n^2,\ \mathbf n\cdot\mathbf n_P)$. Within a class, $r^2=n^2+2s\,(\mathbf n\cdot\mathbf n_P)+s^2n_P^2$ with $s=\frac{\mathbf n\cdot\mathbf n_P}{n_P^2}(\gamma^{-1}-1)-\xi_{jk}\gamma^{-1}$ is constant, so $w(\mathbf r)$ is constant and $\sum_{\rm class}\mathcal Y_a(\mathbf r)\mathcal Y_b(\mathbf r)=$ a quadratic form in the **energy-independent** class moments $C_0=\sum1$, $\mathbf C_1=\sum\mathbf n$, $C_2=\sum\mathbf n\otimes\mathbf n$ ($\mathbf r\otimes\mathbf r=\mathbf n\otimes\mathbf n+s(\mathbf n\otimes\mathbf n_P+\mathbf n_P\otimes\mathbf n)+s^2\mathbf n_P\otimes\mathbf n_P$). Precompute the class table once per (sector, momentum) lazily; per energy loop over classes with $r^2\le R^2$ | the sum is a sum over $O(10^3)$ classes of a scalar weight times a precomputed $4\times4$ block | ≈3–10× more on $F$ (→ $\lesssim1$ ms); memory per sector $\lesssim$ tens of MB if classes are stored only for $r^2$-sorted prefixes | oracle comparison; unit test: class sums equal the direct ellipsoid sum to $10^{-14}$ at 20 random points |
| 6 | $\mathbf P=0,\ \mathbf p=0$ special case (then $\mathbf r=\mathbf n$, $F$ diagonal in $(\ell m)$, sum over $n^2$ shells only) | cubic symmetry | small | oracle comparison |
| 7 | OpenMP over spectator momenta ($F$) and over rows ($G$, $\mathcal K_3$) | 4 cores here | reported separately; **not** the comparison metric | oracle comparison with `OMP_NUM_THREADS=4` |

Not pursued in M1 (recorded for M2+): Woodbury/low-rank updates of $F_3^{-1}$ between nearby energies inside a root refinement; derivative $dQ/dE$ by the same class tables (free with rung 5) for Newton steps (`docs/ROOT_FINDING.md` §3); little-group block structure of $H$ before the solve (the QC must still be built from the full $F_3^{-1}$, invariant 1; a block solve of $H$ is allowed only if proven equivalent).

### 2.4 Lattice data (D1, D3)

Files `data/lattice/L20/20_<ddd>_<irrep>_n<k>.jack`: header `N 1 0 0 1`, then $N=603$ lines `0 x_c` with $x_c=a_tE_{\rm lab}$ on configuration $c$ (rescaled-ensemble values; D1). Recipe (identical to the trusted v32w reader `src/fitter/lattice_data_covariance_cpp.hpp`, research 004 §3.2):
$$
y_c=\frac1{N-1}\sum_{c'\ne c}x_{c'},\qquad y_c\to E^*_c=\sqrt{y_c^2-P^2},\quad P=\frac{2\pi|\mathbf d|}{\xi L/a_s},\qquad
\bar E^*=\frac1N\sum_cE^*_c,\qquad C_{ij}=\frac{N-1}{N}\sum_c(E^{*(i)}_c-\bar E^{*(i)})(E^{*(j)}_c-\bar E^{*(j)}).
$$
For a linear map this equals $C_{ij}=\frac{1}{N(N-1)}\sum_c(x^{(i)}_c-\bar x^{(i)})(x^{(j)}_c-\bar x^{(j)})$, i.e. $\sigma=\mathrm{s.d.}/\sqrt N$. File-name irrep renaming `001→100`, `002→200`, `011→110`. A level is kept iff $\bar E^*\le E_{\rm cut}$. Cross-$L$ blocks are zero. $\chi^2=r^TC^+r$ with $C^+$ the eigen-pseudoinverse dropping $\lambda\le10^{-12}\lambda_{\max}$ (as now). **D3 matching:** sort both lists ascending; with $m$ lattice and $n$ QC levels, $n\ge m$ pairs the lowest $m$; $n<m$ raises `MissingRootError`.

## 3. Milestones and steps

Each step: builder implements → tester (mode `step`) runs the listed ids → commit. Test ids are `tests/acceptance/m1/<file>::<test>`. Commands: fast part `python3 -m pytest tests/acceptance/m1 -m "not slow" -q` (< 2 min on this container once the oracle is bound; ≈10 oracle calls, cached per module), slow part `python3 -m pytest tests/acceptance/m1 -m slow -q` (milestone gate and after steps 6–7), collection check `python3 -m pytest --collect-only -q tests/acceptance` (33 tests). Before any code exists, the suite must fail only with `NotImplementedError` — verified by the planner with a throwaway stub package (`scratch/planner/stubs/qc3`, 33/33 `NotImplementedError`, 0 other errors); the builder's real stubs in `src/python/qc3/` must give the same result in step 1.

### M1 — cache-free single-energy QC (this plan)

| step | what | done-check | acceptance ids to turn green |
|---|---|---|---|
| **1** Build + stubs | `./build.sh` gets targets `core` (compiles `src/qc3_core` as the pybind11 module `qc3._core` into `src/python/qc3/`), `pyext` alias, and `all` no longer requires Minuit2 (fitter target only on demand). `pip install pybind11` allowed. Create the Python package `src/python/qc3/` with every module/function of §5 as stubs raising `NotImplementedError("qc3.<module>.<fn>: not implemented (prompts/001 step N)")`; `qc3.constants` may already hold the real constants. `tests/acceptance/conftest.py` already adds `src/python` to `sys.path`. | `python3 -c "import qc3, qc3.oracle, qc3.fast, qc3.lattice, qc3.qc, qc3.roots, qc3.bench"` works; `pytest --collect-only tests/acceptance` collects without errors; `pytest tests/acceptance/m1 -x -q` fails only with `NotImplementedError`; `./build.sh core` compiles (an empty module is fine at this step). | none green yet; `test_m1_portability.py::test_no_absolute_home_paths_in_m1_code` should pass immediately |
| **2** Oracle wrapper | `qc3.oracle.evaluate_point`: bind the **unchanged** physics of `build_component_entry` (`src/cachegen_cpu/generate_cpu_component_cache.cpp` L223–288: basis, $F$, $\mathcal K_2^{-1}$, $G$, $H^{-1}F$, $F_3$, scaled-LU $F_3^{-1}$, projector cache, $V^\dagger F_3^{-1}V$) plus `K3mat_2plus1` for the four basis matrices ($c=e_i$), with per-stage timings, into `PointResult`; single-threaded; `metadata` per §5.4. Oracle code lives in `src/qc3_core/oracle_point.hpp` and must `#include` the existing headers, not copy formulas. | `pytest tests/acceptance/m1/test_m1_oracle.py -q` (fast part) passes; slow anchors pass (`-m slow`) | `test_m1_oracle.py::test_oracle_anchors_small_points`, `::test_oracle_anchors_k3_baseline_small_points`, `::test_oracle_reality_hermiticity_isometry`, `::test_oracle_metadata_has_production_constants`, `::test_oracle_basis_rule_and_order`, `::test_oracle_timings_present`; slow: `::test_oracle_anchors_all_points` |
| **3** Reference checks | Nothing new in `src/` except `Sector.k2_pwave_form` handling in the oracle (both forms; `coded_v2` default). The reference package is planner-owned. | `pytest tests/acceptance/m1/test_m1_reference.py -q` passes (fast and slow) | `test_m1_reference.py::test_F_G_K2inv_vs_reference_small_points`, `::test_K2inv_pwave_both_forms_vs_reference`, `::test_K3_basis_vs_reference_small_points`, `::test_K3_basis_symmetric_and_linear`; slow: `::test_ingredients_vs_reference_large_points` |
| **4** QC helpers + lattice reader (D1, D3) | `qc3.qc` (`project_Q`, `slogdet`), `qc3.lattice` (§5.6), `qc3.roots.refine_root` (§5.5). Pure NumPy/SciPy. Unit tests in `tests/unit/`. | `pytest tests/acceptance/m1/test_m1_qc_identities.py tests/acceptance/m1/test_m1_lattice.py -q` | `test_m1_qc_identities.py::*` (4 tests), `test_m1_lattice.py::*` (8 tests) |
| **5** Portability | `qc3.lattice.default_data_dir()` honours `QC_DATA_ROOT` (default `<repo>/data`); no `/home/digonto`, `/media/digonto` in `src/python`, `src/qc3_core`, `build.sh`, `configs/qc3/`. Legacy configs/scripts are left untouched (they are not on the M1 path); list them in the step report. | `pytest tests/acceptance/m1/test_m1_portability.py -q` | `test_m1_portability.py::*` (2 tests) |
| **6** Fast path v1 (rungs 1–4) | `qc3.fast.evaluate_point` in `src/qc3_core/fast_point.hpp` (+ bindings): same basis builder as the oracle (shared code), one-pass ellipsoid-truncated $F$, real arithmetic, projector cache, LU-solve projection, symmetric fills. `options` keys at least: `trunc_eps` (default $10^{-16}$), `use_classes` (False until step 7), `threads` (1). Unit test for the truncation bound. | fast part of `test_m1_fast_vs_oracle.py` passes; slow part passes; benchmark report `reports/benchmarks/001-*` written by the tester (protocol §7) showing ≥ 20× single-core at every benchmark point | `test_m1_fast_vs_oracle.py::test_fast_vs_oracle_small_points`, `::test_fast_basis_identical_to_oracle`, `::test_fast_metadata_and_timings`; slow: `::test_fast_vs_oracle_all_anchor_points`, `::test_fast_vs_oracle_random_energies_all_irreps`, `::test_fast_root_vs_oracle_100_A2_windows`, `::test_fast_speedup_single_thread_ge_20x` |
| **7** Fast path v2 (rung 5, optional 6–7) | Class moment tables (`use_classes=True` becomes the default once it passes); then OpenMP as a separate, non-default option. Stop when the per-point single-core time at the largest benchmark point ($N\approx160$) is ≤ 20 ms or a rung gains < 1.5×. | same tests as step 6 (all green with the new default); benchmark `reports/benchmarks/002-*` ≥ 100× single-core at the small points and ≥ 50× at the largest | same ids as step 6 |
| **8** Profile + end-to-end proxy | `qc3.bench.run` + `python -m qc3.bench` implementing §7; tester runs it (mode `benchmark`) after the auditor passes. Includes the **χ² proxy**: time to refine the 12 labelled L20 roots (cutoff 0.335, brackets from `data/windows/L20/*/accepted_windows.csv`) at the baseline $\vec c$ with `qc3.roots.refine_root` on the fast path, 1 core. Compare with 1.4 s (old direct FCN) and 24 ms (hot-window FCN, needed GBs of cache). | `reports/benchmarks/00N-*.md` with the table of §7; `results/bench_*.json` | — (benchmark, not a test) |
| **9** p-wave quantification for D6 | With the fast path: the 12 L20 levels (cutoff 0.335) at the baseline $\vec c$ for `coded_v2` and `q0_consistent` (refine_root in each accepted window; if a root leaves its window, widen by $\pm250$ rows and report); table of $\Delta E$ vs $\sigma_{\rm lat}$ from `qc3.lattice`. Written by the builder as `results/k2_pwave_shift_L20.json` + a short section in its return; the coordinator puts it to Digonto. | file exists with 12 rows, both forms, lattice σ | — |

**Milestone gate (M1):** tester `full` (all m1 tests incl. slow, unit tests, `acceptance-check`), auditor, tester `benchmark`, tag `milestone-1`.

## 4. Acceptance tests and tolerances (all in `tests/acceptance/m1/`, locked)

| test (file::id) | checks | tolerance | justification (research 007 §3.4 unless stated) |
|---|---|---|---|
| `test_m1_oracle.py::test_oracle_anchors_small_points` | 3 anchors (100_A2 0.2915561, 000_A1m 0.26932, 200_A2 0.28664; $\mathcal K_3=0$): $N_1,N_2,N,N_{\rm proj}$ equal; sign equal; $\lvert\Delta\ln\lvert\det Q\rvert\rvert$ | $10^{-10}$ | same code and compiler flags as the anchor run; only the LU library differs (Eigen vs LAPACK) → $\sim10^{-13}$ measured; margin $10^3$ |
| `::test_oracle_anchors_all_points` (slow) | all 12 anchors | $10^{-10}$ | same |
| `::test_oracle_anchors_k3_baseline_small_points` | 3 anchors with $\vec c_{\rm base}$ | $10^{-9}$ | anchor $\mathcal K_3$ comes from the reference, which agrees with the C++ port to $9\times10^{-16}$; sensitivity $\le100\varepsilon$ → $10^{-13}$; margin $10^4$ |
| `::test_oracle_reality_hermiticity_isometry` | max $\lvert\mathrm{Im}\rvert$ of $F,G,\mathcal K_2^{-1},\mathcal K_3$ basis $\le10^{-15}\|\cdot\|_\infty$; $\|F_3^{-1}-F_3^{-1\dagger}\|/\|F_3^{-1}\|\le10^{-14}$; $\|V^\dagger V-1\|\le10^{-12}$; $\lvert\mathrm{Im}\det\rvert$ consistent with sign $\in\{\pm1\}$ | as listed | invariant 5; measured $3\times10^{-19}$, $8\times10^{-16}$ (ingredients exactly real) |
| `::test_oracle_metadata_has_production_constants` | metadata equals invariant-3 constants, `Lbyas`, `irrep`, `k2_pwave_form`, `backend`, `git_commit` | exact | invariant 3 |
| `::test_oracle_basis_rule_and_order` | every basis momentum obeys $E>E^+_{\min}$, none missing, order rule of §2.1 | exact (integers) | definition |
| `test_m1_reference.py::test_F_G_K2inv_vs_reference_small_points` | rel. Frobenius vs `qc3ref` at 3 small points | $10^{-12}$ | measured $10^{-16}$–$10^{-15}$; margin $\ge10^3$ |
| `::test_K2inv_pwave_both_forms_vs_reference` | both forms, 100_A2 point | $10^{-12}$ | closed form; measured $10^{-16}$ |
| `::test_K3_basis_vs_reference_small_points` | 4 basis matrices (iso0, iso1, B, E) vs the angular-projection reference, rel. Frobenius, 2 small points | $10^{-12}$ | quadrature exact for degree ≤ 2 polynomials; measured $8.7\times10^{-16}$ |
| `::test_K3_basis_symmetric_and_linear` | each basis matrix symmetric; $\mathcal K_3(\vec c)=\sum c_iB_i$ | $10^{-12}$ rel. | exact properties |
| `::test_ingredients_vs_reference_large_points` (slow) | 110_A2 0.35000, 111_A2 0.35103 | $10^{-12}$ | as above |
| `test_m1_qc_identities.py` (4 tests) | `project_Q`/`slogdet` vs NumPy on seeded random matrices; invariance under unitary mixing of $V$ columns; $(L\xi)^{-6}$ scaling moves $\ln\lvert\det\rvert$ by $-6N_{\rm proj}\ln(L\xi)$; sign flips with one negative eigenvalue | $10^{-10}$ abs. on $\ln\lvert\det\rvert$ | float64 LU on $N_{\rm proj}\le40$: $\sim10^{-13}$ |
| `test_m1_lattice.py` (8 tests) | synthetic file: mean and $\sigma=\mathrm{s.d.}/\sqrt N$ exact; jackknife formula; $E_{\rm lab}\to E^*$ order; real file `20_000_A1m_n0.jack` vs independent NumPy; $\sigma\in[4,7]\times10^{-4}$; level counts 6/9/12/14 at cutoffs 0.315/0.325/0.335/0.345; covariance symmetric PSD; D3 matching; χ² pseudo-inverse | $10^{-12}$ rel.; counts exact | D1; counts from research 004 Table 1 (recomputed in-test from the files) |
| `test_m1_portability.py` (2 tests) | no absolute `/home/digonto`, `/media/digonto` in M1 code; `QC_DATA_ROOT` honoured | exact | step 5 |
| `test_m1_fast_vs_oracle.py::test_fast_vs_oracle_small_points` | 3 small points, $\mathcal K_3=0$ and $\vec c_{\rm base}$: sign equal; $\lvert\Delta\ln\lvert\det Q\rvert\rvert$; $F,G,\mathcal K_2^{-1},\mathcal K_3$ basis rel. Frobenius; $F_3^{-1}$ rel. Frobenius | $10^{-6}$; $10^{-12}$; $10^{-9}$ | $\ln\det$: LU round-off floor $1.4\times10^{-9}$ at $\kappa(F_3)=4.8\times10^8$, margin 700 (covers $\kappa\lesssim10^{11}$); ingredients: truncation bound $10^{-14}$, summation-order $\lesssim10^{-14}$, margin 100; $F_3^{-1}$: measured response $2\times10^{-11}$ to $10^{-12}$ ingredient noise, margin 50 |
| `::test_fast_basis_identical_to_oracle` | basis arrays equal | exact | §2.1 |
| `::test_fast_vs_oracle_all_anchor_points` (slow) | 12 points, live oracle, as above | as above | as above |
| `::test_fast_vs_oracle_random_energies_all_irreps` (slow) | seed 20261021: 2 uniform energies per irrep in $[0.2631,0.36]$ (re-drawn if within $5\times10^{-5}$ of a dimension jump, detected by $N$ changing within $\pm10^{-5}$), live oracle | as above | generalisation requirement (Digonto) |
| `::test_fast_root_vs_oracle_100_A2_windows` (slow) | `refine_root` on the 2 accepted 100_A2 windows at $\vec c_{\rm base}$ with oracle and fast | $\lvert\Delta E_0\rvert\le10^{-9}$; both roots inside their window | a relative $10^{-10}$ change of $Q$ moves a simple root by $\lesssim10^{-13}$; lattice σ $5\times10^{-4}$; grid spacing $4.8\times10^{-6}$ |
| `::test_fast_speedup_single_thread_ge_20x` (slow) | same process, 1 thread, 110_A2 at 0.35 ($N=108$), median of 3: $t_{\rm oracle}/t_{\rm fast}$ | $\ge20$ | ratio on the same machine; rungs 1–2 alone give $\ge100\times$ on $F$ |
| `::test_fast_metadata_and_timings` | metadata constants as oracle, `backend='fast'`, `trunc_eps`, timings keys present and positive | exact | §5.4 |

Tolerance changes go through the planner only. If a test cannot be met because the **oracle** is found wrong (e.g. an energy where the oracle's cube is not converged matters), stop and report; do not loosen.

### Amendment 2026-10-21 (acceptance-v2): `test_m1_qc_identities.py::test_scale_applies_to_matrix_not_determinant`, last block

Approved by Digonto 2026-10-21; narrow correction of a defective test construction, no tolerance or data change.

- **What was wrong.** The last block means to flip the sign of one eigenvalue of the seeded Hermitian $M$ (`_synthetic(seed=3)`, $N=30$) and check that the LU sign flips. It wrote `w2[0] = -abs(w2[0])`. For seed 3 the eigenvalues are $w_{0..3}=(-6.341,-4.939,-4.616,-4.085)$ (10 of 30 negative, $\min|w|=0.51$), so $w_0<0$ already, `M2 == M` up to rounding, and `assert sa == -sb` was unsatisfiable for *every* correct `slogdet`. The test was wrong, not the code (reproduced by the coordinator and the planner).
- **The fix.** `w2[0] = -w2[0]` (a genuine sign flip; $\|M_2-M\|/\|M\|=0.37$), and the block now also records $\ln|\det|$ of both matrices and asserts $|\ell_a-\ell_b|\le$ `TOL` $=10^{-10}$. Comments in the test explain the construction.
- **Why it does not weaken the test.** (i) The intended check (one sign flip of a full-rank Hermitian matrix flips the LU sign) is now actually exercised; before, it could only fail. (ii) The added $|\det|$-equality check is exact in exact arithmetic ($\det M_2=-\det M$ since $M_2=U\,\mathrm{diag}(w_0\to-w_0)\,U^\dagger$, same scale $(L\xi)^{-6N}$ on both); numerically the reconstruction $U\,\mathrm{diag}(w_2)\,U^\dagger$ and the LU each contribute $\mathcal O(N\varepsilon)\sim10^{-14}$; measured $|\ell_a-\ell_b|=3.4\times10^{-13}$ against `TOL` $=10^{-10}$, margin 300, consistent with the existing §4 justification for this file ("float64 LU on $N_{\rm proj}\le40$: $\sim10^{-13}$"). No tolerance was changed. (iii) Verified against deliberately wrong implementations (throwaway `scratch/planner/verify_sign_flip_amendment*.py`, git-ignored): with NumPy's LU `slogdet` (the builder's `qc3.qc.slogdet`) the test passes; `sign` $\equiv+1$ fails (line 56, the earlier projected-sign check, because the seed-3 projected determinant is negative); `sign` $\equiv-1$, which passes every earlier check of this test, is caught only by the amended `assert sa == -sb` (line 67); `sign = sign(tr Q)` and an $N_{\rm proj}$-rescaled $\ln|\det|$ are caught by the earlier checks (lines 56 and 55). So the amended block adds discriminating power (a constant sign that happens to match the seed-3 projected sign) that the rest of the test does not have.
- **Collection count unchanged:** 33 tests in `tests/acceptance`. No other acceptance test, tolerance or data file was touched. Re-lock as `acceptance-v2`.

## 5. Interface (what the acceptance tests import; builder creates stubs in step 1)

Package root: `src/python/qc3/` (added to `sys.path` by `tests/acceptance/conftest.py`). Compiled extension: `src/python/qc3/_core*.so` from `src/qc3_core/` via `./build.sh core`. All arrays `numpy.ndarray`, float64 or complex128 (tests cast to complex128), C-contiguous, row-major, basis-ordered as §2.1. Energies $E^*$ in $a_t^{-1}$.

### 5.1 `qc3.constants`
```python
PRODUCTION: dict = {"xi": 3.444, "atmK": 0.09698, "atmpi": 0.06906, "eta_1": 1.0, "eta_2": 0.5,
    "alpha": 0.5, "epsilon_H": 0.0, "max_shell_num": 20, "waves_1": (0, 1), "waves_2": (0,),
    "parity": -1, "eig_tol": 0.05, "scatter": {"a0_Kpi": 4.04, "a1_Kpi": -43.2, "a0_KK": 4.12},
    "Q0norm": True, "z_min": 0.02}
C_BASELINE: tuple = (73735.840894011912, -972421.14060757787, 347174.05548116949, -1226756.7068845264)  # label-generation vector (research 004 §3.6)
IRREPS: tuple = ("000_A1m", "100_A2", "110_A2", "111_A2", "200_A2")
NNP: dict = {"000_A1m": (0,0,0), "100_A2": (0,0,1), "110_A2": (1,1,0), "111_A2": (1,1,1), "200_A2": (0,0,2)}
GRID_E0, GRID_DE, GRID_N = 0.26310, 0.0969/19999, 20000
```

### 5.2 `qc3.sector`
```python
@dataclass(frozen=True)
class Sector:
    irrep: str; Lbyas: int = 20; xi: float = 3.444; atmK: float = 0.09698; atmpi: float = 0.06906
    eta_1: float = 1.0; eta_2: float = 0.5; alpha: float = 0.5; epsilon_H: float = 0.0
    max_shell_num: int = 20; waves_1: tuple = (0, 1); waves_2: tuple = (0,)
    parity: int = -1; eig_tol: float = 0.05
    a0_Kpi: float = 4.04; a1_Kpi: float = -43.2; a0_KK: float = 4.12   # r_l = P_l = 0
    k2_pwave_form: str = "coded_v2"          # or "q0_consistent" (§2.2)
    @classmethod
    def production(cls, irrep: str, Lbyas: int = 20, **overrides) -> "Sector": ...
    @property
    def L(self) -> float: ...                 # xi * Lbyas
    @property
    def nnP(self) -> tuple: ...               # integer d vector
    @property
    def P_lab(self) -> np.ndarray: ...        # 2*pi/L * d
```

### 5.3 `qc3.point`
```python
BASIS_DTYPE = np.dtype([("flavour","i4"),("nx","i4"),("ny","i4"),("nz","i4"),("ell","i4"),("m","i4")])
@dataclass
class PointResult:
    sector: Sector; Ecm: float; E_lab: float
    N1: int; N2: int; N: int; Nproj: int
    basis: np.ndarray            # BASIS_DTYPE, length N, order of §2.1
    c: tuple                      # (Kiso0, Kiso1, KB, KE) used for K3, Q, sign, lnabsdet
    F: np.ndarray | None; G: np.ndarray | None; K2inv: np.ndarray | None
    F3: np.ndarray | None; F3inv: np.ndarray | None; Vsel: np.ndarray | None     # Vsel: N x Nproj
    K3_basis: dict | None         # {"iso0","iso1","B","E"} -> N x N, so that K3(c) = sum c_i * K3_basis[i]
    K3: np.ndarray | None; Q: np.ndarray | None                                   # Q: Nproj x Nproj
    sign: int; lnabsdet: float
    metadata: dict; timings: dict  # see 5.4
ALL = ("F","G","K2inv","F3","F3inv","Vsel","K3","Q")
```

### 5.4 `qc3.oracle` and `qc3.fast`
```python
def evaluate_point(sector: Sector, Ecm: float, c=(0.0,0.0,0.0,0.0), want=ALL) -> PointResult: ...        # qc3.oracle
def evaluate_point(sector: Sector, Ecm: float, c=(0.0,0.0,0.0,0.0), want=ALL, options: dict | None = None) -> PointResult: ...  # qc3.fast
```
`sign` and `lnabsdet` are always filled (for the given `c`). `metadata` must contain: every key of `PRODUCTION` with the sector's values, `"Lbyas"`, `"irrep"`, `"nnP"`, `"L"`, `"k2_pwave_form"`, `"backend"` (`"oracle"` or `"fast"`), `"git_commit"` (short hash at build time), `"threads"`, and for `fast` also `"trunc_eps"`, `"use_classes"`. `timings` (seconds, float) keys: `"basis","F","G","K2inv","F3","F3inv","projector","K3","projection","lu","total"` (a stage not run → 0.0). `fast` `options` keys: `trunc_eps` (float, default 1e-16), `use_classes` (bool), `threads` (int, default 1).

### 5.5 `qc3.qc` and `qc3.roots`
```python
def project_Q(F3inv, Vsel, K3, L: float) -> np.ndarray      # (L^-6) V^dag (F3inv + K3) V ; L = xi*Lbyas (so L^6 == (Lbyas*xi)^6)
def slogdet(Q) -> tuple[int, float]                        # partial-pivot LU; sign in {-1,+1}; raises on singular
class MissingRootError(RuntimeError): ...
def refine_root(evaluate, sector: Sector, c, E_left: float, E_right: float, xtol: float = 1e-11, max_iter: int = 200) -> float
    # evaluate = qc3.oracle.evaluate_point or qc3.fast.evaluate_point (called with want=("Q",)).
    # Brent on f(E) = sign * exp(lnabsdet - l0), l0 = max(lnabsdet at the two ends); requires a sign change, else MissingRootError.
```

### 5.6 `qc3.lattice`
```python
def default_data_dir() -> pathlib.Path                      # $QC_DATA_ROOT or <repo>/data
def read_jack(path) -> tuple[dict, np.ndarray]              # header {"N":603,...}, values (N,) float64 (second column)
def jackknife_resample(x: np.ndarray) -> np.ndarray         # y_c = (sum - x_c)/(N-1)
def jackknife_mean_cov(samples: np.ndarray) -> tuple[np.ndarray, np.ndarray]   # samples (N, m) of resampled values -> mean (m,), C (m,m) = (N-1)/N sum (y-ybar)(y-ybar)^T
def elab_to_ecm(E_lab, P: float) -> np.ndarray
@dataclass
class LevelSet: labels: list[tuple[str,int]]; irreps: list[str]; nP: list[tuple]; files: list[str]; Ncfg: int
               Ecm: np.ndarray; cov: np.ndarray; Ecm_cutoff: float; Lbyas: int
def load_levels(data_dir, Lbyas: int, irreps, Ecm_cutoff: float, xi: float = 3.444) -> LevelSet
    # file names 20_<ddd>_<irrep>_n<k>.jack with ddd in {000,001,002,011,111} <-> irrep labels {000_A1m,100_A2,200_A2,110_A2,111_A2};
    # per file: resample, convert each resample to Ecm, keep if mean Ecm <= cutoff; cov over all kept levels (one L).
def match_levels(E_lattice, E_qc) -> list[tuple[int,int]]   # D3; raises MissingRootError if len(E_qc) < len(E_lattice)
def chi2(residual: np.ndarray, cov: np.ndarray, eig_rel_cut: float = 1e-12) -> float
```

### 5.7 `qc3.bench`
```python
def run(irreps=IRREPS, energies=(0.2700, 0.3100, 0.3599), repeats=5, warmup=1, threads=1, backends=("oracle","fast"), c=C_BASELINE, out_json=None) -> dict
# CLI: python -m qc3.bench --out results/bench_<tag>.json [--no-oracle] [--threads 1]
```

Reference package (planner-owned, tests only): `tests/acceptance/m1/qc3ref/` with `kinematics.py`, `fgk.py` (`F_matrix`, `G_matrix`, `K2inv_matrix`), `k3.py` (`K3_basis_matrices`), `basis.py` (`basis_rule_momenta`). Helpers: `tests/acceptance/m1/m1_helpers.py` (cached oracle/fast calls, anchor loading).

## 6. Do-not-retry list

| approach | why not |
|---|---|
| Sphere/cube truncation in $\mathbf n$ with a radius chosen by eye | fails for spectator momenta near the $H$ cutoff ($\gamma$ up to ≈8): measured $3\times10^{-4}$ block errors at $\lvert\mathbf n\rvert^2\le120$. Use the $r^2\le R^2$ rule (§2.3 rung 2) |
| Disk caches of $F_3^{-1}$, hot windows, det grids | D4; the point of M1 |
| Any GPU path (`src/cachegen_gpu`, `src/gpu_det`) | Digonto: no GPU; CPU single-core metric |
| One process per energy (old `generate_eigen_based_F3inv` + Python) | process start-up dominates once a point costs ms |
| Spline/linear interpolation of $F_3^{-1}(E)$ | poles between grid points; replaced by direct evaluation |
| Inverting the projected $F_3$, rescaling by $N_{\rm proj}$, raw $\det$ | invariants 1–2 |
| Tracking "the eigenvalue closest to zero" for roots | branch switching (`docs/ROOT_FINDING.md` §1) — M2 uses $\ln\det$ |
| Changing $\alpha$, $N=20$, $z_{\min}$, $\epsilon_H$ to speed up | invariant 3; they change the oracle |
| Sorted Neumaier summation storing all terms | unnecessary once the kept set is $\sim10^3$ terms: pairwise/Kahan in double meets $10^{-12}$ (test shows it); keep as a debugging option only |
| Treating the `.jack` values as jackknife samples | D1 (factor $602^2$) |

## 7. Benchmark protocol (tester, mode `benchmark`; `reports/benchmarks/NNN-*.md` + `results/bench_*.json`)

- Machine: `lscpu` model, cores, `git rev-parse --short HEAD`, compiler + flags, NumPy/SciPy versions, date.
- Pinning: `OMP_NUM_THREADS=1 taskset -c 0` for the metric; a second run with `OMP_NUM_THREADS=4` reported separately (not the metric).
- Points: the 5 irreps × $E^*\in\{0.2700,0.3100,0.3599\}$ (low $N$, mid, largest $N$); plus the 12 anchor points. $\vec c=\vec c_{\rm base}$. 1 warm-up, 5 repeats, report median and min per stage (`timings`) and per point, for oracle and fast; speed-up = median ratio. The oracle is run once per point (it is slow) unless the tester has time.
- χ² proxy: total time of `refine_root` over the 12 L20 windows at cutoff 0.335 (from `data/windows/L20/*/accepted_windows.csv`, columns `E_left_bracket`, `E_right_bracket`), fast backend, 1 core; number of evaluations per root.
- Memory: peak RSS (`/usr/bin/time -v` if available, else `resource.getrusage`).
- Regression threshold: a slowdown > 10 % at any point versus the previous benchmark report is flagged.

## 8. Risks and stop conditions

- **pybind11/Eigen build problems** in this container: fall back to a C ABI + `ctypes` with caller-allocated buffers (same Python interface). Builder may spend at most 2 attempts before reporting.
- **Dimension jumps:** at energies where $N$ changes, oracle and fast must agree on $N$ (same basis code). The random-energy test avoids $\pm5\times10^{-5}$ of a jump; exact-jump behaviour is M2's concern.
- **Ill-conditioned points** ($\kappa(F_3)\gtrsim10^{11}$) could exceed the $10^{-6}$ $\ln\det$ tolerance while the roots still agree. If that happens, report the point and $\kappa$; the planner decides (do not loosen).
- **Oracle non-convergence:** the oracle's cube is not converged for high-$\gamma$ momenta (research 007 §3.5), but those rows carry $H\lesssim10^{-8}$. The fast path reproduces the oracle by construction (subset of the same cube). If a later physics check shows this matters, it is a change of invariant 3 → Digonto.
- **Stop and ask Digonto:** D6 (p-wave form, after step 9); any need to change a production constant; if oracle and reference disagree beyond tolerance at any point (physics error in the core); if the fast path cannot reach 20× (method failure → planner, then Digonto if a second plan fails).

## 9. Outline of later milestones (not planned in detail here)

- **M2 — roots without labels, smooth root position, fitter:** $g(E)=\mathrm{Tr}[Q^{-1}Q']$ jump classifier scored on the 208 labelled brackets (FP = FN = 0 gate, `docs/ROOT_FINDING.md` §3.3); Newton/secant root refinement inside the FCN using the fast path; D3 matching; `iminuit` (Minuit2) fit of $\vec c$ with the corrected covariance; L20 cutoff scan 0.315–0.345; v32w cross-check ($\chi^2=11.4897$, 12 levels) with the same reader. Minuit2 via `iminuit` 2.33 (pip) — no C++ Minuit2 build.
- **M3 — L24 and combined fits** (needs the L24 `.jack` files from Digonto).
- **M4 — physics output:** uncertainties incl. two-body inputs and $\xi$, spectrum plots, optional connection to $\mathcal M_3$.

## 10. Expected cost

- Builder: step 1 ≈ 1 h (build + stubs); step 2 ≈ 2–3 h (bindings of existing code); steps 3–5 ≈ 1 h each; step 6 ≈ 4–6 h (the core work); step 7 ≈ 4 h; steps 8–9 ≈ 1 h each.
- Tests: fast m1 suite ≤ 2 min (each oracle call 3–6 s, cached per module; ~10 oracle calls); slow suite ≈ 10–15 min (12 anchors + 10 random energies live oracle ≈ 2 min, two oracle root refinements ≈ 2×40 s, speed-up test ≈ 30 s, large reference points ≈ 2 min of Python).
- Memory: < 1 GB everywhere ($N\le200$; class tables ≤ tens of MB per sector).
- Fast path targets: 5–50 ms per point single core after step 7; χ² proxy (12 roots × ≈6 evaluations) ≈ 0.5–3 s single core without any cache (old direct FCN: 1.4 s with a full-grid cache; hot-window FCN: 24 ms with 150 MB caches per sector).

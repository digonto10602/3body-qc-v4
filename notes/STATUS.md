# STATUS

- **Updated:** 2026-10-21 (cloud kickoff)
- **Prompt in force:** `prompts/002-m1-step7-speed-replan-projector-basis-k3-g.md`. It covers step 7 (sub-steps 7.0–7.8). Steps 8, 9 and the milestone gate stay as in `prompts/001`.
- **Acceptance tag:** `acceptance-v4` = `504a689` (47 tests). v4 = v3 plus the 13 step-7 tests (`test_m1_fast_v2.py`, `test_m1_roots.py`), `test_m1_production_default.py`, and the D7 pins of `coded_v2` in the locked comparisons. Earlier tags: `acceptance-v1` = `7801d97`, `-v2` = `bcf12cd`, `-v3` = `e74f94b`. Tags exist **only locally** (tag push gives HTTP 403); recreate them with `git tag -a <tag> <commit> -m ...` if they are missing. Check with `.claude/bin/acceptance-check`.
- **Current milestone:** M1, a cache-free single-energy QC (plus the D1 covariance fix and the D3 level matching).
- **Current step:** M1 step 7.7 (optional $F$ work, in a builder worktree). After it: the milestone gate.
- **Step 8 DONE**, verified by `reports/tests/011-m1-step8-bench-harness.md` (PASS). Code: `qc3.bench` and `python -m qc3.bench`.
  - Smoke-run speed-up vs the oracle: 367–2019× (median 716) including the oracle's projector build, or 348–854× without it.
  - **χ² proxy** (12 labelled L20 roots, `coded_v2`, 1 core, no caches): warm 0.49–0.53 s, cold 0.57–0.63 s.
    - That is 102 evaluations of about 5 ms each.
    - Compare: 1.4 s for the old direct FCN (2.7× slower) and 24 ms for the hot-window FCN (21× faster, but it needed GBs of caches).
    - Next lever, for M2: fewer evaluations per root (warm start from the previous fit iterate, Newton on $\ln|\det Q|$).
  - χ² at the fixed baseline $\vec c$ = 53.890 (12 levels, D1 covariance). This is **not a fit**. An independent recomputation agrees to $1.3\times10^{-12}$.
  - Minor: record the compiler and flags in the bench JSON; use minima for timing comparisons (noise is up to 2×).
- **Steps 7.4–7.6 DONE**, verified by `reports/tests/010-m1-step7-4-5-6-kernels-tail.md` (PASS). **All 47 acceptance tests are green.** Merge commit `3e26d29`.
  - $\mathcal K_{\rm df,3}$ agrees with the oracle to $\le8.9\times10^{-17}$ over 47 energies, including energies at dimension jumps.
  - $G$ is bit-identical to the oracle. The default fills every ordered pair (see the decision below).
  - The tail computes $Q$ by an LU solve; it equals the projection of the explicit inverse to $1.7\times10^{-15}$.
  - Fast vs oracle at the 12 anchors: $\ln|\det Q|$ differs by $\le3.6\times10^{-9}$.
  - Warm time per point: 3.3–13.8 ms. At $N=160$ the speed-up is **581×** single-core, or 554× without counting the oracle's projector build. Cold minus warm is $\le25$ ms.
  - Truly single-threaded: CPU/wall is 0.97–1.00 with BLAS not pinned.
  - **Coordinator decision:** the $G$ fill mode defaults to every ordered pair. That deviates from `prompts/002` §2.3, costs +0.14 ms, and gives bit-equality with the oracle. The symmetric fill is still available as `g_symmetric=True`. The planner is to record this in 002 and METHODS §15 at the gate.
  - **Test-defect candidate (ask Digonto at the gate):** `test_fast_single_thread_cpu_time` fails when run alone, because of OpenBLAS's 0.1–0.3 s thread spin after `import numpy` (verified with a bare probe). It is order-dependent and can't detect real threading at its small point.
  - **N2:** at 111_A2 $E^*=0.3382011860823352$ ($\kappa(Q)=1.5\times10^{13}$), fast vs oracle is $3.06\times10^{-6}$. A 60-digit mpmath check shows the fast path is within $1.2\times10^{-12}$ of exact, and the oracle is the inaccurate side. The point is not in the locked tests; it is for the auditor.
- **Steps 7.2 and 7.3 DONE**, verified by `reports/tests/009-m1-step7-2-3-orbit-projector-basis.md` (PASS). 45 of 47 tests pass; the 2 pending need 7.4–7.6.
  - Orbit projector: $N_{\rm proj}$ agrees with an independent dense projector, with eigen and with the oracle at all 58 production bases. $\|VV^\dagger-P\|/\|P\|\le4\times10^{-16}$.
  - Cold build is at most 1.2 ms over the 58 bases (5.4 s before).
  - The fast basis is bit-identical to the oracle's at 1881 energies, including $\pm$1–3 ulp around every jump.
  - `dimension_jumps` matches the oracle's $N(E)$ jumps to $1.1\times10^{-16}$. The counts are 4/11/15/12/11 jumps, i.e. 5/12/16/13/12 bases.
  - Warm time per point: 4.5–12 ms for $N\le64$, 33–68 ms for $N\ge100$.
  - Note for the auditor/planner: the locked `test_fast_single_thread_cpu_time` uses a small point ($N=35$) where OpenBLAS never threads, so it cannot detect the BLAS threading. Unpinned CPU/wall at $N=160$ is 3.06. The 7.6 done-check is to be measured at a large point.
- **Steps 7.0, 7.1 and 7.8 DONE**, verified by `reports/tests/008-m1-step7-0-1-8-cache-roots-default.md` (PASS). 40 of 47 acceptance tests pass; the 7 pending ones are 7.2/7.3 stubs.
  - 7.0: the cache key is complete. A differential test of 130 inputs found 0 stale hits; a C++ harness checks `norm_tol` and `proj_tol`.
  - 7.1: the discriminant gives zeros $-13.1\ldots-26.6$ and poles $+11.9$, $+12.7$. The derived bound is $|\Delta_h|\ge\ln(h/x_{\rm tol})=11.5$.
    - It is reliable for opposite-type pairs separated by more than about $5.5\times10^{-10}$. Closer pairs come out "ambiguous", never the wrong kind.
    - `kind="any"` is identical to the old Brent on 36 real windows.
  - 7.8: the production default is `q0_consistent` (D7). No test compares across forms.
  - Note: `xtol` must stay below $h e^{-4}\approx1.8\times10^{-8}$. Larger values fail loudly as "ambiguous"; add a docstring note.
- **Step 9 DONE**, verified by `reports/tests/007-m1-step9-pwave-shift.md` (PASS). The table is in `results/k2_pwave_shift_L20.json`.
  - `q0_consistent` lowers all 12 levels, by $0.6$–$3.0\times10^{-4}\,a_t^{-1}$, i.e. $0.11$–$0.39\,\sigma$.
  - The oracle confirms all 24 roots as simple zeros, to $\le2\times10^{-12}$.
  - **Finding A:** `qc3.roots.refine_root` accepts poles of $\det Q$. Two cases are reproducible: 110_A2 n=1 and 111_A2 n=0, with `q0_consistent` and $\pm250$-row windows. The fix goes in `prompts/002`.
  - **Finding B:** a hidden even-order zero, $\det Q\propto(E-E^*)^4$, in 110_A2 at $E^*=0.319468828$, for both forms.
    - It sits exactly at the free level $K(\mathbf 0)K(0,1,0)\pi(1,0,0)$ with $\mathbf P=(1,1,0)$.
    - Sign-change searches cannot see it, and the labels do not contain it.
    - It is probably the known "free solution" artefact of a truncated partial-wave QC (Blanton, Romero-López, Sharpe, arXiv:1901.07095). That is a guess, not yet checked.
    - Question for Digonto: should D3 counting include such levels? See D8.
  - The legacy `lattice_Ecm` in `accepted_windows.csv` differs from the D1 mean by up to $1.2\times10^{-4}$ ($0.19\sigma$). The D1 mean is used.
- **Step 6 DONE**, verified by `reports/tests/006-m1-step6-fast-path-v1.md` (PASS). **All 33 m1 acceptance tests are green, slow ones included.**
  - The oracle refactor is bit-identical to `48d1319` at 24 points.
  - Fast vs oracle: $\max|\Delta\ln|\det Q||=5.1\times10^{-10}$ over the 12 anchors (tolerance $10^{-6}$). $F$ agrees to $\le7.3\times10^{-15}$ relative over 35 random points and at boost factors up to $\gamma\approx10$.
  - Single-core speed-up, warm cache: 63–189×, median 82×. Without the oracle's projector build: 62–91×. The acceptance test measured 92.7×.
  - Fast-path time per point, warm: 26–100 ms. At $N=160$ the basis builder takes 21 ms, the $\mathcal K_{\rm df,3}$ basis 43 ms, $G$ 20 ms and $F$ 10.6 ms.
  - **Cold projector build:** up to 5.4 s (000_A1m 0.34538, $N=103$). Fast cold vs the oracle is then only 1.9×.
  - Tester notes:
    1. The NumPy tail uses OpenBLAS threads unless `OPENBLAS_NUM_THREADS=1`; the benchmark must pin this.
    2. **The projector cache key omits `eig_tol`, `norm_tol` and `proj_tol`.** A stale hit was demonstrated. Production is unaffected, because `eig_tol` is fixed, but it must be fixed.
    3. Rung 4 is only partly done: $G$ and $\mathcal K_{\rm df,3}$ are not filled symmetrically, and the projection uses an explicit inverse.
- **Step 5 DONE**, verified by `reports/tests/005-m1-step5-portability.md` (PASS).
  - `default_data_dir` honours `QC_DATA_ROOT`, otherwise `<repo>/data`, independent of the working directory.
  - The three silent-failure paths now raise.
  - The real L20 numbers are bit-identical to report 004.
  - 147 legacy files outside the M1 path still contain `/home/digonto`-style paths. They are not touched; the counts per directory are in report 005.
  - `accepted_windows.csv` keeps such paths in its provenance column only. That is acceptable (`prompts/001` §4, amendment v3).
  - Minor note: `chi2` on an empty covariance raises with an unclear NumPy message.
- **Step 4 DONE**, verified by `reports/tests/004-m1-step4-qc-lattice.md` (PASS).
  - 12 of 12 step-4 tests and 16 unit tests pass.
  - D1 is confirmed. The legacy reader's variance is $602^2$ times too large: $\sigma=0.3102$ against the correct $5.153\times10^{-4}$ for 20_000_A1m_n0.
  - The 12-level covariance at cutoff 0.335 has condition number $1.8\times10^{3}$. Its correlations run from 0.49 to 0.95.
  - `refine_root` on the oracle converges in 5 calls (13.8 s) at 100_A2 level 0.
- **Step 3 DONE**, verified by `reports/tests/003-m1-step3-reference-pwave.md` (PASS).
  - The `q0_consistent` form is added as a new function. The existing physics code is byte-unchanged.
  - All 5 reference tests and all 7 oracle tests pass.
- **Step 2 DONE**, verified by `reports/tests/002-m1-step2-oracle-wrapper.md` (PASS).
  - All 7 tests in `test_m1_oracle.py` pass.
  - 3 tests in `test_m1_reference.py` already pass.
  - The oracle takes 2.0–10.4 s per point on 1 thread; $F$ is 87–98% of that.
  - The final $Q$ and slogdet are computed in NumPy, so they match the anchors bit for bit. The Eigen LU values are kept in the metadata.
- **Step 1 DONE**, verified by `reports/tests/001` (PASS).
- **Last completed:** kickoff.
  - Planner report: `reports/research/007-single-energy-qc-methods-and-reference.md`.
  - Independent reference: `tests/acceptance/m1/qc3ref/`, NumPy code written from arXiv:2111.12734. It agrees with the oracle dumps to $\le 2.5\times10^{-14}$ (relative).
  - The oracle (existing C++ path, 1 core) takes 2.8–6.5 s per energy point.
- **Plan note:** step 1's last column says `test_no_absolute_home_paths_in_m1_code` "should pass immediately". It cannot: the test calls `qc3.lattice.default_data_dir()`, which step 5 implements. The step table assigns both portability tests to step 5, so step 1 does not have to turn it green. The test itself is fine.

## Decisions from Digonto (answered 2026-10-21)

| id | question | answer | consequence |
|---|---|---|---|
| D1 | What is in the `.jack` files? | "The jack files contain masses for a particular state for configN = (553 or 603 depending on L24 or L20 …) gauge configurations." | They are **per-configuration** values ($N=603$ for L20, $N=553$ for L24), not jackknife samples. `OPEN_ISSUES` #1 is **confirmed**: the multi-L reader inflates the covariance by $(N-1)^2$ and must be fixed (M1 step 4). Every legacy L20 fit is invalid. |
| D2 | Where do the two-body inputs $(4.04,-43.2,4.12)$ come from? | "From my earlier work and from Hadron Spectrum Collaboration papers." | Keep the values. Still open (planner research, not Digonto): cite the exact papers, and check that the coded p-wave form $-1/(a_1q^{*2})$ is the parametrisation in which $a_1=-43.2$ was quoted (`OPEN_ISSUES` #2). |
| D3 | Level matching when there are more QC levels than lattice levels | "$m$ lattice levels and $n\ge m$ QC levels: sort both and use only the lowest ones up to the lattice count." | Rule: sort both ascending and pair them in order; use the lowest $m$ QC levels; ignore the extra QC levels. Unequal counts are no longer an error when $n>m$. ($n<m$ still is: a missing root.) |
| D4 | Missing caches, det grids and hot-window `.bin` files | Not provided. **New direction:** remove the need for cache generation and re-reading altogether. | See below and `prompts/000` "Update 2026-10-21". |
| D5 | Which system? | "KKpi (I=2)", i.e. $K^+K^+\pi^+$. | The 2+1 identical-pair formalism is valid. Uncertainties on $\xi$, $a_tm_K$, $a_tm_\pi$ and the two-body inputs are still not given (M5). |

### Direction for the cloud session (D4, Digonto's words, paraphrased)
- Make the code fast enough that cache generation and saving are not needed.
- Run a loop that takes a **single energy point** and finds what makes solving the QC (or one fit evaluation) take the least time, "using graphical manipulation", so that each point runs fast without caches.

**D6 (answered 2026-10-21, cloud kickoff):** "graphical manipulation" means **symbolic manipulation**, not GPU.
- No GPU. Goal: the lowest possible time to solve the QC at **one energy point** anywhere in $[E_{\rm cm,min},E_{\rm cm,max}]$, on CPU, without caches.
- The existing CPU/GPU cache generators may be read to understand the computation.
- Timings are compared on **1 CPU core**.
- Loop: think → build optimised code → compare with an **oracle** (the last method that worked and has results, or the existing code run on 1 core for 1 energy point in 1 irrep).
- The solution must generalise to all energy points and all irreps. Literature/web techniques or our own are both allowed.

### Decisions answered at the cloud session (2026-10-22)
- **D7, the p-wave $\mathcal K_2^{-1}$ form:** switch production to **`q0_consistent`** ($-1/a_1+q^{*2}|q^*|(1-H)$, consistent with HadSpec $k^3\cot\delta_1=1/a_1$).
  - `coded_v2` stays selectable for comparison and regression.
  - The locked tests built on `coded_v2` data (the oracle anchors, the labelled `previous_model_zero`) pin `coded_v2` explicitly; the planner is amending them, see `prompts/002`.
  - Evidence: `results/k2_pwave_shift_L20.json`, a shift of $0.11$–$0.39\,\sigma$ per level at the baseline $\vec c$.
- **D8, free-level solutions:** exact free-level solutions of the QC are artefacts and are **excluded** from D3 level matching. Example: the even-order zero at $E^*=0.319468828$ in 110_A2. Detecting them belongs to M2 ("detection beyond sign changes").

## Open failures and attempt counters

Counters were reset at `acceptance-v4` (`prompts/002`).

| check | builder attempts | tester rounds | fixer | planner |
|---|---|---|---|---|
| — | 0 | 0 | 0 | 0 |

## Notes for the M1 auditor (from tester 002)
- **Anchor tolerances are below the round-off floor at ill-conditioned points.**
  - Tolerances: oracle vs anchors $10^{-10}$ (with $\vec c=0$) and $10^{-9}$ (with baseline $\vec c$).
  - At 100_A2 0.2915561 with baseline $\vec c$, $\kappa(Q)=3.8\times10^{12}$. At 200_A2 0.34738, $\kappa(F_3)=4.8\times10^{8}$.
  - These tests pass only because the oracle uses the same LAPACK `slogdet` as the anchor script. The margin is 7.3×; it is robust to thread count.
  - So they test reproducibility, not physical accuracy. An independent `eigvalsh` determinant differs by $4\times10^{-7}$ at the worst point.
  - They may fail on a different CPU or BLAS. Changing them is the planner's decision; it is not a blocker now.
- **A tautological assertion.** The last assertion of `test_oracle_reality_hermiticity_isometry` compares NumPy `slogdet` with itself.
- **Commit hash can go stale.** `metadata["git_commit"]` is fixed when `_core` is built. Rebuild `_core` after committing before producing any result that leaves the project.
- **Projector build counted in the oracle time.** The oracle rebuilds the projector on every call: 5.3 s at 000_A1m 0.34538. The speed-up benchmark must report the gain stage by stage.

## Known-unverified

- `build.sh` and `build_gpu.sh` were never run in this layout. Minuit2 (C++) is absent in the cloud container; M2 will use `iminuit` (pip).
- Results in `results/legacy_*` are reference outputs only. Their χ² and parameters are invalid (D1 confirmed issue #1).
- The L24 `.jack` files are not in the repo.
- The planner's prototype and validation scripts (`scratch/planner/*.py`) and its compiled oracle binary are not in git (`scratch/` is gitignored). They are lost when the container ends; report 007 names them.

## Next action

1. Builder: 7.0 (cache key), 7.1 (pole/zero discrimination), 7.8 (D7 default), per `prompts/002` §3. Then the tester.
2. 7.2 orbit projector → 7.3 threshold basis → 7.4 $\mathcal K_{\rm df,3}$ pair kernel → 7.5 $G$ pair kernel → 7.6 C++ tail (7.7 optional).
3. Step 8 (benchmark + χ² proxy). Use the classifying `refine_root`, and pin `coded_v2` for the 12 labelled roots (`prompts/002` §2.5). Pin `OPENBLAS_NUM_THREADS=1`.
4. Milestone gate: tester full, auditor, tester benchmark, tag `milestone-1`.
5. M2 note: with `q0_consistent`, 5 of the 12 labelled windows have no end-point sign change. So M2's own root finder (FP = FN = 0, scored with `coded_v2`) is on the critical path before any `q0_consistent` fit.

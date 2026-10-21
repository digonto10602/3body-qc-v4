# STATUS

- **Updated:** 2026-10-21 (cloud kickoff)
- **Prompt in force:** `prompts/001-m1-cache-free-single-energy-qc.md` (plan §3 steps, §4 tolerances, §5 interface, §7 benchmark protocol).
- **Acceptance tag:** `acceptance-v3` (`tests/acceptance/m1/`, 33 tests, 6 marked `slow`; v2 = v1 plus the sign-flip fix; v3 = v2 plus the portability-test self-scan fix; both in the `prompts/001` §4 amendments, approved by Digonto). Check with `.claude/bin/acceptance-check`. Tags exist **only locally**, because this cloud session cannot push tags (HTTP 403). `acceptance-v1` = `7801d97` , `acceptance-v2` = `bcf12cd` and `acceptance-v3` = `e74f94b`. If they are missing, recreate them with `git tag -a <tag> <commit> -m ...`.
- **Current milestone:** M1, a cache-free single-energy QC (plus the D1 covariance fix and the D3 level matching).
- **Current step:** M1 **step 5**, portability (`default_data_dir`, `QC_DATA_ROOT`, no absolute home paths).
  - It also fixes three silent-failure paths in `lattice.py` that tester 004 found:
    1. `load_levels` drops levels whose $E^*$ is NaN without telling anyone;
    2. `chi2` ignores negative eigenvalues and an all-zero covariance;
    3. file names with a bad `_n<k>` suffix are skipped without a message.
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

### Pending decision for Digonto
- **D7, the p-wave $\mathcal K_2^{-1}$ form.** The planner calls this "D6" in `prompts/001` §2.2 and step 9; it is the same question.
  - HadSpec, arXiv:1411.2004 Sec. VI, quotes $a_1$ with $k^{3}\cot\delta_1=1/a_1$.
  - In the code's $Q_0$ convention, the consistent element is $-1/a_1+q^{*2}|q^*|(1-H)$. The coded form is larger by a factor $1/q^{*2}$.
  - The production default stays `coded_v2` until Digonto decides.
  - Step 9 produces the table of level shifts for the 12 levels, as input to that decision. The planner recommends `q0_consistent`.

## Open failures and attempt counters

M1 counters (reset at kickoff, `acceptance-v1`):

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

2. The tester (mode `step`) checks step 4. Then commit.
3. Steps 5–9 as in `prompts/001` §3.
4. Notes for later steps:
   - Step 8: `qc3.bench.main()` must parse `--out`, `--no-oracle` and `--threads`.

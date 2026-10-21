# STATUS

- **Updated:** 2026-10-21 (cloud kickoff)
- **Prompt in force:** `prompts/001-m1-cache-free-single-energy-qc.md` (plan §3 steps, §4 tolerances, §5 interface, §7 benchmark protocol).
- **Acceptance tag:** `acceptance-v1` (`tests/acceptance/m1/`, 33 tests, 6 marked `slow`). Check with `.claude/bin/acceptance-check`. The tag points at commit `7801d97`. It exists **only locally**: this cloud session cannot push tags (HTTP 403). If it is missing, recreate it with `git tag -a acceptance-v1 7801d97 -m "kickoff: prompts/001-m1-cache-free-single-energy-qc.md"`.
- **Current milestone:** M1, a cache-free single-energy QC (plus the D1 covariance fix and the D3 level matching).
- **Current step:** M1 **step 1**, "Build + stubs".
  - Done: the Python stubs in `src/python/qc3/`. All 33 acceptance tests fail with `NotImplementedError` (checked by the coordinator).
  - Remaining: the `./build.sh core` target (pybind11 module `qc3._core`, an empty module is fine), and `all` no longer requiring Minuit2.
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

## Known-unverified

- `build.sh` and `build_gpu.sh` were never run in this layout. Minuit2 (C++) is absent in the cloud container; M2 will use `iminuit` (pip).
- Results in `results/legacy_*` are reference outputs only. Their χ² and parameters are invalid (D1 confirmed issue #1).
- The L24 `.jack` files are not in the repo.
- The planner's prototype and validation scripts (`scratch/planner/*.py`) and its compiled oracle binary are not in git (`scratch/` is gitignored). They are lost when the container ends; report 007 names them.

## Next action

1. M1 step 1 (rest): builder adds `./build.sh core` (pybind11) and makes `all` Minuit2-optional. Then the tester (mode `step`) runs the step-1 done-checks. Then commit.
2. Step 2: oracle wrapper (`qc3.oracle.evaluate_point`, `src/qc3_core/oracle_point.hpp`).
3. Steps 3–9 as in `prompts/001` §3.

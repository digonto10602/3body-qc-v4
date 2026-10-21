# STATUS

- **Updated:** 2026-10-21
- **Prompt in force:** `prompts/000-consolidation-and-roadmap.md`. Read its "Update 2026-10-21" section first.
- **Current milestone:** M1, the trustworthy and cache-free baseline. Not started.
  - The decisions are answered (below).
  - The next run is the planner kickoff, in a Claude cloud session.
- **Last completed:** M0, consolidation.
  - QC3_release was removed from the repo on 2026-10-21 (not Digonto's code); it is only cited.
  - The physics-loop agents, hooks and rules were copied into `.claude/` so that cloud sessions have them.

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

## Open failures and attempt counters

None yet: no build/test loop has run. Counters start at M1 step 1.

| check | builder attempts | tester rounds | fixer | planner |
|---|---|---|---|---|
| — | 0 | 0 | 0 | 0 |

## Known-unverified

- `build.sh` and `build_gpu.sh` were never run in this layout. Eigen and Minuit2 were absent on the consolidation machine.
- Results in `results/legacy_*` are reference outputs only. Their χ² and parameters are invalid (D1 confirmed issue #1).
- The L24 `.jack` files are not in the repo.

## Next action

1. In the cloud session, run the `planner` (kickoff) on `prompts/000` including its "Update 2026-10-21" section. It writes `prompts/001-*` and `tests/acceptance/m1/`; then run `.claude/bin/acceptance-lock`.
2. Builder loop: M1 step 1 (build), step 2 (portability), step 4 (covariance fix per D1), step 5 (C++ vs an independent reference).
3. Start the single-energy-point timing loop (new step 6).

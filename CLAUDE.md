# CLAUDE.md — project rules for 3body-qc-v4

Digonto's global rules apply: report format, the physics loop, numbering, and the hard correctness rules. `.claude/global-CLAUDE.md` is a **verbatim copy of Digonto's `~/.claude/CLAUDE.md`** (2026-10-21). It is imported here so it also loads in cloud sessions, where `~/.claude/` does not exist:

@.claude/global-CLAUDE.md

**Path note for this repo:** wherever those rules say `~/.claude/bin/<tool>`, use `.claude/bin/<tool>` (run from the repo root) if `~/.claude/bin/` does not exist. The graphify skill they mention (`~/.claude/skills/graphify`) is not in the repo; the `graphify` CLI is (see setup below).

The physics loop is **active**, because `notes/STATUS.md` exists.
- The agents (`builder`, `tester`, `fixer`, `auditor`, `planner`), the guard and session hooks, `/kickoff`, and the helpers `next-num`, `acceptance-lock` and `acceptance-check` are all in `.claude/`.
- Call the helpers as `.claude/bin/<name>` from the repo root.
- The project hooks run only when `~/.claude/hooks/` is absent, so they never run twice on Digonto's machine.

This file adds project-specific facts.

## Setting up a fresh (cloud) checkout
```bash
pip install graphifyy==0.9.53                  # the `graphify` CLI used on the consolidation machine
graphify update .                              # rebuild graphify-out/ for this checkout
# C++: Eigen 3 (header-only, e.g. apt install libeigen3-dev; export EIGEN_INC=/usr/include/eigen3)
#      Minuit2 (standalone or from ROOT; export MINUIT2_ROOT=...) for the fitter only
```
- Python test deps: `pip install numpy scipy pytest`.
- **Build:** `./build.sh all` (CPU tools; fitter needs Minuit2).
- **Fast check (after every edit):** `python3 -m pytest -q tests --ignore=tests/acceptance` — never includes `tests/acceptance`.
- **Acceptance (milestone N):** `python3 -m pytest -q tests/acceptance/mN`, then `.claude/bin/acceptance-check`.
- Cloud sessions normally have **no GPU**. `build_gpu.sh`, `src/cachegen_gpu/` and `src/gpu_det/` cannot be compiled or timed there; CPU work can.
- QC3_release (Blanton, Romero-López, Sharpe; github `ferolo2/QC3_release`) is **not** in this repo: it is not Digonto's code. It is only cited as the published reference implementation. Do not add it back.
- Data in the repo: the L20 `.jack` spectra only. No cache binaries exist anywhere. The L24 `.jack` files are not here.

## Orientation
- Read `notes/STATUS.md` first, then the highest-numbered `prompts/NNN-*.md` (the plan in force), then `docs/OPEN_ISSUES.md`.
- Equations and conventions: `docs/METHODS.md`. Data flow and file formats: `docs/PIPELINE.md`. Root finding: `docs/ROOT_FINDING.md`.
- If `graphify-out/graph.json` exists, use `graphify query/explain/affected` before reading many files. Run `graphify update .` after code changes.
- Every file copied from the old packages is listed with its md5 and original path in `docs/PROVENANCE.tsv`. Do not edit `docs/original_handoffs/` (verbatim history).

## Physics invariants (never change silently)
1. $Q=(L\xi)^{-6}\,V_{\rm sel}^\dagger(F_3^{-1}+\mathcal K_{\rm df,3})V_{\rm sel}$. Project the inverse; scale the matrix, not the determinant; no $N_{\rm proj}$ rescaling.
2. Determinants via LU as (sign, log|det|). Never compare determinants across an $N_{\rm proj}$ change.
3. Production physics constants: $\xi=3.444$, $a_tm_K=0.09698$, $a_tm_\pi=0.06906$, $\eta=(1,\tfrac12)$, $\alpha=0.5$, $\epsilon_H=0$, `max_shell_num`=20, waves $\{0,1\}$ / $\{0\}$, parity $-1$, `eig_tol`=0.05, scatter $(4.04,-43.2,4.12)$. Changing any of them invalidates every cache. Record them in every result's metadata.
4. Grid: $E_i=0.26310+i\cdot0.0969/19999$, $i<20000$. The cache sidecar is the source of truth.
5. Raw-cache complex order is natural (Re, Im) on disk. Any new reader must pass the Hermiticity test ($\|F_3^{-1}-F_3^{-1\dagger}\|/\|F_3^{-1}\|\lesssim10^{-14}$, $\max|\mathrm{Im}\det Q|\approx0$).
6. No root-finding fallback to full-cache reads inside the FCN. Caches load once, outside the FCN.
7. Human-labelled roots are the current ground truth. A new root finder replaces them only after FP = FN = 0 on all labelled sectors (`docs/ROOT_FINDING.md` §3.3).

## Decisions already made by Digonto (2026-10-21; details in `notes/STATUS.md`)
- **D1:** `.jack` files = one value per gauge configuration ($N=603$ L20, $N=553$ L24); build jackknife averages before the covariance.
- **D2:** two-body inputs $(4.04,-43.2,4.12)$ come from Digonto's earlier work and HadSpec papers.
- **D3:** level matching = sort both lists ascending and use the lowest $m$ QC levels for $m$ lattice levels ($n\ge m$ allowed).
- **D4:** cache binaries will not be provided. The priority is a cache-free fast QC (`prompts/000`, "Update 2026-10-21").
- **D5:** the system is KKπ with $I=2$ ($K^+K^+\pi^+$).

## Things that are known to be wrong or unverified (see `docs/OPEN_ISSUES.md`)
- **Do not quote any χ² or parameter from `results/legacy_*`.** The covariance is wrong by $(N-1)^2$ (confirmed by D1).
- The p-wave $\mathcal K_2^{-1}$ form differs from QC3_release. Which form matches the HadSpec definition of $a_1$ is still to be checked; do not change it silently.
- `compare_gpu_cache_vs_cpu_openmp_det_oldscale.cpp` reads $F_3^{-1}$ with a single swap. Do not use its read path.
- `src/gpu_det/` is not validated against the CPU.

## Pushing to GitHub (pre-authorised by Digonto, 2026-10-21)
- Cloud and local sessions working on this repo may `git push` to `origin` (github.com/digonto10602/3body-qc-v4) **without asking**. This overrides the "push" item in the "Stop and ask me" list of `.claude/global-CLAUDE.md`, for this repo only.
- Force-push is allowed only as `git push --force-with-lease=<branch>:<expected-old-sha> origin <branch>`, so that a remote you have not seen is never overwritten. Never use plain `--force` / `-f`.
- Still ask Digonto before rewriting history that is already on GitHub (rebase, amend or squash of pushed commits), deleting branches or tags, or deleting data.
- The coordinator pushes after every commit that passes the tester, so that work survives the cloud session ending.

## Stop and ask Digonto
- Changing any decision D1–D5, or any production constant in invariant 3.
- What "graphical manipulation" means for the speed loop (GPU or something else), if the plan depends on it.
- Any need for cluster time or for data not in the repo (L24 spectra).

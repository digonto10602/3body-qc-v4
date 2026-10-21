# CLAUDE.md — project rules for 3body-qc-v4

The global rules in `~/.claude/CLAUDE.md` apply: the physics loop, report format, numbering, and the hard correctness rules. The physics loop is **active**, because `notes/STATUS.md` exists. This file adds project-specific facts.

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

## Things that are known to be wrong or unverified (see `docs/OPEN_ISSUES.md`)
- **Do not quote any χ² or parameter from `results/legacy_*`.** The covariance convention is suspected wrong by $(N-1)^2$.
- The p-wave $\mathcal K_2^{-1}$ convention differs from QC3_release.
- `compare_gpu_cache_vs_cpu_openmp_det_oldscale.cpp` reads $F_3^{-1}$ with a single swap. Do not use its read path.
- `src/gpu_det/` is not validated against the CPU.

## Stop and ask Digonto
- Anything that touches the decisions D1–D5 listed in `notes/STATUS.md`.
- Any need for the external caches on `/media/digonto/Data` or for cluster time.

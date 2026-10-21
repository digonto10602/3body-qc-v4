# 000 — Consolidation and roadmap (starter plan from chat, 2026-10-21)

**Prompt in force:** this file, until the planner writes `prompts/001-*`.
**Origin:** Digonto's request in chat, paraphrased:
- Consolidate all three-body QC work into one repo.
- The last effort was a faster QC solve that caches the expensive pieces so that $\mathcal K_{\rm df,3}$ fits run quickly.
- Roots are found in several ways; "gradient change" might be the best way.
- **The ultimate goal is to fit three-body K matrices to lattice data.**

## Update 2026-10-21 — Digonto's answers and new priority (read first)

**Decisions D1–D5.** The full text is in `notes/STATUS.md`.
- **D1:** the `.jack` files hold per-configuration values ($N=603$ for L20, $N=553$ for L24). Issue #1 is confirmed.
- **D2:** the two-body inputs come from earlier work and HadSpec papers.
- **D3:** sort both level lists and use only the lowest $m$ QC levels when there are $n\ge m$ of them.
- **D4:** the cache binaries will not be provided.
- **D5:** the system is KKπ with $I=2$.

**New priority (replaces M1 step 6, "Restore inputs"): a cache-free QC.**
- Digonto wants the code fast enough that cache generation and saving (raw $F_3^{-1}$ caches, hot-window binaries, det grids) is no longer needed.
- Method: a loop that takes **one energy point** and measures what makes solving the QC, or one χ² evaluation, cheapest. Each point should then run fast without caches. Digonto said "using graphical manipulation"; the meaning is to be confirmed at kickoff (GPU?). Cloud sessions normally have no GPU.

Constraints for the planner when it refines this into `prompts/001-*`:
1. **Reference value first.** At a fixed point (suggested: L20, irrep 100_A2, $\vec c$ = the v32w best fit, energies from `data/windows/L20`), compute $(\mathrm{sign},\ln|\det Q|)$ in float64 with the existing CPU path. Check it against an independent reference: QC3_release is **not** in this repo and must not be added (not Digonto's code). The planner chooses the reference. Every faster variant must reproduce this reference within a tolerance the planner fixes and justifies in writing.
2. **Profile before optimising.** Time each piece separately: $F$ lattice sums, $G$, $\mathcal K_2^{-1}$, the $F_3$ inversion, the projection with $V_{\rm sel}$, the LU. Report the times in `reports/benchmarks/`.
3. **No disk caches.** Reusing energy-independent pieces in memory within one process is allowed (momentum lists, projectors, $Y_{\ell m}$). Each candidate speed-up must be measured against the reference, never assumed. The invariants in `CLAUDE.md` still hold, in particular "project the inverse".
4. **Report the end-to-end number** that matters: wall time per χ² evaluation for the L20 cutoff-0.335 fit (12 levels). Compare with 1.4 s for the old direct FCN and 24 ms for the hot-window FCN (package logs).

Steps 1–5 below still apply. Step 4 (covariance) is now a confirmed bug fix. Step 3 (raw-cache reader) is low priority, because the caches are not available.

## Goal

A reproducible, fast and correct fit of $\vec c=(K^{\rm iso}_0,K^{\rm iso}_1,K_B,K_E)$ to the KKπ lattice spectrum:
- L20 first, then L20+L24;
- with defensible root identification without human labels;
- with uncertainties.

## Milestone 0 — consolidation (DONE, 2026-10-21)

The repo was assembled from the newest copies, with provenance (`docs/PROVENANCE.tsv`), docs, and six research reports (`reports/research/001`–`006`). Issues found are in `docs/OPEN_ISSUES.md`.

## Milestone 1 — trustworthy baseline (decisions answered 2026-10-21; step 6 replaced, see the update above)

Steps for the planner to refine, each with locked acceptance tests:

1. **Build.** `./build.sh all` compiles in the new layout. Fast unit test: the basis sum $\tilde F+\sum c_i\tilde B^{(i)}$ equals direct assembly on a tiny synthetic system.
2. **Portability.** All scripts and configs take repo-relative paths or environment variables (`QC_CACHE_ROOT`, `QC_DATA_ROOT`). No `/home/digonto` or `/media/digonto` literals remain.
3. **Reader.** One plain raw-cache reader (C++ and NumPy) with a round-trip test (write/read a known non-Hermitian complex matrix) and the Hermiticity checks of `docs/OPEN_ISSUES.md` issue 3. Retire the `variant_04` name.
4. **Covariance.** Fix `load_multil_jack_targets` per decision D1. Unit test: a synthetic per-configuration file gives $\sigma=\mathrm{s.d.}/\sqrt N$; a jackknife-sample file gives $\sqrt{(N-1)/N\sum(\cdot)^2}$. Cross-check against `lattice_data_covariance_cpp.hpp` on `data/lattice/L20`.
5. **Reference physics.** C++ vs an independent reference (QC3_release is cited only, not included) for $\mathcal K_{\rm df,3}$ (`tests/reference/`), $F$, $G$ and $\mathcal K_2^{-1}$ at a few energies on a small box. Resolve the p-wave $\mathcal K_2^{-1}$ convention (D2).
6. ~~**Restore inputs** (D4)~~ **replaced** by the cache-free single-energy timing loop (see "Update 2026-10-21").
7. **Re-run the L20 cutoff fits** (0.315–0.345) with the corrected covariance. Report χ²/ν, parameters, Hesse errors and the correlation matrix. Check them against `tools/audit/chi2_covariance_convention_check.py`.
8. **Cross-check against v32w.** Reproduce the v32w fit ($\chi^2/\nu=1.436$ at 0.335, 12 levels) with the new fitter in single-L mode, within the planner's tolerance.

## Milestone 2 — root finding without human labels ("gradient change")

1. **Offline scoring.** Score the $g(E)=\frac{d}{dE}\ln\det Q$ jump-direction classifier (`docs/ROOT_FINDING.md` §3) against all 208 labelled brackets, using determinant grids. Gate: FP = FN = 0 per sector.
2. **Detection beyond sign changes.** Detect even-order roots and roots at segment edges (local minima of $\ln|\det Q|$).
3. **Smooth root position inside the FCN.** Use a Newton step on $\ln\det$ or the bracket matrix pencil, and re-centre windows on the current model root. This replaces the $10^{100}$ penalty with a smooth treatment.
4. **Level matching** per D3: sort both lists and use the lowest $m$ QC levels. The fitter must accept $n>m$. Root detection must be complete below the cutoff, because a missed low root shifts every pair above it.

## Milestone 3 — L24 and combined fit

1. Get the L24 spectra. Build the L24 det grids and labels, resolving 100_A2 bracket 36.
2. Build the L24 hot windows.
3. Run the L20+L24 fit and study the cutoff dependence.

## Milestone 4 — speed and robustness

1. **FCN speed:** precompute $C^+$; OpenMP over windows; no stdout in the FCN.
2. **Cache format:** indexed, energy-ordered, with full physics metadata including parity, `eig_tol`, scattering parameters and git commit. Fail on singular LU.
3. **GPU determinant backend:** reach CPU parity, then use it for det grids.

## Milestone 5 — physics output

1. **Uncertainties:** parameter errors, including two-body parameter, $\xi$ and mass uncertainties (D5).
2. **Spectrum plots:** spectrum vs $L$ with lattice, model and free levels.
3. **Optional:** connect the fitted $\mathcal K_{\rm df,3}$ to infinite-volume amplitudes (`docs/RELATED_CODES.md`).

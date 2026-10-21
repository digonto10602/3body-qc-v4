# 000 — Consolidation and roadmap (starter plan from chat, 2026-10-21)

**Prompt in force:** this file, until the planner writes `prompts/001-*`.
**Origin:** Digonto's request in chat, paraphrased:
- Consolidate all three-body QC work into one repo.
- The last effort was a faster QC solve that caches the expensive pieces so that $\mathcal K_{\rm df,3}$ fits run quickly.
- Roots are found in several ways; "gradient change" might be the best way.
- **The ultimate goal is to fit three-body K matrices to lattice data.**

## Goal

A reproducible, fast and correct fit of $\vec c=(K^{\rm iso}_0,K^{\rm iso}_1,K_B,K_E)$ to the KKπ lattice spectrum:
- L20 first, then L20+L24;
- with defensible root identification without human labels;
- with uncertainties.

## Milestone 0 — consolidation (DONE, 2026-10-21)

The repo was assembled from the newest copies, with provenance (`docs/PROVENANCE.tsv`), docs, and six research reports (`reports/research/001`–`006`). Issues found are in `docs/OPEN_ISSUES.md`.

## Milestone 1 — trustworthy baseline (blocked on D1, D2, D4 in `notes/STATUS.md`)

Steps for the planner to refine, each with locked acceptance tests:

1. **Build.** `./build.sh all` compiles in the new layout. Fast unit test: the basis sum $\tilde F+\sum c_i\tilde B^{(i)}$ equals direct assembly on a tiny synthetic system.
2. **Portability.** All scripts and configs take repo-relative paths or environment variables (`QC_CACHE_ROOT`, `QC_DATA_ROOT`). No `/home/digonto` or `/media/digonto` literals remain.
3. **Reader.** One plain raw-cache reader (C++ and NumPy) with a round-trip test (write/read a known non-Hermitian complex matrix) and the Hermiticity checks of `docs/OPEN_ISSUES.md` issue 3. Retire the `variant_04` name.
4. **Covariance.** Fix `load_multil_jack_targets` per decision D1. Unit test: a synthetic per-configuration file gives $\sigma=\mathrm{s.d.}/\sqrt N$; a jackknife-sample file gives $\sqrt{(N-1)/N\sum(\cdot)^2}$. Cross-check against `lattice_data_covariance_cpp.hpp` on `data/lattice/L20`.
5. **Reference physics.** C++ vs QC3_release agreement for $\mathcal K_{\rm df,3}$ (`tests/reference/`), $F$, $G$ and $\mathcal K_2^{-1}$ at a few energies on a small box. Resolve the p-wave $\mathcal K_2^{-1}$ convention (D2).
6. **Restore inputs** (D4): raw caches, det grids, hot-window binaries for L20. Validate them against the JSON sidecars (SHA-256, row indices).
7. **Re-run the L20 cutoff fits** (0.315–0.345) with the corrected covariance. Report χ²/ν, parameters, Hesse errors and the correlation matrix. Check them against `tools/audit/chi2_covariance_convention_check.py`.
8. **Cross-check against v32w.** Reproduce the v32w fit ($\chi^2/\nu=1.436$ at 0.335, 12 levels) with the new fitter in single-L mode, within the planner's tolerance.

## Milestone 2 — root finding without human labels ("gradient change")

1. **Offline scoring.** Score the $g(E)=\frac{d}{dE}\ln\det Q$ jump-direction classifier (`docs/ROOT_FINDING.md` §3) against all 208 labelled brackets, using determinant grids. Gate: FP = FN = 0 per sector.
2. **Detection beyond sign changes.** Detect even-order roots and roots at segment edges (local minima of $\ln|\det Q|$).
3. **Smooth root position inside the FCN.** Use a Newton step on $\ln\det$ or the bracket matrix pencil, and re-centre windows on the current model root. This replaces the $10^{100}$ penalty with a smooth treatment.
4. **Level matching rule** beyond sorted order (D3).

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

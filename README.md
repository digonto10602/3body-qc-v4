# 3body-qc-v4 — fitting the three-particle K matrix to lattice spectra

This repo fits the three-body short-distance K matrix $\mathcal K_{\rm df,3}$ of the **KKπ** system (two kaons + one pion, 2+1 RFT formalism) to finite-volume lattice energy levels, using the three-particle **quantization condition** (QC):

$$
\det\Big[(L\xi)^{-6}\,V_{\rm sel}^\dagger\big(F_3^{-1}(E)+\mathcal K_{\rm df,3}(E;\vec c)\big)V_{\rm sel}\Big]=0,\qquad
\vec c=(K^{\rm iso}_0,K^{\rm iso}_1,K_B,K_E).
$$

It consolidates all earlier work: eight packages from June–August 2026 and 71 archived versions back to v31. It keeps only the code needed to go forward, plus enough history to see where each piece came from. **Start with [`notes/STATUS.md`](notes/STATUS.md) and [`docs/OPEN_ISSUES.md`](docs/OPEN_ISSUES.md).**

## The idea that makes fits fast

$F_3^{-1}(E)$ is expensive: an $N\times N$ matrix with $N$ up to about 100, built from lattice sums at each of 20,000 energies. $\mathcal K_{\rm df,3}$, however, is *linear* in the four fit parameters. So everything is cached once, offline:

$$
Q_r(\vec c)=\tilde F_r+\sum_{i=1}^4c_i\,\tilde B^{(i)}_r
$$

Only the grid rows near each accepted root ("hot windows") are cached. A χ² evaluation then needs only small determinants, about **24 ms per call instead of about 1.4 s**, and a full Minuit fit takes 2–4 s.

## The hard part: finding the right roots

About 85 % of the determinant's sign changes are poles or artefacts, not energy levels. The current gate is human labelling of every sign change. The next step is a derivative-based classifier ("gradient change") that removes the human from the loop. See [`docs/ROOT_FINDING.md`](docs/ROOT_FINDING.md).

## Layout

| path | content |
|---|---|
| `src/core/` | physics library: $F$, $G$, $\mathcal K_2^{-1}$, $F_3$, $\mathcal K_{\rm df,3}$, irrep projector, Faddeeva |
| `src/cachegen_gpu/`, `src/cachegen_cpu/` | $F_3^{-1}$ + $V_{\rm sel}$ cache generators (CUDA production, CPU reference) |
| `src/cache_io/` | corrected raw-cache reader and trusted CPU determinant scanner (v33h) |
| `src/cache_build/` | hot-window projected-basis cache builder (v33p); full-grid runtime builder (v33g) |
| `src/fitter/` | multi-$L$ Minuit2 fitter with the fast FCN, χ², covariance |
| `src/gpu_det/` | experimental GPU determinant backend (not validated) |
| `scripts/` | workflow: `classifier/`, `windows/`, `fit/`, `audit/`, `plot/`, `cachegen/` (original run scripts) |
| `tools/` | non-interacting levels, state-line interpolation, cache inspectors, audit checks |
| `configs/` | fitter, cache-build and cache-generation configs |
| `data/` | L20 lattice spectra (`.jack`), human root labels, accepted windows, hot-window metadata |
| `results/` | legacy fit outputs (**χ² values invalid**, see OPEN_ISSUES issue 1) |
| `legacy/` | historical code kept for lineage: v32w single-L fitter, root-finding methods |
| `tests/reference/` | C++ vs QC3_release $\mathcal K_{\rm df,3}$ comparison |
| `third_party/QC3_release/` | fetch instructions and local-edit patch for the reference Python code |
| `docs/` | `METHODS.md`, `PIPELINE.md`, `ROOT_FINDING.md`, `OPEN_ISSUES.md`, `HISTORY.md`, `RELATED_CODES.md`, `PROVENANCE.tsv`, `original_handoffs/` |
| `notes/STATUS.md`, `prompts/`, `reports/` | working state of the build/test loop (see `CLAUDE.md`) |

## Build

```bash
export EIGEN_INC=/path/to/eigen3          # required
export MINUIT2_ROOT=$CONDA_PREFIX         # for the fitter (Minuit2 from ROOT or standalone)
./build.sh all                            # fitter, hot-window builder, v33h scan, v33g builder, CPU cachegen
./build_gpu.sh cachegen                   # CUDA cache generator (nvcc, cuBLAS, cuSOLVER)
```

The compile flags are the ones from the original packages; only the include paths changed. **These build scripts have not yet been run in this layout.**

## Running a fit needs external inputs

The raw caches (several GB per sector), the hot-window binaries, the determinant-grid CSVs and the L24 spectra are not in Git. They are listed in [`docs/OPEN_ISSUES.md`](docs/OPEN_ISSUES.md) issue 4.

## Status (2026-10-21)

Consolidation is done. Before any new physics result:
1. Resolve the lattice covariance convention (`OPEN_ISSUES.md` issue 1: errors likely too large by 602).
2. Resolve the p-wave $\mathcal K_2^{-1}$ convention (issue 2).
3. Restore the external inputs.
4. Re-run the L20 fits.

The roadmap is in [`prompts/000-consolidation-and-roadmap.md`](prompts/000-consolidation-and-roadmap.md).

# RELATED CODES (not part of the build)

## QC3_release (reference Python implementation)

- **Upstream:** https://github.com/ferolo2/QC3_release, commit `9bc4b84`. This is the public code of Blanton, Romero-López and Sharpe (arXiv:2111.12734).
- **Role:** test oracle for $F$, $G$, $\mathcal K_2^{-1}$ and $\mathcal K_{\rm df,3}$. The C++ $\mathcal K_{\rm df,3}$ is a port of `base_code/Kdf3/K3main.py` and `K3E.py`.
- **Not included in this repo.** It is not Digonto's code, so it is cited only. It was removed on 2026-10-21. Earlier work used a local working copy based on commit `9bc4b84` with debug edits.

## Three-body full-matrix HPC/CUDA benchmark (2026-08-21)

- **Repository:** `threebody_fullmatrix_hpc_cuda_bench_v3_codex` (Digonto's GitHub).
- **What it solves:** the **infinite-volume** three-body ladder equation $BD=-G$ in momentum representation ($s=E^2$, $a=16$, $m=1$). It is not the finite-volume QC.
- **What it measured:** at $N=6144$, the GPU-native CUDA path takes 9.2 s per energy on an RTX 3070, versus 555 s serial and 426 s with 8-thread OpenMP. That is about $24\times$ like-for-like; the package's own "60×" figure counts CPU-only diagnostics.
- **Known gaps:** CPU/GPU agreement is $3.4\times10^{-9}$ against a required $10^{-11}$ (blocked). The hybrid scheduler was never benchmarked.
- **Reusable ideas for the QC fitter:**
  - build matrices on the device, one thread per element;
  - keep matrices and solver handles resident on the device;
  - plan memory before launch;
  - keep a separate CPU/GPU parity mode.
- **Not transferable:** the QC fit needs many *small* matrices (batched LU), whereas the benchmark does one large LU.
- **Possible future link:** `copied_reference/integralequation_momrep_based.h:368-412` (`F3infvol`) relates a constant $\mathcal K_{\rm df,3}$ to the infinite-volume amplitude. It would matter after a successful fit. Take it from its upstream V14 package, not from that copy.
- **Full review:** [`reports/research/005-cuda-fullmatrix-benchmark-review.md`](../reports/research/005-cuda-fullmatrix-benchmark-review.md).

## External data locations (source machine)

| what | path |
|---|---|
| raw F3inv/Vsel caches | `/media/digonto/Data/F3inv_cache/` |
| lattice spectra (L20, L24) | `/home/digonto/Codes/KKpi_I2/spectrum/Ecm_data/data` |
| original working tree | `/home/digonto/Codes/Practical_Lattice_v2/3body_quantization_v3/` |

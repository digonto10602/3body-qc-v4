# legacy (not built; kept for lineage)

- `v32w_single_L_fitter/`: the last trusted single-L fitter (v32w; χ²/ν = 1.436). Its `K3df_minuit_fit_v32f` header is the **older** version (no precomputed K3 basis), and it uses the per-configuration covariance reader. It is the reference for M1 step 8.
- `root_finding_lineage/`: representative root finders from v7 to v32. They are explained in docs/ROOT_FINDING.md and reports/research/002 §5:
  - `pole_searching.h` (secant);
  - `file_processor.hpp` (slope-ratio thresholds);
  - GPU v9 (inward shape) and v21 (side-|det| trend);
  - CPU v24 (QC with K3, multi-flip);
  - v30s (eigen-min bisection);
  - v31b (closest-eigenvalue crossing direction, the first "gradient-sign" rule);
  - v32 classifier v2 scanner;
  - version-history notes and the K3iso maker.

The full archive of every version is in the original zips (see docs/HISTORY.md), not here.

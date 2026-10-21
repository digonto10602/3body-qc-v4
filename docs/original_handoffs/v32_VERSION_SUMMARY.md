# Version summary and working status

## Best current versions

| Version | Main use | Status |
|---|---|---|
| v32w | K3df fitting with fixed coarse/refined cache separation | trusted fitting baseline |
| v32y | Memory-bounded QC spectrum vs L using fixed K3df | preferred spectrum-vs-L replacement for v32x |
| v32x | First QC spectrum-vs-L fixed K3df | reference only; memory-heavy |
| v32s/v32t | projected F3inv classifier development | useful classifier baseline |
| v32v | verbose fitter/debug version | useful for diagnosing v32w fixed-coarse bug |

## v32w known good fit

```text
chi2 = 11.489748255611568
ndata = 12
npar = 4
ndof = 8
chi2/dof = 1.436218531951446
model_found = 12/12

K3iso0 = 34470.228772793678
K3iso1 = -1035173.9140440958
K3B    = 253921.14757953485
K3E    = -961680.91711685632
```

## Why v32w is trusted

v32v showed that refined points were being merged back into the initial coarse-grid scan, polluting the sign-flip search. v32w separates the fixed coarse cache from the refined cache, so the initial sign-flip detection remains stable while refinement points are still reused inside the refinement stage.

## Why v32y exists

v32x produced one large all-irrep cache per Lbyas. At high Lbyas, this can fill RAM and swap. v32y processes one irrep at a time, saves/closes the file, and releases memory before moving on. It keeps OpenMP parallelism for expensive calculations but avoids holding all irreps in memory simultaneously.

## Available older packages

The folder `old_versions_archive/` contains all available version/package zip files from this conversation. Some earlier v26-v30 versions discussed in conversation were not available as package zips in the sandbox, so this archive contains the available v31/v32 lineage and noninteracting helper packages.

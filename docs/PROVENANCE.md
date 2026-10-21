# PROVENANCE

`PROVENANCE.tsv` lists every file copied from the original packages: its path in this repo, the md5 of the **original** file, and its path inside the unpacked zips (`<package>/<package>/...`). Use it to trace any file to its source version. The zips themselves are listed in `docs/HISTORY.md`.

Files that were changed after copying (the md5 in the TSV is the pre-change original):

| file | change |
|---|---|
| `src/cachegen_gpu/F3_gpu_cachegen_v32zu.cu` | line 473: `#include "../../k3iso01_swave_cpu_cache_qc_omp_package/source/fv_projector_cartesian_l1_v30m.hpp"` → `#include "fv_projector_cartesian_l1_v30m.hpp"` (the header is in `src/core/`). No other change. |
| `tools/audit/chi2_covariance_convention_check.py` | written during consolidation (not from a package); docstring added and paths pointed at `data/lattice/L20` and `results/legacy_L20_passonly_hotwindow` |

Files written new during consolidation (not in the TSV): `README.md`, `CLAUDE.md`, `build.sh`, `build_gpu.sh`, `.gitignore`, `notes/STATUS.md`, `prompts/000-*`, `docs/*.md` (except `original_handoffs/`), the folder READMEs, `.claude/` (the physics-loop agents, hooks, helpers and a verbatim copy of `~/.claude/CLAUDE.md` as `global-CLAUDE.md`, copied from `~/.claude` on 2026-10-21; helper paths in the agents and hooks made repo-relative), `reports/research/*` (investigation reports; a header was added to each).

# results

Only legacy reference outputs so far. **Every χ², parameter value and parameter error here was computed with the multi-L covariance that is suspected to be too large by $(N-1)^2$ (docs/OPEN_ISSUES.md #1). Do not quote them.** They are kept as regression references (root positions, level assignments, file formats).

| folder | what |
|---|---|
| `legacy_L20_passonly_hotwindow/` | L20 pass-only hot-window fits at cutoffs 0.315/0.325/0.335/0.345 (summaries, per-cutoff levels, model true zeros, logs). Produced 2026-07-15 with `scripts/fit/run_l20_full_cutoff_fits_v2.py`. |
| `legacy_L20_k3b_k3e_zero_tests/` | the same with $K_B$ and/or $K_E$ fixed to 0 |
| `legacy_v32w/` | README of the v32w single-L fit (χ²/ν = 1.436, 12 levels, cutoff 0.335; older covariance reader) |

New results go in new folders with commit, command, parameters, seed and machine, as required by the global rules.

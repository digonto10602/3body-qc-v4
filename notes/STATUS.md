# STATUS

- **Updated:** 2026-10-21
- **Prompt in force:** `prompts/000-consolidation-and-roadmap.md`
- **Current milestone:** M1 (trustworthy baseline). Not started; blocked on decisions D1, D2, D4.
- **Last completed:** M0, the consolidation of 8 packages and 71 archived versions into this repo.

## Decisions needed from Digonto

| id | question | why it matters | where |
|---|---|---|---|
| D1 | Do the `data/lattice/L20/*.jack` files hold **per-configuration** energies or **jackknife samples**? | The multi-L fitter assumes jackknife samples. The evidence says per-configuration, in which case all errors are 602× too large and the L20 fits never really moved. | `docs/OPEN_ISSUES.md` #1 |
| D2 | p-wave $\mathcal K_2^{-1}$: keep the coded $-1/(a_1q^{*2})$ form or use the QC3_release form $-1/a_1$? Where do $a_0^{K\pi}=4.04$, $a_1^{K\pi}=-43.2$, $a_0^{KK}=4.12$ come from? | changes the physics of every cache | #2, #7 |
| D3 | Level-matching rule when a sector has more accepted roots than lattice levels (L20/200_A2 at 0.355; 000_A1m at 0.36) | blocks higher cutoffs | #5 |
| D4 | Access to `/media/digonto/Data/F3inv_cache`, the hot-window `.bin` files, the det-grid CSVs, and the L24 `.jack` files | nothing can be re-run without them | #4 |
| D5 | Uncertainties on $\xi$, $a_tm_K$, $a_tm_\pi$ and the two-body parameters; confirm the isospin channel ($K^+K^+\pi^+$, $I=2$?) | error budget; validity of the 2+1 identical-pair formalism | #7 |

## Open failures and attempt counters

None yet: no build/test loop has run. Counters start at M1 step 1.

| check | builder attempts | tester rounds | fixer | planner |
|---|---|---|---|---|
| — | 0 | 0 | 0 | 0 |

## Known-unverified

- `build.sh` and `build_gpu.sh` were never run in this layout. Eigen and Minuit2 were absent on the consolidation machine.
- Results in `results/legacy_*` are reference outputs only. Their χ² and parameters are invalid under #1.

## Next action

1. Digonto answers D1–D5.
2. Run the `planner` (kickoff) on M1 of `prompts/000`, to write `prompts/001-*` and `tests/acceptance/m1/`, then `acceptance-lock`.
3. Start the build loop at M1 step 1 (build) and step 2 (portability). Neither needs the external data.

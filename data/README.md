# data

| folder | content | source |
|---|---|---|
| `lattice/L20/` | 19 files `20_<irrep>_n<k>.jack`: header `603 1 0 0 1`, then 603 lines `0 <a_t E_lab>`. File irreps `001/002/011` = code irreps `100/200/110`. **Sample convention unresolved** (docs/OPEN_ISSUES.md #1). | L20 hot-window package `input/lattice_spectrum` |
| `labels/` | Human true-zero labels per sign-change bracket: `bracket_id,E_left,E_right,user_label,user_zero_estimate,comment`. L20: 5 sectors, 19 true. L24: 000_A1m, 100_A2 (provisional; 100_A2 bracket 36 unresolved). The L20/100_A2 file uses the old grid start 0.26301. | L20 package `input/` |
| `windows/L20/` | `accepted_windows.csv` per sector, built from the `true` labels | L20 package `output/accepted_windows_L20_user_validated_0p36` |
| `hot_windows_meta/L20/` | Hot-window cache JSON sidecars, validation CSVs and reports. **The `.bin` payloads (about 150 MB each) are not in git.** | L20 package `output/hot_windows_L20_user_validated_0p36` |
| `brackets/` | the only surviving candidate-bracket CSV (L20/100_A2). The other sectors' det grids and bracket files must be regenerated from the raw caches. | v33p package |

Missing (needed for fits): the L24 `.jack` files, the raw caches, the det-grid CSVs and the hot-window `.bin` files. See docs/OPEN_ISSUES.md #4.

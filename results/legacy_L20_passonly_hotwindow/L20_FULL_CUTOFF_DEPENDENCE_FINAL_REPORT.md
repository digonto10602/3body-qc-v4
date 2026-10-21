# L20 full cutoff-dependence final report

## Scope and checkpoint

- Lbyas: `20` only
- Sectors: `000_A1m`, `100_A2`, `110_A2`, `111_A2`, `200_A2`
- Projected-basis cache mode: `hot_windows`
- fallback/full_scan: `0/0` for runnable fits
- L24 status: no active PIDs remained when the requested stop was issued; no L24 work was launched or resumed by this task.

## Prefit coverage

See `PREFIT_COVERAGE_AUDIT.csv/.md`. Cache files required per sector were checked, and requested `max_row_left..max_row_right` sets were compared with JSON `row_indices`.

Exact lattice input directory: `/home/digonto/Codes/KKpi_I2/spectrum/Ecm_data/data`. The fitter parser reads `20_<file_irrep>_n<state>.jack`, skips one header line, reads whitespace column 2 as `En_lab`, converts to Ecm with xi=3.444, filters by mean Ecm cutoff, and sorts by sector order and state.

## Per-cutoff results

### Ecm_cutoff=0.315

- levels/model_found: `6` / `6/6`
- chi2, dof, chi2/dof: `3.262140074829519e-05`, `2`, `1.6310700374147595e-05`
- fallback/full_scan: `0/0`
- parameters: K3iso0=`151422.46172441466` +/- `342579611.35099453`, K3iso1=`-972421.14060757787` +/- `128986425.00538744`, K3B=`347174.05548116949` +/- `134948646.91037357`, K3E=`-1226756.7068845264` +/- `8786443671.2372246`
- model levels: `Ecm_cut_0p315/model_ecm_levels.csv`
- lattice-vs-model: `Ecm_cut_0p315/lattice_vs_model_levels.csv`
- sensitivity/impact: `Ecm_cut_0p315/parameter_sensitivity.csv`, `Ecm_cut_0p315/parameter_impact_summary.csv`

### Ecm_cutoff=0.325

- levels/model_found: `9` / `9/9`
- chi2, dof, chi2/dof: `9.530896622369882e-05`, `5`, `1.9061793244739764e-05`
- fallback/full_scan: `0/0`
- parameters: K3iso0=`125968.22184657381` +/- `1.9999998807907211`, K3iso1=`-972421.14060757787` +/- `1.999999880790718`, K3B=`347174.05548116949` +/- `1.9999998807907189`, K3E=`-1226756.7068845264` +/- `1.9999998807907184`
- model levels: `Ecm_cut_0p325/model_ecm_levels.csv`
- lattice-vs-model: `Ecm_cut_0p325/lattice_vs_model_levels.csv`
- sensitivity/impact: `Ecm_cut_0p325/parameter_sensitivity.csv`, `Ecm_cut_0p325/parameter_impact_summary.csv`

### Ecm_cutoff=0.335

- levels/model_found: `12` / `12/12`
- chi2, dof, chi2/dof: `0.00012832616631556945`, `8`, `1.6040770789446181e-05`
- fallback/full_scan: `0/0`
- parameters: K3iso0=`32040.841899235049` +/- `113211204.38345505`, K3iso1=`-1074716.6512728098` +/- `84176025.162827224`, K3B=`347174.05548116943` +/- `219464496.68809414`, K3E=`-1226756.7068845264` +/- `608399252.10388672`
- model levels: `Ecm_cut_0p335/model_ecm_levels.csv`
- lattice-vs-model: `Ecm_cut_0p335/lattice_vs_model_levels.csv`
- sensitivity/impact: `Ecm_cut_0p335/parameter_sensitivity.csv`, `Ecm_cut_0p335/parameter_impact_summary.csv`

### Ecm_cutoff=0.345

- levels/model_found: `14` / `14/14`
- chi2, dof, chi2/dof: `0.00016711008405490148`, `10`, `1.6711008405490147e-05`
- fallback/full_scan: `0/0`
- parameters: K3iso0=`40607.515418119001` +/- `32989189.926181138`, K3iso1=`-1067567.8661876961` +/- `224279530.06079075`, K3B=`376187.1404411546` +/- `317856751.34917086`, K3E=`-1226756.7068845264` +/- `1116819088.8258309`
- model levels: `Ecm_cut_0p345/model_ecm_levels.csv`
- lattice-vs-model: `Ecm_cut_0p345/lattice_vs_model_levels.csv`
- sensitivity/impact: `Ecm_cut_0p345/parameter_sensitivity.csv`, `Ecm_cut_0p345/parameter_impact_summary.csv`

## Interpretation and caveats

The highest runnable cutoff is operationally preferable only if parameter drift and covariance are acceptable. These remain L20-only fits, not a global production fit. Cutoffs 0.345 and 0.355 require additional accepted classifier windows/cache coverage for the extra lattice levels before they can be fitted.

Next action: review the blocked high-cutoff sectors and decide whether to generate/validate additional accepted windows; L24 should remain inactive until explicitly resumed.

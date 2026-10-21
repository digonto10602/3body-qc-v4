# Lattice input provenance, Ecm_cutoff=0.345

- Directory: `/home/digonto/Codes/KKpi_I2/spectrum/Ecm_data/data`
- Parser: filename `20_<file_irrep>_n<state>.jack`; one header line skipped; whitespace column 2 read as En_lab.
- Ecm conversion: `sqrt(En_lab^2 - P^2)` with xi=3.444 and the irrep momentum shell.
- Levels are filtered by the mean converted Ecm <= 0.34499999999999997, sorted by the fitter order `000_A1m, 100_A2, 110_A2, 111_A2, 200_A2` and state number, then mapped to sorted accepted-window roots.
- Excluded levels are listed in `lattice_levels_excluded.csv` with reason.
- Model_found is compared with the selected lattice count in `fit_metrics.csv`.
- Uncertainties are jackknife covariance-derived Ecm sigmas; correlated chi2 uses the fitter covariance inverse.

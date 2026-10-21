Check these after running:

1. `*_terminal.log`: confirm non-interacting energies/degrees are printed for `110_A2`, `L=24`.
2. `*_projection_matrix_Vsel_audit.dat`: projection residuals should be near machine precision.
3. `*_zero_event_audit.dat`: accepted events should be `BLIND_ZERO_CANDIDATE`; dimension jumps are rejected.
4. `*_final_zero_candidates_vs_nonint.dat`: final clustered blind candidates with post-hoc nearest non-interacting energy comparison.
5. `*_projF3inv_zero_diagnostics.png`: plot of signed log determinant, min eigenvalue, min singular value, and condition number.

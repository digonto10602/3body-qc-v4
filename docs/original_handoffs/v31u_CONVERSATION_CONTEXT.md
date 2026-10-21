# Conversation Context Reconstruction

Generated: 2026-06-02 22:25:16 UTC

## Project purpose

This package preserves the KKpi three-body quantization-condition diagnostic work developed across this conversation. The project studies finite-volume three-body quantization for a two-plus-one system, specifically a `KKπ`-type system, with special attention to how `s`-wave and `p`-wave two-body partial waves enter the finite-volume matrices and the projected quantization condition.

The immediate latest deliverable is the `v31u` diagnostic:

- calculate three-body non-interacting energies over an `Ecm` region,
- calculate group-theory degeneracies of those free shells in a target lattice irrep,
- build `F2` and `G` matrices for `waves_vec_1 = {0}` and `waves_vec_1 = {0,1}`,
- compute `F2iso`, `Giso`, and stable LU-based `slogdet(F2)`, `slogdet(G)`,
- plot the functions and always overlay non-interacting energies on the `Ecm` axis.

## Physics/computational problem

The broader debugging problem was that projected three-body quantization-condition zeros and eigenvalue branches were not behaving as expected when including `s/p` wave mixing. The user expected the determinant/eigenvalue zero structure to be consistent with the physics and with earlier isotropic or `s`-wave-only diagnostics, but branch tracking, closest-to-zero eigenvalues, and determinant sign changes produced unstable or pole-contaminated candidates.

The conversation therefore built a sequence of tests to isolate the source:

1. zero finders based on projected `F3^-1` and `QC = F3^-1 + Kdf,3`,
2. closest-to-zero eigenvalue and tracked branch diagnostics,
3. basis-order and projector consistency checks,
4. finite-volume `F2/G` mixed partial-wave diagnostics,
5. non-interacting shell/degeneracy overlays,
6. stable signed-log-determinant diagnostics after raw determinants underflowed.

## Key definitions

### Momentum/irrep label convention

The working convention used throughout:

```text
000_A1m -> nP = (0,0,0)
100_A2  -> nP = (0,0,1)
110_A2  -> nP = (1,1,0)
111_A2  -> nP = (1,1,1)
200_A2  -> nP = (0,0,2)
```

The latest default is:

```text
list_of_mom = 111_A2
```

### Wave-vector settings

The latest diagnostic compares:

```text
waves_vec_1 = {0}
waves_vec_1 = {0,1}
waves_vec_2 = {0}
```

`waves_vec_1 = {0}` corresponds to the `s`-wave-only first two-body/spectator channel.  
`waves_vec_1 = {0,1}` includes `s` and `p` waves.  
`waves_vec_2 = {0}` remains `s` wave.

### Non-interacting three-body energy

For total integer momentum `nP`, two independent momenta `n1,n2`, and `n3 = nP - n1 - n2`:

```text
p_i = 2π n_i / (ξ L)

E_lab = sqrt(mK^2 + p_1^2)
      + sqrt(mK^2 + p_2^2)
      + sqrt(mπ^2 + p_3^2)

E_cm = sqrt(E_lab^2 - P^2)
```

where `P = 2π nP/(ξ L)`.

For identical kaons, the momentum key is symmetrized by sorting the two kaon momenta.

### Group-theory degeneracy

For each non-interacting energy shell, construct the representation of the little group acting on the shell states. The irrep multiplicity is computed by the character formula:

```text
m_Λ = (1/|G|) Σ_g χ_Λ(g)^* χ_shell(g)
```

where `χ_shell(g)` is the number of symmetrized shell states fixed by the little-group operation.

For one-dimensional irreps like `A2`, `group_theory_degeneracy` is the number of levels in that irrep. The projected dimension is:

```text
projected_irrep_dim = group_theory_degeneracy * dim(Λ)
```

### `F2iso` and `Giso`

The latest diagnostic builds an isotropic vector supported on `ell=0,m=0` rows:

```text
channel-1 ell=0,m=0 weight = 1
channel-2 ell=0,m=0 weight = 1/sqrt(2)
```

and computes:

```text
F2iso = iso^T F2 iso
Giso  = iso^T G  iso
```

### Stable signed log determinant

Raw determinants underflowed for large mixed `s/p` matrices. The latest version uses LU factorization:

```text
det(M) phase = permutation sign * Π_i U_ii/|U_ii|
logabsdet(M) = Σ_i log |U_ii|
slogdet(M) = sign(Re phase(det(M))) * logabsdet(M)
```

This preserves sign-change information without underflowing to zero.

## Important decisions made

- The project must be config-driven. All tweakable parameters belong in config files.
- Generated package layout must use:
  - `source/`
  - `scripts/`
  - `configs/`
  - `bin/`
  - `output/`
- Progress/status prints must appear in run scripts and executables.
- `plot_show` should be configurable.
- Non-interacting energies should always be plotted on the `Ecm` axis for the latest diagnostics.
- Raw `det(F2)`/`det(G)` are not useful for larger `s/p` bases due to numerical underflow; use `slogdet`.
- The old/new config-maker mismatch was found to be an ordering/permutation difference, not a matrix-content error.
- The basis, `F2`, `G`, `K2`, `H`, `F3`, and projectors were found consistent under permutation in `v31q`.
- Therefore the remaining likely issues concern interpretation of zeros/poles, branch tracking through dimension jumps, and whether a sign flip is a true zero or a pole/asymptote.

## Version history summary

### v29-v30 family

Early determinant and projected-F3 inverse scans. These versions explored:
- `det(projected F3^-1)`,
- eigenbranch tracking,
- zero candidates sorted by singular values,
- pole rejection,
- local bisection,
- debug plotting.

They established that simple sign flips could mix true zeros and pole/asymptote artifacts.

### v30z

Focused on minimum/closest eigenvalue local minima with pole rejection for `A2` irreps. It produced grid/eigenvalue plots and highlighted problems with closest-eigenvalue envelope switching.

### v31a-v31b

Moved to simpler sign-flip classification:
- closest-to-zero eigenvalue,
- left positive/right negative = zero,
- left negative/right positive = pole,
- no refinement in some tests,
- heavy coarse scans.

### v31e-v31l

Fitter-family packages using:
- `QC = F3inv + Kdf,3`,
- Minuit2,
- lattice covariance,
- F3inv projection caching,
- fixed-parameter spectrum scans,
- closest-to-zero eigenvalue scans,
- all-eigenvalue and tracked-eigenbranch diagnostics,
- Hungarian/global-overlap assignment.

Known compile issues included Minuit2 include paths and ROOT header requirements. The user environment had Minuit2 at:

```text
/home/digonto/.local/share/mamba/envs/minuit2-fit/include/Minuit2/Minuit2/FCNBase.h
/home/digonto/.local/share/mamba/envs/minuit2-fit/include/ROOT/RSpan.hxx
```

### v31m

Compared different inverse/determinant conventions over the whole energy region with `waves_vec_1 = {0}` and default `111_A2`.

### v31n-v31o

Compared old and new config makers:
- row-by-row flattened basis comparison,
- found that `waves_vec_1={0}` matched directly,
- found that `waves_vec_1={0,1}` did not match row-by-row in `plm_config`.

### v31p

Full mixed `s/p` consistency diagnostics:
- `F2/G` `s/p` mixing,
- `K2` off-ell block checks,
- projector orthogonality/leakage,
- inverse convention comparisons,
- QC closest-eigenvalue sign flips.

### v31q

Permutation consistency test. Key conclusion:
- old/new config makers are physically equivalent,
- complete permutation exists,
- `F2`, `G`, `K2`, `H`, `F3`, and projector transform correctly under permutation,
- the row-order mismatch is not the source of the bad QC zeros.

### v31r-v31s-v31t

Non-interacting shell and `F2/G` diagnostics:
- calculated three-body non-interacting energies,
- calculated group-theory degeneracies,
- plotted `F2iso/Giso`,
- then plotted raw determinants,
- then scaled determinants for visualization.

Problem found:
- raw determinant values for `waves_vec_1={0,1}` were zero due to numerical underflow.

### v31u latest

Current package:
- computes `F2iso/Giso`,
- computes LU-based `slogdet(F2/G)` for both `{0}` and `{0,1}`,
- always plots all non-interacting energies on the `Ecm` axis,
- output is config-driven and organized.

## Known issues and open questions

1. The latest `v31u` diagnostic is a diagnostic tool, not the final quantization-condition solver.
2. It does not resolve which QC sign flips are physical energy levels.
3. It does not perform zero refinement or pole rejection for QC.
4. Non-interacting degeneracy calculation uses a practical little-group character implementation that should be cross-checked against the exact group/projector machinery for all irreps.
5. The current package requires the external codebase headers included under `source/`.
6. Eigen is required at compile time. If `Eigen/Dense` is not found, set:
   ```bash
   export EIGEN3_INCLUDE=/path/to/eigen3
   ```
7. Minuit2 is disabled for the latest diagnostic build via:
   ```bash
   -DV31U_DISABLE_MINUIT -DV31L_DISABLE_MINUIT
   ```
8. The package includes many legacy sources/scripts for context; not every legacy executable is guaranteed to compile without additional dependencies.

## User preferences and constraints

The user explicitly requires:
- always include a config file read by executables,
- all tweakable parameters in configs,
- clear progress/status prints,
- source files under `source/`,
- scripts under `scripts/`,
- configs under `configs/`,
- generated files under `output/`,
- executables under `bin/`,
- self-contained package,
- context summary so future ChatGPT/Codex does not need the original conversation.


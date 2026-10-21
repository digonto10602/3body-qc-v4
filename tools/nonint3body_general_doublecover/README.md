# General three-body non-interacting energies with double-cover support

This package enumerates three-body non-interacting energies and decomposes each degenerate momentum shell into little-group irreps.

It is config driven. For each particle, set

```text
particleN_name
particleN_mass
particleN_J
particleN_P
particleN_C
```

The code supports integer and half-integer spin. Set

```text
double_cover = y
```

to build the double-cover little group. This prints conventional single-cover irreps when present and spinorial irreps when the shell/spin content supports them. For `P=000`, spinorial irreps are labelled `G1p/m`, `G2p/m`, and `Hp/m`. For moving frames the spinorial labels are `E1/E3` for the `Dic4`/`C4v` case and generic `G/G1/G2` labels for lower little groups.

Run:

```bash
bash run_nonint3body_general_doublecover.sh perf
```

or with a custom config:

```bash
bash run_nonint3body_general_doublecover.sh perf configs/my_config.in
```

Outputs are written to `output/`, including one combined file

```text
*_all_irreps_nonint3body.dat
```

and one file per `P + irrep`.

Energy convention:

```text
q = 2*pi/(xi*L)
E_lab = sum_i sqrt(M_i^2 + q^2 |n_i|^2)
E_cm  = sqrt(E_lab^2 - q^2 |P|^2)
n1+n2+n3 = P
```

The irrep content is obtained from the finite-group character projection. Momentum shells are acted on by the little group, and spin/parity factors are included through the configured `J^P` values. `C` is read and printed as metadata but is not used in the generic little-group decomposition.

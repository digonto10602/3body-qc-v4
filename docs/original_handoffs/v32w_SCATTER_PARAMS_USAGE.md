# Where the scatter parameters are used

The requested config values are:

```text
scatter1_00 = 0.01
scatter1_10 = 0.01
scatter2_00 = 0.01
```

They enter in three steps.

## 1. Config parsing

In the included source chain, the config reader loops over `i,j` and reads keys named
`scatter1_ij` and `scatter2_ij`:

```cpp
const std::string k1="scatter1_"+std::to_string(i)+std::to_string(j);
const std::string k2="scatter2_"+std::to_string(i)+std::to_string(j);
par.scatter_params_1[i][j]=comp(get_double(kv,k1,par.scatter_params_1[i][j].real()),0.0);
par.scatter_params_2[i][j]=comp(get_double(kv,k2,par.scatter_params_2[i][j].real()),0.0);
```

This is in:

```text
source/debug_v30r_svd_validated_cartesian_l1_100_A2_zero_scan.cpp
```

## 2. Stored in PhysicsParams

The values are stored as:

```cpp
par.scatter_params_1
par.scatter_params_2
```

## 3. Used in K2 inverse inside F3/F3inv construction

They are passed into the two-body K-matrix inverse builder:

```cpp
K2inv_EREord2_2plus1_mat(
    K2inv,
    par.eta_1,
    par.eta_2,
    par.scatter_params_1,
    par.scatter_params_2,
    En_c,
    plm,
    klm,
    total_P,
    par.atmK,
    par.atmpi,
    par.epsilon_h,
    L
);
```

In the v31zd/v31ze diagnostic executable this call appears in:

```text
source/projected_F3inv_zero_finder_v31z.cpp
```

The same call is also present in the included fitting/projection source chain:

```text
source/K3df_minuit_fit_v31l_lattice_covariance.hpp
source/debug_v30r_svd_validated_cartesian_l1_100_A2_zero_scan.cpp
```

Interpretation in this 2+1 setup:

```text
scatter_params_1 -> channel/pair type 1 used by K2inv_EREord2_2plus1_mat
scatter_params_2 -> channel/pair type 2 used by K2inv_EREord2_2plus1_mat
```

The exact physical mapping of index `(i,j)` follows the convention inside `K2inv_EREord2_2plus1_mat` in `source/K2_functions*.h` and the existing 2+1 implementation.

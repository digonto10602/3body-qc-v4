"""Basis rule (docs/METHODS.md section 2; prompts/001 section 2.1): spectator momenta n in [-N, N]^3 are kept iff
E > E_min^+(k) = omega_k + sqrt(sigma_min + (P-k)^2), sigma_min = |mj^2-mk^2| + 2 z_min (mj+mk) min(mj,mk)/(1+eps_H),
z_min = 0.02; loop order n_x, n_y, n_z from -N to N; each kept momentum carries the (l,m) list
(0,0),(1,-1),(1,0),(1,1) for channel 1 and (0,0) for channel 2; channel 1 rows first.
"""
import numpy as np

from .kinematics import channel_masses, omega

BASIS_DTYPE = np.dtype([("flavour", "i4"), ("nx", "i4"), ("ny", "i4"), ("nz", "i4"), ("ell", "i4"), ("m", "i4")])
LM_LISTS = {0: [(0, 0)], 1: [(1, -1), (1, 0), (1, 1)]}


def kept_momenta(E, P, L, flavour, M1, M2, nmax=20, z_min=0.02, eps_h=0.0):
    """Integer vectors kept for the given channel, in loop order."""
    ms, ma, mb = channel_masses(flavour, M1, M2)
    sig_min = abs(ma * ma - mb * mb) + 2.0 * z_min * (ma + mb) * min(ma, mb) / (1.0 + eps_h)
    P = np.asarray(P, dtype=float)
    out = []
    rng = range(-nmax, nmax + 1)
    for nx in rng:
        for ny in rng:
            for nz in rng:
                k = 2 * np.pi / L * np.array([nx, ny, nz], dtype=float)
                Emin = omega(k, ms) + np.sqrt(sig_min + (P - k) @ (P - k))
                if E > Emin:
                    out.append((nx, ny, nz))
    return out


def expected_basis(E, P, L, M1, M2, waves_1=(0, 1), waves_2=(0,), nmax=20, z_min=0.02, eps_h=0.0):
    """Structured array (BASIS_DTYPE) in the production order."""
    rows = []
    for flavour, waves in ((1, waves_1), (2, waves_2)):
        lms = [lm for ell in waves for lm in LM_LISTS[ell]]
        for (nx, ny, nz) in kept_momenta(E, P, L, flavour, M1, M2, nmax, z_min, eps_h):
            for (ell, m) in lms:
                rows.append((flavour, nx, ny, nz, ell, m))
    return np.array(rows, dtype=BASIS_DTYPE)


def basis_to_list(basis, L):
    """Structured basis array -> list of (flavour, pvec(3,), ell, m) with pvec = 2 pi n / L."""
    out = []
    for row in basis:
        n = np.array([row["nx"], row["ny"], row["nz"]], dtype=float)
        out.append((int(row["flavour"]), 2 * np.pi / L * n, int(row["ell"]), int(row["m"])))
    return out

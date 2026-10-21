"""Reference K_df,3 basis matrices by numerical angular projection (arXiv:2111.12734 Eqs. 28-31).

K_df,3 = Kiso0 + Kiso1*Delta + KB*(Delta_2 + Delta_2') + KE*t22~,
  Delta = (s - M^2)/M^2, Delta_2 = ((p_1 + p_1')^2 - 4 m_1^2)/M^2 (the two kaons of the initial state; primed: final),
  t22~ = (p_2 - p_2')^2/M^2 (pion initial minus pion final), M = 2 m_K + m_pi, s = E^2 - P^2.
Projection in the Q0 convention (Eqs. 20-21, 'delete every |a*|'):
  [K]_{p l'm' i; k l m j} = c_ij/(4 pi (q*_p)^{l'} (q*_k)^{l}) Int dOmega' dOmega Y_{l'm'}(a'^) K Y_{lm}(a^),
  c_11 = 1, c_12 = c_21 = 1/sqrt2, c_22 = 1/2 (Eq. 30); a* = momentum of the kaon of the pair (channel 1) in the
  pair rest frame; q* complex below threshold (the (q*)^l factors cancel analytically).
The integrand is a polynomial of degree <= 2 in each direction, so Gauss-Legendre(6) x uniform phi(12) is exact.
Returns the four basis matrices B_iso0, B_iso1, B_B, B_E with K3(c) = sum_i c_i B_i.
"""
import numpy as np
from numpy.polynomial.legendre import leggauss

from .kinematics import PI, Yreal_unit, boost_from_frame, channel_masses, omega, qstar2, sigma_pair

CIJ = {(1, 1): 1.0, (1, 2): 1 / np.sqrt(2.0), (2, 1): 1 / np.sqrt(2.0), (2, 2): 0.5}
METRIC = np.array([1.0, -1.0, -1.0, -1.0])


def sphere_grid(nth=6, nph=12):
    xg, wg = leggauss(nth)
    ph = 2 * PI * np.arange(nph) / nph
    ct, phg = np.meshgrid(xg, ph, indexing="ij")
    st = np.sqrt(1 - ct ** 2)
    hat = np.stack([st * np.cos(phg), st * np.sin(phg), ct], axis=-1).reshape(-1, 3)
    w = np.outer(wg, np.full(nph, 2 * PI / nph)).reshape(-1)
    return hat, w  # sum(w) = 4 pi


def _state(E, P, k, flav, hat, M1, M2):
    """Four-momenta (n,4) of (kaon_1, kaon_2, pion) for spectator k of flavour flav and pair directions hat."""
    ms, ma, mb = channel_masses(flav, M1, M2)
    oms = omega(k, ms)
    E2, P2 = E - oms, np.asarray(P, dtype=float) - np.asarray(k, dtype=float)
    q = np.sqrt(complex(qstar2(sigma_pair(E, P, k, ms), ma, mb)))
    n = len(hat)
    a = q * hat
    oa = np.sqrt(q * q + ma * ma); ob = np.sqrt(q * q + mb * mb)
    a0, av = boost_from_frame(np.full(n, oa, dtype=complex), a, E2, P2)
    b0, bv = boost_from_frame(np.full(n, ob, dtype=complex), -a, E2, P2)
    spec = np.concatenate([np.full((n, 1), oms, dtype=complex), np.tile(np.asarray(k, dtype=complex), (n, 1))], axis=1)
    pa = np.concatenate([a0[:, None], av], axis=1)
    pb = np.concatenate([b0[:, None], bv], axis=1)
    if flav == 1:
        return spec, pa, pb, q      # kaon (spectator), kaon (pair, +a*), pion (pair, -a*)
    return pa, pb, spec, q          # kaon, kaon, pion (spectator)


def K3_basis_matrices(basis_list, E, P, L, M1, M2):
    """basis_list: [(flavour, pvec, ell, m)], pvec = 2 pi n/L. Returns dict iso0, iso1, B, E -> (N,N) complex arrays."""
    N = len(basis_list)
    hat, w = sphere_grid()
    Mtot = 2 * M1 + M2
    P = np.asarray(P, dtype=float)
    s = E * E - P @ P
    Delta = (s - Mtot ** 2) / Mtot ** 2
    # per-row precomputation (same for rows and columns since the basis is shared)
    pre = []
    for (flav, p, ell, m) in basis_list:
        K1, K2, pion, q = _state(E, P, p, flav, hat, M1, M2)
        D2 = (np.sum((K1 + K2) ** 2 * METRIC, axis=1) - 4 * M1 * M1) / Mtot ** 2  # (n,)
        Y = Yreal_unit(hat, ell, m)                                                 # (n,)
        pre.append((pion, D2, w * Y, q ** ell, flav))
    pion_sq = [np.sum(pi ** 2 * METRIC, axis=1) for (pi, _, _, _, _) in pre]
    out = {k: np.zeros((N, N), dtype=complex) for k in ("iso0", "iso1", "B", "E")}
    for r in range(N):
        pif, D2f, wYf, qf, fi = pre[r]
        for c in range(N):
            pii, D2i, wYi, qi, fj = pre[c]
            norm = CIJ[(fi, fj)] / (4 * PI * qf * qi)
            # angular integrals: sum_f sum_i wY_f * piece(f, i) * wY_i
            sf = wYf.sum(); si = wYi.sum()
            out["iso0"][r, c] = norm * sf * si
            out["iso1"][r, c] = norm * Delta * sf * si
            out["B"][r, c] = norm * ((wYf @ D2f) * si + sf * (wYi @ D2i))
            t22 = (pion_sq[r][:, None] + pion_sq[c][None, :] - 2.0 * (pif * METRIC) @ pii.T) / Mtot ** 2
            out["E"][r, c] = norm * (wYf @ t22 @ wYi)
    return out

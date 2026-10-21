"""Reference F, G, K2^{-1} matrices for the 2+1 system (arXiv:2111.12734 Eqs. 3-15, 20-21, App. A.2).

All matrices are returned in the basis order given (list from qc3ref.basis.basis_to_list). Everything is real;
complex arithmetic is used only for x = q* L/(2 pi) below the pair threshold.
"""
import numpy as np
from scipy.special import erfi

from .kinematics import PI, H_cut, boost_to_frame, calY, channel_masses, omega, qstar2, sigma_pair


def cube_points(nmax):
    rng = np.arange(-nmax, nmax + 1)
    return np.array(np.meshgrid(rng, rng, rng, indexing="ij")).reshape(3, -1).T.astype(float)


def F_block(E, P, p, m_spec, mj, mk, L, alpha, lm_list, n_points, eps_h=0.0):
    """Diagonal-in-momentum block of F^{(i)} for spectator p, rows/cols in lm_list order (Eq. 8, App. A.2)."""
    P = np.asarray(P, dtype=float); p = np.asarray(p, dtype=float)
    om = omega(p, m_spec)
    sigma = sigma_pair(E, P, p, m_spec)
    q2 = qstar2(sigma, mj, mk)
    x = np.sqrt(complex(q2)) * L / (2 * PI)
    gam = (E - om) / np.sqrt(sigma)
    xi_jk = 0.5 * (1.0 + (mj * mj - mk * mk) / sigma)
    nP = (P - p) * L / (2 * PI)
    nP2 = nP @ nP
    n = n_points
    if nP2 == 0.0:
        r = n
    else:
        ndot = n @ nP
        r = n + np.outer((ndot / nP2) * (1.0 / gam - 1.0) - xi_jk / gam, nP)
    r2 = np.einsum("ij,ij->i", r, r)
    w = np.exp(alpha * (x * x - r2)) / (x * x - r2)
    Y = np.array([calY(r, l, m) for (l, m) in lm_list])
    S = (Y * w) @ Y.T
    for a, (la, _) in enumerate(lm_list):
        for b, (lb, _) in enumerate(lm_list):
            S[a, b] *= (2 * PI / L) ** (la + lb)
    ax2 = alpha * x * x
    e = np.exp(ax2)
    er = erfi(np.sqrt(ax2))
    I0 = 4 * PI * gam * (-0.5 * np.sqrt(PI / alpha) * e + 0.5 * PI * x * er)
    I1 = 4 * PI * gam * (-np.sqrt(PI / alpha ** 3) * (1 + 2 * ax2) / 4.0 * e + 0.5 * PI * x ** 3 * er) * (2 * PI / L) ** 2
    I = np.diag([I0 if l == 0 else I1 for (l, _) in lm_list])
    pref = H_cut(sigma, mj, mk, eps_h) / (16 * PI * PI * L ** 4 * om * (E - om))
    return pref * (S - I)


def F_matrix(basis_list, E, P, L, M1, M2, alpha=0.5, nmax=20, eps_h=0.0):
    N = len(basis_list)
    F = np.zeros((N, N), dtype=complex)
    pts = cube_points(nmax)
    a = 0
    while a < N:
        flav, p, _, _ = basis_list[a]
        idx = [a]
        while a + len(idx) < N and basis_list[a + len(idx)][0] == flav and np.array_equal(basis_list[a + len(idx)][1], p):
            idx.append(a + len(idx))
        lm_list = [(basis_list[i][2], basis_list[i][3]) for i in idx]
        ms, mj, mk = channel_masses(flav, M1, M2)
        blk = F_block(E, P, p, ms, mj, mk, L, alpha, lm_list, pts, eps_h)
        for ii, i in enumerate(idx):
            for jj, j in enumerate(idx):
                F[i, j] = blk[ii, jj]
        a += len(idx)
    return F


def G_elem(E, P, p, k, lp, mp, l, m, mi, mj, mk, L, eps_h=0.0):
    """Eq. 12: spectators p (flavour i, mass mi) and k (flavour j, mass mj); exchanged particle mass mk."""
    P = np.asarray(P, dtype=float); p = np.asarray(p, dtype=float); k = np.asarray(k, dtype=float)
    om_p = omega(p, mi); om_k = omega(k, mj)
    kstar = boost_to_frame(om_k, k, E - om_p, P - p)
    pstar = boost_to_frame(om_p, p, E - om_k, P - k)
    H1 = H_cut(sigma_pair(E, P, p, mi), mj, mk, eps_h)
    H2 = H_cut(sigma_pair(E, P, k, mj), mi, mk, eps_h)
    b = P - p - k
    denom = (E - om_p - om_k) ** 2 - b @ b - mk * mk
    return calY(kstar, lp, mp) * H1 * H2 / denom * calY(pstar, l, m) / (2 * om_p * L ** 3) / (2 * om_k * L ** 3)


def G_matrix(basis_list, E, P, L, M1, M2, eps_h=0.0):
    """Flavour matrix of Eq. 5 with the parity factor P_L of Eq. 7: G12 -> sqrt2 (-1)^{l'} G^{(12)}, G21 -> sqrt2 G^{(21)} (-1)^{l}."""
    N = len(basis_list)
    G = np.zeros((N, N), dtype=complex)
    for r, (fi, p, lp, mp) in enumerate(basis_list):
        for c, (fj, k, l, m) in enumerate(basis_list):
            if fi == 1 and fj == 1:
                G[r, c] = G_elem(E, P, p, k, lp, mp, l, m, M1, M1, M2, L, eps_h)
            elif fi == 1 and fj == 2:
                G[r, c] = np.sqrt(2.0) * (-1) ** lp * G_elem(E, P, p, k, lp, mp, l, m, M1, M2, M1, L, eps_h)
            elif fi == 2 and fj == 1:
                G[r, c] = np.sqrt(2.0) * (-1) ** l * G_elem(E, P, p, k, lp, mp, l, m, M2, M1, M1, L, eps_h)
    return G


def K2inv_elem(E, P, k, ell, eta, a, r_eff, P_shape, m_spec, mj, mk, L, pwave_form="coded_v2", eps_h=0.0):
    """Eqs. 14-15 with the 1/(2 omega_k L^3) factor; ERE Eqs. 16-17.
    ell=0: q cot d0 = -1/a0 + r0 q^2/2 + P0 q^4.
    ell=1 'q0_consistent': q^3 cot d1 = -1/a1 + r1 q^2/2 + P1 q^4, threshold term q^2 |q| (1-H)  (Q0 convention of Eqs. 20-21).
    ell=1 'coded_v2'     : -1/(a1 q^2) + r1/2 + P1 q^2, threshold term q^2 |q| (1-H)  (production form, OPEN_ISSUES #2)."""
    P = np.asarray(P, dtype=float); k = np.asarray(k, dtype=float)
    om_k = omega(k, m_spec)
    sigma = sigma_pair(E, P, k, m_spec)
    q2 = qstar2(sigma, mj, mk)
    qabs = abs(np.sqrt(complex(q2)))
    H = H_cut(sigma, mj, mk, eps_h)
    A = eta / (8 * PI * np.sqrt(sigma))
    if ell == 0:
        val = A * (-1.0 / a + 0.5 * r_eff * q2 + P_shape * q2 * q2 + qabs * (1 - H))
    elif ell == 1:
        if pwave_form == "coded_v2":
            val = A * (-1.0 / (a * q2) + 0.5 * r_eff + P_shape * q2 + q2 * qabs * (1 - H))
        elif pwave_form == "q0_consistent":
            val = A * (-1.0 / a + 0.5 * r_eff * q2 + P_shape * q2 * q2 + q2 * qabs * (1 - H))
        else:
            raise ValueError(pwave_form)
    else:
        raise ValueError("only ell <= 1")
    return val / (2 * om_k * L ** 3)


def K2inv_matrix(basis_list, E, P, L, M1, M2, a0_Kpi, a1_Kpi, a0_KK, eta_1=1.0, eta_2=0.5, pwave_form="coded_v2", eps_h=0.0):
    """diag(K2^{(1)-1}, 2 K2^{(2)-1}) (Eq. 6), diagonal in momentum and (l,m)."""
    N = len(basis_list)
    K = np.zeros((N, N), dtype=complex)
    for r, (fi, p, lp, mp) in enumerate(basis_list):
        if fi == 1:
            a = a0_Kpi if lp == 0 else a1_Kpi
            K[r, r] = K2inv_elem(E, P, p, lp, eta_1, a, 0.0, 0.0, M1, M1, M2, L, pwave_form, eps_h)
        else:
            K[r, r] = 2.0 * K2inv_elem(E, P, p, lp, eta_2, a0_KK, 0.0, 0.0, M2, M1, M1, L, pwave_form, eps_h)
    return K

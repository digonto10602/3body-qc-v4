"""Kinematics, cutoff function, harmonics and boosts (arXiv:2111.12734 Eqs. 9-11, 22, 27, App. A.1).

Conventions: units of a_t; L = xi * (L/a_s); P = 2 pi d / L; a 'frame (E2, P2)' is the rest frame of the
four-vector (E2, P2). Masses: M1 = m_K (identical pair), M2 = m_pi.
"""
import numpy as np

PI = np.pi


def J(z):
    """Smooth step J(z) (Eq. 22): 0 for z<=0, exp(-exp(-1/(1-z))/z) for 0<z<1, 1 for z>=1."""
    z = float(np.real(z))
    if z <= 0.0:
        return 0.0
    if z >= 1.0:
        return 1.0
    return float(np.exp(-np.exp(-1.0 / (1.0 - z)) / z))


def H_cut(sigma, mj, mk, eps_h=0.0):
    """Cutoff H = J(z'), z' = (1+eps_H)(sigma - |mj^2-mk^2|)/((mj+mk)^2 - |mj^2-mk^2|) (Eq. 27)."""
    dm2 = abs(mj * mj - mk * mk)
    z = (1.0 + eps_h) * (np.real(sigma) - dm2) / ((mj + mk) ** 2 - dm2)
    return J(z)


def kallen(x, y, z):
    return x * x + y * y + z * z - 2.0 * (x * y + y * z + z * x)


def omega(p, m):
    p = np.asarray(p, dtype=float)
    return np.sqrt(p @ p + m * m)


def sigma_pair(E, P, p, m_spec):
    """sigma_i = (E - omega_p)^2 - (P - p)^2 (Eq. 11)."""
    P = np.asarray(P, dtype=float); p = np.asarray(p, dtype=float)
    return (E - omega(p, m_spec)) ** 2 - (P - p) @ (P - p)


def qstar2(sigma, mj, mk):
    """(q*)^2 = lambda(sigma, mj^2, mk^2)/(4 sigma) (Eq. 11); may be negative below the pair threshold."""
    return kallen(sigma, mj * mj, mk * mk) / (4.0 * sigma)


def calY(r, ell, m):
    """Calligraphic real harmonics sqrt(4pi)|r|^l Y_lm(rhat) for l<=1 (App. A.1, Q0 convention):
    Y_00 = 1, Y_1,-1 = sqrt3 r_y, Y_10 = sqrt3 r_z, Y_11 = sqrt3 r_x. r has shape (..., 3)."""
    r = np.asarray(r)
    if ell == 0:
        return np.ones(r.shape[:-1], dtype=r.dtype)
    if ell == 1:
        idx = {1: 0, 0: 2, -1: 1}[m]
        return np.sqrt(3.0) * r[..., idx]
    raise ValueError("only ell <= 1 is implemented")


def Yreal_unit(hat, ell, m):
    """Orthonormal real spherical harmonics on unit vectors (App. A.1): Y_00 = 1/sqrt(4pi), Y_1m = sqrt(3/4pi) hat_(x,z,y)."""
    hat = np.asarray(hat)
    if ell == 0:
        return np.full(hat.shape[:-1], 1.0 / np.sqrt(4 * PI))
    if ell == 1:
        idx = {1: 0, 0: 2, -1: 1}[m]
        return np.sqrt(3.0 / (4 * PI)) * hat[..., idx]
    raise ValueError("only ell <= 1 is implemented")


def boost_to_frame(p0, p, E2, P2):
    """Spatial part of the four-vector (p0, p) in the rest frame of (E2, P2).
    beta = P2/E2; p* = p + [(gamma-1) p.bhat - gamma beta p0] bhat."""
    p = np.asarray(p, dtype=complex); P2 = np.asarray(P2, dtype=float)
    P2n = np.sqrt(P2 @ P2)
    if P2n == 0.0:
        return p
    bh = P2 / P2n
    beta = P2n / E2
    gam = 1.0 / np.sqrt(1.0 - beta * beta)
    common = (gam - 1.0) * (p @ bh) - gam * beta * p0
    return p + common * bh


def boost_from_frame(q0, q, E2, P2):
    """Inverse of boost_to_frame: four-vector (q0, q) given in the rest frame of (E2, P2) -> lab.
    Vectorised over the leading axis of q (shape (n, 3)); q0 shape (n,). Returns (q0_lab, q_lab)."""
    q = np.asarray(q); q0 = np.asarray(q0); P2 = np.asarray(P2, dtype=float)
    P2n = np.sqrt(P2 @ P2)
    if P2n == 0.0:
        return q0, q
    bh = P2 / P2n
    beta = P2n / E2
    gam = 1.0 / np.sqrt(1.0 - beta * beta)
    qpar = q @ bh
    q0_lab = gam * (q0 + beta * qpar)
    q_lab = q + np.outer((gam - 1.0) * qpar + gam * beta * q0, bh)
    return q0_lab, q_lab


def channel_masses(flavour, M1, M2):
    """(m_spectator, m_pair_a, m_pair_b) for channel 1 (K spectator; pair K pi) and 2 (pi spectator; pair K K).
    The pair particle 'a' carries +a* in the pair rest frame (for channel 1 it is the kaon)."""
    if flavour == 1:
        return M1, M1, M2
    if flavour == 2:
        return M2, M1, M1
    raise ValueError(flavour)

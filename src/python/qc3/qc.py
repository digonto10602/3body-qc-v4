"""QC helpers (prompts/001 step 4): projection of the inverse and (sign, ln|det|) by LU.

CLAUDE.md invariants 1-2: Q = L^-6 V^dag (F3inv + K3) V (project the inverse, scale the matrix, no N_proj
rescaling); the determinant is returned as (sign, ln|det|) from a partial-pivot LU, never as a raw determinant."""
import numpy as np

from .point import MissingRootError  # noqa: F401


def project_Q(F3inv, Vsel, K3, L):
    """(L^-6) V^dag (F3inv + K3) V ; L = xi*Lbyas (so L**6 == (Lbyas*xi)**6).

    The expression order is identical to the oracle's assembly (qc3.oracle) so results are bit-compatible."""
    F3inv = np.asarray(F3inv)
    Vsel = np.asarray(Vsel)
    K3 = np.asarray(K3)
    return Vsel.conj().T @ (F3inv + K3) @ Vsel / L ** 6


def slogdet(Q):
    """(sign, ln|det Q|) by partial-pivot LU (LAPACK via NumPy).  sign is an int in {-1, +1}.

    Raises ValueError if Q is not finite, is singular (det == 0 would give -inf), or if the determinant phase is
    not real (+-1 within 1e-10), i.e. Q is not (numerically) Hermitian/real-determinant."""
    Q = np.asarray(Q)
    if Q.ndim != 2 or Q.shape[0] != Q.shape[1] or Q.shape[0] == 0:
        raise ValueError(f"slogdet: need a non-empty square matrix, got shape {Q.shape}")
    if not np.all(np.isfinite(Q)):
        raise ValueError("slogdet: Q contains non-finite entries")
    sign, lndet = np.linalg.slogdet(Q)
    if sign == 0 or not np.isfinite(lndet):
        raise ValueError("slogdet: Q is singular (zero pivot in LU)")
    if abs(np.imag(sign)) > 1e-10 or abs(abs(sign) - 1.0) > 1e-10:
        raise ValueError(f"slogdet: determinant phase {sign} is not +-1; Q is not Hermitian")
    return (1 if np.real(sign) > 0 else -1), float(lndet)

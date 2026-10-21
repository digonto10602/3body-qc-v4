"""M1 step 4 — qc3.qc implements Q = (L xi)^-6 V^dag (F3inv + K3) V and (sign, ln|det|) by LU with the
project-the-inverse convention (CLAUDE.md invariants 1-2). Synthetic seeded matrices; tolerance 1e-10 on
ln|det| (float64 LU on Nproj <= 40 gives ~1e-13)."""
import numpy as np

TOL = 1e-10
L_BOX = 3.444 * 20


def _synthetic(seed=12345, N=30, Nproj=12):
    rng = np.random.default_rng(seed)
    A = rng.standard_normal((N, N)) + 1j * rng.standard_normal((N, N))
    M = (A + A.conj().T) / 2 + 3.0 * np.eye(N)          # Hermitian
    K = rng.standard_normal((N, N)); K = (K + K.T) / 2 * 0.3
    V, _ = np.linalg.qr(rng.standard_normal((N, Nproj)) + 1j * rng.standard_normal((N, Nproj)))
    return M, K, V


def test_project_Q_matches_definition():
    from qc3.qc import project_Q
    M, K, V = _synthetic()
    Q = np.asarray(project_Q(M, V, K, L_BOX), dtype=complex)
    Q_ref = V.conj().T @ (M + K) @ V / L_BOX ** 6
    assert Q.shape == (V.shape[1], V.shape[1])
    assert np.linalg.norm(Q - Q_ref) <= 1e-14 * np.linalg.norm(Q_ref)


def test_slogdet_matches_numpy_and_sign_is_integer():
    from qc3.qc import project_Q, slogdet
    M, K, V = _synthetic()
    Q = project_Q(M, V, K, L_BOX)
    s, l = slogdet(Q)
    s_ref, l_ref = np.linalg.slogdet(np.asarray(Q, dtype=complex))
    assert s in (-1, 1) and isinstance(int(s), int)
    assert int(np.sign(s_ref.real)) == s and abs(s_ref.imag) <= 1e-10
    assert abs(l - l_ref) <= TOL


def test_lndet_invariant_under_unitary_mixing_of_Vsel():
    from qc3.qc import project_Q, slogdet
    M, K, V = _synthetic(seed=7)
    rng = np.random.default_rng(99)
    U, _ = np.linalg.qr(rng.standard_normal((V.shape[1], V.shape[1])) + 1j * rng.standard_normal((V.shape[1], V.shape[1])))
    s1, l1 = slogdet(project_Q(M, V, K, L_BOX))
    s2, l2 = slogdet(project_Q(M, V @ U, K, L_BOX))
    assert s1 == s2 and abs(l1 - l2) <= TOL


def test_scale_applies_to_matrix_not_determinant():
    from qc3.qc import project_Q, slogdet
    M, K, V = _synthetic(seed=3)
    Nproj = V.shape[1]
    s1, l1 = slogdet(project_Q(M, V, K, L_BOX))
    s0, l0 = np.linalg.slogdet(V.conj().T @ (M + K) @ V)
    assert abs((l1 - l0) + 6 * Nproj * np.log(L_BOX)) <= TOL
    assert s1 == int(np.sign(s0.real))
    # flipping the sign of one eigenvalue of a full-rank Hermitian matrix flips the sign of det and leaves
    # |det| unchanged (M is Hermitian with no zero eigenvalue: seed 3 has min|w| = 0.51, 10 negative of 30).
    # Amendment 2026-10-21 (acceptance-v2): the original wrote w2[0] = -abs(w2[0]); w[0] < 0 for seed 3, so
    # M2 == M and the sign test could never pass.  Now the sign is genuinely flipped.
    w, U = np.linalg.eigh(M)
    w2 = w.copy(); w2[0] = -w2[0]
    M2 = (U * w2) @ U.conj().T
    I = np.eye(M.shape[0])
    sa, la = slogdet(project_Q(M, I, 0 * K, L_BOX))
    sb, lb = slogdet(project_Q(M2, I, 0 * K, L_BOX))
    assert sa == -sb
    assert abs(la - lb) <= TOL

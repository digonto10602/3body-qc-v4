"""Unit tests for qc3.sector and the oracle wrapper (prompts/001 step 2). Small point only (N=56 range)."""
import numpy as np
import pytest

from qc3.constants import NNP, PRODUCTION
from qc3.sector import Sector


def test_sector_properties():
    for irrep, d in NNP.items():
        s = Sector.production(irrep)
        assert s.nnP == d
        assert s.L == pytest.approx(3.444 * 20)
        np.testing.assert_allclose(s.P_lab, 2 * np.pi / s.L * np.array(d))
        assert s.k2_pwave_form == "coded_v2"
        for k in ("xi", "atmK", "atmpi", "eta_1", "eta_2", "alpha", "epsilon_H", "max_shell_num", "parity", "eig_tol"):
            assert getattr(s, k) == PRODUCTION[k]
        assert (s.a0_Kpi, s.a1_Kpi, s.a0_KK) == tuple(PRODUCTION["scatter"].values())


def test_sector_overrides_and_unknown_irrep():
    s = Sector.production("100_A2", Lbyas=24, k2_pwave_form="q0_consistent")
    assert s.Lbyas == 24 and s.L == pytest.approx(3.444 * 24) and s.k2_pwave_form == "q0_consistent"
    with pytest.raises(ValueError):
        Sector.production("123_X").nnP


def test_oracle_rejects_unknown_pwave_form():
    from qc3 import oracle
    with pytest.raises(ValueError):
        oracle.evaluate_point(Sector.production("100_A2", k2_pwave_form="nonsense"), 0.2915561)


def test_oracle_point_shapes_and_consistency():
    from qc3 import oracle
    s = Sector.production("000_A1m")
    c = (1.0e4, -2.0e4, 3.0e3, -1.0e3)
    r = oracle.evaluate_point(s, 0.26932, c)
    N, Np = r.N, r.Nproj
    assert N == r.N1 + r.N2 == len(r.basis) and Np > 0
    for name in ("F", "G", "K2inv", "F3", "F3inv", "K3"):
        assert getattr(r, name).shape == (N, N)
    assert r.Vsel.shape == (N, Np) and r.Q.shape == (Np, Np)
    assert set(r.K3_basis) == {"iso0", "iso1", "B", "E"}
    K3 = sum(ci * r.K3_basis[k] for ci, k in zip(c, ("iso0", "iso1", "B", "E")))
    np.testing.assert_allclose(r.K3, K3, rtol=0, atol=1e-12 * np.abs(K3).max())
    Q = r.Vsel.conj().T @ (r.F3inv + r.K3) @ r.Vsel / s.L ** 6
    np.testing.assert_allclose(r.Q, Q, rtol=1e-12, atol=0)
    sgn, ln = np.linalg.slogdet(r.Q)
    assert int(np.sign(sgn.real)) == r.sign and abs(ln - r.lnabsdet) < 1e-9
    assert r.metadata["backend"] == "oracle" and r.metadata["threads"] == 1
    assert {"basis", "F", "G", "K2inv", "F3", "F3inv", "projector", "K3", "projection", "lu", "total"} <= set(r.timings)
    # F3 F3inv = 1
    np.testing.assert_allclose(r.F3 @ r.F3inv, np.eye(N), atol=1e-6)


def test_oracle_want_subset():
    from qc3 import oracle
    r = oracle.evaluate_point(Sector.production("000_A1m"), 0.26932, want=("Q",))
    assert r.Q is not None and r.F is None and r.K3_basis is None and np.isfinite(r.lnabsdet)


def test_k2_pwave_forms_differ_only_in_channel1_pwave_diagonal():
    """prompts/001 section 2.2: q0_consistent replaces  -1/(a1 q^2)  by  -1/a1  in the l=1 channel-1 diagonal
    entries of K2^{-1}; everything else (basis, F, G, s-wave, channel 2) is identical."""
    from qc3 import oracle
    irrep, Ecm = "100_A2", 0.2915561
    s1 = Sector.production(irrep)
    s2 = Sector.production(irrep, k2_pwave_form="q0_consistent")
    r1 = oracle.evaluate_point(s1, Ecm)
    r2 = oracle.evaluate_point(s2, Ecm)
    assert r2.metadata["k2_pwave_form"] == "q0_consistent" and r1.metadata["k2_pwave_form"] == "coded_v2"
    assert np.array_equal(r1.basis, r2.basis)
    assert np.array_equal(r1.F, r2.F) and np.array_equal(r1.G, r2.G)
    ell, flav = np.asarray(r1.basis["ell"]), np.asarray(r1.basis["flavour"])
    pw = (ell == 1) & (flav == 1)
    assert pw.sum() > 0
    D = np.abs(r1.K2inv - r2.K2inv)
    off = ~np.eye(r1.N, dtype=bool)
    assert D[off].max() == 0.0 and (r1.K2inv[off] == 0).all()
    assert D[np.diag(~pw)].max() == 0.0 if (~pw).any() else True
    assert (np.diag(D)[pw] > 0).all()
    # difference = 1/(2 w_k L^3) * eta_1/(8 pi sqrt(sigma)) * (-1/a1) * (1 - 1/q^2)
    mK, mpi, L = s1.atmK, s1.atmpi, s1.L
    P = s1.P_lab
    for idx in np.flatnonzero(pw):
        row = r1.basis[idx]
        k = 2 * np.pi / L * np.array([row["nx"], row["ny"], row["nz"]], dtype=float)
        wk = np.sqrt(mK ** 2 + k @ k)
        sigma = (r1.E_lab - wk) ** 2 - (P - k) @ (P - k)
        q2 = (sigma - (mK + mpi) ** 2) * (sigma - (mK - mpi) ** 2) / (4 * sigma)
        expect = 1 / (2 * wk * L ** 3) * s1.eta_1 / (8 * np.pi * np.sqrt(sigma)) * (-1 / s1.a1_Kpi) * (1 - 1 / q2)
        got = (r2.K2inv[idx, idx] - r1.K2inv[idx, idx]).real
        assert abs(got - expect) <= 1e-12 * abs(expect), (idx, got, expect)

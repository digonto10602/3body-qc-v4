"""M1 step 2 — the oracle wrapper (existing C++ path, unchanged physics) reproduces its stored anchors
and satisfies the reality / Hermiticity / isometry identities. Tolerances: prompts/001 section 4."""
import numpy as np
import pytest

from m1_helpers import (C_BASELINE, C_ZERO, PRODUCTION, anchor_points, as_c, basis_list, E_lab_of, L_BOX,
                        LBYAS, NNP, oracle_point, rel_frob)

TOL_LNDET_ANCHOR = 1e-10      # same code, same flags; only the LU library differs (~1e-13 measured)
TOL_LNDET_BASELINE = 1e-9     # anchor K3 from the reference (9e-16 agreement with the port); sensitivity <= 100 eps


def _check_anchor(p, c, tol):
    r = oracle_point(p["irrep"], p["Ecm"], c)
    assert (r.N1, r.N2, r.N, r.Nproj) == (p["N1"], p["N2"], p["N"], p["Nproj"]), (p["irrep"], p["Ecm"])
    ref = p["K3_zero"] if tuple(c) == C_ZERO else p["K3_baseline"]
    assert int(r.sign) == ref["sign"], (p["irrep"], p["Ecm"], r.sign, ref["sign"])
    assert abs(r.lnabsdet - ref["lnabsdet"]) <= tol, (p["irrep"], p["Ecm"], r.lnabsdet, ref["lnabsdet"])


def test_oracle_anchors_small_points():
    pts = anchor_points(small_only=True)
    assert len(pts) >= 3
    for p in pts:
        _check_anchor(p, C_ZERO, TOL_LNDET_ANCHOR)


@pytest.mark.slow
def test_oracle_anchors_all_points():
    pts = anchor_points()
    assert len(pts) == 12
    for p in pts:
        _check_anchor(p, C_ZERO, TOL_LNDET_ANCHOR)


def test_oracle_anchors_k3_baseline_small_points():
    pts = anchor_points(small_only=True, with_baseline=True)
    assert len(pts) >= 3
    for p in pts:
        assert tuple(p["K3_baseline"]["c"]) == C_BASELINE
        _check_anchor(p, C_BASELINE, TOL_LNDET_BASELINE)


def test_oracle_reality_hermiticity_isometry():
    p = anchor_points(small_only=True)[0]
    r = oracle_point(p["irrep"], p["Ecm"], C_BASELINE)
    for name in ("F", "G", "K2inv"):
        M = as_c(getattr(r, name))
        assert np.abs(M.imag).max() <= 1e-15 * np.abs(M).max(), name
        assert rel_frob(M, M.T) <= 1e-12, f"{name} not symmetric"
    for name, M in r.K3_basis.items():
        M = as_c(M)
        assert np.abs(M.imag).max() <= 1e-15 * max(np.abs(M).max(), 1e-300), name
        assert rel_frob(M, M.T) <= 1e-12, f"K3 basis {name} not symmetric"
    Fi = as_c(r.F3inv)
    assert rel_frob(Fi, Fi.conj().T) <= 1e-14
    V = as_c(r.Vsel)
    assert V.shape == (r.N, r.Nproj)
    assert np.linalg.norm(V.conj().T @ V - np.eye(r.Nproj)) <= 1e-12
    Q = as_c(r.Q)
    assert Q.shape == (r.Nproj, r.Nproj)
    s, l = np.linalg.slogdet(Q)
    assert abs(s.imag) <= 1e-10 and int(np.sign(s.real)) == int(r.sign)
    assert abs(l - r.lnabsdet) <= 1e-10
    assert r.sign in (-1, 1)


def test_oracle_metadata_has_production_constants():
    p = anchor_points(small_only=True)[0]
    r = oracle_point(p["irrep"], p["Ecm"], C_ZERO)
    md = r.metadata
    for k, v in PRODUCTION.items():
        if k == "scatter":
            for kk, vv in v.items():
                assert md["scatter"][kk] == vv, kk
        elif isinstance(v, tuple):
            assert tuple(md[k]) == v, k
        else:
            assert md[k] == v, k
    assert md["Lbyas"] == LBYAS and md["irrep"] == p["irrep"] and tuple(md["nnP"]) == NNP[p["irrep"]]
    assert abs(md["L"] - L_BOX) <= 1e-12
    # "coded_v2" is the form PINNED by m1_helpers.sector (amendment acceptance-v4, D7): the anchors were made with it
    assert md["backend"] == "oracle" and md["k2_pwave_form"] == "coded_v2"
    assert isinstance(md["git_commit"], str) and len(md["git_commit"]) >= 7
    assert md["threads"] == 1
    assert tuple(r.c) == C_ZERO and abs(r.Ecm - p["Ecm"]) <= 1e-15
    E_lab, _ = E_lab_of(p["irrep"], p["Ecm"])
    assert abs(r.E_lab - E_lab) <= 1e-13


def test_oracle_basis_rule_and_order():
    from qc3ref.basis import expected_basis
    for p in anchor_points(small_only=True):
        r = oracle_point(p["irrep"], p["Ecm"], C_ZERO)
        E_lab, P = E_lab_of(p["irrep"], p["Ecm"])
        exp = expected_basis(E_lab, P, L_BOX, PRODUCTION["atmK"], PRODUCTION["atmpi"])
        assert len(r.basis) == len(exp) == r.N
        for name in exp.dtype.names:
            assert np.array_equal(np.asarray(r.basis[name]), exp[name]), (p["irrep"], name)
        assert r.N1 == int(np.sum(exp["flavour"] == 1)) and r.N2 == int(np.sum(exp["flavour"] == 2))


def test_oracle_timings_present():
    p = anchor_points(small_only=True)[0]
    r = oracle_point(p["irrep"], p["Ecm"], C_ZERO)
    keys = {"basis", "F", "G", "K2inv", "F3", "F3inv", "projector", "K3", "projection", "lu", "total"}
    assert keys <= set(r.timings)
    assert r.timings["total"] > 0 and r.timings["F"] > 0
    assert all(float(v) >= 0 for v in r.timings.values())

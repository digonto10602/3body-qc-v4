"""M1 step 2 — the oracle wrapper (existing C++ path, unchanged physics) reproduces its stored anchors
and satisfies the reality / Hermiticity / isometry identities. Tolerances: prompts/001 section 4."""
import numpy as np
import pytest

from m1_helpers import (C_BASELINE, C_ZERO, PRODUCTION, anchor_points, as_c, basis_list, E_lab_of, L_BOX,
                        LBYAS, NNP, oracle_point, rel_frob)

TOL_LNDET_ANCHOR = 1e-10      # FLOOR since acceptance-v8 (see anchor_tolerance); same code, same flags; only the LU library differs (~1e-13 measured)
TOL_LNDET_BASELINE = 1e-9     # FLOOR since acceptance-v8; anchor K3 from the reference (9e-16 agreement with the port); sensitivity <= 100 eps

# Amendment acceptance-v5 (audit 001 ISSUE D): NumPy slogdet(Q) vs the oracle's INDEPENDENT Eigen partial-pivot LU
# (metadata "lnabsdet_eigen_lu", "sign_eigen_lu"; C++ assembles its own Q from the same F3inv, K3, Vsel).  The two
# values differ by the summation order of the projection and by the two LU codes' rounding; a backward error dQ with
# |dQ| ~ u|Q| (u = 2^-53) moves ln|det Q| by tr(Q^-1 dQ), which for random-sign rounding over O(N_proj N) terms is
# ~ sqrt(N_proj N) kappa(Q) u <~ 30 kappa u.  Bound: 50 kappa(Q) u, floored at 1e-12 (>= 10 ulps of |ln|det|| <~ 1e3).
# Measured 2026-10-23 at the 4 small baseline anchors: kappa 8.9e4 / 2.2e6 / 3.6e6 / 3.8e12, gaps 4.3e-14 / 0 / 2.8e-14 /
# 1.1e-9, bounds 9.9e-10 / 2.4e-8 / 4.0e-8 / 4.2e-2 (margin >= 2e4 everywhere).  The sign is exact: a flip needs a
# relative perturbation >= 1/kappa >= 2.6e-13 >> u.
EIGEN_LU_KAPPA_FACTOR = 50.0
EIGEN_LU_FLOOR = 1e-12


def eigen_lu_tolerance(Q):
    """kappa(Q)-aware bound on |ln|det Q|_numpy - ln|det Q|_eigen|: max(1e-12, 50 kappa_2(Q) u)."""
    kappa = float(np.linalg.cond(np.asarray(Q, dtype=np.complex128)))
    return max(EIGEN_LU_FLOOR, EIGEN_LU_KAPPA_FACTOR * kappa * np.finfo(float).eps), kappa


# Amendment acceptance-v8 (2026-10-27, prompts/001 section 4): the anchor comparisons below used FIXED bounds
# (1e-10 for c = 0, 1e-9 for the baseline c), which are under the BLAS round-off floor at ill-conditioned points.  The
# anchors were produced on another CPU; the oracle's ln|det Q| is NumPy slogdet(Q) = sum ln|diag LU|, and a different
# OpenBLAS kernel changes the summation order of F, G, F3^{-1}, the projection V^dagger (.) V and the LU.  First-order
# analysis: Q = s V^dagger (F3^{-1} + K3) V; a perturbation dQ with ||dQ|| ~ u ||Q|| (u = 2^-53 ~ 1.1e-16, eps = 2.2e-16)
# moves ln|det Q| by Re tr(Q^{-1} dQ), |.| <= N_proj ||Q^{-1}|| ||dQ||, i.e. <~ sqrt(N_proj N) kappa_2(Q) u for
# random-sign rounding (same estimate as test D above).  The operative condition number is kappa_2(Q) of the matrix
# whose LU is taken, NOT kappa(F3): Q adds K3 and projects, and at 100_A2 0.2915561 with the baseline c, kappa(F3) = 9.0e6
# but kappa(Q) = 3.8e12 and the observed deviation is 1.1e-9 (Sandybridge kernel), which kappa(F3) would not cover and
# kappa(Q) does.  At c = 0 the two coincide or differ by <= 5 (measured).  Bound: max(fixed, 50 kappa_2(Q) eps), K = 50
# as in test D; measured deviation/(kappa eps) <= 0.2 at all 16 (point, c) pairs and both kernels tested (Haswell,
# Sandybridge), so the margin is >= 250 at every point.  The sign must match exactly (a flip needs a relative
# perturbation >= 1/kappa(Q) >> u).  The fixed bounds are kept as floors, so nothing is tighter than before anywhere.
ANCHOR_KAPPA_FACTOR = 50.0


def anchor_tolerance(Q, fixed):
    """kappa-aware bound on |ln|det Q|_oracle - ln|det Q|_anchor|: max(fixed, 50 kappa_2(Q) eps); returns (tol, kappa)."""
    kappa = float(np.linalg.cond(np.asarray(Q, dtype=np.complex128)))
    return max(fixed, ANCHOR_KAPPA_FACTOR * kappa * np.finfo(float).eps), kappa


def _check_anchor(p, c, fixed_tol):
    r = oracle_point(p["irrep"], p["Ecm"], c)
    assert (r.N1, r.N2, r.N, r.Nproj) == (p["N1"], p["N2"], p["N"], p["Nproj"]), (p["irrep"], p["Ecm"])
    ref = p["K3_zero"] if tuple(c) == C_ZERO else p["K3_baseline"]
    tol, kappa = anchor_tolerance(as_c(r.Q), fixed_tol)
    tag = (p["irrep"], p["Ecm"], "kappa", kappa, "tol", tol)
    assert int(r.sign) == ref["sign"], (tag, r.sign, ref["sign"])
    assert abs(r.lnabsdet - ref["lnabsdet"]) <= tol, (tag, r.lnabsdet, ref["lnabsdet"])


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
    assert r.sign in (-1, 1)
    # Amendment acceptance-v5 (audit 001 ISSUE D).  Before: `abs(l - r.lnabsdet) <= 1e-10`, but r.lnabsdet IS
    # np.linalg.slogdet(Q) (qc3.oracle), so that compared NumPy with itself.  Now: NumPy's LU vs the oracle's independent
    # Eigen LU, at every small baseline anchor (the first one, kappa = 3.8e12, gets a loose bound; the other three,
    # kappa <= 3.6e6, get bounds <= 4e-8), with the kappa-aware tolerance defined above.
    for q in anchor_points(small_only=True, with_baseline=True):
        rq = oracle_point(q["irrep"], q["Ecm"], C_BASELINE)
        tol, kappa = eigen_lu_tolerance(rq.Q)
        tag = (q["irrep"], q["Ecm"], "kappa", kappa, "tol", tol)
        assert int(rq.metadata["sign_eigen_lu"]) == int(rq.sign), (tag, rq.metadata["sign_eigen_lu"], rq.sign)
        gap = abs(float(rq.metadata["lnabsdet_eigen_lu"]) - rq.lnabsdet)
        assert gap <= tol, (tag, "gap", gap, rq.lnabsdet, rq.metadata["lnabsdet_eigen_lu"])


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

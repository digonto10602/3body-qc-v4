"""Unit tests for qc3.fit.QCChi2 / make_trackers (M2 step 5, prompts/003 section 2.6).  Synthetic trackers whose roots are
affine in c (E_k(c) = a_k + J_k . c, exact Jacobian), so the gradient, D3 with n > m, the penalty path, D9 and the call counting
are tested without the C++ core; one real-core test checks D8 (a free-level solution is never matched).  No randomness."""
from types import SimpleNamespace

import numpy as np
import pytest

from qc3 import fit, lattice
from qc3.point import LevelSet
from qc3.roots import AmbiguousRootError
from qc3.tracking import RootSetChangedError

C0 = np.array([1.0, -2.0, 0.5, 0.25])


class FakeTracker:
    """Affine roots E_k(c) = a_k + J_k . (c - C0); `slopes` < 0 marks wrong-sign zeros; `drop_below` removes roots at the top
    (to mimic a sector that loses levels)."""

    def __init__(self, irrep, a, J, slopes=None, n_keep=None, Lbyas=20, E_max=0.36):
        self.sector = SimpleNamespace(irrep=irrep, Lbyas=Lbyas)
        self.a, self.J = np.asarray(a, float), np.asarray(J, float)
        self.slopes = np.ones(len(a)) if slopes is None else np.asarray(slopes, float)
        self.n_keep = len(a) if n_keep is None else n_keep
        self.c, self.E_max, self.n_updates, self.n_evaluations, self.n_rescans = C0.copy(), E_max, 0, 0, 0

    def update(self, c):
        self.c = np.asarray(c, float).copy()
        self.n_updates += 1
        return self.energies()

    def energies(self):
        return (self.a + self.J @ (self.c - C0))[: self.n_keep]

    def jacobian(self):
        return self.J[: self.n_keep].copy()

    def wrong_sign_roots(self):
        E = self.energies()
        return [SimpleNamespace(E=float(E[k]), slope=float(self.slopes[k])) for k in range(self.n_keep) if self.slopes[k] < 0]

    def rescan(self):
        self.n_rescans += 1
        return None


def make_levels(labels, E, sigma_rel=1e-3):
    m = len(labels)
    E = np.asarray(E, float)
    s = sigma_rel * np.ones(m)
    cov = np.outer(s, s) * (0.3 + 0.7 * np.eye(m))        # correlated, SPD
    irreps = sorted({l[0] for l in labels})
    return LevelSet(labels=labels, irreps=irreps, nP=[(0, 0, 0)] * m, files=[""] * m, Ncfg=10, Ecm=E, cov=cov, Ecm_cutoff=0.36, Lbyas=20)


def setup(n_keep=None, slopes=None):
    # sector A: three QC roots, two lattice levels (D3 with n > m); sector B: two roots, two levels
    JA = np.array([[-1e-6, 2e-7, 0, 1e-7], [-2e-6, 0, 3e-7, 0], [-3e-6, 1e-7, 1e-7, 5e-7]])
    JB = np.array([[-1.5e-6, 1e-7, 2e-7, 0], [-2.5e-6, 3e-7, 0, 2e-7]])
    tA = FakeTracker("000_A1m", [0.280, 0.300, 0.330], JA, slopes=slopes, n_keep=n_keep)
    tB = FakeTracker("100_A2", [0.290, 0.310], JB)
    labels = [("000_A1m", 0), ("000_A1m", 1), ("100_A2", 0), ("100_A2", 1)]
    lv = make_levels(labels, [0.2803, 0.3004, 0.2899, 0.3103])
    return lv, {"000_A1m": tA, "100_A2": tB}


def test_chi2_equals_lattice_chi2_of_the_matched_roots_and_d3_uses_lowest():
    lv, tr = setup()
    f = fit.QCChi2(lv, tr)
    c = C0 + np.array([10.0, -5.0, 3.0, 1.0]) * 1e3
    E, J = f.model(c)
    # sector A: n=3 > m=2 -> the lowest two roots (0.28.., 0.30..), not the 0.33 one
    assert E[0] == tr["000_A1m"].energies()[0] and E[1] == tr["000_A1m"].energies()[1]
    assert E[2] == tr["100_A2"].energies()[0] and E[3] == tr["100_A2"].energies()[1]
    assert np.array_equal(J[0], tr["000_A1m"].J[0]) and np.array_equal(J[3], tr["100_A2"].J[1])
    assert f(c) == lattice.chi2(lv.Ecm - E, lv.cov)
    assert f.n_penalised_calls == 0 and f.errordef == 1.0


def test_gradient_matches_finite_differences_of_chi2():
    lv, tr = setup()
    f = fit.QCChi2(lv, tr)
    c = C0 + np.array([3.0, -2.0, 1.0, 4.0]) * 1e3
    g = f.grad(c)
    g_fd = np.empty(4)
    for i in range(4):
        d = np.zeros(4); d[i] = 1e-3 * max(abs(c[i]), 1.0)
        g_fd[i] = (f(c + d) - f(c - d)) / (2 * d[i])
    assert np.linalg.norm(g - g_fd) <= 1e-8 * np.linalg.norm(g_fd)      # affine roots: chi^2 is quadratic, central FD exact
    assert f.n_grad_calls == 1 and f.n_calls == 8


def test_sign_convention_residual_is_lattice_minus_model():
    lv, tr = setup()
    f = fit.QCChi2(lv, tr)
    E, _ = f.model(C0)
    assert np.allclose(f.residual(C0), lv.Ecm - E, rtol=0, atol=0)


def test_penalty_path_is_finite_counted_and_has_zero_gradient_rows():
    lv, tr = setup(n_keep=1)                                           # sector A keeps 1 root for 2 levels
    f = fit.QCChi2(lv, tr, penalty_sigmas=100.0)
    E, J = f.model(C0)
    assert np.all(np.isfinite(E)) and np.all(np.isfinite(J))
    assert f.penalised.tolist() == [False, True, False, False]         # the highest lattice level of the sector is penalised
    sigma = np.sqrt(np.diag(lv.cov))
    assert f.residual(C0)[1] == pytest.approx(100.0 * sigma[1], rel=0, abs=1e-15)
    assert np.all(J[1] == 0.0)
    chi2 = f(C0)
    assert np.isfinite(chi2) and 1e3 < chi2 < 1e100
    assert f.n_penalised_calls == 1                                    # f(C0) after model(C0): cached, counted once
    f(C0 + 1e3)                                                        # new c -> new model, counted again
    assert f.n_penalised_calls == 2
    g = f.grad(C0)
    assert np.all(np.isfinite(g))
    # the penalised level contributes no gradient of its own
    r = f.residual(C0)
    assert np.allclose(g, -2.0 * J.T @ (f.Cinv @ r), rtol=1e-12, atol=0)


def test_penalty_is_a_failure_of_the_whole_sector_not_an_exception():
    lv, tr = setup(n_keep=0)                                           # sector A has no root at all
    f = fit.QCChi2(lv, tr)
    assert f(C0) > 1e4 and f.penalised.tolist() == [True, True, False, False]


def test_d9_wrong_sign_zeros_are_flagged_and_stay_in_the_matching():
    lv, tr = setup(slopes=[1.0, -1.0, -1.0])                           # 2nd (matched) and 3rd (unmatched, n > m) are wrong-sign
    f = fit.QCChi2(lv, tr)
    c = C0 + 1e3
    E, _ = f.model(c)
    assert E[1] == tr["000_A1m"].energies()[1]                          # the wrong-sign zero is still matched (D3 unchanged)
    w = f.wrong_sign_zeros
    assert [x["matched"] for x in w] == [True, False] and [x["level_index"] for x in w] == [1, None]
    assert all(x["irrep"] == "000_A1m" and x["slope"] < 0 for x in w)
    assert w[0]["E"] == E[1] and f.n_calls_with_wrong_sign == 1
    f(c); f.grad(c)                                                     # same c: cached, not counted again
    assert f.n_calls_with_wrong_sign == 1
    f(c + 1.0)
    assert f.n_calls_with_wrong_sign == 2
    lv2, tr2 = setup()
    f2 = fit.QCChi2(lv2, tr2); f2(C0)
    assert f2.wrong_sign_zeros == [] and f2.n_calls_with_wrong_sign == 0


def test_tracker_errors_propagate_and_verify_roots_calls_rescan():
    lv, tr = setup()
    f = fit.QCChi2(lv, tr)
    f.verify_roots()
    assert tr["000_A1m"].n_rescans == 1 and tr["100_A2"].n_rescans == 1

    class Boom(FakeTracker):
        def update(self, c):
            raise AmbiguousRootError("loud")

    tr["000_A1m"].__class__ = Boom
    with pytest.raises(AmbiguousRootError):
        fit.QCChi2(lv, tr)(C0 + 1.0)

    class Changed(FakeTracker):
        def rescan(self):
            raise RootSetChangedError("a zero entered from above E_max")

    tr["100_A2"].__class__ = Changed
    with pytest.raises(RootSetChangedError):
        fit.QCChi2(lv, tr).verify_roots()


def test_constructor_validation():
    lv, tr = setup()
    with pytest.raises(ValueError):
        fit.QCChi2(lv, {"000_A1m": tr["000_A1m"]})                      # missing tracker
    with pytest.raises(ValueError):
        fit.QCChi2(lv, tr, penalty_sigmas=0.0)
    with pytest.raises(ValueError):
        fit.QCChi2(lv, tr)(np.zeros(3))
    bad = make_levels([("000_A1m", 0)], [0.28])
    bad.cov = np.array([[-1.0]])
    with pytest.raises(ValueError):
        fit.QCChi2(bad, tr)


def test_cov_eig_pinv_agrees_with_lattice_chi2():
    lv, _ = setup()
    r = np.array([1e-4, -2e-4, 3e-4, 5e-5])
    w, U = lattice.cov_eig_pinv(lv.cov)
    assert float(np.sum((U.T @ r) ** 2 / w)) == lattice.chi2(r, lv.cov)
    assert np.allclose((U / w) @ U.T, np.linalg.inv(lv.cov), rtol=1e-9)


def test_real_trackers_d8_free_level_solution_is_never_matched():
    """110_A2, coded_v2, baseline, window [0.3190, 0.3200]: one genuine zero (0.31915494) and a free-level solution (0.319468828,
    D8) in the same window.  A lattice level next to the free level must be matched to the genuine zero."""
    from qc3.constants import C_BASELINE
    tr = fit.make_trackers(["110_A2"], C_BASELINE, Lbyas=20, k2_pwave_form="coded_v2", E_min=0.3190, E_max=0.3200)
    assert len(tr["110_A2"].free_level_solutions) == 1
    lv = make_levels([("110_A2", 0)], [0.31947])
    f = fit.QCChi2(lv, tr)
    E, J = f.model(np.asarray(C_BASELINE))
    assert abs(E[0] - 0.3191549427) < 1e-8 and J.shape == (1, 4) and J[0, 0] < 0
    assert f.n_penalised_calls == 0
    # a window holding only the free level: the sector has no zero -> penalty, not a match to the free level
    tr2 = fit.make_trackers(["110_A2"], C_BASELINE, Lbyas=20, k2_pwave_form="coded_v2", E_min=0.31930, E_max=0.31990)
    assert len(tr2["110_A2"].energies()) == 0 and len(tr2["110_A2"].free_level_solutions) == 1
    f2 = fit.QCChi2(lv, tr2)
    assert f2(np.asarray(C_BASELINE)) > 1e4 and f2.n_penalised_calls == 1

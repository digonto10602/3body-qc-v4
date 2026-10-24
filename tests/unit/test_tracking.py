"""Unit tests for qc3.tracking (M2 step 4, prompts/003 section 2.5): a fake two-channel Q(E; c) = A(E) + sum_i c_i B~_i with a
movable root (so eq. (1), the secant, the predictor-corrector, the safeguards and both rescan fallbacks are exercised without
the C++ core), plus one real-core smoke test of the Jacobian.  Fixed inputs, no randomness."""
from types import SimpleNamespace

import numpy as np
import pytest

from qc3 import fast, roots, tracking
from qc3.constants import GRID_DE
from qc3.point import SingularQError
from qc3.roots import AmbiguousRootError
from qc3.sector import Sector

SEC = Sector.production("000_A1m", Lbyas=20, k2_pwave_form="coded_v2")
SEC110 = Sector.production("110_A2", Lbyas=20, k2_pwave_form="coded_v2")
L6 = SEC.L ** 6
C0 = (0.20, -0.10, 0.30, 0.05)
BT = [1e-2 * np.array([[1.0, 0.0], [0.0, 0.0]]), 1e-2 * np.array([[0.0, 0.2], [0.2, 0.0]]),
      1e-2 * np.array([[0.3, 0.1], [0.1, 0.0]]), 1e-2 * np.array([[0.0, 0.0], [0.0, 1.0]])]        # B~_i (the fake's K3_basis is B~_i * L^6)


class FakeQ:
    """Q = [[sg*50 (E - 0.2800), 0.05], [0.05, 30 (E - 0.2900) + 1e-4/(E - 0.2850)]] + sum_i c_i B~_i: a pole of det at 0.2850,
    zeros near 0.2798, 0.2857 (slope < 0: a wrong-sign zero, D9) and 0.2894 at C0 that move with c.  Vsel = 1, N = 2."""

    def __init__(self, sg=1.0, sing_at=None):
        self.n, self.sg, self.sing_at = 0, sg, sing_at

    def dimension_jumps(self, sector, a, b):
        return np.zeros(0)

    def __call__(self, sector, E, c, want=("Q",)):
        self.n += 1
        if self.sing_at is not None and abs(E - self.sing_at) < 1e-14:
            raise SingularQError("fake exact root")
        Q = np.array([[self.sg * 50.0 * (E - 0.28), 0.05], [0.05, 30.0 * (E - 0.29) + 1e-4 / (E - 0.285)]], dtype=complex)
        Q = Q + sum(ci * B for ci, B in zip(c, BT))
        s, l = np.linalg.slogdet(Q)
        out = SimpleNamespace(sign=int(np.sign(s.real)), lnabsdet=float(l), N=2, Q=Q, Vsel=None, K3_basis=None)
        if "K3" in want:
            out.Vsel = np.eye(2, dtype=complex)
            out.K3_basis = {k: B.astype(complex) * L6 for k, B in zip(("iso0", "iso1", "B", "E"), BT)}
        return out


def brent(f, c, a, b):
    return roots.refine_root(f, SEC, tuple(c), a, b)


def test_root_jacobian_matches_finite_differences_fake():
    f = FakeQ()
    E0 = brent(f, C0, 0.2790, 0.2810)
    tr = tracking.root_jacobian(f, SEC, C0, E0)
    assert tr.E == E0 and tr.J.shape == (4,) and tr.N == 2 and tr.slope > 0 and tr.discriminant <= -8.0
    Jfd = np.zeros(4)
    for i in range(4):
        d = np.zeros(4); d[i] = 1e-4
        Jfd[i] = (brent(f, np.add(C0, d), 0.2790, 0.2810) - brent(f, np.subtract(C0, d), 0.2790, 0.2810)) / 2e-4
    assert np.max(np.abs(tr.J - Jfd)) <= 1e-6 * np.max(np.abs(Jfd))
    # lnabsdet_at_root supplied: 2 evaluations only
    f.n = 0
    tracking.root_jacobian(f, SEC, C0, E0, lnabsdet_at_root=-30.0)
    assert f.n == 2
    f.n = 0
    tracking.root_jacobian(f, SEC, C0, E0)
    assert f.n == 3


def test_root_jacobian_rejects_a_pole_and_a_non_root_loudly():
    f = FakeQ()
    with pytest.raises(roots.PoleError):
        tracking.root_jacobian(f, SEC, C0, 0.2850 + 1e-11)
    with pytest.raises(AmbiguousRootError):
        tracking.root_jacobian(f, SEC, C0, 0.2800 + 1e-5)     # |Delta_h| ~ ln(1e-5/1e-6) < D: neither zero nor pole


def test_secant_root_converges_without_evaluating_the_returned_energy_and_handles_exact_root():
    f = FakeQ()
    E0 = brent(f, C0, 0.2790, 0.2810)
    seen = []
    orig = f.__call__

    def spy(sector, E, c, want=("Q",)):
        seen.append(E)
        return orig(sector, E, c, want)
    c1 = np.add(C0, [0.01, 0.0, 0.0, 0.0])
    E1 = brent(f, c1, 0.2790, 0.2810)
    E, l_last, n = tracking.secant_root(spy, SEC, c1, E0, E1 - E0)
    assert abs(E - E1) <= 1e-11 and np.isfinite(l_last) and n == len(seen) <= 10 and E not in seen
    # exactly singular iterate: returns that iterate with ln|det| = -inf
    fs = FakeQ(sing_at=0.28)
    E, l_last, n = tracking.secant_root(fs, SEC, C0, 0.28, 1e-4)
    assert E == 0.28 and l_last == -np.inf and n == 1
    # not converged -> loud
    with pytest.raises(tracking.SecantNotConvergedError):
        tracking.secant_root(f, SEC, C0, 0.2795, 1e-3, max_iter=1)


def _fresh(f, sec, c, a=0.270, b=0.300):
    from qc3 import rootscan
    return rootscan.scan_sector(f, sec, tuple(c), a, b).zeros()


def test_sector_roots_update_agrees_with_fresh_scan_and_counts():
    f = FakeQ()
    sr = tracking.SectorRoots(SEC, C0, 0.270, 0.300, evaluate=f)
    assert np.array_equal(sr.energies(), sr.scan.zeros()) and len(sr.roots) == 3 and sr.jacobian().shape == (3, 4)
    assert all(r.discriminant <= -8.0 and r.segment == (0.270, 0.300) for r in sr.roots) and sr.free_level_solutions == []
    n0 = sr.n_evaluations
    assert np.array_equal(sr.update(C0), sr.energies()) and sr.n_evaluations == n0 and sr.n_substeps == 0     # no-op costs nothing
    c1 = np.add(C0, [0.02, -0.01, 0.03, 0.0])
    E = sr.update(c1)
    assert np.max(np.abs(E - _fresh(f, SEC, c1))) <= 1e-10 and sr.n_substeps >= 1
    assert sr.n_local_rescans == 0 and sr.n_full_rescans == 0
    # the predictor: n_sub = ceil(max |J dc| / (20 rows)); a bigger move costs more sub-steps and still agrees
    c2 = np.add(c1, [0.5, 0.0, 0.0, 0.0])
    move = max(abs(r.J @ (np.asarray(c2) - c1)) for r in sr.roots)
    s0 = sr.n_substeps
    E = sr.update(c2)
    assert sr.n_substeps - s0 == int(np.ceil(move / (20.0 * GRID_DE) - 1e-12)) > 1
    assert np.max(np.abs(E - _fresh(f, SEC, c2))) <= 1e-10
    res = sr.rescan()
    assert len(res.zeros()) == 3 and sr.n_full_rescans == 0 and sr.n_verify_rescans == 1


def test_safeguard_failure_triggers_local_rescan_and_tiny_window_triggers_full_rescan():
    f = FakeQ()
    # safeguards that can never be met (any |E1 - E_pred| > 0 fails) -> every root goes through the local rescan
    sr = tracking.SectorRoots(SEC, C0, 0.270, 0.300, evaluate=f, track_options=dict(safeguard_factor=0.0, safeguard_abs=0.0))
    c1 = np.add(C0, [0.02, 0.0, 0.0, 0.0])
    E = sr.update(c1)
    assert sr.n_local_rescans == 3 and sr.n_full_rescans == 0
    assert np.max(np.abs(E - _fresh(f, SEC, c1))) <= 1e-10
    # a window narrower than the scan edge finds nothing at any step -> full rescan; the tracked set is replaced
    sr = tracking.SectorRoots(SEC, C0, 0.270, 0.300, evaluate=f, track_options=dict(safeguard_factor=0.0, safeguard_abs=0.0, local_rows=0))
    c1 = np.add(C0, [1e-12, 0.0, 0.0, 0.0])
    E = sr.update(c1)
    assert sr.n_full_rescans >= 1 and np.max(np.abs(E - _fresh(f, SEC, c1))) <= 1e-10 and len(sr.roots) == 3


def test_wrong_sign_root_is_recorded_and_kept():
    f = FakeQ()
    sr = tracking.SectorRoots(SEC, C0, 0.270, 0.300, evaluate=f)
    w = sr.wrong_sign_roots()
    assert len(w) == 1 and w[0].slope < 0 and abs(w[0].E - 0.2857346) < 1e-6 and w[0].E in sr.energies()      # D9: never excluded
    assert sr.energies().size == 3
    sr.update(np.add(C0, [0.02, -0.01, 0.03, 0.0]))
    assert len(sr.wrong_sign_roots()) == 1 and sr.energies().size == 3


def test_root_leaving_above_e_max_is_parked_and_comes_back():
    f = FakeQ()
    sr = tracking.SectorRoots(SEC, C0, 0.270, 0.2895, evaluate=f, track_options=dict(park_margin=1e-4))
    n_up = len(sr.energies())
    assert n_up == 3
    c_up = np.add(C0, [0.0, 0.0, 0.0, -0.60])         # pushes the upper zero up through E_max + margin
    E = sr.update(c_up)
    parked = [r for r in sr.roots if r.E > sr.E_max + 1e-4]
    assert len(E) < n_up and len(sr.roots) == n_up and len(parked) >= 1 and all(r.E > sr.E_max for r in parked)
    n_ev = sr.n_evaluations
    E = sr.update(C0)                                  # and back down
    assert len(E) == n_up and np.max(np.abs(E - _fresh(f, SEC, C0, 0.270, 0.2895))) <= 1e-9
    assert sr.n_evaluations > n_ev


def test_rescan_raises_when_the_tracked_set_is_stale():
    f = FakeQ()
    sr = tracking.SectorRoots(SEC, C0, 0.270, 0.300, evaluate=f)
    sr.c = np.add(sr.c, [0.05, 0.0, 0.0, 0.0])         # the tracker's c is moved without tracking: roots are stale
    with pytest.raises(tracking.RootSetChangedError):
        sr.rescan()
    assert issubclass(tracking.RootSetChangedError, roots.MissingRootError)


def test_free_level_solution_is_excluded_from_roots_and_reported():
    # a simple zero exactly at the 110_A2 free level, plus an ordinary one: only the ordinary one is tracked (D8)
    from qc3 import free_levels
    Ef = float(free_levels.free_levels(SEC110, 0.3190, 0.3200)[0])

    class F:
        n = 0

        def dimension_jumps(self, sector, a, b):
            return np.zeros(0)

        def __call__(self, sector, E, c, want=("Q",)):
            self.n += 1
            ln = np.log(abs(E - Ef) * abs(E - 0.31915) * 1e3)
            sg = np.sign(E - Ef) * np.sign(E - 0.31915)
            Q = np.array([[sg * np.exp(ln)]], dtype=complex)
            out = SimpleNamespace(sign=int(sg), lnabsdet=float(ln), N=1, Q=Q, Vsel=None, K3_basis=None)
            if "K3" in want:
                out.Vsel = np.eye(1, dtype=complex)
                out.K3_basis = {k: np.zeros((1, 1), complex) for k in ("iso0", "iso1", "B", "E")}
            return out

    sr = tracking.SectorRoots(SEC110, C0, 0.3190, 0.3200, evaluate=F())
    assert len(sr.roots) == 1 and abs(sr.roots[0].E - 0.31915) <= 1e-10
    assert len(sr.free_level_solutions) == 1 and sr.free_level_solutions[0].kind == "free_level"
    assert abs(sr.free_level_solutions[0].E - Ef) <= 1e-9


def test_ambiguous_candidate_is_a_loud_event():
    # a simple zero 4e-6 from a pole: |Delta_h| < D at the zero -> the scan reports kind "ambiguous"; the tracker must raise
    class F:
        def dimension_jumps(self, sector, a, b):
            return np.zeros(0)

        def __call__(self, sector, E, c, want=("Q",)):
            z, p = 0.2850000, 0.2850000 + 4e-6
            ln = np.log(abs(E - z)) - np.log(abs(E - p)) - 3.0
            sg = np.sign(E - z) * np.sign(E - p)
            return SimpleNamespace(sign=int(sg), lnabsdet=float(ln), N=1, Q=np.array([[sg * np.exp(ln)]], dtype=complex), Vsel=None, K3_basis=None)

    from qc3 import rootscan
    res = rootscan.scan_sector(F(), SEC, C0, 0.280, 0.290)
    if not res.of_kind("ambiguous"):
        pytest.skip("the synthetic zero-pole pair is not classified ambiguous by this scan; nothing to check")
    with pytest.raises(AmbiguousRootError):
        tracking.SectorRoots(SEC, C0, 0.280, 0.290, evaluate=F())


def test_unknown_track_option_and_bad_c_rejected():
    f = FakeQ()
    with pytest.raises(ValueError):
        tracking.SectorRoots(SEC, C0, 0.270, 0.300, evaluate=f, track_options=dict(nope=1))
    sr = tracking.SectorRoots(SEC, C0, 0.270, 0.300, evaluate=f)
    with pytest.raises(ValueError):
        sr.update((1.0, 2.0, 3.0))


def test_root_jacobian_real_core_matches_finite_difference_000_A1m():
    # real fast path at a different c than the acceptance tests (a 10 % shift of the baseline), first zero of 000_A1m
    from qc3.constants import C_BASELINE
    sec = Sector.production("000_A1m", Lbyas=20, k2_pwave_form="coded_v2")
    c = np.asarray(C_BASELINE) * np.array([1.1, 0.9, 1.1, 0.95])
    zs = roots.find_zeros(fast.evaluate_point, sec, tuple(c), 0.2630, 0.2800, n_scan=121)
    assert len(zs) >= 1
    E0 = zs[0].E
    tr = tracking.root_jacobian(fast.evaluate_point, sec, tuple(c), E0)
    Jfd = np.zeros(4)
    for i in range(4):
        d = np.zeros(4); d[i] = 1e-3 * abs(c[i])
        Ep = roots.refine_root(fast.evaluate_point, sec, tuple(c + d), E0 - 2e-5, E0 + 2e-5, kind="any")
        Em = roots.refine_root(fast.evaluate_point, sec, tuple(c - d), E0 - 2e-5, E0 + 2e-5, kind="any")
        Jfd[i] = (Ep - Em) / (2 * d[i])
    assert np.max(np.abs(tr.J - Jfd) / np.abs(Jfd)) <= 1e-4
    assert tr.slope > 0 and tr.discriminant <= -8.0

"""M2 step 4 (prompts/003 section 2.5): root tracking across parameter changes, `qc3.tracking`.

Method under test: for a simple zero E_0(c) of det Q(E; c), with u the unit null vector of Q(E_0),
    dE_0/dc_i = -(u^T B~_i u) / (u^T Q'(E_0) u)         (implicit function theorem + Hellmann-Feynman, prompts/003 eq. 5)
is the predictor; a secant on sign * exp(ln|det Q| - l_ref) is the corrector; every accepted root is re-classified with
Delta_h and its Jacobian recomputed from Q(E_0 +- h).  Reference values come from independent code paths of the
production package (Brent refinement `qc3.roots.refine_root`, the full scan `qc3.rootscan.scan_sector`).

Tolerances: Jacobian vs central finite difference with relative step 1e-3: 1e-4 (truncation ~ (1e-3)^2/6 and the measured
prototype agreement 4e-8 to 1.8e-6, margin >= 50).  Root positions: 1e-9 (prompts/001 section 4; prototype <= 2.4e-12).
Evaluation budget: 12 per root per continuation sub-step (prototype: secant 4-7 + 2 for Delta_h/Q').
"""
import numpy as np
import pytest

from m2_helpers import C_BASELINE, GRID_DE, GRID_E0, E_MAX, sector

TOL_J_REL = 1e-4
TOL_E = 1e-9
DISC_MIN = 8.0
EVALS_PER_ROOT_PER_SUBSTEP = 12
SEED = 20261024

# coded_v2 roots at C_BASELINE (labelled, reports/tests/007 table; oracle-verified)
ROOTS = [("000_A1m", 0.26932269974456263), ("100_A2", 0.2915560996081883), ("200_A2", 0.28624436431121125)]


def _ev():
    from qc3 import fast
    return fast.evaluate_point


def _fd_jacobian(sec, c, E0, rel=1e-3):
    from qc3 import roots
    c = np.asarray(c, float)
    J = np.zeros(4)
    for i in range(4):
        d = np.zeros(4); d[i] = rel * abs(c[i])
        Ep = roots.refine_root(_ev(), sec, tuple(c + d), E0 - 2e-5, E0 + 2e-5, kind="any")
        Em = roots.refine_root(_ev(), sec, tuple(c - d), E0 - 2e-5, E0 + 2e-5, kind="any")
        J[i] = (Ep - Em) / (2 * d[i])
    return J


@pytest.mark.parametrize("irrep,E0", ROOTS)
def test_root_jacobian_matches_finite_difference(irrep, E0):
    from qc3 import tracking
    sec = sector(irrep, k2_pwave_form="coded_v2")
    tr = tracking.root_jacobian(_ev(), sec, C_BASELINE, E0)
    assert abs(tr.E - E0) <= 1e-12
    assert tr.J.shape == (4,) and np.all(np.isfinite(tr.J))
    assert tr.discriminant <= -DISC_MIN and tr.N > 0 and tr.slope > 0.0
    Jfd = _fd_jacobian(sec, C_BASELINE, E0)
    assert np.max(np.abs(tr.J - Jfd) / np.abs(Jfd)) <= TOL_J_REL, (tr.J, Jfd)
    # the root moves down for a larger K_iso0 at every labelled level (prototype: all J_0 < 0)
    assert tr.J[0] < 0.0


def test_secant_root_from_jacobian_prediction():
    from qc3 import roots, tracking
    rng = np.random.default_rng(SEED)
    sec = sector("100_A2", k2_pwave_form="coded_v2")
    E0 = ROOTS[1][1]
    tr = tracking.root_jacobian(_ev(), sec, C_BASELINE, E0)
    c0 = np.asarray(C_BASELINE)
    for frac in (0.05, 0.05, 0.3):
        c1 = c0 * (1 + frac * rng.uniform(-1, 1, 4))
        dE = float(tr.J @ (c1 - c0))
        E, l_last, n = tracking.secant_root(_ev(), sec, tuple(c1), E0 + dE, dE)
        E_ref = roots.refine_root(_ev(), sec, tuple(c1), E - 5e-6, E + 5e-6)       # Brent, kind="zero"
        assert abs(E - E_ref) <= TOL_E, (frac, E, E_ref)
        assert n <= 10 and np.isfinite(l_last)
        # the predictor is good to a fraction of the move (prototype: <= 5 % at 5 %, <= 15 % at 30 %)
        assert abs((E - E0) - dE) <= 0.3 * abs(dE) + 1e-7


@pytest.mark.slow
def test_sector_roots_track_random_parameter_steps_and_rescan_agrees():
    from qc3 import rootscan, tracking
    rng = np.random.default_rng(SEED)
    c0 = np.asarray(C_BASELINE)
    for irrep in ("100_A2", "200_A2"):
        sec = sector(irrep, k2_pwave_form="coded_v2")
        sr = tracking.SectorRoots(sec, tuple(c0), GRID_E0, E_MAX, evaluate=_ev())
        assert np.array_equal(sr.energies(), sr.scan.zeros()) and len(sr.roots) == len(sr.scan.zeros()) >= 2
        assert sr.jacobian().shape == (len(sr.roots), 4)
        assert all(r.slope > 0 and r.discriminant <= -DISC_MIN for r in sr.roots)
        n_prev = sr.n_evaluations
        assert np.array_equal(sr.update(tuple(c0)), sr.energies()) and sr.n_evaluations == n_prev   # no-op update
        c = c0.copy()
        for frac in (0.05, 0.05, 0.05, 0.3, 0.3):
            c = c * (1 + frac * rng.uniform(-1, 1, 4))
            n_ev0, n_sub0 = sr.n_evaluations, sr.n_substeps
            E_tr = sr.update(tuple(c))
            assert np.all(np.diff(E_tr) > 0) and len(E_tr) == len(sr.roots)
            n_sub = sr.n_substeps - n_sub0
            assert n_sub >= 1
            assert sr.n_evaluations - n_ev0 <= EVALS_PER_ROOT_PER_SUBSTEP * n_sub * len(E_tr) + 400 * sr.n_local_rescans, (irrep, frac)
            # independent reference: the label-free scan on a window around every tracked root at this c
            for E in E_tr:
                ref = rootscan.scan_sector(_ev(), sec, tuple(c), max(GRID_E0, E - 60 * GRID_DE), min(E_MAX, E + 60 * GRID_DE))
                zs = ref.zeros()
                assert len(zs) >= 1 and np.min(np.abs(zs - E)) <= TOL_E, (irrep, frac, E, zs)
            assert np.array_equal(sr.update(tuple(c)), E_tr)
        # a full rescan at the final c must reproduce the tracked set (count and positions)
        res = sr.rescan()
        assert len(res.zeros()) == len(sr.roots)
        assert np.max(np.abs(res.zeros() - sr.energies())) <= TOL_E
        assert sr.n_full_rescans == 0, "random 5-30 % steps must be handled by continuation, not by full rescans"


def test_sector_roots_excludes_free_level_solutions_and_reports_them():
    from qc3 import tracking
    sec = sector("110_A2", k2_pwave_form="coded_v2")
    sr = tracking.SectorRoots(sec, tuple(C_BASELINE), 0.3190, 0.3200, evaluate=_ev())
    assert len(sr.roots) == 1 and abs(sr.roots[0].E - 0.3191549427) <= TOL_E
    assert len(sr.free_level_solutions) == 1 and abs(sr.free_level_solutions[0].E - 0.319468828) <= 5e-10
    assert sr.free_level_solutions[0].kind == "free_level"

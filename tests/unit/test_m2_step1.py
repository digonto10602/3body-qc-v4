"""Unit tests for M2 step 1 (prompts/003): qc3.free_levels, the interface stubs, SingularQError.  Fixed inputs only."""
import itertools

import numpy as np
import pytest

from qc3 import free_levels as fl
from qc3.point import SingularQError
from qc3.sector import Sector

MK, MPI, XI = 0.09698, 0.06906, 3.444


def _naive(d, Lbyas, E_min, E_max, n=3):
    """Loop-only enumeration of unordered-kaon states; returns {rounded E*: [states]} (independent of qc3.free_levels)."""
    q = 2 * np.pi / (XI * Lbyas)
    out = {}
    pts = list(itertools.product(range(-n, n + 1), repeat=3))
    for a in pts:
        for b in pts:
            if b < a:
                continue
            c = tuple(d[i] - a[i] - b[i] for i in range(3))
            El = sum(np.sqrt(m * m + q * q * sum(x * x for x in v)) for m, v in ((MK, a), (MK, b), (MPI, c)))
            Es = np.sqrt(El * El - q * q * sum(x * x for x in d))
            if E_min <= Es <= E_max:
                out.setdefault(round(Es, 9), []).append((a, b, c))
    return out


@pytest.mark.parametrize("irrep,d", [("000_A1m", (0, 0, 0)), ("100_A2", (0, 0, 1)), ("110_A2", (1, 1, 0)), ("111_A2", (1, 1, 1)),
                                     ("200_A2", (0, 0, 2))])
def test_small_box_matches_naive_enumeration(irrep, d):
    Lbyas = 8                                     # large q: few levels, small enumeration cube
    sec = Sector.production(irrep, Lbyas=Lbyas)
    lv = fl.free_level_states(sec, 0.0, 0.70)
    ref = _naive(d, Lbyas, 0.0, 0.70)
    assert len(lv) == len(ref) >= 3
    for level, key in zip(lv, sorted(ref)):
        assert abs(level.E - key) < 2e-9
        assert level.multiplicity == len(ref[key]) == len(level.states)
        assert set(level.states) == set(ref[key])
    assert np.all(np.diff(fl.free_levels(sec, 0.0, 0.70)) > 1e-11)


def test_threshold_and_first_moving_level_closed_form():
    q = 2 * np.pi / (XI * 20)
    s0 = Sector.production("000_A1m", Lbyas=20)
    lv0 = fl.free_level_states(s0, 0.0, 0.40)
    assert lv0[0].E == pytest.approx(2 * MK + MPI, abs=1e-15)
    # next rest-frame levels with |n|^2 = 1: K(0)K(e)pi(-e) (6 choices of e) has E = mK + w_K + w_pi;
    # K(e)K(-e)pi(0) (3 axes, unordered pair) has E = 2 w_K + mpi
    E_a = MK + np.sqrt(MK ** 2 + q * q) + np.sqrt(MPI ** 2 + q * q)
    E_b = 2 * np.sqrt(MK ** 2 + q * q) + MPI
    levels = {round(l.E, 12): l for l in lv0}
    assert levels[round(E_a, 12)].multiplicity == 6
    assert levels[round(E_b, 12)].multiplicity == 3
    # one-axis boost: d = (0,0,1); lowest level has one kaon carrying the momentum (the lighter pion would cost more), E* = sqrt(E_lab^2 - q^2)
    s1 = Sector.production("100_A2", Lbyas=20)
    E_lab = MK + np.sqrt(MK ** 2 + q * q) + MPI
    assert fl.free_levels(s1, 0.0, 0.5)[0] == pytest.approx(np.sqrt(E_lab ** 2 - q * q), abs=1e-15)


def test_range_is_inclusive_and_degenerate_ranges_are_empty():
    s = Sector.production("110_A2", Lbyas=20)
    E = fl.free_levels(s, 0.2, 0.36)
    assert len(E) > 2
    assert np.array_equal(fl.free_levels(s, E[1], E[2]), E[1:3])           # both ends included
    assert len(fl.free_levels(s, E[1] + 1e-9, E[2] - 1e-9)) == 0
    assert len(fl.free_levels(s, 0.36, 0.2)) == 0                         # reversed range
    assert len(fl.free_levels(s, 0.0, 0.1)) == 0                          # below threshold
    assert fl.free_levels(s, 0.2, 0.36).dtype == np.float64


def test_states_are_canonical_sums_and_energies_recomputed():
    s = Sector.production("111_A2", Lbyas=24)
    q = 2 * np.pi / (XI * 24)
    d = np.array(s.nnP)
    for lv in fl.free_level_states(s, 0.2, 0.36):
        assert lv.multiplicity == len(lv.states) >= 1 and list(lv.states) == sorted(lv.states)
        for n1, n2, n3 in lv.states:
            assert n1 <= n2
            assert tuple(np.array(n1) + np.array(n2) + np.array(n3)) == tuple(d)
            El = sum(np.sqrt(m * m + q * q * float(np.dot(v, v))) for m, v in ((MK, n1), (MK, n2), (MPI, n3)))
            assert np.sqrt(El ** 2 - q * q * float(d @ d)) == pytest.approx(lv.E, abs=1e-11)


def test_identical_kaons_are_not_double_counted():
    # rest frame, K(0)K(0)pi(0) has multiplicity 1; K(e1)K(e2)pi(-e1-e2) with e1 != e2 appears once, not twice
    s = Sector.production("000_A1m", Lbyas=20)
    st = fl.free_level_states(s, 0.0, 0.40)
    n_states = sum(l.multiplicity for l in st)
    all_states = [x for l in st for x in l.states]
    assert len(all_states) == len(set(all_states)) == n_states
    for n1, n2, _ in all_states:
        assert n1 <= n2


def test_interface_stubs_raise_not_implemented_and_modules_import():
    # rootscan.scan_sector / read_labels / score_labels: step 2 (tests/unit/test_rootscan.py); qc3.tracking: step 4 (tests/unit/test_tracking.py);
    # fit.QCChi2 / make_trackers: step 5 (tests/unit/test_m2_step5_chi2.py); the rest of qc3.fit stays a stub until step 6
    from qc3 import fit, rootscan, tracking
    s = Sector.production("000_A1m")
    calls = [
        (fit.fit_cutoff, (None, {}, (0, 0, 0, 0))),
        (fit.verify_fit, (None, {}, None)),
        (fit.cutoff_scan, (None, 20, (0.335,), (0, 0, 0, 0), "q0_consistent")),
        (fit.save_results, ([], "x.json")),
        (fit.main, ([],)),
    ]
    for fn, args in calls:
        with pytest.raises(NotImplementedError, match=r"prompts/003 step \d"):
            fn(*args)
    assert issubclass(tracking.RootSetChangedError, RuntimeError)
    assert set(rootscan.KINDS) == {"zero", "pole", "even_zero", "free_level", "ambiguous"}
    assert rootscan.SCAN_DEFAULTS["free_tol"] == 1e-7 and rootscan.SCAN_DEFAULTS["h2"] == 3e-6
    assert tracking.TRACK_DEFAULTS["local_steps_rows"] == (4.0, 1.0, 0.25)
    assert fit.QCChi2.errordef == 1.0 and len(fit.C_V32W) == 4


def test_label_files_exist_and_are_relative_to_the_data_root():
    from qc3 import lattice, rootscan
    assert len(rootscan.LABEL_FILES) == 7
    for (L, irrep), rel in rootscan.LABEL_FILES.items():
        p = lattice.default_data_dir() / rel
        assert p.is_file(), p
        assert f"L{L}_{irrep}" in rel


def test_fast_maps_exact_zero_pivot_to_SingularQError(monkeypatch):
    from qc3 import _core, fast

    def boom(*a, **k):
        raise ValueError("fast: Q is singular (zero pivot in LU)")
    monkeypatch.setattr(_core, "fast_point", boom)
    with pytest.raises(SingularQError, match="Q is singular"):
        fast.evaluate_point(Sector.production("000_A1m"), 0.3, (0.0, 0.0, 0.0, 0.0), want=("Q",))
    assert issubclass(SingularQError, ValueError)

    def other(*a, **k):
        raise ValueError("some other failure")
    monkeypatch.setattr(_core, "fast_point", other)
    with pytest.raises(ValueError, match="some other failure") as ei:
        fast.evaluate_point(Sector.production("000_A1m"), 0.3, (0.0, 0.0, 0.0, 0.0), want=("Q",))
    assert not isinstance(ei.value, SingularQError)

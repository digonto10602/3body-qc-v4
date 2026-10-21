"""M1 steps 6-7 — the cache-free fast path reproduces the oracle at every energy and irrep, and is faster.
Tolerances (prompts/001 section 4, reports/research/007 section 3.4): ln|det Q| 1e-6 absolute (LU round-off floor
1.4e-9 at kappa(F3)=4.8e8, margin 700); ingredients 1e-12 relative Frobenius (truncation bound 1e-14, margin 100);
F3inv 1e-9 (measured response 2e-11 to 1e-12 ingredient noise, margin 50); roots 1e-9 in a_t^{-1}."""
import time

import numpy as np
import pytest

from m1_helpers import (C_BASELINE, C_ZERO, IRREPS, PRODUCTION, accepted_windows, anchor_points, as_c, fast_point,
                        grid_energy, oracle_point, rel_frob, sector)

TOL_LNDET = 1e-6
TOL_ING = 1e-12
TOL_F3INV = 1e-9
TOL_ROOT = 1e-9
E_MIN, E_MAX = 0.26310, 0.36000


def _compare(irrep, Ecm, c):
    o = oracle_point(irrep, Ecm, c)
    f = fast_point(irrep, Ecm, c)
    tag = (irrep, Ecm, tuple(c))
    assert (f.N1, f.N2, f.N, f.Nproj) == (o.N1, o.N2, o.N, o.Nproj), tag
    for name in o.basis.dtype.names:
        assert np.array_equal(np.asarray(f.basis[name]), np.asarray(o.basis[name])), (tag, name)
    assert int(f.sign) == int(o.sign), (tag, f.sign, o.sign, f.lnabsdet, o.lnabsdet)
    assert abs(f.lnabsdet - o.lnabsdet) <= TOL_LNDET, (tag, f.lnabsdet, o.lnabsdet)
    for name in ("F", "G", "K2inv"):
        d = rel_frob(getattr(f, name), getattr(o, name))
        assert d <= TOL_ING, (tag, name, d)
    for name in ("iso0", "iso1", "B", "E"):
        d = rel_frob(f.K3_basis[name], o.K3_basis[name])
        assert d <= TOL_ING, (tag, "K3_" + name, d)
    d = rel_frob(f.F3inv, o.F3inv)
    assert d <= TOL_F3INV, (tag, "F3inv", d)
    # the projected matrices must agree too (V may differ by a unitary mixing inside degenerate subspaces,
    # so compare the eigenvalues of Q, which are invariant)
    eo = np.sort(np.linalg.eigvals(as_c(o.Q)).real)
    ef = np.sort(np.linalg.eigvals(as_c(f.Q)).real)
    assert np.allclose(ef, eo, rtol=1e-6, atol=1e-6 * np.abs(eo).max()), tag
    return o, f


def test_fast_vs_oracle_small_points():
    pts = anchor_points(small_only=True)
    assert len(pts) >= 3
    for p in pts:
        _compare(p["irrep"], p["Ecm"], C_ZERO)
        _compare(p["irrep"], p["Ecm"], C_BASELINE)


def test_fast_basis_identical_to_oracle():
    for p in anchor_points(small_only=True):
        o = oracle_point(p["irrep"], p["Ecm"], C_ZERO)
        f = fast_point(p["irrep"], p["Ecm"], C_ZERO)
        assert f.basis.dtype.names == o.basis.dtype.names
        assert len(f.basis) == len(o.basis) == o.N
        for name in o.basis.dtype.names:
            assert np.array_equal(np.asarray(f.basis[name]), np.asarray(o.basis[name])), name


def test_fast_metadata_and_timings():
    p = anchor_points(small_only=True)[0]
    o = oracle_point(p["irrep"], p["Ecm"], C_ZERO)
    f = fast_point(p["irrep"], p["Ecm"], C_ZERO)
    for k in PRODUCTION:
        assert f.metadata[k] == o.metadata[k], k
    assert f.metadata["backend"] == "fast" and f.metadata["threads"] == 1
    assert "trunc_eps" in f.metadata and 0 < float(f.metadata["trunc_eps"]) <= 1e-14
    assert "use_classes" in f.metadata
    keys = {"basis", "F", "G", "K2inv", "F3", "F3inv", "projector", "K3", "projection", "lu", "total"}
    assert keys <= set(f.timings) and f.timings["total"] > 0


@pytest.mark.slow
def test_fast_vs_oracle_all_anchor_points():
    pts = anchor_points()
    assert len(pts) == 12
    for p in pts:
        _compare(p["irrep"], p["Ecm"], C_BASELINE)


@pytest.mark.slow
def test_fast_vs_oracle_random_energies_all_irreps():
    """Generalisation check: seeded random energies over the whole production range, every irrep, live oracle."""
    rng = np.random.default_rng(20261021)
    for irrep in IRREPS:
        drawn = 0
        tries = 0
        while drawn < 2 and tries < 20:
            tries += 1
            Ecm = float(rng.uniform(E_MIN, E_MAX))
            # avoid the immediate neighbourhood of a dimension jump (N changes within +-1e-5)
            Nm = fast_point(irrep, Ecm - 1e-5, C_ZERO).N
            Np = fast_point(irrep, Ecm + 1e-5, C_ZERO).N
            if Nm != Np:
                continue
            _compare(irrep, Ecm, C_BASELINE)
            drawn += 1
        assert drawn == 2, irrep


@pytest.mark.slow
def test_fast_root_vs_oracle_100_A2_windows():
    from qc3 import fast, oracle
    from qc3.roots import refine_root
    rows = accepted_windows("100_A2")
    assert len(rows) == 2
    sec = sector("100_A2")
    for row in rows:
        E_left, E_right = grid_energy(row["row_left"]), grid_energy(row["row_right"])
        label = float(row["zero_estimate_initial"])
        assert E_left < label < E_right
        E_o = refine_root(oracle.evaluate_point, sec, C_BASELINE, E_left, E_right, xtol=1e-11)
        E_f = refine_root(fast.evaluate_point, sec, C_BASELINE, E_left, E_right, xtol=1e-11)
        assert E_left < E_o < E_right and E_left < E_f < E_right
        assert abs(E_f - E_o) <= TOL_ROOT, (row["bracket_id"], E_f, E_o)
        # the human label was linear interpolation on the 4.845e-6 grid: agreement to 5e-6
        assert abs(E_o - label) <= 5e-6, (row["bracket_id"], E_o, label)


@pytest.mark.slow
def test_fast_speedup_single_thread_ge_20x():
    from qc3 import fast, oracle
    sec = sector("110_A2")
    Ecm = 0.35
    fast.evaluate_point(sec, Ecm, C_BASELINE, options={"threads": 1})  # warm-up (fills in-memory caches)
    t_f = []
    for _ in range(3):
        t0 = time.perf_counter(); fast.evaluate_point(sec, Ecm, C_BASELINE, options={"threads": 1}); t_f.append(time.perf_counter() - t0)
    t_o = []
    for _ in range(3):
        t0 = time.perf_counter(); oracle.evaluate_point(sec, Ecm, C_BASELINE); t_o.append(time.perf_counter() - t0)
    ratio = np.median(t_o) / np.median(t_f)
    assert ratio >= 20.0, (np.median(t_o), np.median(t_f), ratio)

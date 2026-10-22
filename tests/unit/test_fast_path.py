"""Unit tests for the fast path v1 (prompts/001 step 6). Small points only."""
import numpy as np
import pytest

from qc3 import fast, oracle
from qc3.sector import Sector

POINTS = [("100_A2", 0.2915561), ("000_A1m", 0.26932), ("111_A2", 0.32036)]


def _rel(a, b):
    return np.linalg.norm(a - b) / np.linalg.norm(b)


@pytest.mark.parametrize("irrep,Ecm", POINTS)
def test_truncation_error_below_bound(irrep, Ecm):
    """Ellipsoid truncation at trunc_eps=1e-16 vs an (almost) untruncated sum (trunc_eps=1e-40, R^2 ~ x^2+190):
    the dropped terms change F by <= 1e-14 relative (plan section 2.3 rung 2)."""
    s = Sector.production(irrep)
    f_def = fast.evaluate_point(s, Ecm, want=("F",))
    f_ref = fast.evaluate_point(s, Ecm, want=("F",), options={"trunc_eps": 1e-40})
    assert f_ref.metadata["lattice_terms"] > f_def.metadata["lattice_terms"]
    assert _rel(f_def.F, f_ref.F) <= 1e-14
    # looser tolerance really truncates more and the error grows accordingly (bound is not vacuous)
    f_loose = fast.evaluate_point(s, Ecm, want=("F",), options={"trunc_eps": 1e-8})
    assert f_loose.metadata["lattice_terms"] < f_def.metadata["lattice_terms"]
    assert 1e-14 < _rel(f_loose.F, f_ref.F) <= 1e-6


def test_term_count_is_much_smaller_than_cube():
    f = fast.evaluate_point(Sector.production("100_A2"), 0.2915561, want=("F",))
    n_momenta = f.metadata["lattice_cube_terms_per_entry"] // 41 ** 3
    assert f.metadata["lattice_terms"] < 0.2 * n_momenta * 41 ** 3


def _proj(V):
    V = np.asarray(V)
    return V @ V.conj().T


@pytest.mark.parametrize("method", ["orbit", "eigen"])
def test_cached_projector_identical_to_cold_build(method):
    """The cache returns the stored V_sel unchanged (bit-identical hit vs miss).  Against the oracle's cold build the
    comparison is of the subspace V V^dag (prompts/002 7.2): the orbit construction differs from the eigenvectors by a
    unitary mixing inside each (flavour, orbit, l) block."""
    s = Sector.production("110_A2")
    Ecm = 0.30733
    opt = {"projector": method}
    o = oracle.evaluate_point(s, Ecm, want=("Vsel",))          # oracle clears the cache first: cold build
    fast.clear_caches()
    f1 = fast.evaluate_point(s, Ecm, want=("Vsel",), options=opt)   # miss (cold build inside the fast path)
    f2 = fast.evaluate_point(s, Ecm, want=("Vsel",), options=opt)   # hit
    assert not f1.metadata["projector_cache_hit"] and f2.metadata["projector_cache_hit"]
    assert f1.metadata["projector_method"] == method
    assert np.array_equal(f1.Vsel, f2.Vsel)
    assert f1.Nproj == o.Nproj
    assert _rel(_proj(f1.Vsel), _proj(o.Vsel)) <= 1e-12
    if method == "eigen":
        assert np.array_equal(o.Vsel, f1.Vsel)                 # same builder, same bits
    # a different energy with the same basis signature also hits and gives the same V_sel
    f3 = fast.evaluate_point(s, Ecm + 1e-5, want=("Vsel",), options=opt)
    if f3.N == o.N:
        assert f3.metadata["projector_cache_hit"] and np.array_equal(f1.Vsel, f3.Vsel)


@pytest.mark.parametrize("irrep,Ecm", POINTS)
def test_orbit_and_eigen_projectors_agree(irrep, Ecm):
    """prompts/002 7.2: equal N_proj and equal V V^dag to 1e-12 at the three unit points; V orthonormal."""
    s = Sector.production(irrep)
    fo = fast.evaluate_point(s, Ecm, want=("Vsel",), options={"projector": "orbit"})
    fe = fast.evaluate_point(s, Ecm, want=("Vsel",), options={"projector": "eigen"})
    assert fo.Nproj == fe.Nproj and fo.Vsel.shape == (fo.N, fo.Nproj)
    assert _rel(_proj(fo.Vsel), _proj(fe.Vsel)) <= 1e-12
    assert np.linalg.norm(fo.Vsel.conj().T @ fo.Vsel - np.eye(fo.Nproj)) <= 1e-12


def test_orbit_projector_independent_of_eig_tol():
    """The orbit construction does not use eig_tol (eigenvalues are exactly 0/1): 17 columns at 100_A2 0.2915561 for
    eig_tol 0.05 and 1e-30 (the eigen method gives 17 and 5; 5 is a round-off artefact), in either order."""
    Ecm = 0.2915561
    fast.clear_caches()
    n = [fast.evaluate_point(Sector.production("100_A2", eig_tol=t), Ecm, want=("Vsel",)).Nproj for t in (0.05, 1e-30, 0.3)]
    assert n == [17, 17, 17]
    assert fast.cache_info()["projector_entries"] == 3     # eig_tol stays in the key


def test_dimension_jumps_basic_properties():
    for irrep, n_expected_min in (("000_A1m", 3), ("100_A2", 8), ("110_A2", 8), ("111_A2", 8), ("200_A2", 8)):
        s = Sector.production(irrep)
        j = fast.dimension_jumps(s, 0.2631, 0.36)
        assert j.ndim == 1 and len(j) >= n_expected_min and np.all(np.diff(j) > 0) and np.all((j > 0.2631) & (j < 0.36))
        # a sub-interval gives the matching subset
        sub = fast.dimension_jumps(s, 0.29, 0.33)
        assert np.array_equal(sub, j[(j > 0.29) & (j < 0.33)])
        # N increases across every jump
        for e in j[:: max(1, len(j) // 3)]:
            lo = fast.evaluate_point(s, e - 2e-9, want=("Q",)).N
            hi = fast.evaluate_point(s, e + 2e-9, want=("Q",)).N
            assert hi > lo


@pytest.mark.parametrize("irrep", ["000_A1m", "100_A2", "110_A2", "111_A2", "200_A2"])
def test_threshold_basis_equals_oracle_basis(irrep):
    """prompts/002 7.3: the threshold builder's Kin (momenta as complex doubles, n, ell, m, N1, N2) is BIT-identical to
    the oracle's make_kinematics (config_maker_4_momentum_first) at 10 seeded energies per irrep over the whole
    production range, and at the exact jump energies +- 1 ulp."""
    from qc3 import _core
    from qc3.fast import _params
    rng = np.random.default_rng(20261022)
    s = Sector.production(irrep)
    jumps = fast.dimension_jumps(s, 0.2631, 0.36)
    energies = list(rng.uniform(0.2631, 0.36, size=10)) + list(jumps[:2]) + [float(np.nextafter(jumps[0], 1.0)), float(np.nextafter(jumps[0], 0.0))]
    for Ecm in energies:
        r = _core.compare_kinematics(_params(s, float(Ecm), (0.0, 0.0, 0.0, 0.0), False))
        assert r["identical"], (irrep, Ecm, r)
        assert r["N"] > 0


def test_options_and_metadata():
    s = Sector.production("000_A1m")
    r = fast.evaluate_point(s, 0.26932, options={"trunc_eps": 1e-15})
    assert r.metadata["backend"] == "fast" and r.metadata["trunc_eps"] == 1e-15
    assert r.metadata["use_classes"] is False and r.metadata["threads"] == 1
    with pytest.raises(NotImplementedError):
        fast.evaluate_point(s, 0.26932, options={"use_classes": True})
    with pytest.raises(NotImplementedError):
        fast.evaluate_point(s, 0.26932, options={"threads": 4})
    with pytest.raises(ValueError):
        fast.evaluate_point(s, 0.26932, options={"bogus": 1})
    with pytest.raises(ValueError):
        fast.evaluate_point(s, 0.26932, options={"trunc_eps": 0.0})


def test_both_pwave_forms_follow_oracle():
    for form in ("coded_v2", "q0_consistent"):
        s = Sector.production("100_A2", k2_pwave_form=form)
        o = oracle.evaluate_point(s, 0.2915561)
        f = fast.evaluate_point(s, 0.2915561)
        assert f.metadata["k2_pwave_form"] == form
        assert np.array_equal(f.K2inv, o.K2inv) and _rel(f.F, o.F) < 1e-12
        assert f.sign == o.sign and abs(f.lnabsdet - o.lnabsdet) < 1e-6


def test_projector_cache_key_includes_eig_tol():
    """prompts/002 7.0 (report 006 note 2): eig_tol=0.05 then 1e-30 in one process gives Nproj 17 then 5 at 100_A2
    0.2915561, in either order, and the fast-owned cache holds two entries."""
    Ecm = 0.2915561
    fast.clear_caches()

    def nproj(tol):
        return fast.evaluate_point(Sector.production("100_A2", eig_tol=tol), Ecm, want=("Vsel",),
                                   options={"projector": "eigen"}).Nproj

    base = fast.cache_info()                       # counters are cumulative; clear_caches only drops the entries
    assert nproj(0.05) == 17 and nproj(1e-30) == 5
    info = fast.cache_info()
    assert info["projector_entries"] == 2 and info["projector_misses"] - base["projector_misses"] == 2
    assert nproj(0.05) == 17 and fast.cache_info()["projector_hits"] - base["projector_hits"] == 1
    fast.clear_caches()
    assert nproj(1e-30) == 5 and nproj(0.05) == 17 and fast.cache_info()["projector_entries"] == 2
    fast.clear_caches()
    assert fast.cache_info()["projector_entries"] == 0


def test_projector_option_values():
    s = Sector.production("000_A1m")
    assert fast.evaluate_point(s, 0.26932, want=("Q",)).metadata["projector_method"] == "orbit"   # default
    with pytest.raises(ValueError):
        fast.evaluate_point(s, 0.26932, options={"projector": "bogus"})

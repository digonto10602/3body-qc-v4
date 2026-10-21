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


def test_cached_projector_identical_to_cold_build():
    s = Sector.production("110_A2")
    Ecm = 0.30733
    o = oracle.evaluate_point(s, Ecm, want=("Vsel",))          # oracle clears the cache first: cold build
    fast.clear_caches()
    f1 = fast.evaluate_point(s, Ecm, want=("Vsel",))           # miss (cold build inside the fast path)
    f2 = fast.evaluate_point(s, Ecm, want=("Vsel",))           # hit
    assert not f1.metadata["projector_cache_hit"] and f2.metadata["projector_cache_hit"]
    assert np.array_equal(o.Vsel, f1.Vsel) and np.array_equal(o.Vsel, f2.Vsel)
    # a different energy with the same basis signature also hits and gives the same V_sel
    f3 = fast.evaluate_point(s, Ecm + 1e-5, want=("Vsel",))
    if f3.N == o.N:
        assert f3.metadata["projector_cache_hit"] and np.array_equal(o.Vsel, f3.Vsel)


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

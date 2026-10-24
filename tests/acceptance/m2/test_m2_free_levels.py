"""M2 step 1 (prompts/003 section 2.1, D8): non-interacting K K pi levels per sector.

The production enumerator `qc3.free_levels` is compared with an independent brute-force enumeration written in
m2_helpers (same kinematics: E_lab = sum_a sqrt(m_a^2 + q^2 n_a^2), sum n_a = d, E* = sqrt(E_lab^2 - P^2); identical
kaons -> unordered pairs).  The agreement is pure float64 arithmetic on three square roots, so 1e-12 is 10^3 above
the round-off (a few 1e-16 relative on E* ~ 0.3).
"""
import numpy as np
import pytest

from m2_helpers import ATMK, ATMPI, E_MAX, FREE_LEVEL_110_A2, GRID_E0, IRREPS, brute_force_free_levels, sector

TOL_E = 1e-12


def _states_as_set(states):
    """states: iterable of (n_K1, n_K2, n_pi) triples of int 3-tuples -> frozenset for order-free comparison."""
    return frozenset(tuple(tuple(int(x) for x in v) for v in st) for st in states)


@pytest.mark.parametrize("Lbyas", [20, 24])
@pytest.mark.parametrize("irrep", IRREPS)
def test_free_levels_match_brute_force(irrep, Lbyas):
    from qc3 import free_levels as fl
    sec = sector(irrep, Lbyas=Lbyas, k2_pwave_form="q0_consistent")
    E_prod = fl.free_levels(sec, GRID_E0, E_MAX)
    E_ref, states_ref = brute_force_free_levels(irrep, Lbyas, GRID_E0, E_MAX)
    assert E_prod.dtype == np.float64 and E_prod.ndim == 1
    assert np.all(np.diff(E_prod) > 0), "free_levels must be sorted and distinct"
    assert len(E_prod) == len(E_ref), (irrep, Lbyas, E_prod, E_ref)
    assert np.max(np.abs(E_prod - E_ref)) <= TOL_E
    levels = fl.free_level_states(sec, GRID_E0, E_MAX)
    assert [lv.E for lv in levels] == pytest.approx(list(E_ref), abs=TOL_E)
    for lv, st_ref in zip(levels, states_ref):
        assert lv.multiplicity == len(st_ref) == len(lv.states), (irrep, Lbyas, lv.E)
        assert _states_as_set(lv.states) == _states_as_set(st_ref), (irrep, Lbyas, lv.E)
        for n1, n2, n3 in lv.states:
            assert tuple(n1) <= tuple(n2), "identical kaons: canonical order n_K1 <= n_K2 (lexicographic)"
            assert tuple(int(a) + int(b) + int(c) for a, b, c in zip(n1, n2, n3)) == tuple(sec.nnP)


def test_threshold_and_the_110_A2_quartic_free_level():
    from qc3 import free_levels as fl
    # rest frame: the lowest free level is the three-particle threshold 2 m_K + m_pi (all momenta zero)
    sec0 = sector("000_A1m", k2_pwave_form="q0_consistent")
    E0 = fl.free_levels(sec0, 0.2, 0.30)
    assert abs(E0[0] - (2 * ATMK + ATMPI)) <= 1e-14
    lv0 = fl.free_level_states(sec0, 0.2, 0.30)[0]
    assert lv0.multiplicity == 1 and _states_as_set(lv0.states) == _states_as_set([((0, 0, 0), (0, 0, 0), (0, 0, 0))])
    # 110_A2: the level at which reports/tests/007 found |det Q| ~ (E - E*)^4 without a sign change (D8)
    sec = sector("110_A2", k2_pwave_form="q0_consistent")
    levels = fl.free_level_states(sec, 0.3194, 0.3195)
    assert len(levels) == 1
    lv = levels[0]
    assert abs(lv.E - FREE_LEVEL_110_A2) <= 5e-10          # the report quotes 9 digits
    assert lv.multiplicity == 2
    assert _states_as_set(lv.states) == _states_as_set([((0, 0, 0), (0, 1, 0), (1, 0, 0)), ((0, 0, 0), (1, 0, 0), (0, 1, 0))])
    # the range is inclusive and empty ranges are allowed
    assert len(fl.free_levels(sec, lv.E, lv.E)) == 1
    assert len(fl.free_levels(sec, 0.3195, 0.3196)) == 0

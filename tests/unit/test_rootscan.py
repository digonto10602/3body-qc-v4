"""Unit tests for qc3.rootscan (M2 step 2, prompts/003 section 2.2): synthetic det Q(E) with known zeros, poles, even-order
zeros, smooth extrema, a zero-pole pair inside one coarse step and a dimension jump.  Fixed inputs, no randomness."""
from types import SimpleNamespace

import numpy as np
import pytest

from qc3 import free_levels, rootscan
from qc3.constants import GRID_DE
from qc3.roots import DimensionJumpError
from qc3.sector import Sector

SEC = Sector.production("000_A1m", Lbyas=20, k2_pwave_form="coded_v2")
ROWS = GRID_DE


class Fake:
    """det Q(E) = s0 prod (E-z)^m / prod (E-p) * exp(bg(E)) with Q = diag(det, 1) (so slope = d det/dE at a zero, u = e1).
    Across each jump the dimension N increases by one and the sign flips (a dimension change, not a root)."""

    def __init__(self, zeros=(), poles=(), s0=1.0, bg=None, jumps=(), offset=-3.0):
        self.zeros, self.poles, self.s0, self.bg, self.jumps, self.offset = list(zeros), list(poles), s0, bg, list(jumps), offset
        self.n = 0

    def dimension_jumps(self, sector, a, b):
        return np.array([j for j in self.jumps if a < j < b])

    def __call__(self, sector, E, c, want=("Q",)):
        self.n += 1
        ln = self.offset + (self.bg(E) if self.bg else 0.0)
        sg = self.s0
        for z, m in self.zeros:
            ln += m * np.log(max(abs(E - z), 1e-300))
            sg *= np.sign(E - z) ** m
        for p in self.poles:
            ln -= np.log(max(abs(E - p), 1e-300))
            sg *= np.sign(E - p)
        k = sum(E > j for j in self.jumps)
        sg *= (-1) ** k
        g = sg * np.exp(min(ln, 700.0))
        return SimpleNamespace(sign=int(sg), lnabsdet=float(ln), N=30 + k, Q=np.diag([g, 1.0]).astype(complex))


def scan(f, a=0.27, b=0.30, **opts):
    return rootscan.scan_sector(f, SEC, (0.0, 0.0, 0.0, 0.0), a, b, options=opts or None)


def test_simple_zeros_poles_and_orders():
    f = Fake(zeros=[(0.2800123, 1), (0.2911111, 1)], poles=[0.2851234], s0=1.0)
    r = scan(f)
    zs = r.zeros()
    assert len(zs) == 2 and np.allclose(zs, [0.2800123, 0.2911111], atol=1e-10)
    assert np.allclose(r.poles(), [0.2851234], atol=1e-10)
    for c in r.of_kind("zero"):
        assert c.order_int == 1 and abs(c.order - 1.0) < 1e-3 and c.sign_change and c.via == "sign"
        assert c.discriminant <= -8.0 and c.N == 30 and c.segment == (0.27, 0.30)
    p = r.of_kind("pole")[0]
    assert p.discriminant >= 8.0 and p.slope is None and p.order is None
    assert r.n_evaluations == f.n and r.n_grid == len(r.grid) and r.n_grid < r.n_evaluations
    assert np.all(np.diff(r.grid) > 0) and np.max(np.diff(r.grid)) <= 4.0 * ROWS * (1 + 1e-9)
    assert [c.E for c in r.candidates] == sorted(c.E for c in r.candidates)


def test_slope_sign_is_recorded_not_used_to_drop():
    # det = +(E - z) rises through z (slope > 0); det = -(E - z) falls through z (slope < 0): both are kept as zeros (D9)
    r = scan(Fake(zeros=[(0.2812345, 1)], s0=1.0))
    r2 = scan(Fake(zeros=[(0.2812345, 1)], s0=-1.0))
    assert len(r.zeros()) == 1 and r.of_kind("zero")[0].slope > 0
    assert len(r2.zeros()) == 1 and r2.of_kind("zero")[0].slope < 0
    s = r.of_kind("zero")[0].slope
    assert abs(s - np.exp(-3.0)) / np.exp(-3.0) < 1e-3       # d/dE of exp(offset) (E - z) = exp(offset)


@pytest.mark.parametrize("m,expected_kind", [(2, "even_zero"), (4, "even_zero")])
def test_even_order_zero_without_sign_change_is_found_by_the_minimum_rule(m, expected_kind):
    f = Fake(zeros=[(0.2833337, m)], s0=1.0)
    r = scan(f)
    assert len(r.candidates) == 1 and len(r.zeros()) == 0
    c = r.candidates[0]
    assert c.kind == expected_kind and c.order_int == m and abs(c.order - m) < 1e-3 and not c.sign_change and c.via == "min"
    assert abs(c.E - 0.2833337) < 1e-9 and abs(c.slope) < 1e-9 and c.discriminant < -8.0


def test_even_zero_at_a_free_level_is_a_free_level_solution_and_not_a_zero():
    Ef = float(free_levels.free_levels(SEC, 0.27, 0.36)[0])
    f = Fake(zeros=[(Ef + 3e-9, 2), (Ef - 2.5e-3, 1)], s0=1.0)
    r = scan(f, 0.27, 0.36)
    fl = r.of_kind("free_level")
    assert len(fl) == 1 and fl[0].free_level == Ef and fl[0].E == Ef and fl[0].order_int == 2 and not fl[0].sign_change
    assert len(r.zeros()) == 1 and abs(r.zeros()[0] - (Ef - 2.5e-3)) < 1e-10      # the genuine zero 2.5e-3 away stays a zero
    assert not r.of_kind("even_zero") and not r.of_kind("ambiguous")


def test_free_level_probe_finds_a_double_zero_the_grid_rules_cannot_see():
    # a double zero 1 row inside a smooth plateau of the coarse grid is still found by the minimum rule; with the grid
    # minimum rule switched off by an equal-ln background the probe at the free level must find it anyway
    Ef = float(free_levels.free_levels(SEC, 0.27, 0.36)[1])
    f = Fake(zeros=[(Ef, 2)], s0=1.0)
    r = scan(f, 0.27, 0.36)
    assert [c.kind for c in r.candidates] == ["free_level"] and r.candidates[0].via in ("min", "free")


def test_smooth_minimum_and_maximum_yield_no_candidate_but_the_maximum_triggers_a_fine_rescan():
    bg_min = lambda E: 3.0 * ((E - 0.2871) / 0.004) ** 2
    bg_max = lambda E: -3.0 * ((E - 0.2871) / 0.004) ** 2
    base = scan(Fake(zeros=[]), 0.28, 0.294)
    rmin = scan(Fake(bg=bg_min), 0.28, 0.294)
    rmax = scan(Fake(bg=bg_max), 0.28, 0.294)
    assert rmin.candidates == [] and rmax.candidates == []
    fine = 2 * 4 / 0.25 + 1     # the maximum's window: +-1 coarse step at the fine step
    assert rmax.n_evaluations - base.n_evaluations == pytest.approx(fine, abs=2)
    assert rmin.n_evaluations > base.n_evaluations      # the bounded minimiser ran and the smooth minimum was discarded


@pytest.mark.parametrize("shift", [0.0, 0.7, 1.3, 2.1, 3.0, 3.6])
def test_zero_pole_pair_inside_one_coarse_step(shift):
    # pole at p, zero at p + 1.7 rows (s = |V^T u|^2 small): with a 4-row grid most phases see no sign change at all
    a = 0.27 + shift * ROWS
    p = 0.2850000
    z = p + 1.7 * ROWS
    f = Fake(zeros=[(z, 1)], poles=[p], s0=1.0, bg=lambda E: 0.0)
    r = scan(f, a, 0.30)
    near = [c for c in r.candidates if 0.2849 < c.E < 0.2852]
    assert sorted(c.kind for c in near) == ["pole", "zero"], [(c.kind, c.E, c.via) for c in near]
    assert abs(r.zeros()[0] - z) < 1e-9 and abs(r.poles()[0] - p) < 1e-9
    straddles = not np.any((r.grid > p) & (r.grid < z))
    if straddles:
        assert {c.via for c in near} & {"min", "max"}, "a straddled pair can only be found through the extrema rules"


def test_dimension_jump_is_a_segment_boundary_not_a_root():
    f = Fake(zeros=[(0.2750001, 1), (0.2950001, 1)], jumps=[0.2850001], s0=1.0)
    r = scan(f)
    assert np.array_equal(r.jumps, [0.2850001])
    assert r.segments == [(0.27, 0.2850001), (0.2850001, 0.30)]
    assert len(r.zeros()) == 2 and not r.poles().size                       # the sign flip across the jump is not a candidate
    assert {c.N for c in r.candidates} == {30, 31}
    for c in r.candidates:
        assert c.segment[0] <= c.E <= c.segment[1]
    assert np.all(np.abs(r.grid - 0.2850001) >= 2e-6 * (1 - 1e-9))        # edge guard


def test_zero_next_to_a_pole_is_resolved_by_the_fine_rescan_and_a_very_close_one_is_loud():
    p = 0.2850004
    z = p - 8.0e-6                       # 1.65 rows below the pole: resolved, order estimate unpolluted (h2 = 3e-6)
    r = scan(Fake(zeros=[(z, 1)], poles=[p]), 0.27, 0.30)
    assert abs(r.zeros()[0] - z) < 1e-9 and abs(r.poles()[0] - p) < 1e-9 and r.of_kind("zero")[0].order_int == 1
    # 3.8e-6 below the pole (the 200_A2 configuration): the pole pollutes the order estimate (m = 1.41, prompts/003 item 5);
    # the candidate is "ambiguous" (loud), never a silent zero
    z2 = p - 3.8e-6
    r2 = scan(Fake(zeros=[(z2, 1)], poles=[p]), 0.27, 0.30)
    assert [c.kind for c in r2.candidates] == ["ambiguous", "pole"] and not r2.zeros().size
    assert abs(r2.candidates[0].E - z2) < 1e-9 and abs(r2.candidates[1].E - p) < 1e-9


def test_options_are_validated():
    with pytest.raises(ValueError, match="unknown options"):
        scan(Fake(), nonsense=1)
    with pytest.raises(ValueError, match="xtol"):
        scan(Fake(), xtol=1e-6)
    with pytest.raises(ValueError):
        rootscan.scan_sector(Fake(), SEC, (0, 0, 0, 0), 0.3, 0.27)
    r = scan(Fake(), step_rows=8.0)
    assert r.options["step_rows"] == 8.0 and r.options["fine_step_rows"] == 0.25
    with pytest.raises(ValueError):
        r.of_kind("nonsense")


def test_nonconstant_dimension_inside_a_segment_is_loud():
    f = Fake(zeros=[(0.2750001, 1)], jumps=[0.2850001])
    f.dimension_jumps = lambda sector, a, b: np.array([])              # hide the jump: N changes inside one segment
    with pytest.raises(DimensionJumpError):
        scan(f)


def _cand(E, kind="zero"):
    return rootscan.RootCandidate(E=E, kind=kind, order=1.0, order_int=1, sign_change=True, discriminant=-20.0, slope=1.0, N=1,
                                  segment=(0, 1), bracket=(E, E), via="sign", free_level=None, n_evaluations=1)


def test_score_labels_counts_and_exclusions():
    labels = [dict(bracket_id=1, E_left=0.30, E_right=0.30 + 4 * ROWS, label="true"),
              dict(bracket_id=2, E_left=0.31, E_right=0.31 + 4 * ROWS, label="false"),
              dict(bracket_id=3, E_left=0.32, E_right=0.32 + 4 * ROWS, label="true"),
              dict(bracket_id=4, E_left=0.33, E_right=0.33 + 4 * ROWS, label="false"),
              dict(bracket_id=5, E_left=0.34, E_right=0.34 + 4 * ROWS, label="true")]
    cands = [_cand(0.30 + 2 * ROWS), _cand(0.33 + ROWS), _cand(0.35), _cand(0.30 - 0.2 * ROWS), _cand(0.31, "pole"),
             _cand(0.34 + 2 * ROWS, "free_level")]
    res = SimpleNamespace(candidates=cands)
    sc = rootscan.score_labels(res, labels, slack_rows=0.25)
    # bracket 1 holds two zeros (one in the slack margin) -> FN; 2 TN; 3 empty -> FN; 4 holds one -> FP; 5 only a free_level -> FN
    assert (sc["tp"], sc["fp"], sc["fn"], sc["tn"]) == (0, 1, 3, 1)
    assert [m["type"] for m in sc["mismatches"]] == ["FN", "FN", "FP", "FN"] or len(sc["mismatches"]) == 4
    assert [z.E for z in sc["unlabelled_zeros"]] == [0.35]
    sc2 = rootscan.score_labels(res, labels, slack_rows=0.1, exclude_ids=(3,))
    assert sc2["excluded"] == [3] and sc2["tp"] == 1 and sc2["fn"] == 1
    with pytest.raises(ValueError, match="unlabelled"):
        rootscan.score_labels(res, [dict(bracket_id=9, E_left=0.3, E_right=0.31, label="")])


def test_read_labels_relative_and_absolute_paths():
    from qc3.lattice import default_data_dir
    rows = rootscan.read_labels(rootscan.LABEL_FILES[(20, "000_A1m")])
    assert len(rows) == 10 and sum(r["label"] == "true" for r in rows) == 3
    assert rows[0].keys() == {"bracket_id", "E_left", "E_right", "label"} and all(r["E_left"] < r["E_right"] for r in rows)
    rows2 = rootscan.read_labels(str(default_data_dir() / rootscan.LABEL_FILES[(20, "000_A1m")]))
    assert rows == rows2
    assert set(rootscan.LABEL_FILES) == {(20, "000_A1m"), (20, "100_A2"), (20, "110_A2"), (20, "111_A2"), (20, "200_A2"),
                                         (24, "000_A1m"), (24, "100_A2")}
    for key in rootscan.LABEL_FILES:
        assert (default_data_dir() / rootscan.LABEL_FILES[key]).is_file()

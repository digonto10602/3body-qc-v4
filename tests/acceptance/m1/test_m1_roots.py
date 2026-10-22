"""M1 step 7.1 (prompts/002 section 2.5): the root refiner must return ZEROS of det Q, not poles.

Facts used (results/k2_pwave_shift_L20.json, step 9, fast backend, q0_consistent, c = baseline):
- Brent on the sign of det Q converged to a POLE in the widened windows of 110_A2 level 1 and 111_A2 level 0;
- 100_A2 level 0, 110_A2 level 0 and 200_A2 level 2 have a pole+zero pair inside the widened window, so the window
  ends have equal signs and the zero is invisible to a single bracket.
Discriminant: Delta_h = l(E*) - [l(E*-h) + l(E*+h)]/2 with l = ln|det Q|, h = 1e-6: zero <=> Delta_h <= -4, pole <=> >= +4.
Exact-arithmetic magnitude ln(h/xtol) = 11.5 for simple zeros/poles; measured 11.6-21.6 (zeros), 11.9-12.7 (poles).
Energies in a_t^{-1}; windows are derived from data/windows/L20/*/accepted_windows.csv (provenance in docs/PROVENANCE.tsv).
"""
import numpy as np
import pytest

from m1_helpers import C_BASELINE, GRID_DE, accepted_windows, sector

WIDEN_ROWS = 250                 # as in scripts/m1/k2_pwave_shift.py
DISC_MIN = 8.0                   # |Delta_h| at a clean zero/pole; exact bound 11.5, measured >= 11.6 (margin 1.45x on the bound)
TOL_E = 1e-8                     # same backend, Brent xtol 1e-11; 500x below the grid step 4.845e-6
TOL_LABEL = 5e-6                 # label resolution (prompts/001 section 4)

# step-9 values (results/k2_pwave_shift_L20.json, fast backend, q0_consistent, c = baseline)
STEP9_POINTS = [
    # irrep, level, kind, E
    ("110_A2", 1, "pole", 0.31612911112995534),
    ("110_A2", 1, "zero", 0.31564501983275506),
    ("111_A2", 0, "pole", 0.3195543911709816),
    ("111_A2", 0, "zero", 0.3192435345105379),
    ("100_A2", 0, "zero", 0.2914252169900899),
    ("110_A2", 0, "zero", 0.3066350722070048),
    ("200_A2", 2, "zero", 0.33077631011278064),
]
HIDDEN_ZERO_WINDOWS = [("100_A2", 0, 0.2914252169900899), ("110_A2", 0, 0.3066350722070048), ("200_A2", 2, 0.33077631011278064)]


def _window(irrep, level, widen_rows=WIDEN_ROWS):
    rows = [r for r in accepted_windows(irrep) if int(r["lattice_level_index"]) == level]
    assert len(rows) == 1, (irrep, level)
    r = rows[0]
    left, right = float(r["E_left_bracket"]), float(r["E_right_bracket"])
    return left - widen_rows * GRID_DE, right + widen_rows * GRID_DE, r


def _evaluate():
    from qc3 import fast
    return fast.evaluate_point


def test_classify_sign_change_on_step9_cases():
    from qc3 import roots
    ev = _evaluate()
    for irrep, level, kind, E in STEP9_POINTS:
        sec = sector(irrep, k2_pwave_form="q0_consistent")
        info = roots.classify_sign_change(ev, sec, C_BASELINE, E)
        tag = (irrep, level, kind, E)
        assert info.kind == kind, (tag, info)
        assert abs(info.discriminant) >= DISC_MIN, (tag, info.discriminant)
        assert (info.discriminant < 0) == (kind == "zero"), tag
        assert info.N > 0 and np.isfinite(info.lnabsdet), tag
        # the two side values agree (the local picture around a simple zero/pole is symmetric)
        assert abs(info.lnabsdet_minus - info.lnabsdet_plus) <= 0.5, (tag, info.lnabsdet_minus, info.lnabsdet_plus)
    # a zero of ln|E - E0| type is a local minimum, a pole a local maximum: kinds are exclusive
    assert {k for (_, _, k, _) in STEP9_POINTS} == {"zero", "pole"}


def test_refine_root_rejects_pole_110_A2_n1_q0():
    from qc3 import roots
    ev = _evaluate()
    sec = sector("110_A2", k2_pwave_form="q0_consistent")
    lo, hi, _ = _window("110_A2", 1)
    info_any = roots.refine_root_info(ev, sec, C_BASELINE, lo, hi, kind="any")
    assert info_any.kind == "pole", info_any
    assert abs(info_any.E - 0.31612911112995534) <= TOL_E, info_any.E
    assert info_any.discriminant >= DISC_MIN, info_any.discriminant
    with pytest.raises(roots.PoleError):
        roots.refine_root(ev, sec, C_BASELINE, lo, hi)          # default kind="zero"
    with pytest.raises(roots.MissingRootError):                 # PoleError is a MissingRootError
        roots.refine_root(ev, sec, C_BASELINE, lo, hi, kind="zero")


def test_refine_root_classifies_12_coded_v2_zeros():
    """Precondition of the step-8 chi^2 proxy: every labelled window inside the cutoff holds a ZERO for coded_v2."""
    from qc3 import roots
    ev = _evaluate()
    n = 0
    for irrep in ("000_A1m", "100_A2", "110_A2", "111_A2", "200_A2"):
        for r in accepted_windows(irrep):
            if r["inside_Ecm_cutoff"].strip().lower() != "true":
                continue
            lo, hi = float(r["E_left_bracket"]), float(r["E_right_bracket"])
            sec = sector(irrep, k2_pwave_form="coded_v2")
            info = roots.refine_root_info(ev, sec, C_BASELINE, lo, hi)
            tag = (irrep, r["lattice_level_index"])
            assert info.kind == "zero", (tag, info)
            assert info.discriminant <= -DISC_MIN, (tag, info.discriminant)
            assert lo < info.E < hi, (tag, info.E, lo, hi)
            assert abs(info.E - float(r["previous_model_zero"])) <= TOL_LABEL, (tag, info.E, r["previous_model_zero"])
            E_plain = roots.refine_root(ev, sec, C_BASELINE, lo, hi)
            assert abs(E_plain - info.E) <= TOL_E, tag
            n += 1
    assert n == 12, n


@pytest.mark.slow
def test_find_zeros_recovers_hidden_zeros():
    from qc3 import roots
    ev = _evaluate()
    for irrep, level, E_zero in HIDDEN_ZERO_WINDOWS:
        sec = sector(irrep, k2_pwave_form="q0_consistent")
        lo, hi, _ = _window(irrep, level)
        changes = roots.find_sign_changes(ev, sec, C_BASELINE, lo, hi, n_scan=241)
        with pytest.raises(roots.MissingRootError):           # equal signs at both ends (pole + zero inside)
            roots.refine_root(ev, sec, C_BASELINE, lo, hi, kind="any")
        kinds = [c.kind for c in changes]
        assert "ambiguous" not in kinds, (irrep, level, changes)
        assert kinds.count("pole") >= 1, (irrep, level, kinds)
        zeros = roots.find_zeros(ev, sec, C_BASELINE, lo, hi, n_scan=241)
        assert len(zeros) == 1, (irrep, level, [z.E for z in zeros])
        assert zeros[0].kind == "zero" and abs(zeros[0].E - E_zero) <= TOL_E, (irrep, level, zeros[0].E, E_zero)
        assert zeros[0].discriminant <= -DISC_MIN
        assert all(lo <= c.E <= hi for c in changes)
        assert [c.E for c in changes] == sorted(c.E for c in changes)

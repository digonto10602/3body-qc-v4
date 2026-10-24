"""M2 steps 2-3 (prompts/003 sections 2.2-2.4): the label-free root finder `qc3.rootscan.scan_sector` and its gate
(CLAUDE.md invariant 7): FP = FN = 0 against every human-labelled bracket, scored with coded_v2 at the label-generation
vector C_BASELINE, in all 7 labelled sectors (5 x L20, 2 x L24; 208 brackets, 31 true).

Scoring (written here independently of qc3.rootscan.score_labels):
  * a labelled TRUE bracket must contain exactly one candidate of kind "zero" (widened by SLACK_ROWS on each side);
  * a labelled FALSE bracket must contain no "zero";
  * every "zero" must lie inside some labelled bracket (no unlabelled zeros);
  * every "free_level" candidate must be within FREE_TOL of a brute-force free level (D8), and no "even_zero" or
    "ambiguous" kinds may remain (all even-order zeros found at the baseline are free-level solutions).
SLACK_ROWS = 0.25 (1.2e-6): the prototype found labelled zeros as close as 1.5e-8 to a bracket edge (200_A2 bracket 18)
while the fast path and the legacy linear interpolation differ by <= 4.4e-8 (reports/tests/011 3d); neighbouring
brackets are >= 7 rows apart, so a quarter row is unambiguous.  DISC_MIN = 8 as in tests/acceptance/m1/test_m1_roots.py
(exact bound 11.5, prototype minimum 11.60).  Evaluation cap 10000 = half the legacy 20000-row grid per sector.
"""
import numpy as np
import pytest

from m2_helpers import (C_BASELINE, E_MAX, FREE_LEVEL_110_A2, GRID_DE, GRID_E0, IRREPS, LABEL_COUNTS, LABEL_EXCLUDE,
                        brute_force_free_levels, read_label_rows, sector)

SLACK_ROWS = 0.25
DISC_MIN = 8.0
FREE_TOL = 1e-7
TOL_E = 1e-9              # root positions, same backend (prompts/001 section 4: Brent xtol 1e-11, measured <= 2e-12)
MAX_EVALUATIONS = 10000
ZERO_LIKE = ("zero", "even_zero", "free_level", "ambiguous")

_SCANS = {}


def _scan(irrep, Lbyas=20, form="coded_v2", E_min=GRID_E0, E_max=E_MAX, c=C_BASELINE):
    key = (irrep, Lbyas, form, E_min, E_max, tuple(c))
    if key not in _SCANS:
        from qc3 import fast, rootscan
        _SCANS[key] = rootscan.scan_sector(fast.evaluate_point, sector(irrep, Lbyas=Lbyas, k2_pwave_form=form), tuple(c),
                                           E_min, E_max)
    return _SCANS[key]


def _score(result, rows, exclude=()):
    zeros = [cnd for cnd in result.candidates if cnd.kind == "zero"]
    tp = fp = fn = tn = 0
    mismatches, covered = [], set()
    for r in rows:
        if r["bracket_id"] in exclude:
            continue
        lo, hi = r["E_left"] - SLACK_ROWS * GRID_DE, r["E_right"] + SLACK_ROWS * GRID_DE
        inside = [z for z in zeros if lo <= z.E <= hi]
        covered.update(id(z) for z in inside)
        if r["label"] == "true":
            if len(inside) == 1:
                tp += 1
            else:
                fn += 1; mismatches.append(("FN", r["bracket_id"], [z.E for z in inside]))
        elif r["label"] == "false":
            if len(inside) == 0:
                tn += 1
            else:
                fp += 1; mismatches.append(("FP", r["bracket_id"], [z.E for z in inside]))
        else:
            raise AssertionError(f"unlabelled bracket {r['bracket_id']}")
    unlabelled = [z.E for z in zeros if id(z) not in covered]
    return dict(tp=tp, fp=fp, fn=fn, tn=tn, mismatches=mismatches, unlabelled_zeros=unlabelled)


def _common_checks(result, irrep, Lbyas):
    from qc3 import fast
    sec = result.sector
    assert sec.k2_pwave_form == "coded_v2" and tuple(result.c) == tuple(C_BASELINE)
    # segments: the scan splits the range at the basis-dimension jumps of the fast path, exactly
    jumps = fast.dimension_jumps(sec, result.E_min, result.E_max)
    assert np.array_equal(np.asarray(result.jumps), jumps)
    edges = [result.E_min] + list(jumps) + [result.E_max]
    assert [tuple(s) for s in result.segments] == [(a, b) for a, b in zip(edges[:-1], edges[1:])]
    assert result.n_evaluations <= MAX_EVALUATIONS, result.n_evaluations
    assert result.n_grid < result.n_evaluations
    Es = [cnd.E for cnd in result.candidates]
    assert Es == sorted(Es)
    E_free, _ = brute_force_free_levels(irrep, Lbyas, result.E_min - 1e-3, result.E_max + 1e-3)
    for cnd in result.candidates:
        assert cnd.kind in ("zero", "pole", "free_level", "even_zero", "ambiguous"), cnd
        assert cnd.segment[0] <= cnd.E <= cnd.segment[1] and cnd.N > 0, cnd
        assert np.isfinite(cnd.discriminant), cnd
        if cnd.kind == "zero":
            assert cnd.discriminant <= -DISC_MIN, cnd
            assert cnd.sign_change is True and cnd.order_int == 1, cnd
            assert cnd.slope is not None and cnd.slope > 0.0, ("physical zeros have d lambda/dE > 0 (arXiv:1901.07095 Eq. 64)", cnd)
            assert np.min(np.abs(E_free - cnd.E)) > FREE_TOL, cnd
        elif cnd.kind == "pole":
            assert cnd.discriminant >= DISC_MIN, cnd
        elif cnd.kind == "free_level":
            assert cnd.free_level is not None and abs(cnd.E - cnd.free_level) <= FREE_TOL, cnd
            assert np.min(np.abs(E_free - cnd.E)) <= FREE_TOL, cnd
            assert cnd.discriminant <= -DISC_MIN, cnd
    kinds = {cnd.kind for cnd in result.candidates}
    assert "ambiguous" not in kinds and "even_zero" not in kinds, [c for c in result.candidates if c.kind in ("ambiguous", "even_zero")]


def _gate(irrep, Lbyas):
    from qc3 import rootscan
    res = _scan(irrep, Lbyas)
    _common_checks(res, irrep, Lbyas)
    rows = read_label_rows(Lbyas, irrep)
    n_br, n_true = LABEL_COUNTS[(Lbyas, irrep)]
    assert len(rows) == n_br and sum(r["label"] == "true" for r in rows) == n_true
    exclude = LABEL_EXCLUDE.get((Lbyas, irrep), ())
    sc = _score(res, rows, exclude)
    assert sc["fp"] == 0 and sc["fn"] == 0, (irrep, Lbyas, sc)
    assert sc["tp"] == n_true - sum(1 for r in rows if r["bracket_id"] in exclude and r["label"] == "true")
    assert sc["unlabelled_zeros"] == [], (irrep, Lbyas, sc)
    # the production scorer must agree with this independent scoring
    prod = rootscan.score_labels(res, rootscan.read_labels(rootscan.LABEL_FILES[(Lbyas, irrep)]), slack_rows=SLACK_ROWS,
                                 exclude_ids=exclude)
    assert (prod["tp"], prod["fp"], prod["fn"], prod["tn"]) == (sc["tp"], sc["fp"], sc["fn"], sc["tn"])
    assert prod["unlabelled_zeros"] == []
    return res, sc


def test_scan_000_A1m_coded_v2_reproduces_labels():
    res, sc = _gate("000_A1m", 20)
    assert (sc["tp"], sc["tn"]) == (3, 7)
    # the two free-level solutions of this sector (orders 2 and 4 at 0.33534 and 0.34453, prototype) are flagged, not zeros
    fls = [c for c in res.candidates if c.kind == "free_level"]
    assert len(fls) == 2 and sorted(c.order_int for c in fls) == [2, 4]
    assert all(c.sign_change is False for c in fls)
    assert len(res.zeros()) == 3 and len(res.poles()) == 7
    assert res.seconds > 0 and isinstance(res.options, dict) and res.options["step_rows"] == 4.0


@pytest.mark.slow
@pytest.mark.parametrize("irrep", [ir for ir in IRREPS if ir != "000_A1m"])
def test_scan_L20_sector_FP_FN_zero(irrep):
    res, sc = _gate(irrep, 20)
    n_br, n_true = LABEL_COUNTS[(20, irrep)]
    assert sc["tp"] == n_true and sc["tn"] == n_br - n_true


@pytest.mark.slow
@pytest.mark.parametrize("irrep", ["000_A1m", "100_A2"])
def test_scan_L24_sector_FP_FN_zero(irrep):
    res, sc = _gate(irrep, 24)
    n_br, n_true = LABEL_COUNTS[(24, irrep)]
    excl = LABEL_EXCLUDE.get((24, irrep), ())
    assert sc["tp"] + sc["fn"] == n_true and sc["tp"] + sc["tn"] == n_br - len(excl)


def test_free_level_solutions_detected_with_order_110_A2_and_200_A2():
    # 110_A2: quartic even-order zero exactly at the free level (reports/tests/007 finding 2) next to the labelled zero 0.3191549
    res = _scan("110_A2", E_min=0.3190, E_max=0.3200)
    fls = [c for c in res.candidates if c.kind == "free_level"]
    assert len(fls) == 1
    fl = fls[0]
    assert abs(fl.E - FREE_LEVEL_110_A2) <= 5e-10 and abs(fl.E - fl.free_level) <= FREE_TOL
    assert fl.sign_change is False and fl.order_int == 4 and abs(fl.order - 4.0) <= 0.05
    assert abs(fl.slope) < 1e-6, ("an even-order zero has d lambda/dE ~ 0 at the root", fl.slope)
    zs = res.zeros()
    assert len(zs) == 1 and abs(zs[0] - 0.3191549427) <= 1e-9        # labelled bracket 17 (coded_v2 root, reports/tests/007)
    # 200_A2: a free-level solution 3.8e-6 below a pole (prototype); it must be flagged free_level, never counted as a zero
    E_free, _ = brute_force_free_levels("200_A2", 20, 0.3465, 0.3470)
    assert len(E_free) == 1
    res2 = _scan("200_A2", E_min=0.3465, E_max=0.3470)
    fls2 = [c for c in res2.candidates if c.kind == "free_level"]
    assert len(fls2) == 1 and abs(fls2[0].E - E_free[0]) <= FREE_TOL
    assert len(res2.zeros()) == 0
    assert any(c.kind == "pole" and 0.34678 < c.E < 0.34679 for c in res2.candidates)


# step-9 q0_consistent roots at C_BASELINE, oracle-verified to <= 2e-12 (reports/tests/007 check 2b; results/k2_pwave_shift_L20.json)
Q0_ROOTS = [("000_A1m", 0.2693226784857993), ("100_A2", 0.2914252169900899), ("100_A2", 0.3003162581565104),
            ("110_A2", 0.30663507220700487), ("110_A2", 0.31564501983275506), ("110_A2", 0.3189395550629675),
            ("110_A2", 0.3286459066230467), ("111_A2", 0.3192435345105379), ("111_A2", 0.33403517451625425),
            ("200_A2", 0.28605537552657084), ("200_A2", 0.2991156800151639), ("200_A2", 0.3307763101127805)]


@pytest.mark.slow
def test_scan_q0_consistent_finds_the_oracle_verified_roots_without_labels():
    """With the production form the labelled windows no longer bracket the roots (5 of 12, reports/tests/007); the
    label-free scan must find all 12 step-9 roots and exactly 12 zeros below the 0.335 cutoff (1/2/4/2/3 per sector)."""
    counts = {"000_A1m": 1, "100_A2": 2, "110_A2": 4, "111_A2": 2, "200_A2": 3}
    for irrep in IRREPS:
        res = _scan(irrep, form="q0_consistent")
        zs = res.zeros()
        assert sum(zs <= 0.335) == counts[irrep], (irrep, zs)
        for ir, E in Q0_ROOTS:
            if ir == irrep:
                assert np.min(np.abs(zs - E)) <= TOL_E, (irrep, E, zs)
        assert not any(c.kind in ("ambiguous", "even_zero") for c in res.candidates)
        for c in res.candidates:
            if c.kind == "zero":
                assert c.slope > 0.0

"""Unit tests for the audit-002 follow-ups (reports/audits/002-m2-milestone-audit.md N1-N5): nothing is dropped silently.
N2 rootscan relabelling + dead bands, N3 full-rescan changes, N4 invalid fits, N1 caveats + opt-in multi-start, N5 Jacobian near a
pole on the real core.  Synthetic inputs are fixed; the multi-start seed is fit.MULTISTART_SEED."""
import json
import pathlib
import warnings
from types import SimpleNamespace

import numpy as np
import pytest

from qc3 import fit, free_levels, rootscan, tracking
from qc3.constants import GRID_DE, GRID_E0
from qc3.sector import Sector

from test_m2_step6_fit import C0 as FIT_C0, setup as fit_setup
from test_rootscan import Fake, SEC, scan
from test_tracking import C0 as TR_C0, FakeQ, SEC as TR_SEC

ROOT = pathlib.Path(__file__).resolve().parents[2]


# ----------------------------------------------------------------------------------------------------------------- N2a
def test_sign_change_zero_at_a_free_level_is_recorded_and_warned():
    Ef = float(free_levels.free_levels(SEC, 0.27, 0.36)[0])
    f = Fake(zeros=[(Ef + 3e-9, 1), (Ef - 2.5e-3, 1)], s0=1.0)
    with pytest.warns(RuntimeWarning, match="relabelled free_level"):
        r = scan(f, 0.27, 0.36)
    assert len(r.relabelled_free_levels) == 1 and len(r.relabelled_sign_change_zeros()) == 1
    rec = r.relabelled_free_levels[0]
    assert rec["kind_orig"] == "zero" and rec["sign_change"] and rec["free_level"] == Ef and rec["E"] == Ef
    assert abs(rec["E_orig"] - (Ef + 3e-9)) < 1e-9 and rec["slope"] is not None
    assert len(r.zeros()) == 1 and abs(r.zeros()[0] - (Ef - 2.5e-3)) < 1e-10          # D8 behaviour itself is unchanged
    assert len(r.of_kind("free_level")) == 1


def test_even_zero_at_a_free_level_is_recorded_without_a_warning():
    Ef = float(free_levels.free_levels(SEC, 0.27, 0.36)[0])
    f = Fake(zeros=[(Ef + 3e-9, 2)], s0=1.0)
    with warnings.catch_warnings():
        warnings.simplefilter("error", RuntimeWarning)
        r = scan(f, 0.27, 0.36)
    assert len(r.relabelled_free_levels) == 1 and r.relabelled_free_levels[0]["kind_orig"] == "even_zero"
    assert not r.relabelled_free_levels[0]["sign_change"] and r.relabelled_sign_change_zeros() == []


# ----------------------------------------------------------------------------------------------------------------- N2b
def test_dead_bands_are_recorded_and_quiet_when_nothing_is_in_them():
    f = Fake(zeros=[(0.2750001, 1), (0.2950001, 1)], jumps=[0.2850001], s0=1.0)
    with warnings.catch_warnings():
        warnings.simplefilter("error", RuntimeWarning)
        r = scan(f)
    assert len(r.dead_bands) == 4 and all(b["probed"] and b["sign_change"] is False and b["error"] is None for b in r.dead_bands)
    assert [b["kind"] for b in r.dead_bands] == ["low_end", "high_end", "low_end", "high_end"]
    b = r.dead_bands[2]
    assert b["lo"] == 0.2850001 and abs(b["hi"] - 0.2850021) < 1e-12 and b["segment"] == (0.2850001, 0.30)
    assert r.dead_band_sign_changes() == [] and len(r.zeros()) == 2


@pytest.mark.parametrize("where,lo_expected", [("above_jump", 0.2850001), ("below_jump", 0.2850001 - 2e-6)])
def test_zero_inside_a_dead_band_is_detected_by_the_probe_and_warned(where, lo_expected):
    z = 0.2850001 + 1.0e-6 if where == "above_jump" else 0.2850001 - 1.0e-6
    f = Fake(zeros=[(z, 1)], jumps=[0.2850001], s0=1.0)
    with pytest.warns(RuntimeWarning, match="changes sign across the unscanned band"):
        r = scan(f)
    sc = r.dead_band_sign_changes()
    assert len(sc) == 1 and abs(sc[0]["lo"] - lo_expected) < 1e-12 and sc[0]["lo"] <= z <= sc[0]["hi"]
    assert len(r.zeros()) == 0            # the dead band zero is NOT reported as a root, only flagged (never silent)


def test_segment_shorter_than_two_edges_is_recorded_and_warned():
    f = Fake(zeros=[(0.2750001, 1)], jumps=[0.2850001, 0.2850031], s0=1.0)
    with pytest.warns(RuntimeWarning, match="not scanned at all"):
        r = scan(f)
    whole = [b for b in r.dead_bands if b["kind"] == "whole_segment"]
    assert len(whole) == 1 and whole[0]["probed"] and whole[0]["segment"] == (0.2850001, 0.2850031)


def test_free_level_next_to_a_segment_end_is_listed_as_skipped():
    Ef = float(free_levels.free_levels(SEC, 0.27, 0.36)[0])
    f = Fake(jumps=[Ef + 1e-6], s0=1.0)
    r = scan(f, 0.27, 0.36)
    assert any(abs(s["E"] - Ef) < 1e-15 and "not probed" in s["reason"] for s in r.skipped_free_levels)


# ----------------------------------------------------------------------------------------------------------------- N3
def test_full_rescan_change_is_counted_listed_and_warned():
    f = FakeQ()
    sr = tracking.SectorRoots(TR_SEC, TR_C0, 0.270, 0.300, evaluate=f)
    assert (sr.n_full_rescan_changes, sr.full_rescan_events) == (0, [])
    dc = np.array([0.02, 0.0, 0.0, 0.0])
    cs = tuple(float(x) for x in np.add(TR_C0, dc))
    # the scan at the displaced c, as _load_full_scan would have produced it
    sr._load_full_scan(cs)
    new = [r.E for r in sr.roots]
    with warnings.catch_warnings():
        warnings.simplefilter("error", RuntimeWarning)
        sr._note_full_rescan(cs, [(E - 1e-9, E) for E in new])       # prediction == new root: no event
    assert sr.n_full_rescan_changes == 0
    far = [(E - 1e-3, E - 1e-3 + 1e-8) for E in new]                          # prediction 1e-3 away from the root found, tolerance ~1e-6
    with pytest.warns(RuntimeWarning, match="changed the tracked set"):
        sr._note_full_rescan(cs, far)
    assert sr.n_full_rescan_changes == 1 and sr.full_rescan_events[0]["max_deviation"] == pytest.approx(1e-3, abs=1e-7)
    assert not sr.full_rescan_events[0]["count_changed"]
    with pytest.warns(RuntimeWarning, match="3 -> 2|2 -> 3|changed the tracked set"):
        sr._note_full_rescan(cs, [(E - 1e-9, E) for E in new[:2]])   # one fewer prediction than roots found
    assert sr.n_full_rescan_changes == 2 and sr.full_rescan_events[1]["count_changed"]


def test_real_full_rescan_through_update_records_deviation_or_nothing():
    # safeguards that never accept + a window too narrow for the local rescan -> full rescan inside update(); the counter and the
    # event list must be consistent whatever the outcome (an event is a warning, never an exception)
    f = FakeQ()
    sr = tracking.SectorRoots(TR_SEC, TR_C0, 0.270, 0.300, evaluate=f,
                              track_options=dict(safeguard_factor=0.0, safeguard_abs=0.0, local_rows=0))
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        sr.update(np.add(TR_C0, [1e-12, 0.0, 0.0, 0.0]))
    assert sr.n_full_rescans >= 1 and sr.n_full_rescan_changes == len(sr.full_rescan_events) <= sr.n_full_rescans


# ----------------------------------------------------------------------------------------------------------------- N4
def _bad_verification():
    return {"ok": False, "counts_equal": False, "max_root_diff": float("inf"), "chi2_fresh": 5.0, "chi2_diff": 1.5, "chi2_bound": 1e-6,
            "per_sector": {}, "wrong_sign_zeros_fresh": [], "n_evaluations": 0, "seconds": 0.0, "margin": 0.002}


def test_failed_verification_marks_the_fit_invalid_and_warns(monkeypatch, capsys):
    lv, tr = fit_setup()
    monkeypatch.setattr(fit, "verify_fit", lambda *a, **k: _bad_verification())
    r = fit.fit_cutoff(lv, tr, FIT_C0 + np.array([2e3, -1e3, 5e2, 3e3]))
    assert r.valid is False and any(x.startswith("verification_failed") for x in r.invalid_reasons)
    assert "WARNING" in capsys.readouterr().err
    text = fit.format_fit(r)
    assert "INVALID FIT" in text and "valid=False" in text
    assert fit._jsonable({"v": r.valid, "r": r.invalid_reasons})["v"] is False


def test_passing_verification_keeps_the_fit_valid(monkeypatch):
    lv, tr = fit_setup()
    ok = dict(_bad_verification(), ok=True, counts_equal=True, max_root_diff=0.0, chi2_diff=0.0)
    monkeypatch.setattr(fit, "verify_fit", lambda *a, **k: ok)
    r = fit.fit_cutoff(lv, tr, FIT_C0 + np.array([2e3, -1e3, 5e2, 3e3]))
    assert r.valid and r.invalid_reasons == [] and "INVALID" not in fit.format_fit(r)


# ----------------------------------------------------------------------------------------------------------------- N1
START = FIT_C0 + np.array([2e3, -1e3, 5e2, 3e3])


def test_caveats_default_list_and_json_field():
    lv, tr = fit_setup()
    r = fit.fit_cutoff(lv, tr, START, verify=False)
    txt = " | ".join(r.caveats)
    assert "NOT rescaled by chi2/ndof" in txt and "no multi-start" in txt and "not verified" in txt
    assert r.multistart == 0 and r.multistart_results == [] and r.multistart_seed is None
    assert "CAVEAT:" in fit.format_fit(r)
    d = {f.name for f in fit.fields(fit.FitResult)}
    assert {"caveats", "invalid_reasons", "multistart", "multistart_results", "n_full_rescan_changes"} <= d


def test_multistart_is_opt_in_seeded_reproducible_and_does_not_change_the_fit():
    lv, tr = fit_setup()
    r0 = fit.fit_cutoff(lv, tr, START, verify=False)
    lv, tr = fit_setup()
    r1 = fit.fit_cutoff(lv, tr, START, verify=False, multistart=3)
    lv, tr = fit_setup()
    r2 = fit.fit_cutoff(lv, tr, START, verify=False, multistart=3)
    assert r1.multistart == 3 and len(r1.multistart_results) == 3 and r1.multistart_seed == fit.MULTISTART_SEED == 20261027
    assert np.array_equal(r1.c, r0.c) and r1.chi2 == r0.chi2                              # the reported fit is untouched
    assert [e["c_start"] for e in r1.multistart_results] == [e["c_start"] for e in r2.multistart_results]   # seeded
    for e in r1.multistart_results:                                                      # exactly quadratic chi2: same minimum
        assert e["valid"] and abs(e["delta_chi2"]) < fit.MULTISTART_CHI2_TOL
        rel = np.abs(np.array(e["c_start"]) - r1.c) / np.abs(r1.c)
        assert 0 < rel.max() < 1e-4
    assert any("multi-start check passed" in c for c in r1.caveats) and not any("no multi-start" in c for c in r1.caveats)


def test_multistart_finding_a_different_minimum_adds_the_not_robust_caveat(monkeypatch):
    lv, tr = fit_setup()
    fake = [{"index": 0, "valid": True, "c": [0.0] * 4, "chi2": 99.0, "delta_chi2": 50.0, "edm": 1e-5, "c_start": [0.0] * 4},
            {"index": 1, "valid": False, "error": "AmbiguousRootError: boom"}]
    monkeypatch.setattr(fit, "_multistart", lambda *a, **k: fake)
    r = fit.fit_cutoff(lv, tr, START, verify=False, multistart=2)
    txt = " | ".join(r.caveats)
    assert "minimum not robust: 1 of 2" in txt and "AmbiguousRootError: boom" in txt


def test_errors_disagree_and_piecewise_smooth_caveats():
    lv, tr = fit_setup()
    r = fit.fit_cutoff(lv, tr, START, verify=False)
    fcn = fit.QCChi2(lv, tr)
    # a pole within 1 sigma of the first model level and a dimension jump within 1 sigma of the second
    tr["000_A1m"]._poles = [float(r.E_model[0]) + 0.5e-3]
    tr["000_A1m"].scan = SimpleNamespace(jumps=np.array([float(r.E_model[1]) - 0.5e-3]))
    tr["100_A2"]._poles = [float(r.E_model[3]) + 5e-3]                                  # farther than 1 sigma (1e-3): not reported
    tr["100_A2"].scan = SimpleNamespace(jumps=np.zeros(0))
    items = fit._piecewise_smooth_items(fcn, tr, r.E_model)
    assert [i["level"] for i in items] == [["000_A1m", 0], ["000_A1m", 1]] and items[0]["poles"] and items[1]["jumps"]
    r.errors_disagree, r.errors_disagree_params = True, ["Kiso0"]
    txt = " | ".join(fit._make_caveats(r, fcn, tr, [0, 1, 2, 3]))
    assert "Hesse errors unreliable" in txt and "only piecewise smooth" in txt


def test_near_threshold_error_gap_is_reported():
    lv, tr = fit_setup()
    r = fit.fit_cutoff(lv, tr, START, verify=False)
    fcn = fit.QCChi2(lv, tr)
    r.error_gap = np.array([1.15, 1.0, 0.95, np.nan])
    r.errors_disagree = False
    txt = " | ".join(fit._make_caveats(r, fcn, tr, [0, 1, 2]))
    assert "marginal" in txt and "Kiso0" in txt and "Kiso1" not in txt.split("marginal")[0].split("for")[-1]


def test_cli_has_multistart_option():
    a = fit.build_parser().parse_args(["--multistart", "3"])
    assert a.multistart == 3 and fit.build_parser().parse_args([]).multistart == 0


# ----------------------------------------------------------------------------------------------------------------- N5
def test_root_jacobian_vs_central_finite_differences_110_A2_near_a_pole(capsys):
    """Real core, 110_A2, q0_consistent, c = the A/0.335 minimum of results/fit_L20_q0_consistent_m2.json.  The analytic eq. (1)
    Jacobian (Q' by a two-point central difference, h = 1e-6) against central differences of the tracked roots (tracker.update to
    c +- d_i, d_i set so that the root moves by ~3e-7).  The relative error is printed and bounded."""
    from qc3 import fast
    doc = json.loads((ROOT / "results" / "fit_L20_q0_consistent_m2.json").read_text())
    c = np.asarray([f for f in doc["fits"] if abs(f["cutoff"] - 0.335) < 1e-9][0]["c"], dtype=float)
    sec = Sector.production("110_A2", Lbyas=20, k2_pwave_form="q0_consistent")
    sr = tracking.SectorRoots(sec, c, GRID_E0, 0.340, evaluate=fast.evaluate_point)
    E0, J = sr.energies(), sr.jacobian()
    poles = np.asarray(sr._poles)
    d_pole = np.array([np.min(np.abs(poles - e)) if len(poles) else np.inf for e in E0])
    Jfd = np.zeros_like(J)
    for i in range(4):
        d = np.zeros(4)
        scale = np.max(np.abs(J[:, i]))
        d[i] = 3e-7 / scale
        Ep = sr.update(c + d).copy()
        Em = sr.update(c - d).copy()
        sr.update(c)
        assert len(Ep) == len(Em) == len(E0)
        Jfd[:, i] = (Ep - Em) / (2.0 * d[i])
    big = np.abs(J) > 1e-3 * np.max(np.abs(J), axis=1, keepdims=True)           # components that matter for the root's movement
    rel = np.abs(J - Jfd) / np.abs(Jfd)
    with capsys.disabled():
        print("\n[N5] 110_A2 A/0.335: root E, distance to nearest pole [rows], max rel. error (components > 1e-3 of the largest)")
        for k in range(len(E0)):
            print(f"[N5]   E={E0[k]:.9f}  d_pole={d_pole[k] / GRID_DE:7.2f} rows  max rel err = {np.max(np.where(big[k], rel[k], 0.0)):.2e}"
                  f"  per-component {np.array2string(rel[k], precision=1)}")
    assert np.all(rel[big] < 1e-3)         # observed <= 2.4e-5 (the audit estimated ~2e-3 for the nearest-pole zero)

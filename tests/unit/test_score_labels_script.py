"""M2 step 3: scoring, JSON record and exit code of scripts/m2/score_labels.py (fake scans; no QC evaluation)."""
import importlib.util
import json
import pathlib

import numpy as np
import pytest

from qc3 import rootscan
from qc3.constants import GRID_DE

ROOT = pathlib.Path(__file__).resolve().parents[2]
_spec = importlib.util.spec_from_file_location("score_labels_script", ROOT / "scripts" / "m2" / "score_labels.py")
SL = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(SL)

KEY = (20, "000_A1m")
LABELS = rootscan.read_labels(rootscan.LABEL_FILES[KEY])          # 10 brackets, 3 true


def cand(E, kind="zero", free_level=None):
    return rootscan.RootCandidate(E=E, kind=kind, order=1.0, order_int=1, sign_change=kind != "free_level", discriminant=-20.0,
                                  slope=1.0, N=10, segment=(0.0, 1.0), bracket=(E, E), via="sign", free_level=free_level,
                                  n_evaluations=5)


def fake_scan(mutate=None):
    """A scan that reproduces the labels exactly (a zero in the middle of each TRUE bracket), then applies mutate(cands)."""
    def scan(sec, c, a, b):
        cands = [cand(0.5 * (r["E_left"] + r["E_right"])) for r in LABELS if r["label"] == "true"]
        if mutate:
            cands = mutate(cands)
        return rootscan.ScanResult(sector=sec, c=tuple(c), E_min=a, E_max=b, candidates=sorted(cands, key=lambda t: t.E),
                                   jumps=np.array([0.3]), segments=[(a, 0.3), (0.3, b)], n_evaluations=123, n_grid=100,
                                   seconds=0.5, options={}, grid=np.zeros(100))
    return scan


def run(tmp_path, scan, only=("20:000_A1m",)):
    out = tmp_path / "gate.json"
    code = SL.main(["--out", str(out), "--only", *only], scan=scan, meta={"git_commit": "test"})
    return code, json.loads(out.read_text())


def test_perfect_scan_passes_with_exit_code_zero_and_writes_the_record(tmp_path):
    code, doc = run(tmp_path, fake_scan())
    assert code == 0 and doc["gate_passed"] is True and doc["git_commit"] == "test"
    s = doc["sectors"]["L20/000_A1m"]
    assert (s["tp"], s["fp"], s["fn"], s["tn"]) == (3, 0, 0, 7) and s["n_scored"] == 10 and s["passed"]
    assert s["n_evaluations"] == 123 and s["seconds"] == 0.5 and len(s["candidates"]) == 3 and s["kinds"] == {"zero": 3}
    assert {"E", "kind", "order_int", "discriminant", "slope", "via", "free_level"} <= set(s["candidates"][0])


def test_missing_zero_is_a_false_negative_and_exits_nonzero(tmp_path):
    code, doc = run(tmp_path, fake_scan(lambda cs: cs[1:]))
    s = doc["sectors"]["L20/000_A1m"]
    assert code == 1 and doc["gate_passed"] is False and s["fn"] == 1 and s["fp"] == 0 and not s["passed"]
    assert s["mismatches"][0]["type"] == "FN" and any("false negative" in p for p in s["problems"])


def test_zero_in_a_false_bracket_is_a_false_positive_and_exits_nonzero(tmp_path):
    fb = next(r for r in LABELS if r["label"] == "false")
    code, doc = run(tmp_path, fake_scan(lambda cs: cs + [cand(0.5 * (fb["E_left"] + fb["E_right"]))]))
    s = doc["sectors"]["L20/000_A1m"]
    assert code == 1 and s["fp"] == 1 and s["fn"] == 0 and s["mismatches"][0]["bracket_id"] == fb["bracket_id"]


def test_unlabelled_zero_ambiguous_and_stray_free_level_fail_the_gate(tmp_path):
    assert run(tmp_path, fake_scan(lambda cs: cs + [cand(0.2)]))[0] == 1                              # outside every bracket
    assert run(tmp_path, fake_scan(lambda cs: cs + [cand(0.2, "ambiguous")]))[0] == 1
    assert run(tmp_path, fake_scan(lambda cs: cs + [cand(0.2, "even_zero")]))[0] == 1
    assert run(tmp_path, fake_scan(lambda cs: cs + [cand(0.2, "free_level", 0.2)]))[0] == 1            # 0.2 is not a free level
    assert run(tmp_path, fake_scan(lambda cs: cs + [cand(0.2, "pole")]))[0] == 0                       # poles are fine


def test_free_level_candidate_at_a_free_level_is_accepted(tmp_path):
    from qc3.free_levels import free_levels
    from qc3.sector import Sector
    Ef = float(free_levels(Sector.production("000_A1m", Lbyas=20, k2_pwave_form="coded_v2"), 0.2631, 0.36)[0])
    code, doc = run(tmp_path, fake_scan(lambda cs: cs + [cand(Ef, "free_level", Ef)]))
    assert code == 0 and doc["sectors"]["L20/000_A1m"]["kinds"]["free_level"] == 1


def test_exclusion_is_applied_and_reported_for_L24_100_A2():
    labels = rootscan.read_labels(rootscan.LABEL_FILES[(24, "100_A2")])
    assert SL.EXCLUDE[(24, "100_A2")] == (36,) and any(r["bracket_id"] == 36 for r in labels)
    zeros = [cand(0.5 * (r["E_left"] + r["E_right"])) for r in labels if r["label"] == "true" and r["bracket_id"] != 36]
    res = rootscan.ScanResult(sector=None, c=(0,) * 4, E_min=0.0, E_max=1.0, candidates=zeros, jumps=np.zeros(0), segments=[],
                              n_evaluations=1, n_grid=1, seconds=0.0, options={}, grid=np.zeros(1))
    rec = SL.score_sector(res, labels, exclude_ids=SL.EXCLUDE[(24, "100_A2")])
    assert rec["excluded_bracket_ids"] == [36] and rec["n_scored"] == len(labels) - 1 and rec["passed"]


def test_exit_code_rules_and_bad_arguments(tmp_path):
    ok = {"a": {"passed": True, "fp": 0, "fn": 0}}
    assert SL.exit_code(ok) == 0
    assert SL.exit_code({"a": {"passed": True, "fp": 1, "fn": 0}}) == 1
    assert SL.exit_code({"a": {"passed": True, "fp": 0, "fn": 2}}) == 1
    assert SL.exit_code({"a": {"passed": False, "fp": 0, "fn": 0}}) == 1
    assert SL.exit_code({}) == 1                                    # nothing scored is not a pass
    with pytest.raises(ValueError, match="unknown sector"):
        SL.run_gate(["20:999_A9"], scan=fake_scan())


def test_metadata_records_commit_constants_and_form():
    m = SL.metadata()
    assert m["k2_pwave_form"] == "coded_v2" and m["production_constants"]["xi"] == 3.444 and len(m["c_baseline"]) == 4
    assert m["git_commit"] and "core_built_at_commit" in m and m["excluded_brackets"] == {"L24/100_A2": [36]}

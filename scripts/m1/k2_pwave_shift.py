"""M1 step 9: p-wave K2^{-1} form quantification for decision D7 (prompts/001 sec 2.2 and sec 3 row 9).

For the 12 L20 levels at Ecm cutoff 0.335 and the baseline c, refine the zero of det Q with the fast path
for k2_pwave_form = coded_v2 and q0_consistent; report the level shift against the lattice sigma (D1 covariance).

Run (single thread):
  OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 PYTHONPATH=src/python python3 scripts/m1/k2_pwave_shift.py
"""
import csv
import json
import pathlib
import platform
import subprocess
import sys
import time

import numpy as np

from qc3 import constants, fast, lattice
from qc3.point import MissingRootError
from qc3.roots import find_sign_changes, refine_root_info
from qc3.sector import Sector

ROOT = pathlib.Path(__file__).resolve().parents[2]
CUTOFF = 0.335
LBYAS = 20
WIDEN_ROWS = 250
FORMS = ("coded_v2", "q0_consistent")
SCAN_PTS = 241


def read_windows():
    out = {}
    for irrep in constants.IRREPS:
        p = ROOT / "data" / "windows" / f"L{LBYAS}" / f"L{LBYAS}_{irrep}" / "accepted_windows.csv"
        with open(p) as fh:
            for r in csv.DictReader(fh):
                if r["inside_Ecm_cutoff"].strip().lower() == "true":
                    out[(irrep, int(r["lattice_level_index"]))] = r
    return out


def _ev(sec, c, E, counter):
    counter[0] += 1
    return fast.evaluate_point(sec, float(E), c, want=("Q",))


def find_root(sec, c, left, right, counter):
    """Plan rule: refine_root_info (default kind="zero": poles raise PoleError) in [left,right]; if there is no sign change
    or the sign change is a pole/ambiguous, widen by +-WIDEN_ROWS grid rows.  Returns dict(E, window, note, ln_root, ends,
    kind, discriminant)."""
    def ev(sector, E, cc, want=("Q",)):
        counter[0] += 1
        return fast.evaluate_point(sector, E, cc, want=want)
    wl, wr = left - WIDEN_ROWS * constants.GRID_DE, right + WIDEN_ROWS * constants.GRID_DE
    notes = []
    for tag, (lo, hi) in (("accepted window", (left, right)), (f"widened +-{WIDEN_ROWS} rows", (wl, wr))):
        try:
            info = refine_root_info(ev, sec, c, lo, hi)
        except MissingRootError as e:      # includes PoleError / AmbiguousRootError
            notes.append(f"{tag}: {type(e).__name__}: {e}")
            continue
        ends = [float(_ev(sec, c, x, counter).lnabsdet) for x in (lo, hi)]
        return {"E": info.E, "window": (lo, hi), "note": tag + ("; earlier: " + " | ".join(notes) if notes else ""),
                "ln_root": info.lnabsdet, "ends": ends, "kind": info.kind, "discriminant": info.discriminant}
    return {"E": None, "window": (wl, wr), "note": "MissingRootError (plan rule): " + " | ".join(notes), "ln_root": None,
            "ends": None, "kind": None, "discriminant": None}


def scan_roots(sec, c, lo, hi, npts, counter):
    """Supplementary (not the plan rule): qc3.roots.find_sign_changes on npts points; keep the zeros.
    Returns (list of dict(E, ln_root, bracket, discriminant) sorted by E, number of pole sign changes)."""
    def ev(sector, E, cc, want=("Q",)):
        counter[0] += 1
        return fast.evaluate_point(sector, E, cc, want=want)
    changes = find_sign_changes(ev, sec, c, lo, hi, n_scan=npts)
    zs = [{"E": r.E, "ln_root": r.lnabsdet, "bracket": list(r.bracket), "discriminant": r.discriminant}
          for r in changes if r.kind == "zero"]
    return zs, sum(1 for r in changes if r.kind != "zero")


def main():
    t_start = time.perf_counter()
    wins = read_windows()
    ls = lattice.load_levels(ROOT / "data", LBYAS, list(constants.IRREPS), CUTOFF)
    labels = [tuple(l) for l in ls.labels]
    assert len(wins) == 12, f"expected 12 windows inside the cutoff, got {len(wins)}"
    assert len(labels) == 12, f"expected 12 lattice levels, got {len(labels)}"
    assert set(wins) == set(labels), f"window/level mismatch: {sorted(set(wins) ^ set(labels))}"
    sigma = np.sqrt(np.diag(ls.cov))
    c = tuple(constants.C_BASELINE)
    rows, timings = [], {}
    for i, (irrep, n) in enumerate(labels):
        w = wins[(irrep, n)]
        left, right = float(w["E_left_bracket"]), float(w["E_right_bracket"])
        row = {"irrep": irrep, "level_index": n, "lattice_Ecm_window_csv": float(w["lattice_Ecm"]),
               "lattice_Ecm_mean": float(ls.Ecm[i]), "lattice_sigma": float(sigma[i]),
               "accepted_window": [left, right], "previous_model_zero": float(w["previous_model_zero"])}
        for form in FORMS:
            sec = Sector.production(irrep, LBYAS, k2_pwave_form=form)
            cnt = [0]
            t0 = time.perf_counter()
            fr = find_root(sec, c, left, right, cnt)
            dt = time.perf_counter() - t0
            entry = {"E_root": fr["E"], "window_used": list(fr["window"]), "note": fr["note"],
                     "kind": fr["kind"], "discriminant": fr["discriminant"], "lnabsdet_at_root": fr["ln_root"], "lnabsdet_at_window_ends": fr["ends"],
                     "n_evaluations": cnt[0], "seconds": dt}
            if form == "q0_consistent":
                # supplementary scan in the +-WIDEN_ROWS window: all validated zeros there (not part of the plan rule)
                cnt2 = [0]
                t0 = time.perf_counter()
                zs, npol = scan_roots(sec, c, left - WIDEN_ROWS * constants.GRID_DE, right + WIDEN_ROWS * constants.GRID_DE,
                                      SCAN_PTS, cnt2)
                entry["scan"] = {"npts": SCAN_PTS, "zeros_in_widened_window": zs, "n_pole_sign_changes": npol,
                                 "n_evaluations": cnt2[0], "seconds": time.perf_counter() - t0}
                ref = row["coded_v2"]["E_root"]
                entry["E_root_scan_nearest_to_coded_v2"] = (min(zs, key=lambda z: abs(z["E"] - ref))["E"] if zs else None)
            row[form] = entry
            timings[f"{irrep}_{n}_{form}"] = dt
            print(f"{irrep} n={n} {form}: E={fr['E']} evals={cnt[0]} t={dt:.1f}s {fr['note']}", flush=True)
            if form == "q0_consistent":
                print(f"    scan: zeros {[z['E'] for z in zs]} poles {npol}", flush=True)
        a, b = row["coded_v2"]["E_root"], row["q0_consistent"]["E_root"]
        row["dE"] = None if (a is None or b is None) else b - a
        row["dE_over_sigma"] = None if row["dE"] is None else row["dE"] / row["lattice_sigma"]
        bs = row["q0_consistent"]["E_root_scan_nearest_to_coded_v2"]
        row["dE_scan_assisted"] = None if (a is None or bs is None) else bs - a
        row["dE_scan_assisted_over_sigma"] = None if row["dE_scan_assisted"] is None else row["dE_scan_assisted"] / row["lattice_sigma"]
        row["coded_v2_minus_previous_model_zero"] = None if a is None else a - row["previous_model_zero"]
        rows.append(row)
    commit = subprocess.run(["git", "rev-parse", "--short", "HEAD"], cwd=ROOT, capture_output=True, text=True).stdout.strip()
    dirty = bool(subprocess.run(["git", "status", "--porcelain"], cwd=ROOT, capture_output=True, text=True).stdout.strip())
    diffs = [abs(r["coded_v2_minus_previous_model_zero"]) for r in rows if r["coded_v2_minus_previous_model_zero"] is not None]
    out = {
        "description": "K2^{-1} p-wave form level shifts, L20, cutoff 0.335, baseline c (prompts/001 step 9, decision D7)",
        "git_commit": commit, "git_tree_dirty": dirty,
        "command": "OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 PYTHONPATH=src/python python3 scripts/m1/k2_pwave_shift.py",
        "machine": platform.platform(), "python": sys.version.split()[0],
        "production_constants": constants.PRODUCTION, "c_baseline": list(c), "Lbyas": LBYAS, "Ecm_cutoff": CUTOFF,
        "backend": "fast", "covariance": "D1 (jackknife of per-configuration values, qc3.lattice.load_levels)",
        "Ncfg": ls.Ncfg, "widen_rows": WIDEN_ROWS, "grid_dE": constants.GRID_DE,
        "max_abs_coded_v2_minus_previous_model_zero": max(diffs) if diffs else None,
        "levels": rows, "timings_seconds": timings, "total_seconds": time.perf_counter() - t_start,
    }
    (ROOT / "results").mkdir(exist_ok=True)
    (ROOT / "results" / "k2_pwave_shift_L20.json").write_text(json.dumps(out, indent=2))
    print(f"total {out['total_seconds']:.1f}s; max|coded_v2 - previous_model_zero| = {out['max_abs_coded_v2_minus_previous_model_zero']}")


if __name__ == "__main__":
    main()

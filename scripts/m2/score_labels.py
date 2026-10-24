"""M2 step 3 (prompts/003 section 3): the label gate.  Runs the label-free scan `qc3.rootscan.scan_sector` on the 7 labelled
sectors (rootscan.LABEL_FILES: 5 x L20, 2 x L24), coded_v2 at the label-generation vector C_BASELINE, scores every sector
against the human labels with `rootscan.score_labels` (CLAUDE.md invariant 7: FP = FN = 0) and writes
results/m2_label_gate.json (commit, environment, constants, per-sector counts, evaluations, seconds, every candidate).

Gate (exit code 1 if violated, 0 otherwise):
  * FP = FN = 0 in every sector, and no kind-"zero" outside every non-excluded labelled bracket,
  * no "ambiguous" / "even_zero" candidate (every even-order zero must be a D8 free-level solution),
  * every "free_level" candidate within FREE_TOL of a non-interacting level (qc3.free_levels).
L24/100_A2 bracket 36 is excluded and reported (research 004 section 3.5, unresolved label).

Run (single thread, about 2 min):
  OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 PYTHONPATH=src/python python3 scripts/m2/score_labels.py [--out PATH] [--only 20:100_A2 ...]
"""
from __future__ import annotations

import argparse
import json
import pathlib
import platform
import subprocess
import sys
import time

ROOT = pathlib.Path(__file__).resolve().parents[2]
SLACK_ROWS = 0.25            # as tests/acceptance/m2/test_m2_rootscan.py
FREE_TOL = 1e-7
E_MAX = 0.36
FORM = "coded_v2"            # the labels were made with coded_v2 at C_BASELINE
EXCLUDE = {(24, "100_A2"): (36,)}
DEFAULT_OUT = ROOT / "results" / "m2_label_gate.json"


def _jsonable(x):
    import numpy as np
    if isinstance(x, dict):
        return {str(k): _jsonable(v) for k, v in x.items()}
    if isinstance(x, (list, tuple)):
        return [_jsonable(v) for v in x]
    if isinstance(x, np.ndarray):
        return [_jsonable(v) for v in x.tolist()]
    if isinstance(x, (np.floating, float)):
        return float(x) if np.isfinite(x) else repr(float(x))
    if isinstance(x, np.integer):
        return int(x)
    return x


def candidate_record(c) -> dict:
    return {"E": c.E, "kind": c.kind, "order": c.order, "order_int": c.order_int, "sign_change": c.sign_change,
            "discriminant": c.discriminant, "slope": c.slope, "N": c.N, "segment": list(c.segment), "bracket": list(c.bracket),
            "via": c.via, "free_level": c.free_level, "n_evaluations": c.n_evaluations}


def score_sector(result, labels, exclude_ids=(), slack_rows=SLACK_ROWS, free_levels=None, free_tol=FREE_TOL) -> dict:
    """Per-sector record for the JSON: counts (qc3.rootscan.score_labels), the problems that fail the gate, every candidate.
    `free_levels`: array of the sector's non-interacting levels (None: skip the free-level check)."""
    from qc3 import rootscan
    sc = rootscan.score_labels(result, labels, slack_rows=slack_rows, exclude_ids=exclude_ids)
    problems = []
    if sc["fp"]:
        problems.append(f"{sc['fp']} false positive(s)")
    if sc["fn"]:
        problems.append(f"{sc['fn']} false negative(s)")
    if sc["unlabelled_zeros"]:
        problems.append(f"{len(sc['unlabelled_zeros'])} zero(s) outside every scored bracket")
    for kind in ("ambiguous", "even_zero"):
        n = len(result.of_kind(kind))
        if n:
            problems.append(f"{n} candidate(s) of kind {kind!r} (every even-order zero must be a free-level solution)")
    if free_levels is not None:
        import numpy as np
        fl = np.asarray(free_levels, dtype=float)
        for c in result.of_kind("free_level"):
            if fl.size == 0 or float(np.min(np.abs(fl - c.E))) > free_tol:
                problems.append(f"free_level candidate at E={c.E!r} is not within {free_tol} of a free level")
    kinds = {}
    for c in result.candidates:
        kinds[c.kind] = kinds.get(c.kind, 0) + 1
    return {"tp": sc["tp"], "fp": sc["fp"], "fn": sc["fn"], "tn": sc["tn"], "n_scored": sc["tp"] + sc["fp"] + sc["fn"] + sc["tn"],
            "n_labelled": len(labels), "excluded_bracket_ids": sc["excluded"], "mismatches": sc["mismatches"],
            "unlabelled_zeros": [candidate_record(z) for z in sc["unlabelled_zeros"]],
            "problems": problems, "passed": not problems, "kinds": kinds,
            "n_evaluations": result.n_evaluations, "n_grid": result.n_grid, "seconds": result.seconds,
            "jumps": result.jumps, "segments": result.segments, "candidates": [candidate_record(c) for c in result.candidates]}


def gate_failed(sectors: dict) -> bool:
    """True iff any sector has FP or FN or another gate problem (sectors: {name: score_sector record})."""
    return any((not s["passed"]) or s["fp"] > 0 or s["fn"] > 0 for s in sectors.values())


def exit_code(sectors: dict) -> int:
    return 1 if (not sectors or gate_failed(sectors)) else 0


def summary_lines(sectors: dict) -> list:
    out = [f"{'sector':<14}{'TP':>4}{'FP':>4}{'FN':>4}{'TN':>5}{'scored':>8}{'evals':>8}{'sec':>8}  kinds / problems"]
    for name, s in sectors.items():
        kinds = ",".join(f"{k}={v}" for k, v in sorted(s["kinds"].items()))
        out.append(f"{name:<14}{s['tp']:>4}{s['fp']:>4}{s['fn']:>4}{s['tn']:>5}{s['n_scored']:>8}{s['n_evaluations']:>8}"
                   f"{s['seconds']:>8.1f}  {kinds}" + ("" if s["passed"] else "  PROBLEMS: " + "; ".join(s["problems"])))
    out.append(f"{'total':<14}{sum(s['tp'] for s in sectors.values()):>4}{sum(s['fp'] for s in sectors.values()):>4}"
               f"{sum(s['fn'] for s in sectors.values()):>4}{sum(s['tn'] for s in sectors.values()):>5}"
               f"{sum(s['n_scored'] for s in sectors.values()):>8}{sum(s['n_evaluations'] for s in sectors.values()):>8}"
               f"{sum(s['seconds'] for s in sectors.values()):>8.1f}")
    return out


def _git(*args):
    r = subprocess.run(["git", *args], cwd=ROOT, capture_output=True, text=True)
    return r.stdout.strip()


def metadata() -> dict:
    import numpy
    import scipy
    from qc3 import _core, constants
    return {"git_commit": _git("rev-parse", "--short", "HEAD"), "git_tree_dirty": bool(_git("status", "--porcelain")),
            "core_built_at_commit": _core.git_commit(), "machine": platform.platform(), "python": sys.version.split()[0],
            "numpy": numpy.__version__, "scipy": scipy.__version__,
            "production_constants": constants.PRODUCTION, "k2_pwave_form": FORM, "c_baseline": list(constants.C_BASELINE),
            "grid_E0": constants.GRID_E0, "grid_dE": constants.GRID_DE, "E_max": E_MAX, "backend": "fast", "slack_rows": SLACK_ROWS,
            "free_tol": FREE_TOL, "excluded_brackets": {f"L{L}/{ir}": list(v) for (L, ir), v in EXCLUDE.items()},
            "command": "OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 PYTHONPATH=src/python python3 scripts/m2/score_labels.py " + " ".join(sys.argv[1:])}


def run_gate(only=None, scan=None, progress=print) -> dict:
    """Scan and score the labelled sectors; returns {"L20/000_A1m": record, ...}.  `scan(sector, c, E_min, E_max)` defaults to
    rootscan.scan_sector with the fast evaluator; `only` is a list of "L:irrep" strings."""
    from qc3 import constants, fast, rootscan
    from qc3.free_levels import free_levels
    from qc3.sector import Sector
    if scan is None:
        def scan(sec, c, a, b):
            return rootscan.scan_sector(fast.evaluate_point, sec, c, a, b)
    keys = list(rootscan.LABEL_FILES)
    if only:
        wanted = {tuple(o.split(":", 1)) for o in only}
        keys = [k for k in keys if (str(k[0]), k[1]) in wanted]
        if len(keys) != len(wanted):
            raise ValueError(f"--only: unknown sector in {only}; allowed {[f'{L}:{ir}' for L, ir in rootscan.LABEL_FILES]}")
    out = {}
    for Lbyas, irrep in keys:
        sec = Sector.production(irrep, Lbyas=Lbyas, k2_pwave_form=FORM)
        res = scan(sec, tuple(constants.C_BASELINE), constants.GRID_E0, E_MAX)
        rec = score_sector(res, rootscan.read_labels(rootscan.LABEL_FILES[(Lbyas, irrep)]), exclude_ids=EXCLUDE.get((Lbyas, irrep), ()),
                           free_levels=free_levels(sec, constants.GRID_E0 - 1e-3, E_MAX + 1e-3))
        out[f"L{Lbyas}/{irrep}"] = rec
        progress(f"L{Lbyas}/{irrep}: TP={rec['tp']} FP={rec['fp']} FN={rec['fn']} TN={rec['tn']} evals={rec['n_evaluations']} "
                 f"t={rec['seconds']:.1f}s {'OK' if rec['passed'] else 'PROBLEMS: ' + '; '.join(rec['problems'])}")
    return out


def main(argv=None, scan=None, meta=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--out", default=str(DEFAULT_OUT), help="result JSON (default results/m2_label_gate.json)")
    ap.add_argument("--only", nargs="*", default=None, help='subset of sectors, e.g. "20:100_A2 24:000_A1m"')
    args = ap.parse_args(argv)
    t0 = time.perf_counter()
    sectors = run_gate(args.only, scan=scan)
    code = exit_code(sectors)
    doc = {"description": "M2 step 3 label gate: label-free scan_sector vs human labels, coded_v2, C_BASELINE (prompts/003 section 3)",
           **(metadata() if meta is None else meta), "subset": args.only, "gate_passed": code == 0,
           "total_seconds": time.perf_counter() - t0, "sectors": sectors}
    out = pathlib.Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(_jsonable(doc), indent=1))
    print("\n".join(summary_lines(sectors)))
    print(f"gate {'PASSED' if code == 0 else 'FAILED'}; wrote {out}")
    return code


if __name__ == "__main__":
    sys.exit(main())

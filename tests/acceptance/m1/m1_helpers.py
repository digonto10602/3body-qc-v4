"""Shared helpers for the M1 acceptance tests (planner-owned, locked).

Nothing here calls the production interface at import time. Oracle/fast results are cached per process
because one oracle point costs 3-6 s on one core.
"""
import json
import pathlib

import numpy as np

HERE = pathlib.Path(__file__).resolve().parent
ROOT = HERE.parents[2]
ANCHOR_FILE = HERE / "data" / "oracle_anchors_L20.json"

# Production constants (CLAUDE.md invariant 3) — written out here on purpose, not imported from qc3.
PRODUCTION = {
    "xi": 3.444, "atmK": 0.09698, "atmpi": 0.06906, "eta_1": 1.0, "eta_2": 0.5,
    "alpha": 0.5, "epsilon_H": 0.0, "max_shell_num": 20, "waves_1": (0, 1), "waves_2": (0,),
    "parity": -1, "eig_tol": 0.05, "scatter": {"a0_Kpi": 4.04, "a1_Kpi": -43.2, "a0_KK": 4.12},
}
C_BASELINE = (73735.840894011912, -972421.14060757787, 347174.05548116949, -1226756.7068845264)
C_ZERO = (0.0, 0.0, 0.0, 0.0)
IRREPS = ("000_A1m", "100_A2", "110_A2", "111_A2", "200_A2")
NNP = {"000_A1m": (0, 0, 0), "100_A2": (0, 0, 1), "110_A2": (1, 1, 0), "111_A2": (1, 1, 1), "200_A2": (0, 0, 2)}
LBYAS = 20
L_BOX = PRODUCTION["xi"] * LBYAS
GRID_E0, GRID_DE = 0.26310, 0.0969 / 19999

_CACHE = {}


def load_anchors():
    with open(ANCHOR_FILE) as f:
        return json.load(f)


def anchor_points(small_only=False, with_baseline=False):
    pts = load_anchors()["points"]
    if small_only:
        pts = [p for p in pts if p["N"] <= 60]
    if with_baseline:
        pts = [p for p in pts if "K3_baseline" in p]
    return pts


def sector(irrep, **overrides):
    """Production sector for the tests, with the p-wave K2^{-1} form PINNED to "coded_v2" unless overridden.

    Amendment acceptance-v4 (2026-10-22, prompts/001 section 4, prompts/002 section 4): Digonto's decision D7 switches
    the production default to "q0_consistent".  The anchors (data/oracle_anchors_L20.json), the labelled windows and
    every comparison in the locked tests were produced with "coded_v2"; pinning it here keeps those tests testing
    exactly what they tested before.  The production default itself is checked in test_m1_production_default.py,
    which builds its Sector directly (not through this helper)."""
    from qc3.sector import Sector
    overrides.setdefault("k2_pwave_form", "coded_v2")
    return Sector.production(irrep, Lbyas=LBYAS, **overrides)


def _key(backend, irrep, Ecm, c, overrides, options):
    return (backend, irrep, float(Ecm), tuple(float(x) for x in c), tuple(sorted(overrides.items())),
            tuple(sorted((options or {}).items())))


def oracle_point(irrep, Ecm, c=C_ZERO, **overrides):
    key = _key("oracle", irrep, Ecm, c, overrides, None)
    if key not in _CACHE:
        from qc3 import oracle
        _CACHE[key] = oracle.evaluate_point(sector(irrep, **overrides), float(Ecm), tuple(c))
    return _CACHE[key]


def fast_point(irrep, Ecm, c=C_ZERO, options=None, **overrides):
    key = _key("fast", irrep, Ecm, c, overrides, options)
    if key not in _CACHE:
        from qc3 import fast
        _CACHE[key] = fast.evaluate_point(sector(irrep, **overrides), float(Ecm), tuple(c), options=options)
    return _CACHE[key]


def as_c(a):
    return np.asarray(a, dtype=np.complex128)


def rel_frob(a, b):
    """||a-b||_F / ||b||_F (b is the reference)."""
    a = as_c(a); b = as_c(b)
    nb = np.linalg.norm(b)
    return float(np.linalg.norm(a - b) / nb) if nb > 0 else float(np.linalg.norm(a - b))


def E_lab_of(irrep, Ecm):
    P = 2 * np.pi / L_BOX * np.array(NNP[irrep], dtype=float)
    return float(np.sqrt(Ecm ** 2 + P @ P)), P


def basis_list(result):
    """PointResult.basis (structured) -> [(flavour, pvec, ell, m)] for qc3ref."""
    from qc3ref.basis import basis_to_list
    return basis_to_list(result.basis, L_BOX)


def accepted_windows(irrep):
    """Rows of data/windows/L20/L20_<irrep>/accepted_windows.csv as dicts (strings)."""
    import csv
    path = ROOT / "data" / "windows" / "L20" / f"L20_{irrep}" / "accepted_windows.csv"
    with open(path) as f:
        return list(csv.DictReader(f))


def grid_energy(row):
    return GRID_E0 + int(row) * GRID_DE

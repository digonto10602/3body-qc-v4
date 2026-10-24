"""Shared helpers for the M2 acceptance tests (planner-owned, locked with the tests; prompts/003).

Nothing here calls the production interface at import time.  Everything that the tests compare against is
written out here independently of `qc3` (constants, label files, a brute-force free-level enumerator), so that
the tests do not become tautological.
"""
import csv
import itertools
import pathlib

import numpy as np

HERE = pathlib.Path(__file__).resolve().parent
ROOT = HERE.parents[2]
DATA = ROOT / "data"

# Production constants (CLAUDE.md invariant 3), written out on purpose, not imported from qc3.
XI, ATMK, ATMPI = 3.444, 0.09698, 0.06906
IRREPS = ("000_A1m", "100_A2", "110_A2", "111_A2", "200_A2")
NNP = {"000_A1m": (0, 0, 0), "100_A2": (0, 0, 1), "110_A2": (1, 1, 0), "111_A2": (1, 1, 1), "200_A2": (0, 0, 2)}
GRID_E0, GRID_DE, GRID_N = 0.26310, 0.0969 / 19999, 20000
E_MAX = 0.36
# label-generation vector (research 004 section 3.6); the human labels were made with it and with coded_v2
C_BASELINE = (73735.840894011912, -972421.14060757787, 347174.05548116949, -1226756.7068845264)
# v32w best fit (reports/research/003 section 7.1; chi2 = 11.489748255611568 with v32w's own inputs)
C_V32W = (34470.228772793678, -1035173.9140440958, 253921.14757953485, -961680.91711685632)
V32W_CHI2 = 11.489748255611568

# Human true-zero labels (data/README.md; research 004 section 3.5).  208 brackets, 31 true.
LABEL_FILES = {
    (20, "000_A1m"): "v33p_user_truezero_labels_L20_000_A1m_E026310_0360.csv",
    (20, "100_A2"): "v33p_user_truezero_labels_L20_100_A2_E026301_0360.csv",   # old grid start 0.26301; brackets are energies, unaffected
    (20, "110_A2"): "v33p_gpu_user_truezero_labels_L20_110_A2.csv",
    (20, "111_A2"): "v33p_gpu_user_truezero_labels_L20_111_A2.csv",
    (20, "200_A2"): "v33p_gpu_user_truezero_labels_L20_200_A2.csv",
    (24, "000_A1m"): "v33p_user_truezero_labels_L24_000_A1m_E026310_0360.csv",
    (24, "100_A2"): "v33p_user_truezero_labels_L24_100_A2_E026310_0360.csv",
}
LABEL_COUNTS = {   # (brackets, true) per sector, research 002 section 5.2
    (20, "000_A1m"): (10, 3), (20, "100_A2"): (25, 2), (20, "110_A2"): (36, 4), (20, "111_A2"): (30, 5), (20, "200_A2"): (35, 5),
    (24, "000_A1m"): (20, 3), (24, "100_A2"): (52, 9),
}
# L24/100_A2 bracket 36 is documented as unresolved (research 004 section 3.5; prompts/000 M3 step 1): not scored.
LABEL_EXCLUDE = {(24, "100_A2"): (36,)}

# Free-level solutions of the truncated QC seen at the baseline (D8; reports/tests/007 finding 2 for 110_A2)
FREE_LEVEL_110_A2 = 0.319468828          # K(0,0,0) K(0,1,0) pi(1,0,0), P = (1,1,0); quoted to 9 digits in report 007

# Legacy L20 hot-window fit at cutoff 0.335 (results/legacy_L20_passonly_hotwindow, 2026-07-15, coded_v2): its parameter
# vector and its 12 model levels are the only legacy model levels in the repo at a known c (the v32w levels are not in
# the repo; research 003 section 7.1 has only its parameters and chi2).  Used for the code-equivalence check of
# prompts/003 section 4 amendment 2026-10-24 (acceptance-v7).  Its chi2 and errors are invalid (OPEN_ISSUES #1) and
# are never used; the model levels were linear interpolations of det Q sign changes on the 1-row grid.
C_LEGACY_0335 = (32040.841899235049, -1074716.6512728098, 347174.05548116943, -1226756.7068845264)
LEGACY_0335_FIT_LEVELS = ROOT / "results" / "legacy_L20_passonly_hotwindow" / "Ecm_cut_0p335" / "l20_full_cut_0p335_fit_levels_L20.dat"

# A genuine zero of det Q 8.15e-6 (1.7 grid rows) above a pole, 110_A2 at C_LEGACY_0335, coded_v2 (planner scratch
# 2026-10-24, Brent xtol 1e-11 + Delta_h: pole +12.16, zero -21.42).  On the production 4-row coarse grid both fall
# inside one step, so a sign-change scan alone sees neither (prompts/003 section 2.2 item 4b).
PAIR_110_A2_POLE = 0.315947156958689
PAIR_110_A2_ZERO = 0.315955307088961


def read_legacy_0335_fit_levels():
    """[(irrep, level_index, lattice_Ecm, model_Ecm)] from the legacy fit-levels table (12 rows, internal irrep names)."""
    rows = []
    with open(LEGACY_0335_FIT_LEVELS) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            p = line.split()
            rows.append((p[3], int(p[5]), float(p[6]), float(p[8])))
    assert len(rows) == 12
    return rows


def label_path(Lbyas, irrep):
    return DATA / "labels" / LABEL_FILES[(Lbyas, irrep)]


def read_label_rows(Lbyas, irrep):
    """Independent reader of the label CSV: [{"bracket_id": int, "E_left": float, "E_right": float, "label": str}]."""
    with open(label_path(Lbyas, irrep)) as fh:
        rows = []
        for r in csv.DictReader(fh):
            rows.append({"bracket_id": int(r["bracket_id"]), "E_left": float(r["E_left"]), "E_right": float(r["E_right"]),
                         "label": r["user_label"].strip().lower()})
    return rows


def sector(irrep, Lbyas=20, **overrides):
    """Production sector; the p-wave form must be given explicitly by every test (no silent default)."""
    from qc3.sector import Sector
    if "k2_pwave_form" not in overrides:
        raise ValueError("m2 tests must pin k2_pwave_form explicitly")
    return Sector.production(irrep, Lbyas=Lbyas, **overrides)


FREE_MERGE_TOL = 1e-11        # states whose E* agree within this are one level (prompts/003 section 2.1 convention)
FREE_DEGENERACY_SPREAD = 1e-13  # members of one level must agree to this (symmetry-degenerate states are bit-identical)


def brute_force_free_levels(irrep, Lbyas, E_min, E_max, nmax=6):
    """Non-interacting K K pi energies with total momentum d = NNP[irrep], any irrep, sorted distinct E* in [E_min, E_max].
    E_lab = sum_a sqrt(m_a^2 + q^2 n_a^2), n_1 + n_2 + n_3 = d, q = 2 pi / (xi Lbyas); E* = sqrt(E_lab^2 - q^2 d^2).
    Identical kaons: unordered pair {n_1, n_2}.  Independent of qc3.free_levels.

    Grouping (acceptance-v7, prompts/003 section 4 amendment 2026-10-24): the states are sorted by E* and consecutive
    energies closer than FREE_MERGE_TOL form one level.  The reported E* of a level is the exact float64 energy of its
    first member -- NOT a rounded value (acceptance-v6 rounded the key to 1e-11, which quantised the reference by up to
    5e-12 and made the 1e-12 comparison unsatisfiable for an exact implementation).  Symmetry-degenerate states give
    bit-identical E* (same multiset of |n_a|^2, commutative float additions), which is asserted (FREE_DEGENERACY_SPREAD)."""
    q = 2.0 * np.pi / (XI * Lbyas)
    d = np.asarray(NNP[irrep], dtype=int)
    rng = range(-nmax, nmax + 1)
    rows = []
    for n1 in itertools.product(rng, repeat=3):
        n1 = np.asarray(n1)
        w1 = np.sqrt(ATMK ** 2 + q * q * float(n1 @ n1))
        if w1 > E_max:            # the other two particles add at least their masses
            continue
        for n2 in itertools.product(rng, repeat=3):
            if tuple(n2) < tuple(n1):
                continue
            n2 = np.asarray(n2)
            n3 = d - n1 - n2
            E_lab = w1 + np.sqrt(ATMK ** 2 + q * q * float(n2 @ n2)) + np.sqrt(ATMPI ** 2 + q * q * float(n3 @ n3))
            Es = np.sqrt(E_lab ** 2 - q * q * float(d @ d))
            if E_min <= Es <= E_max:
                rows.append((float(Es), (tuple(int(x) for x in n1), tuple(int(x) for x in n2), tuple(int(x) for x in n3))))
    rows.sort(key=lambda t: t[0])
    groups = []                   # [(list of exact E*, list of states)]
    for E, st in rows:
        if groups and E - groups[-1][0][-1] <= FREE_MERGE_TOL:
            groups[-1][0].append(E)
            groups[-1][1].append(st)
        else:
            groups.append(([E], [st]))
    for g_E, _ in groups:
        assert g_E[-1] - g_E[0] <= FREE_DEGENERACY_SPREAD, ("free-level reference: a merged group is not exactly degenerate", g_E)
    Es = np.array([g_E[0] for g_E, _ in groups], dtype=float)
    states = [st for _, st in groups]
    return Es, states


def grid_energy(row):
    return GRID_E0 + int(row) * GRID_DE


def lambda_min_cov(cov):
    return float(np.linalg.eigvalsh(np.asarray(cov)).min())


def chi2_bound_from_root_error(chi2, cov, dE):
    """Rigorous bound on |delta chi2| when every model energy moves by at most dE (Cauchy-Schwarz in the C^-1 metric):
    |delta chi2| <= 2 sqrt(chi2) sqrt(dr^T C^-1 dr) + dr^T C^-1 dr, with dr^T C^-1 dr <= m dE^2 / lambda_min(C)."""
    m = np.asarray(cov).shape[0]
    q = m * dE ** 2 / lambda_min_cov(cov)
    return 2.0 * np.sqrt(max(chi2, 0.0)) * np.sqrt(q) + q

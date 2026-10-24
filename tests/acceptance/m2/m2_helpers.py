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


def brute_force_free_levels(irrep, Lbyas, E_min, E_max, nmax=6):
    """Non-interacting K K pi energies with total momentum d = NNP[irrep], any irrep, sorted distinct E* in [E_min, E_max].
    E_lab = sum_a sqrt(m_a^2 + q^2 n_a^2), n_1 + n_2 + n_3 = d, q = 2 pi / (xi Lbyas); E* = sqrt(E_lab^2 - q^2 d^2).
    Identical kaons: unordered pair {n_1, n_2}.  Independent of qc3.free_levels."""
    q = 2.0 * np.pi / (XI * Lbyas)
    d = np.asarray(NNP[irrep], dtype=int)
    rng = range(-nmax, nmax + 1)
    vals = {}
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
                key = round(Es, 11)
                vals.setdefault(key, []).append((tuple(int(x) for x in n1), tuple(int(x) for x in n2), tuple(int(x) for x in n3)))
    Es = np.array(sorted(vals.keys()), dtype=float)
    states = [vals[k] for k in sorted(vals.keys())]
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

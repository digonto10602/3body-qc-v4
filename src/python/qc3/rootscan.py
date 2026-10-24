"""Label-free, complete root finder per sector (prompts/003 section 2.2; D8; docs/METHODS.md section 12).

STEP 1 STUBS: the data definitions below are real, every function raises NotImplementedError until prompts/003 step 2
(scan_sector) / step 3 (score_labels as used by the gate script) implements it.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .sector import Sector

SCAN_DEFAULTS = dict(step_rows=4.0, fine_step_rows=0.25, edge=2e-6, h=1e-6, h2=3e-6, D=4.0, xtol=1e-11, min_xatol=1e-8,
                     free_tol=1e-7, order_tol=0.3)
KINDS = ("zero", "pole", "even_zero", "free_level", "ambiguous")

# (Lbyas, irrep) -> path relative to qc3.lattice.default_data_dir(): the 7 human true-zero label files of data/labels
LABEL_FILES: dict = {
    (20, "000_A1m"): "labels/v33p_user_truezero_labels_L20_000_A1m_E026310_0360.csv",
    (20, "100_A2"): "labels/v33p_user_truezero_labels_L20_100_A2_E026301_0360.csv",
    (20, "110_A2"): "labels/v33p_gpu_user_truezero_labels_L20_110_A2.csv",
    (20, "111_A2"): "labels/v33p_gpu_user_truezero_labels_L20_111_A2.csv",
    (20, "200_A2"): "labels/v33p_gpu_user_truezero_labels_L20_200_A2.csv",
    (24, "000_A1m"): "labels/v33p_user_truezero_labels_L24_000_A1m_E026310_0360.csv",
    (24, "100_A2"): "labels/v33p_user_truezero_labels_L24_100_A2_E026310_0360.csv",
}


@dataclass
class RootCandidate:
    E: float
    kind: str                  # one of KINDS
    order: float | None        # estimated order of the zero (m of section 2.2 item 5); None for poles
    order_int: int | None
    sign_change: bool
    discriminant: float        # Delta_h
    slope: float | None        # u^T Q'(E) u for zero-like kinds, None for poles
    N: int
    segment: tuple
    bracket: tuple
    via: str                   # "sign" | "min" | "free"
    free_level: float | None
    n_evaluations: int


@dataclass
class ScanResult:
    sector: Sector
    c: tuple
    E_min: float
    E_max: float
    candidates: list           # list[RootCandidate], sorted by E
    jumps: np.ndarray
    segments: list
    n_evaluations: int
    n_grid: int
    seconds: float
    options: dict

    def zeros(self) -> np.ndarray:
        """E of kind == "zero", sorted."""
        raise NotImplementedError("qc3.rootscan.ScanResult.zeros: prompts/003 step 2")

    def poles(self) -> np.ndarray:
        raise NotImplementedError("qc3.rootscan.ScanResult.poles: prompts/003 step 2")

    def of_kind(self, kind: str) -> list:
        raise NotImplementedError("qc3.rootscan.ScanResult.of_kind: prompts/003 step 2")


def scan_sector(evaluate, sector: Sector, c, E_min: float, E_max: float, options: dict | None = None) -> ScanResult:
    raise NotImplementedError("qc3.rootscan.scan_sector: prompts/003 step 2")


def read_labels(path) -> list:
    """Rows {"bracket_id": int, "E_left": float, "E_right": float, "label": "true"|"false"|""}; `path` absolute or
    relative to qc3.lattice.default_data_dir() (as the LABEL_FILES values)."""
    raise NotImplementedError("qc3.rootscan.read_labels: prompts/003 step 2")


def score_labels(result: ScanResult, labels: list, slack_rows: float = 0.25, exclude_ids=()) -> dict:
    """{"tp","fp","fn","tn": int, "mismatches": list[dict], "unlabelled_zeros": list[RootCandidate], "excluded": list[int]}"""
    raise NotImplementedError("qc3.rootscan.score_labels: prompts/003 step 2")

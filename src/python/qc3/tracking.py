"""Root tracking across fit-parameter changes (prompts/003 section 2.5; docs/METHODS.md section 13).

STEP 1 STUBS: the data definitions below are real; every function and method raises NotImplementedError until
prompts/003 step 4 implements it.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from . import fast as _fast
from .point import MissingRootError
from .rootscan import RootCandidate, ScanResult  # noqa: F401
from .sector import Sector

TRACK_DEFAULTS = dict(h=1e-6, D=4.0, xtol=1e-11, max_iter=12, max_move_rows=20.0, safeguard_factor=10.0, safeguard_abs=1e-6,
                      local_rows=60, local_steps_rows=(4.0, 1.0, 0.25), park_margin=0.005)


class RootSetChangedError(MissingRootError):
    """A full scan at the current c disagrees with the tracked root set (count, or a root differs by > 1e-9)."""


@dataclass
class TrackedRoot:
    E: float
    J: np.ndarray              # dE0/dc_i, shape (4,)
    slope: float
    discriminant: float
    N: int
    segment: tuple


def root_jacobian(evaluate, sector: Sector, c, E0: float, h=1e-6, D=4.0, lnabsdet_at_root: float | None = None) -> TrackedRoot:
    """Eq. (1): J_i = -u^T B_i u / u^T Q'(E0) u; 2 evaluations at E0 +- h (+1 at E0 if lnabsdet_at_root is None)."""
    raise NotImplementedError("qc3.tracking.root_jacobian: prompts/003 step 4")


def secant_root(evaluate, sector: Sector, c, E_pred: float, scale: float, xtol=1e-11, max_iter=12) -> tuple:
    """(E, lnabsdet of the last evaluated iterate, n_evaluations); never evaluates at the returned E."""
    raise NotImplementedError("qc3.tracking.secant_root: prompts/003 step 4")


class SectorRoots:
    """Tracked zeros of one sector.  Attributes: sector, c, E_min, E_max, scan, roots (incl. parked), free_level_solutions,
    n_evaluations, n_local_rescans, n_full_rescans, n_substeps."""

    def __init__(self, sector: Sector, c0, E_min: float, E_max: float, evaluate=_fast.evaluate_point,
                 scan_options: dict | None = None, track_options: dict | None = None):
        raise NotImplementedError("qc3.tracking.SectorRoots.__init__: prompts/003 step 4")

    def update(self, c) -> np.ndarray:
        raise NotImplementedError("qc3.tracking.SectorRoots.update: prompts/003 step 4")

    def energies(self) -> np.ndarray:
        raise NotImplementedError("qc3.tracking.SectorRoots.energies: prompts/003 step 4")

    def jacobian(self) -> np.ndarray:
        raise NotImplementedError("qc3.tracking.SectorRoots.jacobian: prompts/003 step 4")

    def rescan(self) -> ScanResult:
        raise NotImplementedError("qc3.tracking.SectorRoots.rescan: prompts/003 step 4")

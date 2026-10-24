"""Point result container (prompts/001 section 5.3)."""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .sector import Sector

BASIS_DTYPE = np.dtype([("flavour", "i4"), ("nx", "i4"), ("ny", "i4"), ("nz", "i4"),
                        ("ell", "i4"), ("m", "i4")])
ALL = ("F", "G", "K2inv", "F3", "F3inv", "Vsel", "K3", "Q")


class MissingRootError(RuntimeError):
    """Raised when a root is not bracketed / fewer QC levels than lattice levels."""


class SingularQError(ValueError):
    """The fast path's LU of Q met an exactly zero pivot (C++ "fast: Q is singular"): E is an exact root of det Q
    (prompts/003 section 2.5; never evaluate exactly at a converged root, but if it happens it is a documented case)."""


@dataclass
class PointResult:
    sector: Sector
    Ecm: float
    E_lab: float
    N1: int
    N2: int
    N: int
    Nproj: int
    basis: np.ndarray            # BASIS_DTYPE, length N, order of plan section 2.1
    c: tuple                     # (Kiso0, Kiso1, KB, KE)
    F: np.ndarray | None
    G: np.ndarray | None
    K2inv: np.ndarray | None
    F3: np.ndarray | None
    F3inv: np.ndarray | None
    Vsel: np.ndarray | None      # N x Nproj
    K3_basis: dict | None        # {"iso0","iso1","B","E"} -> N x N
    K3: np.ndarray | None
    Q: np.ndarray | None         # Nproj x Nproj
    sign: int
    lnabsdet: float
    metadata: dict
    timings: dict


@dataclass
class LevelSet:
    labels: list      # list[tuple[str, int]]
    irreps: list
    nP: list          # list[tuple]
    files: list
    Ncfg: int
    Ecm: np.ndarray
    cov: np.ndarray
    Ecm_cutoff: float
    Lbyas: int

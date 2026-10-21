"""Sector definition (prompts/001 section 5.2)."""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .constants import NNP


@dataclass(frozen=True)
class Sector:
    irrep: str
    Lbyas: int = 20
    xi: float = 3.444
    atmK: float = 0.09698
    atmpi: float = 0.06906
    eta_1: float = 1.0
    eta_2: float = 0.5
    alpha: float = 0.5
    epsilon_H: float = 0.0
    max_shell_num: int = 20
    waves_1: tuple = (0, 1)
    waves_2: tuple = (0,)
    parity: int = -1
    eig_tol: float = 0.05
    a0_Kpi: float = 4.04
    a1_Kpi: float = -43.2
    a0_KK: float = 4.12            # r_l = P_l = 0
    k2_pwave_form: str = "coded_v2"  # or "q0_consistent" (plan section 2.2)

    @classmethod
    def production(cls, irrep: str, Lbyas: int = 20, **overrides) -> "Sector":
        """Production sector (CLAUDE.md invariant 3); `overrides` replace individual fields."""
        return cls(irrep=irrep, Lbyas=Lbyas, **overrides)

    @property
    def L(self) -> float:  # xi * Lbyas  (box size in units of a_t)
        return self.xi * self.Lbyas

    @property
    def nnP(self) -> tuple:  # integer d vector
        try:
            return tuple(NNP[self.irrep])
        except KeyError:
            raise ValueError(f"unknown irrep {self.irrep!r}; known: {sorted(NNP)}") from None

    @property
    def P_lab(self) -> np.ndarray:  # 2*pi/L * d
        return 2.0 * np.pi / self.L * np.asarray(self.nnP, dtype=float)

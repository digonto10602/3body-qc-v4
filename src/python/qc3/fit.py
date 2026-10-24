"""chi^2 of the tracked QC levels and the two-stage L20 fit (prompts/003 sections 2.6-2.7; docs/METHODS.md section 14).

STEP 1 STUBS: the data definitions below are real; every function and method raises NotImplementedError until
prompts/003 step 5 (QCChi2, make_trackers) / step 6 (fit_cutoff, verify_fit, cutoff_scan, save_results, CLI) implements it.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from . import fast as _fast
from .constants import GRID_E0, IRREPS, K2_PWAVE_FORM_DEFAULT
from .point import LevelSet

# v32w best fit (reports/research/003 section 7.1; chi2 = 11.489748255611568 with v32w's own inputs, not reproducible here)
C_V32W = (34470.228772793678, -1035173.9140440958, 253921.14757953485, -961680.91711685632)
V32W_CHI2 = 11.489748255611568


class QCChi2:
    errordef = 1.0

    def __init__(self, levels: LevelSet, trackers: dict, penalty_sigmas: float = 100.0):
        raise NotImplementedError("qc3.fit.QCChi2.__init__: prompts/003 step 5")

    def __call__(self, c) -> float:
        raise NotImplementedError("qc3.fit.QCChi2.__call__: prompts/003 step 5")

    def grad(self, c) -> np.ndarray:
        """Eq. (2), shape (4,)."""
        raise NotImplementedError("qc3.fit.QCChi2.grad: prompts/003 step 5")

    def model(self, c) -> tuple:
        """(E_model aligned with levels.labels, J_model (m, 4))."""
        raise NotImplementedError("qc3.fit.QCChi2.model: prompts/003 step 5")


@dataclass
class FitResult:
    cutoff: float
    k2_pwave_form: str
    Lbyas: int
    labels: list
    E_lat: np.ndarray
    sigma: np.ndarray
    cov_lattice: np.ndarray
    E_model: np.ndarray
    residual_over_sigma: np.ndarray
    c_start: tuple
    c: np.ndarray
    errors: np.ndarray
    errors_gauss_newton: np.ndarray
    covariance: np.ndarray
    correlation: np.ndarray
    chi2: float
    chi2_start: float
    chi2_lm: float
    ndof: int
    chi2_per_dof: float
    edm: float
    valid: bool
    accurate: bool
    nfcn: int
    ngrad: int
    n_evaluations: int
    n_local_rescans: int
    n_full_rescans: int
    n_penalised_calls: int
    n_penalised_calls_final_stage: int
    seconds: float
    minimizer: str
    verification: dict | None


def make_trackers(irreps, c0, Lbyas=20, k2_pwave_form=K2_PWAVE_FORM_DEFAULT, E_min=GRID_E0, E_max=0.36,
                  evaluate=_fast.evaluate_point) -> dict:
    raise NotImplementedError("qc3.fit.make_trackers: prompts/003 step 5")


def fit_cutoff(levels: LevelSet, trackers: dict, c_start, *, minimizer="lm+migrad", strategy=1, ncall=2000,
               lm_options=None, fixed=(), verify=True) -> FitResult:
    raise NotImplementedError("qc3.fit.fit_cutoff: prompts/003 step 6")


def verify_fit(fit: FitResult, trackers: dict, levels: LevelSet) -> dict:
    """{"counts_equal","max_root_diff","chi2_fresh","chi2_diff","n_evaluations","seconds","per_sector"}"""
    raise NotImplementedError("qc3.fit.verify_fit: prompts/003 step 6")


def cutoff_scan(data_dir, Lbyas, cutoffs, c_start, k2_pwave_form, irreps=IRREPS, E_min=GRID_E0, E_max=0.36, out_json=None,
                verify=True, fixed=(), evaluate=_fast.evaluate_point) -> list:
    raise NotImplementedError("qc3.fit.cutoff_scan: prompts/003 step 6")


def save_results(fits: list, path, extra: dict | None = None) -> None:
    raise NotImplementedError("qc3.fit.save_results: prompts/003 step 6")


def main(argv=None) -> int:
    """CLI: python -m qc3.fit --cutoffs 0.315 0.325 0.335 0.345 --form q0_consistent --out results/fit_L20_..._<tag>.json
    [--fix KB KE] [--c-start baseline|v32w]"""
    raise NotImplementedError("qc3.fit.main: prompts/003 step 6")


if __name__ == "__main__":
    raise SystemExit(main())

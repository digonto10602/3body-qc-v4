"""chi^2 of the tracked QC levels and the two-stage L20 fit (prompts/003 sections 2.6-2.7; docs/METHODS.md section 14).

Step 5 (done): QCChi2, make_trackers.  STEP 6 STUBS: fit_cutoff, verify_fit, cutoff_scan, save_results and main raise
NotImplementedError until prompts/003 step 6 implements them.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from . import fast as _fast
from . import lattice as _lattice
from .constants import GRID_E0, IRREPS, K2_PWAVE_FORM_DEFAULT
from .point import LevelSet
from .sector import Sector
from .tracking import SectorRoots as _SectorRoots

# v32w best fit (reports/research/003 section 7.1; chi2 = 11.489748255611568 with v32w's own inputs, not reproducible here)
C_V32W = (34470.228772793678, -1035173.9140440958, 253921.14757953485, -961680.91711685632)
V32W_CHI2 = 11.489748255611568


class QCChi2:
    """chi^2(c) = r^T C^+ r, r_i = E_lat_i - E_model_i(c) (prompts/003 section 2.6; docs/METHODS.md section 14).

    * covariance: the D1 jackknife covariance of `levels` and its eigenvalue pseudo-inverse, exactly as `qc3.lattice.chi2`;
    * model levels: per sector `trackers[irrep].update(c)` -> `energies()` (D8: free-level solutions are never in it) -> D3
      (`qc3.lattice.match_levels`: sort both lists, the lowest m QC levels for the m lattice levels, n >= m allowed);
    * a sector with n < m QC zeros does not abort and does not return 1e100: its (m - n) highest lattice levels get the
      finite residual `penalty_sigmas * sigma_k` (E_model = E_lat - penalty_sigmas * sigma_k) and a zero gradient row;
      `n_penalised_calls` counts the model evaluations in which this happened (prompts/003 section 8);
    * D9 (flag only): tracked zeros with slope < 0 stay in the matching and are listed in `wrong_sign_zeros`
      ({"irrep","E","slope","matched","level_index"}; level_index indexes `levels.labels`), `n_calls_with_wrong_sign` counts the
      model evaluations in which any occurred;
    * gradient, eq. (2): dchi^2/dc = -2 J_model^T C^+ r, J_model = dE_model/dc from the trackers (penalised rows are 0).

    The model of the most recent parameter vector is cached (so __call__ followed by grad at the same c, as Migrad does, updates
    the trackers and the counters once).  Errors of the trackers (AmbiguousRootError, RootSetChangedError, ...) propagate."""
    errordef = 1.0

    def __init__(self, levels: LevelSet, trackers: dict, penalty_sigmas: float = 100.0, eig_rel_cut: float = 1e-12):
        from . import lattice
        self.levels = levels
        self.trackers = trackers
        self.penalty_sigmas = float(penalty_sigmas)
        if not self.penalty_sigmas > 0.0:
            raise ValueError("penalty_sigmas must be > 0")
        self.labels = [tuple(l) for l in levels.labels]
        self.m = len(self.labels)
        if self.m == 0:
            raise ValueError("QCChi2: the LevelSet has no levels")
        self.E_lat = np.asarray(levels.Ecm, dtype=np.float64)
        cov = np.asarray(levels.cov, dtype=np.float64)
        if cov.shape != (self.m, self.m) or self.E_lat.shape != (self.m,):
            raise ValueError(f"QCChi2: levels inconsistent: {self.m} labels, Ecm {self.E_lat.shape}, cov {cov.shape}")
        self.cov = cov
        self._w, self._U = lattice.cov_eig_pinv(cov, eig_rel_cut)       # chi^2 = sum((U^T r)^2 / w), same as lattice.chi2
        self.sigma = np.sqrt(np.diag(cov))
        self.Cinv = (self._U / self._w) @ self._U.T
        self.irreps = []                                              # sectors that have levels, in first-appearance order
        for irrep, _ in self.labels:
            if irrep not in self.irreps:
                self.irreps.append(irrep)
        self.index = {irrep: [i for i, l in enumerate(self.labels) if l[0] == irrep] for irrep in self.irreps}
        for irrep in self.irreps:
            if irrep not in trackers:
                raise ValueError(f"QCChi2: no tracker for sector {irrep!r} (trackers: {sorted(trackers)})")
            if int(trackers[irrep].sector.Lbyas) != int(levels.Lbyas):
                raise ValueError(f"QCChi2: tracker {irrep} has L/a_s={trackers[irrep].sector.Lbyas}, levels have {levels.Lbyas}")
        self.n_calls = self.n_grad_calls = self.n_penalised_calls = self.n_calls_with_wrong_sign = 0
        self.wrong_sign_zeros: list = []
        self.penalised: np.ndarray = np.zeros(self.m, dtype=bool)     # levels penalised in the last model evaluation
        self._cache_key = None
        self._cache = None

    # ------------------------------------------------------------ model
    def _model(self, c):
        c = np.asarray([float(x) for x in c], dtype=np.float64)
        if c.shape != (4,):
            raise ValueError("c must have 4 components (Kiso0, Kiso1, KB, KE)")
        key = c.tobytes()
        if self._cache_key == key and all(np.array_equal(self.trackers[ir].c, c) for ir in self.irreps):
            return self._cache
        E_model = np.empty(self.m)
        J_model = np.zeros((self.m, 4))
        penalised = np.zeros(self.m, dtype=bool)
        wrong: list = []
        for irrep in self.irreps:
            tr = self.trackers[irrep]
            E_qc = np.asarray(tr.update(c), dtype=np.float64)
            J_qc = np.asarray(tr.jacobian(), dtype=np.float64).reshape(len(E_qc), 4)
            idx = np.asarray(self.index[irrep])
            order = np.argsort(self.E_lat[idx], kind="stable")
            n_match = min(len(idx), len(E_qc))
            matched = idx[order[:n_match]]                     # D3: the lowest n_match lattice levels get QC levels
            pairs = _lattice.match_levels(self.E_lat[matched], E_qc)
            matched_qc = {}
            for a, b in pairs:
                i = int(matched[a])
                E_model[i], J_model[i] = E_qc[b], J_qc[b]
                matched_qc[b] = i
            for i in idx[order[n_match:]]:                     # fewer zeros than levels: finite penalty, zero gradient
                E_model[i] = self.E_lat[i] - self.penalty_sigmas * self.sigma[i]
                penalised[i] = True
            for r in tr.wrong_sign_roots():
                b = int(np.flatnonzero(E_qc == r.E)[0])        # r.E is an entry of energies() (same floats)
                wrong.append({"irrep": irrep, "E": float(r.E), "slope": float(r.slope), "matched": b in matched_qc,
                              "level_index": matched_qc.get(b)})
        self.penalised = penalised
        self.wrong_sign_zeros = wrong
        if penalised.any():
            self.n_penalised_calls += 1
        if wrong:
            self.n_calls_with_wrong_sign += 1
        self._cache_key, self._cache = key, (E_model, J_model)
        return self._cache

    def model(self, c) -> tuple:
        """(E_model aligned with levels.labels, J_model (m, 4)); penalised levels: E_model = E_lat - penalty*sigma, J row 0."""
        E, J = self._model(c)
        return E.copy(), J.copy()

    def residual(self, c) -> np.ndarray:
        """r = E_lat - E_model(c) (penalised levels: +penalty_sigmas * sigma)."""
        return self.E_lat - self._model(c)[0]

    def __call__(self, c) -> float:
        self.n_calls += 1
        r = self.residual(c)
        return float(np.sum((self._U.T @ r) ** 2 / self._w))

    def grad(self, c) -> np.ndarray:
        """Eq. (2): dchi^2/dc = -2 J_model^T C^+ r, shape (4,)."""
        self.n_grad_calls += 1
        _, J = self._model(c)
        return -2.0 * (J.T @ (self.Cinv @ self.residual(c)))

    def verify_roots(self) -> dict:
        """Explicit verification (prompts/003 section 2.5 item 6): a fresh full scan of every sector at its current c
        (`SectorRoots.rescan`).  A zero that entered the range from above E_max is caught only here.  Raises
        RootSetChangedError / AmbiguousRootError; returns {irrep: ScanResult}."""
        return {irrep: self.trackers[irrep].rescan() for irrep in self.irreps}

    @property
    def n_evaluations(self) -> int:
        """QC evaluations made by the trackers of the sectors in use (all calls)."""
        return sum(self.trackers[ir].n_evaluations for ir in self.irreps)


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
    """{irrep: SectorRoots}: production sectors (invariant 3) with the given p-wave form, initial label-free scan at c0 over
    [E_min, E_max]."""
    return {irrep: _SectorRoots(Sector.production(irrep, Lbyas, k2_pwave_form=k2_pwave_form), c0, E_min, E_max, evaluate=evaluate)
            for irrep in irreps}


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

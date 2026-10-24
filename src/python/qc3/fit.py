"""chi^2 of the tracked QC levels and the two-stage L20 fit (prompts/003 sections 2.6-2.7; docs/METHODS.md section 14).

* `QCChi2` (step 5): chi^2 = r^T C^+ r of the tracked, D3-matched, D8-filtered QC zeros, analytic gradient, D9 wrong-sign flags.
* `fit_cutoff` (step 6): stage 1 Levenberg-Marquardt (scipy `lm`, whitened residual, analytic Jacobian), stage 2 iminuit Migrad +
  Hesse from the LM minimum; reports the `unconstrained` flags.  Never Migrad first (prompts/003 section 8).
* `verify_fit`: a fresh label-free scan of every sector at the fitted parameters; compares tracked and fresh zeros up to the highest
  matched root plus VERIFY_MARGIN (a zero can enter the tracked range from above E_max; it matters only if it lands below the
  highest matched level).
* `cutoff_scan`, `save_results`, `main` (`python -m qc3.fit`), `code_equivalence` (`--code-equivalence`).

Tracked range: for a cutoff E_c the trackers cover [E_min, E_c + TRACK_MARGIN].  Why 0.005 (= ~10 sigma of a level, sigma ~ 5e-4):
a root pushed above E_c + TRACK_MARGIN (+ the tracker's own park margin, where it is still followed by the tangent) while it is one
of the m matched ones costs the fit a residual of >= TRACK_MARGIN ~ 10 sigma on the highest level, i.e. chi^2 >= ~100, while the
finite penalty for a missing level is 100 sigma (+1e4); so the penalty only acts where the fit is already far from any minimum.
The cost of a chi^2 call grows with E_max (82-109 ms at 0.34, 183-255 ms at 0.36; reports/tests/017).
"""
from __future__ import annotations

import argparse
import json
import sys
import time
from dataclasses import dataclass, field, fields

import numpy as np

from . import fast as _fast
from . import lattice as _lattice
from .constants import C_BASELINE, GRID_DE, GRID_E0, IRREPS, K2_PWAVE_FORM_DEFAULT
from .point import LevelSet, MissingRootError
from .roots import AmbiguousRootError
from .sector import Sector
from .tracking import SectorRoots as _SectorRoots

# v32w best fit (reports/research/003 section 7.1; chi2 = 11.489748255611568 with v32w's own inputs, not reproducible here)
C_V32W = (34470.228772793678, -1035173.9140440958, 253921.14757953485, -961680.91711685632)
V32W_CHI2 = 11.489748255611568
# legacy hot-window fit at cutoff 0.335 (results/legacy_L20_passonly_hotwindow; its chi^2 and errors are invalid, OPEN_ISSUES #1)
C_LEGACY_0335 = (32040.841899235049, -1074716.6512728098, 347174.05548116943, -1226756.7068845264)
PARAM_NAMES = ("Kiso0", "Kiso1", "KB", "KE")
TRACK_MARGIN = 0.005       # trackers cover [E_min, cutoff + TRACK_MARGIN]; justification in the module docstring
VERIFY_MARGIN = 0.002      # verify_fit compares tracked and fresh zeros up to (highest matched root + VERIFY_MARGIN)
EDM_MAX = 2e-3             # Migrad: EDM < 10 x 0.002 * tol * errordef (iminuit documents that Migrad may exceed edm_max by 10x)
VERIFY_ROOT_TOL = 1e-9


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

    def _match_sector(self, irrep, E_qc):
        """D3 for one sector: [(level index i into levels.labels, index b into E_qc)]; the lowest min(m_irrep, n) lattice levels
        get the lowest QC levels (n >= m allowed); the remaining (highest) lattice levels are left to the penalty."""
        idx = np.asarray(self.index[irrep])
        order = np.argsort(self.E_lat[idx], kind="stable")
        matched = idx[order[:min(len(idx), len(E_qc))]]
        return [(int(matched[a]), int(b)) for a, b in _lattice.match_levels(self.E_lat[matched], E_qc)]

    def chi2_from_roots(self, roots: dict) -> tuple:
        """(chi^2, E_model) for given QC zeros {irrep: array} with the same D3 matching and penalty as the model (no trackers
        involved; used by `verify_fit` for the fresh-scan roots)."""
        E_model = np.empty(self.m)
        for irrep in self.irreps:
            E_qc = np.asarray(roots[irrep], dtype=np.float64)
            idx = np.asarray(self.index[irrep])
            order = np.argsort(self.E_lat[idx], kind="stable")
            for i, b in self._match_sector(irrep, E_qc):
                E_model[i] = E_qc[b]
            for i in idx[order[min(len(idx), len(E_qc)):]]:
                E_model[i] = self.E_lat[i] - self.penalty_sigmas * self.sigma[i]
        r = self.E_lat - E_model
        return float(np.sum((self._U.T @ r) ** 2 / self._w)), E_model

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
            matched_qc = {}
            for i, b in self._match_sector(irrep, E_qc):
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
    """One fit at one cutoff.  `n_evaluations` = all QC evaluations of the trackers (initial scans included); `seconds` = wall time
    of the LM + Migrad + Hesse stages (the initial scans are `seconds_scan`, the verification is `verification["seconds"]`).
    Fixed parameters: c[i] == c_start[i] exactly; errors / errors_gauss_newton entries 0; covariance and correlation rows and
    columns 0 (correlation[i, i] == 0)."""
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
    # amendment v7 (prompts/003 section 2.6 D9, section 2.7 item 4)
    fixed: tuple = ()
    n_free: int = 4
    unconstrained: bool = False
    unconstrained_reasons: list = field(default_factory=list)
    wrong_sign_zeros: list = field(default_factory=list)
    n_wrong_sign_zeros: int = 0
    # additions (step 6): cost and diagnostics
    n_chi2_calls: int = 0               # QCChi2.__call__ count over the whole fit (LM's own residual evaluations are lm_nfev)
    n_calls_with_wrong_sign: int = 0    # model evaluations in which any tracked zero had slope < 0 (D9)
    lm_status: int = 0                  # scipy least_squares status (1-4 converged, 0 max_nfev, -1 failure)
    lm_nfev: int = 0
    lm_njev: int = 0
    seconds_scan: float = 0.0           # building the trackers (initial label-free scans), set by cutoff_scan
    seconds_lm: float = 0.0
    seconds_migrad: float = 0.0         # Migrad + Hesse
    E_max_tracked: float = 0.0
    migrad_calls: int = 1               # 2 if Migrad had to be restarted once because fmin.is_valid was False
    errors_fd_hessian: np.ndarray = field(default_factory=lambda: np.zeros(4))   # sqrt(diag(2 H^-1)), H = central differences of the analytic gradient (0.02 sigma steps): independent check of the Hesse errors
    fd_hessian_rel: float | None = None       # step (in units of sigma) that gave errors_fd_hessian; None if all failed
    fd_hessian_note: str = ""                 # why a larger step failed (e.g. a level moved across a basis-dimension jump)
    ambiguous_trials_lm: list = field(default_factory=list)   # LM trial points rejected because a scan candidate was ambiguous (see fit_cutoff)


def make_trackers(irreps, c0, Lbyas=20, k2_pwave_form=K2_PWAVE_FORM_DEFAULT, E_min=GRID_E0, E_max=0.36,
                  evaluate=_fast.evaluate_point) -> dict:
    """{irrep: SectorRoots}: production sectors (invariant 3) with the given p-wave form, initial label-free scan at c0 over
    [E_min, E_max]."""
    return {irrep: _SectorRoots(Sector.production(irrep, Lbyas, k2_pwave_form=k2_pwave_form), c0, E_min, E_max, evaluate=evaluate)
            for irrep in irreps}


def _free_indices(fixed) -> tuple:
    fixed = tuple(fixed)
    unknown = [n for n in fixed if n not in PARAM_NAMES]
    if unknown:
        raise ValueError(f"fixed: unknown parameter name(s) {unknown}; allowed: {PARAM_NAMES}")
    if len(set(fixed)) != len(fixed):
        raise ValueError(f"fixed: duplicate names in {fixed}")
    free = tuple(i for i, n in enumerate(PARAM_NAMES) if n not in fixed)
    if not free:
        raise ValueError("fixed: all four parameters are fixed")
    return fixed, free


def _gn_errors(fcn: QCChi2, J: np.ndarray, free) -> np.ndarray:
    """Gauss-Newton errors sqrt(diag((J^T C^+ J)^-1)) = sqrt(diag(2 H_GN^-1)) of the free parameters (length 4, zeros for fixed)."""
    Jf = J[:, list(free)]
    A = Jf.T @ fcn.Cinv @ Jf
    out = np.zeros(4)
    out[list(free)] = np.sqrt(np.diag(np.linalg.inv(A)))
    return out


def _fd_hessian_errors(fcn: QCChi2, c, errors, free, rels=(0.02, 0.002)) -> tuple:
    """Independent cross-check of the Hesse errors (a diagnostic, not a fit input): central differences of the ANALYTIC gradient
    with steps rel * sigma_i, H symmetrised, errors sqrt(diag(2 H^-1)) of the free parameters.  Returns (errors (length 4, zeros for
    fixed), rel used, message).  A displaced point can move a level across a basis-dimension jump (seen: 200_A2 level 7e-6 above
    a jump at cutoff 0.345), which makes the tracker raise MissingRootError; the next smaller step is then tried, and if all fail
    the errors are NaN and the message says why (recorded in the result, never silent).  Leaves the trackers displaced."""
    fr = list(free)
    last = ""
    for rel in rels:
        try:
            H = np.zeros((len(fr), len(fr)))
            for a, i in enumerate(fr):
                d = np.zeros(4)
                d[i] = rel * errors[i]
                H[:, a] = (fcn.grad(c + d) - fcn.grad(c - d))[fr] / (2.0 * d[i])
        except MissingRootError as e:
            last = f"rel {rel}: {type(e).__name__}: {str(e)[:200]}"
            continue
        out = np.zeros(4)
        out[fr] = np.sqrt(np.diag(2.0 * np.linalg.inv(0.5 * (H + H.T))))
        return out, rel, last
    return np.full(4, np.nan), None, last


def _unconstrained(names_free, c, errors, correlation, accurate) -> list:
    """Reporting flags only (prompts/003 section 2.7 item 4): not pass/fail criteria."""
    reasons = []
    if not accurate:
        reasons.append("hesse_not_accurate")
    for i in names_free:
        if errors[i] >= abs(c[i]):
            reasons.append(f"error_exceeds_value:{PARAM_NAMES[i]}")
    for a in range(len(names_free)):
        for b in range(a + 1, len(names_free)):
            i, j = names_free[a], names_free[b]
            if abs(correlation[i, j]) >= 0.99:
                reasons.append(f"correlation_ge_0.99:{PARAM_NAMES[i]},{PARAM_NAMES[j]}")
    return reasons


def fit_cutoff(levels: LevelSet, trackers: dict, c_start, *, minimizer="lm+migrad", strategy=1, ncall=2000,
               lm_options=None, fixed=(), verify=True, reject_ambiguous_trials=True) -> FitResult:
    """Two-stage fit at one cutoff (prompts/003 section 2.7): LM from `c_start` on the whitened residual with the analytic
    Jacobian (the fixed parameters removed from the LM vector), then iminuit Migrad + Hesse from the LM minimum (`fixed` fixed in
    Migrad).  `trackers` must sit at `c_start` (as built by `make_trackers(..., c_start, ...)`).  Loud: a start vector with fewer
    zeros than lattice levels raises MissingRootError; Migrad invalid after one restart raises RuntimeError; Hesse failure raises
    RuntimeError; tracker errors propagate.  With verify=True the result carries `verify_fit(...)`.

    Stage-1 deviation from "ambiguous stops the fit" (prompts/003 section 7.5), made loud and countable: LM's first trial steps can
    be very large (MINPACK step bound 100 x |x/x_scale|) and land on a parameter vector where a zero sits within ~1 row of a pole,
    so that the scan's order estimate is polluted and the candidate is `ambiguous` (seen: 100_A2, Table B, cutoff 0.315, trial
    c = (4.8e5, -1.25e6, 0, 0): zero 2.8e-6 below a pole).  With `reject_ambiguous_trials` (default) such a trial point is given the
    all-levels-penalised residual (a finite, huge chi^2), so LM shrinks its trust region; every event is listed in
    `ambiguous_trials_lm` (c and the message), counted, and printed as a WARNING.  Only the LM stage does this; Migrad, the start
    vector, the minimum and the verification still raise AmbiguousRootError."""
    import scipy.optimize
    from iminuit import Minuit
    if minimizer != "lm+migrad":
        raise ValueError(f"minimizer {minimizer!r}: only 'lm+migrad' is implemented (Migrad alone is on the do-not-retry list)")
    fixed, free = _free_indices(fixed)
    c_start = tuple(float(x) for x in c_start)
    if len(c_start) != 4:
        raise ValueError("c_start must have 4 components")
    cs = np.asarray(c_start)
    fcn = QCChi2(levels, trackers)
    n_free, ndof = len(free), fcn.m - len(free)
    if ndof <= 0:
        raise ValueError(f"{fcn.m} levels for {n_free} free parameters")
    t_fit = time.perf_counter()
    forms = {trackers[ir].sector.k2_pwave_form for ir in fcn.irreps}
    if len(forms) != 1:
        raise ValueError(f"trackers disagree on the p-wave form: {forms}")
    E_max_tracked = max(trackers[ir].E_max for ir in fcn.irreps)

    chi2_start = fcn(cs)
    if fcn.penalised.any():
        bad = [fcn.labels[i] for i in np.flatnonzero(fcn.penalised)]
        raise MissingRootError(f"start vector {c_start}: fewer QC zeros than lattice levels for {bad}")
    _, J0 = fcn.model(cs)
    sig_gn = _gn_errors(fcn, J0, free)

    # ---- stage 1: Levenberg-Marquardt (MINPACK lmder) on rho = L^T r, C^+ = L L^T
    Lc = np.linalg.cholesky(0.5 * (fcn.Cinv + fcn.Cinv.T))
    fr = list(free)

    def pack(x):
        c = cs.copy()
        c[fr] = x
        return c

    rho_reject = Lc.T @ (fcn.penalty_sigmas * fcn.sigma)       # the all-levels-penalised residual
    ambiguous_trials: list = []

    def rho(x):
        try:
            return Lc.T @ fcn.residual(pack(x))
        except AmbiguousRootError as e:
            if not reject_ambiguous_trials:
                raise
            ambiguous_trials.append({"c": pack(x).tolist(), "message": str(e)})
            return rho_reject

    def jac(x):
        return -Lc.T @ fcn.model(pack(x))[1][:, fr]

    opts = dict(ftol=1e-12, xtol=1e-10, gtol=1e-10, max_nfev=200)
    if lm_options:
        opts.update(lm_options)
    t0 = time.perf_counter()
    lm = scipy.optimize.least_squares(rho, cs[fr], jac=jac, method="lm", x_scale=sig_gn[fr], **opts)
    c_lm = pack(lm.x)
    chi2_lm = fcn(c_lm)
    seconds_lm = time.perf_counter() - t0

    # ---- stage 2: Migrad + Hesse from the LM minimum
    pen0 = fcn.n_penalised_calls
    t1 = time.perf_counter()
    m = Minuit(fcn, c_lm, grad=fcn.grad, name=PARAM_NAMES)      # fcn.errordef == 1 (least squares)
    m.errors = np.where(sig_gn > 0, sig_gn, 1.0)
    m.strategy = int(strategy)
    for n in fixed:
        m.fixed[n] = True
    m.migrad(ncall=int(ncall))
    migrad_calls = 1
    if not m.fmin.is_valid:
        m.migrad(ncall=int(ncall))                              # one restart, recorded
        migrad_calls = 2
    if not m.fmin.is_valid:
        raise RuntimeError(f"Migrad not valid after a restart at cutoff {levels.Ecm_cutoff}: {m.fmin}")
    m.hesse()
    if m.covariance is None:
        raise RuntimeError(f"Hesse returned no covariance at cutoff {levels.Ecm_cutoff}: {m.fmin}")
    seconds_migrad = time.perf_counter() - t1
    n_pen_final = fcn.n_penalised_calls - pen0
    fm = m.fmin
    edm = float(fm.edm)
    accurate = bool(fm.has_accurate_covar)
    c = np.asarray(m.values, dtype=float)
    c[[i for i in range(4) if i not in free]] = cs[[i for i in range(4) if i not in free]]
    errors = np.zeros(4)
    errors[fr] = np.asarray(m.errors, dtype=float)[fr]
    covariance = np.asarray(m.covariance, dtype=float).copy()
    correlation = np.asarray(m.covariance.correlation(), dtype=float).copy()
    seconds = time.perf_counter() - t_fit

    errors_fd, fd_rel, fd_msg = _fd_hessian_errors(fcn, c, errors, free)
    # ---- the model at the minimum (re-positions the trackers at c: Hesse and the FD Hessian end at displaced points)
    chi2 = fcn(c)
    E_model, J_fit = fcn.model(c)
    if fcn.penalised.any():
        raise MissingRootError(f"the minimum {c.tolist()} has fewer QC zeros than lattice levels")
    wrong = [dict(w) for w in fcn.wrong_sign_zeros]
    gn = _gn_errors(fcn, J_fit, free)
    res = FitResult(
        cutoff=float(levels.Ecm_cutoff), k2_pwave_form=forms.pop(), Lbyas=int(levels.Lbyas), labels=[tuple(l) for l in levels.labels],
        E_lat=fcn.E_lat.copy(), sigma=fcn.sigma.copy(), cov_lattice=fcn.cov.copy(), E_model=E_model,
        residual_over_sigma=(fcn.E_lat - E_model) / fcn.sigma, c_start=c_start, c=c, errors=errors, errors_gauss_newton=gn,
        covariance=covariance, correlation=correlation, chi2=chi2, chi2_start=chi2_start, chi2_lm=chi2_lm, ndof=ndof,
        chi2_per_dof=chi2 / ndof, edm=edm, valid=bool(fm.is_valid and np.isfinite(edm) and edm < EDM_MAX), accurate=accurate,
        nfcn=int(m.nfcn), ngrad=int(m.ngrad), n_evaluations=fcn.n_evaluations,
        n_local_rescans=sum(trackers[ir].n_local_rescans for ir in fcn.irreps),
        n_full_rescans=sum(trackers[ir].n_full_rescans for ir in fcn.irreps),
        n_penalised_calls=fcn.n_penalised_calls, n_penalised_calls_final_stage=n_pen_final, seconds=seconds,
        minimizer=minimizer, verification=None, fixed=fixed, n_free=n_free,
        unconstrained=False, unconstrained_reasons=[], wrong_sign_zeros=wrong, n_wrong_sign_zeros=len(wrong),
        n_chi2_calls=fcn.n_calls, n_calls_with_wrong_sign=fcn.n_calls_with_wrong_sign, lm_status=int(lm.status),
        lm_nfev=int(lm.nfev), lm_njev=int(lm.njev or 0), seconds_lm=seconds_lm, seconds_migrad=seconds_migrad,
        E_max_tracked=E_max_tracked, migrad_calls=migrad_calls, ambiguous_trials_lm=ambiguous_trials,
        errors_fd_hessian=errors_fd, fd_hessian_rel=fd_rel, fd_hessian_note=fd_msg)
    res.unconstrained_reasons = _unconstrained(list(free), c, errors, correlation, accurate)
    res.unconstrained = len(res.unconstrained_reasons) > 0
    if verify:
        res.verification = verify_fit(res, trackers, levels)
    return res


def verify_fit(fit: FitResult, trackers: dict, levels: LevelSet, margin: float = VERIFY_MARGIN) -> dict:
    """Fresh label-free scan of every sector at the fitted parameters (the trackers must sit at fit.c).

    Compared per sector: the tracked zeros and the fresh zeros up to E_cmp = min(E_max, E_top + margin), E_top the highest tracked
    root that D3 matches to a lattice level (a zero can enter the tracked range from above E_max, 111_A2; it matters only if it
    lands below the highest matched level, so roots above E_cmp are not compared).  chi^2_fresh = chi^2 of the fresh zeros (all
    those <= E_max) with the same D3 matching.  An ambiguous / even-order candidate <= E_cmp raises AmbiguousRootError (loud).
    Returns {"counts_equal", "max_root_diff" (inf if the counts differ), "chi2_fresh", "chi2_diff" (absolute), "chi2_bound"
    (rigorous bound for root differences of VERIFY_ROOT_TOL), "ok", "n_evaluations", "seconds", "margin", "per_sector":
    {irrep: {"tracked", "fresh", "E_compare_max", "n_fresh_above_compare", "counts_equal"}}, "wrong_sign_zeros_fresh"}."""
    from . import rootscan as _rootscan
    fcn = QCChi2(levels, trackers)
    t0 = time.perf_counter()
    n0 = fcn.n_evaluations
    per_sector, fresh_all, wrong_fresh = {}, {}, []
    counts_equal, max_diff = True, 0.0
    for irrep in fcn.irreps:
        tr = trackers[irrep]
        if not np.array_equal(tr.c, fit.c):
            raise ValueError(f"verify_fit: tracker {irrep} sits at {tr.c.tolist()}, not at the fit {np.asarray(fit.c).tolist()}")
        tracked = tr.energies()
        n_match = min(len(fcn.index[irrep]), len(tracked))
        E_top = tracked[n_match - 1] if n_match else tr.E_min
        E_cmp = min(tr.E_max, float(E_top) + float(margin))
        scan = _rootscan.scan_sector(tr._ev, tr.sector, tuple(tr.c), tr.E_min, tr.E_max, tr.scan_options)
        tr.n_verify_rescans += 1
        bad = [c for c in scan.candidates if c.kind in ("ambiguous", "even_zero") and c.E <= E_cmp]
        if bad:
            raise AmbiguousRootError(f"verify_fit {irrep}: {len(bad)} ambiguous/even-order candidate(s) <= {E_cmp:.6f}, first at "
                                     f"E={bad[0].E!r} (kind {bad[0].kind}); a physics decision, not dropped silently")
        fz = scan.zeros()
        fresh_all[irrep] = fz
        tk, fr = tracked[tracked <= E_cmp], fz[fz <= E_cmp]
        eq = len(tk) == len(fr)
        counts_equal &= eq
        if eq and len(tk):
            max_diff = max(max_diff, float(np.max(np.abs(tk - fr))))
        per_sector[irrep] = {"tracked": tracked.tolist(), "fresh": fz.tolist(), "E_compare_max": E_cmp,
                             "n_fresh_above_compare": int(np.sum(fz > E_cmp)), "counts_equal": bool(eq)}
        wrong_fresh += [{"irrep": irrep, "E": float(c.E), "slope": float(c.slope)} for c in scan.of_kind("zero")
                        if c.slope is not None and c.slope < 0.0]
    chi2_fresh, _ = fcn.chi2_from_roots(fresh_all)
    q = fcn.m * VERIFY_ROOT_TOL ** 2 / float(np.linalg.eigvalsh(fcn.cov).min())
    bound = 2.0 * np.sqrt(max(fit.chi2, 0.0)) * np.sqrt(q) + q
    diff = abs(chi2_fresh - fit.chi2)
    return {"counts_equal": bool(counts_equal), "max_root_diff": max_diff if counts_equal else float("inf"),
            "chi2_fresh": chi2_fresh, "chi2_diff": diff, "chi2_bound": float(bound),
            "ok": bool(counts_equal and max_diff <= VERIFY_ROOT_TOL and diff <= bound),
            "n_evaluations": fcn.n_evaluations - n0, "seconds": time.perf_counter() - t0, "margin": float(margin),
            "per_sector": per_sector, "wrong_sign_zeros_fresh": wrong_fresh}


def _log(msg):
    print(msg, file=sys.stderr, flush=True)


def cutoff_scan(data_dir, Lbyas, cutoffs, c_start, k2_pwave_form, irreps=IRREPS, E_min=GRID_E0, E_max=0.36, out_json=None,
                verify=True, fixed=(), evaluate=_fast.evaluate_point) -> list:
    """One `fit_cutoff` per cutoff, every one from `c_start` (not warm-started: prompts/003 section 2.7 item 4 and section 8).
    Fresh trackers per cutoff over [E_min, min(E_max, cutoff + TRACK_MARGIN)] for the sectors that have levels at that cutoff.
    Writes `save_results` if `out_json` is given."""
    fits = []
    for cutoff in cutoffs:
        lv = _lattice.load_levels(data_dir, Lbyas, list(irreps), float(cutoff))
        used = [ir for ir in irreps if any(l[0] == ir for l in lv.labels)]
        e_hi = min(float(E_max), float(cutoff) + TRACK_MARGIN)
        t0 = time.perf_counter()
        trackers = make_trackers(used, c_start, Lbyas=Lbyas, k2_pwave_form=k2_pwave_form, E_min=E_min, E_max=e_hi, evaluate=evaluate)
        t_scan = time.perf_counter() - t0
        _log(f"[cutoff {cutoff}] {len(lv.labels)} levels in {used}; scans to E_max={e_hi:.4f} took {t_scan:.1f} s")
        r = fit_cutoff(lv, trackers, c_start, fixed=fixed, verify=verify)
        r.seconds_scan = t_scan
        fits.append(r)
        _log(f"[cutoff {cutoff}] chi2={r.chi2:.4f} ndof={r.ndof} valid={r.valid} edm={r.edm:.2e} fit {r.seconds:.1f} s"
             + (f" verify ok={r.verification['ok']} {r.verification['seconds']:.1f} s" if r.verification else ""))
        if out_json is not None:
            save_results(fits, out_json)           # after every cutoff, so that a long run keeps what it has
    return fits


# ------------------------------------------------------------------------------------------------ output
def _jsonable(x):
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
    if isinstance(x, np.bool_):
        return bool(x)
    return x


def _document(fits: list, extra: dict | None) -> dict:
    import os
    import iminuit
    from . import bench
    from .constants import GRID_N, PRODUCTION
    prov = bench.provenance()
    prov["iminuit"] = iminuit.__version__
    prov["machine"] = bench.machine_info()
    prov["command"] = " ".join(sys.argv)
    prov["threads_env"] = {k: os.environ.get(k) for k in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS")}
    f0 = fits[0] if fits else None
    doc = {
        "kind": "qc3.fit", "is_a_fit": True, "provenance": prov,
        "constants": {
            "production": PRODUCTION, "k2_pwave_form": f0.k2_pwave_form if f0 else None,
            "c_start": list(f0.c_start) if f0 else None, "fixed": list(f0.fixed) if f0 else None,
            "covariance": "D1: .jack files = one value per gauge configuration (N=603 at L20); jackknife covariance of the mean "
                          "(qc3.lattice.load_levels), eigenvalue pseudo-inverse (eig_rel_cut 1e-12)",
            "matching": "D3: sort both lists ascending, the lowest m QC levels for the m lattice levels (n >= m allowed); "
                        "n < m: finite penalty 100 sigma per missing level, counted",
            "exclusion": "D8: exact free-level solutions of the QC are excluded (qc3.rootscan kind 'free_level')",
            "residue_sign": "D9 flag only: zeros with u^T Q' u < 0 stay in the matching and are listed (wrong_sign_zeros)",
            "grid": {"E0": GRID_E0, "dE": GRID_DE, "N": GRID_N},
            "tracking": {"track_margin": TRACK_MARGIN, "verify_margin": VERIFY_MARGIN, "edm_max": EDM_MAX,
                         "verify_root_tol": VERIFY_ROOT_TOL},
        },
        "seeds": {"stochastic_inputs": None, "note": "deterministic: no random numbers; start vectors are inputs"},
        "fits": [_jsonable({f.name: getattr(fit, f.name) for f in fields(fit)}) for fit in fits],
    }
    if extra:
        doc.update(_jsonable(extra))
    return doc


def save_results(fits: list, path, extra: dict | None = None) -> None:
    """JSON: kind, is_a_fit, provenance (commit, dirty flag, core commit/version, numpy/scipy/iminuit, date, machine, command),
    constants, seeds, fits (every FitResult field), plus `extra` keys at the top level."""
    import pathlib
    p = pathlib.Path(path)
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text(json.dumps(_document(fits, extra), indent=1) + "\n")


def format_fit(r: FitResult) -> str:
    """Plain-text summary of one fit (chi2, chi2/nu, values +- Hesse errors, correlation, flags, D9 warnings)."""
    free = [i for i, n in enumerate(PARAM_NAMES) if n not in r.fixed]
    lines = [f"cutoff {r.cutoff:.3f}  form {r.k2_pwave_form}  m={len(r.labels)} levels  n_free={r.n_free}  ndof={r.ndof}  "
             f"chi2={r.chi2:.4f}  chi2/nu={r.chi2_per_dof:.4f}  (start {r.chi2_start:.3f}, LM {r.chi2_lm:.4f})  valid={r.valid}  "
             f"EDM={r.edm:.2e}  accurate_cov={r.accurate}"]
    for i in range(4):
        tag = ("fixed" if PARAM_NAMES[i] in r.fixed else
               f"+- {r.errors[i]:.4g} (GN {r.errors_gauss_newton[i]:.4g}, FD-Hessian {r.errors_fd_hessian[i]:.4g})")
        lines.append(f"    {PARAM_NAMES[i]:6s} = {r.c[i]: .6e} {tag}")
    if len(free) > 1:
        lines.append("    correlation (free): " + "; ".join(f"{PARAM_NAMES[i]},{PARAM_NAMES[j]}={r.correlation[i, j]:+.3f}"
                                                            for a, i in enumerate(free) for j in free[a + 1:]))
    lines.append(f"    unconstrained={r.unconstrained} {r.unconstrained_reasons}")
    lines.append(f"    max|r/sigma|={np.max(np.abs(r.residual_over_sigma)):.2f}  nfcn={r.nfcn} ngrad={r.ngrad} chi2_calls={r.n_chi2_calls} "
                 f"QC_evals={r.n_evaluations} local/full rescans={r.n_local_rescans}/{r.n_full_rescans} "
                 f"penalised(final)={r.n_penalised_calls_final_stage}  time scan/LM/Migrad+Hesse="
                 f"{r.seconds_scan:.1f}/{r.seconds_lm:.1f}/{r.seconds_migrad:.1f} s")
    if r.verification:
        v = r.verification
        lines.append(f"    verification: ok={v['ok']} counts_equal={v['counts_equal']} max_root_diff={v['max_root_diff']:.2e} "
                     f"|chi2_fresh-chi2|={v['chi2_diff']:.2e} (bound {v['chi2_bound']:.2e}) {v['seconds']:.1f} s")
    if r.fd_hessian_note:
        lines.append(f"    NOTE: FD-Hessian cross-check step 0.02 sigma failed ({r.fd_hessian_note[:140]}); used rel step {r.fd_hessian_rel}")
    for a in r.ambiguous_trials_lm:
        lines.append(f"    WARNING: LM trial point rejected (ambiguous candidate) c={[float(f'{x:.6g}') for x in a['c']]}: {a['message'][:160]}")
    for w in r.wrong_sign_zeros:
        lines.append(f"    WARNING: wrong-sign zero {w['irrep']} E={w['E']:.9f} slope={w['slope']:.3e} matched={w['matched']}")
    if r.verification:
        for w in r.verification["wrong_sign_zeros_fresh"]:
            lines.append(f"    WARNING: wrong-sign zero (fresh scan) {w['irrep']} E={w['E']:.9f} slope={w['slope']:.3e}")
    return "\n".join(lines)


# ------------------------------------------------------------------------------------------------ code equivalence
def _read_legacy_fit_levels():
    from .bench import ROOT
    path = ROOT / "results" / "legacy_L20_passonly_hotwindow" / "Ecm_cut_0p335" / "l20_full_cut_0p335_fit_levels_L20.dat"
    rows = []
    with open(path) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            p = line.split()
            rows.append((p[3], int(p[5]), float(p[6]), float(p[8])))
    return str(path.relative_to(ROOT)), rows


def code_equivalence(data_dir, out_json=None, Lbyas=20, cutoff=0.335, irreps=IRREPS, E_max=0.36, evaluate=_fast.evaluate_point,
                     verify=True) -> dict:
    """prompts/003 section 2.7 item 4 (Digonto 2026-10-24: a code-equivalence check, not a reproduction of v32w's 11.4897).
    coded_v2, 12 levels at 0.335, D1 covariance: (i) our model levels and corrected chi^2 at C_V32W; (ii) the refit from C_V32W;
    (iii) our zeros at C_LEGACY_0335 against the 12 legacy `model_Ecm` of the legacy hot-window fit (legacy chi^2 / errors invalid)."""
    lv = _lattice.load_levels(data_dir, Lbyas, list(irreps), float(cutoff))
    out = {"form": "coded_v2", "cutoff": float(cutoff), "labels": [list(l) for l in lv.labels], "E_lat": lv.Ecm.tolist()}
    # (i) C_V32W
    tr = make_trackers(irreps, C_V32W, Lbyas=Lbyas, k2_pwave_form="coded_v2", E_min=GRID_E0, E_max=E_max, evaluate=evaluate)
    f = QCChi2(lv, tr)
    chi2 = f(np.asarray(C_V32W))
    E_model, _ = f.model(np.asarray(C_V32W))
    out["c_v32w"] = {"c": list(C_V32W), "chi2_corrected_D1": chi2, "ndof_if_4_free": len(lv.labels) - 4,
                     "v32w_reported_chi2_not_reproduced_not_compared": V32W_CHI2, "n_penalised_calls": f.n_penalised_calls,
                     "E_model": E_model.tolist(), "residual_over_sigma": ((lv.Ecm - E_model) / f.sigma).tolist(),
                     "tracked_zeros_by_sector": {ir: tr[ir].energies().tolist() for ir in f.irreps},
                     "n_qc_zeros_vs_levels": {ir: [len(tr[ir].energies()), len(f.index[ir])] for ir in f.irreps},
                     "wrong_sign_zeros": f.wrong_sign_zeros}
    # (iii) C_LEGACY_0335
    trl = make_trackers(irreps, C_LEGACY_0335, Lbyas=Lbyas, k2_pwave_form="coded_v2", E_min=GRID_E0, E_max=E_max, evaluate=evaluate)
    fl = QCChi2(lv, trl)
    chi2_l = fl(np.asarray(C_LEGACY_0335))
    E_model_l, _ = fl.model(np.asarray(C_LEGACY_0335))
    path, legacy = _read_legacy_fit_levels()
    rows = []
    for i, (ir, n, E_lat_leg, E_leg) in enumerate(legacy):
        zs, poles = trl[ir].energies(), np.asarray(trl[ir].scan.poles())
        near = float(zs[np.argmin(np.abs(zs - E_leg))])
        d_pole = float(np.min(np.abs(poles - E_leg))) if len(poles) else float("inf")
        rows.append({"irrep": ir, "n": n, "legacy_model_Ecm": E_leg, "our_nearest_zero": near, "diff": near - E_leg,
                     "model_D3": float(E_model_l[i]), "nearest_pole_distance": d_pole,
                     "nearest_pole_distance_rows": d_pole / GRID_DE, "legacy_lattice_Ecm": E_lat_leg,
                     "labels_match": [ir, n] == list(lv.labels[i])})
    out["c_legacy_0335"] = {"c": list(C_LEGACY_0335), "legacy_table": path, "chi2_corrected_D1": chi2_l,
                            "n_penalised_calls": fl.n_penalised_calls, "rows": rows,
                            "max_abs_diff_all": max(abs(r["diff"]) for r in rows),
                            "max_abs_diff_no_pole_within_10_rows": max([abs(r["diff"]) for r in rows if r["nearest_pole_distance_rows"] > 10] or [0.0]),
                            "note": "legacy model levels are linear interpolations of det Q sign changes on the 1-row grid; the "
                                    "legacy chi^2 and errors are invalid (OPEN_ISSUES #1) and are not used"}
    # (ii) refit from C_V32W (coded_v2, 4 free)
    fits = cutoff_scan(data_dir, Lbyas, (cutoff,), C_V32W, "coded_v2", irreps=irreps, E_max=E_max, out_json=None, verify=verify,
                       evaluate=evaluate)
    out["refit_from_c_v32w"] = "fits[0] of this file"
    if out_json is not None:
        save_results(fits, out_json, extra={"code_equivalence": out})
    return {"code_equivalence": out, "fits": fits}


# ------------------------------------------------------------------------------------------------ CLI
def build_parser() -> argparse.ArgumentParser:
    ap = argparse.ArgumentParser(prog="python -m qc3.fit", description="L20 two-stage K_df,3 fit (prompts/003 section 2.7). "
                                 "Run single-threaded: OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1.")
    ap.add_argument("--cutoffs", type=float, nargs="+", default=[0.315, 0.325, 0.335, 0.345])
    ap.add_argument("--form", default=K2_PWAVE_FORM_DEFAULT, choices=("q0_consistent", "coded_v2"))
    ap.add_argument("--out", default=None, help="result JSON")
    ap.add_argument("--fix", nargs="+", default=[], metavar="NAME", help=f"parameters to fix, from {PARAM_NAMES}")
    ap.add_argument("--c-start", nargs="+", default=["baseline"], help="baseline | v32w | legacy0335 | four floats")
    ap.add_argument("--no-verify", action="store_true")
    ap.add_argument("--code-equivalence", action="store_true",
                    help="coded_v2 code-equivalence record at 0.335 (c_v32w and c_legacy0335); ignores --cutoffs/--form/--fix")
    ap.add_argument("--data-dir", default=None)
    ap.add_argument("--Lbyas", type=int, default=20)
    ap.add_argument("--e-max", type=float, default=0.36, help="upper cap of the tracked range (default 0.36; trackers use cutoff + margin)")
    return ap


def _start_vector(tokens, fixed) -> tuple:
    presets = {"baseline": C_BASELINE, "v32w": C_V32W, "legacy0335": C_LEGACY_0335}
    if len(tokens) == 1 and tokens[0] in presets:
        c = list(presets[tokens[0]])
        for n in fixed:
            c[PARAM_NAMES.index(n)] = 0.0         # --fix: fixed parameters start (and stay) at 0 unless given explicitly
        return tuple(c)
    if len(tokens) == 4:
        return tuple(float(t) for t in tokens)
    raise SystemExit(f"--c-start: expected baseline|v32w|legacy0335 or four floats, got {tokens}")


def main(argv=None) -> int:
    """CLI: python -m qc3.fit --cutoffs 0.315 0.325 0.335 0.345 --form q0_consistent --out results/fit_L20_..._<tag>.json
    [--fix KB KE] [--c-start baseline|v32w|legacy0335|<4 floats>] [--code-equivalence]"""
    args = build_parser().parse_args(argv)
    data_dir = args.data_dir or _lattice.default_data_dir()
    if args.code_equivalence:
        res = code_equivalence(data_dir, args.out, Lbyas=args.Lbyas, E_max=args.e_max, verify=not args.no_verify)
        ce = res["code_equivalence"]
        print(f"code equivalence (coded_v2, cutoff {ce['cutoff']}, D1 covariance)")
        print(f"  chi2 at c_v32w = {ce['c_v32w']['chi2_corrected_D1']:.4f} (v32w quoted {V32W_CHI2:.4f}: different inputs, not compared)")
        cl = ce["c_legacy_0335"]
        print(f"  chi2 at c_legacy0335 = {cl['chi2_corrected_D1']:.4f}; max |our zero - legacy model level| = "
              f"{cl['max_abs_diff_all']:.3e} (no pole within 10 rows: {cl['max_abs_diff_no_pole_within_10_rows']:.3e})")
        print(format_fit(res["fits"][0]))
        return 0
    fixed, _ = _free_indices(args.fix) if args.fix else ((), None)
    c_start = _start_vector(args.c_start, fixed)
    fits = cutoff_scan(data_dir, args.Lbyas, args.cutoffs, c_start, args.form, E_max=args.e_max, out_json=args.out,
                       verify=not args.no_verify, fixed=fixed)
    print(f"c_start={c_start} fixed={fixed} form={args.form}")
    for r in fits:
        print(format_fit(r))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

"""Root tracking across fit-parameter changes (prompts/003 section 2.5; docs/METHODS.md section 13).

With Q(E; c) = F~(E) + sum_i c_i B~_i(E) (exact; B~_i = V^dagger B_i V / (L xi)^6 from PointResult.K3_basis and .Vsel) and u the
unit null vector of the real symmetric Q(E0), the Hellmann-Feynman theorem and the implicit function theorem give

        J_i = dE0/dc_i = - (u^T B~_i u) / (u^T Q'(E0) u),      Q' = [Q(E0+h) - Q(E0-h)] / 2h,  u from Qbar = [Q(E0+h) + Q(E0-h)] / 2.   (1)

`root_jacobian` evaluates (1) from the two points E0 +- h (the same two points give the Delta_h pole/zero re-classification
and the slope u^T Q' u).  `secant_root` is the corrector (a secant on g(E) = sign * exp(l(E) - l_ref), which is proportional to
det Q and so linear at a simple zero) and never evaluates at the returned energy (an exactly zero pivot in the LU is
`SingularQError`).  `SectorRoots` is the predictor-corrector continuation of all zeros of one sector with the safeguards, the
local and the full rescan fallbacks of section 2.5.

Loud events (never a silent drop): a candidate of kind "ambiguous" or an even-order zero away from a free level raises
`AmbiguousRootError`; a full scan that disagrees with the tracked set raises `RootSetChangedError`.  Fallbacks that are part of
the method (secant not converged, landed on a pole, dimension jump on the path, predictor/safeguard violation) are counted in
`n_local_rescans` / `n_full_rescans`.  D8: free-level solutions are kept apart in `free_level_solutions`.  D9: the sign of
`slope` is recorded, `wrong_sign_roots()` lists the tracked zeros with slope < 0, nothing is ever excluded.
"""
from __future__ import annotations

import bisect
from dataclasses import dataclass

import numpy as np

from . import fast as _fast
from . import rootscan as _rootscan
from . import roots as _roots
from .constants import GRID_DE
from .point import MissingRootError, SingularQError
from .rootscan import RootCandidate, ScanResult  # noqa: F401
from .roots import AmbiguousRootError, DimensionJumpError, PoleError  # noqa: F401
from .sector import Sector

TRACK_DEFAULTS = dict(h=1e-6, D=4.0, xtol=1e-11, max_iter=12, max_move_rows=20.0, safeguard_factor=10.0, safeguard_abs=1e-6,
                      local_rows=60, local_steps_rows=(4.0, 1.0, 0.25), park_margin=0.005)
_DISTINCT = 1e-9                 # two tracked roots closer than this are "merged/crossed"
_IMAG_TOL = 1e-10                # relative bound on Im Q, Im B~ (Q is real symmetric)
_KEYS = ("iso0", "iso1", "B", "E")


class RootSetChangedError(MissingRootError):
    """A full scan at the current c disagrees with the tracked root set (count, or a root differs by > 1e-9)."""


class SecantNotConvergedError(MissingRootError):
    """The secant corrector did not converge within max_iter (or met a flat / non-finite g); the caller falls back to a rescan."""


@dataclass
class TrackedRoot:
    E: float
    J: np.ndarray              # dE0/dc_i, shape (4,)
    slope: float               # u^T Q'(E0) u  (D9: its sign is recorded, > 0 physical)
    discriminant: float        # Delta_h(E0) = l(E0) - [l(E0-h) + l(E0+h)]/2
    N: int
    segment: tuple             # (E_lo, E_hi) between basis-dimension jumps (root_jacobian alone: (E0-h, E0+h), the checked window)


def _real(A, what):
    A = np.asarray(A)
    if np.iscomplexobj(A):
        scale = float(np.max(np.abs(A))) or 1.0
        im = float(np.max(np.abs(A.imag)))
        if im > _IMAG_TOL * scale:
            raise ValueError(f"root_jacobian: {what} is not real (max|Im| = {im:.3e}, max|.| = {scale:.3e})")
        A = A.real
    return np.asarray(A, dtype=float)


def root_jacobian(evaluate, sector: Sector, c, E0: float, h=1e-6, D=4.0, lnabsdet_at_root: float | None = None) -> TrackedRoot:
    """Eq. (1): J_i = -u^T B_i u / u^T Q'(E0) u; 2 evaluations at E0 +- h (+1 at E0 if lnabsdet_at_root is None).

    Also returns slope = u^T Q' u, the discriminant Delta_h and N.  Raises PoleError / AmbiguousRootError when Delta_h > -D
    (the point is not a zero), DimensionJumpError if the basis dimension differs among E0 - h, E0, E0 + h.  If the evaluation
    at E0 itself raises SingularQError, E0 is an exact root: l(E0) = -inf (documented case, not a silent catch)."""
    E0, h, c = float(E0), float(h), tuple(float(x) for x in c)
    N0 = None
    if lnabsdet_at_root is None:
        try:
            r0 = evaluate(sector, E0, c, want=("Q",))
            l0, N0 = float(r0.lnabsdet), int(r0.N)
        except SingularQError:
            l0 = -np.inf
    else:
        l0 = float(lnabsdet_at_root)
    want = ("Q", "Vsel", "K3")
    pm = evaluate(sector, E0 - h, c, want=want)
    pp = evaluate(sector, E0 + h, c, want=want)
    Ns = {"E0-h": int(pm.N), "E0+h": int(pp.N)}
    if N0 is not None:
        Ns["E0"] = N0
    if len(set(Ns.values())) != 1:
        raise DimensionJumpError(f"basis dimension changes near E0={E0!r}: {Ns}")
    if np.isnan(l0) or not (np.isfinite(pm.lnabsdet) and np.isfinite(pp.lnabsdet)):
        raise ValueError(f"root_jacobian: non-finite ln|det Q| near E0={E0!r}: {pm.lnabsdet}, {l0}, {pp.lnabsdet}")
    disc = l0 - 0.5 * (float(pm.lnabsdet) + float(pp.lnabsdet))
    if disc > -D:
        info = f"E0={E0!r} is not a zero (Delta_h = {disc:.3f}, need <= {-D})"
        raise (PoleError if disc >= D else AmbiguousRootError)(info)
    Vm, Vp = np.asarray(pm.Vsel), np.asarray(pp.Vsel)
    if Vm.shape != Vp.shape or not np.array_equal(Vm, Vp):
        raise ValueError(f"root_jacobian: the projector V_sel differs between E0-h and E0+h at E0={E0!r}")
    Qm, Qp = _real(pm.Q, "Q(E0-h)"), _real(pp.Q, "Q(E0+h)")
    Qbar = 0.5 * (Qp + Qm)
    w, v = np.linalg.eigh(0.5 * (Qbar + Qbar.T))
    u = v[:, int(np.argmin(np.abs(w)))]
    denom = float(u @ ((Qp - Qm) / (2.0 * h)) @ u)
    L6 = float(sector.L) ** 6
    Vh = Vp.conj().T
    J = np.empty(4)
    for i, k in enumerate(_KEYS):
        # B~_i(E0) = average of the two points (error O(h^2)); same V at both points (checked above)
        Bi = 0.5 * (_real(Vh @ np.asarray(pm.K3_basis[k]) @ Vp, f"B~_{k}(E0-h)") + _real(Vh @ np.asarray(pp.K3_basis[k]) @ Vp, f"B~_{k}(E0+h)")) / L6
        J[i] = -float(u @ Bi @ u) / denom
    return TrackedRoot(E=E0, J=J, slope=denom, discriminant=float(disc), N=int(pp.N), segment=(E0 - h, E0 + h))


def secant_root(evaluate, sector: Sector, c, E_pred: float, scale: float, xtol=1e-11, max_iter=12) -> tuple:
    """Secant corrector.  g(E) = sign * exp(l(E) - l_ref) (proportional to det Q) is iterated from E_pred and
    E_pred + max(1e-8, 0.1|scale|) until |E_{k+1} - E_k| < xtol.  Returns (E, lnabsdet of the last evaluated iterate,
    n_evaluations); the returned E is never evaluated.  SingularQError at an iterate (exactly zero pivot) returns that
    iterate with ln|det| = -inf.  SecantNotConvergedError if not converged in max_iter secant steps or if g is flat/non-finite."""
    c = tuple(float(x) for x in c)
    state = {"n": 0, "lref": None}

    def g(E):
        state["n"] += 1
        r = evaluate(sector, float(E), c, want=("Q",))
        s, l = int(r.sign), float(r.lnabsdet)
        if state["lref"] is None:
            if not np.isfinite(l):
                raise SecantNotConvergedError(f"secant: non-finite ln|det Q| = {l} at the starting point E={E!r}")
            state["lref"] = l
        with np.errstate(over="ignore"):
            val = s * float(np.exp(l - state["lref"]))
        if not np.isfinite(val):
            raise SecantNotConvergedError(f"secant: g(E) is not finite at E={E!r} (ln|det Q| = {l}, l_ref = {state['lref']})")
        return val, l

    eps = max(1e-8, 0.1 * abs(float(scale)))
    Ea, Eb = float(E_pred), float(E_pred) + eps
    try:
        ga, la = g(Ea)
        gb, lb = g(Eb)
    except SingularQError:
        return (Ea if state["n"] <= 1 else Eb), -np.inf, state["n"]
    for _ in range(int(max_iter)):
        if gb == ga:
            raise SecantNotConvergedError(f"secant: flat g (g(Ea) == g(Eb)) at E={Eb!r}")
        Ec = Eb - gb * (Eb - Ea) / (gb - ga)
        if not np.isfinite(Ec):
            raise SecantNotConvergedError(f"secant: non-finite iterate from E={Eb!r}")
        if abs(Ec - Eb) < xtol:
            return float(Ec), lb, state["n"]
        try:
            gc, lc = g(Ec)
        except SingularQError:
            return float(Ec), -np.inf, state["n"]
        Ea, ga, Eb, gb, lb = Eb, gb, Ec, gc, lc
    raise SecantNotConvergedError(f"secant: not converged in {max_iter} steps (last iterates {Ea!r}, {Eb!r})")


class _CountingEvaluate:
    """evaluate(sector, E, c, want) with a call counter; forwards an evaluator's own `dimension_jumps` attribute (rootscan uses it)."""

    def __init__(self, evaluate):
        self.evaluate, self.n = evaluate, 0
        dj = getattr(evaluate, "dimension_jumps", None)
        if dj is not None:
            self.dimension_jumps = dj

    def __call__(self, sector, E, c, want=("Q",)):
        self.n += 1
        return self.evaluate(sector, E, c, want=want)


class SectorRoots:
    """Tracked zeros of one sector.  Attributes: sector, c, E_min, E_max, scan (the initial scan, kept), roots (all tracked zeros
    including parked ones, sorted), free_level_solutions (D8, kept apart, never tracked or matched), n_evaluations,
    n_local_rescans, n_full_rescans (continuation fallbacks), n_verify_rescans (explicit `rescan()` calls), n_substeps."""

    def __init__(self, sector: Sector, c0, E_min: float, E_max: float, evaluate=_fast.evaluate_point,
                 scan_options: dict | None = None, track_options: dict | None = None):
        self.track = dict(TRACK_DEFAULTS)
        if track_options:
            unknown = set(track_options) - set(TRACK_DEFAULTS)
            if unknown:
                raise ValueError(f"SectorRoots: unknown track_options {sorted(unknown)}; allowed: {sorted(TRACK_DEFAULTS)}")
            self.track.update(track_options)
        t = self.track
        _roots._check_xtol(t["xtol"], t["h"], t["D"])
        self.sector = sector
        self.c = np.asarray([float(x) for x in c0], dtype=float)
        if self.c.shape != (4,):
            raise ValueError("c0 must have 4 components (Kiso0, Kiso1, KB, KE)")
        self.E_min, self.E_max = float(E_min), float(E_max)
        self.scan_options = dict(scan_options) if scan_options else None
        self._ev = _CountingEvaluate(evaluate)
        self.n_local_rescans = self.n_full_rescans = self.n_verify_rescans = self.n_substeps = 0
        self.roots: list = []
        self.free_level_solutions: list = []
        self._poles: list = []
        self._load_full_scan(tuple(self.c))

    # ---------------------------------------------------------------- bookkeeping
    @property
    def n_evaluations(self) -> int:
        return self._ev.n

    def _segment_of(self, E):
        edges = [self.E_min] + [float(x) for x in self.scan.jumps] + [self.E_max]
        k = min(max(bisect.bisect_right(edges, E) - 1, 0), len(edges) - 2)
        return (edges[k], edges[k + 1])

    def _check_loud(self, scan: ScanResult, where: str):
        bad = scan.of_kind("ambiguous") + scan.of_kind("even_zero")
        if bad:
            b = bad[0]
            raise AmbiguousRootError(f"{where}: {len(bad)} candidate(s) of kind ambiguous/even_zero in sector {self.sector.irrep}, first at "
                                     f"E={b.E!r} (kind {b.kind}, order {b.order}, Delta_h {b.discriminant:.3f}); a physics decision, "
                                     f"not dropped silently (prompts/003 section 2.6)")

    def _make_root(self, E, c, lnabsdet_at_root=None):
        t = self.track
        tr = root_jacobian(self._ev, self.sector, c, E, t["h"], t["D"], lnabsdet_at_root)
        tr.segment = self._segment_of(tr.E)
        return tr

    def _load_full_scan(self, c):
        """Full scan at c; the tracked set is replaced by its zeros (also the constructor's initial scan)."""
        scan = _rootscan.scan_sector(self._ev, self.sector, c, self.E_min, self.E_max, self.scan_options)
        self._check_loud(scan, "full scan")
        self.scan = scan                       # its jumps/segments are used by _segment_of and the safeguards
        self.roots = sorted((self._make_root(z, c) for z in scan.zeros()), key=lambda r: r.E)
        self.free_level_solutions = list(scan.of_kind("free_level"))
        self._poles = [float(p) for p in scan.poles()]
        return scan

    # ---------------------------------------------------------------- public
    def energies(self) -> np.ndarray:
        """Tracked zeros with E <= E_max, sorted."""
        return np.array([r.E for r in self.roots if r.E <= self.E_max], dtype=float)

    def jacobian(self) -> np.ndarray:
        """(len(energies()), 4): dE/dc_i of the tracked zeros with E <= E_max, same order."""
        rows = [r.J for r in self.roots if r.E <= self.E_max]
        return np.array(rows, dtype=float).reshape(len(rows), 4)

    def wrong_sign_roots(self) -> list:
        """D9 flag only: tracked zeros with E <= E_max and slope < 0.  They stay in energies()."""
        return [r for r in self.roots if r.E <= self.E_max and r.slope < 0.0]

    def rescan(self) -> ScanResult:
        """Full scan at the current c (explicit verification; the tracked set is unchanged).  AmbiguousRootError for an ambiguous/
        even-order candidate; RootSetChangedError if its zeros() differ from the tracked set in count or by > 1e-9."""
        before = self._ev.n
        scan = _rootscan.scan_sector(self._ev, self.sector, tuple(self.c), self.E_min, self.E_max, self.scan_options)
        self.n_verify_rescans += 1
        self._check_loud(scan, "rescan")
        fresh, tracked = scan.zeros(), self.energies()
        if len(fresh) != len(tracked):
            raise RootSetChangedError(f"rescan {self.sector.irrep}: {len(fresh)} zeros in a fresh scan vs {len(tracked)} tracked "
                                      f"(fresh {fresh.tolist()}, tracked {tracked.tolist()}); {self._ev.n - before} evaluations")
        if len(fresh) and np.max(np.abs(fresh - tracked)) > _DISTINCT:
            raise RootSetChangedError(f"rescan {self.sector.irrep}: max |fresh - tracked| = {np.max(np.abs(fresh - tracked)):.3e} > {_DISTINCT}")
        return scan

    def update(self, c) -> np.ndarray:
        """Move the tracked zeros to the parameters c (section 2.5): sub-steps of at most max_move_rows grid rows, tangent
        predictor, secant corrector, Delta_h re-check, safeguards, local then full rescan fallbacks.  Returns energies()."""
        c_new = np.asarray([float(x) for x in c], dtype=float)
        if c_new.shape != (4,):
            raise ValueError("c must have 4 components (Kiso0, Kiso1, KB, KE)")
        c_old = self.c.copy()
        dc = c_new - c_old
        if not np.any(dc != 0.0):
            return self.energies()
        t = self.track
        active = [r for r in self.roots if r.E <= self.E_max + t["park_margin"]]
        move = max((abs(float(r.J @ dc)) for r in active), default=0.0)
        n_sub = max(1, int(np.ceil(move / (t["max_move_rows"] * GRID_DE) - 1e-12)))
        for s in range(1, n_sub + 1):
            c_s = c_new if s == n_sub else c_old + (s / n_sub) * dc
            self._substep(c_s)
            self.n_substeps += 1
        return self.energies()

    # ---------------------------------------------------------------- continuation internals
    def _substep(self, c_s):
        t = self.track
        dc_s = c_s - self.c
        cs = tuple(float(x) for x in c_s)
        new: list = []
        for k, r in enumerate(self.roots):
            scale = float(r.J @ dc_s)
            E_pred = r.E + scale
            if E_pred > self.E_max + t["park_margin"]:
                # parked: the tangent predictor moves it, no evaluation (re-tracked when it comes back below E_max + margin)
                new.append(TrackedRoot(E=E_pred, J=r.J.copy(), slope=r.slope, discriminant=r.discriminant, N=r.N, segment=r.segment))
                continue
            tr = self._correct(cs, r, E_pred, scale, new)
            if tr is None:
                tr = self._local_rescan(cs, E_pred, scale, new)
            if tr is None:
                self.n_full_rescans += 1
                self._load_full_scan(cs)
                self.c = c_s.copy()
                return
            new.append(tr)
        new.sort(key=lambda x: x.E)
        for a, b in zip(new[:-1], new[1:]):
            if b.E - a.E < _DISTINCT:     # cannot happen after the per-root distinctness check; loud if it ever does
                raise RootSetChangedError(f"tracked roots merged: {a.E!r}, {b.E!r} (sector {self.sector.irrep})")
        self.roots = new
        self.c = c_s.copy()

    def _distinct(self, E, taken):
        return all(abs(E - x.E) >= _DISTINCT for x in taken)

    def _correct(self, cs, r, E_pred, scale, taken):
        """Secant corrector + Delta_h re-check + safeguards; None -> local rescan.  AmbiguousRootError propagates (loud)."""
        t = self.track
        try:
            E1, l_last, _ = secant_root(self._ev, self.sector, cs, E_pred, scale, t["xtol"], t["max_iter"])
        except SecantNotConvergedError:
            return None
        if abs(E1 - E_pred) > t["safeguard_factor"] * abs(scale) + t["safeguard_abs"]:
            return None
        lo, hi = min(r.E, E_pred, E1), max(r.E, E_pred, E1)
        if any(lo < p < hi for p in self._poles):
            return None
        if any(lo < j < hi for j in self.scan.jumps):
            return None
        if not self._distinct(E1, taken):
            return None
        try:
            tr = self._make_root(E1, cs, lnabsdet_at_root=l_last)
        except (PoleError, AmbiguousRootError, DimensionJumpError):
            # the corrector's iterate is not a zero (Delta_h > -D): a false secant convergence next to a pole (111_A2 q0_consistent,
            # reports 2026-10-24).  Not dropped: the local rescan re-classifies the neighbourhood and its _check_loud still raises
            # AmbiguousRootError for a genuinely ambiguous / even-order candidate.
            return None
        return tr

    def _local_rescan(self, cs, E_pred, scale, taken):
        """Windowed label-free scan around E_pred at successively finer steps; the zero nearest E_pred.  None -> full rescan."""
        t = self.track
        w = max(t["local_rows"] * GRID_DE, 3.0 * abs(scale))
        lo, hi = max(self.E_min, E_pred - w), E_pred + w
        self.n_local_rescans += 1
        for step in t["local_steps_rows"]:
            opts = dict(self.scan_options or {})
            opts["step_rows"] = float(step)
            opts["fine_step_rows"] = min(float(step), opts.get("fine_step_rows", _rootscan.SCAN_DEFAULTS["fine_step_rows"]))
            res = _rootscan.scan_sector(self._ev, self.sector, cs, lo, hi, opts)
            self._check_loud(res, f"local rescan [{lo:.7f}, {hi:.7f}]")
            self._poles = [p for p in self._poles if not (lo <= p <= hi)] + [float(p) for p in res.poles()]
            zs = [z for z in res.zeros() if self._distinct(z, taken)]
            if zs:
                E = min(zs, key=lambda z: abs(z - E_pred))
                return self._make_root(float(E), cs)
        return None

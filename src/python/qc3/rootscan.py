"""Label-free, complete root finder per sector (prompts/003 section 2.2; D8; docs/METHODS.md section 12).

scan_sector(evaluate, sector, c, E_min, E_max) finds every zero and pole of det Q(E) in [E_min, E_max] without labels:

  1. split at the exact basis-dimension jumps (qc3.fast.dimension_jumps; an evaluator may carry its own
     `dimension_jumps(sector, E_min, E_max)` attribute, used instead -- this is how the unit tests inject a jump);
  2. coarse grid of spacing <= step_rows rows inside each segment (edge `edge` away from the ends);
  3. sign changes -> Brent (qc3.roots) -> Delta_h classification -> order estimate -> slope;
  4. grid-local minima of ln|det Q| with equal signs -> bounded minimisation (even-order zeros);
  4b. grid-local maxima with equal signs -> fine rescan of the 2-step window (zero-pole pair inside one coarse step);
  7. free-level probe (D8): zero-like candidates within free_tol of a non-interacting level become kind "free_level";
  8. fine rescans (fine_step_rows) around every candidate, and the 4b windows.

Kinds: "zero" (odd order, sign change; the physical roots), "pole", "even_zero", "free_level" (D8 artefact: excluded from
zeros()), "ambiguous".  The slope u^T Q' u is recorded for every zero-like candidate and never used to drop a root (D9).
SingularQError at the centre of a classification (E exactly on a root) is the documented exact-root case: ln|det| = -inf.
"""
from __future__ import annotations

import csv
import pathlib
import time
from dataclasses import dataclass

import numpy as np
from scipy.optimize import minimize_scalar

from . import fast as _fast
from . import roots as _roots
from .constants import GRID_DE
from .free_levels import free_levels as _free_levels
from .point import SingularQError
from .sector import Sector

SCAN_DEFAULTS = dict(step_rows=4.0, fine_step_rows=0.25, edge=2e-6, h=1e-6, h2=3e-6, D=4.0, xtol=1e-11, min_xatol=1e-8,
                     free_tol=1e-7, order_tol=0.3)
KINDS = ("zero", "pole", "even_zero", "free_level", "ambiguous")
_ZERO_LIKE = ("zero", "even_zero", "ambiguous")
_POLISH_HALF_WIDTH = 2e-8       # item 4: polish the minimum on +-2e-8
_POLISH_XATOL = 1e-12
_PROBE_OFFSET = 2.5e-9          # free-level probe: centre point offset from E_f (an even zero is within ~6e-9 of it)

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
    order: float | None        # estimated order of the zero (m of section 2.2 item 5); None for poles / unclassifiable
    order_int: int | None
    sign_change: bool
    discriminant: float        # Delta_h
    slope: float | None        # u^T Q'(E) u for zero-like kinds, None for poles
    N: int
    segment: tuple
    bracket: tuple
    via: str                   # "sign" | "min" | "max" | "free"  ("max": found in the fine rescan of an item-4b window)
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
    grid: np.ndarray           # the coarse grid energies of all segments, sorted; n_grid == len(grid)

    def zeros(self) -> np.ndarray:
        """E of kind == "zero", sorted."""
        return np.array(sorted(c.E for c in self.candidates if c.kind == "zero"), dtype=float)

    def poles(self) -> np.ndarray:
        return np.array(sorted(c.E for c in self.candidates if c.kind == "pole"), dtype=float)

    def of_kind(self, kind: str) -> list:
        if kind not in KINDS:
            raise ValueError(f"unknown kind {kind!r}; allowed: {KINDS}")
        return [c for c in self.candidates if c.kind == kind]


# ---------------------------------------------------------------------------------------------------------------------
# evaluation wrapper
# ---------------------------------------------------------------------------------------------------------------------
class _Ev:
    """Counts evaluate() calls.  Calling it returns (sign, ln|det Q|, N) as qc3.roots expects; .full() the PointResult."""

    def __init__(self, evaluate, sector, c):
        self.evaluate, self.sector, self.c, self.n = evaluate, sector, tuple(float(x) for x in c), 0

    def full(self, E):
        self.n += 1
        return self.evaluate(self.sector, float(E), self.c, want=("Q",))

    def __call__(self, E):
        r = self.full(E)
        return int(r.sign), float(r.lnabsdet), int(r.N)


def _check_n(E, **Ns):
    if len(set(Ns.values())) != 1:
        raise _roots.DimensionJumpError(f"basis dimension changes near E={E!r}: {Ns}")


def _slope(pm, pp, h):
    """u^T Q'(E) u with Q' = [Q(E+h) - Q(E-h)]/2h and u the unit eigenvector of the smallest-|lambda| eigenvalue of
    Qbar = [Q(E+h) + Q(E-h)]/2 (prompts/003 section 2.2 item 6, eq. (1)); Q is real symmetric (imaginary part ~ 1e-16)."""
    Qm, Qp = np.asarray(pm.Q), np.asarray(pp.Q)
    Qbar = 0.5 * (Qp + Qm)
    scale = float(np.max(np.abs(Qbar))) or 1.0
    if np.max(np.abs(Qbar.imag)) > 1e-10 * scale:
        raise ValueError(f"slope: Q is not real (max|Im Q| = {np.max(np.abs(Qbar.imag)):.3e})")
    A = Qbar.real
    w, v = np.linalg.eigh(0.5 * (A + A.T))
    u = v[:, int(np.argmin(np.abs(w)))]
    dQ = (Qp.real - Qm.real) / (2.0 * h)
    return float(u @ dQ @ u)


def _characterize(ev, E, bracket, segment, via, o):
    """Classify a candidate at E (items 3, 5, 6): Delta_h at h and h2 -> zero / pole / even_zero / ambiguous, order, slope."""
    n0 = ev.n
    h1, h2, D = o["h"], o["h2"], o["D"]
    pm, pp = ev.full(E - h1), ev.full(E + h1)
    try:
        s0, l0, N0 = ev(E)
    except SingularQError:
        # E is exactly a root of det Q (zero pivot): ln|det| = -inf there (documented exact-root case, prompts/003 section 2.5)
        s0, l0, N0 = 0, -np.inf, int(pm.N)
    sm, lm, Nm = int(pm.sign), float(pm.lnabsdet), int(pm.N)
    sp, lp, Np = int(pp.sign), float(pp.lnabsdet), int(pp.N)
    _check_n(E, minus=Nm, centre=N0, plus=Np)
    if not (np.isfinite(lm) and np.isfinite(lp)) or np.isnan(l0) or l0 == np.inf:
        raise ValueError(f"scan: non-finite ln|det Q| around E={E!r}: {lm}, {l0}, {lp}")
    disc = float(l0 - 0.5 * (lm + lp))
    sign_change = bool(sm * sp < 0)
    kw = dict(E=float(E), sign_change=sign_change, discriminant=disc, N=N0, segment=tuple(segment),
              bracket=tuple(bracket), via=via, free_level=None)
    if disc >= D:
        return RootCandidate(kind="pole", order=None, order_int=None, slope=None, n_evaluations=ev.n - n0, **kw)
    slope = _slope(pm, pp, h1)
    if disc > -D:
        return RootCandidate(kind="ambiguous", order=None, order_int=None, slope=slope, n_evaluations=ev.n - n0, **kw)
    # zero-like: order m = (lbar(h2) - lbar(h1)) / ln(h2/h1)  (= (Delta_h1 - Delta_h2)/ln(h2/h1); the centre drops out)
    s2m, l2m, N2m = ev(E - h2)
    s2p, l2p, N2p = ev(E + h2)
    _check_n(E, minus2=N2m, centre=N0, plus2=N2p)
    m = (0.5 * (l2m + l2p) - 0.5 * (lm + lp)) / np.log(h2 / h1)
    mi = int(round(m))
    if abs(m - mi) > o["order_tol"] or mi < 1:
        kind = "ambiguous"
    elif mi % 2 == 1 and sign_change:
        kind = "zero"
    elif mi % 2 == 0 and not sign_change:
        kind = "even_zero"
    else:
        kind = "ambiguous"
    return RootCandidate(kind=kind, order=float(m), order_int=mi, slope=slope, n_evaluations=ev.n - n0, **kw)


def _disc_at(ev, E, h, E_centre=None):
    """Delta_h(E) (3 evaluations); -inf if E is exactly a root.  E_centre: evaluate the centre there instead of at E (used by
    the free-level probe: F has a pole exactly at a free level, the fast path then raises NONFINITE_COMPONENT)."""
    sm, lm, Nm = ev(E - h)
    sp, lp, Np = ev(E + h)
    try:
        s0, l0, N0 = ev(E if E_centre is None else E_centre)
    except SingularQError:
        return -np.inf
    _check_n(E, minus=Nm, centre=N0, plus=Np)
    return float(l0 - 0.5 * (lm + lp))


def _minimise_ell(ev, lo, hi, xatol):
    """argmin of ln|det Q| on [lo, hi] (scipy 'bounded' = Brent's fmin).  The search runs on t = (x - lo)/(hi - lo) in [0, 1]:
    scipy's tolerance sqrt(eps)|x| + xatol/3 would otherwise be ~4e-9 at x = 0.28, regardless of xatol."""
    w = float(hi) - float(lo)
    res = minimize_scalar(lambda t: ev(lo + w * t)[1], bounds=(0.0, 1.0), method="bounded",
                          options={"xatol": xatol / w, "maxiter": 500})
    if not res.success:
        raise RuntimeError(f"scan: bounded minimisation of ln|det Q| on [{lo!r}, {hi!r}] did not converge: {res.message}")
    return float(lo + w * res.x)


def _polished_minimum_candidate(ev, E0, lo, hi, segment, bracket, via, o):
    """Polish a minimum of ln|det Q| near E0 (+-2e-8, xatol 1e-12) and classify it.  A candidate that turns out to have a
    sign change (a simple zero found by the minimiser) is re-refined by Brent between E-h and E+h."""
    a, b = max(lo, E0 - _POLISH_HALF_WIDTH), min(hi, E0 + _POLISH_HALF_WIDTH)
    E = _minimise_ell(ev, a, b, _POLISH_XATOL)
    cand = _characterize(ev, E, bracket, segment, via, o)
    if cand.sign_change:
        n0 = ev.n
        E = _roots._brent(ev, E - o["h"], E + o["h"], o["xtol"], 200)
        cand = _characterize(ev, E, bracket, segment, via, o)
        cand.n_evaluations += ev.n - n0 - cand.n_evaluations
    return cand


# ---------------------------------------------------------------------------------------------------------------------
# the scan
# ---------------------------------------------------------------------------------------------------------------------
class _Scan:
    def __init__(self, ev, o, step, fine_step, segs):
        self.ev, self.o, self.step, self.fine_step = ev, o, step, fine_step
        self.segs = segs                 # list of (a, b, lo, hi) with hi > lo
        self.found = []                  # insertion order = priority for duplicate merging
        self.max_windows = []            # item-4b windows (lo, hi, seg)

    def known_in(self, lo, hi, sign_change_only=False):
        tol = 10.0 * self.o["xtol"]
        return any(lo - tol <= c.E <= hi + tol and (c.sign_change or not sign_change_only) for c in self.found)

    def block(self, Es, s, l, N, seg, fine):
        """Items 3, 4, 4b on one grid block (coarse or fine)."""
        ev, o = self.ev, self.o
        a, b = seg[0], seg[1]
        if len(set(N.tolist())) != 1:
            raise _roots.DimensionJumpError(f"basis dimension is not constant on the segment ({a!r}, {b!r}): N = {sorted(set(N.tolist()))}")
        n = len(Es)
        for k in range(n - 1):
            if s[k] * s[k + 1] >= 0:
                continue
            if fine and self.known_in(Es[k], Es[k + 1], sign_change_only=True):
                continue
            n0 = ev.n
            E = _roots._brent(ev, float(Es[k]), float(Es[k + 1]), o["xtol"], 200, ends=((int(s[k]), float(l[k])), (int(s[k + 1]), float(l[k + 1]))))
            via = "sign"
            if fine and any(w[0] <= Es[k] and Es[k + 1] <= w[1] + 1e-15 for w in self.max_windows):
                via = "max"
            cand = _characterize(ev, E, (float(Es[k]), float(Es[k + 1])), (a, b), via, o)
            cand.n_evaluations = ev.n - n0
            self.found.append(cand)
        for k in range(1, n - 1):
            if not (s[k - 1] == s[k] == s[k + 1]):
                continue
            if l[k] < l[k - 1] and l[k] <= l[k + 1]:
                lo, hi = float(Es[k - 1]), float(Es[k + 1])
                if fine and self.known_in(lo, hi):
                    continue
                n0 = ev.n
                Em = _minimise_ell(ev, lo, hi, o["min_xatol"])
                if _disc_at(ev, Em, o["h"]) > -o["D"]:
                    continue                                   # smooth minimum
                cand = _polished_minimum_candidate(ev, Em, lo, hi, (a, b), (lo, hi), "min", o)
                cand.n_evaluations = ev.n - n0
                self.found.append(cand)
            elif (not fine) and l[k] > l[k - 1] and l[k] >= l[k + 1]:
                self.max_windows.append((float(Es[k - 1]), float(Es[k + 1]), seg))   # item 4b

    def fine_rescan(self):
        """Item 8: +-1 coarse step around every candidate and over every item-4b window, at the fine step."""
        ivs = []
        for c in self.found:
            seg = next(sg for sg in self.segs if sg[0] == c.segment[0] and sg[1] == c.segment[1])
            ivs.append((max(seg[2], c.E - self.step), min(seg[3], c.E + self.step), seg))
        ivs += [(w[0], w[1], w[2]) for w in self.max_windows]
        ivs.sort(key=lambda t: (t[2][0], t[0]))
        merged = []
        for lo, hi, seg in ivs:
            if merged and merged[-1][2][0] == seg[0] and lo <= merged[-1][1] + 1e-15:
                merged[-1][1] = max(merged[-1][1], hi)
            else:
                merged.append([lo, hi, seg])
        n_fine = 0
        for lo, hi, seg in merged:
            n = int(np.ceil((hi - lo) / self.fine_step - 1e-9)) + 1
            Es = np.linspace(lo, hi, max(n, 2))
            pts = [self.ev(E) for E in Es]
            s = np.array([p[0] for p in pts]); l = np.array([p[1] for p in pts]); N = np.array([p[2] for p in pts])
            self.block(Es, s, l, N, seg, fine=True)
            n_fine += len(Es)
        return n_fine


def _dedupe(found, tol):
    out = []
    for c in found:
        if any(abs(c.E - k.E) <= tol for k in out):
            continue
        out.append(c)
    return out


def scan_sector(evaluate, sector: Sector, c, E_min: float, E_max: float, options: dict | None = None) -> ScanResult:
    """Find every zero and pole of det Q(E; c) in [E_min, E_max] without labels (prompts/003 section 2.2).
    `evaluate(sector, E, c, want=("Q",)) -> PointResult` (qc3.fast.evaluate_point).  Single thread, no disk cache."""
    t0 = time.perf_counter()
    o = dict(SCAN_DEFAULTS)
    if options:
        unknown = set(options) - set(SCAN_DEFAULTS)
        if unknown:
            raise ValueError(f"scan_sector: unknown options {sorted(unknown)}; allowed: {sorted(SCAN_DEFAULTS)}")
        o.update(options)
    _roots._check_xtol(o["xtol"], o["h"], o["D"])
    if not (o["h"] < o["h2"]):
        raise ValueError("scan_sector: need h < h2")
    E_min, E_max = float(E_min), float(E_max)
    if not E_min < E_max:
        raise ValueError("scan_sector: need E_min < E_max")
    c = tuple(float(x) for x in c)
    ev = _Ev(evaluate, sector, c)
    jumps_fn = getattr(evaluate, "dimension_jumps", None) or _fast.dimension_jumps
    jumps = np.asarray(jumps_fn(sector, E_min, E_max), dtype=float)
    edges = [E_min] + [float(x) for x in jumps] + [E_max]
    segments = [(a, b) for a, b in zip(edges[:-1], edges[1:])]
    step, fine_step = o["step_rows"] * GRID_DE, o["fine_step_rows"] * GRID_DE
    segs = [(a, b, a + o["edge"], b - o["edge"]) for a, b in segments if b - o["edge"] > a + o["edge"]]
    scan = _Scan(ev, o, step, fine_step, segs)

    grids = []
    for seg in segs:
        lo, hi = seg[2], seg[3]
        n = int(np.ceil((hi - lo) / step - 1e-9)) + 1
        Es = np.linspace(lo, hi, max(n, 2))
        pts = [ev(E) for E in Es]
        s = np.array([p[0] for p in pts]); l = np.array([p[1] for p in pts]); N = np.array([p[2] for p in pts])
        scan.block(Es, s, l, N, seg, fine=False)
        grids.append(Es)
    grid = np.concatenate(grids) if grids else np.zeros(0)
    n_after_coarse = ev.n

    # item 7: free-level probe -- a free-level solution that the grid rules missed
    free_E = _free_levels(sector, E_min, E_max)
    for Ef in free_E:
        Ef = float(Ef)
        seg = next((sg for sg in segs if sg[2] + o["h2"] <= Ef <= sg[3] - o["h2"]), None)
        if seg is None:
            continue
        if any(abs(k.E - Ef) <= o["free_tol"] and k.kind in _ZERO_LIKE for k in scan.found):
            continue
        n0 = ev.n
        if _disc_at(ev, Ef, o["h"], E_centre=Ef + _PROBE_OFFSET) > -o["D"]:
            continue
        cand = _polished_minimum_candidate(ev, Ef, seg[2], seg[3], (seg[0], seg[1]), (Ef - _POLISH_HALF_WIDTH, Ef + _POLISH_HALF_WIDTH), "free", o)
        cand.n_evaluations = ev.n - n0
        scan.found.append(cand)

    # item 8: fine rescans
    scan.fine_rescan()

    found = _dedupe(scan.found, 10.0 * o["xtol"])
    for cnd in found:                                             # item 7: the D8 labelling
        if cnd.kind in _ZERO_LIKE and free_E.size:
            j = int(np.argmin(np.abs(free_E - cnd.E)))
            if abs(free_E[j] - cnd.E) <= o["free_tol"]:
                # An even-order zero is only located to ~3e-10 (ln|det Q| has a round-off floor at ~x^m; the minimum is a
                # plateau), while the free level is exact to ~1e-16 and the solution sits exactly on it (arXiv:1901.07095
                # Sec. 4.4.3): E is set to the free-level energy.  The other fields (order, slope, discriminant) are kept.
                cnd.kind, cnd.free_level, cnd.E = "free_level", float(free_E[j]), float(free_E[j])
    found.sort(key=lambda t: t.E)
    return ScanResult(sector=sector, c=c, E_min=E_min, E_max=E_max, candidates=found, jumps=jumps, segments=segments,
                      n_evaluations=ev.n, n_grid=int(len(grid)), seconds=time.perf_counter() - t0, options=o, grid=grid)


# ---------------------------------------------------------------------------------------------------------------------
# labels
# ---------------------------------------------------------------------------------------------------------------------
def read_labels(path) -> list:
    """Rows {"bracket_id": int, "E_left": float, "E_right": float, "label": "true"|"false"|""}; `path` absolute or
    relative to qc3.lattice.default_data_dir() (as the LABEL_FILES values)."""
    p = pathlib.Path(path)
    if not p.is_absolute():
        from .lattice import default_data_dir
        p = pathlib.Path(default_data_dir()) / p
    with open(p) as fh:
        return [{"bracket_id": int(r["bracket_id"]), "E_left": float(r["E_left"]), "E_right": float(r["E_right"]),
                 "label": r["user_label"].strip().lower()} for r in csv.DictReader(fh)]


def score_labels(result: ScanResult, labels: list, slack_rows: float = 0.25, exclude_ids=()) -> dict:
    """{"tp","fp","fn","tn": int, "mismatches": list[dict], "unlabelled_zeros": list[RootCandidate], "excluded": list[int]}

    A TRUE bracket (widened by slack_rows grid rows on each side) must contain exactly one kind-"zero" candidate, a FALSE
    bracket none, and no zero may lie outside every non-excluded bracket (invariant 7; prompts/003 section 2.2 gate)."""
    exclude = set(int(i) for i in exclude_ids)
    zeros = [c for c in result.candidates if c.kind == "zero"]
    tp = fp = fn = tn = 0
    mismatches, covered, excluded = [], set(), []
    for r in labels:
        if r["bracket_id"] in exclude:
            excluded.append(r["bracket_id"])
            continue
        lo, hi = r["E_left"] - slack_rows * GRID_DE, r["E_right"] + slack_rows * GRID_DE
        inside = [z for z in zeros if lo <= z.E <= hi]
        covered.update(id(z) for z in inside)
        if r["label"] == "true":
            if len(inside) == 1:
                tp += 1
            else:
                fn += 1
                mismatches.append({"type": "FN", "bracket_id": r["bracket_id"], "E": [z.E for z in inside]})
        elif r["label"] == "false":
            if not inside:
                tn += 1
            else:
                fp += 1
                mismatches.append({"type": "FP", "bracket_id": r["bracket_id"], "E": [z.E for z in inside]})
        else:
            raise ValueError(f"unlabelled bracket {r['bracket_id']} (label {r['label']!r})")
    return {"tp": tp, "fp": fp, "fn": fn, "tn": tn, "mismatches": mismatches,
            "unlabelled_zeros": [z for z in zeros if id(z) not in covered], "excluded": excluded}

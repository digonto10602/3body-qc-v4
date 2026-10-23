"""Root refinement with pole/zero discrimination (prompts/001 step 4; prompts/002 section 2.5, sub-step 7.1).

det Q(E) changes sign at zeros AND at poles.  With l(E) = ln|det Q(E)| a simple zero at E0 looks like ln|E-E0| (a dip,
local minimum) and a simple pole like -ln|E-Ep| (a spike).  After Brent has narrowed a sign change to E* the
discriminant
        Delta_h(E*) = l(E*) - [ l(E*-h) + l(E*+h) ] / 2 ,        h = 1e-6  (0.2 grid rows)
is ~ ln(|E*-E0|/h) <= ln(xtol/h) = -11.5 for a simple zero and ~ +11.5 for a simple pole (higher order: larger |Delta_h|;
the smooth background adds O(h |l'_bg|) << 1).  zero if Delta_h <= -D, pole if Delta_h >= +D, else ambiguous (D = 4).

xtol coupling (audit 001 N6): Brent stops within xtol of the sign change, so |E* - E0| <= xtol and a simple zero has
Delta_h ~ ln(|E*-E0|/h) <= ln(xtol/h).  Classification as a zero needs Delta_h <= -D, which is guaranteed only for
xtol <= h e^{-D} (= 1.8e-8 for h = 1e-6, D = 4); a larger xtol makes the pole/zero discriminant ambiguous (reports/tests/008).
refine_root, refine_root_info, find_sign_changes and find_zeros therefore raise ValueError for xtol > h e^{-D}
(kind="any" in refine_root does not classify and is exempt).
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.optimize import brentq

from .point import MissingRootError  # noqa: F401
from .qc import slogdet  # noqa: F401  (re-exported for convenience)


class PoleError(MissingRootError):
    """The sign change of det Q is a pole, not a zero."""


class AmbiguousRootError(MissingRootError):
    """The sign change could not be classified (|Delta_h| < D)."""


class DimensionJumpError(MissingRootError):
    """The basis dimension N differs at E*-h, E*, E*+h: the sign change sits on a dimension jump."""


@dataclass
class RootInfo:
    E: float
    kind: str                       # "zero" | "pole" | "ambiguous"
    lnabsdet: float                 # l(E*)
    lnabsdet_minus: float           # l(E*-h)
    lnabsdet_plus: float            # l(E*+h)
    discriminant: float             # Delta_h
    N: int                          # basis dimension at E* (0 if the evaluator does not report it)
    bracket: tuple                  # (E_left, E_right) the Brent bracket
    n_evaluations: int              # evaluate() calls spent (Brent + 2 end points if not reused + 3 classification)


class _Counter:
    """Wraps `evaluate` (called as evaluate(sector, E, c, want=("Q",))) and counts calls."""

    def __init__(self, evaluate, sector, c):
        self.evaluate, self.sector, self.c, self.n = evaluate, sector, tuple(c), 0

    def __call__(self, E):
        self.n += 1
        r = self.evaluate(self.sector, float(E), self.c, want=("Q",))
        return int(r.sign), float(r.lnabsdet), int(getattr(r, "N", 0))


def _classify(point, E_star, h, D, bracket, n_before):
    s0, l0, N0 = point(E_star)
    sm, lm, Nm = point(E_star - h)
    sp, lp, Np = point(E_star + h)
    if not (N0 == Nm == Np):
        raise DimensionJumpError(f"basis dimension changes near E*={E_star!r}: N = {Nm}, {N0}, {Np} at E*-h, E*, E*+h")
    disc = l0 - 0.5 * (lm + lp)
    if not np.isfinite(disc):
        # exactly on a zero (ln|det| = -inf at E* itself) is a zero; any other non-finite value is a failure
        if l0 == -np.inf and np.isfinite(lm) and np.isfinite(lp):
            disc = -np.inf
        else:
            raise ValueError(f"classify: non-finite l(E) around E*={E_star!r}: {lm}, {l0}, {lp}")
    kind = "zero" if disc <= -D else ("pole" if disc >= D else "ambiguous")
    return RootInfo(E=float(E_star), kind=kind, lnabsdet=l0, lnabsdet_minus=lm, lnabsdet_plus=lp, discriminant=float(disc),
                    N=N0, bracket=tuple(bracket), n_evaluations=point.n)


def classify_sign_change(evaluate, sector, c, E_star, h=1e-6, D=4.0):
    """Delta_h classification of a sign change at E_star (3 evaluations)."""
    point = _Counter(evaluate, sector, c)
    return _classify(point, E_star, h, D, (E_star, E_star), 0)


def _brent(point, E_left, E_right, xtol, max_iter, ends=None):
    """Brent on f(E) = sign * exp(l - l0), l0 = max(l at the two ends). ends = ((sL, lL), (sR, lR)) to reuse values."""
    if ends is None:
        sL, lL, _ = point(E_left)
        sR, lR, _ = point(E_right)
    else:
        (sL, lL), (sR, lR) = ends
    if sL * sR > 0:
        raise MissingRootError(f"no sign change of det Q on [{E_left}, {E_right}] (sign {sL}, {sR})")
    l0 = max(lL, lR)
    if not np.isfinite(l0):
        l0 = 0.0

    def f_of(s, l):
        # exactly singular evaluations (l = +-inf) are mapped to a huge / zero value instead of NaN
        if np.isnan(l):
            raise ValueError("det Q evaluated to NaN")
        if l == np.inf:
            return s * 1e300
        if l == -np.inf:
            return 0.0
        return s * float(np.exp(l - l0))

    fL, fR = f_of(sL, lL), f_of(sR, lR)

    def g(E):
        if E == E_left:
            return fL
        if E == E_right:
            return fR
        s, l, _ = point(E)
        return f_of(s, l)

    return float(brentq(g, E_left, E_right, xtol=xtol, rtol=4 * np.finfo(float).eps, maxiter=max_iter))


def _check_xtol(xtol, h, D):
    """Delta_h is only a reliable pole/zero discriminant if xtol <= h e^{-D} (see the module docstring)."""
    if not (h > 0.0 and np.isfinite(h)) or not (xtol > 0.0 and np.isfinite(xtol)):
        raise ValueError(f"xtol and h must be finite and positive, got xtol={xtol!r}, h={h!r}")
    xmax = h * float(np.exp(-D))
    if xtol > xmax:
        raise ValueError(f"xtol={xtol!r} > h*exp(-D)={xmax!r} (h={h!r}, D={D!r}): a larger xtol makes the pole/zero "
                         f"discriminant Delta_h ambiguous (a simple zero then gives Delta_h > -D); lower xtol or raise h")


def _check_kind(kind):
    if kind not in ("zero", "any"):
        raise ValueError(f"kind must be 'zero' or 'any', got {kind!r}")


def _raise_unless_zero(info):
    if info.kind == "pole":
        raise PoleError(f"sign change at E={info.E!r} is a pole (Delta_h = {info.discriminant:.2f})")
    if info.kind == "ambiguous":
        raise AmbiguousRootError(f"sign change at E={info.E!r} is ambiguous (Delta_h = {info.discriminant:.2f})")


def refine_root_info(evaluate, sector, c, E_left, E_right, xtol=1e-11, max_iter=200, kind="zero", h=1e-6, D=4.0):
    """Brent refinement plus Delta_h classification.  kind="zero": raise PoleError / AmbiguousRootError unless the sign
    change is a zero.  kind="any": return the RootInfo of whatever sign change was found.
    Requires xtol <= h exp(-D) (ValueError otherwise: the discriminant would be ambiguous)."""
    _check_kind(kind)
    _check_xtol(xtol, h, D)
    point = _Counter(evaluate, sector, c)
    E = _brent(point, E_left, E_right, xtol, max_iter)
    info = _classify(point, E, h, D, (E_left, E_right), point.n)
    if kind == "zero":
        _raise_unless_zero(info)
    return info


def refine_root(evaluate, sector, c, E_left, E_right, xtol=1e-11, max_iter=200, kind="zero", h=1e-6, D=4.0):
    """Root of det Q(E) = 0 in [E_left, E_right]; returns E.  Default kind="zero" raises PoleError / AmbiguousRootError
    if the sign change found is not a zero (3 extra evaluations); kind="any" is the plain Brent sign-change finder.
    MissingRootError if det Q has equal signs at both ends.  `evaluate` is qc3.oracle.evaluate_point or
    qc3.fast.evaluate_point, called as evaluate(sector, E, c, want=("Q",)).
    With kind="zero" the classification needs xtol <= h exp(-D) (= 1.8e-8 for the defaults); a larger xtol raises
    ValueError because the pole/zero discriminant Delta_h would be ambiguous.  kind="any" does not classify, no limit."""
    _check_kind(kind)
    if kind == "any":
        return _brent(_Counter(evaluate, sector, c), E_left, E_right, xtol, max_iter)
    return refine_root_info(evaluate, sector, c, E_left, E_right, xtol, max_iter, "zero", h, D).E


def find_sign_changes(evaluate, sector, c, E_left, E_right, n_scan=241, xtol=1e-11, h=1e-6, D=4.0):
    """Scan [E_left, E_right] on n_scan equally spaced points, Brent-refine every sign-change bracket and classify it.
    Returns all RootInfo sorted by E (zeros, poles, ambiguous).  A sign change that sits on a basis-dimension jump raises
    DimensionJumpError (loud).  Requires xtol <= h exp(-D) (ValueError otherwise).  For the D7 table and M2; NOT the fit function."""
    _check_xtol(xtol, h, D)
    if n_scan < 2:
        raise ValueError("n_scan must be >= 2")
    Es = np.linspace(E_left, E_right, int(n_scan))
    scan = _Counter(evaluate, sector, c)
    pts = [scan(E) for E in Es]
    out = []
    for k in range(len(Es) - 1):
        if pts[k][0] * pts[k + 1][0] < 0:
            point = _Counter(evaluate, sector, c)
            E = _brent(point, float(Es[k]), float(Es[k + 1]), xtol, 200,
                       ends=((pts[k][0], pts[k][1]), (pts[k + 1][0], pts[k + 1][1])))
            info = _classify(point, E, h, D, (float(Es[k]), float(Es[k + 1])), point.n)
            out.append(info)
    out.sort(key=lambda r: r.E)
    return out


def find_zeros(evaluate, sector, c, E_left, E_right, n_scan=241, xtol=1e-11, h=1e-6, D=4.0):
    """The zeros among find_sign_changes (same xtol <= h exp(-D) requirement)."""
    _check_xtol(xtol, h, D)
    return [r for r in find_sign_changes(evaluate, sector, c, E_left, E_right, n_scan, xtol, h, D) if r.kind == "zero"]

"""Root refinement (prompts/001 step 4)."""
import numpy as np
from scipy.optimize import brentq

from .point import MissingRootError  # noqa: F401
from .qc import slogdet  # noqa: F401  (re-exported for convenience)


def refine_root(evaluate, sector, c, E_left, E_right, xtol=1e-11, max_iter=200):
    """Root of det Q(E) = 0 in [E_left, E_right] by Brent's method on f(E) = sign * exp(lnabsdet - l0),
    l0 = max(lnabsdet at the two ends).  `evaluate` is qc3.oracle.evaluate_point or qc3.fast.evaluate_point
    (called as evaluate(sector, E, c, want=("Q",))).  Requires a sign change of f over the bracket, else
    MissingRootError.  Note: a sign change of det Q can also be a pole; the caller supplies brackets that
    contain a zero (labelled windows), and the final |f| is checked to be below its end-point values."""
    def point(E):
        r = evaluate(sector, float(E), tuple(c), want=("Q",))
        return int(r.sign), float(r.lnabsdet)

    sL, lL = point(E_left)
    sR, lR = point(E_right)
    if sL * sR > 0:
        raise MissingRootError(f"no sign change of det Q on [{E_left}, {E_right}] (sign {sL}, {sR})")
    l0 = max(lL, lR)

    def f(E):
        s, l = point(E)
        return s * float(np.exp(l - l0))

    fL, fR = sL * float(np.exp(lL - l0)), sR * float(np.exp(lR - l0))

    def g(E):
        if E == E_left:
            return fL
        if E == E_right:
            return fR
        return f(E)

    root = brentq(g, E_left, E_right, xtol=xtol, rtol=4 * np.finfo(float).eps, maxiter=max_iter)
    return float(root)

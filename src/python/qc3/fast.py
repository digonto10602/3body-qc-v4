"""Fast backend, cache-free, single core (prompts/001 step 6 rungs 1-4, prompts/002 steps 7.0-7.6).

Implemented in src/qc3_core/fast_point.hpp: the threshold basis builder (fast_basis.hpp, bit-identical to the oracle's),
the K2^{-1} of the oracle, the pair kernels for K_df,3 (fast_k3.hpp) and G (fast_g.hpp: same formulas as
K3_functions_2plus1.hpp / G_functions_v2.h, one kinematics evaluation per unordered momentum pair), the
orbit/character projector (fast_projector.hpp); the F lattice sum (prompts/002 7.7, fast_fclass.hpp) uses class moment
tables: per canonical nP (cubic-group orbit of the pair momentum) the cube points are grouped into classes of equal
(n^2, n.nP), the ellipsoid r^2 <= R^2 (per-term bound <= trunc_eps, same point set as the direct sum) is selected per class,
momenta with equal (canonical nP, |p|^2) share one sum which is rotated back; real arithmetic.  options["use_classes"]=False
selects the direct term-by-term sum (cross-check).
Energy-independent pieces (the projector V_sel per basis signature, the basis threshold tables) are cached in memory
per process.  The tail is in C++ (Eigen, one thread, no BLAS): F3 from the block structure of F,
Q = L^-6 V^dag (X/s + K3 V) with X = (F3/s)^-1 V by one LU solve with N_proj right-hand sides (mathematically the
projection of the inverse, invariant 1), (sign, ln|det Q|) by partial-pivot LU (invariant 2).  The full inverse F3^-1
is formed only if "F3inv" is in `want`.  This module does not import qc3.qc (no BLAS call, hence no BLAS threads)."""
import numpy as np

from .oracle import K2_PWAVE_FORMS, _metadata
from .point import ALL, BASIS_DTYPE, PointResult
from .sector import Sector  # noqa: F401

# "projector": "orbit" (default, prompts/002 section 2.1) builds V_sel from the little-group orbits and characters (no dense
# projector, no eigen-solve, mandatory self-check P V = V); "eigen" = the existing eigen-decomposition builder, kept for
# cross-checks.  Both go through the fast-owned in-memory cache with a complete key (step 7.0).
PROJECTOR_METHODS = ("orbit", "eigen")
DEFAULT_OPTIONS = {"trunc_eps": 1e-16, "use_classes": True, "threads": 1, "projector": "orbit"}


def _params(sector, Ecm, c, compute_k3):
    return dict(
        irrep=sector.irrep, Lbyas=float(sector.Lbyas), xi=sector.xi, atmK=sector.atmK, atmpi=sector.atmpi,
        eta_1=sector.eta_1, eta_2=sector.eta_2, alpha=sector.alpha, epsilon_H=sector.epsilon_H,
        max_shell_num=int(sector.max_shell_num), waves_1=list(sector.waves_1), waves_2=list(sector.waves_2),
        parity=int(sector.parity), eig_tol=sector.eig_tol, a0_Kpi=sector.a0_Kpi, a1_Kpi=sector.a1_Kpi,
        a0_KK=sector.a0_KK, k2_pwave_form=sector.k2_pwave_form, Ecm=float(Ecm), c=list(c),
        compute_k3=bool(compute_k3),
    )


def evaluate_point(sector, Ecm, c=(0.0, 0.0, 0.0, 0.0), want=ALL, options=None):
    opts = dict(DEFAULT_OPTIONS)
    if options:
        unknown = set(options) - set(DEFAULT_OPTIONS)
        if unknown:
            raise ValueError(f"fast: unknown options {sorted(unknown)}")
        opts.update(options)
    if int(opts["threads"]) != 1:
        raise NotImplementedError("qc3.fast: only threads=1 is implemented (OpenMP is an optional later rung)")
    if opts["projector"] not in PROJECTOR_METHODS:
        raise ValueError(f"fast: projector must be one of {PROJECTOR_METHODS}, got {opts['projector']!r}")
    trunc_eps = float(opts["trunc_eps"])
    if not (0.0 < trunc_eps < 1.0):
        raise ValueError("fast: trunc_eps must be in (0, 1)")
    if sector.k2_pwave_form not in K2_PWAVE_FORMS:
        raise ValueError(f"fast: unknown k2_pwave_form {sector.k2_pwave_form!r}; allowed: {K2_PWAVE_FORMS}")
    from . import _core as core  # lazy: the compiled extension is never imported at package import time

    c = tuple(float(x) for x in c)
    if len(c) != 4:
        raise ValueError("c must be (Kiso0, Kiso1, KB, KE)")
    want = set(want)
    params = _params(sector, Ecm, c, bool("K3" in want))
    d = core.fast_point(params, trunc_eps, opts["projector"], sorted(want), use_classes=bool(opts["use_classes"]))
    if tuple(d["nnP"]) != tuple(sector.nnP):
        raise RuntimeError(f"nnP mismatch between Sector {sector.nnP} and C++ parse_label {d['nnP']}")
    if abs(d["sign_imag"]) > 1e-10:   # same check as qc3.qc.slogdet: the determinant phase must be +-1
        raise ValueError(f"fast: determinant phase has imaginary part {d['sign_imag']:.3e}; Q is not Hermitian")
    b = d["basis"]
    basis = np.empty(b.shape[0], dtype=BASIS_DTYPE)
    for j, name in enumerate(BASIS_DTYPE.names):
        basis[name] = b[:, j]

    timings = dict(d["timings"])   # "projection" and "lu" are measured in C++ (prompts/002 step 7.6)

    md = _metadata(sector, core)
    md["backend"] = "fast"
    md["trunc_eps"] = trunc_eps
    md["use_classes"] = bool(opts["use_classes"])
    md["projector_method"] = opts["projector"]
    md["tail"] = "cpp"   # prompts/002 step 7.6: F3, Q and the LU in C++ (Eigen, one thread); no NumPy/BLAS call
    md["projector_cache_hit"] = bool(d["stats"]["projector_cache_hit"])
    md["lattice_terms"] = int(d["stats"]["lattice_terms"])
    md["lattice_cube_terms_per_entry"] = int(d["stats"]["lattice_cube_per_entry"])
    pick = lambda name: d[name] if name in want else None  # noqa: E731
    return PointResult(
        sector=sector, Ecm=float(Ecm), E_lab=d["E_lab"], N1=d["N1"], N2=d["N2"], N=d["N"], Nproj=d["Nproj"],
        basis=basis, c=c,
        F=pick("F"), G=pick("G"), K2inv=pick("K2inv"), F3=pick("F3"), F3inv=pick("F3inv"), Vsel=pick("Vsel"),
        K3_basis=d["K3_basis"], K3=pick("K3"), Q=d["Q"] if "Q" in want else None,
        sign=int(d["sign"]), lnabsdet=float(d["lnabsdet"]), metadata=md, timings=timings,
    )


def clear_caches():
    """Drop the fast path's in-memory caches (projector V_sel; threshold tables from step 7.3): next call is cold."""
    from . import _core as core
    core.clear_caches()


def cache_info():
    """{"projector_entries", "projector_hits", "projector_misses", "threshold_entries"} of the fast path's caches."""
    from . import _core as core
    return dict(core.cache_info())


def dimension_jumps(sector, E_min, E_max):
    """Sorted E_cm strictly inside (E_min, E_max) where the basis dimension N changes (either channel).

    The basis keeps a spectator momentum n iff E_lab > t(n) (docs/METHODS.md section 2, the oracle's own function
    En_min_plus_for_config with zval 0.02); the jumps are the distinct thresholds t converted to E* = sqrt(t^2 - P^2).
    A sign change of det Q across a jump is a change of dimension, not a root (M2's classifier uses this)."""
    from . import _core as core
    E_min, E_max = float(E_min), float(E_max)
    if not E_min < E_max:
        raise ValueError("dimension_jumps: need E_min < E_max")
    return np.asarray(core.dimension_jumps(_params(sector, 0.5 * (E_min + E_max), (0.0, 0.0, 0.0, 0.0), False), E_min, E_max), dtype=float)

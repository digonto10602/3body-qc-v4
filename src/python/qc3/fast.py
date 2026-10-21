"""Fast backend, cache-free, single core (prompts/001 step 6: rungs 1-4).

Implemented in src/qc3_core/fast_point.hpp: the oracle's basis builder, G, K2^{-1}, K3 and projector, but the
F lattice sum is done in one pass per spectator momentum, truncated to the ellipsoid r^2 <= R^2 inside the oracle's
cube (per-term bound <= trunc_eps), in real arithmetic.  Energy-independent pieces (the projector V_sel per basis
signature) are cached in memory per process.  Final Q and (sign, ln|det|) via qc3.qc (same expressions as the
oracle)."""
import time

import numpy as np

from .oracle import K2_PWAVE_FORMS, _metadata
from .point import ALL, BASIS_DTYPE, PointResult
from .qc import project_Q, slogdet
from .sector import Sector  # noqa: F401

DEFAULT_OPTIONS = {"trunc_eps": 1e-16, "use_classes": False, "threads": 1}


def evaluate_point(sector, Ecm, c=(0.0, 0.0, 0.0, 0.0), want=ALL, options=None):
    opts = dict(DEFAULT_OPTIONS)
    if options:
        unknown = set(options) - set(DEFAULT_OPTIONS)
        if unknown:
            raise ValueError(f"fast: unknown options {sorted(unknown)}")
        opts.update(options)
    if opts["use_classes"]:
        raise NotImplementedError("qc3.fast: use_classes=True is prompts/001 step 7")
    if int(opts["threads"]) != 1:
        raise NotImplementedError("qc3.fast: only threads=1 is implemented (OpenMP is an optional later rung)")
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
    params = dict(
        irrep=sector.irrep, Lbyas=float(sector.Lbyas), xi=sector.xi, atmK=sector.atmK, atmpi=sector.atmpi,
        eta_1=sector.eta_1, eta_2=sector.eta_2, alpha=sector.alpha, epsilon_H=sector.epsilon_H,
        max_shell_num=int(sector.max_shell_num), waves_1=list(sector.waves_1), waves_2=list(sector.waves_2),
        parity=int(sector.parity), eig_tol=sector.eig_tol, a0_Kpi=sector.a0_Kpi, a1_Kpi=sector.a1_Kpi,
        a0_KK=sector.a0_KK, k2_pwave_form=sector.k2_pwave_form, Ecm=float(Ecm), c=list(c),
        compute_k3=bool("K3" in want),
    )
    d = core.fast_point(params, trunc_eps)
    if tuple(d["nnP"]) != tuple(sector.nnP):
        raise RuntimeError(f"nnP mismatch between Sector {sector.nnP} and C++ parse_label {d['nnP']}")
    b = d["basis"]
    basis = np.empty(b.shape[0], dtype=BASIS_DTYPE)
    for j, name in enumerate(BASIS_DTYPE.names):
        basis[name] = b[:, j]

    # final assembly: same expressions and LAPACK calls as the oracle (qc3.qc.project_Q / slogdet)
    t0 = time.perf_counter()
    Q = project_Q(d["F3inv"], d["Vsel"], d["K3"], sector.L)
    t_proj = time.perf_counter() - t0
    t0 = time.perf_counter()
    sgn, lndet = slogdet(Q)
    t_lu = time.perf_counter() - t0
    timings = dict(d["timings"])
    timings["projection"] = t_proj
    timings["lu"] = t_lu
    timings["total"] += t_proj + t_lu

    md = _metadata(sector, core)
    md["backend"] = "fast"
    md["trunc_eps"] = trunc_eps
    md["use_classes"] = False
    md["projector_cache_hit"] = bool(d["stats"]["projector_cache_hit"])
    md["lattice_terms"] = int(d["stats"]["lattice_terms"])
    md["lattice_cube_terms_per_entry"] = int(d["stats"]["lattice_cube_per_entry"])
    pick = lambda name: d[name] if name in want else None  # noqa: E731
    return PointResult(
        sector=sector, Ecm=float(Ecm), E_lab=d["E_lab"], N1=d["N1"], N2=d["N2"], N=d["N"], Nproj=d["Nproj"],
        basis=basis, c=c,
        F=pick("F"), G=pick("G"), K2inv=pick("K2inv"), F3=pick("F3"), F3inv=pick("F3inv"), Vsel=pick("Vsel"),
        K3_basis=d["K3_basis"], K3=pick("K3"), Q=Q if "Q" in want else None,
        sign=int(sgn), lnabsdet=float(lndet), metadata=md, timings=timings,
    )


def clear_caches():
    """Drop the in-memory per-process caches (projector V_sel per basis signature): next call is a cold build."""
    from . import _core as core
    core.clear_caches()

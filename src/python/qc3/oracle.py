"""Oracle backend: the existing C++ physics path (build_component_entry + K3mat_2plus1), unchanged,
single-threaded, one energy, no cache (prompts/001 step 2).  Implemented in src/qc3_core/oracle_point.hpp."""
import time

import numpy as np

from .constants import PRODUCTION
from .point import ALL, BASIS_DTYPE, PointResult
from .sector import Sector  # noqa: F401


def _metadata(sector, core):
    md = {k: v for k, v in PRODUCTION.items() if k not in ("scatter",)}
    # sector values for every production key
    for k in ("xi", "atmK", "atmpi", "eta_1", "eta_2", "alpha", "epsilon_H", "max_shell_num", "parity", "eig_tol"):
        md[k] = getattr(sector, k)
    md["waves_1"] = tuple(sector.waves_1)
    md["waves_2"] = tuple(sector.waves_2)
    md["scatter"] = {"a0_Kpi": sector.a0_Kpi, "a1_Kpi": sector.a1_Kpi, "a0_KK": sector.a0_KK}
    md["Lbyas"] = sector.Lbyas
    md["irrep"] = sector.irrep
    md["nnP"] = tuple(sector.nnP)
    md["L"] = sector.L
    md["k2_pwave_form"] = sector.k2_pwave_form
    md["backend"] = "oracle"
    md["git_commit"] = core.git_commit()
    md["threads"] = 1
    return md


def evaluate_point(sector, Ecm, c=(0.0, 0.0, 0.0, 0.0), want=ALL):
    """Oracle evaluation of one point.  `sign`/`lnabsdet` (for `c`) are always filled; `want` selects which
    matrices are returned (others are None).  Q = L^-6 V^dag (F3inv + K3(c)) V, no N_proj rescaling."""
    if sector.k2_pwave_form != "coded_v2":
        raise ValueError(f"oracle: k2_pwave_form {sector.k2_pwave_form!r} not available yet "
                         "(only 'coded_v2'; 'q0_consistent' is prompts/001 step 3)")
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
        a0_KK=sector.a0_KK, Ecm=float(Ecm), c=list(c), compute_k3=bool("K3" in want),
    )
    d = core.oracle_point(params)
    if tuple(d["nnP"]) != tuple(sector.nnP):
        raise RuntimeError(f"nnP mismatch between Sector {sector.nnP} and C++ parse_label {d['nnP']}")
    b = d["basis"]
    basis = np.empty(b.shape[0], dtype=BASIS_DTYPE)
    for j, name in enumerate(BASIS_DTYPE.names):
        basis[name] = b[:, j]
    # Final assembly in NumPy/LAPACK, exactly as the stored anchors were produced (tests/acceptance/m1/data):
    #   Q = V^dag (F3inv + K3) V / L^6 ;  (sign, ln|det Q|) = LAPACK partial-pivot LU of the complex128 Q.
    # Invariants 1-2: the inverse is projected, L^-6 scales the matrix, no N_proj rescaling.  The C++ (Eigen)
    # versions of the same two lines are returned only as a cross-check (metadata "*_eigen*"): at
    # ill-conditioned points (kappa(Q) ~ 5e8..4e12) a 1e-16 difference in summation order changes ln|det| by
    # 1e-10..1e-9, which is the size of the acceptance tolerance.
    L = sector.L
    Vsel, F3inv, K3 = d["Vsel"], d["F3inv"], d["K3"]
    t0 = time.perf_counter()
    Q = Vsel.conj().T @ (F3inv + K3) @ Vsel / L ** 6
    t_proj = time.perf_counter() - t0
    t0 = time.perf_counter()
    sgn, lndet = np.linalg.slogdet(Q)
    t_lu = time.perf_counter() - t0
    if abs(sgn.imag) > 1e-10 or abs(abs(sgn) - 1.0) > 1e-10:
        raise RuntimeError(f"oracle: det(Q) phase {sgn} is not +-1 (Q not Hermitian?)")
    timings = dict(d["timings"])
    timings["projection_eigen"] = timings["projection"]
    timings["lu_eigen"] = timings["lu"]
    timings["projection"] = t_proj
    timings["lu"] = t_lu
    timings["total"] += t_proj + t_lu - timings["projection_eigen"] - timings["lu_eigen"]
    md = _metadata(sector, core)
    md["sign_phase_imag"] = d["sign_imag"]
    md["lnabsdet_eigen_lu"] = d["lnabsdet"]
    md["sign_eigen_lu"] = int(d["sign"])
    pick = lambda name: d[name] if name in want else None  # noqa: E731
    return PointResult(
        sector=sector, Ecm=float(Ecm), E_lab=d["E_lab"], N1=d["N1"], N2=d["N2"], N=d["N"], Nproj=d["Nproj"],
        basis=basis, c=c,
        F=pick("F"), G=pick("G"), K2inv=pick("K2inv"), F3=pick("F3"), F3inv=pick("F3inv"), Vsel=pick("Vsel"),
        K3_basis=d["K3_basis"], K3=pick("K3"), Q=Q if "Q" in want else None,
        sign=int(np.sign(sgn.real)), lnabsdet=float(lndet), metadata=md, timings=timings,
    )

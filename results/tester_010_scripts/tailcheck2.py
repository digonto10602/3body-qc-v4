"""Tail: Q = L^-6 V^T (F3^-1 + K3) V, recomputed from fast's F, G, K2inv with an EXPLICIT np.linalg.inv, vs fast (LU solve) and the oracle."""
import sys, json, subprocess
import numpy as np
sys.path.insert(0, "/home/user/3body-qc-v4/src/python"); sys.path.insert(0, "/home/user/3body-qc-v4/tests/acceptance/m1")
from qc3 import fast, oracle
from m1_helpers import sector, load_anchors, C_BASELINE
commit = subprocess.check_output(["git", "-C", "/home/user/3body-qc-v4", "rev-parse", "--short", "HEAD"]).decode().strip()
rows = []
for p in load_anchors()["points"]:
    s = sector(p["irrep"]); E = p["Ecm"]
    f = fast.evaluate_point(s, E, C_BASELINE)
    F, G, K2 = np.asarray(f.F).real, np.asarray(f.G).real, np.asarray(f.K2inv)
    H = K2 + F + G
    F3 = F / 3.0 - F @ np.linalg.solve(H, F)
    F3inv = np.linalg.inv(F3)
    V = np.asarray(f.Vsel); K3 = np.asarray(f.K3)
    L6 = f.sector.L ** 6
    Q = V.conj().T @ (F3inv + K3) @ V / L6
    sg, ld = np.linalg.slogdet(Q)
    ev = np.linalg.eigvalsh((Q + Q.conj().T) / 2); ld_eig = float(np.log(np.abs(ev)).sum()); sg_eig = int(np.prod(np.sign(ev)))
    herm = float(np.linalg.norm(Q - Q.conj().T) / np.linalg.norm(Q))
    rows.append(dict(irrep=p["irrep"], Ecm=E, N=f.N, Nproj=f.Nproj,
                     Q_rel=float(np.linalg.norm(f.Q - Q) / np.linalg.norm(Q)),
                     F3_rel=float(np.linalg.norm(np.asarray(f.F3).real - F3) / np.linalg.norm(F3)),
                     F3inv_rel=float(np.linalg.norm(np.asarray(f.F3inv).real - F3inv) / np.linalg.norm(F3inv)),
                     dln_numpy_inv=abs(ld - f.lnabsdet), dln_eigvalsh=abs(ld_eig - f.lnabsdet), sign_ok=(int(sg.real) == f.sign == sg_eig),
                     herm_fast_Q=float(np.linalg.norm(f.Q - np.asarray(f.Q).conj().T) / np.linalg.norm(f.Q)),
                     max_imag_Q=float(np.abs(np.asarray(f.Q).imag).max()), V_real=bool(np.abs(V.imag).max() == 0.0),
                     kappa_Q=float(np.linalg.cond(f.Q)), kappa_F3=float(np.linalg.cond(F3))))
    print(rows[-1], flush=True)
# eigen projector: complex V? tail complex branch vs orbit
cplx = []
for irrep, E in (("000_A1m", 0.30), ("100_A2", 0.2915561), ("110_A2", 0.30733), ("111_A2", 0.32036), ("200_A2", 0.28664)):
    s = sector(irrep)
    fo = fast.evaluate_point(s, E, C_BASELINE, want=("Q", "Vsel"), options={"projector": "orbit"})
    fe = fast.evaluate_point(s, E, C_BASELINE, want=("Q", "Vsel"), options={"projector": "eigen"})
    Ve = np.asarray(fe.Vsel)
    # independent Q for the eigen V, from the orbit run's F3inv/K3 (same physics)
    fa = fast.evaluate_point(s, E, C_BASELINE, options={"projector": "eigen"})
    Qe = Ve.conj().T @ (np.asarray(fa.F3inv) + np.asarray(fa.K3)) @ Ve / fa.sector.L ** 6
    cplx.append(dict(irrep=irrep, Ecm=E, Nproj_orbit=fo.Nproj, Nproj_eigen=fe.Nproj, eigen_V_max_imag=float(np.abs(Ve.imag).max()),
                     dln_orbit_eigen=abs(fo.lnabsdet - fe.lnabsdet), sign_equal=(fo.sign == fe.sign), Q_eigen_vs_numpy=float(np.linalg.norm(fe.Q - Qe) / np.linalg.norm(Qe)),
                     herm_eigen_Q=float(np.linalg.norm(fe.Q - np.asarray(fe.Q).conj().T) / np.linalg.norm(fe.Q)), md_tail=fe.metadata["tail"]))
    print(cplx[-1], flush=True)
# error paths through the Python interface
errs = {}
s = sector("100_A2")
for name, c in (("nan_c", (float("nan"), 0.0, 0.0, 0.0)), ("inf_c", (float("inf"), 0.0, 0.0, 0.0))):
    try:
        r = fast.evaluate_point(s, 0.2915561, c, want=("Q",)); errs[name] = "NO RAISE: lnabsdet=%r sign=%r" % (r.lnabsdet, r.sign)
    except Exception as e: errs[name] = "%s: %s" % (type(e).__name__, str(e)[:120])
print(errs)
json.dump({"commit": commit, "c": list(C_BASELINE), "k2_pwave_form": "coded_v2", "anchors": rows, "eigen_complex_branch": cplx, "error_paths": errs}, open("/home/user/3body-qc-v4/results/tester_010_tail.json", "w"), indent=1)

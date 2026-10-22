"""Independent projector check (own NumPy, from scratch).  For every one of the production bases: build the dense
P = (1/|G|) sum_M chi(M) Pi(M) U_orb(M) from signed-permutation matrices (A2 / A1u character = det M, Pi = det M for parity -1),
and compare V V^dag of the fast orbit projector, the fast eigen projector and (optionally) the oracle."""
import sys, json, itertools, subprocess, time, numpy as np
sys.path.insert(0, "/home/user/3body-qc-v4/src/python"); sys.path.insert(0, "/home/user/3body-qc-v4/tests/acceptance/m1")
from qc3 import fast, oracle
from m1_helpers import sector, IRREPS, NNP, C_ZERO
commit = subprocess.check_output(["git", "-C", "/home/user/3body-qc-v4", "rev-parse", "--short", "HEAD"]).decode().strip()
WITH_ORACLE = int(sys.argv[1]) if len(sys.argv) > 1 else 12     # number of bases (largest first) that also get an oracle call; 0 = none, -1 = all
mats = []
for perm in itertools.permutations(range(3)):
    for sg in itertools.product([1, -1], repeat=3):
        M = np.zeros((3, 3), int)
        for r in range(3): M[r, perm[r]] = sg[r]
        mats.append(M)
assert len(mats) == 48
AX = {1: 0, 0: 2, -1: 1}   # m -> cartesian axis (x,z,y)
def dense_P(basis, nnP):
    N = len(basis); idx = {tuple(b): j for j, b in enumerate(basis)}
    LG = [M for M in mats if np.array_equal(M @ np.array(nnP), np.array(nnP))]
    P = np.zeros((N, N))
    for M in LG:
        det = int(round(np.linalg.det(M)))
        chi = det; Pi = det           # A2 (C4v,C2v,C3v) and A1u (Oh): chi = det M ; intrinsic parity of three pseudoscalars
        U = np.zeros((N, N))
        for j, (fl, nx, ny, nz, l, m) in enumerate(basis):
            kp = tuple(int(x) for x in M @ np.array([nx, ny, nz]))
            for mp in range(-l, l + 1):
                d = 1.0 if l == 0 else float(M[AX[mp], AX[m]])
                if d == 0: continue
                U[idx[(fl,) + kp + (l, mp)], j] += Pi * d
        P += chi * U
    return P / len(LG), len(LG)
res = {"commit": commit, "bases": []}
allb = []
for irrep in IRREPS:
    s = sector(irrep)
    jumps = np.asarray(fast.dimension_jumps(s, 0.2631, 0.36))
    edges = np.concatenate([[0.2631], jumps, [0.36]]); mids = 0.5*(edges[1:]+edges[:-1])
    for E in mids:
        allb.append((irrep, float(E)))
print("bases:", len(allb))
rows = []
for irrep, E in allb:
    s = sector(irrep)
    fo = fast.evaluate_point(s, E, C_ZERO, want=("Vsel",), options={"projector": "orbit"})
    fe = fast.evaluate_point(s, E, C_ZERO, want=("Vsel",), options={"projector": "eigen"})
    basis = [tuple(int(x) for x in b) for b in fo.basis]
    t0 = time.perf_counter(); P, nlg = dense_P(basis, NNP[irrep]); tP = time.perf_counter()-t0
    idem = np.linalg.norm(P@P-P)/np.linalg.norm(P); herm = np.linalg.norm(P-P.T)/np.linalg.norm(P)
    nprojP = int(round(np.trace(P)))
    def rel(V): V = np.asarray(V); return float(np.linalg.norm(V@V.conj().T - P)/np.linalg.norm(P))
    Vo = np.asarray(fo.Vsel)
    row = {"irrep": irrep, "Ecm": E, "N": fo.N, "LG": nlg, "Nproj_mine": nprojP, "Nproj_orbit": fo.Nproj, "Nproj_eigen": fe.Nproj,
           "P_idem": float(idem), "P_herm": float(herm), "relVV_orbit": rel(Vo), "relVV_eigen": rel(fe.Vsel),
           "orth_orbit": float(np.linalg.norm(Vo.conj().T@Vo-np.eye(fo.Nproj))), "orth_eigen": float(np.linalg.norm(np.asarray(fe.Vsel).conj().T@np.asarray(fe.Vsel)-np.eye(fe.Nproj))),
           "imag_Vorbit": float(np.abs(Vo.imag).max()), "t_orbit_proj_ms": 1e3*fo.timings["projector"], "t_eigen_proj_ms": 1e3*fe.timings["projector"]}
    rows.append((row, P))
    print(irrep, f"{E:.6f}", "N", fo.N, "Np", nprojP, fo.Nproj, fe.Nproj, "idem %.1e"%idem, "orbit %.1e eigen %.1e"%(row["relVV_orbit"], row["relVV_eigen"]), "orth %.1e"%row["orth_orbit"], flush=True)
# oracle
order = sorted(range(len(rows)), key=lambda i: -rows[i][0]["N"])
sel = order if WITH_ORACLE < 0 else order[:WITH_ORACLE]
# also make sure every irrep has an oracle point
for irrep in IRREPS:
    c = [i for i in order if rows[i][0]["irrep"] == irrep]
    for i in c[:2]:
        if i not in sel: sel.append(i)
for i in sel:
    row, P = rows[i]
    s = sector(row["irrep"])
    t0 = time.perf_counter()
    o = oracle.evaluate_point(s, row["Ecm"], C_ZERO, want=("Vsel",))
    Vor = np.asarray(o.Vsel)
    row["Nproj_oracle"] = o.Nproj; row["relVV_oracle"] = float(np.linalg.norm(Vor@Vor.conj().T-P)/np.linalg.norm(P))
    row["oracle_N"] = o.N
    print("oracle", row["irrep"], f"{row['Ecm']:.6f}", "N", o.N, "Np", o.Nproj, "relVV_oracle %.1e" % row["relVV_oracle"], "%.1fs" % (time.perf_counter()-t0), flush=True)
res["bases"] = [r for r, _ in rows]
json.dump(res, open("/home/user/3body-qc-v4/results/tester_009_projector_check.json", "w"), indent=1)

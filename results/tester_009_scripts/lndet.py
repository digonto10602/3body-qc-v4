"""ln|det Q| orbit vs eigen vs oracle (and an independent eigenvalue-based value with MY dense projector basis), c = baseline."""
import sys, json, subprocess, time, numpy as np
sys.path.insert(0, "/home/user/3body-qc-v4/src/python"); sys.path.insert(0, "/home/user/3body-qc-v4/tests/acceptance/m1")
sys.path.insert(0, "/tmp/claude-0/-home-user-3body-qc-v4/756827e6-de12-5b97-9946-026a972abb7d/scratchpad/t")
from qc3 import fast, oracle
from m1_helpers import sector, IRREPS, NNP, C_BASELINE, L_BOX
import importlib.util
commit = subprocess.check_output(["git", "-C", "/home/user/3body-qc-v4", "rev-parse", "--short", "HEAD"]).decode().strip()
# my dense P, copied from projcheck (same code)
import itertools
mats = []
for perm in itertools.permutations(range(3)):
    for sg in itertools.product([1, -1], repeat=3):
        M = np.zeros((3, 3), int)
        for r in range(3): M[r, perm[r]] = sg[r]
        mats.append(M)
AX = {1: 0, 0: 2, -1: 1}
def dense_P(basis, nnP):
    N = len(basis); idx = {tuple(b): j for j, b in enumerate(basis)}
    LG = [M for M in mats if np.array_equal(M @ np.array(nnP), np.array(nnP))]
    P = np.zeros((N, N))
    for M in LG:
        det = int(round(np.linalg.det(M)))
        U = np.zeros((N, N))
        for j, (fl, nx, ny, nz, l, m) in enumerate(basis):
            kp = tuple(int(x) for x in M @ np.array([nx, ny, nz]))
            for mp in range(-l, l + 1):
                d = 1.0 if l == 0 else float(M[AX[mp], AX[m]])
                if d == 0: continue
                U[idx[(fl,) + kp + (l, mp)], j] += det * det * d
        P += det * 0 + U * 1.0  # chi*Pi = det^2 = 1 already folded in (U carries Pi=det, chi=det)
    return P / len(LG)
# energies: lowest, mid-range, and the largest anchors; 2 per irrep, plus largest N
PTS = {"000_A1m": [0.30000, 0.34538], "100_A2": [0.2915561, 0.33500], "110_A2": [0.30733, 0.34000],
       "111_A2": [0.31000, 0.35103], "200_A2": [0.30500, 0.34738]}
rows = []
for irrep in IRREPS:
    s = sector(irrep)
    for E in PTS[irrep]:
        fo = fast.evaluate_point(s, E, C_BASELINE, want=("Q", "Vsel"), options={"projector": "orbit"})
        fe = fast.evaluate_point(s, E, C_BASELINE, want=("Q", "Vsel"), options={"projector": "eigen"})
        o = oracle.evaluate_point(s, E, C_BASELINE)
        basis = [tuple(int(x) for x in b) for b in o.basis]
        P = dense_P(basis, NNP[irrep])
        w, U = np.linalg.eigh(P)
        Vm = U[:, w > 0.5]
        assert abs(np.linalg.norm(P@P-P)) < 1e-12 and Vm.shape[1] == o.Nproj, (irrep, E, Vm.shape, o.Nproj)
        A = np.asarray(o.F3inv) + np.asarray(o.K3)
        Vo = np.asarray(o.Vsel)
        Qo_raw = Vo.conj().T @ A @ Vo
        scale = np.linalg.norm(Qo_raw) / np.linalg.norm(np.asarray(o.Q))
        Qm = (Vm.conj().T @ A @ Vm) / scale
        ev_m = np.linalg.eigvalsh(0.5*(Qm+Qm.conj().T)); ev_o = np.linalg.eigvalsh(0.5*(np.asarray(o.Q)+np.asarray(o.Q).conj().T))
        ln_m = float(np.sum(np.log(np.abs(ev_m)))); sg_m = int(np.prod(np.sign(ev_m)))
        ln_eig_o = float(np.sum(np.log(np.abs(ev_o))))
        row = {"irrep": irrep, "Ecm": E, "N": o.N, "Nproj": o.Nproj, "Nproj_orbit": fo.Nproj, "Nproj_eigen": fe.Nproj,
               "sign_orbit": int(fo.sign), "sign_eigen": int(fe.sign), "sign_oracle": int(o.sign), "sign_mine": sg_m,
               "ln_orbit": fo.lnabsdet, "ln_eigen": fe.lnabsdet, "ln_oracle": o.lnabsdet, "ln_mine_eigvalsh_myV": ln_m, "ln_oracle_eigvalsh": ln_eig_o,
               "d_orbit_eigen": abs(fo.lnabsdet-fe.lnabsdet), "d_orbit_oracle": abs(fo.lnabsdet-o.lnabsdet), "d_eigen_oracle": abs(fe.lnabsdet-o.lnabsdet),
               "d_orbit_mine": abs(fo.lnabsdet-ln_m), "scale_over_L6": float(scale/L_BOX**6), "cond_Q": float(np.max(np.abs(ev_o))/np.min(np.abs(ev_o)))}
        rows.append(row)
        print(irrep, E, "N", o.N, "Np", o.Nproj, "orb-eig %.1e orb-orc %.1e eig-orc %.1e orb-mine %.1e" % (row["d_orbit_eigen"], row["d_orbit_oracle"], row["d_eigen_oracle"], row["d_orbit_mine"]),
              "scale/L^6 %.6f" % row["scale_over_L6"], "cond %.1e" % row["cond_Q"], flush=True)
json.dump({"commit": commit, "c": list(C_BASELINE), "k2_pwave_form": "coded_v2", "rows": rows}, open("/home/user/3body-qc-v4/results/tester_009_lndet.json", "w"), indent=1)

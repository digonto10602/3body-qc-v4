"""111_A2 at 0.3382011860823352 (kappa(Q)=1.5e13): who is right, fast or oracle?  High-precision (mpmath) Q from the oracle's own F,G,K2inv,K3 and V."""
import sys, json
import numpy as np, mpmath as mp
sys.path.insert(0, "/home/user/3body-qc-v4/src/python"); sys.path.insert(0, "/home/user/3body-qc-v4/tests/acceptance/m1")
from qc3 import fast, oracle
from m1_helpers import sector, C_BASELINE
mp.mp.dps = 60
s = sector("111_A2"); E = 0.3382011860823352
o = oracle.evaluate_point(s, E, C_BASELINE)
f = fast.evaluate_point(s, E, C_BASELINE)
fe = fast.evaluate_point(s, E, C_BASELINE, options={"projector": "eigen"})
print("N", o.N, f.N, "Nproj", o.Nproj, f.Nproj, "lnabs oracle %.10f fast-orbit %.10f fast-eigen %.10f" % (o.lnabsdet, f.lnabsdet, fe.lnabsdet))
def hp_ln(F, G, K2, K3, V, L6):
    F = mp.matrix(np.asarray(F).real.tolist()); G = mp.matrix(np.asarray(G).real.tolist()); K2 = mp.matrix(np.asarray(K2).real.tolist()); K3m = mp.matrix(np.asarray(K3).real.tolist())
    Vm = mp.matrix(np.asarray(V).real.tolist())
    H = K2 + F + G
    F3 = F / 3 - F * (mp.inverse(H) * F)
    X = mp.inverse(F3) * Vm
    Q = (Vm.T * (X + K3m * Vm)) / L6
    return mp.log(abs(mp.det(Q)))
L6 = s.L ** 6
res = {}
# (a) oracle matrices, oracle V ; (b) oracle matrices, orbit V ; (c) fast matrices, orbit V
res["hp_oracle_mats_oracleV"] = float(hp_ln(o.F, o.G, o.K2inv, o.K3, o.Vsel, L6))
res["hp_oracle_mats_orbitV"] = float(hp_ln(o.F, o.G, o.K2inv, o.K3, f.Vsel, L6))
res["hp_fast_mats_orbitV"] = float(hp_ln(f.F, f.G, f.K2inv, f.K3, f.Vsel, L6))
res["double_oracle"] = o.lnabsdet; res["double_fast_orbit"] = f.lnabsdet; res["double_fast_eigen"] = fe.lnabsdet
res["kappa_Q_oracle"] = float(np.linalg.cond(np.asarray(o.Q))); res["kappa_Q_fast"] = float(np.linalg.cond(np.asarray(f.Q)))
res["F_rel_fast_vs_oracle"] = float(np.linalg.norm(np.asarray(f.F) - np.asarray(o.F)) / np.linalg.norm(np.asarray(o.F)))
res["F3_rel_fast_vs_oracle"] = float(np.linalg.norm(np.asarray(f.F3) - np.asarray(o.F3)) / np.linalg.norm(np.asarray(o.F3)))
for k, v in res.items(): print(k, v)
print("fast-orbit - HP: %.3e   oracle - HP: %.3e   fast-eigen - HP: %.3e" % (res["double_fast_orbit"] - res["hp_fast_mats_orbitV"], res["double_oracle"] - res["hp_oracle_mats_oracleV"], res["double_fast_eigen"] - res["hp_fast_mats_orbitV"]))
json.dump(res, open("/home/user/3body-qc-v4/results/tester_010_illcond.json", "w"), indent=1)

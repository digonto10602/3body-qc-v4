"""Which G fill mode is closer to the exact ln|det Q| (60-digit) of ITS OWN matrices, and how far apart are the two exact values?"""
import sys, json
import numpy as np, mpmath as mp
sys.path.insert(0, "/home/user/3body-qc-v4/src/python"); sys.path.insert(0, "/home/user/3body-qc-v4/tests/acceptance/m1")
from qc3 import fast, oracle, _core
from qc3.fast import _params
from m1_helpers import sector, C_BASELINE
mp.mp.dps = 60
def hp(F, G, K2, K3, V, L6):
    M = lambda a: mp.matrix(np.asarray(a).real.tolist())
    F, G, K2, K3, V = M(F), M(G), M(K2), M(K3), M(V)
    F3 = F / 3 - F * (mp.inverse(K2 + F + G) * F)
    Q = (V.T * (mp.inverse(F3) * V + K3 * V)) / L6
    return mp.log(abs(mp.det(Q)))
out = []
for irrep, E in (("100_A2", 0.291556099879), ("200_A2", 0.34738), ("111_A2", 0.32036), ("100_A2", 0.291798362)):
    s = sector(irrep); L6 = s.L ** 6
    o = oracle.evaluate_point(s, E, C_BASELINE)
    fo = fast.evaluate_point(s, E, C_BASELINE)
    d = _core.fast_point(_params(s, E, C_BASELINE, True), 1e-16, "orbit", ["G", "K3", "Vsel", "F", "K2inv"], True)
    r = dict(irrep=irrep, Ecm=E, N=o.N,
             hp_ordered=float(hp(fo.F, fo.G, fo.K2inv, fo.K3, fo.Vsel, L6)), hp_symmetric=float(hp(d["F"], d["G"], d["K2inv"], d["K3"], d["Vsel"], L6)),
             double_ordered=fo.lnabsdet, double_symmetric=float(d["lnabsdet"]), double_oracle=o.lnabsdet)
    r["exact_diff_sym_minus_ord"] = r["hp_symmetric"] - r["hp_ordered"]
    r["err_ordered_vs_own_hp"] = r["double_ordered"] - r["hp_ordered"]; r["err_symmetric_vs_own_hp"] = r["double_symmetric"] - r["hp_symmetric"]
    r["err_oracle_vs_ordered_hp"] = r["double_oracle"] - r["hp_ordered"]
    out.append(r); print(r, flush=True)
json.dump(out, open("/home/user/3body-qc-v4/results/tester_010_gmode_hp.json", "w"), indent=1)

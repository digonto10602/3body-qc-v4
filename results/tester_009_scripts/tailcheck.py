import sys, numpy as np
sys.path.insert(0, "/home/user/3body-qc-v4/src/python"); sys.path.insert(0, "/home/user/3body-qc-v4/tests/acceptance/m1")
from qc3 import fast
from qc3.qc import slogdet
from m1_helpers import sector, anchor_points, C_BASELINE, as_c
for p in anchor_points(small_only=True):
    s = sector(p["irrep"])
    full = fast.evaluate_point(s, p["Ecm"], C_BASELINE, want=("F","G","K2inv","F3","F3inv","Vsel","K3","Q"), options={"projector": "orbit"})
    sub = fast.evaluate_point(s, p["Ecm"], C_BASELINE, want=("Q",), options={"projector": "orbit"})
    s_np, l_np = slogdet(as_c(full.Q))
    print(p["irrep"], p["Ecm"], "tail:", full.metadata["tail"], sub.metadata["tail"], "| F3inv None for subset:", sub.F3inv is None and full.F3inv is not None,
          "| sign==np:", int(full.sign) == s_np, "dln(full,np) %.1e" % abs(full.lnabsdet-l_np), "dln(sub,full) %.1e" % abs(sub.lnabsdet-full.lnabsdet), "sign sub==full", int(sub.sign)==int(full.sign),
          "relQ %.1e" % (np.linalg.norm(as_c(sub.Q)-as_c(full.Q))/np.linalg.norm(as_c(full.Q))))

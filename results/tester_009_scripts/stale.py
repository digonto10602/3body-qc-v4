import sys, json, numpy as np, subprocess
sys.path.insert(0, "/home/user/3body-qc-v4/src/python"); sys.path.insert(0, "/home/user/3body-qc-v4/tests/acceptance/m1")
from qc3 import _core, fast
from qc3.fast import _params
from qc3.sector import Sector
rng = np.random.default_rng(20261022)
bad = 0; n = 0
base = dict(Lbyas=20, xi=3.444, atmK=0.09698, atmpi=0.06906, epsilon_H=0.0, max_shell_num=20, eta_1=1.0, eta_2=0.5, alpha=0.5, eig_tol=0.05)
variants = []
for i in range(60):
    v = dict(base)
    key = ["Lbyas", "xi", "atmK", "atmpi", "epsilon_H", "max_shell_num", "alpha", "eta_2", "waves"][i % 9]
    if key == "Lbyas": v["Lbyas"] = float(rng.choice([16, 18, 20, 24]))
    elif key == "xi": v["xi"] = float(rng.uniform(3.0, 4.0))
    elif key == "atmK": v["atmK"] = float(rng.uniform(0.09, 0.10))
    elif key == "atmpi": v["atmpi"] = float(rng.uniform(0.06, 0.075))
    elif key == "epsilon_H": v["epsilon_H"] = float(rng.uniform(0, 0.3))
    elif key == "max_shell_num": v["max_shell_num"] = int(rng.choice([8, 12, 16, 20]))
    elif key == "alpha": v["alpha"] = float(rng.uniform(0.3, 0.7))
    elif key == "eta_2": v["eta_2"] = float(rng.uniform(0.3, 0.7))
    variants.append((key, v))
for irrep in ("100_A2", "111_A2"):
    for key, v in variants:
        for wv in ((0, 1), (0,)) if key == "waves" else ((0, 1),):
            s = Sector.production(irrep, k2_pwave_form="coded_v2", **{k: x for k, x in v.items()})
            p = _params(s, 0.31, (0., 0., 0., 0.), False)
            if key == "waves": p["waves_1"] = list(wv)
            r = _core.compare_kinematics(p); n += 1
            if not r["identical"]: bad += 1; print("STALE/MISMATCH", irrep, key, v)
print("cases", n, "mismatches", bad, "cache", fast.cache_info())

import sys, json, subprocess, pickle, time
import numpy as np
which = sys.argv[1]  # 'dbl' or 'ld'
pkg = "/home/user/3body-qc-v4/src/python" if which == "dbl" else "/tmp/claude-0/-home-user-3body-qc-v4/756827e6-de12-5b97-9946-026a972abb7d/scratchpad/ld"
sys.path.insert(0, pkg); sys.path.insert(0, "/home/user/3body-qc-v4/tests/acceptance/m1")
import qc3; assert (which == "dbl") == ("scratchpad" not in qc3.__file__)
from qc3 import fast
from m1_helpers import sector, load_anchors, IRREPS, C_BASELINE
pts = [(p["irrep"], p["Ecm"]) for p in load_anchors()["points"]]
pts += [("200_A2",0.346802591729),("111_A2",0.337572574251),("111_A2",0.330752399270),("111_A2",0.337780731301)]
for r in json.load(open("/home/user/3body-qc-v4/results/m1_step6_gamma_stress.json")): pts.append((r["irrep"], r["Ecm"]))
out = {}
for irrep, E in pts:
    s = sector(irrep)
    for mode in (True, False):
        r = fast.evaluate_point(s, E, C_BASELINE, want=("F", "Q"), options={"use_classes": mode})
        ts = []
        for _ in range(11):
            t = fast.evaluate_point(s, E, C_BASELINE, want=("Q",), options={"use_classes": mode}).timings
            ts.append(t["F"])
        out[(irrep, E, mode)] = (np.asarray(r.F), r.lnabsdet, r.sign, 1e3*float(np.median(ts)))
pickle.dump(out, open(f"{sys.argv[2]}", "wb"))

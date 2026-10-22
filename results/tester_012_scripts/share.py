import sys, json
import numpy as np
sys.path.insert(0, "/home/user/3body-qc-v4/src/python"); sys.path.insert(0, "/home/user/3body-qc-v4/tests/acceptance/m1")
from qc3 import fast, _core
from qc3.fast import _params
from m1_helpers import sector, load_anchors, C_BASELINE
rows = []
for p in load_anchors()["points"]:
    s = sector(p["irrep"]); E = p["Ecm"]; P = np.array(s.nnP)
    r = fast.evaluate_point(s, E, C_BASELINE, want=("Q",))
    keys = set(); nblocks = 0; seen = set()
    for x in r.basis:
        blk = (int(x['flavour']), int(x['nx']), int(x['ny']), int(x['nz']))
        if blk in seen: continue
        seen.add(blk); nblocks += 1
        nP = P - np.array(blk[1:]); c = tuple(int(v) for v in sorted(np.abs(nP), reverse=True))
        keys.add((blk[0], c, blk[1]**2 + blk[2]**2 + blk[3]**2))
    d = _core.fast_point(_params(s, E, C_BASELINE, False), 1e-16, "orbit", ["Q"])["stats"]
    rows.append(dict(irrep=p["irrep"], Ecm=E, blocks_basis=nblocks, f_blocks=int(d["f_blocks"]), distinct_flavour_c_p2=len(keys), f_sums=int(d["f_sums"])))
    print(rows[-1])
print("f_sums == distinct (flavour, canonical c, |p|^2) at all anchors:", all(r["f_sums"] == r["distinct_flavour_c_p2"] for r in rows))
json.dump(rows, open("/home/user/3body-qc-v4/results/tester_012_share_counts.json", "w"), indent=1)

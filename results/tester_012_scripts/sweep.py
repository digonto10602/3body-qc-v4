import sys, json
import numpy as np
sys.path.insert(0, "/home/user/3body-qc-v4/src/python"); sys.path.insert(0, "/home/user/3body-qc-v4/tests/acceptance/m1")
from qc3 import fast
from m1_helpers import sector, IRREPS, C_BASELINE
rng = np.random.default_rng(20261022)
N = int(sys.argv[1]); rows = []; bad = 0
worst = {}
for irrep in IRREPS:
    s = sector(irrep)
    E = np.concatenate([rng.uniform(0.26310, 0.36000, N), np.asarray(fast.dimension_jumps(s, 0.2631, 0.36)) + rng.choice([-1, 1], 11 if irrep!="000_A1m" else 11)[:0].sum()])
    jumps = np.asarray(fast.dimension_jumps(s, 0.2631, 0.36))
    E = np.concatenate([rng.uniform(0.26310, 0.36000, N), jumps + 3e-9, jumps - 3e-9])
    nterm_bad = 0; mxF = 0; mxd = 0
    for e in E:
        fc = fast.evaluate_point(s, float(e), C_BASELINE, want=("F", "Q"))
        fd = fast.evaluate_point(s, float(e), C_BASELINE, want=("F", "Q"), options={"use_classes": False})
        r = np.linalg.norm(fc.F - fd.F) / np.linalg.norm(fd.F)
        ok = fc.metadata["lattice_terms"] == fd.metadata["lattice_terms"]
        nterm_bad += (not ok); mxF = max(mxF, r); mxd = max(mxd, abs(fc.lnabsdet - fd.lnabsdet) if fc.sign == fd.sign else 99)
        rows.append(dict(irrep=irrep, Ecm=float(e), rel=float(r), terms_equal=bool(ok), N=fc.N))
    worst[irrep] = dict(n=len(E), term_mismatch=int(nterm_bad), max_F_rel=mxF, max_dln=mxd)
    print(irrep, worst[irrep], flush=True)
json.dump(dict(seed=20261022, c=list(C_BASELINE), n_per_irrep_random=N, summary=worst, rows=rows), open("/home/user/3body-qc-v4/results/tester_012_sweep.json", "w"))

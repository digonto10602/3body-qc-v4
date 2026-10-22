import sys, json, statistics as st
sys.path.insert(0, "/home/user/3body-qc-v4/src/python"); sys.path.insert(0, "/home/user/3body-qc-v4/tests/acceptance/m1")
from qc3 import fast, _core
from qc3.fast import _params
from m1_helpers import sector, C_BASELINE
out = {}
for irrep, E in (("111_A2", 0.35103), ("000_A1m", 0.34538), ("100_A2", 0.2915561), ("200_A2", 0.34738)):
    par = _params(sector(irrep), E, C_BASELINE, False)
    g = {False: [], True: []}; k3 = []
    for sym in (False, True): _core.fast_point(par, 1e-16, "orbit", ["Q"], sym)
    for _ in range(41):
        for sym in (False, True):
            d = _core.fast_point(par, 1e-16, "orbit", ["Q"], sym); g[sym].append(d["timings"]["G"])
    N = d["N"]
    out["%s_%s" % (irrep, E)] = dict(N=N, G_ordered_ms=1e3 * st.median(g[False]), G_symmetric_ms=1e3 * st.median(g[True]), diff_ms=1e3 * (st.median(g[False]) - st.median(g[True])))
    print(irrep, E, out["%s_%s" % (irrep, E)], flush=True)
json.dump(out, open("/home/user/3body-qc-v4/results/tester_010_gmodes.json", "w"), indent=1)

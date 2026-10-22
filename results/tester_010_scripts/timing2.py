"""Speed-up vs oracle at 111_A2 0.35103 (taskset -c 0, BLAS not pinned) + G fill-mode cost.  Medians of 7."""
import sys, json, time, os, statistics as st, subprocess
import numpy as np
sys.path.insert(0, "/home/user/3body-qc-v4/src/python"); sys.path.insert(0, "/home/user/3body-qc-v4/tests/acceptance/m1")
from qc3 import fast, oracle, _core
from qc3.fast import _params
from m1_helpers import sector, C_BASELINE
commit = subprocess.check_output(["git", "-C", "/home/user/3body-qc-v4", "rev-parse", "--short", "HEAD"]).decode().strip()
s = sector("111_A2"); E = 0.35103
fast.evaluate_point(s, E, C_BASELINE, want=("Q",))
fw, ft = [], []
for _ in range(15):
    t0 = time.perf_counter(); r = fast.evaluate_point(s, E, C_BASELINE, want=("Q",)); fw.append(time.perf_counter() - t0); ft.append(r.timings["total"])
ow, op, ot = [], [], []
for _ in range(7):
    t0 = time.perf_counter(); o = oracle.evaluate_point(s, E, C_BASELINE); ow.append(time.perf_counter() - t0); op.append(float(o.timings.get("projector", 0.0))); ot.append(o.timings.get("total", 0.0))
res = dict(commit=commit, affinity=sorted(os.sched_getaffinity(0)), OPENBLAS=os.environ.get("OPENBLAS_NUM_THREADS"), OMP=os.environ.get("OMP_NUM_THREADS"),
           fast_wall_median_ms=1e3 * st.median(fw), fast_wall_min_ms=1e3 * min(fw), fast_cpp_total_median_ms=1e3 * st.median(ft), fast_wall_all_ms=[1e3 * x for x in fw],
           oracle_wall_median_s=st.median(ow), oracle_wall_all_s=ow, oracle_projector_median_s=st.median(op), oracle_total_median_s=st.median(ot))
res["speedup_with_oracle_projector"] = res["oracle_wall_median_s"] / (res["fast_wall_median_ms"] * 1e-3)
res["speedup_without_oracle_projector"] = (res["oracle_wall_median_s"] - res["oracle_projector_median_s"]) / (res["fast_wall_median_ms"] * 1e-3)
res["speedup_cpp_total_with"] = res["oracle_wall_median_s"] / (res["fast_cpp_total_median_ms"] * 1e-3)
# G fill modes through the core: stage timings, N=160
par = _params(s, E, C_BASELINE, False)
modes = {}
for sym in (False, True):
    gs, tt, lns = [], [], []
    for _ in range(15):
        d = _core.fast_point(par, 1e-16, "orbit", ["Q"], sym); gs.append(d["timings"]["G"]); tt.append(d["timings"]["total"]); lns.append(d["lnabsdet"])
    modes["symmetric" if sym else "ordered"] = dict(G_median_ms=1e3 * st.median(gs), total_median_ms=1e3 * st.median(tt), lnabsdet=lns[-1])
res["g_modes_N160"] = modes
# K3 stage, F stage etc. at N=160
stg = {}
rs = [fast.evaluate_point(s, E, C_BASELINE, want=("Q",)).timings for _ in range(15)]
for k in rs[0]: stg[k] = 1e3 * st.median(r[k] for r in rs)
res["stage_ms_N160"] = stg
print(json.dumps(res, indent=1))
json.dump(res, open("/home/user/3body-qc-v4/results/tester_010_speedup.json", "w"), indent=1)

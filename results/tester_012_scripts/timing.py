import sys, json, time, os, subprocess, statistics as st
import numpy as np
sys.path.insert(0, "/home/user/3body-qc-v4/src/python"); sys.path.insert(0, "/home/user/3body-qc-v4/tests/acceptance/m1")
from qc3 import fast, _core
from m1_helpers import sector, load_anchors, C_BASELINE
commit = subprocess.check_output(["git", "-C", "/home/user/3body-qc-v4", "rev-parse", "--short", "HEAD"]).decode().strip()
NW, NC = 21, 9
res = []
for p in load_anchors()["points"]:
    s = sector(p["irrep"]); E = p["Ecm"]; row = dict(irrep=p["irrep"], Ecm=E, N=p["N"], Nproj=p["Nproj"])
    for mode in (True, False):
        opt = {"use_classes": mode}
        cold_w, cold_t = [], []
        for _ in range(NC):
            fast.clear_caches()
            t0 = time.perf_counter(); r = fast.evaluate_point(s, E, C_BASELINE, want=("Q",), options=opt); cold_w.append(time.perf_counter() - t0); cold_t.append(r.timings)
        w, tt = [], []
        fast.evaluate_point(s, E, C_BASELINE, want=("Q",), options=opt)
        for _ in range(NW):
            t0 = time.perf_counter(); r = fast.evaluate_point(s, E, C_BASELINE, want=("Q",), options=opt); w.append(time.perf_counter() - t0); tt.append(r.timings)
        key = "classes" if mode else "direct"
        row[key] = dict(warm_wall_median_ms=1e3 * st.median(w), warm_wall_min_ms=1e3 * min(w), warm_total_median_ms=1e3 * st.median(t["total"] for t in tt),
                        warm_F_median_ms=1e3 * st.median(t["F"] for t in tt), warm_F_min_ms=1e3 * min(t["F"] for t in tt),
                        warm_stage_median_ms={k: 1e3 * st.median(t[k] for t in tt) for k in tt[0]},
                        cold_wall_median_ms=1e3 * st.median(cold_w), cold_wall_min_ms=1e3 * min(cold_w),
                        cold_F_median_ms=1e3 * st.median(t["F"] for t in cold_t),
                        cold_minus_warm_median_ms=1e3 * (st.median(cold_w) - st.median(w)), cold_minus_warm_min_ms=1e3 * (min(cold_w) - min(w)))
    print("%s %.7f N=%3d  classes warm wall med %.2f min %.2f (cpp %.2f, F %.3f)  direct warm wall med %.2f (F %.2f)  speedup med %.2fx min %.2fx | cold classes %.1f ms cold-warm %.1f/%.1f ms" % (
        row["irrep"], E, row["N"], row["classes"]["warm_wall_median_ms"], row["classes"]["warm_wall_min_ms"], row["classes"]["warm_total_median_ms"], row["classes"]["warm_F_median_ms"],
        row["direct"]["warm_wall_median_ms"], row["direct"]["warm_F_median_ms"], row["direct"]["warm_wall_median_ms"] / row["classes"]["warm_wall_median_ms"],
        row["direct"]["warm_wall_min_ms"] / row["classes"]["warm_wall_min_ms"], row["classes"]["cold_wall_median_ms"], row["classes"]["cold_minus_warm_median_ms"], row["classes"]["cold_minus_warm_min_ms"]), flush=True)
    res.append(row)
json.dump(dict(commit=commit, affinity=sorted(os.sched_getaffinity(0)), OMP=os.environ.get("OMP_NUM_THREADS"), OPENBLAS=os.environ.get("OPENBLAS_NUM_THREADS"), n_warm=NW, n_cold=NC, c=list(C_BASELINE),
               k2_pwave_form="coded_v2 (m1_helpers.sector)", rows=res), open("/home/user/3body-qc-v4/results/tester_012_timing.json", "w"), indent=1)

import sys, json, subprocess, time, statistics as st, numpy as np
sys.path.insert(0, "/home/user/3body-qc-v4/src/python"); sys.path.insert(0, "/home/user/3body-qc-v4/tests/acceptance/m1")
from qc3 import fast, _core
from qc3.fast import _params
from m1_helpers import sector, load_anchors, IRREPS, C_BASELINE, C_ZERO
commit = subprocess.check_output(["git", "-C", "/home/user/3body-qc-v4", "rev-parse", "--short", "HEAD"]).decode().strip()
out = {"commit": commit, "cpu_affinity": sorted(__import__("os").sched_getaffinity(0)), "c": list(C_BASELINE)}
KEYS = ["basis", "projector", "F", "G", "K2inv", "K3", "projection", "lu", "total"]
def ev(s, E, c=C_BASELINE, opt=None):
    t0 = time.perf_counter(); r = fast.evaluate_point(s, E, c, want=("Q",), options=opt); w = time.perf_counter() - t0
    return r, w
# ---- 1. 12 anchors: cold (clear_caches) and warm
pts = load_anchors()["points"]
rows = []
for p in pts:
    s = sector(p["irrep"]); E = p["Ecm"]
    fast.clear_caches()
    rc, wc = ev(s, E)
    cold_t = dict(rc.timings); assert rc.metadata["projector_cache_hit"] is False
    ws = []; tt = []
    for _ in range(7):
        rw, w = ev(s, E); ws.append(w); tt.append(dict(rw.timings))
    med = {k: st.median(t.get(k, 0.0) for t in tt) for k in tt[0]}
    row = {"irrep": p["irrep"], "Ecm": E, "N": rc.N, "Nproj": rc.Nproj, "wall_cold_ms": 1e3*wc, "wall_warm_ms": 1e3*st.median(ws),
           "cold_stage_ms": {k: 1e3*v for k, v in cold_t.items()}, "warm_stage_ms": {k: 1e3*v for k, v in med.items()}}
    row["cold_minus_warm_ms"] = row["wall_cold_ms"] - row["wall_warm_ms"]
    row["cold_minus_warm_total_ms"] = 1e3*(cold_t["total"] - med["total"])
    rows.append(row)
    print(f'{p["irrep"]} {E} N={rc.N} Np={rc.Nproj} cold {1e3*wc:7.1f} warm {1e3*st.median(ws):6.1f} cold-warm {row["cold_minus_warm_ms"]:6.1f} | cold basis {1e3*cold_t["basis"]:.2f} proj {1e3*cold_t["projector"]:.3f} | warm basis {1e3*med["basis"]:.3f} proj {1e3*med["projector"]:.3f} F {1e3*med["F"]:.1f} G {1e3*med["G"]:.1f} K3 {1e3*med["K3"]:.1f} tot {1e3*med["total"]:.1f}', flush=True)
out["anchors"] = rows
# ---- 2. cold projector at all production bases
bases = []
for irrep in IRREPS:
    s = sector(irrep)
    j = np.asarray(fast.dimension_jumps(s, 0.2631, 0.36)); edges = np.concatenate([[0.2631], j, [0.36]]); mids = 0.5*(edges[1:]+edges[:-1])
    for E in mids:
        res = []
        for rep in range(5):
            fast.clear_caches(); r, _ = ev(s, float(E), C_ZERO)
            res.append((r.timings["projector"], r.timings["basis"]))
        r_warm, _ = ev(s, float(E), C_ZERO)
        bases.append({"irrep": irrep, "Ecm": float(E), "N": r.N, "Nproj": r.Nproj, "cold_projector_ms_min": 1e3*min(a for a, b in res), "cold_projector_ms_max": 1e3*max(a for a, b in res),
                      "cold_basis_ms_median": 1e3*st.median(b for a, b in res), "warm_basis_ms": 1e3*r_warm.timings["basis"]})
mx = max(bases, key=lambda b: b["cold_projector_ms_max"])
print("bases", len(bases), "max cold projector (ms) over 58 bases x5 repeats:", mx)
out["bases"] = bases
# ---- 3. builder: threshold table build vs per-point basis; the 111_A2 0.35103 case
s = sector("111_A2"); E = 0.35103
tab = []
for rep in range(7):
    fast.clear_caches()
    t0 = time.perf_counter(); fast.dimension_jumps(s, 0.2631, 0.36); t_tab = time.perf_counter() - t0   # builds both tables (cold)
    ci = fast.cache_info()
    t0 = time.perf_counter(); r1, w1 = ev(s, E); t1 = time.perf_counter() - t0                          # tables present, projector cold
    fast.clear_caches()
    r2, w2 = ev(s, E)                                                                                   # everything cold
    r3, w3 = ev(s, E)                                                                                   # warm
    tab.append({"table_build_ms_both_channels": 1e3*t_tab, "threshold_entries": ci["threshold_entries"], "projcold_basis_ms": 1e3*r1.timings["basis"], "projcold_proj_ms": 1e3*r1.timings["projector"], "projcold_total_ms": 1e3*r1.timings["total"],
                "allcold_basis_ms": 1e3*r2.timings["basis"], "allcold_proj_ms": 1e3*r2.timings["projector"], "allcold_total_ms": 1e3*r2.timings["total"],
                "warm_basis_ms": 1e3*r3.timings["basis"], "warm_total_ms": 1e3*r3.timings["total"], "allcold_wall_ms": 1e3*w2, "warm_wall_ms": 1e3*w3})
med = {k: st.median(t[k] for t in tab) for k in tab[0]}
print("111_A2 0.35103 builder breakdown (median of 7):", json.dumps(med, indent=1))
out["builder_111_A2_0.35103"] = {"reps": tab, "median": med}
json.dump(out, open("/home/user/3body-qc-v4/results/tester_009_timing.json", "w"), indent=1)

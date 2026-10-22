import sys, time, json, statistics as st
sys.path.insert(0, "/home/user/3body-qc-v4/src/python"); sys.path.insert(0, "/home/user/3body-qc-v4/tests/acceptance/m1")
from qc3 import fast
from m1_helpers import sector, C_BASELINE
s = sector("111_A2"); E = 0.35103
keys = ["basis", "projector", "F", "G", "K2inv", "K3", "projection", "lu", "total"]
def ev():
    r = fast.evaluate_point(s, E, C_BASELINE, want=("Q",)); return {k: 1e3*r.timings.get(k, 0.0) for k in keys}
ev()
cold, warm = [], []
for rep in range(30):
    fast.clear_caches(); cold.append(ev())
    w = [ev() for _ in range(3)]; warm.append({k: st.median(x[k] for x in w) for k in keys})
out = {}
print("stage: median cold, median warm, median (cold-warm) paired, IQR of paired diff")
for k in keys:
    d = sorted(c[k]-w[k] for c, w in zip(cold, warm))
    out[k] = {"cold_med": st.median(c[k] for c in cold), "warm_med": st.median(w[k] for w in warm), "diff_med": st.median(d), "diff_q1": d[7], "diff_q3": d[22], "diff_max": d[-1]}
    print(f"{k:11s} cold {out[k]['cold_med']:7.2f} warm {out[k]['warm_med']:7.2f} diff med {out[k]['diff_med']:6.2f} [{d[7]:6.2f},{d[22]:6.2f}] max {d[-1]:6.2f}")
json.dump({"reps": 30, "stages": out}, open("/home/user/3body-qc-v4/results/tester_009_coldwarm_111A2.json", "w"), indent=1)

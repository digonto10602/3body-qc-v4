import sys, time, json
sys.path.insert(0, "/home/user/3body-qc-v4/src/python"); sys.path.insert(0, "/home/user/3body-qc-v4/tests/acceptance/m1")
t0 = time.perf_counter()
from qc3 import fast
from m1_helpers import sector, C_BASELINE
t_imp = time.perf_counter() - t0
s = sector("111_A2"); E = 0.35103
def ev():
    t = time.perf_counter(); r = fast.evaluate_point(s, E, C_BASELINE, want=("Q",)); w = time.perf_counter() - t
    return {"wall": 1e3*w, **{k: 1e3*v for k, v in r.timings.items() if k in ("basis", "projector", "F", "G", "K3", "total")}}
first = ev(); fast.clear_caches(); cold2 = ev(); warm = [ev() for _ in range(3)]
print(json.dumps({"import_ms": 1e3*t_imp, "first": first, "cold2": cold2, "warm": warm[-1]}))

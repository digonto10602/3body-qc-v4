import sys, time, os, json
sys.path.insert(0, "/home/user/3body-qc-v4/src/python"); sys.path.insert(0, "/home/user/3body-qc-v4/tests/acceptance/m1")
from qc3 import fast
from m1_helpers import sector, C_BASELINE
out = {}
for irrep, E in (("000_A1m", 0.26932), ("100_A2", 0.2915561), ("111_A2", 0.35103), ("000_A1m", 0.34538)):
    s = sector(irrep); fast.evaluate_point(s, E, C_BASELINE, want=("Q",))
    t0c, t0w = time.process_time(), time.perf_counter()
    for _ in range(20): fast.evaluate_point(s, E, C_BASELINE, want=("Q",))
    cpu, wall = time.process_time()-t0c, time.perf_counter()-t0w
    out[f"{irrep}_{E}"] = {"cpu": cpu, "wall": wall, "ratio": cpu/wall}
print(os.environ.get("OPENBLAS_NUM_THREADS"), json.dumps(out))

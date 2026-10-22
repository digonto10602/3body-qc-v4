"""CPU/wall of 20 warm calls (the locked-test protocol) for a given point, FIRST measured window of the process, BLAS unpinned."""
import sys, time, os, json
sys.path.insert(0, "/home/user/3body-qc-v4/src/python"); sys.path.insert(0, "/home/user/3body-qc-v4/tests/acceptance/m1")
irrep, E = sys.argv[1], float(sys.argv[2])
from qc3 import fast
from m1_helpers import sector, C_BASELINE
delay = float(sys.argv[3]) if len(sys.argv) > 3 else 0.0
s = sector(irrep); fast.evaluate_point(s, E, C_BASELINE, want=("Q",))
if delay: time.sleep(delay)
t0c, t0w = time.process_time(), time.perf_counter()
for _ in range(20): fast.evaluate_point(s, E, C_BASELINE, want=("Q",))
cpu, wall = time.process_time() - t0c, time.perf_counter() - t0w
nthr = len(os.listdir("/proc/self/task"))
print(json.dumps({"point": f"{irrep}_{E}", "cpu": cpu, "wall": wall, "ratio": cpu / wall, "threads_in_process": nthr}))

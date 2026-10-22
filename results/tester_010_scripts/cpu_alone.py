"""The locked test body, run in a fresh process, optionally after sleeping `delay` seconds (to let OpenBLAS idle)."""
import sys, time, os, json
sys.path.insert(0, "/home/user/3body-qc-v4/src/python"); sys.path.insert(0, "/home/user/3body-qc-v4/tests/acceptance/m1")
delay = float(sys.argv[1]) if len(sys.argv) > 1 else 0.0
from qc3 import fast
from m1_helpers import sector, C_BASELINE
s = sector("000_A1m"); E = 0.26932
fast.evaluate_point(s, E, C_BASELINE, want=("Q",))
if delay: time.sleep(delay)
t0c, t0w = time.process_time(), time.perf_counter()
for _ in range(20): fast.evaluate_point(s, E, C_BASELINE, want=("Q",))
cpu, wall = time.process_time() - t0c, time.perf_counter() - t0w
print(json.dumps({"delay": delay, "OPENBLAS": os.environ.get("OPENBLAS_NUM_THREADS"), "cpu": cpu, "wall": wall, "ratio": cpu / wall}))

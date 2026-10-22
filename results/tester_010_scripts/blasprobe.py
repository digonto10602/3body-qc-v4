"""Bare numpy probe: CPU/wall of a pure sleep right after `import numpy` (no qc3, no BLAS call)."""
import time, os, json, sys
t0c, t0w = time.process_time(), time.perf_counter()   # note: process_time already includes import cost; reset below
import numpy as np
t1c, t1w = time.process_time(), time.perf_counter()
res = {"import_cpu": t1c - t0c, "import_wall": t1w - t0w}
w = []
for k in range(6):
    a, b = time.process_time(), time.perf_counter()
    time.sleep(0.1)
    w.append((time.process_time() - a, time.perf_counter() - b))
res["sleep_windows_cpu_over_wall"] = [c / ww for c, ww in w]
res["OPENBLAS_NUM_THREADS"] = os.environ.get("OPENBLAS_NUM_THREADS")
print(json.dumps(res))

import sys, os
sys.path.insert(0, "/home/user/3body-qc-v4/src/python"); sys.path.insert(0, "/home/user/3body-qc-v4/tests/acceptance/m1")
from qc3 import fast
from m1_helpers import sector, C_BASELINE
s = sector("111_A2")
fast.evaluate_point(s, 0.35103, C_BASELINE, want=("Q",))   # imports done
try: open("/MARK_START_COLD")
except OSError: pass
fast.clear_caches()
fast.evaluate_point(s, 0.35103, C_BASELINE, want=("Q",))   # cold: builds F tables, projector, thresholds
try: open("/MARK_START_WARM")
except OSError: pass
for _ in range(5): fast.evaluate_point(s, 0.35103, C_BASELINE, want=("Q",))
try: open("/MARK_END")
except OSError: pass
print(fast.cache_info())

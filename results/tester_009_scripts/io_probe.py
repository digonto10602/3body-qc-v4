import sys, os
sys.path.insert(0, "/home/user/3body-qc-v4/src/python"); sys.path.insert(0, "/home/user/3body-qc-v4/tests/acceptance/m1")
from qc3 import fast
from m1_helpers import sector
s = sector("111_A2"); fast.evaluate_point(s, 0.30, (0,0,0,0), want=("Q",)); fast.clear_caches()   # import/lazy init done before the markers
os.write(2, b"\nMARK_BEGIN\n")
fast.dimension_jumps(s, 0.2631, 0.36)
r1 = fast.evaluate_point(s, 0.35103, (0,0,0,0), want=("Q",))
r2 = fast.evaluate_point(sector("000_A1m"), 0.34538, (0,0,0,0), want=("Q",))
r3 = fast.evaluate_point(s, 0.35103, (0,0,0,0), want=("Q",))
os.write(2, b"\nMARK_END\n")
print(fast.cache_info())

import sys, json, resource, os
import numpy as np
sys.path.insert(0, "/home/user/3body-qc-v4/src/python"); sys.path.insert(0, "/home/user/3body-qc-v4/tests/acceptance/m1")
from qc3 import fast
from m1_helpers import sector, IRREPS, C_BASELINE
def rss(): 
    for l in open('/proc/self/status'):
        if l.startswith('VmRSS'): return int(l.split()[1])/1024
cs = {}; per_irrep = {}
r0 = rss(); fast.clear_caches()
nbases = 0; stats = []
for irrep in IRREPS:
    s = sector(irrep); P = np.array(s.nnP)
    j = np.concatenate([[0.2631], np.asarray(fast.dimension_jumps(s, 0.2631, 0.36)), [0.36]])
    mids = 0.5 * (j[1:] + j[:-1])
    per_irrep[irrep] = len(mids)
    for E in mids:
        f = fast.evaluate_point(s, float(E), C_BASELINE, want=("Q",)); nbases += 1
        b = f.basis
        for row in set((int(x['flavour']), int(x['nx']), int(x['ny']), int(x['nz'])) for x in b):
            nP = P - np.array(row[1:]); c = tuple(sorted(np.abs(nP), reverse=True))
            cs.setdefault(tuple(int(v) for v in c), set()).add(irrep)
info = fast.cache_info(); r1 = rss()
print("bases", nbases, per_irrep); print("f_table_entries", info["f_table_entries"], "distinct canonical c from basis enumeration", len(cs))
print("RSS MB before/after", r0, r1)
print(sorted(cs))
json.dump(dict(bases=nbases, per_irrep=per_irrep, f_table_entries=info["f_table_entries"], c_list=sorted(cs), rss_before_MB=r0, rss_after_MB=r1), open("/home/user/3body-qc-v4/results/tester_012_tables.json", "w"))

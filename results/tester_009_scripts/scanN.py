"""Independent scan of the basis dimension N(E) with the ORACLE's config_maker (via compare_kinematics hook) by bisection
on the 20000-point production grid, then refinement of each jump to double precision in E_cm.  Compare to dimension_jumps
and to the threshold formula (own NumPy)."""
import sys, json, subprocess, numpy as np
sys.path.insert(0, "/home/user/3body-qc-v4/src/python")
sys.path.insert(0, "/home/user/3body-qc-v4/tests/acceptance/m1")
from qc3 import _core, fast
from qc3.fast import _params
from m1_helpers import sector, IRREPS, NNP, PRODUCTION, L_BOX
commit = subprocess.check_output(["git", "-C", "/home/user/3body-qc-v4", "rev-parse", "--short", "HEAD"]).decode().strip()
E0, DE, NG = 0.26310, 0.0969 / 19999, 20000
def grid(i): return 0.26310 + i * 0.0969 / 19999
out = {"commit": commit, "grid": "E_i=0.26310+i*0.0969/19999, i<20000", "irreps": {}}
def thresholds(irrep):
    mK, mpi = PRODUCTION["atmK"], PRODUCTION["atmpi"]
    P = 2*np.pi/L_BOX*np.array(NNP[irrep], float)
    g = np.arange(-20, 21); nx, ny, nz = np.meshgrid(g, g, g, indexing="ij")
    k = 2*np.pi/L_BOX*np.stack([nx.ravel(), ny.ravel(), nz.ravel()], 1).astype(float)
    res = []
    for (mi, mj, mk) in ((mK, mK, mpi), (mpi, mK, mK)):
        sig = abs(mj*mj-mk*mk) + 2*0.02*(mj+mk)*min(mj, mk)
        res.append(np.sqrt((k*k).sum(1)+mi*mi) + np.sqrt(sig + ((P-k)**2).sum(1)))
    return np.concatenate(res), float(P@P)
for irrep in IRREPS:
    s = sector(irrep)
    calls = [0]
    def N(E, which="oracle"):
        calls[0] += 1
        r = _core.compare_kinematics(_params(s, float(E), (0., 0., 0., 0.), False))
        assert r["identical"], (irrep, E)
        return r["N"]
    cache = {}
    def Ng(i):
        if i not in cache: cache[i] = N(grid(i))
        return cache[i]
    cells = []
    def rec(a, b):
        if Ng(a) == Ng(b): return
        if b - a == 1: cells.append((a, b)); return
        m = (a + b) // 2
        rec(a, m); rec(m, b)
    rec(0, NG - 1)
    # jumps in a cell may be several at the same cell: count N change
    bases = sorted(set(cache.values()))
    # refine
    jpos = []
    for (a, b) in cells:
        lo, hi = grid(a), grid(b); nlo, nhi = Ng(a), Ng(b)
        for _ in range(80):
            mid = 0.5*(lo+hi)
            if mid <= lo or mid >= hi: break
            if N(mid) == nlo: lo = mid
            else: hi = mid
        jpos.append({"cell": [a, b], "N_lo": nlo, "N_hi": nhi, "E_lo": lo, "E_hi": hi})
    # distinct bases along the grid: count distinct N at all grid points via cells
    nb = 1 + len(cells)  # if each cell has one jump
    dj = fast.dimension_jumps(s, 0.2631, 0.36)
    t, Psq = thresholds(irrep)
    tin = np.sort(t[(np.sqrt(t**2-Psq) > 0.2631) & (np.sqrt(t**2-Psq) < 0.36)])
    merged = []
    for v in tin:
        if not merged or abs(v-merged[-1]) > 1e-12*v: merged.append(v)
    rule = np.sqrt(np.array(merged)**2 - Psq)
    # match each scanned cell to dimension_jumps
    in_cell = [int(np.sum((dj > grid(a)) & (dj <= grid(b)))) for (a, b) in cells]
    # number of dj not in any cell
    allin = int(np.sum([(dj > grid(a)) & (dj <= grid(b)) for a, b in cells], axis=0).sum()) if cells else 0
    dmax = max(abs(j["E_hi"]-j["E_lo"]) for j in jpos)
    # nearest dj for each oracle jump
    dd = []
    for j in jpos:
        Ej = 0.5*(j["E_lo"]+j["E_hi"])
        k = int(np.argmin(np.abs(dj-Ej))); dd.append(abs(dj[k]-Ej))
    out["irreps"][irrep] = {
        "n_grid_cells_with_change": len(cells), "n_distinct_N_on_grid": len(bases), "n_dimension_jumps": int(len(dj)),
        "n_rule_jumps": int(len(rule)), "dj_per_cell": in_cell, "dj_in_cells_total": allin,
        "max_|dj-rule|/rule": float(np.max(np.abs(dj-rule)/rule)) if len(dj)==len(rule) else None,
        "max_|oracle_jump-dj|": float(max(dd)) if dd else 0.0, "max_bisect_bracket": float(dmax),
        "N_values_on_grid": bases, "oracle_N_calls": calls[0], "jumps_dj": [float(x) for x in dj],
        "jumps_oracle_bisect": [0.5*(j["E_lo"]+j["E_hi"]) for j in jpos], "cells": jpos,
    }
    print(irrep, "distinct N", len(bases), "cells", len(cells), "dj", len(dj), "rule", len(rule), "cell counts", set(in_cell),
          "max|dj-rule|/rule", out["irreps"][irrep]["max_|dj-rule|/rule"], "max|orc-dj|", max(dd) if dd else 0, "calls", calls[0], flush=True)
json.dump(out, open("/home/user/3body-qc-v4/results/tester_009_basis_scan.json", "w"), indent=1)

import sys, json, subprocess, numpy as np
sys.path.insert(0, "/home/user/3body-qc-v4/src/python"); sys.path.insert(0, "/home/user/3body-qc-v4/tests/acceptance/m1")
from qc3 import _core, fast
from qc3.fast import _params
from m1_helpers import sector, IRREPS
commit = subprocess.check_output(["git", "-C", "/home/user/3body-qc-v4", "rev-parse", "--short", "HEAD"]).decode().strip()
rng = np.random.default_rng(20261022)
res = {"commit": commit, "seed": 20261022, "irreps": {}}
tot = 0; bad = 0
for irrep in IRREPS:
    s = sector(irrep)
    jumps = fast.dimension_jumps(s, 0.2631, 0.36)
    E = list(rng.uniform(0.2631, 0.36, 300)) + [0.2631, 0.36]
    ulp = []
    for j in jumps:
        x = float(j)
        for k in range(-3, 4):
            y = x
            for _ in range(abs(k)): y = float(np.nextafter(y, np.inf if k > 0 else -np.inf))
            ulp.append(y)
    n_ident = 0
    jump_side = []
    for e in E + ulp:
        r = _core.compare_kinematics(_params(s, float(e), (0., 0., 0., 0.), False))
        tot += 1
        if not (r["identical"] and r["N"] == r["N_fast"]): bad += 1; print("MISMATCH", irrep, repr(e), r)
    # N change across the jump at +-1ulp: find the ulp at which N flips
    flips = []
    for j in jumps:
        x = float(j); Ns = []
        for k in range(-3, 4):
            y = x
            for _ in range(abs(k)): y = float(np.nextafter(y, np.inf if k > 0 else -np.inf))
            Ns.append(_core.compare_kinematics(_params(s, y, (0., 0., 0., 0.), False))["N"])
        flips.append(Ns)
    # location of the flip in ulps relative to dimension_jumps value
    offs = [int(np.argmax(np.diff(Ns) != 0)) - 2 if np.any(np.diff(Ns) != 0) else None for Ns in flips]
    res["irreps"][irrep] = {"n_jumps": len(jumps), "n_energies": len(E) + len(ulp), "flip_position_ulps": offs}
    print(irrep, "jumps", len(jumps), "energies", len(E) + len(ulp), "flip position (ulp offset of the last point below the flip; 0 = at dimension_jumps value):", sorted(set(offs)), flush=True)
res["total_energies"] = tot; res["mismatches"] = bad
print("total", tot, "mismatches", bad)
json.dump(res, open("/home/user/3body-qc-v4/results/tester_009_basis_bitident.json", "w"), indent=1)

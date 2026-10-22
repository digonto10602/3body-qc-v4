"""Tester 010: K3 basis, G, tail vs oracle at 12 anchors + random energies (seed 20261022). coded_v2 pinned."""
import sys, json, subprocess, time
import numpy as np
sys.path.insert(0, "/home/user/3body-qc-v4/src/python"); sys.path.insert(0, "/home/user/3body-qc-v4/tests/acceptance/m1")
from qc3 import fast, oracle, _core
from qc3.fast import _params
from m1_helpers import sector, load_anchors, IRREPS, C_BASELINE, C_ZERO
commit = subprocess.check_output(["git", "-C", "/home/user/3body-qc-v4", "rev-parse", "--short", "HEAD"]).decode().strip()
rng = np.random.default_rng(20261022)
NAMES = ("iso0", "iso1", "B", "E")
def rel(a, b):
    a = np.asarray(a, dtype=complex); b = np.asarray(b, dtype=complex)
    nb = np.linalg.norm(b); return float(np.linalg.norm(a - b) / nb) if nb > 0 else float(np.linalg.norm(a - b))
pts = [(p["irrep"], p["Ecm"], "anchor") for p in load_anchors()["points"]]
# random energies: 3 per irrep uniform + 2 adjacent to dimension jumps (jump +- 1e-9 ... in nearest grid-like offsets)
for irrep in IRREPS:
    s = sector(irrep)
    for E in rng.uniform(0.26310, 0.36000, 3): pts.append((irrep, float(E), "random"))
    j = np.asarray(fast.dimension_jumps(s, 0.2631, 0.36))
    for jj in rng.choice(j, 2, replace=False):
        eps = 10.0 ** rng.uniform(-9, -6)
        pts.append((irrep, float(jj + eps), "jump+%.1e" % eps)); pts.append((irrep, float(jj - eps), "jump-%.1e" % eps))
rows = []
for irrep, E, tag in pts:
    s = sector(irrep)
    o = oracle.evaluate_point(s, E, C_BASELINE)
    f = fast.evaluate_point(s, E, C_BASELINE)
    par = _params(s, E, C_BASELINE, False)
    Gs = np.asarray(_core.fast_point(par, 1e-16, "orbit", ["G", "Q"], True)["G"])
    fs = _core.fast_point(par, 1e-16, "orbit", ["Q"], True)
    row = dict(irrep=irrep, Ecm=E, tag=tag, N=o.N, Nproj=o.Nproj, N_fast=f.N, Nproj_fast=f.Nproj)
    for n in NAMES:
        row["K3_" + n] = rel(f.K3_basis[n], o.K3_basis[n])
        row["K3sym_" + n] = rel(np.asarray(o.K3_basis[n]), np.asarray(o.K3_basis[n]).T)
        row["K3fsym_" + n] = rel(np.asarray(f.K3_basis[n]), np.asarray(f.K3_basis[n]).T)
    row["K3"] = rel(f.K3, o.K3)
    row["G"] = rel(f.G, o.G); row["Gsym"] = rel(Gs, o.G); row["G_oracle_asym"] = rel(np.asarray(o.G).T, o.G)
    row["G_max_abs_diff"] = float(np.abs(np.asarray(f.G) - np.asarray(o.G)).max())
    row["Gs_exact_sym"] = bool(np.array_equal(Gs, Gs.T))
    row["lnabs_oracle"] = o.lnabsdet; row["lnabs_fast"] = f.lnabsdet; row["lnabs_fast_Gsym"] = float(fs["lnabsdet"])
    row["dln"] = abs(f.lnabsdet - o.lnabsdet); row["dln_Gsym"] = abs(float(fs["lnabsdet"]) - o.lnabsdet)
    row["sign_ok"] = (f.sign == o.sign); row["sign_Gsym_ok"] = (int(fs["sign"]) == o.sign)
    row["Q_rel"] = rel(f.Q, o.Q)
    # NumPy tail from fast's own F3inv, K3, V
    V = np.asarray(f.Vsel); Qn = V.conj().T @ (np.asarray(f.F3inv) + np.asarray(f.K3)) @ V / f.sector.L ** 6
    row["Q_vs_numpy_proj_inverse"] = rel(f.Q, Qn)
    sg, ld = np.linalg.slogdet(Qn); row["dln_vs_numpy_slogdet_own"] = abs(ld - f.lnabsdet); row["sign_numpy"] = int(sg.real)
    # kappa(Q)
    row["kappa_Q"] = float(np.linalg.cond(np.asarray(o.Q)))
    rows.append(row)
    print(f'{irrep} {E:.9f} {tag:14s} N={o.N:3d}/{f.N:3d} K3max={max(row["K3_"+n] for n in NAMES):.1e} K3={row["K3"]:.1e} G={row["G"]:.1e} Gsym={row["Gsym"]:.1e} dln={row["dln"]:.1e} dlnGs={row["dln_Gsym"]:.1e} sign={row["sign_ok"]} Qnp={row["Q_vs_numpy_proj_inverse"]:.1e} kQ={row["kappa_Q"]:.1e}', flush=True)
json.dump({"commit": commit, "seed": 20261022, "c": list(C_BASELINE), "k2_pwave_form": "coded_v2", "rows": rows}, open("/home/user/3body-qc-v4/results/tester_010_compare.json", "w"), indent=1)

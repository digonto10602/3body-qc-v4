"""Tester 012: F (classes) vs oracle vs use_classes=False; per spectator block; ln|det Q| vs oracle. coded_v2 pinned."""
import sys, json, subprocess
import numpy as np
sys.path.insert(0, "/home/user/3body-qc-v4/src/python"); sys.path.insert(0, "/home/user/3body-qc-v4/tests/acceptance/m1")
from qc3 import fast, oracle, _core
from m1_helpers import sector, load_anchors, IRREPS, C_BASELINE
commit = subprocess.check_output(["git", "-C", "/home/user/3body-qc-v4", "rev-parse", "--short", "HEAD"]).decode().strip()
rng = np.random.default_rng(20261022)
def rel(a, b):
    nb = np.linalg.norm(b); return float(np.linalg.norm(a - b) / nb) if nb > 0 else float(np.linalg.norm(a - b))
def blocks(basis):
    out = []; i = 0; n = len(basis)
    while i < n:
        j = i + 1
        while j < n and all(basis[k][j] == basis[k][i] for k in ("flavour", "nx", "ny", "nz")): j += 1
        out.append((i, j)); i = j
    return out
pts = [(p["irrep"], p["Ecm"], "anchor") for p in load_anchors()["points"]]
for irrep in IRREPS:
    for E in rng.uniform(0.26310, 0.36000, 3): pts.append((irrep, float(E), "random"))
for r in json.load(open("/home/user/3body-qc-v4/results/m1_step6_gamma_stress.json")): pts.append((r["irrep"], r["Ecm"], "gamma%.2f" % r["gamma_max"]))
rows = []
for irrep, E, tag in pts:
    s = sector(irrep)
    o = oracle.evaluate_point(s, E, C_BASELINE)
    fc = fast.evaluate_point(s, E, C_BASELINE, want=("F", "Q"))
    fd = fast.evaluate_point(s, E, C_BASELINE, want=("F", "Q"), options={"use_classes": False})
    Fo, Fc, Fd = (np.asarray(x.F) for x in (o, fc, fd))
    bl = blocks(o.basis)
    def maxblock(A, B):
        m = 0.0
        for (i, j) in bl:
            m = max(m, rel(A[i:j, i:j], B[i:j, i:j]))
        return m
    offmask = np.ones(Fo.shape, bool)
    for (i, j) in bl: offmask[i:j, i:j] = False
    row = dict(irrep=irrep, Ecm=E, tag=tag, N=o.N, nblocks=len(bl),
               F_c_vs_o=rel(Fc, Fo), F_d_vs_o=rel(Fd, Fo), F_c_vs_d=rel(Fc, Fd),
               blk_c_vs_o=maxblock(Fc, Fo), blk_c_vs_d=maxblock(Fc, Fd), blk_d_vs_o=maxblock(Fd, Fo),
               offblock_max=float(np.abs(Fc[offmask]).max()) if offmask.any() else 0.0, sym_exact=bool(np.array_equal(Fc, Fc.T)),
               terms_c=fc.metadata["lattice_terms"], terms_d=fd.metadata["lattice_terms"],
               dln_c=abs(fc.lnabsdet - o.lnabsdet), dln_d=abs(fd.lnabsdet - o.lnabsdet), dln_cd=abs(fc.lnabsdet - fd.lnabsdet),
               sign_ok=(fc.sign == o.sign), Nproj=o.Nproj, kappa=float(np.linalg.cond(np.asarray(o.Q))))
    rows.append(row)
    print(f'{irrep} {E:.9f} {tag:10s} N={o.N:3d} Fc/o={row["F_c_vs_o"]:.1e} Fd/o={row["F_d_vs_o"]:.1e} Fc/d={row["F_c_vs_d"]:.1e} blk c/o={row["blk_c_vs_o"]:.1e} c/d={row["blk_c_vs_d"]:.1e} terms {row["terms_c"]}=={row["terms_d"]} dln={row["dln_c"]:.1e}/{row["dln_d"]:.1e} sign={row["sign_ok"]} off={row["offblock_max"]:.0e} k={row["kappa"]:.1e}', flush=True)
json.dump({"commit": commit, "seed": 20261022, "c": list(C_BASELINE), "k2_pwave_form": "coded_v2", "rows": rows}, open("/home/user/3body-qc-v4/results/tester_012_accuracy.json", "w"), indent=1)

"""fast K3 basis and G vs the planner's independent NumPy reference qc3ref (not the oracle)."""
import sys, json, subprocess
import numpy as np
sys.path.insert(0, "/home/user/3body-qc-v4/src/python"); sys.path.insert(0, "/home/user/3body-qc-v4/tests/acceptance/m1")
from qc3 import fast
from m1_helpers import sector, basis_list, E_lab_of, L_BOX, PRODUCTION, C_ZERO, rel_frob
from qc3ref import fgk, k3
MK, MPI = PRODUCTION["atmK"], PRODUCTION["atmpi"]
commit = subprocess.check_output(["git", "-C", "/home/user/3body-qc-v4", "rev-parse", "--short", "HEAD"]).decode().strip()
rows = []
for irrep, E in (("100_A2", 0.2915561), ("000_A1m", 0.26932), ("110_A2", 0.30733), ("111_A2", 0.32036), ("200_A2", 0.28664), ("110_A2", 0.35), ("111_A2", 0.35103), ("000_A1m", 0.34538)):
    r = fast.evaluate_point(sector(irrep), E, C_ZERO, want=("G", "K3"))
    Elab, P = E_lab_of(irrep, E); bl = basis_list(r)
    G_ref = fgk.G_matrix(bl, Elab, P, L_BOX, MK, MPI)
    kb = k3.K3_basis_matrices(bl, Elab, P, L_BOX, MK, MPI)
    row = {"irrep": irrep, "Ecm": E, "N": r.N, "G_vs_qc3ref": rel_frob(r.G, G_ref)}
    for n in ("iso0", "iso1", "B", "E"): row["K3_%s_vs_qc3ref" % n] = rel_frob(r.K3_basis[n], kb[n])
    rows.append(row); print(row, flush=True)
json.dump({"commit": commit, "rows": rows}, open("/home/user/3body-qc-v4/results/tester_010_refcmp.json", "w"), indent=1)

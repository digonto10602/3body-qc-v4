"""Non-production parameters (eps_H != 0, other xi, masses, eta, q0_consistent): fast K3/G/ln|det Q| vs oracle.  Not part of acceptance; generalisation check."""
import sys, json
import numpy as np
sys.path.insert(0, "/home/user/3body-qc-v4/src/python"); sys.path.insert(0, "/home/user/3body-qc-v4/tests/acceptance/m1")
from qc3 import fast, oracle
from qc3.sector import Sector
from m1_helpers import C_BASELINE
rel = lambda a, b: float(np.linalg.norm(np.asarray(a) - np.asarray(b)) / np.linalg.norm(np.asarray(b)))
rows = []
cases = [("100_A2", 0.30, dict(epsilon_H=0.25, xi=3.2, atmK=0.10, atmpi=0.075, eta_2=0.7)),
         ("111_A2", 0.335, dict(epsilon_H=0.5, atmK=0.11, atmpi=0.06906, k2_pwave_form="q0_consistent")),
         ("110_A2", 0.31, dict(epsilon_H=0.0, xi=3.0, alpha=0.4, max_shell_num=12)),
         ("000_A1m", 0.30, dict(epsilon_H=0.1, atmK=0.09698, atmpi=0.12))]
for irrep, E, kw in cases:
    s = Sector.production(irrep, **kw)
    try:
        o = oracle.evaluate_point(s, E, C_BASELINE); f = fast.evaluate_point(s, E, C_BASELINE)
        r = dict(irrep=irrep, Ecm=E, params=kw, N=o.N, Nf=f.N, G=rel(f.G, o.G), K3=rel(f.K3, o.K3), **{"K3_" + n: rel(f.K3_basis[n], o.K3_basis[n]) for n in ("iso0", "iso1", "B", "E")},
                 F=rel(f.F, o.F), dln=abs(f.lnabsdet - o.lnabsdet), sign_ok=(f.sign == o.sign))
    except Exception as e:
        r = dict(irrep=irrep, Ecm=E, params=kw, error=f"{type(e).__name__}: {str(e)[:200]}")
    rows.append(r); print(r, flush=True)
json.dump(rows, open("/home/user/3body-qc-v4/results/tester_010_offprod.json", "w"), indent=1)

import sys, json
import numpy as np
sys.path.insert(0, "/home/user/3body-qc-v4/src/python"); sys.path.insert(0, "/home/user/3body-qc-v4/tests/acceptance/m1")
from qc3 import fast, oracle
from m1_helpers import sector, C_BASELINE
def rel(a,b): return float(np.linalg.norm(a-b)/np.linalg.norm(b))
out=[]
for irrep,E in [("200_A2",0.346802591729),("111_A2",0.337572574251),("111_A2",0.330752399270),("111_A2",0.337780731301)]:
    s=sector(irrep); o=oracle.evaluate_point(s,E,C_BASELINE)
    fc=fast.evaluate_point(s,E,C_BASELINE,want=("F","Q")); fd=fast.evaluate_point(s,E,C_BASELINE,want=("F","Q"),options={"use_classes":False})
    ld=fast.evaluate_point(s,E,C_BASELINE,want=("F","Q"),options={"trunc_eps":1e-40})
    row=dict(irrep=irrep,E=E,Fc_o=rel(fc.F,o.F),Fd_o=rel(fd.F,o.F),Fc_d=rel(fc.F,fd.F),Fd40_o=rel(ld.F,o.F),dln_c=abs(fc.lnabsdet-o.lnabsdet),dln_d=abs(fd.lnabsdet-o.lnabsdet),kappa=float(np.linalg.cond(np.asarray(o.Q))), maxF=float(np.abs(o.F).max()), medF=float(np.median(np.abs(o.F[o.F!=0]))))
    print(row,flush=True); out.append(row)
json.dump(out,open("/home/user/3body-qc-v4/results/tester_012_outliers.json","w"),indent=1)

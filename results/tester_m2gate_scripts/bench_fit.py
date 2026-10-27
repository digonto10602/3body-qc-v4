import time, json, subprocess, numpy as np
from qc3 import fit, lattice, constants, _core
t0=time.perf_counter()
fits=fit.cutoff_scan(lattice.default_data_dir(),20,[0.335],constants.C_BASELINE,'q0_consistent',verify=True,multistart=0)
wall=time.perf_counter()-t0
r=fits[0]
ver=r.verification['seconds'] if r.verification else None
d=dict(commit=subprocess.run(['git','rev-parse','--short','HEAD'],capture_output=True,text=True).stdout.strip(),core_commit=_core.git_commit(),cutoff=0.335,table='A (4 free parameters)',form='q0_consistent',c_start=list(constants.C_BASELINE),
  wall_total_s=wall,seconds_scan=r.seconds_scan,seconds_fit_cutoff=r.seconds,seconds_lm=r.seconds_lm,seconds_migrad=r.seconds_migrad,seconds_verify=ver,
  seconds_other_in_fit=r.seconds-r.seconds_lm-r.seconds_migrad, n_chi2_calls=r.n_chi2_calls,n_evaluations=r.n_evaluations,nfcn=r.nfcn,lm_nfev=r.lm_nfev,lm_njev=r.lm_njev,
  chi2=r.chi2,ndof=r.ndof,edm=r.edm,valid=r.valid,accurate=r.accurate,c=list(r.c),n_local_rescans=r.n_local_rescans,n_full_rescans=r.n_full_rescans,verification_ok=r.verification['ok'] if r.verification else None,seed=None)
json.dump(d,open('results/bench_m2_gate_fit0335.json','w'),indent=1)
for k,v in d.items(): print(k,v)
import resource; print("peak_rss_MB",resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024)

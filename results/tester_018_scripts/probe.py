import json
from multiprocessing import Pool
from indep import *
from run_min import job
HOME="/home/user/3body-qc-v4/results/"
d=json.load(open(HOME+"fit_L20_q0_consistent_KBKE0_m2.json")); f=d['fits'][2]
labels=[tuple(l) for l in f['labels']]; form=f['k2_pwave_form']
c0=np.array(f['c']); err=np.array(f['errors'])
Emax={ir:max(f['E_model'][i] for i,l in enumerate(labels) if l[0]==ir)+0.004 for ir in set(l[0] for l in labels)}
ss=[-0.3,-0.6,-0.8,-0.9,-1.0,-1.1,-1.3,-1.6,-2.0,-3.0]
reqs=[]
for s in ss:
    v=c0.copy(); v[1]+=s*err[1]; reqs.append(v)
if __name__=="__main__":
    with Pool(4) as p: vals=p.map(job,[(v,labels,form,Emax) for v in reqs])
    for s,v in zip(ss,vals): print("dKiso1 = %.1f sigma_H: chi2 - chi2_min ="%s, (v-f['chi2']) if not isinstance(v,tuple) else v)

import json
from multiprocessing import Pool
from indep import *
from run_min import job
HOME="/home/user/3body-qc-v4/results/"
d=json.load(open(HOME+"fit_L20_q0_consistent_m2.json")); f=d['fits'][2]
labels=[tuple(l) for l in f['labels']]; form=f['k2_pwave_form']
c0=np.array(f['c']); err=np.array(f['errors'])
Emax={ir:max(f['E_model'][i] for i,l in enumerate(labels) if l[0]==ir)+0.004 for ir in set(l[0] for l in labels)}
reqs={}
for i in range(4):
    for sc in (0.003,0.01,0.03):
        for s in (1,-1):
            v=c0.copy(); v[i]+=s*sc*err[i]; reqs[(i,sc,s)]=v
keys=list(reqs)
if __name__=="__main__":
    with Pool(2) as p: vals=p.map(job,[(reqs[k],labels,form,Emax) for k in keys])
    F=dict(zip(keys,vals))
    f0=f['chi2']
    for i in range(4):
        for sc in (0.003,0.01,0.03):
            fp,fm=F[(i,sc,1)],F[(i,sc,-1)]
            print(i,sc,"fp-f0 %.6e fm-f0 %.6e  g*sigma %.4e  curv H*sig^2 %.4f"%(fp-f0,fm-f0,(fp-fm)/(2*sc),(fp+fm-2*f0)/sc**2))

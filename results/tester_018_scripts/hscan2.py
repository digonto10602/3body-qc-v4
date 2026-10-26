import json,sys
from multiprocessing import Pool
from indep import *
from run_min import job
HOME="/home/user/3body-qc-v4/results/"
def go(fn,idx,tag):
    d=json.load(open(HOME+fn)); f=d['fits'][idx]
    labels=[tuple(l) for l in f['labels']]; form=f['k2_pwave_form']
    free=[i for i in range(4) if f['errors'][i]>0]
    c0=np.array(f['c']); err=np.array(f['errors_fd_hessian'])
    Emax={ir:max(f['E_model'][i] for i,l in enumerate(labels) if l[0]==ir)+0.004 for ir in set(l[0] for l in labels)}
    reqs={}
    for i in free:
        for sc in (0.003,0.01,0.03):
            for s in (1,-1):
                v=c0.copy(); v[i]+=s*sc*err[i]; reqs[(i,sc,s)]=v
    keys=list(reqs)
    with Pool(4) as p: vals=p.map(job,[(reqs[k],labels,form,Emax) for k in keys])
    F=dict(zip(keys,vals)); f0=f['chi2']
    print(tag,"(steps in units of the FD-Hessian sigma)")
    for i in free:
        for sc in (0.003,0.01,0.03):
            fp,fm=F[(i,sc,1)],F[(i,sc,-1)]
            if isinstance(fp,tuple) or isinstance(fm,tuple): print(i,sc,"ERR",fp if isinstance(fp,tuple) else "", fm if isinstance(fm,tuple) else ""); continue
            print(i,sc,"g*sigma %.4e  Hsig2/2 %.4f  fp-f0 %.3e fm-f0 %.3e"%((fp-fm)/(2*sc),(fp+fm-2*f0)/sc**2/2,fp-f0,fm-f0),flush=True)
if __name__=="__main__":
    go("fit_L20_q0_consistent_KBKE0_m2.json",2,"B0335")
    go("fit_L20_q0_consistent_KBKE0_m2.json",3,"B0345")
    go("fit_L20_q0_consistent_m2.json",3,"A0345")

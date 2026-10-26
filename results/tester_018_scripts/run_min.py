import json,sys,time,itertools
from multiprocessing import Pool
from indep import *
HOME="/home/user/3body-qc-v4/results/"
def job(args):
    c,labels,form,Emax=args
    mean,cov,Cinv=make_chi(labels)
    try:
        return chi2_at(c,labels,mean,Cinv,form,Emax)[0]
    except Exception as e:
        return ("ERR",type(e).__name__+": "+str(e)[:150])
def run(tag,fname,idx):
    d=json.load(open(HOME+fname)); f=d['fits'][idx]
    labels=[tuple(l) for l in f['labels']]; form=f['k2_pwave_form']
    free=[i for i in range(4) if i not in [ ["Kiso0","Kiso1","KB","KE"].index(n) for n in f['fixed']]]
    c0=np.array(f['c']); err=np.array(f['errors']); cov=np.array(f['covariance'])
    Emax={ir:max(f['E_model'][i] for i,l in enumerate(labels) if l[0]==ir)+0.004 for ir in set(l[0] for l in labels)}
    h={i:0.03*err[i] for i in free}
    pts={}
    def key(v): return tuple(np.round(v,6))
    reqs={"0":c0.copy()}
    for i in free:
        for s in (+1,-1):
            v=c0.copy(); v[i]+=s*h[i]; reqs[f"{s}{i}"]=v
    for i,j in itertools.combinations(free,2):
        for s in (+1,-1):
            v=c0.copy(); v[i]+=s*h[i]; v[j]+=s*h[j]; reqs[f"{s}{i}{j}"]=v
    # conditional 1-sigma steps using reported Hesse covariance
    for i in free:
        sh=cov[:,i]/err[i]
        for sc in (0.1,0.3,1.0):
            for s in (+1,-1):
                reqs[f"cond{s*sc}_{i}"]=c0+s*sc*sh
        for s in (+1,-1):
            v=c0.copy(); v[i]+=s*err[i]; reqs[f"ax{s}_{i}"]=v
    keys=list(reqs)
    t=time.time()
    with Pool(4) as p:
        vals=p.map(job,[(reqs[k],labels,form,Emax) for k in keys])
    F=dict(zip(keys,vals))
    print('ERRS',{k:v for k,v in F.items() if isinstance(v,tuple)},flush=True)
    json.dump({k:v for k,v in F.items()},open(f'raw_{tag}.json','w'))
    errs={k:v for k,v in F.items() if isinstance(v,(tuple,list))}
    F={k:(float('nan') if isinstance(v,(tuple,list)) else v) for k,v in F.items()}
    f0=F["0"]
    out={"tag":tag,"chi2_reported":f['chi2'],"chi2_mine":f0,"rel_diff":abs(f0/f['chi2']-1),"free":free,"c":c0.tolist(),"hesse_err":err.tolist(),"fd_analytic_err":f['errors_fd_hessian']}
    fr=free; n=len(fr)
    H=np.zeros((n,n)); g=np.zeros(n)
    for a,i in enumerate(fr):
        g[a]=(F[f"1{i}"]-F[f"-1{i}"])/(2*h[i])
        H[a,a]=(F[f"1{i}"]+F[f"-1{i}"]-2*f0)/h[i]**2
    for a,i in enumerate(fr):
        for b,j in enumerate(fr):
            if a<b:
                k1=f"1{i}{j}";k2=f"-1{i}{j}"
                H[a,b]=H[b,a]=(F[k1]+F[k2]-F[f"1{i}"]-F[f"1{j}"]+2*f0-F[f"-1{i}"]-F[f"-1{j}"])/(2*h[i]*h[j])
    covm=2*np.linalg.inv(H)
    errm=np.sqrt(np.diag(covm))
    out["grad_times_sigma"]=[g[a]*err[i] for a,i in enumerate(fr)]
    out["grad_times_mysigma"]=[g[a]*errm[a] for a in range(n)]
    out["H_eig"]=np.linalg.eigvalsh(H).tolist()
    out["my_fd_hessian_errors"]=errm.tolist()
    def dv(k):
        v=F[k]; return None if v!=v else v-f0
    out["cond_dchi2"]={str(i):{sc:[dv(f"cond{sg*sc}_{i}") for sg in (1,-1)] for sc in (0.1,0.3,1.0)} for i in fr}
    out["axis_dchi2"]={str(i):[dv(f"ax1_{i}"),dv(f"ax-1_{i}")] for i in fr}
    # conditional steps using MY covariance
    out["axis_expected_1_over_1mR2"]={str(i):float(1/(1-(covm[a,[b for b in range(n) if b!=a]]**2).sum()/1)) for a,i in enumerate(fr)} if False else None
    # reported Hesse covariance implied axis dchi2 = 0.5 H_ii sigma^2 (my H)
    out["axis_dchi2_from_myH"]={str(i):0.5*H[a,a]*err[i]**2 for a,i in enumerate(fr)}
    # predicted cond step Delta chi2 using my H and the Hesse cov direction
    out["pred_cond_dchi2_myH"]={str(i):0.5*(cov[:,i]/err[i])[fr]@H@(cov[:,i]/err[i])[fr] for i in fr}
    out["seconds"]=time.time()-t
    out["F"]=F; out["errs"]=errs
    json.dump(out,open(f"min_{tag}.json","w"),indent=1)
    print(tag,{k:v for k,v in out.items() if k!="F"},flush=True)
if __name__=="__main__":
    run("A0345","fit_L20_q0_consistent_m2.json",3)
    run("B0335","fit_L20_q0_consistent_KBKE0_m2.json",2)
    run("B0345","fit_L20_q0_consistent_KBKE0_m2.json",3)

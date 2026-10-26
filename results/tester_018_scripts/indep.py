import sys, json, glob, os, re, time
import numpy as np
sys.path.insert(0, "/home/user/3body-qc-v4/src/python")
os.environ.setdefault("OMP_NUM_THREADS","1"); os.environ.setdefault("OPENBLAS_NUM_THREADS","1")
from qc3.sector import Sector
from qc3 import rootscan, fast
DATA="/home/user/3body-qc-v4/data/lattice/L20"
XI=3.444; L=20
DDD={"000_A1m":("000",0.0),"100_A2":("001",1.0),"200_A2":("002",2.0),"110_A2":("011",np.sqrt(2)),"111_A2":("111",np.sqrt(3))}

def D1_levels(labels):
    """labels list of (irrep,n). Own numpy implementation of the D1 recipe."""
    cols=[];means=[]
    for irrep,n in labels:
        ddd,dn=DDD[irrep]
        fn=f"{DATA}/20_{ddd}_{'A1m' if irrep.endswith('A1m') else 'A2'}_n{n}.jack"
        lines=open(fn).read().split("\n")
        N=int(lines[0].split()[0])
        x=np.array([float(t.split()[1]) for t in lines[1:] if t.strip()])
        assert len(x)==N
        y=(x.sum()-x)/(N-1)
        P=2*np.pi*dn/(XI*L)
        e=np.sqrt(y**2-P**2)
        cols.append(e)
    Y=np.column_stack(cols); N=Y.shape[0]
    mean=Y.mean(0); d=Y-mean
    cov=(N-1)/N*d.T@d
    return mean,0.5*(cov+cov.T),N

def make_chi(labels):
    mean,cov,N=D1_levels(labels)
    w,U=np.linalg.eigh(cov)
    keep=w>1e-12*w.max()
    assert keep.all(), w
    Cinv=(U/w)@U.T
    return mean,cov,Cinv

def zeros_fresh(irrep,c,Emax,form):
    sec=Sector.production(irrep,20,k2_pwave_form=form)
    sr=rootscan.scan_sector(fast.evaluate_point,sec,tuple(c),0.2631,Emax)
    bad=[x for x in sr.candidates if x.kind in("ambiguous","even_zero") and x.E<=Emax]
    assert not bad, (irrep,[(b.kind,b.E) for b in bad])
    return sr.zeros(), [x for x in sr.of_kind("zero") if x.slope is not None and x.slope<0]

def chi2_at(c,labels,mean,Cinv,form,Emax_by_irrep):
    """D3: lowest m QC zeros (kind zero only) matched to sorted lattice levels per irrep."""
    Em=np.empty(len(labels))
    for irrep in sorted(set(l[0] for l in labels)):
        idx=[i for i,l in enumerate(labels) if l[0]==irrep]
        z,_=zeros_fresh(irrep,c,Emax_by_irrep[irrep],form)
        order=sorted(idx,key=lambda i:mean[i])
        assert len(z)>=len(idx),(irrep,len(z),len(idx))
        for k,i in enumerate(order): Em[i]=z[k]
    r=mean-Em
    return float(r@Cinv@r),Em

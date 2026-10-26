import json, math, numpy as np
HOME="/home/user/3body-qc-v4/"
xi,L=3.444,20.0
DDD={"000_A1m":("000","A1m",0),"100_A2":("001","A2",1),"200_A2":("002","A2",4),"110_A2":("011","A2",2),"111_A2":("111","A2",3)}
def chi(r,C):
    w,v=np.linalg.eigh(C); keep=np.abs(w)>1e-12*np.abs(w).max()
    Ci=(v[:,keep]/w[keep])@v[:,keep].T
    return float(r@Ci@r), int(keep.sum())
for fname in ["fit_L20_q0_consistent_m2","fit_L20_q0_consistent_KBKE0_m2"]:
  d=json.load(open(HOME+"results/"+fname+".json"))
  for f in d['fits']:
    S1=[];S2=[]
    for irrep,n in f['labels']:
        ddd,nm,n2=DDD[irrep]
        x=np.loadtxt(f"{HOME}data/lattice/L20/20_{ddd}_{nm}_n{n}.jack",skiprows=1)[:,1]
        P=2*math.pi*math.sqrt(n2)/(xi*L)
        S1.append(np.sqrt(x**2-P**2))                          # tool: convert each configuration
        y=(x.sum()-x)/(len(x)-1); S2.append(np.sqrt(y**2-P**2))  # D1: convert each jackknife resample
    S1=np.array(S1);S2=np.array(S2);N=S1.shape[1]
    Em=np.array(f['E_model'])
    out=[]
    # tool recipe
    dev=S1-S1.mean(1,keepdims=True); Ct=(N-1)/N*dev@dev.T/(N-1)**2; rt=Em-S1.mean(1)
    c_tool,_=chi(rt,Ct)
    # D1 pipeline: cov (N-1)/N sum dev^2 of resamples, mean of resamples
    dev2=S2-S2.mean(1,keepdims=True); C2=(N-1)/N*dev2@dev2.T; r2=Em-S2.mean(1)
    c_d1,_=chi(r2,C2)
    # hybrids
    c_h1,_=chi(Em-S2.mean(1),Ct)   # tool cov with D1 mean
    c_h2,_=chi(rt,C2)              # D1 cov with tool mean
    print(fname[-12:],f['cutoff'],"reported %.8f  D1(recomputed) %.8f  tool %.6f  diff(tool-ours) %.4f | tool-cov+D1-mean %.6f  D1-cov+tool-mean %.6f"%(f['chi2'],c_d1,c_tool,c_tool-f['chi2'],c_h1,c_h2))
    print("     cov rel diff Ct vs C2: %.3e; mean shift max %.2e (in sigma %.3f)"%(np.abs(Ct-C2).max()/np.abs(C2).max(),np.abs(S1.mean(1)-S2.mean(1)).max(),np.max(np.abs(S1.mean(1)-S2.mean(1))/np.sqrt(np.diag(C2)))))

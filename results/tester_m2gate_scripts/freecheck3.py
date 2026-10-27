import json, numpy as np
from qc3 import fast, constants, free_levels
from qc3.sector import Sector
res = json.load(open('results/fit_L20_q0_consistent_m2.json'))
cA = np.array([f for f in res['fits'] if f['cutoff']==0.335][0]['c']); cb=np.array(constants.C_BASELINE)
rng = np.random.default_rng(20261027)
cs = {'baseline': cb, 'A0.335_min': cA}
for i in (1,2): cs[f'rand{i}_pm10%'] = cb*(1+rng.uniform(-0.1,0.1,4))
cs['c=0']=np.zeros(4)
ds=np.geomspace(2e-8,3e-5,25)
out={}
for form in ('coded_v2','q0_consistent'):
    sec = Sector.production('110_A2', Lbyas=20, k2_pwave_form=form)
    Ef = float(free_levels.free_levels(sec, 0.3190, 0.3200)[0])
    for name,c in cs.items():
        c=tuple(float(x) for x in c); r={}
        for s,lab in ((1,'+'),(-1,'-')):
            y=np.array([np.exp(fast.evaluate_point(sec,Ef+s*d,c).lnabsdet/4) for d in ds])
            # y = |p1 (E-E0)| (1 + b (E-E0)): fit y = p0 + p1*x + p2*x^2 with x=s*d (signed offset from E_f); E0 = root of p0+p1 x
            X=s*ds; P=np.polyfit(X,y,2)
            roots=np.roots(P); r[lab]=P[2]/(-P[1]) if False else float(min(roots,key=abs).real)
        out[f'{form}|{name}']=r
        print(f'{form:14s} {name:12s} E0-E_f from +side fit: {r["+"]:+.2e}  from -side fit: {r["-"]:+.2e}',flush=True)
json.dump(out,open('results/tester_m2gate_freelevel_quartic_fit.json','w'),indent=1)

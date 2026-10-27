import json, numpy as np
from qc3 import fast, rootscan, constants, free_levels
from qc3.sector import Sector
out=json.load(open('results/tester_m2gate_freelevel_c_independence.json'))
res = json.load(open('results/fit_L20_q0_consistent_m2.json'))
cA = np.array([f for f in res['fits'] if f['cutoff']==0.335][0]['c']); cb=np.array(constants.C_BASELINE)
rng = np.random.default_rng(20261027)
cs = {'baseline': cb, 'A0.335_min': cA}
for i in (1,2): cs[f'rand{i}_pm10%_of_baseline'] = cb*(1+rng.uniform(-0.1,0.1,4))
cs['c=0']=np.zeros(4)
fine={}
for form in ('coded_v2','q0_consistent'):
    sec = Sector.production('110_A2', Lbyas=20, k2_pwave_form=form)
    Ef = float(free_levels.free_levels(sec, 0.3190, 0.3200)[0])
    for name,c in cs.items():
        c=tuple(float(x) for x in c)
        rr=out[f'{form}|{name}']['relabelled'][0]
        # fine sign pattern in +-2e-9
        d=np.linspace(-2e-9,2e-9,2001)
        sg=[];ln=[]
        for x in d:
            r=fast.evaluate_point(sec,Ef+x,c); sg.append(r.sign); ln.append(r.lnabsdet)
        sg=np.array(sg); ch=np.where(sg[1:]!=sg[:-1])[0]
        # exponent: local slope of ln|det| vs ln|d| at 1e-8..1e-6 (both sides)
        def l(x): return fast.evaluate_point(sec,Ef+x,c).lnabsdet
        ex=[(l(s*1e-6)-l(s*1e-7))/np.log(10) for s in (+1,-1)]
        ex2=[(l(s*1e-7)-l(s*1e-8))/np.log(10) for s in (+1,-1)]
        fine[f'{form}|{name}']=dict(E_orig_minus_Ef=rr['E_orig']-Ef, order=rr['order'], slope=rr['slope'], n_sign_flips_pm2e9=len(ch),
             flip_positions_first=[float(d[i]) for i in ch[:6]], exp_1e7_1e6=ex, exp_1e8_1e7=ex2)
        print(form,name,'E_orig-Ef=%.3e'%(rr['E_orig']-Ef),'order=%.5f'%rr['order'],'nflips(+-2e-9,2001pts)=',len(ch),'exponent(1e-6/1e-7; +,-)=',np.round(ex,3),' (1e-7/1e-8)',np.round(ex2,3),flush=True)
json.dump(fine,open('results/tester_m2gate_freelevel_fine.json','w'),indent=1,default=float)

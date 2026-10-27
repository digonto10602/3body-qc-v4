import json, sys, numpy as np, time
from qc3 import fast, rootscan, constants, free_levels
from qc3.sector import Sector
from scipy.optimize import brentq
res = json.load(open('results/fit_L20_q0_consistent_m2.json'))
cA = [f for f in res['fits'] if f['cutoff']==0.335][0]['c']
cb = np.array(constants.C_BASELINE)
rng = np.random.default_rng(20261027)
cs = {'baseline': cb, 'A0.335_min': np.array(cA)}
for i in (1,2):
    cs[f'rand{i}_pm10%_of_baseline'] = cb*(1+rng.uniform(-0.1,0.1,4))
cs['c=0'] = np.zeros(4)
out = {}
for form in ('coded_v2','q0_consistent'):
    sec = Sector.production('110_A2', Lbyas=20, k2_pwave_form=form)
    Ef = float(free_levels.free_levels(sec, 0.3190, 0.3200)[0])
    print(form, 'E_f=%.14f'%Ef, flush=True)
    for name,c in cs.items():
        c=tuple(float(x) for x in c)
        t=time.time()
        sr = rootscan.scan_sector(fast.evaluate_point, sec, c, 0.3170, 0.3220)
        rec = {'c':list(c)}
        near=[x for x in sr.candidates if abs(x.E-Ef)<2e-3]
        rec['candidates_near']=[dict(E=x.E,dE=x.E-Ef,kind=x.kind,order=x.order,slope=x.slope,sign_change=x.sign_change,free_level=x.free_level) for x in near]
        rec['relabelled']=sr.relabelled_free_levels
        rec['skipped']=sr.skipped_free_levels
        # direct probe of det Q around E_f
        def f(E):
            r=fast.evaluate_point(sec,E,c); return r.sign, r.lnabsdet, r.N
        probe=[]
        for d in (1e-3,1e-4,1e-6,1e-8,1e-9,1e-10,1e-11,1e-12,1e-13):
            sm=f(Ef-d); sp=f(Ef+d)
            probe.append(dict(d=d,sign_m=sm[0],ln_m=sm[1],sign_p=sp[0],ln_p=sp[1],N=sm[2]))
        rec['probe']=probe
        # sign-change bisection of nearest sign-change within 1e-8 of E_f
        sg=lambda E: f(E)[0]
        a,b=Ef-1e-8,Ef+1e-8
        rec['sign_a_b']=[sg(a),sg(b)]
        if sg(a)!=sg(b):
            lo,hi=a,b
            for _ in range(80):
                m=0.5*(lo+hi)
                if sg(m)==sg(lo): lo=m
                else: hi=m
            rec['bisect_zero_minus_Ef']=0.5*(lo+hi)-Ef
        # nearest ordinary zeros / levels in window
        rec['zeros_in_window']=[float(z) for z in sr.zeros()]
        rec['seconds']=time.time()-t
        out[f'{form}|{name}']=rec
        print(form,name,'c=',[f'{x:.4g}' for x in c], flush=True)
        for x in rec['candidates_near']: print('   cand',x)
        print('   relabelled',[ (r['E']-Ef, r['kind_orig'], r['order'], r['slope']) for r in rec['relabelled']],' signs a,b',rec['sign_a_b'],' bisect',rec.get('bisect_zero_minus_Ef'), flush=True)
        print('   zeros in window',[round(z,8) for z in rec['zeros_in_window']], flush=True)
        print('   probe (d: sign-,ln-, sign+,ln+):',[(p['d'],p['sign_m'],round(p['ln_m'],3),p['sign_p'],round(p['ln_p'],3)) for p in probe],flush=True)
json.dump(out,open('results/tester_m2gate_freelevel_c_independence.json','w'),indent=1,default=float)

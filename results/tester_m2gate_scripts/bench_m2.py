import sys, time, json, platform, subprocess, resource, numpy as np
from qc3 import fast, rootscan, constants, fit, lattice, _core
from qc3.sector import Sector
IR=list(constants.IRREPS); cb=np.array(constants.C_BASELINE)
out={'commit':subprocess.run(['git','rev-parse','--short','HEAD'],capture_output=True,text=True).stdout.strip(),'core_commit':_core.git_commit(),
     'form':'q0_consistent','seed':20261027,'c_baseline':list(cb),'machine':platform.platform(),'pinned':'taskset core 0, OMP/OPENBLAS=1'}
# (b) scan_sector per sector
b={}
for Emax in (0.34,0.36):
    for ir in IR:
        sec=Sector.production(ir,Lbyas=20,k2_pwave_form='q0_consistent')
        ts=[];ne=None
        for rep in range(2):
            t=time.perf_counter(); r=rootscan.scan_sector(fast.evaluate_point,sec,tuple(cb),constants.GRID_E0,Emax); ts.append(time.perf_counter()-t); ne=r.n_evaluations
        kinds={}
        for c in r.candidates: kinds[c.kind]=kinds.get(c.kind,0)+1
        b[f'{ir}|{Emax}']=dict(evals=ne,wall_s=ts,kinds=kinds,n_grid=r.n_grid,internal_seconds=r.seconds)
        print(f'(b) {ir} E_max={Emax}: evals={ne} wall min {min(ts):.2f} s  (reps {ts[0]:.2f},{ts[1]:.2f}) kinds {kinds}',flush=True)
out['b_scan_sector']=b
# (c) chi2 call cost
lv=lattice.load_levels(lattice.default_data_dir(),20,IR,0.335)
c={}
for rel in (1e-3,2e-2):
    trk=fit.make_trackers(IR,tuple(cb),Lbyas=20,k2_pwave_form='q0_consistent',E_min=constants.GRID_E0,E_max=0.335+0.005)
    f=fit.QCChi2(lv,trk); f(cb)
    rng=np.random.default_rng(20261027); ts=[];ne=[];chi=[]
    for k in range(60):
        c_new=cb*(1+rel*rng.choice([-1,1],4)*rng.uniform(0.5,1,4))
        e0=f.n_evaluations; t=time.perf_counter(); v=f(c_new); ts.append(time.perf_counter()-t); ne.append(f.n_evaluations-e0); chi.append(v)
    ts=np.array(ts)*1e3
    c[str(rel)]=dict(ms=ts.tolist(),min=float(ts.min()),median=float(np.median(ts)),max=float(ts.max()),p90=float(np.percentile(ts,90)),evals_median=float(np.median(ne)),evals_min=int(min(ne)),evals_max=int(max(ne)),
        local=sum(t.n_local_rescans for t in trk.values()),full=sum(t.n_full_rescans for t in trk.values()),pen=f.n_penalised_calls,chi2_median=float(np.median(chi)))
    print(f'(c) rel={rel}: n=60 ms min/med/p90/max = {ts.min():.1f}/{np.median(ts):.1f}/{np.percentile(ts,90):.1f}/{ts.max():.1f}; evals/call med {np.median(ne)} [{min(ne)},{max(ne)}]; local {c[str(rel)]["local"]} full {c[str(rel)]["full"]} pen {f.n_penalised_calls}',flush=True)
out['c_chi2_cost']=c
json.dump(out,open('results/bench_m2_gate_scan_chi2.json','w'),indent=1)

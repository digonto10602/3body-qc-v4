import sys, json, numpy as np
from qc3 import fit, lattice
from qc3.constants import C_BASELINE
tab, cut = sys.argv[1], float(sys.argv[2]); out=sys.argv[3]
fixed = (2,3) if tab=='B' else ()
cs0 = np.array(C_BASELINE, float)
if fixed: cs0[list(fixed)] = 0.0
res=[]
for eps, seed in [(0.0,0),(1e-10,1),(1e-10,2),(1e-8,1),(1e-8,2),(1e-6,1)]:
    rng=np.random.default_rng(seed)
    cs = cs0*(1+eps*rng.standard_normal(4)) if eps>0 else cs0.copy()
    if fixed: cs[list(fixed)]=0.0
    lv = lattice.load_levels(lattice.default_data_dir(), 20, list(fit.IRREPS), cut)
    used=[ir for ir in fit.IRREPS if any(l[0]==ir for l in lv.labels)]
    tr = fit.make_trackers(used, tuple(cs), Lbyas=20, k2_pwave_form='q0_consistent', E_min=fit.GRID_E0, E_max=min(0.36,cut+fit.TRACK_MARGIN))
    r = fit.fit_cutoff(lv, tr, tuple(cs), fixed=tuple(fit.PARAM_NAMES[i] for i in fixed), verify=False)
    d=dict(eps=eps,seed=seed,chi2=r.chi2,edm=r.edm,c=list(map(float,r.c)),hesse=list(map(float,r.errors)),fd=list(map(float,r.errors_fd_hessian)),
           gap=[float(x) for x in r.error_gap],accurate=bool(r.accurate),reasons=list(r.unconstrained_reasons),disagree=list(r.errors_disagree_params),nfcn=r.nfcn)
    res.append(d); print(json.dumps(d),flush=True)
    json.dump(res,open(out,'w'))

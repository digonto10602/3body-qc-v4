import json
from multiprocessing import Pool
from indep import *
from run_min import job
d=json.load(open("/home/user/3body-qc-v4/results/fit_L20_coded_v2_code_equivalence_m2.json")); ce=d['code_equivalence']; f=d['fits'][0]
labels=[tuple(l) for l in ce['labels']]
Emax={ir:0.345 for ir in set(l[0] for l in labels)}
cs=[ce['c_v32w']['c'],ce['c_legacy_0335']['c'],f['c']]
if __name__=="__main__":
    with Pool(3) as p: vals=p.map(job,[(c,labels,"coded_v2",Emax) for c in cs])
    print("mine",vals); print("reported",ce['c_v32w']['chi2_corrected_D1'],ce['c_legacy_0335']['chi2_corrected_D1'],f['chi2'])

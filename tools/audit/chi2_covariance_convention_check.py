"""Recompute the legacy L20 pass-only chi^2 under two readings of the .jack files.

asis      : file values are jackknife samples (what qc_fitter_norm_refine_v2_multiL.cpp assumes)
corrected : file values are per-configuration values (what lattice_data_covariance_cpp.hpp assumes)
The two covariances differ by exactly (N-1)^2. See docs/OPEN_ISSUES.md, issue 1.
Usage (from repo root): python3 tools/audit/chi2_covariance_convention_check.py [repo_root]
"""
import sys, csv, os, math, numpy as np
root = sys.argv[1] if len(sys.argv) > 1 else "."
jackdir = os.path.join(root, "data/lattice/L20")
xi, L = 3.444, 20.0
for cut in ["0p315","0p325","0p335","0p345"]:
    f = os.path.join(root, f"results/legacy_L20_passonly_hotwindow/Ecm_cut_{cut}/lattice_vs_model_levels.csv")
    rows = list(csv.DictReader(open(f)))
    samples, model, lat = [], [], []
    for r in rows:
        fn = os.path.basename(r["source_lattice_file"])
        d = np.loadtxt(os.path.join(jackdir, fn), skiprows=1)[:,1]
        n = fn.split("_")[1]; n2 = sum(int(c)**2 for c in n)
        P = 2*math.pi*math.sqrt(n2)/(xi*L)
        samples.append(np.sqrt(d**2 - P**2)); model.append(float(r["model_Ecm"])); lat.append(float(r["lattice_Ecm"]))
    S = np.array(samples); N = S.shape[1]
    dev = S - S.mean(axis=1, keepdims=True)
    cov_asis = (N-1)/N * dev @ dev.T            # fitter convention (files treated as JK samples)
    cov_raw = cov_asis / (N-1)**2               # files treated as rescaled ensemble -> JK cov
    r = np.array(model) - S.mean(axis=1)
    def chi2(C):
        w, v = np.linalg.eigh(C); keep = np.abs(w) > 1e-12*np.abs(w).max()
        Ci = (v[:,keep]/w[keep]) @ v[:,keep].T
        return float(r @ Ci @ r)
    npar = 4; dof = len(r)-npar
    print(cut, "nlev=%d" % len(r), "chi2_asis=%.4e" % chi2(cov_asis), "chi2_corrected=%.3f" % chi2(cov_raw),
          "dof=%d" % dof, "chi2/dof_corrected=%.3f" % (chi2(cov_raw)/dof),
          "cond(corr)=%.2e" % np.linalg.cond(cov_raw/np.sqrt(np.outer(np.diag(cov_raw),np.diag(cov_raw)))),
          "max|r|/sigma=%.2f" % np.max(np.abs(r)/np.sqrt(np.diag(cov_raw))))

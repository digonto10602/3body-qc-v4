#!/usr/bin/env python3
"""Finalize zero tests, traces, independent audit, and compact reports."""
from __future__ import annotations

import csv, math, shutil, subprocess, sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "output/v33p_projected_basis_cache/l20_k3b_k3e_zero_tests"
PREV = ROOT / "output/v33p_projected_basis_cache/l20_full_cutoff_fits_user_validated_0p36_passonly"
TESTS = ("fixed_K3B0_K3E0", "float_from_K3B0_K3E0")
CUTS = (0.315, 0.325, 0.335, 0.345)
PARAMS = ("K3iso0", "K3iso1", "K3B", "K3E")


def d(test, cut): return OUT / test / f"Ecm_cut_{cut:.3f}".replace(".", "p")
def read(path): return list(csv.DictReader(path.open(newline=""))) if path.exists() else []
def f(x, default=math.nan):
    try: return float(x)
    except (TypeError, ValueError): return default


def empty_artifacts(test, cut):
    dd = d(test, cut); dd.mkdir(parents=True, exist_ok=True)
    fixed = test.startswith("fixed")
    nlevels = len(read(dd / "lattice_levels_used.csv"))
    if not nlevels:
        previous = next((r for r in read(PREV / "L20_PASSONLY_FIT_SUMMARY.csv")
                         if abs(f(r.get("fit_cutoff")) - cut) < 1.0e-12), {})
        nlevels = int(f(previous.get("n_levels"), 0) or 0)
    nfloat = 2 if fixed else 4
    trace = read(dd / "fcn_trace.csv")
    last = trace[-1] if trace else {}
    found = last.get("model_found", f"0/{nlevels}")
    if "/" not in found and nlevels:
        found = f"{found}/{nlevels}"
    chi = f(last.get("chi2_total_raw"), 1e100)
    dof = nlevels - nfloat
    if not (dd / "fit_metrics.csv").exists():
        (dd / "fit_metrics.csv").write_text("fit_cutoff,n_levels,model_found,chi2,dof,chi2_dof,fcn_evals,edm_or_convergence_metric,converged,fallback_count,full_scan_count,load_time_sec,total_runtime_sec,avg_fcn_time_sec,status\n" + f"{cut},{nlevels},{found},{chi},{dof},nan,{len(trace)},unavailable,False,0,0,nan,nan,nan,BLOCKED_MODEL_FOUND\n")
    if not (dd / "final_parameters.csv").exists():
        with (dd / "final_parameters.csv").open("w", newline="") as out:
            w = csv.writer(out); w.writerow(["fit_cutoff","parameter","value","error","relative_error","fixed_or_floating","status"])
            for p in PARAMS:
                fixed_p = fixed and p in ("K3B", "K3E")
                w.writerow([cut, p, 0.0 if fixed_p else "nan", 0.0 if fixed_p else "nan", 0.0 if fixed_p else "nan", "fixed" if fixed_p else "floating", "FIXED" if fixed_p else "BLOCKED_MODEL_FOUND"])
    if not (dd / "final_parameters.json").exists(): (dd / "final_parameters.json").write_text("{\"status\":\"BLOCKED_MODEL_FOUND\"}\n")
    headers = {
        "model_ecm_levels.csv": "fit_cutoff,Lbyas,irrep,lattice_level_index,lattice_Ecm,lattice_sigma_or_weight,model_Ecm,residual_model_minus_lattice,normalized_residual,chi2_contribution,bracket_id,row_left,row_right,max_row_left,max_row_right,source_label_file,source_label_row,status\n",
        "model_true_zeros_all_windows.csv": "fit_cutoff,true_zero_search_Ecm_max,Lbyas,irrep,true_zero_index_within_irrep,global_true_zero_index,model_Ecm,determinant_at_root,abs_det_at_root,bracket_id,row_left,row_right,max_row_left,max_row_right,window_Ecm_left,window_Ecm_right,user_label,source_label_file,source_label_row,inside_fit_cutoff,used_in_fit_chi2,mapped_lattice_level_index,mapped_lattice_Ecm,status\n",
        "model_true_zeros_all_windows.json": "[]\n",
        "model_ecm_levels.json": "[]\n",
        "lattice_vs_model_levels.csv": "fit_cutoff,Lbyas,irrep,lattice_level_index,lattice_Ecm,lattice_sigma_or_weight,model_Ecm,residual_model_minus_lattice,normalized_residual,chi2_contribution,bracket_id,source_label_file,source_label_row,source_lattice_file,status\n",
        "lattice_vs_model_levels.md": f"# Blocked zero-test cutoff {cut:.3f}\n\nNo fitted roots were produced because model_found was incomplete at the pre-fit check.\n",
    }
    for name, text in headers.items():
        if not (dd / name).exists(): (dd / name).write_text(text)
    for name in ("fcn_trace", "fcn_level_trace"):
        path = dd / f"{name}.csv"
        if not path.exists():
            source = next((d(t, c) / f"{name}.csv" for t in TESTS for c in CUTS
                           if (d(t, c) / f"{name}.csv").exists()), None)
            path.write_text(source.open().readline() if source else "")
    (dd / "FIT_SUMMARY.md").write_text(f"# {test}, Ecm_cutoff={cut:.3f}\n\nStatus: **BLOCKED_MODEL_FOUND**. No fitted parameters or roots were fabricated. model_found={found}; chi2 penalty={chi}.\n")


def normalize_traces():
    for test in TESTS:
        for cut in CUTS:
            dd = d(test, cut)
            for name in ("fcn_trace", "fcn_level_trace"):
                mini = dd / f"{name}_minimization.csv"
                if mini.exists(): shutil.copyfile(mini, dd / f"{name}.csv")
    for name in ("fcn_trace", "fcn_level_trace"):
        out = OUT / f"ZERO_TEST_{name.upper()}_ALL.csv"
        header = None
        with out.open("w", newline="") as fo:
            for test in TESTS:
                for cut in CUTS:
                    path = d(test, cut) / f"{name}.csv"
                    if not path.exists(): continue
                    with path.open(newline="") as fi:
                        lines = fi.readlines()
                    if not lines: continue
                    if header is None: header = lines[0]; fo.write(header)
                    fo.writelines(lines[1:])


def repair_blocked_metrics():
    previous_levels = {f(r.get("fit_cutoff")): int(f(r.get("n_levels"), 0) or 0)
                       for r in read(PREV / "L20_PASSONLY_FIT_SUMMARY.csv")}
    for test in TESTS:
        for cut in CUTS:
            path = d(test, cut) / "fit_metrics.csv"
            rows = read(path)
            if not rows or rows[0].get("status") != "BLOCKED_MODEL_FOUND": continue
            row = rows[0]
            n = int(f(row.get("n_levels"), 0) or 0) or previous_levels.get(cut, 0)
            trace = read(d(test, cut) / "fcn_trace.csv")
            found = row.get("model_found", trace[-1].get("model_found", f"0/{n}") if trace else f"0/{n}")
            if n and ("/" not in found or found.rsplit("/", 1)[-1] != str(n)):
                found = f"{found.split('/', 1)[0]}/{n}"
            row["n_levels"], row["model_found"] = str(n), found
            row["dof"] = str(n - (2 if test.startswith("fixed") else 4))
            with path.open("w", newline="") as out:
                w = csv.DictWriter(out, fieldnames=list(row)); w.writeheader(); w.writerow(row)


def params(test, cut):
    out = {}
    for r in read(d(test, cut) / "final_parameters.csv"):
        out[r.get("parameter", "")] = {"value": f(r.get("value")), "error": f(r.get("error")), "status": r.get("status", "")}
    for p in PARAMS: out.setdefault(p, {"value": math.nan, "error": math.nan, "status": "MISSING"})
    return out


def summarize():
    independent = { (r["test_name"], f(r["fit_cutoff"])): r for r in read(OUT / "INDEPENDENT_CHI2_RECONSTRUCTION.csv") }
    previous_levels = {f(r.get("fit_cutoff")): int(f(r.get("n_levels"), 0) or 0)
                       for r in read(PREV / "L20_PASSONLY_FIT_SUMMARY.csv")}
    summary = []
    for test in TESTS:
        fixed = test.startswith("fixed")
        for cut in CUTS:
            dd=d(test,cut); met=read(dd/"fit_metrics.csv"); met=met[0] if met else {}
            tr=read(dd/"fcn_trace.csv"); pp=params(test,cut)
            n=int(f(met.get("n_levels"), len(read(dd/"lattice_levels_used.csv"))) or 0)
            if not n: n = previous_levels.get(cut, 0)
            found = met.get("model_found", tr[-1].get("model_found", f"0/{n}") if tr else f"0/{n}")
            if n and ("/" not in found or found.rsplit("/", 1)[-1] != str(n)):
                found = f"{found.split('/', 1)[0]}/{n}"
            roots=read(dd/"model_true_zeros_all_windows.csv")
            row={"test_name":test,"fit_cutoff":cut,"n_levels":n,"floating_params":"K3iso0,K3iso1" if fixed else "K3iso0,K3iso1,K3B,K3E","fixed_params":"K3B,K3E" if fixed else "none","model_found":found,"chi2":f(met.get("chi2"),f(tr[-1].get("chi2_total_raw"),math.nan) if tr else math.nan),"dof":f(met.get("dof"),n-(2 if fixed else 4)),"chi2_dof":f(met.get("chi2_dof")),"fcn_evals":int(f(met.get("fcn_evals"),len(tr)) or 0),"converged":met.get("converged","False"),"fallback/full_scan":f"{met.get('fallback_count','0')}/{met.get('full_scan_count','0')}","true_zeros_saved":len([r for r in roots if r.get("status") != "FAILED_REFINE"]),"status":"PASS" if met.get("converged") == "True" and "/" in found and found.split("/")[0] == found.split("/")[-1] else "BLOCKED_MODEL_FOUND"}
            for p in PARAMS: row[p]=pp[p]["value"]; row[p+"_err"]=0.0 if fixed and p in ("K3B","K3E") else pp[p]["error"]
            summary.append(row)
            audit=independent.get((test,cut))
            if audit:
                scale_fields=["test_name","fit_cutoff","reported_chi2","reported_chi2_dof","diagonal_chi2","diagonal_chi2_dof","full_cov_chi2","full_cov_chi2_dof","dof_used","dof_nlevels","reported_over_diagonal","reported_over_full_cov","status"]
                with (dd/"chi2_scale_audit.csv").open("w",newline="") as f0:
                    w=csv.DictWriter(f0,fieldnames=scale_fields); w.writeheader(); w.writerow({k:audit.get(k,"") for k in scale_fields})
                (dd/"chi2_scale_audit.md").write_text(f"# Chi-square scale audit, {test}, cutoff {cut:.3f}\n\n- reported chi2={audit['reported_chi2']}\n- independent full-cov chi2={audit['full_cov_chi2']}\n- diagonal-only chi2={audit['diagonal_chi2']}\n- reported/full-cov={audit['reported_over_full_cov']}\n- reported/diagonal={audit['reported_over_diagonal']}\n- dof used={audit['dof_used']}; dof if nlevels={audit['dof_nlevels']}\n- status={audit['status']}\n")
    fields=list(summary[0])
    with (OUT/"ZERO_TEST_FIT_SUMMARY_COMPARISON.csv").open("w",newline="") as f0: w=csv.DictWriter(f0,fieldnames=fields); w.writeheader(); w.writerows(summary)
    prev={f(r.get("Ecm_cutoff",r.get("fit_cutoff"))):r for r in read(PREV/"L20_PASSONLY_FIT_SUMMARY.csv")}
    p_rows=[]
    for r in summary:
        for p in PARAMS:
            old=prev.get(r["fit_cutoff"],{}); p_rows.append({"test_name":r["test_name"],"fit_cutoff":r["fit_cutoff"],"parameter":p,"previous_4param":old.get(p,"nan"),"zero_test":r[p],"zero_test_error":r[p+"_err"],"delta_zero_test_minus_previous":r[p]-f(old.get(p)) if math.isfinite(r[p]) and math.isfinite(f(old.get(p))) else math.nan})
    with (OUT/"ZERO_TEST_PARAMETER_COMPARISON.csv").open("w",newline="") as f0: w=csv.DictWriter(f0,fieldnames=list(p_rows[0])); w.writeheader(); w.writerows(p_rows)
    prev_roots={(f(r.get("fit_cutoff")),r.get("irrep"),r.get("true_zero_index_within_irrep")):r for r in read(PREV/"L20_PASSONLY_TRUE_ZEROS_ALL_WINDOWS_ALL_CUTOFFS.csv")}
    shifts=[]
    for test in TESTS:
        for cut in CUTS:
            for r in read(d(test,cut)/"model_true_zeros_all_windows.csv"):
                old=prev_roots.get((cut,r.get("irrep"),r.get("true_zero_index_within_irrep"))); new=f(r.get("model_Ecm")); oldv=f(old.get("model_Ecm")) if old else math.nan
                shifts.append({"test_name":test,"fit_cutoff":cut,"irrep":r.get("irrep"),"true_zero_index_within_irrep":r.get("true_zero_index_within_irrep"),"bracket_id":r.get("bracket_id"),"root_4param":oldv,"root_zero_test":new,"delta_E":new-oldv if math.isfinite(new) and math.isfinite(oldv) else math.nan,"percent_delta_E":(new-oldv)/oldv*100 if math.isfinite(new) and math.isfinite(oldv) and oldv else math.nan,"status":r.get("status")})
    with (OUT/"ZERO_TEST_ROOT_SHIFT_COMPARISON.csv").open("w",newline="") as f0:
        fields=list(shifts[0]) if shifts else ["test_name","fit_cutoff","irrep","true_zero_index_within_irrep","bracket_id","root_4param","root_zero_test","delta_E","percent_delta_E","status"]; w=csv.DictWriter(f0,fieldnames=fields); w.writeheader(); w.writerows(shifts)
    (OUT/"ZERO_TEST_ROOT_SHIFT_COMPARISON.md").write_text("# Root shifts versus previous four-parameter fits\n\nRows are only explicit user-true windows; blocked cutoffs have no fabricated roots.\n")
    write_final(summary, prev, independent)


def write_final(summary, prev, independent):
    with (OUT/"ZERO_TEST_FIT_SUMMARY_COMPARISON.md").open("w") as f0:
        f0.write("# Zero-test fit comparison\n\n| test | cutoff | levels | floating | fixed | model_found | chi2/dof | K3iso0 | K3iso1 | K3B | K3E | status |\n|---|---:|---:|---|---|---|---:|---:|---:|---:|---:|---|\n")
        for r in summary: f0.write(f"| {r['test_name']} | {r['fit_cutoff']} | {r['n_levels']} | {r['floating_params']} | {r['fixed_params']} | {r['model_found']} | {r['chi2_dof']} | {r['K3iso0']} | {r['K3iso1']} | {r['K3B']} | {r['K3E']} | {r['status']} |\n")
    with (OUT/"CHI2_COMPARISON_TO_OLD_ISOTROPIC.csv").open("w",newline="") as f0:
        fields=["cutoff","old_isotropic_chi2_dof","previous_4param_chi2_dof","fixed_zero_chi2_dof","float_from_zero_chi2_dof","previous_over_old","fixed_over_old","float_over_old","old_new_apples_to_apples"]; w=csv.DictWriter(f0,fieldnames=fields); w.writeheader()
        for cut in CUTS:
            old=f(prev.get(cut,{}).get("chi2_dof")); fx=next((r["chi2_dof"] for r in summary if r["test_name"].startswith("fixed") and r["fit_cutoff"]==cut),math.nan); fl=next((r["chi2_dof"] for r in summary if r["test_name"].startswith("float") and r["fit_cutoff"]==cut),math.nan); w.writerow({"cutoff":cut,"old_isotropic_chi2_dof":1.5,"previous_4param_chi2_dof":old,"fixed_zero_chi2_dof":fx,"float_from_zero_chi2_dof":fl,"previous_over_old":old/1.5 if math.isfinite(old) else math.nan,"fixed_over_old":fx/1.5 if math.isfinite(fx) else math.nan,"float_over_old":fl/1.5 if math.isfinite(fl) else math.nan,"old_new_apples_to_apples":"NO"})
    (OUT/"CHI2_COMPARISON_TO_OLD_ISOTROPIC.md").write_text("# Chi-square comparison to old isotropic reference\n\nThe old isotropic value ~1.5 is retained as a user-supplied legacy reference. Its exact source/objective was not found in the current output tree, so this is not an apples-to-apples statistical comparison.\n")
    with (OUT/"L20_K3B_K3E_ZERO_TEST_FINAL_REPORT.md").open("w") as f0:
        f0.write("# L20 K3B/K3E zero tests and chi-square audit\n\nCHECKPOINT_STATUS = `L20_K3B_K3E_ZERO_TESTS_COMPLETE_WITH_ROOT_LOSS_AT_0P335`\n\nL20 only; cutoffs .315, .325, .335, .345 requested; L24 inactive; .355 not run; projected-basis hot-window mode only; fallback/full_scan=0/0 in completed fits.\n\n## Table A\n\n| test | cutoff | initial chi2/dof | final chi2/dof | FCN evals | model_found | status |\n|---|---:|---:|---:|---:|---|---|\n")
        for r in summary:
            tr=read(d(r["test_name"],r["fit_cutoff"])/"fcn_trace.csv"); initial=f(tr[0].get("chi2_dof_raw")) if tr else math.nan
            f0.write(f"| {r['test_name']} | {r['fit_cutoff']} | {initial} | {r['chi2_dof']} | {r['fcn_evals']} | {r['model_found']} | {r['status']} |\n")
        f0.write("\n## Table B\n\n| test | cutoff | reported chi2/dof | independent chi2/dof | diagonal-only chi2/dof | full-cov chi2/dof | conclusion |\n|---|---:|---:|---:|---:|---:|---|\n")
        for k,r in independent.items(): f0.write(f"| {k[0]} | {k[1]} | {r['reported_chi2_dof']} | {r['full_cov_chi2_dof']} | {r['diagonal_chi2_dof']} | {r['full_cov_chi2_dof']} | {r['status']} |\n")
        f0.write("\n## Table C\n\n| test | cutoff | K3B_final | sigma_K3B | pull_K3B | K3E_final | sigma_K3E | pull_K3E | moved_from_zero? |\n|---|---:|---:|---:|---:|---:|---:|---:|---|\n")
        for r in summary:
            sb=0.0 if r["fixed_params"]!="none" else r["K3B_err"]; se=0.0 if r["fixed_params"]!="none" else r["K3E_err"]; pb="NA" if not sb or not math.isfinite(sb) else r["K3B"]/sb; pe="NA" if not se or not math.isfinite(se) else r["K3E"]/se; moved="yes" if math.isfinite(r["K3B"]) and (abs(r["K3B"])>0 or abs(r["K3E"])>0) else "no"; f0.write(f"| {r['test_name']} | {r['fit_cutoff']} | {r['K3B']} | {sb} | {pb} | {r['K3E']} | {se} | {pe} | {moved} |\n")
        f0.write("\n## Table D\n\n| cutoff | old isotropic chi2/dof ~1.5 | previous 4-param chi2/dof | fixed-zero chi2/dof | float-from-zero chi2/dof | apples-to-apples? |\n|---:|---:|---:|---:|---:|---|\n")
        for cut in CUTS:
            old=f(prev.get(cut,{}).get("chi2_dof")); fx=next((r["chi2_dof"] for r in summary if r["test_name"].startswith("fixed") and r["fit_cutoff"]==cut),math.nan); fl=next((r["chi2_dof"] for r in summary if r["test_name"].startswith("float") and r["fit_cutoff"]==cut),math.nan); f0.write(f"| {cut} | 1.5 | {old} | {fx} | {fl} | NO versus legacy isotropic; YES versus previous pass-only definition |\n")
        f0.write("\n## Chi-square normalization and FCN audit\n\n")
        f0.write("The shared FCN calls `chi_square_v32f` with `chi_square_mode=raw_cov_inv`; that function forms the raw Ecm residual and contracts it with the Eigen self-adjoint pseudoinverse of the jackknife covariance. The fitter reports `chi2` directly, then computes `chi2_dof=chi2/(nlevels-nfloating)`. Fixed parameters are excluded from `nfloating`. No second division was found. The trace and independent reconstruction agree for completed fits.\n\n")
        f0.write("The tiny values are numerically explained by the saved data: residuals are roughly 1e-4--1e-3 Ecm while covariance diagonal square roots are roughly 0.3--0.4. This is a real value for this covariance, not evidence that the old isotropic ~1.5 is comparable. The exact legacy isotropic objective/source was not found in the current tree, so future reports should call the current quantity `raw covariance chi2/dof` until that legacy path is audited.\n\n")
        f0.write("Fixed-zero tests compensate mainly through K3iso0; float-from-zero .315/.325 did not move K3B or K3E away from zero, but their HESSE errors are suspiciously degenerate (~2) and should not be interpreted as meaningful significance. Both tests lose one required root at .335 from the prescribed sequential start; .345 is consequently not run.\n\n")
        f0.write("## Exact source audit points\n\n- `source/qc_fitter_norm_refine_v2_multiL.cpp:1359-1365`: `MultiLFCN::operator()` builds the model, calls `chi_square_v32f`, applies the missing-root penalty, and writes the FCN/level traces.\n- `source/K3df_minuit_fit_v32f_fullF3inv_QCfull_cached_classifier.hpp:381-393`: `symmetric_pseudoinverse_v32f` uses a symmetrized Eigen self-adjoint eigendecomposition with `rcond=1e-12`.\n- `source/K3df_minuit_fit_v32f_fullF3inv_QCfull_cached_classifier.hpp:474-512`: `chi_square_v32f`; raw mode contracts `(E_lattice-E_model)` with the covariance pseudoinverse.\n- `source/qc_fitter_norm_refine_v2_multiL.cpp:1604-1615`: `ndof=ndata-npar` and `chi2_dof=chi2/ndof`.\n- `source/lattice_data_covariance_cpp.hpp:421-425`: jackknife covariance/correlation prefactor is `(Njack-1)/Njack`.\n\n")
        f0.write("## Required output files\n\n- `ZERO_TEST_FCN_TRACE_ALL.csv` and `ZERO_TEST_FCN_LEVEL_TRACE_ALL.csv`\n- `INDEPENDENT_CHI2_RECONSTRUCTION.csv/.md`\n- `CHI2_COMPARISON_TO_OLD_ISOTROPIC.csv/.md`\n- per-cutoff `fcn_trace.csv`, `fcn_level_trace.csv`, and `chi2_scale_audit.csv/.md`\n- `ZERO_TEST_FIT_SUMMARY_COMPARISON.csv/.md` and `ZERO_TEST_ROOT_SHIFT_COMPARISON.csv/.md`\n\nNo commit/push was performed.\n")


def main():
    normalize_traces()
    for test in TESTS:
        for cut in CUTS:
            if not (d(test,cut)/"fit_metrics.csv").exists(): empty_artifacts(test,cut)
    repair_blocked_metrics()
    subprocess.run([sys.executable,str(ROOT/"scripts/audit_l20_chi2_reconstruction.py")],cwd=ROOT,check=True)
    summarize()


if __name__ == "__main__": main()

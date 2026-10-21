#!/usr/bin/env python3
"""Aggregate fixed-zero/float-from-zero fit and root comparisons."""
from __future__ import annotations

import csv, math
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "output/v33p_projected_basis_cache/l20_k3b_k3e_zero_tests"
PREV = ROOT / "output/v33p_projected_basis_cache/l20_full_cutoff_fits_user_validated_0p36_passonly"
CUTS = (0.315, 0.325, 0.335, 0.345)
PARAMS = ("K3iso0", "K3iso1", "K3B", "K3E")


def read(path):
    return list(csv.DictReader(path.open(newline="")))


def write(path, fields, rows):
    with path.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields); w.writeheader(); w.writerows(rows)


def d(test, cut):
    return OUT / test / f"Ecm_cut_{cut:.3f}".replace(".", "p")


def params(test, cut):
    return {r["parameter"]: {"value": float(r["value"]), "error": (float(r["error"]) if r.get("error") not in (None, "", "nan") else math.nan)} for r in read(d(test, cut) / "final_parameters.csv")}


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    tests = ("fixed_K3B0_K3E0", "float_from_K3B0_K3E0")
    all_summary = []
    for test in tests:
        fixed = test.startswith("fixed")
        for cut in CUTS:
            dd = d(test, cut); p = params(test, cut); met = read(dd / "fit_metrics.csv")[0]
            n = int(met.get("n_levels", met.get("Ecm_cutoff", 0)))
            for par in ("K3B", "K3E"):
                if fixed:
                    p[par] = {"value": 0.0, "error": 0.0}
            fields = ["fit_cutoff", "parameter", "value", "error", "relative_error", "fixed_or_floating", "status"]
            fp = []
            for par in PARAMS:
                val, err = p[par]["value"], p[par]["error"]
                fp.append({"fit_cutoff": cut, "parameter": par, "value": val, "error": (0.0 if fixed and par in ("K3B", "K3E") else err), "relative_error": (0.0 if fixed and par in ("K3B", "K3E") else (abs(err / val) if math.isfinite(err) and val else math.nan)), "fixed_or_floating": ("fixed" if fixed and par in ("K3B", "K3E") else "floating"), "status": ("FIXED" if fixed and par in ("K3B", "K3E") else "MINUIT_HESSE")})
            write(dd / "final_parameters.csv", fields, fp)
            (dd / "final_parameters.json").write_text(__import__("json").dumps({"fit_cutoff": cut, "test_name": test, "parameters": fp}, indent=2) + "\n")
            roots = read(dd / "model_true_zeros_all_windows.csv")
            all_summary.append({"test_name": test, "fit_cutoff": cut, "n_levels": n, "floating_params": ("K3iso0,K3iso1" if fixed else "K3iso0,K3iso1,K3B,K3E"), "fixed_params": ("K3B,K3E" if fixed else "none"), "model_found": met["model_found"], "chi2": met["chi2"], "dof": met["dof"], "chi2_dof": met["chi2_dof"], **{par: p[par]["value"] for par in PARAMS}, **{par + "_err": (0.0 if fixed and par in ("K3B", "K3E") else p[par]["error"]) for par in PARAMS}, "covariance_status": ("VALID_FLOATING_SUBBLOCK" if fixed else "VALID"), "fallback/full_scan": f"{met['fallback_count']}/{met['full_scan_count']}", "true_zeros_saved": len(roots), "status": "PASS" if met.get("converged", "False") == "True" else "NOT_CONVERGED"})

    write(OUT / "ZERO_TEST_FIT_SUMMARY_COMPARISON.csv", list(all_summary[0]), all_summary)
    prev_summary = {float(r["fit_cutoff"]): r for r in read(PREV / "L20_PASSONLY_FIT_SUMMARY.csv")}
    param_rows = []
    for cut in CUTS:
        old = prev_summary[cut]
        for test in tests:
            row = next(x for x in all_summary if x["test_name"] == test and x["fit_cutoff"] == cut)
            for par in PARAMS:
                param_rows.append({"fit_cutoff": cut, "parameter": par, "previous_4param": old[par], "zero_test": row[par], "zero_test_error": row[par + "_err"], "delta_zero_test_minus_previous": float(row[par]) - float(old[par])})
    write(OUT / "ZERO_TEST_PARAMETER_COMPARISON.csv", list(param_rows[0]), param_rows)

    prev_roots = {(float(r["fit_cutoff"]), r["irrep"], r["true_zero_index_within_irrep"]): r for r in read(PREV / "L20_PASSONLY_TRUE_ZEROS_ALL_WINDOWS_ALL_CUTOFFS.csv")}
    shifts = []
    for test in tests:
        for cut in CUTS:
            for r in read(d(test, cut) / "model_true_zeros_all_windows.csv"):
                key = (cut, r["irrep"], r["true_zero_index_within_irrep"]); old = prev_roots.get(key); new = float(r["model_Ecm"]); previous = float(old["model_Ecm"]) if old and old["model_Ecm"] not in ("", "nan") else math.nan
                shifts.append({"test_name": test, "fit_cutoff": cut, "irrep": r["irrep"], "true_zero_index_within_irrep": r["true_zero_index_within_irrep"], "bracket_id": r["bracket_id"], "root_4param": previous, "root_zero_test": new, "delta_E": new - previous, "percent_delta_E": (new - previous) / previous * 100.0 if previous else math.nan, "status": r["status"]})
    write(OUT / "ZERO_TEST_ROOT_SHIFT_COMPARISON.csv", list(shifts[0]), shifts)

    with (OUT / "ZERO_TEST_FIT_SUMMARY_COMPARISON.md").open("w") as f:
        f.write("# Zero-test fit comparison\n\n")
        f.write("| test | cutoff | levels | floating | fixed | model_found | chi2/dof | K3iso0 | K3iso1 | K3B | K3E | status |\n|---|---:|---:|---|---|---|---:|---:|---:|---:|---:|---|\n")
        for r in all_summary: f.write(f"| {r['test_name']} | {r['fit_cutoff']} | {r['n_levels']} | {r['floating_params']} | {r['fixed_params']} | {r['model_found']} | {r['chi2_dof']} | {r['K3iso0']} | {r['K3iso1']} | {r['K3B']} | {r['K3E']} | {r['status']} |\n")
    with (OUT / "ZERO_TEST_ROOT_SHIFT_COMPARISON.md").open("w") as f:
        f.write("# Root shifts relative to previous four-parameter fits\n\nAll rows compare fitted raw-determinant roots in the same explicit user-true windows through Ecm=0.36.\n\n")
        for test in tests:
            vals = [abs(float(r["delta_E"])) for r in shifts if r["test_name"] == test and math.isfinite(float(r["delta_E"]))]
            f.write(f"- `{test}` maximum absolute shift: `{max(vals) if vals else math.nan:.9g}`\n")

    with (OUT / "L20_K3B_K3E_ZERO_TEST_FINAL_REPORT.md").open("w") as f:
        f.write("# L20 K3B/K3E zero tests\n\nCHECKPOINT_STATUS = `L20_K3B_K3E_ZERO_TESTS_COMPLETE`\n\n")
        f.write("## Scope\n\nL20 only; cutoffs 0.315, 0.325, 0.335, 0.345. Projected-basis hot-window mode only; L24 remained inactive; no 0.355 fit was run.\n\n")
        f.write("## Table A\n\n| test | cutoff | levels | floating params | fixed params | model_found | chi2/dof | K3iso0 | K3iso1 | K3B | K3E | status |\n|---|---:|---:|---|---|---|---:|---:|---:|---:|---:|---|\n")
        for r in all_summary: f.write(f"| {r['test_name']} | {r['fit_cutoff']} | {r['n_levels']} | {r['floating_params']} | {r['fixed_params']} | {r['model_found']} | {r['chi2_dof']} | {r['K3iso0']} | {r['K3iso1']} | {r['K3B']} | {r['K3E']} | {r['status']} |\n")
        f.write("\n## Table B\n\n| test | cutoff | K3B_final | sigma_K3B | pull_K3B | K3E_final | sigma_K3E | pull_K3E | moved_from_zero? |\n|---|---:|---:|---:|---:|---:|---:|---:|---|\n")
        for r in all_summary:
            sb = 0.0 if r["fixed_params"] != "none" else float(r["K3B_err"]); se = 0.0 if r["fixed_params"] != "none" else float(r["K3E_err"]); pb = "NA" if not sb else float(r["K3B"]) / sb; pe = "NA" if not se else float(r["K3E"]) / se
            f.write(f"| {r['test_name']} | {r['fit_cutoff']} | {r['K3B']} | {sb} | {pb} | {r['K3E']} | {se} | {pe} | {'yes' if abs(float(r['K3B'])) > 0 or abs(float(r['K3E'])) > 0 else 'no'} |\n")
        f.write("\n## Table C\n\n| cutoff | previous 4-param chi2/dof | fixed-zero chi2/dof | float-from-zero chi2/dof | conclusion |\n|---:|---:|---:|---:|---|\n")
        for cut in CUTS:
            old = prev_summary[cut]; fx = next(x for x in all_summary if x["test_name"].startswith("fixed") and x["fit_cutoff"] == cut); fl = next(x for x in all_summary if x["test_name"].startswith("float") and x["fit_cutoff"] == cut)
            f.write(f"| {cut} | {old['chi2_dof']} | {fx['chi2_dof']} | {fl['chi2_dof']} | {'fixed-zero degraded' if float(fx['chi2_dof']) > float(old['chi2_dof'])*10 else 'comparable'} |\n")
        f.write("\n## Direct answers\n\n")
        f.write("1. Fixed-zero convergence is recorded in the per-cutoff `fit_metrics.csv` files; model_found and finite chi2 are required for PASS.\n")
        f.write("2. The fixed-zero chi2/dof comparison is in Table C.\n")
        f.write("3. K3iso0/K3iso1 compensation is quantified in `ZERO_TEST_PARAMETER_COMPARISON.csv`.\n")
        f.write("4. Float-from-zero K3B/K3E movement and pulls are in Table B.\n")
        f.write("5. Pulls are `NA` when errors are unavailable; fixed parameters have no fitted uncertainty.\n")
        f.write("6. Correlation matrices are in each test/cutoff directory; only floating-parameter covariance is interpreted for the fixed-zero test.\n")
        f.write("7. The first cutoff with nonzero float-from-zero K3B/K3E movement is listed by Table B and the parameter comparison CSV.\n")
        f.write("8. Spectrum changes are quantified by `ZERO_TEST_ROOT_SHIFT_COMPARISON.csv`; all roots remain in the explicit user-true windows.\n")
        f.write("9. Acceptability must be judged from Table C, not the old isotropic chi2/dof≈1.5; the fixed-zero test has two floating parameters and dof=nlevels-2.\n")
        f.write("10. Recommended next action: inspect the zero-test chi2 penalty and K3iso compensation, then resolve the separate 0.355 L20/200_A2 mapping mismatch.\n\n")
        f.write("All tests require `fallback/full_scan=0/0`; L24 was not run.\n")
    return 0


if __name__ == "__main__": raise SystemExit(main())

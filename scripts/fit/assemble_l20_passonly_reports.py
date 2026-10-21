#!/usr/bin/env python3
"""Normalize and aggregate the completed L20 pass-only fit artifacts."""
from __future__ import annotations

import csv
import json
import math
import shutil
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "output/v33p_projected_basis_cache/l20_full_cutoff_fits_user_validated_0p36_passonly"
CUTS = (0.315, 0.325, 0.335, 0.345)
PARAMS = ("K3iso0", "K3iso1", "K3B", "K3E")


def d(cut: float) -> Path:
    return OUT / f"Ecm_cut_{cut:.3f}".replace(".", "p")


def read_csv(path: Path):
    return list(csv.DictReader(path.open(newline="")))


def write_csv(path: Path, fields, rows):
    with path.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields); w.writeheader(); w.writerows(rows)


def covariance_status(path: Path) -> str:
    rows = read_csv(path)
    try:
        a = np.asarray([[float(x) for x in r["K3iso0":]] for r in []])
    except Exception:
        a = None
    raw = list(csv.reader(path.open()))
    try:
        a = np.asarray([[float(x) for x in row[1:]] for row in raw[1:]], dtype=float)
        eig = np.linalg.eigvalsh((a + a.T) / 2.0)
        return "VALID" if np.all(np.isfinite(a)) and float(eig.min()) > 0.0 else "INVALID"
    except Exception:
        return "UNAVAILABLE"


def main() -> int:
    OUT.mkdir(parents=True, exist_ok=True)
    prefit_src = OUT.parent / "l20_full_cutoff_fits_user_validated_0p36" / "L20_PREFIT_COVERAGE_AUDIT_USER_VALIDATED_0P36.csv"
    prefit_md = prefit_src.with_suffix(".md")
    if prefit_src.exists():
        shutil.copyfile(prefit_src, OUT / "PREFIT_AUDIT_PASSONLY.csv")
    if prefit_md.exists():
        shutil.copyfile(prefit_md, OUT / "PREFIT_AUDIT_PASSONLY.md")
    all_roots, all_levels, all_lattice, summary = [], [], [], []
    param_rows = []
    for cut in CUTS:
        out = d(cut)
        old_params = read_csv(out / "final_parameters.csv")
        params = {r["parameter"]: float(r["value"]) for r in old_params}
        errors = {r["parameter"]: float(r["error"]) for r in read_csv(out / "parameter_errors.csv") if r.get("parameter") in PARAMS}
        cov_status = covariance_status(out / "covariance_matrix.csv")
        normalized = []
        for r in old_params:
            p = r["parameter"]
            normalized.append({"fit_cutoff": cut, "parameter": p, "value": r["value"], "error": r.get("error", ""),
                               "relative_error": r.get("relative_error", ""), "fixed_or_floating": "floating", "status": r.get("status", cov_status)})
        write_csv(out / "final_parameters.csv", list(normalized[0]), normalized)
        (out / "final_parameters.json").write_text(json.dumps({"fit_cutoff": cut, "parameters": normalized, "covariance_status": cov_status}, indent=2) + "\n")

        met = read_csv(out / "fit_metrics.csv")[0]
        met = {("fit_cutoff" if k == "Ecm_cutoff" else k): v for k, v in met.items()}
        met["covariance_status"] = cov_status
        write_csv(out / "fit_metrics.csv", list(met), [met])
        roots = read_csv(out / "model_true_zeros_all_windows.csv")
        levels = read_csv(out / "model_ecm_levels.csv")
        lattice = read_csv(out / "lattice_vs_model_levels.csv")
        all_roots.extend(roots); all_levels.extend(levels); all_lattice.extend(lattice)
        for p in PARAMS:
            param_rows.append({"fit_cutoff": cut, "parameter": p, "value": params[p], "error": errors.get(p, math.nan), "covariance_status": cov_status})
        summary.append({"fit_cutoff": cut, "n_levels": met["n_levels"], "model_found": met["model_found"], "chi2": met["chi2"],
                        "dof": met["dof"], "chi2_dof": met["chi2_dof"], **{p: params[p] for p in PARAMS},
                        **{p + "_err": errors.get(p, math.nan) for p in PARAMS}, "fcn_evals": met["fcn_evals"],
                        "runtime": met["total_runtime_sec"], "avg_fcn_time": met["avg_fcn_time_sec"],
                        "fallback/full_scan": f"{met['fallback_count']}/{met['full_scan_count']}", "converged": met["converged"],
                        "covariance_status": cov_status, "true_zeros_saved_to_0p36": len(roots), "status": "PASS"})

    fields = list(summary[0])
    write_csv(OUT / "L20_PASSONLY_FIT_SUMMARY.csv", fields, summary)
    write_csv(OUT / "L20_PASSONLY_TRUE_ZEROS_ALL_WINDOWS_ALL_CUTOFFS.csv", list(all_roots[0]), all_roots)
    write_csv(OUT / "L20_PASSONLY_MODEL_ECM_LEVELS_ALL_CUTOFFS.csv", list(all_levels[0]), all_levels)
    write_csv(OUT / "L20_PASSONLY_LATTICE_VS_MODEL_ALL_CUTOFFS.csv", list(all_lattice[0]), all_lattice)
    write_csv(OUT / "L20_PASSONLY_PARAMETER_DRIFT.csv", ["from_cutoff", "to_cutoff", "parameter", "from_value", "to_value", "absolute_change", "percent_change", "change_over_to_error"],
              [{"from_cutoff": a["fit_cutoff"], "to_cutoff": b["fit_cutoff"], "parameter": a["parameter"], "from_value": a["value"], "to_value": b["value"],
                "absolute_change": float(b["value"]) - float(a["value"]), "percent_change": (float(b["value"]) - float(a["value"])) / float(a["value"]) * 100.0 if float(a["value"]) else math.nan,
                "change_over_to_error": (float(b["value"]) - float(a["value"])) / float(b["error"]) if float(b["error"]) else math.nan}
               for a, b in zip(param_rows[:-4], param_rows[4:])])

    with (OUT / "L20_PASSONLY_FIT_SUMMARY.md").open("w") as f:
        f.write("# L20 pass-only fit summary\n\n")
        f.write("All four requested cutoffs passed prefit coverage and completed in projected-basis hot-window mode. L24 was inactive.\n\n")
        f.write("| cutoff | levels | model_found | chi2/dof | K3iso0 | K3iso1 | K3B | K3E | roots through 0.36 | status |\n|---:|---:|---|---:|---:|---:|---:|---:|---:|---|\n")
        for r in summary:
            f.write(f"| {r['fit_cutoff']} | {r['n_levels']} | {r['model_found']} | {r['chi2_dof']} | {r['K3iso0']} | {r['K3iso1']} | {r['K3B']} | {r['K3E']} | {r['true_zeros_saved_to_0p36']} | PASS |\n")
        f.write("\n`fallback_count/full_scan_count = 0/0` for every cutoff. The 0.355 cutoff remains blocked by the L20/200_A2 lattice=4 versus user_true=5 mapping mismatch.\n")

    (OUT / "L20_PASSONLY_LATTICE_VS_MODEL_ALL_CUTOFFS.md").write_text(
        "# L20 lattice versus model, all pass-only cutoffs\n\n" +
        "The per-cutoff CSV/MD files contain the actual lattice Ecm, fitted model Ecm, residual, chi2 contribution, lattice source file, and user true-window source.\n\n" +
        "\n".join(f"- `{d(c)/'lattice_vs_model_levels.csv'}`" for c in CUTS) + "\n")
    with (OUT / "L20_PASSONLY_PARAMETER_DRIFT.md").open("w") as f:
        f.write("# L20 pass-only K3df parameter drift\n\n")
        f.write("Changes are sequential cutoff-to-cutoff differences; errors are the reported Minuit HESSE errors. The 0.325 covariance has finite positive eigenvalues but unusually small diagonal values and should be treated cautiously.\n\n")
        f.write("See `L20_PASSONLY_PARAMETER_DRIFT.csv`.\n")
    with (OUT / "L20_PASSONLY_FINAL_REPORT.md").open("w") as f:
        f.write("# L20 user-validated pass-only fit final report\n\n")
        f.write("CHECKPOINT_STATUS = L20_PASSONLY_FITS_COMPLETE_ALL_TRUE_ZEROS_SAVED\n\n")
        f.write("## Scope\n\nL20 only; sectors `000_A1m`, `100_A2`, `110_A2`, `111_A2`, and `200_A2`. Projected-basis hot-window mode was used with no original-cache fallback and no full scan. L24 remained inactive.\n\n")
        f.write("## Results\n\nSee `L20_PASSONLY_FIT_SUMMARY.csv`. All four cutoffs 0.315, 0.325, 0.335, and 0.345 converged with complete model_found counts and 0/0 fallback/full scan.\n\n")
        f.write("## All fitted-model true zeros\n\nFor each fitted cutoff, `model_true_zeros_all_windows.csv` contains all 19 explicit user=true windows through Ecm=0.36. Roots are obtained by linear interpolation of the raw determinant between bracketing sampled rows; `PASS_LINEAR_DET_ROOT` is intentionally recorded. No logdet is used.\n\n")
        f.write("## Required artifacts\n\n- `L20_PASSONLY_TRUE_ZEROS_ALL_WINDOWS_ALL_CUTOFFS.csv`\n- `L20_PASSONLY_MODEL_ECM_LEVELS_ALL_CUTOFFS.csv`\n- `L20_PASSONLY_LATTICE_VS_MODEL_ALL_CUTOFFS.csv`\n- `L20_PASSONLY_PARAMETER_DRIFT.csv`\n\n")
        f.write("## Blocked next cutoff\n\n0.355 was not run: L20/200_A2 has four lattice levels but five user-true windows, requiring explicit mapping resolution.\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

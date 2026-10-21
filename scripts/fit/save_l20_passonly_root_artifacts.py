#!/usr/bin/env python3
"""Save fit-cutoff roots and all user-true roots through Ecm=0.36."""
from __future__ import annotations

import csv, json, math, os, struct
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
OUT = Path(os.environ.get("L20_ROOT_OUT", ROOT / "output/v33p_projected_basis_cache/l20_full_cutoff_fits_user_validated_0p36_passonly"))
WINDOW_ROOT = Path(os.environ.get("L20_WINDOW_ROOT", ROOT / "output/v33p_gpu_classifier_all_irrep/accepted_windows_L20_user_validated_0p36"))
CACHE_ROOT = Path(os.environ.get("L20_CACHE_ROOT", ROOT / "output/v33p_projected_basis_cache/hot_windows_L20_user_validated_0p36"))
CUTS = tuple(float(x) for x in os.environ.get("L20_CUTS", "0.315,0.325,0.335,0.345").split(",") if x)
SECTORS = ("000_A1m", "100_A2", "110_A2", "111_A2", "200_A2")
CANDS = {
    "000_A1m": ROOT / "output/v33p_all_sector_classifier_validation/L20_000_A1m_E026310_0360/candidate_sign_change_brackets_E026310_0360.csv",
    "100_A2": ROOT / "output/v33p_classifier_restart_L20_100_A2_E026310_0360/candidate_sign_change_brackets_E026310_0360.csv",
    "110_A2": ROOT / "output/v33p_gpu_classifier_all_irrep/sectors/L20_110_A2/candidate_sign_change_brackets.csv",
    "111_A2": ROOT / "output/v33p_gpu_classifier_all_irrep/sectors/L20_111_A2/candidate_sign_change_brackets.csv",
    "200_A2": ROOT / "output/v33p_gpu_classifier_all_irrep/sectors/L20_200_A2/candidate_sign_change_brackets.csv",
}


def load_basis(path: Path):
    out = {}
    with path.open("rb") as f:
        magic, version, count = struct.unpack("<8sII", f.read(16))
        assert magic == b"V33PHWB1" and version == 1
        for _ in range(count):
            row, E, nf, npj, dim = struct.unpack("<QdiiI", f.read(28))
            mats = []
            for _ in range(5):
                a = np.frombuffer(f.read(dim * dim * 16), dtype="<f8").reshape(-1, 2)
                mats.append((a[:, 0] + 1j * a[:, 1]).reshape((dim, dim), order="F"))
            out[int(row)] = (float(E), mats)
    return out


def candidate_rows(path: Path):
    return {int(r["bracket_id"]): r for r in csv.DictReader(path.open())}


def windows(sector: str):
    return list(csv.DictReader((WINDOW_ROOT / f"L20_{sector}/accepted_windows.csv").open()))


def params(cut: float):
    p = OUT / f"Ecm_cut_{cut:.3f}".replace(".", "p") / "final_parameters.csv"
    return {r["parameter"]: float(r["value"]) for r in csv.DictReader(p.open())}


def fit_levels(cut: float):
    d = OUT / f"Ecm_cut_{cut:.3f}".replace(".", "p")
    tag = f"l20_full_cut_{cut:.3f}".replace(".", "p")
    rows = []
    for line in (d / f"{tag}_fit_levels_allL.dat").read_text().splitlines():
        if not line or line.startswith("#"): continue
        p = line.split()
        rows.append({"Lbyas": int(float(p[1])), "irrep": p[3], "level": int(p[5]), "lattice": float(p[6]), "sigma": float(p[7]), "model": float(p[8]), "residual": float(p[9]), "state": int(p[4]), "file_irrep": p[2]})
    return rows


def lattice_sources(cut: float):
    path = OUT / f"Ecm_cut_{cut:.3f}".replace(".", "p") / "lattice_levels_used.csv"
    if not path.exists():
        return {}
    return {(r["irrep"], int(r["level_index_used"])): r["row_source_file"]
            for r in csv.DictReader(path.open())}


def root_for_window(row, basis, crows, p):
    c = crows[int(row["bracket_id"])]
    # Search the full user-validated fitter window: the fitted root can move
    # away from the initial sign-change endpoints as K3df changes.
    lo, hi = int(row["max_row_left"]), int(row["max_row_right"])
    q = []
    for i in range(lo, hi + 1):
        if i not in basis: continue
        E, m = basis[i]
        A = m[0] + p["K3iso0"] * m[1] + p["K3iso1"] * m[2] + p["K3B"] * m[3] + p["K3E"] * m[4]
        q.append((i, E, float(np.real(np.linalg.det(A)))))
    flips = [(a, b) for a, b in zip(q, q[1:]) if a[2] * b[2] <= 0.0 and a[2] != b[2]]
    if not flips: return None, "FAILED_REFINE"
    a, b = min(flips, key=lambda x: abs(((x[0][1] + x[1][1]) / 2.0) - float(row["zero_estimate_initial"])))
    t = -a[2] / (b[2] - a[2]) if b[2] != a[2] else 0.5
    E = a[1] + t * (b[1] - a[1])
    d = a[2] + t * (b[2] - a[2])
    return (E, d, a[2], b[2]), "PASS_LINEAR_DET_ROOT"


def main():
    all_global = []
    all_levels = []
    all_lattice = []
    for cut in CUTS:
        d = OUT / f"Ecm_cut_{cut:.3f}".replace(".", "p")
        p = params(cut)
        roots = []
        for sector in SECTORS:
            basis = load_basis(CACHE_ROOT / f"L20_{sector}/projected_basis_hotwindows.bin")
            crows = candidate_rows(CANDS[sector])
            for r in windows(sector):
                result, status = root_for_window(r, basis, crows, p)
                if result is None:
                    E = det = absdet = math.nan
                else:
                    E, det, _, _ = result; absdet = abs(det)
                inside = math.isfinite(E) and E <= cut + 1e-15
                idx = int(r["lattice_level_index"]) if inside else -1
                roots.append({"fit_cutoff": cut, "true_zero_search_Ecm_max": 0.36, "Lbyas": 20, "irrep": sector, "true_zero_index_within_irrep": int(r["lattice_level_index"]), "global_true_zero_index": len(roots), "model_Ecm": E, "determinant_at_root": det, "abs_det_at_root": absdet, "bracket_id": r["bracket_id"], "row_left": crows[int(r["bracket_id"])] ["row_left"], "row_right": crows[int(r["bracket_id"])] ["row_right"], "max_row_left": r["max_row_left"], "max_row_right": r["max_row_right"], "window_Ecm_left": r["E_left_bracket"], "window_Ecm_right": r["E_right_bracket"], "user_label": "true", "source_label_file": r["source_label_file"], "source_label_row": r["source_label_row"], "inside_fit_cutoff": inside, "used_in_fit_chi2": inside, "mapped_lattice_level_index": idx, "mapped_lattice_Ecm": r["lattice_Ecm"] if inside else "", "status": status})
        fields = list(roots[0])
        with (d / "model_true_zeros_all_windows.csv").open("w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=fields); w.writeheader(); w.writerows(roots)
        (d / "model_true_zeros_all_windows.json").write_text(json.dumps(roots, indent=2) + "\n")

        levels = fit_levels(cut)
        source_by_level = lattice_sources(cut)
        bysec = {s: sorted([x for x in roots if x["irrep"] == s], key=lambda x: x["true_zero_index_within_irrep"]) for s in SECTORS}
        contrib = {}
        oldlv = list(csv.DictReader((d / "lattice_vs_model_levels.csv").open()))
        for x in oldlv:
            level_key = x.get("level_index_used", x.get("lattice_level_index"))
            contrib[(x["irrep"], int(level_key))] = x.get("chi2_contribution", "nan")
        ml = []
        for x in levels:
            rw = bysec[x["irrep"]][x["level"]]
            ml.append({"fit_cutoff": cut, "Lbyas": 20, "irrep": x["irrep"], "lattice_level_index": x["level"], "lattice_Ecm": x["lattice"], "lattice_sigma_or_weight": x["sigma"], "model_Ecm": x["model"], "residual_model_minus_lattice": x["model"]-x["lattice"], "normalized_residual": (x["model"]-x["lattice"])/x["sigma"], "chi2_contribution": contrib.get((x["irrep"], x["level"]), "nan"), "bracket_id": rw["bracket_id"], "row_left": rw["row_left"], "row_right": rw["row_right"], "max_row_left": rw["max_row_left"], "max_row_right": rw["max_row_right"], "source_label_file": rw["source_label_file"], "source_label_row": rw["source_label_row"], "status": "MAPPED"})
        lf = ["fit_cutoff","Lbyas","irrep","lattice_level_index","lattice_Ecm","lattice_sigma_or_weight","model_Ecm","residual_model_minus_lattice","normalized_residual","chi2_contribution","bracket_id","source_label_file","source_label_row","source_lattice_file","status"]
        with (d / "model_ecm_levels.csv").open("w", newline="") as f:
            w=csv.DictWriter(f,fieldnames=["fit_cutoff","Lbyas","irrep","lattice_level_index","lattice_Ecm","lattice_sigma_or_weight","model_Ecm","residual_model_minus_lattice","normalized_residual","chi2_contribution","bracket_id","row_left","row_right","max_row_left","max_row_right","source_label_file","source_label_row","status"]); w.writeheader(); w.writerows(ml)
        (d / "model_ecm_levels.json").write_text(json.dumps(ml, indent=2) + "\n")
        with (d / "lattice_vs_model_levels.csv").open("w", newline="") as f:
            w=csv.DictWriter(f,fieldnames=lf); w.writeheader()
            for x in ml:
                src = source_by_level.get((x["irrep"], int(x["lattice_level_index"])), "")
                w.writerow({**{k:x[k] for k in lf if k in x}, "source_lattice_file": src})
        with (d / "lattice_vs_model_levels.md").open("w") as f:
            f.write(f"# Lattice versus model, fit cutoff {cut:.3f}\n\n")
            f.write("| irrep | level | lattice Ecm | model Ecm | model-lattice | sigma | chi2 contribution | source lattice file | source true-window file |\n|---|---:|---:|---:|---:|---:|---:|---|---|\n")
            for x in ml:
                src = source_by_level.get((x["irrep"], int(x["lattice_level_index"])), "")
                f.write(f"| {x['irrep']} | {x['lattice_level_index']} | {x['lattice_Ecm']:.17g} | {x['model_Ecm']:.17g} | {x['residual_model_minus_lattice']:.17g} | {x['lattice_sigma_or_weight']:.17g} | {x['chi2_contribution']} | {src} | {x['source_label_file']} |\n")
        all_global.extend(roots); all_levels.extend(ml)

    fields=list(all_global[0]);
    with (OUT/"L20_PASSONLY_TRUE_ZEROS_ALL_WINDOWS_ALL_CUTOFFS.csv").open("w",newline="") as f:
        w=csv.DictWriter(f,fieldnames=fields); w.writeheader(); w.writerows(all_global)
    with (OUT/"L20_PASSONLY_MODEL_ECM_LEVELS_ALL_CUTOFFS.csv").open("w",newline="") as f:
        w=csv.DictWriter(f,fieldnames=list(all_levels[0])); w.writeheader(); w.writerows(all_levels)
    (OUT/"L20_PASSONLY_TRUE_ZEROS_ALL_WINDOWS_ALL_CUTOFFS.md").write_text("# All user-true fitted-model roots\n\nEvery fitted cutoff includes every explicit user=true window through Ecm=0.36.\n")
    return 0


if __name__ == "__main__": raise SystemExit(main())

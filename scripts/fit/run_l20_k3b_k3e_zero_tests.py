#!/usr/bin/env python3
"""Run fixed-zero and float-from-zero L20 cutoff tests."""
from __future__ import annotations

import csv, hashlib, os, shutil, subprocess, sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
import run_l20_full_cutoff_fits_v2 as suite

WINDOW_ROOT = ROOT / "output/v33p_gpu_classifier_all_irrep/accepted_windows_L20_user_validated_0p36"
CACHE_ROOT = ROOT / "output/v33p_projected_basis_cache/hot_windows_L20_user_validated_0p36"
OUT = ROOT / "output/v33p_projected_basis_cache/l20_k3b_k3e_zero_tests"
SECTORS = ("000_A1m", "100_A2", "110_A2", "111_A2", "200_A2")
CUTS = (0.315, 0.325, 0.335, 0.345)
START = (73735.840894011912, -972421.14060757787, 0.0, 0.0)
BASE_TEMPLATE = ROOT / "output/v33p_projected_basis_cache/fitter_projected_basis_hot_windows_all.in"


def set_key(text: str, key: str, value: str) -> str:
    import re
    pat = re.compile(rf"(?m)^\s*{re.escape(key)}\s*=.*$")
    line = f"{key} = {value}"
    return pat.sub(line, text, count=1) if pat.search(text) else text + "\n" + line + "\n"


def zero_sources(test_out: Path) -> dict[str, Path]:
    root = test_out / "normalized_user_validated_windows"
    root.mkdir(parents=True, exist_ok=True)
    fields = ["bracket_id", "E_left", "E_right", "zero_estimate", "user_label", "inside_Ecm_cutoff", "row_left", "row_right", "max_row_left", "max_row_right", "lattice_level_index", "source_label_file", "source_label_row"]
    out = {}
    for sector in SECTORS:
        src = WINDOW_ROOT / f"L20_{sector}/accepted_windows.csv"
        dst = root / f"accepted_truezeros_L20_{sector}.csv"
        with src.open(newline="") as fi, dst.open("w", newline="") as fo:
            reader = csv.DictReader(fi); writer = csv.DictWriter(fo, fieldnames=fields); writer.writeheader()
            for r in reader:
                writer.writerow({"bracket_id": r["bracket_id"], "E_left": r["E_left_bracket"], "E_right": r["E_right_bracket"], "zero_estimate": r["zero_estimate_initial"], "user_label": "true", "inside_Ecm_cutoff": "true", "row_left": r["row_left"], "row_right": r["row_right"], "max_row_left": r["max_row_left"], "max_row_right": r["max_row_right"], "lattice_level_index": r["lattice_level_index"], "source_label_file": r["source_label_file"], "source_label_row": r["source_label_row"]})
        out[sector] = dst
    return out


def prefit_audit() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    rows = []
    for cut in CUTS:
        used, _ = suite.lattice_records(cut)
        for sector in SECTORS:
            wpath = WINDOW_ROOT / f"L20_{sector}/accepted_windows.csv"
            cdir = CACHE_ROOT / f"L20_{sector}"
            windows = [r for r in csv.DictReader(wpath.open(newline="")) if float(r["zero_estimate_initial"]) <= cut + 1e-15]
            requested = set()
            for r in windows:
                requested.update(range(int(r["max_row_left"]), int(r["max_row_right"]) + 1))
            meta = __import__("json").loads((cdir / "projected_basis_hotwindows.json").read_text())
            saved = set(meta.get("row_indices", [])); missing = sorted(requested - saved)
            hash_ok = meta.get("accepted_windows_sha256") == hashlib.sha256(wpath.read_bytes()).hexdigest()
            files_ok = all((cdir / x).exists() for x in ("projected_basis_hotwindows.bin", "projected_basis_hotwindows.json", "projected_basis_hotwindows_validation.csv", "PROJECTED_BASIS_HOTWINDOW_CACHE_REPORT.md"))
            lattice_n = sum(r["irrep"] == sector for r in used)
            status = "PASS" if files_ok and hash_ok and not missing and len(windows) == lattice_n else "BLOCKED"
            rows.append({"cutoff": cut, "sector": sector, "lattice_levels": lattice_n, "user_true_windows": len(windows), "missing_rows": len(missing), "hash_match": hash_ok, "cache_files": files_ok, "status": status})
    fields = list(rows[0]);
    with (OUT / "PREFIT_AUDIT_ZERO_TESTS.csv").open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields); w.writeheader(); w.writerows(rows)
    with (OUT / "PREFIT_AUDIT_ZERO_TESTS.md").open("w") as f:
        f.write("# Zero-test prefit audit\n\nOnly L20 user-true windows and projected-basis hot-window caches were checked.\n\n")
        f.write("| cutoff | sector | lattice | user true | missing rows | hash | files | status |\n|---:|---|---:|---:|---:|---|---|---|\n")
        for r in rows: f.write(f"| {r['cutoff']} | {r['sector']} | {r['lattice_levels']} | {r['user_true_windows']} | {r['missing_rows']} | {r['hash_match']} | {r['cache_files']} | {r['status']} |\n")
    if any(r["status"] != "PASS" for r in rows): raise SystemExit("zero-test prefit audit failed")


def run_test(name: str, fixed_zero: bool) -> None:
    test_out = OUT / name
    test_out.mkdir(parents=True, exist_ok=True)
    paths = zero_sources(test_out)
    base = BASE_TEMPLATE.read_text()
    for sector, path in paths.items(): base = set_key(base, f"accepted_zeros_file_L20_{sector}", str(path))
    base = set_key(base, "projected_basis_cache_root", str(CACHE_ROOT))
    base_path = test_out / "fitter_input_base.in"; base_path.write_text(base)

    def make_config(cut, guess, outdir, zero_paths, tag):
        text = base
        values = {"Lbyas_values": "20", "irreps_L24": "", "list_of_mom": "000_A1m 100_A2 110_A2 111_A2 200_A2", "Ecm_cutoff": f"{cut:.17g}", "energy_cutoff": f"{cut:.17g}", "one_fcn_only": "0", "max_fcn": "400", "max_fcn_evals": "500", "output_dir": str(outdir), "output_tag": tag, "K3iso0_guess": f"{guess[0]:.17g}", "K3iso1_guess": f"{guess[1]:.17g}", "K3B_guess": "0", "K3E_guess": "0", "K3iso0_step": "1000", "K3iso1_step": "10000", "K3B_step": "1000", "K3E_step": "10000", "float_K3iso0": "1", "float_K3iso1": "1", "float_K3B": "0" if fixed_zero else "1", "float_K3E": "0" if fixed_zero else "1", "chi2_trace_path": str(outdir / "fcn_trace.csv"), "chi2_level_trace_path": str(outdir / "fcn_level_trace.csv"), "chi2_trace_test_name": ("fixed_K3B0_K3E0" if fixed_zero else "float_from_K3B0_K3E0")}
        for k, v in values.items(): text = set_key(text, k, v)
        for k, v in zero_paths.items(): text = set_key(text, k, str(v))
        p = outdir / "fitter_input.txt"; p.write_text(text); return p

    old = (suite.BASE, suite.OUT, suite.WINDOW_ROOT, suite.CACHE_ROOT, suite.CUTS, suite.START, suite.make_config)
    suite.BASE, suite.OUT, suite.WINDOW_ROOT, suite.CACHE_ROOT, suite.CUTS, suite.START, suite.make_config = base_path, test_out, WINDOW_ROOT, CACHE_ROOT, CUTS, START, make_config
    failed = None
    try:
        suite.main()
    except (SystemExit, RuntimeError) as exc:
        failed = str(exc)
    finally:
        suite.BASE, suite.OUT, suite.WINDOW_ROOT, suite.CACHE_ROOT, suite.CUTS, suite.START, suite.make_config = old
    env = os.environ.copy(); env.update({"L20_ROOT_OUT": str(test_out), "L20_WINDOW_ROOT": str(WINDOW_ROOT), "L20_CACHE_ROOT": str(CACHE_ROOT)})
    completed = [f"{cut:.3f}" for cut in CUTS if (test_out / f"Ecm_cut_{cut:.3f}".replace(".", "p") / "final_parameters.csv").exists()]
    if completed:
        env["L20_CUTS"] = ",".join(completed)
        subprocess.run([sys.executable, str(ROOT / "scripts/save_l20_passonly_root_artifacts.py")], cwd=ROOT, env=env, check=True)
    if failed:
        (test_out / "TEST_ABORTED.md").write_text(f"# {name}\n\nStatus: BLOCKED after the first failed cutoff.\n\nReason: `{failed}`\n\nNo later cutoff was fitted and no roots were fabricated.\n")


def main() -> int:
    prefit_audit()
    run_test("fixed_K3B0_K3E0", True)
    run_test("float_from_K3B0_K3E0", False)
    subprocess.run([sys.executable, str(ROOT / "scripts/assemble_l20_zero_test_reports.py")], cwd=ROOT, check=True)
    return 0


if __name__ == "__main__": raise SystemExit(main())

#!/usr/bin/env python3
"""Run only the four L20 cutoffs that pass the user-validated prefit gate."""
from __future__ import annotations

import csv
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(Path(__file__).parent))
import run_l20_full_cutoff_fits_v2 as suite

WINDOW_ROOT = ROOT / "output/v33p_gpu_classifier_all_irrep/accepted_windows_L20_user_validated_0p36"
CACHE_ROOT = ROOT / "output/v33p_projected_basis_cache/hot_windows_L20_user_validated_0p36"
OUT = ROOT / "output/v33p_projected_basis_cache/l20_full_cutoff_fits_user_validated_0p36_passonly"
SECTORS = ("000_A1m", "100_A2", "110_A2", "111_A2", "200_A2")


def write_fitter_zero_sources() -> dict[str, Path]:
    root = OUT / "normalized_user_validated_windows"
    root.mkdir(parents=True, exist_ok=True)
    paths = {}
    fields = ["bracket_id", "E_left", "E_right", "zero_estimate", "user_label",
              "inside_Ecm_cutoff", "row_left", "row_right", "max_row_left",
              "max_row_right", "lattice_level_index", "source_label_file", "source_label_row"]
    for sector in SECTORS:
        source = WINDOW_ROOT / f"L20_{sector}/accepted_windows.csv"
        target = root / f"accepted_truezeros_L20_{sector}.csv"
        with source.open(newline="") as src, target.open("w", newline="") as dst:
            reader = csv.DictReader(src)
            writer = csv.DictWriter(dst, fieldnames=fields)
            writer.writeheader()
            for row in reader:
                writer.writerow({
                    "bracket_id": row["bracket_id"],
                    "E_left": row["E_left_bracket"],
                    "E_right": row["E_right_bracket"],
                    "zero_estimate": row["zero_estimate_initial"],
                    "user_label": "true",
                    "inside_Ecm_cutoff": "true",
                    "row_left": row["row_left"],
                    "row_right": row["row_right"],
                    "max_row_left": row["max_row_left"],
                    "max_row_right": row["max_row_right"],
                    "lattice_level_index": row["lattice_level_index"],
                    "source_label_file": row["source_label_file"],
                    "source_label_row": row["source_label_row"],
                })
        paths[sector] = target
    return paths


def main() -> int:
    OUT.mkdir(parents=True, exist_ok=True)
    zero_paths = write_fitter_zero_sources()
    base_text = suite.BASE.read_text()
    for sector, path in zero_paths.items():
        base_text = suite.set_key(base_text, f"accepted_zeros_file_L20_{sector}", str(path))
    base_text = suite.set_key(base_text, "projected_basis_cache_root", str(CACHE_ROOT))
    base_path = OUT / "fitter_input_base_user_validated.in"
    base_path.write_text(base_text)
    suite.BASE = base_path
    suite.OUT = OUT
    suite.WINDOW_ROOT = WINDOW_ROOT
    suite.CACHE_ROOT = CACHE_ROOT
    suite.CUTS = (0.315, 0.325, 0.335, 0.345)
    return suite.main()


if __name__ == "__main__":
    raise SystemExit(main())

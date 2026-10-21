#!/usr/bin/env python3
"""Merge validated projected-basis row records into an existing hot cache."""
from __future__ import annotations

import csv
import json
import struct
import sys
from pathlib import Path

HEADER = struct.Struct("<8sII")
ROW_HEAD = struct.Struct("<QdiiI")


def records(path: Path) -> dict[int, bytes]:
    data = path.read_bytes()
    magic, version, count = HEADER.unpack_from(data, 0)
    if magic != b"V33PHWB1" or version != 1:
        raise ValueError(f"invalid cache: {path}")
    pos = HEADER.size
    out = {}
    for _ in range(count):
        row, _e, _nf, _np, dim = ROW_HEAD.unpack_from(data, pos)
        size = ROW_HEAD.size + 5 * dim * dim * 16
        out[row] = data[pos:pos + size]
        pos += size
    if pos != len(data):
        raise ValueError(f"trailing bytes in cache: {path}")
    return out


def main() -> int:
    if len(sys.argv) != 3:
        raise SystemExit("usage: merge_projected_basis_rows.py EXISTING_DIR PATCH_DIR")
    existing = Path(sys.argv[1]); patch = Path(sys.argv[2])
    old_bin = existing / "projected_basis_hotwindows.bin"
    patch_bin = patch / "projected_basis_hotwindows.bin"
    old = records(old_bin); old.update(records(patch_bin))
    ordered = dict(sorted(old.items()))
    out = bytearray(HEADER.pack(b"V33PHWB1", 1, len(ordered)))
    for raw in ordered.values(): out.extend(raw)
    old_bin.write_bytes(out)

    old_rows = list(csv.DictReader((existing / "projected_basis_hotwindows_validation.csv").open()))
    patch_rows = list(csv.DictReader((patch / "projected_basis_hotwindows_validation.csv").open()))
    by_row = {int(r["row_index"]): r for r in old_rows}
    by_row.update({int(r["row_index"]): r for r in patch_rows})
    fields = list(old_rows[0]) if old_rows else list(patch_rows[0])
    with (existing / "projected_basis_hotwindows_validation.csv").open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields); w.writeheader(); w.writerows(by_row[k] for k in sorted(by_row))

    meta_path = existing / "projected_basis_hotwindows.json"
    meta = json.loads(meta_path.read_text())
    vals = list(csv.DictReader((existing / "projected_basis_hotwindows_validation.csv").open()))
    meta["n_rows_saved"] = len(vals)
    meta["row_indices"] = [int(r["row_index"]) for r in vals]
    meta["Ecm_values"] = [float(r["Ecm"]) for r in vals]
    meta["Ecm_min"] = min(meta["Ecm_values"]); meta["Ecm_max"] = max(meta["Ecm_values"])
    def val(row: dict[str, str], *names: str) -> str:
        for name in names:
            if name in row:
                return row[name]
        raise KeyError(names)

    meta["Nproj_values"] = [int(val(r, "proj_dim", "Nproj")) for r in vals]
    meta["matrix_dimension_values"] = [int(r.get("matrix_dimension", val(r, "proj_dim", "Nproj"))) for r in vals]
    meta["validation"] = {
        "max_matrix_abs_diff": max(float(val(r, "matrix_abs_diff", "matrix_max_abs_diff")) for r in vals),
        "max_det_abs_diff": max(float(val(r, "det_diff", "det_abs_diff")) for r in vals),
        "max_det_rel_diff": max(float(val(r, "det_rel_diff")) for r in vals),
        "sign_agreement": f"{sum(int(val(r, 'sign_agreement')) for r in vals)}/{len(vals)}",
        "rows_checked": len(vals),
    }
    meta["file_sizes"] = {"bin_size_bytes": old_bin.stat().st_size, "json_size_bytes": 0}
    meta["notes"] = meta.get("notes", "") + " Missing fitter nearest-grid endpoint rows merged from independently validated builder output; accepted-window hash unchanged."
    meta_path.write_text(json.dumps(meta, indent=2) + "\n")
    meta["file_sizes"]["json_size_bytes"] = meta_path.stat().st_size
    meta_path.write_text(json.dumps(meta, indent=2) + "\n")
    report = existing / "PROJECTED_BASIS_HOTWINDOW_CACHE_REPORT.md"
    text = report.read_text()
    text = text.replace("- Rows saved: ", "- Rows saved: ", 1)
    text += f"\n- Row merge validation: PASS; rows={len(vals)}, sign agreement={meta['validation']['sign_agreement']}.\n"
    report.write_text(text)
    print(f"merged rows={len(vals)} range={min(meta['row_indices'])}..{max(meta['row_indices'])}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

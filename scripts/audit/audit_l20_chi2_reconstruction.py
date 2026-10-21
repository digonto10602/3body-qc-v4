#!/usr/bin/env python3
"""Independent chi-square reconstruction for the L20 zero-structure tests."""
from __future__ import annotations

import csv
import math
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "output/v33p_projected_basis_cache/l20_k3b_k3e_zero_tests"
TESTS = ("fixed_K3B0_K3E0", "float_from_K3B0_K3E0")
CUTS = (0.315, 0.325, 0.335, 0.345)


def read_csv(path: Path):
    return list(csv.DictReader(path.open(newline="")))


def parse_matrix(path: Path) -> np.ndarray | None:
    if not path.exists():
        return None
    rows = []
    for line in path.read_text().splitlines():
        if line.strip() and not line.lstrip().startswith("#"):
            rows.append([float(x) for x in line.split()])
    return np.asarray(rows, dtype=float) if rows else None


def fit_dir(test: str, cut: float) -> Path:
    return OUT / test / f"Ecm_cut_{cut:.3f}".replace(".", "p")


def main() -> int:
    rows_out = []
    md = ["# Independent L20 chi-square reconstruction", "", "The script recomputes residual contractions from the saved lattice/model tables.", ""]
    for test in TESTS:
        for cut in CUTS:
            d = fit_dir(test, cut)
            metrics = read_csv(d / "fit_metrics.csv") if (d / "fit_metrics.csv").exists() else []
            levels = read_csv(d / "lattice_vs_model_levels.csv") if (d / "lattice_vs_model_levels.csv").exists() else []
            if not metrics or not levels:
                status = "BLOCKED_OR_NOT_RUN"
                rows_out.append({"test_name": test, "fit_cutoff": cut, "n_levels": len(levels), "reported_chi2": math.nan, "reported_chi2_dof": math.nan, "dof_used": math.nan, "dof_nlevels": len(levels), "dof_nlevels_minus_nfloat": math.nan, "diagonal_chi2": math.nan, "diagonal_chi2_dof": math.nan, "full_cov_chi2": math.nan, "full_cov_chi2_dof": math.nan, "reported_over_diagonal": math.nan, "reported_over_full_cov": math.nan, "covariance_file": "", "status": status})
                continue
            m = metrics[0]
            residual = np.asarray([float(r["residual_model_minus_lattice"]) for r in levels], dtype=float)
            sigma = np.asarray([float(r.get("lattice_Ecm_err_or_sigma", r.get("lattice_sigma_or_weight", "nan"))) for r in levels], dtype=float)
            diagonal = float(np.sum((residual / sigma) ** 2))
            cov_candidates = sorted(d.glob("*_covariance_allL.dat"))
            cov_path = cov_candidates[0] if cov_candidates else None
            cov = parse_matrix(cov_path) if cov_path else None
            full = math.nan
            if cov is not None and cov.shape == (len(residual), len(residual)):
                full = float(residual @ np.linalg.pinv(cov, rcond=1e-12) @ residual)
            reported = float(m["chi2"])
            dof = float(m["dof"])
            nfloat = len(residual) - dof
            rows_out.append({"test_name": test, "fit_cutoff": cut, "n_levels": len(residual), "reported_chi2": reported, "reported_chi2_dof": float(m["chi2_dof"]), "dof_used": dof, "dof_nlevels": len(residual), "dof_nlevels_minus_nfloat": dof, "diagonal_chi2": diagonal, "diagonal_chi2_dof": diagonal / dof if dof > 0 else math.nan, "full_cov_chi2": full, "full_cov_chi2_dof": full / dof if dof > 0 and math.isfinite(full) else math.nan, "reported_over_diagonal": reported / diagonal if diagonal else math.nan, "reported_over_full_cov": reported / full if full else math.nan, "covariance_file": str(cov_path) if cov_path else "", "status": "PASS" if math.isfinite(full) and abs(full - reported) <= max(1e-10, 1e-7 * max(1.0, abs(reported))) else "RECONSTRUCTION_MISMATCH"})
            rr = rows_out[-1]
            md.append(f"## {test}, cutoff {cut:.3f}")
            md.append("")
            md.append(f"- levels={rr['n_levels']}; dof used={rr['dof_used']}; inferred floating parameters={nfloat:g}")
            md.append(f"- reported chi2={reported:.17g}; reconstructed full-cov chi2={full:.17g}; diagonal-only chi2={diagonal:.17g}")
            md.append(f"- reported/full-cov ratio={rr['reported_over_full_cov']:.17g}; reported/diagonal ratio={rr['reported_over_diagonal']:.17g}; status={rr['status']}")
            md.append(f"- covariance source: `{cov_path}`")
            md.append("")
    fields = list(rows_out[0])
    with (OUT / "INDEPENDENT_CHI2_RECONSTRUCTION.csv").open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields); w.writeheader(); w.writerows(rows_out)
    (OUT / "INDEPENDENT_CHI2_RECONSTRUCTION.md").write_text("\n".join(md) + "\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

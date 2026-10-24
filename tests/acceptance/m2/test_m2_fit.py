"""M2 steps 5-6 (prompts/003 sections 2.6-2.7): the chi^2 of the tracked QC levels (`qc3.fit.QCChi2`) with the D1 covariance,
D3 matching and D8 exclusion, its analytic gradient, and the two-stage fit (Levenberg-Marquardt, then iminuit Migrad + Hesse)
with a fresh-scan verification of the roots at the minimum.

Tolerances:
  * model energies vs independent Brent roots: 1e-9 (prompts/001 section 4; same backend);
  * chi^2 vs an independent recomputation from the same roots: 1e-4 (rigorous bound 1.3e-5 for root differences of 1e-11,
    m2_helpers.chi2_bound_from_root_error; margin 8);
  * analytic gradient vs central differences (relative step 1e-4): 1e-3 of |grad| (truncation ~1e-8, round-off ~1e-6);
  * tracked vs fresh-scan roots at the minimum: 1e-9, hence |chi2 difference| <= bound(1e-9) <= 1.4e-3 at every cutoff
    (lambda_min(C) >= 1.6e-9, m <= 14, chi2 <= 50; the test evaluates the bound from the actual covariance);
  * Migrad convergence: EDM < 10 x 0.002*tol*errordef = 2e-3 (iminuit documents that Migrad may exceed edm_max by 10x);
  * Gauss-Newton vs Hesse errors within a factor 3 (prototype: 15-60 % apart; GN drops the residual x curvature term);
  * legacy code-equivalence (amendment 2026-10-24, acceptance-v7): our coded_v2 zeros at the legacy 0.335 fit's c vs the
    legacy model levels, which were linear interpolations of det Q sign changes on the 1-row grid (Delta = 4.845e-6):
    2e-7 when no pole lies within 10 rows (interpolation error ~ Delta^2/(8 d_pole) <= 6e-8, measured <= 4.0e-8 on 11
    levels, margin 5x) and one grid row otherwise (both codes bracket the zero in the same row pair; measured 6.1e-7 for
    110_A2 n=1 with a pole 1.7 rows away).
Reported-only quantities (Digonto 2026-10-24, flag-only): wrong-sign zeros (D9) and the `unconstrained` flag are asserted
to be present and self-consistent, never to have a particular value.
"""
import json

import numpy as np
import pytest

from m2_helpers import (C_BASELINE, C_LEGACY_0335, C_V32W, GRID_DE, GRID_E0, IRREPS, chi2_bound_from_root_error,
                        read_legacy_0335_fit_levels, sector)

TOL_E = 1e-9
TOL_CHI2_SAME_ROOTS = 1e-4
TOL_GRAD_REL = 1e-3
EDM_MAX = 2e-3
LBYAS = 20


def _ev():
    from qc3 import fast
    return fast.evaluate_point


def _levels(cutoff):
    from qc3 import lattice
    return lattice.load_levels(lattice.default_data_dir(), LBYAS, list(IRREPS), cutoff)


def _window_roots(levels, c, form):
    """Independent roots: Brent (qc3.roots.refine_root, kind='zero') on the labelled accepted windows for coded_v2 at
    C_BASELINE; the pairing to levels is by (irrep, level index)."""
    import csv
    import pathlib
    from qc3 import roots
    root = pathlib.Path(__file__).resolve().parents[3]
    out = np.full(len(levels.labels), np.nan)
    for i, (irrep, n) in enumerate(levels.labels):
        with open(root / "data" / "windows" / f"L{LBYAS}" / f"L{LBYAS}_{irrep}" / "accepted_windows.csv") as fh:
            rows = [r for r in csv.DictReader(fh) if int(r["lattice_level_index"]) == n]
        assert len(rows) == 1
        lo, hi = float(rows[0]["E_left_bracket"]), float(rows[0]["E_right_bracket"])
        out[i] = roots.refine_root(_ev(), sector(irrep, k2_pwave_form=form), tuple(c), lo, hi)
    return out


def test_qcchi2_at_baseline_matches_independent_roots_gradient_and_chi2():
    """Cutoff 0.315 (6 levels in 4 sectors), coded_v2, c = C_BASELINE, tracked range limited to [E0, 0.32]."""
    from qc3 import fit, lattice
    lv = _levels(0.315)
    assert len(lv.labels) == 6
    irreps = sorted({l[0] for l in lv.labels})
    trackers = fit.make_trackers(irreps, C_BASELINE, Lbyas=LBYAS, k2_pwave_form="coded_v2", E_min=GRID_E0, E_max=0.32, evaluate=_ev())
    f = fit.QCChi2(lv, trackers)
    assert f.errordef == 1.0
    chi2 = f(np.asarray(C_BASELINE))
    E_model, J_model = f.model(np.asarray(C_BASELINE))
    E_ref = _window_roots(lv, C_BASELINE, "coded_v2")
    assert np.max(np.abs(E_model - E_ref)) <= TOL_E
    chi2_ref = lattice.chi2(lv.Ecm - E_ref, lv.cov)
    assert abs(chi2 - chi2_ref) <= TOL_CHI2_SAME_ROOTS, (chi2, chi2_ref)
    assert chi2_bound_from_root_error(chi2, lv.cov, 1e-11) <= TOL_CHI2_SAME_ROOTS
    assert J_model.shape == (6, 4) and np.all(J_model[:, 0] < 0)
    # analytic gradient vs central finite differences of the chi^2 itself
    c0 = np.asarray(C_BASELINE)
    g = f.grad(c0)
    g_fd = np.zeros(4)
    for i in range(4):
        d = np.zeros(4); d[i] = 1e-4 * abs(c0[i])
        g_fd[i] = (f(c0 + d) - f(c0 - d)) / (2 * d[i])
    assert np.linalg.norm(g - g_fd) <= TOL_GRAD_REL * np.linalg.norm(g_fd), (g, g_fd)
    assert f.n_calls >= 9 and f.n_grad_calls >= 1
    # D3 with n > m: a sector keeps only its lowest m roots; nothing is penalised at the baseline
    assert f.n_penalised_calls == 0


@pytest.mark.slow
def test_cutoff_scan_q0_consistent_0335_two_stage_fit_verifies_and_records(tmp_path):
    """Cutoff 0.335 (12 levels, the best-constrained of the four production cutoffs), all four parameters free, from
    C_BASELINE.  Prototype: LM 19 residual evaluations, Migrad EDM 1e-11, 12 tracked roots vs fresh scan <= 1.6e-12."""
    from qc3 import fit, lattice
    out = tmp_path / "fit_L20_q0_test.json"
    fits = fit.cutoff_scan(lattice.default_data_dir(), LBYAS, (0.335,), C_BASELINE, "q0_consistent", out_json=out, verify=True,
                           evaluate=_ev())
    assert len(fits) == 1
    r = fits[0]
    assert r.cutoff == 0.335 and r.k2_pwave_form == "q0_consistent" and len(r.labels) == 12 and r.ndof == 8
    assert r.valid and np.isfinite(r.edm) and r.edm < EDM_MAX
    assert r.chi2 <= r.chi2_start and r.chi2 <= r.chi2_lm + 1e-6 and r.chi2_lm - r.chi2 <= EDM_MAX
    assert abs(r.chi2_per_dof * r.ndof - r.chi2) <= 1e-12
    assert np.all(np.asarray(r.errors) > 0) and np.all(np.asarray(r.errors_gauss_newton) > 0)
    assert np.all(np.asarray(r.errors) / np.asarray(r.errors_gauss_newton) < 3.0)
    assert np.all(np.asarray(r.errors_gauss_newton) / np.asarray(r.errors) < 3.0)
    cov = np.asarray(r.covariance); corr = np.asarray(r.correlation)
    assert cov.shape == (4, 4) and np.allclose(cov, cov.T) and np.linalg.eigvalsh(cov).min() > 0
    assert np.allclose(np.diag(corr), 1.0) and np.all(np.abs(corr) <= 1.0 + 1e-12)
    assert np.allclose(np.sqrt(np.diag(cov)), r.errors)
    assert r.n_penalised_calls_final_stage == 0 and r.n_full_rescans >= 0
    # amendment 2026-10-24 (acceptance-v7): all four parameters free; flag-only fields present and self-consistent
    assert tuple(r.fixed) == () and r.n_free == 4
    assert isinstance(r.unconstrained, bool) and isinstance(r.unconstrained_reasons, list)
    assert r.unconstrained == (len(r.unconstrained_reasons) > 0)
    assert isinstance(r.wrong_sign_zeros, list) and r.n_wrong_sign_zeros == len(r.wrong_sign_zeros)
    for w in r.wrong_sign_zeros:                              # D9: recorded and reported, never dropped from matching
        assert w["irrep"] in IRREPS and np.isfinite(w["E"]) and w["slope"] < 0.0 and isinstance(w["matched"], bool)
    # verification at the minimum: a fresh label-free scan of every sector agrees with the tracked roots
    v = r.verification
    assert v["counts_equal"] is True
    assert v["max_root_diff"] <= TOL_E, v
    bound = chi2_bound_from_root_error(r.chi2, np.asarray(r.cov_lattice), TOL_E)
    assert abs(v["chi2_fresh"] - r.chi2) <= bound, (v, bound)
    assert np.max(np.abs(np.asarray(r.residual_over_sigma))) < 10.0
    # provenance in the written file
    d = json.loads(out.read_text())
    assert d["is_a_fit"] is True
    prov = d["provenance"]
    for key in ("git_commit", "core_built_at_commit", "numpy", "scipy", "iminuit", "date", "machine"):
        assert key in prov, key
    assert d["constants"]["k2_pwave_form"] == "q0_consistent" and d["constants"]["production"]["xi"] == 3.444
    assert d["constants"]["covariance"].startswith("D1") and d["constants"]["matching"].startswith("D3")
    assert d["seeds"]["stochastic_inputs"] is None
    assert len(d["fits"]) == 1 and d["fits"][0]["cutoff"] == 0.335 and d["fits"][0]["chi2"] == r.chi2
    assert d["fits"][0]["verification"]["max_root_diff"] == v["max_root_diff"]
    assert d["fits"][0]["fixed"] == [] and d["fits"][0]["unconstrained"] == r.unconstrained
    assert d["fits"][0]["wrong_sign_zeros"] == r.wrong_sign_zeros and "wrong_sign_zeros_fresh" in d["fits"][0]["verification"]


@pytest.mark.slow
def test_v32w_code_equivalence_chi2_at_c_v32w_is_recorded_not_reproduced():
    """Same form (coded_v2), same recipe (D1 covariance, D3 matching, 12 levels at 0.335) as the v32w fit, but the repo's
    .jack inputs: the chi^2 at c_v32w is computed by the production path and checked against an independent recomputation.
    It is NOT compared with v32w's 11.4897 (the prototype gives 55.04: the inputs differ; prompts/003 section 7)."""
    from qc3 import fit, lattice, roots
    lv = _levels(0.335)
    assert len(lv.labels) == 12
    trackers = fit.make_trackers(IRREPS, C_V32W, Lbyas=LBYAS, k2_pwave_form="coded_v2", E_min=GRID_E0, E_max=0.36, evaluate=_ev())
    f = fit.QCChi2(lv, trackers)
    chi2 = f(np.asarray(C_V32W))
    E_model, _ = f.model(np.asarray(C_V32W))
    assert np.all(np.isfinite(E_model)) and f.n_penalised_calls == 0
    # D3: 000_A1m and 111_A2 have more QC zeros than lattice levels at this c; the lowest ones are used
    for irrep in ("000_A1m", "111_A2"):
        n_lat = sum(1 for l in lv.labels if l[0] == irrep)
        assert len(trackers[irrep].energies()) > n_lat
        lowest = np.sort(trackers[irrep].energies())[:n_lat]
        assert np.allclose(E_model[[i for i, l in enumerate(lv.labels) if l[0] == irrep]], lowest)
    # independent Brent refinement around every model energy
    E_ref = np.array([roots.refine_root(_ev(), sector(ir, k2_pwave_form="coded_v2"), C_V32W, E - 2 * GRID_DE, E + 2 * GRID_DE)
                      for (ir, _), E in zip(lv.labels, E_model)])
    assert np.max(np.abs(E_model - E_ref)) <= TOL_E
    assert abs(chi2 - lattice.chi2(lv.Ecm - E_ref, lv.cov)) <= TOL_CHI2_SAME_ROOTS
    assert chi2 > 0.0


@pytest.mark.slow
def test_cutoff_scan_q0_consistent_0335_fixed_KB_KE_two_parameter_fit_records_flags():
    """Amendment 2026-10-24 (acceptance-v7; Digonto: report BOTH tables).  The two-parameter table fixes K_B = K_E = 0
    (start (K_iso0, K_iso1) = baseline).  At this start vector every sector has exactly as many q0_consistent zeros
    <= 0.335 as lattice levels (planner scratch 2026-10-24: 1/2/4/2/3), so the fit is well posed.  Same acceptance as the
    four-parameter fit plus: the fixed parameters stay exactly 0 with zero errors and zero covariance rows/columns, the
    free 2x2 block is SPD with unit-diagonal correlation, ndof = 12 - 2."""
    from qc3 import fit, lattice
    c_start = (C_BASELINE[0], C_BASELINE[1], 0.0, 0.0)
    fits = fit.cutoff_scan(lattice.default_data_dir(), LBYAS, (0.335,), c_start, "q0_consistent", out_json=None, verify=True,
                           fixed=("KB", "KE"), evaluate=_ev())
    assert len(fits) == 1
    r = fits[0]
    assert tuple(r.fixed) == ("KB", "KE") and r.n_free == 2 and r.ndof == 10 and len(r.labels) == 12
    assert tuple(r.c_start) == c_start
    c = np.asarray(r.c); err = np.asarray(r.errors); cov = np.asarray(r.covariance); corr = np.asarray(r.correlation)
    assert c[2] == 0.0 and c[3] == 0.0
    assert err.shape == (4,) and err[2] == 0.0 and err[3] == 0.0 and np.all(err[:2] > 0)
    assert cov.shape == (4, 4) and np.all(cov[2:, :] == 0.0) and np.all(cov[:, 2:] == 0.0)
    assert np.allclose(cov[:2, :2], cov[:2, :2].T) and np.linalg.eigvalsh(cov[:2, :2]).min() > 0
    assert np.allclose(np.sqrt(np.diag(cov)[:2]), err[:2])
    assert np.allclose(np.diag(corr)[:2], 1.0) and abs(corr[0, 1]) <= 1.0 + 1e-12 and np.all(corr[2:, :] == 0.0)
    assert np.all(np.asarray(r.errors_gauss_newton)[:2] > 0) and np.all(np.asarray(r.errors_gauss_newton)[2:] == 0.0)
    assert r.valid and np.isfinite(r.edm) and r.edm < EDM_MAX
    assert r.chi2 <= r.chi2_start and r.chi2 <= r.chi2_lm + 1e-6 and r.chi2_lm - r.chi2 <= EDM_MAX
    assert abs(r.chi2_per_dof * r.ndof - r.chi2) <= 1e-12
    assert r.n_penalised_calls_final_stage == 0
    assert isinstance(r.unconstrained, bool) and r.unconstrained == (len(r.unconstrained_reasons) > 0)
    assert isinstance(r.wrong_sign_zeros, list) and r.n_wrong_sign_zeros == len(r.wrong_sign_zeros)
    v = r.verification
    assert v["counts_equal"] is True and v["max_root_diff"] <= TOL_E, v
    assert abs(v["chi2_fresh"] - r.chi2) <= chi2_bound_from_root_error(r.chi2, np.asarray(r.cov_lattice), TOL_E), v


@pytest.mark.slow
def test_legacy_0335_fit_levels_code_equivalence_at_c_legacy():
    """Amendment 2026-10-24 (acceptance-v7; Digonto: v32w = code-equivalence check).  The v32w model levels are not in the
    repo; the only legacy model levels at a known parameter vector are the 12 levels of the legacy hot-window fit at
    cutoff 0.335 (results/legacy_L20_passonly_hotwindow, coded_v2).  Our label-free zeros at that vector must agree with
    them to the legacy interpolation accuracy (module docstring), the D3 pairing of the lowest zeros must reproduce the
    legacy level assignment, and the chi^2 there (D1 covariance) must be finite without penalties.  Independent of the
    legacy chi^2 and errors, which are invalid (OPEN_ISSUES #1)."""
    from qc3 import fit, lattice
    lv = _levels(0.335)
    legacy = read_legacy_0335_fit_levels()
    assert [(ir, n) for ir, n, _, _ in legacy] == list(lv.labels)
    # the legacy table's lattice values differ from the D1 means by up to 1.2e-4 (0.19 sigma; STATUS 2026-10-23: the legacy
    # reader converted the mean E_lab, D1 converts per sample).  Only the legacy MODEL levels are used below.
    assert np.max(np.abs(np.array([e for _, _, e, _ in legacy]) - lv.Ecm)) <= 2e-4
    trackers = fit.make_trackers(IRREPS, C_LEGACY_0335, Lbyas=LBYAS, k2_pwave_form="coded_v2", E_min=GRID_E0, E_max=0.36,
                                 evaluate=_ev())
    f = fit.QCChi2(lv, trackers)
    chi2 = f(np.asarray(C_LEGACY_0335))
    E_model, _ = f.model(np.asarray(C_LEGACY_0335))
    assert np.isfinite(chi2) and chi2 > 0.0 and f.n_penalised_calls == 0
    nearest = np.zeros(len(legacy)); tol = np.zeros(len(legacy))
    for i, (irrep, n, _, E_leg) in enumerate(legacy):
        zs = trackers[irrep].energies(); poles = np.asarray(trackers[irrep].scan.poles())
        nearest[i] = zs[np.argmin(np.abs(zs - E_leg))]
        pole_within_10_rows = len(poles) and np.min(np.abs(poles - E_leg)) <= 10 * GRID_DE
        tol[i] = GRID_DE if pole_within_10_rows else 2e-7
        assert abs(nearest[i] - E_leg) <= tol[i], (irrep, n, E_leg, nearest[i], tol[i])
    assert np.allclose(E_model, nearest, rtol=0, atol=1e-12), "D3 pairing must reproduce the legacy level assignment"
    assert len(trackers["110_A2"].energies()[trackers["110_A2"].energies() <= 0.335]) == 4   # includes the pole-paired zero

"""M1 step 3 — the production ingredients (through the oracle) agree with the planner's independent NumPy
reference written from arXiv:2111.12734 (qc3ref). Tolerance 1e-12 relative Frobenius: measured agreement is
1e-16 to 2.5e-14 (reports/research/007 section 3.2-3.3), margin >= 40."""
import numpy as np
import pytest

from m1_helpers import C_BASELINE, C_ZERO, PRODUCTION, anchor_points, as_c, basis_list, E_lab_of, L_BOX, oracle_point, rel_frob
from qc3ref import fgk, k3

TOL = 1e-12
MK, MPI = PRODUCTION["atmK"], PRODUCTION["atmpi"]
SC = PRODUCTION["scatter"]


def _ingredients_vs_reference(irrep, Ecm, pwave_form="coded_v2", check_k3=True):
    kw = {"k2_pwave_form": pwave_form}   # always explicit (amendment acceptance-v4, D7: the default is no longer coded_v2)
    r = oracle_point(irrep, Ecm, C_ZERO, **kw)
    E, P = E_lab_of(irrep, Ecm)
    bl = basis_list(r)
    F_ref = fgk.F_matrix(bl, E, P, L_BOX, MK, MPI, alpha=PRODUCTION["alpha"], nmax=PRODUCTION["max_shell_num"])
    G_ref = fgk.G_matrix(bl, E, P, L_BOX, MK, MPI)
    K_ref = fgk.K2inv_matrix(bl, E, P, L_BOX, MK, MPI, SC["a0_Kpi"], SC["a1_Kpi"], SC["a0_KK"],
                             PRODUCTION["eta_1"], PRODUCTION["eta_2"], pwave_form)
    out = {"F": rel_frob(r.F, F_ref), "G": rel_frob(r.G, G_ref), "K2inv": rel_frob(r.K2inv, K_ref)}
    if check_k3:
        kb = k3.K3_basis_matrices(bl, E, P, L_BOX, MK, MPI)
        for name in ("iso0", "iso1", "B", "E"):
            out["K3_" + name] = rel_frob(r.K3_basis[name], kb[name])
    for name, d in out.items():
        assert d <= TOL, (irrep, Ecm, pwave_form, name, d)
    return out


def test_F_G_K2inv_vs_reference_small_points():
    pts = anchor_points(small_only=True)
    assert len(pts) >= 3
    for p in pts:
        _ingredients_vs_reference(p["irrep"], p["Ecm"], check_k3=False)


def test_K2inv_pwave_both_forms_vs_reference():
    p = [q for q in anchor_points(small_only=True) if q["irrep"] == "100_A2"][0]
    _ingredients_vs_reference(p["irrep"], p["Ecm"], "coded_v2", check_k3=False)
    _ingredients_vs_reference(p["irrep"], p["Ecm"], "q0_consistent", check_k3=False)
    # the switch must act: the two forms differ in the p-wave rows
    r1 = oracle_point(p["irrep"], p["Ecm"], C_ZERO)
    r2 = oracle_point(p["irrep"], p["Ecm"], C_ZERO, k2_pwave_form="q0_consistent")
    assert r2.metadata["k2_pwave_form"] == "q0_consistent"
    assert rel_frob(r1.K2inv, r2.K2inv) > 1e-3
    pw = np.asarray(r1.basis["ell"]) == 1
    d = np.abs(as_c(r1.K2inv) - as_c(r2.K2inv))
    assert d[np.ix_(~pw, ~pw)].max() == 0.0 and d[np.ix_(pw, pw)].max() > 0.0


def test_K3_basis_vs_reference_small_points():
    pts = [p for p in anchor_points(small_only=True) if p["irrep"] in ("100_A2", "000_A1m")]
    assert len(pts) >= 2
    for p in pts[:2]:
        _ingredients_vs_reference(p["irrep"], p["Ecm"], check_k3=True)


def test_K3_basis_symmetric_and_linear():
    p = anchor_points(small_only=True)[0]
    r0 = oracle_point(p["irrep"], p["Ecm"], C_ZERO)
    rc = oracle_point(p["irrep"], p["Ecm"], C_BASELINE)
    for name, M in r0.K3_basis.items():
        assert rel_frob(M, as_c(M).T) <= TOL, name
    K3_sum = sum(ci * as_c(r0.K3_basis[n]) for ci, n in zip(C_BASELINE, ("iso0", "iso1", "B", "E")))
    assert rel_frob(rc.K3, K3_sum) <= TOL
    assert np.abs(as_c(r0.K3)).max() == 0.0


@pytest.mark.slow
def test_ingredients_vs_reference_large_points():
    for irrep, Ecm in (("110_A2", 0.35), ("111_A2", 0.35103)):
        _ingredients_vs_reference(irrep, Ecm, check_k3=True)

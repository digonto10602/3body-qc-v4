"""Unit tests for prompts/002 sub-steps 7.4 (K3 pair kernel), 7.5 (G pair kernel) and 7.6 (C++ tail).

Small systems only (N <= 68), one oracle call per point (cached in the module).  Tolerances: the pair kernels are
checked against the oracle's K3mat_2plus1 / G_2plus1_mat to 1e-12 relative Frobenius (the locked acceptance tolerance);
the measured deviations are far smaller (K3 <= 4e-17, G full 0.0, G symmetric <= the oracle's own asymmetry)."""
import numpy as np
import pytest

from qc3 import _core, fast, oracle
from qc3.fast import _params
from qc3.sector import Sector

C_BASE = (73735.840894011912, -972421.14060757787, 347174.05548116949, -1226756.7068845264)
POINTS = [("100_A2", 0.2915561), ("000_A1m", 0.26932), ("111_A2", 0.32036), ("110_A2", 0.30733), ("200_A2", 0.28664)]
NAMES = ("iso0", "iso1", "B", "E")
_CACHE = {}


def _sec(irrep):
    return Sector.production(irrep, k2_pwave_form="coded_v2")


def _oracle(irrep, Ecm):
    if (irrep, Ecm) not in _CACHE:
        _CACHE[(irrep, Ecm)] = oracle.evaluate_point(_sec(irrep), Ecm, C_BASE)
    return _CACHE[(irrep, Ecm)]


def _rel(a, b):
    a = np.asarray(a, dtype=complex); b = np.asarray(b, dtype=complex)
    nb = np.linalg.norm(b)
    return float(np.linalg.norm(a - b) / nb) if nb > 0 else float(np.linalg.norm(a - b))


@pytest.mark.parametrize("irrep,Ecm", POINTS)
def test_k3_pair_kernel_matches_oracle_and_is_symmetric(irrep, Ecm):
    o = _oracle(irrep, Ecm)
    f = fast.evaluate_point(_sec(irrep), Ecm, C_BASE)
    assert f.N == o.N
    for n in NAMES:
        assert _rel(f.K3_basis[n], o.K3_basis[n]) <= 1e-12, (irrep, n)
        assert np.array_equal(f.K3_basis[n], np.asarray(f.K3_basis[n]).T), (irrep, n)   # exact symmetry by construction
        assert np.all(np.asarray(f.K3_basis[n]).imag == 0.0)
    assert _rel(f.K3, o.K3) <= 1e-12
    # K3(c) = sum_i c_i B_i, same order of accumulation as the oracle
    assert _rel(f.K3, sum(ci * np.asarray(f.K3_basis[n]) for ci, n in zip(C_BASE, NAMES))) <= 1e-15


def test_k3_basis_matrices_are_independent_of_which_are_needed():
    """Only the matrices with c_i != 0 are built when K3 is not requested: K3(c) must be the same either way."""
    irrep, Ecm = "110_A2", 0.30733
    full = fast.evaluate_point(_sec(irrep), Ecm, C_BASE, want=("K3", "Q"))
    for c in [(0.0, 0.0, C_BASE[2], 0.0), (0.0, 0.0, 0.0, C_BASE[3]), (C_BASE[0], C_BASE[1], 0.0, 0.0)]:
        sub = fast.evaluate_point(_sec(irrep), Ecm, c, want=("Q",))
        ref = fast.evaluate_point(_sec(irrep), Ecm, c, want=("K3", "Q"))
        assert sub.K3_basis is None and ref.K3_basis is not None
        assert sub.sign == ref.sign and abs(sub.lnabsdet - ref.lnabsdet) <= 1e-12 * max(1.0, abs(ref.lnabsdet))
        assert np.array_equal(sub.Q, ref.Q)
    assert full.K3_basis is not None


@pytest.mark.parametrize("irrep,Ecm", POINTS)
def test_g_pair_kernel_matches_oracle(irrep, Ecm):
    """Default fill (every ordered pair): same operation order as G_2plus1_mat, so G agrees to round-off or better.
    Symmetric fill (plan 2.3): exactly symmetric, and it deviates from the oracle by at most the oracle's own asymmetry."""
    o = _oracle(irrep, Ecm)
    f = fast.evaluate_point(_sec(irrep), Ecm, C_BASE)
    assert _rel(f.G, o.G) <= 1e-14, (irrep, _rel(f.G, o.G))
    par = _params(_sec(irrep), Ecm, C_BASE, False)
    Gs = np.asarray(_core.fast_point(par, 1e-16, "orbit", ["G"], True)["G"])
    asym = _rel(np.asarray(o.G).T, o.G)
    assert np.array_equal(Gs, Gs.T)
    assert _rel(Gs, o.G) <= 1e-12
    assert _rel(Gs, o.G) <= asym + 1e-15, (irrep, _rel(Gs, o.G), asym)
    Gf = np.asarray(_core.fast_point(par, 1e-16, "orbit", ["G"], False)["G"])
    assert np.array_equal(Gf, np.asarray(f.G))
    # the (2,2) block (pion-pion) is exactly zero
    n1 = f.N1
    assert not np.any(np.asarray(f.G)[n1:, n1:])


@pytest.mark.parametrize("irrep,Ecm", POINTS)
def test_cpp_tail_matches_oracle_projection(irrep, Ecm):
    """sign and ln|det Q| equal the oracle's (1e-6, the locked tolerance); Q equals NumPy's V^dag (F3inv + K3) V / L^6
    built from the fast path's OWN F3inv, K3, V (same projector), i.e. LU-solve with N_proj right-hand sides == projecting
    the inverse."""
    o = _oracle(irrep, Ecm)
    f = fast.evaluate_point(_sec(irrep), Ecm, C_BASE)
    assert f.metadata["tail"] == "cpp"
    assert f.sign == o.sign and abs(f.lnabsdet - o.lnabsdet) <= 1e-6, (f.lnabsdet, o.lnabsdet)
    V = np.asarray(f.Vsel)
    Q_np = V.conj().T @ (np.asarray(f.F3inv) + np.asarray(f.K3)) @ V / f.sector.L ** 6
    assert _rel(f.Q, Q_np) <= 1e-8, (irrep, _rel(f.Q, Q_np))
    assert f.Q.shape == (f.Nproj, f.Nproj) and f.Q.dtype == np.complex128


def test_want_subsets_and_lazy_conversion():
    irrep, Ecm = "100_A2", 0.2915561
    sub = fast.evaluate_point(_sec(irrep), Ecm, C_BASE, want=("Q",))
    full = fast.evaluate_point(_sec(irrep), Ecm, C_BASE)
    for name in ("F", "G", "K2inv", "F3", "F3inv", "Vsel", "K3", "K3_basis"):
        assert getattr(sub, name) is None, name
    assert sub.Q is not None and sub.sign == full.sign and abs(sub.lnabsdet - full.lnabsdet) <= 1e-9
    # the full inverse is formed only when asked: its stage is not part of the timed projection path
    only_inv = fast.evaluate_point(_sec(irrep), Ecm, C_BASE, want=("F3inv",))
    assert only_inv.F3inv is not None and only_inv.Q is None and only_inv.F is None
    assert np.array_equal(only_inv.F3inv, full.F3inv)
    for k in ("basis", "F", "K2inv", "G", "F3", "projector", "F3inv", "K3", "projection", "lu", "total"):
        assert k in full.timings


def test_fast_module_has_no_numpy_tail():
    """prompts/002 7.6: qc3.fast must not call into qc3.qc (NumPy/BLAS) any more."""
    assert not hasattr(fast, "project_Q") and not hasattr(fast, "slogdet")
    import ast
    tree = ast.parse(open(fast.__file__).read())
    imported = set()
    for node in ast.walk(tree):
        if isinstance(node, ast.ImportFrom):
            imported.add((node.module or "", node.level))
        elif isinstance(node, ast.Import):
            imported.update((a.name, 0) for a in node.names)
    assert ("qc", 1) not in imported and all(name.split(".")[-1] != "qc" for name, _ in imported), imported


def test_complex_projector_branch_agrees_with_real_branch():
    """The eigen projector may return complex columns; the tail then solves with [Re V | Im V] right-hand sides.  Whichever
    branch runs, Q must be the Hermitian projection of the inverse: compare orbit (real V) and eigen at the same point."""
    irrep, Ecm = "111_A2", 0.32036
    fo = fast.evaluate_point(_sec(irrep), Ecm, C_BASE, options={"projector": "orbit"})
    fe = fast.evaluate_point(_sec(irrep), Ecm, C_BASE, options={"projector": "eigen"})
    assert fo.sign == fe.sign and abs(fo.lnabsdet - fe.lnabsdet) <= 1e-8
    eo = np.sort(np.linalg.eigvals(np.asarray(fo.Q)).real)
    ee = np.sort(np.linalg.eigvals(np.asarray(fe.Q)).real)
    assert np.allclose(eo, ee, rtol=1e-6, atol=1e-6 * np.abs(eo).max())
    # Q is Hermitian
    for f in (fo, fe):
        Q = np.asarray(f.Q)
        assert np.linalg.norm(Q - Q.conj().T) <= 1e-9 * np.linalg.norm(Q)

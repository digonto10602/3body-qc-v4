"""Unit tests for qc3.qc, qc3.roots, qc3.lattice (prompts/001 step 4). Synthetic data, fixed seeds."""
import types

import numpy as np
import pytest

from qc3.lattice import (chi2, elab_to_ecm, jackknife_mean_cov, jackknife_resample, load_levels, match_levels,
                         read_jack)
from qc3.point import MissingRootError
from qc3.qc import project_Q, slogdet
from qc3.roots import refine_root


def _write_jack(path, values, header_n=None):
    with open(path, "w") as f:
        f.write(f"{len(values) if header_n is None else header_n} 1 0 0 1\n")
        for v in values:
            f.write(f" 0 {v:.17g}\n")


# ---------------------------------------------------------------- qc
def test_project_Q_expression_and_real_inputs():
    rng = np.random.default_rng(0)
    F = rng.standard_normal((6, 6)); F = F + F.T
    K = rng.standard_normal((6, 6)); K = K + K.T
    V = np.linalg.qr(rng.standard_normal((6, 3)))[0]
    L = 68.88
    Q = project_Q(F, V, K, L)
    assert np.array_equal(Q, V.conj().T @ (F + K) @ V / L ** 6)
    assert Q.shape == (3, 3)


def test_slogdet_values_and_errors():
    s, l = slogdet(np.diag([2.0, -3.0, 0.5]))
    assert s == -1 and isinstance(s, int) and abs(l - np.log(3.0)) < 1e-14
    with pytest.raises(ValueError):
        slogdet(np.diag([1.0, 0.0]))                  # singular: never -inf silently
    with pytest.raises(ValueError):
        slogdet(np.array([[1.0, np.nan], [0.0, 1.0]]))
    with pytest.raises(ValueError):
        slogdet(np.zeros((2, 3)))
    with pytest.raises(ValueError):
        slogdet(np.array([[1j, 0], [0, 1.0]]))        # phase i: not +-1


# ---------------------------------------------------------------- roots
def _fake_eval(E0, sing_lndet=0.0):
    def evaluate(sector, E, c, want=("Q",)):
        d = E - E0
        return types.SimpleNamespace(sign=1 if d > 0 else -1, lnabsdet=float(np.log(abs(d)) if d != 0 else -np.inf) + sing_lndet)
    return evaluate


def test_refine_root_finds_simple_zero_and_raises_without_sign_change():
    ev = _fake_eval(0.3012345678901)
    r = refine_root(ev, None, (0, 0, 0, 0), 0.30, 0.31, xtol=1e-13)
    assert abs(r - 0.3012345678901) <= 1e-12
    with pytest.raises(MissingRootError):
        refine_root(ev, None, (0, 0, 0, 0), 0.31, 0.32)


# ---------------------------------------------------------------- lattice
def test_read_jack_header_mismatch_raises(tmp_path):
    p = tmp_path / "bad.jack"
    _write_jack(p, [0.1, 0.2, 0.3], header_n=4)
    with pytest.raises(ValueError):
        read_jack(p)


def test_sigma_is_sd_over_sqrt_N_for_per_config_file(tmp_path):
    rng = np.random.default_rng(11)
    N = 400
    x = 0.30 + 0.01 * rng.standard_normal(N)
    d = tmp_path / "L20"; d.mkdir()
    _write_jack(d / "20_000_A1m_n0.jack", x)
    ls = load_levels(tmp_path, 20, ["000_A1m"], 1.0)
    assert ls.Ncfg == N
    assert abs(np.sqrt(ls.cov[0, 0]) - x.std(ddof=1) / np.sqrt(N)) <= 1e-10 * x.std(ddof=1) / np.sqrt(N)
    assert abs(ls.Ecm[0] - x.mean()) <= 1e-13


def test_jackknife_samples_input_gives_sqrt_N_minus_1_over_N_sum(tmp_path):
    rng = np.random.default_rng(12)
    y = 0.3 + 0.001 * rng.standard_normal((50, 1))
    mean, cov = jackknife_mean_cov(y)
    N = 50
    assert abs(np.sqrt(cov[0, 0]) - np.sqrt((N - 1) / N * np.sum((y[:, 0] - y.mean()) ** 2))) <= 1e-15


def test_cutoff_selection_and_cov_shape(tmp_path):
    rng = np.random.default_rng(13)
    d = tmp_path / "lattice" / "L20"; d.mkdir(parents=True)
    for n, m0 in enumerate((0.30, 0.34, 0.31)):
        _write_jack(d / f"20_000_A1m_n{n}.jack", m0 + 0.002 * rng.standard_normal(30))
    ls = load_levels(tmp_path, 20, ["000_A1m"], 0.32)
    assert ls.labels == [("000_A1m", 0), ("000_A1m", 2)] and ls.cov.shape == (2, 2)
    assert np.all(ls.Ecm <= 0.32)
    ls0 = load_levels(tmp_path, 20, ["000_A1m"], 0.1)
    assert len(ls0.Ecm) == 0 and ls0.cov.shape == (0, 0)


def test_elab_to_ecm():
    assert elab_to_ecm(np.array([5.0]), 3.0)[0] == pytest.approx(4.0)


def test_match_levels_edge_cases():
    assert match_levels(np.array([]), np.array([0.3])) == []
    # extra QC levels are ignored (lowest m used)
    assert match_levels(np.array([0.30]), np.array([0.5, 0.31, 0.4])) == [(0, 1)]
    # ties/identical energies: stable ordering, deterministic
    assert match_levels(np.array([0.3, 0.3]), np.array([0.3, 0.3])) == [(0, 0), (1, 1)]
    with pytest.raises(MissingRootError):
        match_levels(np.array([0.3]), np.array([]))


def test_chi2_matches_solve_and_drops_null_direction():
    C = np.diag([1.0, 4.0, 0.0])
    r = np.array([1.0, 2.0, 5.0])
    assert chi2(r, C) == pytest.approx(1.0 + 1.0)

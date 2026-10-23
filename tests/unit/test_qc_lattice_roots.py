"""Unit tests for qc3.qc, qc3.roots, qc3.lattice (prompts/001 step 4). Synthetic data, fixed seeds."""
import pathlib
import types

import numpy as np
import pytest

from qc3.lattice import (chi2, elab_to_ecm, jackknife_mean_cov, jackknife_resample, load_levels, match_levels,
                         read_jack)
from qc3.point import MissingRootError
from qc3.qc import project_Q, slogdet
from qc3 import roots as R
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


def test_jackknife_mean_cov_vs_explicit_loops():
    """Independent check: explicit double loops over samples and levels (no vectorised formula shared with the code)."""
    rng = np.random.default_rng(12)
    N, m = 40, 3
    y = 0.3 + 0.001 * rng.standard_normal((N, m))
    mean, cov = jackknife_mean_cov(y)
    mean_ref = [sum(y[c][i] for c in range(N)) / N for i in range(m)]
    cov_ref = np.zeros((m, m))
    for i in range(m):
        for j in range(m):
            acc = 0.0
            for c in range(N):
                acc += (y[c][i] - mean_ref[i]) * (y[c][j] - mean_ref[j])
            cov_ref[i, j] = (N - 1) / N * acc
    assert np.allclose(mean, mean_ref, rtol=0, atol=1e-15)
    assert np.allclose(cov, cov_ref, rtol=1e-12, atol=0)
    # and the leave-one-out resample, loop version
    x = rng.standard_normal(N)
    yr = jackknife_resample(x)
    for c in range(N):
        assert abs(yr[c] - sum(x[k] for k in range(N) if k != c) / (N - 1)) < 1e-14


def test_load_levels_nonfinite_ecm_raises_naming_file(tmp_path):
    d = tmp_path / "L20"; d.mkdir()
    # P for 100_A2 at L20 is 2 pi/68.88 = 0.0912; E_lab below it -> E_lab^2 < P^2
    _write_jack(d / "20_001_A2_n0.jack", [0.05, 0.06, 0.05, 0.06, 0.05])
    with pytest.raises(ValueError, match="20_001_A2_n0.jack"):
        load_levels(tmp_path, 20, ["100_A2"], 1.0)


def test_chi2_rejects_asymmetric_non_psd_and_nonpositive():
    r = np.array([1.0, 2.0])
    with pytest.raises(ValueError):
        chi2(r, np.array([[1.0, 0.5], [0.0, 1.0]]))               # not symmetric
    with pytest.raises(ValueError):
        chi2(r, np.array([[1.0, 2.0], [2.0, 1.0]]))               # eigenvalue -1 < 0
    with pytest.raises(ValueError):
        chi2(r, np.zeros((2, 2)))                                 # lambda_max = 0
    with pytest.raises(ValueError):
        chi2(r, -np.eye(2))                                       # lambda_max < 0
    with pytest.raises(ValueError):
        chi2(r, np.eye(3))                                        # shape mismatch
    assert chi2(r, np.eye(2)) == pytest.approx(5.0)               # fine case unchanged


def test_load_levels_malformed_level_file_name_raises(tmp_path):
    d = tmp_path / "L20"; d.mkdir()
    _write_jack(d / "20_001_A2_n0.jack", [0.30, 0.31, 0.29])
    _write_jack(d / "20_001_A2_nX.jack", [0.30, 0.31, 0.29])      # matches the prefix, bad suffix
    with pytest.raises(ValueError, match="nX"):
        load_levels(tmp_path, 20, ["100_A2"], 1.0)
    (d / "20_001_A2_nX.jack").unlink()
    _write_jack(d / "notes_other.jack", [1.0, 2.0, 3.0])          # does not match the pattern at all: ignored
    _write_jack(d / "24_001_A2_n0.jack", [0.30, 0.31, 0.29])      # other L: ignored
    ls = load_levels(tmp_path, 20, ["100_A2"], 1.0)
    assert ls.labels == [("100_A2", 0)]


def test_default_data_dir_env_and_repo_root(tmp_path, monkeypatch):
    from qc3.lattice import default_data_dir
    monkeypatch.setenv("QC_DATA_ROOT", "/nonexistent/qc_data")
    assert default_data_dir() == pathlib.Path("/nonexistent/qc_data")
    with pytest.raises(FileNotFoundError):                         # loud failure, not an empty LevelSet
        load_levels(default_data_dir(), 20, ["000_A1m"], 1.0)
    monkeypatch.delenv("QC_DATA_ROOT")
    monkeypatch.chdir(tmp_path)                                    # independent of cwd
    d = default_data_dir()
    assert d == pathlib.Path(__file__).resolve().parents[2] / "data"
    assert (d / "lattice" / "L20" / "20_000_A1m_n0.jack").exists()


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


# ---------------------------------------------------------------- pole/zero discrimination (prompts/002 7.1)
def _synthetic(kind, E0, amp=1.0, bg=0.0, jump_at=None):
    """ln|det| = +-amp*ln|E-E0| + smooth background; sign flips at E0.  kind 'zero' (dip) or 'pole' (spike)."""
    def evaluate(sector, E, c, want=("Q",)):
        d = E - E0
        ln = (amp if kind == "zero" else -amp) * (np.log(abs(d)) if d != 0 else -np.inf) + bg * (E - E0)
        N = 50 if (jump_at is None or E < jump_at) else 52
        return types.SimpleNamespace(sign=1 if d > 0 else -1, lnabsdet=float(ln), N=N)
    return evaluate


def test_classify_synthetic_zero_pole_ambiguous():
    E0 = 0.3012345678
    z = R.classify_sign_change(_synthetic("zero", E0), None, (0,) * 4, E0 + 1e-12)
    p = R.classify_sign_change(_synthetic("pole", E0), None, (0,) * 4, E0 + 1e-12)
    assert z.kind == "zero" and z.discriminant == pytest.approx(np.log(1e-12 / 1e-6), abs=1e-3) and z.N == 50
    assert p.kind == "pole" and p.discriminant == pytest.approx(-np.log(1e-12 / 1e-6), abs=1e-3)
    assert z.n_evaluations == 3 and z.lnabsdet_minus == pytest.approx(z.lnabsdet_plus, abs=1e-3)
    # a sign change with a very weak singularity (amp=0.05): |Delta_h| < D -> ambiguous
    a = R.classify_sign_change(_synthetic("zero", E0, amp=0.05), None, (0,) * 4, E0 + 1e-12)
    assert a.kind == "ambiguous" and abs(a.discriminant) < 4.0
    # smooth linear background does not flip the verdict (O(h * slope))
    zb = R.classify_sign_change(_synthetic("zero", E0, bg=50.0), None, (0,) * 4, E0 + 1e-12)
    assert zb.kind == "zero"


def test_refine_root_kind_zero_pole_any():
    E0 = 0.3012345678
    z, p = _synthetic("zero", E0), _synthetic("pole", E0)
    c = (0,) * 4
    assert abs(R.refine_root(z, None, c, 0.30, 0.31) - E0) < 1e-10
    with pytest.raises(R.PoleError):
        R.refine_root(p, None, c, 0.30, 0.31)
    with pytest.raises(MissingRootError):
        R.refine_root(p, None, c, 0.30, 0.31, kind="zero")
    assert abs(R.refine_root(p, None, c, 0.30, 0.31, kind="any") - E0) < 1e-10      # old behaviour
    info = R.refine_root_info(p, None, c, 0.30, 0.31, kind="any")
    assert info.kind == "pole" and info.discriminant > 4 and info.bracket == (0.30, 0.31)
    with pytest.raises(R.AmbiguousRootError):
        R.refine_root(_synthetic("zero", E0, amp=0.05), None, c, 0.30, 0.31)
    with pytest.raises(ValueError):
        R.refine_root(z, None, c, 0.30, 0.31, kind="bogus")
    with pytest.raises(MissingRootError):
        R.refine_root(z, None, c, 0.31, 0.32)


def test_classify_dimension_jump_raises():
    E0 = 0.3012345678
    ev = _synthetic("zero", E0, jump_at=E0 + 5e-7)       # N changes between E*-h.. E*+h
    with pytest.raises(R.DimensionJumpError):
        R.classify_sign_change(ev, None, (0,) * 4, E0 + 1e-12)


def test_find_sign_changes_and_zeros_synthetic():
    E1, E2 = 0.3010, 0.3020           # a zero and a pole inside one window
    def evaluate(sector, E, c, want=("Q",)):
        z, p = E - E1, E - E2
        with np.errstate(divide="ignore"):
            ln = np.log(abs(z)) - np.log(abs(p))
        return types.SimpleNamespace(sign=(1 if z > 0 else -1) * (1 if p > 0 else -1), lnabsdet=float(ln), N=50)
    ch = R.find_sign_changes(evaluate, None, (0,) * 4, 0.2995, 0.3035, n_scan=41)
    assert [c.kind for c in ch] == ["zero", "pole"] and [c.E for c in ch] == sorted(c.E for c in ch)
    assert abs(ch[0].E - E1) < 1e-9 and abs(ch[1].E - E2) < 1e-9
    zs = R.find_zeros(evaluate, None, (0,) * 4, 0.2995, 0.3035, n_scan=41)
    assert len(zs) == 1 and zs[0].kind == "zero"
    # equal signs at both ends (zero+pole inside) -> plain refine_root finds nothing, as in step 9
    with pytest.raises(MissingRootError):
        R.refine_root(evaluate, None, (0,) * 4, 0.2995, 0.3035, kind="any")


def test_xtol_guard_rejects_ambiguous_discriminant_settings():
    """audit 001 N6: xtol > h*exp(-D) makes Delta_h of a simple zero ambiguous -> ValueError naming the reason (loud)."""
    E0 = 0.3012345678
    z = _synthetic("zero", E0)
    c = (0,) * 4
    h, D = 1e-6, 4.0
    xmax = h * np.exp(-D)
    # at the limit it still works (and classifies a zero)
    assert abs(R.refine_root(z, None, c, 0.30, 0.31, xtol=xmax) - E0) <= xmax
    for xtol in (2 * xmax, 1e-6, 1e-3):
        with pytest.raises(ValueError, match="discriminant"):
            R.refine_root(z, None, c, 0.30, 0.31, xtol=xtol)
        with pytest.raises(ValueError, match="discriminant"):
            R.refine_root_info(z, None, c, 0.30, 0.31, xtol=xtol, kind="any")
        with pytest.raises(ValueError, match="discriminant"):
            R.find_zeros(z, None, c, 0.30, 0.31, n_scan=11, xtol=xtol)
        with pytest.raises(ValueError, match="discriminant"):
            R.find_sign_changes(z, None, c, 0.30, 0.31, n_scan=11, xtol=xtol)
    # the limit follows h and D
    with pytest.raises(ValueError, match="discriminant"):
        R.refine_root(z, None, c, 0.30, 0.31, xtol=1e-9, h=1e-8)
    assert abs(R.refine_root(z, None, c, 0.30, 0.31, xtol=1e-9, h=1e-6) - E0) < 1e-9
    # kind="any" does not classify: no limit
    assert abs(R.refine_root(z, None, c, 0.30, 0.31, xtol=1e-6, kind="any") - E0) < 1e-6
    for bad in (0.0, -1e-11, float("nan")):
        with pytest.raises(ValueError):
            R.refine_root(z, None, c, 0.30, 0.31, xtol=bad)

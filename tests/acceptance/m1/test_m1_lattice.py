"""M1 step 4 — lattice data per decision D1 (per-configuration values: resample, convert, jackknife), D3
matching, chi^2 pseudo-inverse. Tolerances 1e-12 relative (closed-form statistics), counts exact."""
import os
import pathlib

import numpy as np
import pytest

from m1_helpers import ROOT, LBYAS, NNP, PRODUCTION

XI = PRODUCTION["xi"]
DATA = ROOT / "data"


def _write_jack(path, values):
    with open(path, "w") as f:
        f.write(f"{len(values)} 1 0 0 1\n")
        for v in values:
            f.write(f" 0 {v:.17g}\n")


def _P_of(irrep):
    return 2 * np.pi * np.sqrt(sum(d * d for d in NNP[irrep])) / (XI * LBYAS)


def test_read_jack_and_resample_synthetic(tmp_path):
    from qc3.lattice import jackknife_resample, read_jack
    x = np.array([0.30, 0.31, 0.29, 0.305, 0.295])
    p = tmp_path / "20_000_A1m_n0.jack"
    _write_jack(p, x)
    header, vals = read_jack(p)
    assert header["N"] == 5 and vals.dtype == np.float64 and np.array_equal(vals, x)
    y = jackknife_resample(vals)
    y_ref = (x.sum() - x) / (len(x) - 1)
    assert np.allclose(y, y_ref, rtol=0, atol=1e-15)


def test_jackknife_mean_cov_equals_sample_cov_over_N():
    from qc3.lattice import jackknife_mean_cov, jackknife_resample
    rng = np.random.default_rng(603)
    N, m = 603, 3
    A = rng.standard_normal((m, m))
    x = rng.standard_normal((N, m)) @ A + np.array([0.27, 0.30, 0.33])
    y = np.column_stack([jackknife_resample(x[:, j]) for j in range(m)])
    mean, cov = jackknife_mean_cov(y)
    assert np.allclose(mean, x.mean(axis=0), rtol=1e-12, atol=0)
    cov_ref = np.cov(x, rowvar=False, ddof=1) / N          # sigma = s.d./sqrt(N) for a linear map (D1)
    assert np.allclose(cov, cov_ref, rtol=1e-10, atol=0)
    assert np.allclose(np.sqrt(np.diag(cov)), x.std(axis=0, ddof=1) / np.sqrt(N), rtol=1e-10)


def test_load_levels_conversion_order_matches_legacy_reader(tmp_path):
    """Legacy (trusted v32w) order: resample E_lab, convert each resample to E*, then jackknife statistics."""
    from qc3.lattice import load_levels
    rng = np.random.default_rng(1)
    N = 50
    x0 = 0.31 + 0.01 * rng.standard_normal(N)
    x1 = 0.32 + 0.01 * rng.standard_normal(N) + 0.5 * (x0 - 0.31)
    d = tmp_path / "L20"; d.mkdir()
    _write_jack(d / "20_001_A2_n0.jack", x0)
    _write_jack(d / "20_001_A2_n1.jack", x1)
    P = _P_of("100_A2")
    ls = load_levels(tmp_path, LBYAS, ["100_A2"], Ecm_cutoff=1.0, xi=XI)
    assert ls.Ncfg == N and len(ls.Ecm) == 2 and ls.labels == [("100_A2", 0), ("100_A2", 1)]
    Y = np.column_stack([(x.sum() - x) / (N - 1) for x in (x0, x1)])
    Ecm_c = np.sqrt(Y ** 2 - P ** 2)
    mean_ref = Ecm_c.mean(axis=0)
    cov_ref = (N - 1) / N * (Ecm_c - mean_ref).T @ (Ecm_c - mean_ref)
    assert np.allclose(ls.Ecm, mean_ref, rtol=1e-12, atol=0)
    assert np.allclose(ls.cov, cov_ref, rtol=1e-10, atol=0)
    assert list(ls.nP[0]) == [0, 0, 1] and ls.Lbyas == LBYAS


def test_real_file_000_A1m_n0_vs_independent_numpy():
    from qc3.lattice import load_levels
    path = DATA / "lattice" / "L20" / "20_000_A1m_n0.jack"
    raw = np.loadtxt(path, skiprows=1)[:, 1]
    assert len(raw) == 603
    mean_ref = raw.mean()
    sigma_ref = raw.std(ddof=1) / np.sqrt(len(raw))
    ls = load_levels(DATA, LBYAS, ["000_A1m"], Ecm_cutoff=0.30, xi=XI)
    assert len(ls.Ecm) == 1 and ls.Ncfg == 603
    assert abs(ls.Ecm[0] - mean_ref) <= 1e-12 * mean_ref
    assert abs(np.sqrt(ls.cov[0, 0]) - sigma_ref) <= 1e-10 * sigma_ref
    assert 4e-4 <= sigma_ref <= 7e-4            # the (N-1)^2 bug would give ~0.31


def test_level_counts_per_cutoff_L20():
    """Counts recomputed from the files; must also equal the table of reports/research/004 (6/9/12/14/17)."""
    from qc3.lattice import load_levels
    irreps = ["000_A1m", "100_A2", "110_A2", "111_A2", "200_A2"]
    expected = {0.315: 6, 0.325: 9, 0.335: 12, 0.345: 14, 0.355: 17}
    # independent recount
    files = sorted((DATA / "lattice" / "L20").glob("20_*.jack"))
    assert len(files) == 19
    means = []
    for f in files:
        ddd = f.name.split("_")[1]
        irrep = {"000": "000_A1m", "001": "100_A2", "002": "200_A2", "011": "110_A2", "111": "111_A2"}[ddd]
        raw = np.loadtxt(f, skiprows=1)[:, 1]
        y = (raw.sum() - raw) / (len(raw) - 1)
        means.append(np.mean(np.sqrt(y ** 2 - _P_of(irrep) ** 2)))
    means = np.array(means)
    for cut, n in expected.items():
        assert int(np.sum(means <= cut)) == n, cut
        ls = load_levels(DATA, LBYAS, irreps, Ecm_cutoff=cut, xi=XI)
        assert len(ls.Ecm) == n == ls.cov.shape[0] == ls.cov.shape[1], cut
        assert np.all(ls.Ecm <= cut)


def test_covariance_symmetric_psd_and_correlated():
    from qc3.lattice import load_levels
    ls = load_levels(DATA, LBYAS, ["000_A1m", "100_A2", "110_A2", "111_A2", "200_A2"], Ecm_cutoff=0.335, xi=XI)
    C = ls.cov
    assert np.allclose(C, C.T, rtol=0, atol=1e-18)
    w = np.linalg.eigvalsh(C)
    assert w.min() >= -1e-12 * w.max()
    corr = C / np.sqrt(np.outer(np.diag(C), np.diag(C)))
    off = corr[~np.eye(len(C), dtype=bool)]
    assert np.abs(off).max() < 1.0 and np.abs(off).max() > 0.0     # same ensemble: levels are correlated


def test_match_levels_D3():
    from qc3.lattice import MissingRootError, match_levels
    lat = np.array([0.30, 0.28, 0.33])
    qc = np.array([0.31, 0.275, 0.34, 0.36, 0.29])
    pairs = match_levels(lat, qc)
    # lowest 3 QC levels (0.275, 0.29, 0.31) paired with sorted lattice (0.28, 0.30, 0.33)
    assert sorted(pairs) == sorted([(1, 1), (0, 4), (2, 0)])
    assert match_levels(np.array([0.3]), np.array([0.3])) == [(0, 0)]
    with pytest.raises(MissingRootError):
        match_levels(np.array([0.3, 0.31]), np.array([0.3]))


def test_chi2_pseudoinverse():
    from qc3.lattice import chi2
    rng = np.random.default_rng(5)
    A = rng.standard_normal((4, 4)); C = A @ A.T + np.eye(4)
    r = rng.standard_normal(4)
    assert abs(chi2(r, C) - r @ np.linalg.solve(C, r)) <= 1e-10 * abs(r @ np.linalg.solve(C, r))
    # rank-deficient: one exactly zero direction is dropped (lambda <= 1e-12 lambda_max)
    w, U = np.linalg.eigh(C)
    w[0] = 0.0
    C0 = (U * w) @ U.T
    ref = sum((U[:, i] @ r) ** 2 / w[i] for i in range(1, 4))
    assert abs(chi2(r, C0) - ref) <= 1e-10 * abs(ref)

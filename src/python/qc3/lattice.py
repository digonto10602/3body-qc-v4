"""Lattice .jack reader and level matching (prompts/001 step 4, 5).

Decision D1: a .jack file holds ONE VALUE PER GAUGE CONFIGURATION (N=603 for L20), not jackknife samples.  Recipe
(identical to the trusted v32w reader src/fitter/lattice_data_covariance_cpp.hpp):
    y_c = (sum_{c'} x_{c'} - x_c)/(N-1)                      (leave-one-out resample of E_lab)
    E*_c = sqrt(y_c^2 - P^2),  P = 2 pi |d| / (xi L/a_s)     (each resample converted to the cm frame)
    mean = (1/N) sum_c E*_c,  C_ij = (N-1)/N sum_c (E*_i,c - mean_i)(E*_j,c - mean_j)
Decision D3: levels are matched by sorting both lists and pairing the lowest m QC levels with the m lattice levels."""
from __future__ import annotations

import pathlib

import numpy as np

from .constants import NNP
from .point import LevelSet, MissingRootError  # noqa: F401

# file name code ddd  <->  irrep label
_DDD_TO_IRREP = {"000": "000_A1m", "001": "100_A2", "002": "200_A2", "011": "110_A2", "111": "111_A2"}
_IRREP_TO_DDD = {v: k for k, v in _DDD_TO_IRREP.items()}


def default_data_dir():
    raise NotImplementedError("qc3.lattice.default_data_dir: prompts/001 step 5")


def read_jack(path):
    """Read a .jack file: header line `N 1 0 0 1`, then N lines `0 value`.  Returns (header dict, values (N,))."""
    path = pathlib.Path(path)
    with open(path) as f:
        head = f.readline().split()
        if len(head) < 1:
            raise ValueError(f"{path}: empty header")
        header = {"N": int(head[0]), "fields": tuple(int(h) for h in head[1:])}
        vals = []
        for lineno, line in enumerate(f, start=2):
            tok = line.split()
            if not tok:
                continue
            if len(tok) < 2:
                raise ValueError(f"{path}:{lineno}: expected `0 value`, got {line!r}")
            vals.append(float(tok[1]))
    vals = np.asarray(vals, dtype=np.float64)
    if len(vals) != header["N"]:
        raise ValueError(f"{path}: header says N={header['N']} but {len(vals)} values found")
    return header, vals


def jackknife_resample(x):
    """Leave-one-out means: y_c = (sum(x) - x_c)/(N-1)."""
    x = np.asarray(x, dtype=np.float64)
    if x.ndim != 1 or len(x) < 2:
        raise ValueError("jackknife_resample needs a 1-d array with at least 2 entries")
    return (x.sum() - x) / (len(x) - 1)


def jackknife_mean_cov(samples):
    """samples (N, m) of resampled values -> mean (m,), C (m,m) = (N-1)/N sum_c (y_c - ybar)(y_c - ybar)^T."""
    y = np.asarray(samples, dtype=np.float64)
    if y.ndim == 1:
        y = y[:, None]
    N = y.shape[0]
    if N < 2:
        raise ValueError("jackknife_mean_cov needs at least 2 samples")
    mean = y.mean(axis=0)
    d = y - mean
    cov = (N - 1) / N * (d.T @ d)
    return mean, 0.5 * (cov + cov.T)


def elab_to_ecm(E_lab, P):
    """E* = sqrt(E_lab^2 - P^2) (elementwise)."""
    E_lab = np.asarray(E_lab, dtype=np.float64)
    return np.sqrt(E_lab ** 2 - P ** 2)


def _lattice_dir(data_dir, Lbyas):
    data_dir = pathlib.Path(data_dir)
    for cand in (data_dir / "lattice" / f"L{Lbyas}", data_dir / f"L{Lbyas}"):
        if cand.is_dir():
            return cand
    raise FileNotFoundError(f"no lattice directory L{Lbyas} under {data_dir} (tried lattice/L{Lbyas}, L{Lbyas})")


def load_levels(data_dir, Lbyas, irreps, Ecm_cutoff, xi=3.444):
    """All levels of the given irreps with mean E* <= Ecm_cutoff at one L, in (irrep order, n index) order,
    with the jackknife covariance over all kept levels (one L; same ensemble so levels are correlated)."""
    d = _lattice_dir(data_dir, Lbyas)
    labels, nps, files, cols, means = [], [], [], [], []
    Ncfg = None
    for irrep in irreps:
        if irrep not in _IRREP_TO_DDD:
            raise ValueError(f"unknown irrep {irrep!r}")
        ddd = _IRREP_TO_DDD[irrep]
        nP = tuple(NNP[irrep])
        P = 2.0 * np.pi * np.sqrt(sum(v * v for v in nP)) / (xi * Lbyas)
        found = []
        for f in d.glob(f"{Lbyas}_{ddd}_*_n*.jack"):
            stem = f.stem  # e.g. 20_001_A2_n3
            try:
                found.append((int(stem.rsplit("_n", 1)[1]), f))
            except ValueError:
                continue
        for n, f in sorted(found):
            header, x = read_jack(f)
            if Ncfg is None:
                Ncfg = header["N"]
            elif header["N"] != Ncfg:
                raise ValueError(f"{f}: N={header['N']} differs from {Ncfg}")
            ecm = elab_to_ecm(jackknife_resample(x), P)
            mean = float(ecm.mean())
            if mean <= Ecm_cutoff:
                labels.append((irrep, n)); nps.append(nP); files.append(str(f)); cols.append(ecm); means.append(mean)
    if cols:
        _, cov = jackknife_mean_cov(np.column_stack(cols))
        Ecm = np.asarray(means, dtype=np.float64)
        # mean over the resamples (same numbers as `means`, exactly)
    else:
        Ecm, cov = np.zeros(0), np.zeros((0, 0))
    return LevelSet(labels=labels, irreps=list(irreps), nP=nps, files=files, Ncfg=int(Ncfg or 0), Ecm=Ecm,
                    cov=cov, Ecm_cutoff=float(Ecm_cutoff), Lbyas=int(Lbyas))


def match_levels(E_lattice, E_qc):
    """D3: sort both lists ascending; pair the lowest m QC levels with the m lattice levels (n >= m allowed).
    Returns [(i_lattice, j_qc)] indices into the ORIGINAL arrays, ordered by lattice energy.
    Raises MissingRootError if there are fewer QC levels than lattice levels."""
    E_lattice = np.asarray(E_lattice, dtype=np.float64)
    E_qc = np.asarray(E_qc, dtype=np.float64)
    m, n = len(E_lattice), len(E_qc)
    if n < m:
        raise MissingRootError(f"only {n} QC levels for {m} lattice levels")
    il = np.argsort(E_lattice, kind="stable")
    iq = np.argsort(E_qc, kind="stable")[:m]
    return [(int(a), int(b)) for a, b in zip(il, iq)]


def chi2(residual, cov, eig_rel_cut=1e-12):
    """chi^2 = r^T C^+ r with the eigenvalue pseudo-inverse dropping eigenvalues <= eig_rel_cut * lambda_max."""
    r = np.asarray(residual, dtype=np.float64)
    C = np.asarray(cov, dtype=np.float64)
    w, U = np.linalg.eigh(0.5 * (C + C.T))
    keep = w > eig_rel_cut * w.max()
    proj = U[:, keep].T @ r
    return float(np.sum(proj ** 2 / w[keep]))

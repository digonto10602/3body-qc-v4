"""Lattice .jack reader and level matching (prompts/001 step 4, 5)."""
from .point import LevelSet, MissingRootError  # noqa: F401


def default_data_dir():
    raise NotImplementedError("qc3.lattice.default_data_dir: prompts/001 step 5")


def read_jack(path):
    raise NotImplementedError("qc3.lattice.read_jack: prompts/001 step 4")


def jackknife_resample(x):
    raise NotImplementedError("qc3.lattice.jackknife_resample: prompts/001 step 4")


def jackknife_mean_cov(samples):
    raise NotImplementedError("qc3.lattice.jackknife_mean_cov: prompts/001 step 4")


def elab_to_ecm(E_lab, P):
    raise NotImplementedError("qc3.lattice.elab_to_ecm: prompts/001 step 4")


def load_levels(data_dir, Lbyas, irreps, Ecm_cutoff, xi=3.444):
    raise NotImplementedError("qc3.lattice.load_levels: prompts/001 step 4")


def match_levels(E_lattice, E_qc):
    raise NotImplementedError("qc3.lattice.match_levels: prompts/001 step 4")


def chi2(residual, cov, eig_rel_cut=1e-12):
    raise NotImplementedError("qc3.lattice.chi2: prompts/001 step 4")

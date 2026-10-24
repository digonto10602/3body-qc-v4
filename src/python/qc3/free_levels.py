"""Non-interacting K K pi levels per sector (prompts/003 section 2.1; D8 input; docs/METHODS.md section 14).

With lattice momentum q = 2 pi / (xi L/a_s) and total momentum d = sector.nnP (any irrep: a conservative superset),

    E_lab = sqrt(mK^2 + q^2 n1^2) + sqrt(mK^2 + q^2 n2^2) + sqrt(mpi^2 + q^2 n3^2),   n1 + n2 + n3 = d,
    E*    = sqrt(E_lab^2 - q^2 d^2)                                               (energies in a_t^-1).

The two kaons are identical: an unordered pair, stored as n1 <= n2 (lexicographic on the integer triples).  Levels whose
E* agree within MERGE_TOL = 1e-11 are merged (consecutive sorted values closer than MERGE_TOL form one level); the
multiplicity is the number of distinct (n1, n2, n3) assignments.  The reported E is the mean of the merged energies.
These are the energies at which the truncated QC has "free-level solutions" (arXiv:1901.07095 Sec. 4.4.3), which D8
excludes from the level matching.

Enumeration bound: a particle cannot carry more momentum than the total energy allows.  With
E_lab,max = sqrt(E_max^2 + q^2 d^2), every kaon has omega_K <= E_lab,max - mK - mpi, hence
|n| <= sqrt((E_lab,max - mK - mpi)^2 - mK^2) / q; the third particle is fixed by n3 = d - n1 - n2 and is
filtered by the final energy window.  A small margin (1e-9 in energy) guards the window edges.
"""
from __future__ import annotations

import itertools
from dataclasses import dataclass

import numpy as np

from .sector import Sector

MERGE_TOL = 1e-11          # levels closer than this are one level (prompts/003 section 2.1)
_EDGE_MARGIN = 1e-9        # enumeration margin around [E_min, E_max]; the final filter is exact and inclusive


@dataclass(frozen=True)
class FreeLevel:
    E: float
    multiplicity: int
    states: tuple      # tuple of (n_K1, n_K2, n_pi), each an int 3-tuple, n_K1 <= n_K2 lexicographic


def _kaon_momenta(q, mK, mpi, E_lab_max):
    """All integer vectors n (lexicographic order) that a kaon can carry below E_lab_max."""
    omega_max = E_lab_max - mK - mpi
    if omega_max < mK:
        return np.zeros((0, 3), dtype=np.int64)
    n2max = (omega_max ** 2 - mK ** 2) / (q * q)
    nmax = int(np.floor(np.sqrt(n2max))) + 1
    rng = range(-nmax, nmax + 1)
    cand = np.array(list(itertools.product(rng, repeat=3)), dtype=np.int64)    # itertools order == lexicographic
    n2 = np.einsum("ij,ij->i", cand, cand).astype(float)
    return cand[n2 <= n2max * (1.0 + 1e-12) + 1e-9]


def _raw_states(sector: Sector, E_min: float, E_max: float):
    """[(E*, (n1, n2, n3))] for every state with E_min - margin <= E* <= E_max + margin (unsorted, unmerged)."""
    if not (E_max >= E_min):
        return []
    q = 2.0 * np.pi / (sector.xi * sector.Lbyas)
    d = np.asarray(sector.nnP, dtype=np.int64)
    P2 = q * q * float(d @ d)
    mK, mpi = sector.atmK, sector.atmpi
    E_hi = E_max + _EDGE_MARGIN
    E_lab_max = np.sqrt(E_hi ** 2 + P2)
    kaons = _kaon_momenta(q, mK, mpi, E_lab_max)
    if len(kaons) == 0:
        return []
    wK = np.sqrt(mK ** 2 + q * q * np.einsum("ij,ij->i", kaons, kaons))
    out = []
    for i in range(len(kaons)):
        # unordered pair: j >= i in the lexicographic order of `kaons` <=> n1 <= n2
        n1 = kaons[i]
        n2 = kaons[i:]
        n3 = d - n1 - n2
        E_lab = wK[i] + wK[i:] + np.sqrt(mpi ** 2 + q * q * np.einsum("ij,ij->i", n3, n3))
        Es = np.sqrt(np.maximum(E_lab ** 2 - P2, 0.0))
        keep = np.nonzero((Es >= E_min - _EDGE_MARGIN) & (Es <= E_hi))[0]
        for k in keep:
            out.append((float(Es[k]), (tuple(int(x) for x in n1), tuple(int(x) for x in n2[k]), tuple(int(x) for x in n3[k]))))
    return out


def free_level_states(sector: Sector, E_min: float, E_max: float) -> list:
    """Distinct free levels with E_min <= E* <= E_max (inclusive), sorted ascending, with multiplicities and states."""
    raw = sorted(_raw_states(sector, float(E_min), float(E_max)))
    levels = []
    start = 0
    for k in range(1, len(raw) + 1):
        if k == len(raw) or raw[k][0] - raw[k - 1][0] > MERGE_TOL:
            group = raw[start:k]
            E = float(np.mean([g[0] for g in group]))
            if E_min <= E <= E_max:
                states = tuple(sorted(g[1] for g in group))
                levels.append(FreeLevel(E=E, multiplicity=len(states), states=states))
            start = k
    return levels


def free_levels(sector: Sector, E_min: float, E_max: float) -> np.ndarray:
    """Sorted distinct free-level energies E* in [E_min, E_max] (inclusive), float64."""
    return np.array([lv.E for lv in free_level_states(sector, E_min, E_max)], dtype=np.float64)

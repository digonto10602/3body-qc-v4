"""qc3ref — the planner's independent NumPy reference for the 2+1 (KKpi) three-particle QC ingredients.

Written from the published formulas (Blanton, Romero-Lopez, Sharpe, arXiv:2111.12734, JHEP 02 (2022) 098:
Eqs. 3-15, 20-21, 28-31, App. A), not from the C++ code. Used only by tests/acceptance; never by production code.
See reports/research/007-single-energy-qc-methods-and-reference.md for the validation numbers.
"""
from . import kinematics, basis, fgk, k3  # noqa: F401

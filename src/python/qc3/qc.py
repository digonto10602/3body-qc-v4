"""QC helpers (prompts/001 step 4)."""
from .point import MissingRootError  # noqa: F401


def project_Q(F3inv, Vsel, K3, L):
    """(L^-6) V^dag (F3inv + K3) V ; L = xi*Lbyas."""
    raise NotImplementedError("qc3.qc.project_Q: prompts/001 step 4")


def slogdet(Q):
    """Partial-pivot LU; sign in {-1,+1}; raises on singular."""
    raise NotImplementedError("qc3.qc.slogdet: prompts/001 step 4")

"""Fast backend (prompts/001 steps 6-7)."""
from .point import ALL, PointResult  # noqa: F401
from .sector import Sector  # noqa: F401


def evaluate_point(sector, Ecm, c=(0.0, 0.0, 0.0, 0.0), want=ALL, options=None):
    raise NotImplementedError("qc3.fast.evaluate_point: prompts/001 step 6")

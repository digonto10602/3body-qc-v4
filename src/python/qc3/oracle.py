"""Oracle backend: unchanged physics of build_component_entry (prompts/001 step 2)."""
from .point import ALL, PointResult  # noqa: F401
from .sector import Sector  # noqa: F401


def evaluate_point(sector, Ecm, c=(0.0, 0.0, 0.0, 0.0), want=ALL):
    raise NotImplementedError("qc3.oracle.evaluate_point: prompts/001 step 2")

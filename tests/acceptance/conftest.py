"""Acceptance-test configuration (planner-owned, locked with the acceptance tests).

Makes the production Python package (src/python/qc3) and the planner reference package
(tests/acceptance/m1/qc3ref) importable, and registers the `slow` marker.
No call into the production interface happens at import time.
"""
import pathlib
import sys

ROOT = pathlib.Path(__file__).resolve().parents[2]
for _p in (ROOT / "src" / "python", ROOT / "tests" / "acceptance" / "m1"):
    if str(_p) not in sys.path:
        sys.path.insert(0, str(_p))


def pytest_configure(config):
    config.addinivalue_line("markers", "slow: long-running acceptance test (run at milestone gates with -m slow)")

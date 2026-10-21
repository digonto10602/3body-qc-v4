"""Unit-test configuration: make src/python (package qc3) importable."""
import pathlib
import sys

ROOT = pathlib.Path(__file__).resolve().parents[2]
_p = ROOT / "src" / "python"
if str(_p) not in sys.path:
    sys.path.insert(0, str(_p))

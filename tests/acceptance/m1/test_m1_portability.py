"""M1 step 5 — no machine-specific absolute paths on the M1 production path; QC_DATA_ROOT honoured."""
import os
import pathlib

from m1_helpers import ROOT

FORBIDDEN = ("/home/digonto", "/media/digonto", "/home/digimonk")
SCAN = ["src/python", "src/qc3_core", "build.sh", "configs/qc3", "tests/acceptance"]
EXT = {".py", ".cpp", ".hpp", ".h", ".cc", ".sh", ".in", ".json", ".txt", ".md", ".toml", ".cfg"}


def _files():
    for s in SCAN:
        p = ROOT / s
        if p.is_file():
            yield p
        elif p.is_dir():
            for f in p.rglob("*"):
                if f.is_file() and f.suffix in EXT and "__pycache__" not in f.parts:
                    yield f


def test_no_absolute_home_paths_in_m1_code():
    from qc3.lattice import default_data_dir
    default_data_dir()
    hits = []
    for f in _files():
        try:
            text = f.read_text(errors="ignore")
        except OSError:
            continue
        for bad in FORBIDDEN:
            if bad in text:
                hits.append((str(f.relative_to(ROOT)), bad))
    assert not hits, hits


def test_data_root_env_var(tmp_path, monkeypatch):
    from qc3.lattice import default_data_dir
    monkeypatch.setenv("QC_DATA_ROOT", str(tmp_path))
    assert pathlib.Path(default_data_dir()).resolve() == tmp_path.resolve()
    monkeypatch.delenv("QC_DATA_ROOT", raising=False)
    d = pathlib.Path(default_data_dir()).resolve()
    assert d == (ROOT / "data").resolve()
    assert (d / "lattice" / "L20" / "20_000_A1m_n0.jack").exists()

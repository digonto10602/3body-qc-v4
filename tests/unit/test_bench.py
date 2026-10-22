"""Unit tests for qc3.bench (prompts/001 step 8): CLI parsing, JSON schema, chi^2 proxy = 12 zeros."""
import json
import os

import pytest

from qc3 import bench
from qc3.constants import K2_PWAVE_FORM_DEFAULT

TOP_KEYS = {"provenance", "machine", "env", "constants", "k2_pwave_form", "settings", "seeds", "peak_rss_mb", "summary",
            "points", "chi2_proxy"}
STAGE_KEYS = set(bench.STAGES)


@pytest.fixture
def restore_env():
    saved_env = {v: os.environ.get(v) for v in bench.THREAD_VARS}
    saved_aff = os.sched_getaffinity(0) if hasattr(os, "sched_getaffinity") else None
    yield
    for v, x in saved_env.items():
        if x is None:
            os.environ.pop(v, None)
        else:
            os.environ[v] = x
    if saved_aff is not None:
        os.sched_setaffinity(0, saved_aff)


def test_cli_parses_defaults_and_options():
    ap = bench.build_parser()
    a = ap.parse_args([])
    assert a.out is None and not a.no_oracle and a.threads == 1 and a.k2_pwave_form == K2_PWAVE_FORM_DEFAULT
    a = ap.parse_args(["--out", "x.json", "--no-oracle", "--threads", "4", "--k2-pwave-form", "coded_v2", "--repeats", "2"])
    assert a.out == "x.json" and a.no_oracle and a.threads == 4 and a.k2_pwave_form == "coded_v2" and a.repeats == 2
    with pytest.raises(SystemExit):
        ap.parse_args(["--k2-pwave-form", "nonsense"])


def test_set_thread_env(restore_env):
    env = bench.set_thread_env(1)
    assert env["OMP_NUM_THREADS"] == "1" and env["OPENBLAS_NUM_THREADS"] == "1" and os.environ["OMP_NUM_THREADS"] == "1"
    with pytest.raises(ValueError):
        bench.set_thread_env(0)


def test_json_schema_keys(tmp_path, restore_env):
    out = tmp_path / "b.json"
    res = bench.run(irreps=("100_A2",), energies=(0.2900,), repeats=2, warmup=1, backends=("fast",), out_json=out,
                    include_anchors=False, include_proxy=False, verbose=False)
    on_disk = json.loads(out.read_text())
    assert TOP_KEYS <= set(res) and TOP_KEYS <= set(on_disk)
    assert res["k2_pwave_form"] == K2_PWAVE_FORM_DEFAULT
    assert {"git_commit", "core_built_at_commit", "numpy", "scipy", "date"} <= set(res["provenance"])
    assert {"cpu_model", "cores_online"} <= set(res["machine"])
    assert res["env"]["OMP_NUM_THREADS"] == "1" and res["env"]["OPENBLAS_NUM_THREADS"] == "1"
    assert res["constants"]["production"]["xi"] == 3.444 and len(res["constants"]["c"]) == 4
    (pt,) = res["points"]
    assert pt["irrep"] == "100_A2" and pt["tag"] == "grid"
    f = pt["fast"]
    assert STAGE_KEYS <= set(f["stages"]) and {"median", "min", "max", "n"} <= set(f["stages"]["total"])
    assert f["stages"]["total"]["n"] == 2 and f["wall"]["median"] > 0 and f["N"] > 0
    assert f["cold"]["wall"] > 0 and STAGE_KEYS <= set(f["cold"]["stages"])
    assert res["chi2_proxy"] is None


def test_chi2_proxy_returns_12_zeros():
    p = bench.chi2_proxy(repeats=1)
    assert p["k2_pwave_form"] == "coded_v2" and p["is_a_fit"] is False
    assert p["n_roots"] == 12 and p["n_zero"] == 12 and p["all_kind_zero"]
    assert all(r["kind"] == "zero" for r in p["levels"]) and len(p["levels"]) == 12
    assert all(r["n_evaluations"] >= 5 for r in p["levels"])
    assert p["evaluations_total"] == sum(r["n_evaluations"] for r in p["levels"])
    # the roots lie inside their labelled windows and within 1e-6 (0.2 grid rows) of the labelled coded_v2 zeros
    assert all(r["bracket"][0] <= r["E_root"] <= r["bracket"][1] for r in p["levels"])
    assert max(abs(r["E_root"] - r["previous_model_zero"]) for r in p["levels"]) < 1e-6
    assert p["chi2_at_baseline_c"] > 0 and p["warm"]["wall_s"]["median"] > 0
    assert set(p["reference_times_s"]) == {"old_direct_fcn", "hot_window_fcn"}

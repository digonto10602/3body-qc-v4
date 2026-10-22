"""Benchmark driver and chi^2 proxy (prompts/001 step 8, protocol section 7; prompts/002 section 2.5 and sub-step 7.8).

    OMP_NUM_THREADS=1 python -m qc3.bench --out results/bench_<tag>.json [--no-oracle] [--threads 1]
                                          [--k2-pwave-form q0_consistent|coded_v2]

What is measured (everything single core unless --threads says otherwise):
  * points: 5 irreps x E* in {0.2700, 0.3100, 0.3599} plus the 12 oracle-anchor points (L20), c = C_BASELINE;
  * per point and backend: warm-up calls, then `repeats` timed calls; median / min / max of the wall time and of every
    per-stage entry of PointResult.timings; for the fast backend also one cold call (after qc3.fast.clear_caches());
  * speed-up = median(oracle total) / median(fast warm total);
  * chi^2 proxy: wall time and number of evaluations of qc3.roots.refine_root_info(kind="zero") (= refine_root with the
    pole/zero classifier) over the 12 labelled L20 roots (cutoff 0.335, brackets from data/windows/L20/*/accepted_windows.csv)
    on the fast backend, k2_pwave_form pinned to "coded_v2" (the labels are coded_v2 zeros), c = C_BASELINE, plus the
    chi^2 of those 12 levels against the lattice (qc3.lattice, D1 covariance, D3 matching).
    THE chi^2 NUMBER IS NOT A FIT: c is the fixed baseline label-generation vector, nothing is minimised.

Thread pinning: the BLAS/OpenMP thread variables are read when numpy is imported, so this module imports numpy lazily and
`main()` sets the variables before anything imports numpy.  `run()` called from a process that already imported numpy cannot
change them any more; the JSON records `numpy_imported_before_env` and the measured CPU/wall ratio of the proxy to show it.
"""
from __future__ import annotations

import argparse
import csv
import datetime
import json
import os
import pathlib
import platform
import subprocess
import sys
import time

from .constants import C_BASELINE, IRREPS, K2_PWAVE_FORM_DEFAULT, PRODUCTION

ROOT = pathlib.Path(__file__).resolve().parents[3]
ANCHOR_FILE = ROOT / "tests" / "acceptance" / "m1" / "data" / "oracle_anchors_L20.json"
STAGES = ("basis", "F", "G", "K2inv", "F3", "F3inv", "projector", "K3", "projection", "lu", "total")
THREAD_VARS = ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS", "VECLIB_MAXIMUM_THREADS")
LBYAS = 20
CUTOFF = 0.335
PROXY_K2_FORM = "coded_v2"          # the 12 labels are coded_v2 zeros (prompts/002 section 2.5 / 7.8)
REFERENCE_TIMES_S = {"old_direct_fcn": 1.4, "hot_window_fcn": 0.024}   # prompts/001 step 8; hot window needed GBs of caches
SEEDS = {"stochastic_inputs": None, "note": "deterministic: fixed energies and windows, no random numbers"}


# ---------------------------------------------------------------- environment / provenance

def set_thread_env(threads: int) -> dict:
    """Force every BLAS/OpenMP thread variable to `threads` in os.environ; returns what was set."""
    threads = int(threads)
    if threads < 1:
        raise ValueError("threads must be >= 1")
    for v in THREAD_VARS:
        os.environ[v] = str(threads)
    return {v: os.environ[v] for v in THREAD_VARS}


def pin_cores(threads: int, core: int | None = None) -> dict:
    """Pin the process to `threads` cores (os.sched_setaffinity, like taskset -c).  Returns before/after sets."""
    if not hasattr(os, "sched_setaffinity"):
        return {"available": False, "before": None, "after": None}
    before = sorted(os.sched_getaffinity(0))
    if core is None:
        chosen = before[: max(1, int(threads))]
    else:
        chosen = [core + k for k in range(max(1, int(threads)))]
    if not set(chosen) <= set(before):
        raise ValueError(f"requested cores {chosen} are not in the allowed set {before}")
    os.sched_setaffinity(0, set(chosen))
    return {"available": True, "before": before, "after": sorted(os.sched_getaffinity(0))}


def _cpu_model() -> str:
    try:
        with open("/proc/cpuinfo") as fh:
            for line in fh:
                if line.startswith("model name"):
                    return line.split(":", 1)[1].strip()
    except OSError:
        pass
    return platform.processor() or "unknown"


def _git(*args) -> str:
    try:
        return subprocess.run(["git", "-C", str(ROOT), *args], capture_output=True, text=True, check=False).stdout.strip()
    except OSError:
        return ""


def machine_info() -> dict:
    return {"cpu_model": _cpu_model(), "cores_online": os.cpu_count(),
            "cores_allowed": len(os.sched_getaffinity(0)) if hasattr(os, "sched_getaffinity") else None,
            "platform": platform.platform(), "python": sys.version.split()[0]}


def provenance() -> dict:
    import numpy
    try:
        import scipy
        scipy_v = scipy.__version__
    except ImportError:   # scipy is a hard dependency of qc3.roots, so this cannot normally happen
        scipy_v = None
    from . import _core
    return {"git_commit": _git("rev-parse", "--short", "HEAD"),
            "git_tree_dirty": bool(_git("status", "--porcelain", "--", "src", "tests", "scripts", "data")),
            "core_built_at_commit": _core.git_commit(),   # fixed when _core was compiled (may differ from HEAD)
            "core_version": str(_core.version()),
            "numpy": numpy.__version__, "scipy": scipy_v,
            "date": datetime.datetime.now(datetime.timezone.utc).isoformat(timespec="seconds"),
            "build_flags": "see build.sh COMMON (-std=c++17 -O3 -fopenmp); not introspectable from the module"}


def _peak_rss_mb() -> float:
    import resource
    return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024.0   # Linux: kB


# ---------------------------------------------------------------- statistics

def _stats(xs) -> dict:
    import numpy as np
    a = np.asarray(xs, dtype=float)
    return {"median": float(np.median(a)), "min": float(a.min()), "max": float(a.max()), "n": int(a.size)}


def _time_point(backend_name, evaluate, sector, E, c, repeats, warmup):
    """Warm-up + repeats of one point; returns {"wall": stats, "stages": {stage: stats}, "N", "Nproj", "metadata_subset"}."""
    for _ in range(int(warmup)):
        evaluate(sector, E, c, want=("Q",))
    walls, stages = [], {s: [] for s in STAGES}
    last = None
    for _ in range(int(repeats)):
        t0 = time.perf_counter()
        r = evaluate(sector, E, c, want=("Q",))
        walls.append(time.perf_counter() - t0)
        for s in STAGES:
            stages[s].append(float(r.timings.get(s, 0.0)))
        last = r
    return {"wall": _stats(walls), "stages": {s: _stats(v) for s, v in stages.items()},
            "N": int(last.N), "Nproj": int(last.Nproj), "sign": int(last.sign), "lnabsdet": float(last.lnabsdet)}


# ---------------------------------------------------------------- chi^2 proxy

def read_windows(data_dir=None, Lbyas=LBYAS) -> dict:
    """{(irrep, lattice_level_index): csv row} for the rows with inside_Ecm_cutoff == true (cutoff 0.335 labels)."""
    data_dir = pathlib.Path(data_dir) if data_dir is not None else ROOT / "data"
    out = {}
    for irrep in IRREPS:
        p = data_dir / "windows" / f"L{Lbyas}" / f"L{Lbyas}_{irrep}" / "accepted_windows.csv"
        with open(p) as fh:
            for r in csv.DictReader(fh):
                if r["inside_Ecm_cutoff"].strip().lower() == "true":
                    out[(irrep, int(r["lattice_level_index"]))] = r
    return out


def _proxy_pass(sectors, windows, labels, c, evaluate):
    """One pass over the 12 roots; returns (rows, wall, cpu)."""
    from .roots import refine_root_info
    rows = []
    t0, p0 = time.perf_counter(), time.process_time()
    for key in labels:
        w = windows[key]
        left, right = float(w["E_left_bracket"]), float(w["E_right_bracket"])
        info = refine_root_info(evaluate, sectors[key[0]], c, left, right, kind="zero")
        rows.append({"irrep": key[0], "level_index": key[1], "E_root": info.E, "kind": info.kind,
                     "discriminant": info.discriminant, "n_evaluations": info.n_evaluations, "N": info.N,
                     "bracket": [left, right], "previous_model_zero": float(w["previous_model_zero"])})
    return rows, time.perf_counter() - t0, time.process_time() - p0


def chi2_proxy(c=C_BASELINE, repeats=3, cutoff=CUTOFF, data_dir=None, Lbyas=LBYAS, k2_pwave_form=PROXY_K2_FORM) -> dict:
    """Time the refinement of the 12 labelled L20 roots on the fast backend (kind="zero"), cold once then `repeats` warm
    passes, and compute the chi^2 of the 12 levels at the fixed baseline c.  NOT a fit (nothing is minimised).

    A root that fails (PoleError, AmbiguousRootError, MissingRootError, ...) is not caught: it propagates, the run stops."""
    import numpy as np

    from . import fast, lattice
    from .sector import Sector

    data_dir = pathlib.Path(data_dir) if data_dir is not None else ROOT / "data"
    c = tuple(float(x) for x in c)
    windows = read_windows(data_dir, Lbyas)
    ls = lattice.load_levels(data_dir, Lbyas, list(IRREPS), cutoff)
    labels = [tuple(l) for l in ls.labels]
    if set(labels) != set(windows) or len(labels) != 12:
        raise RuntimeError(f"windows/levels mismatch at cutoff {cutoff}: {sorted(set(labels) ^ set(windows))}, n={len(labels)}")
    sectors = {irrep: Sector.production(irrep, Lbyas, k2_pwave_form=k2_pwave_form) for irrep in IRREPS}

    fast.clear_caches()
    rows_cold, wall_cold, cpu_cold = _proxy_pass(sectors, windows, labels, c, fast.evaluate_point)
    warm = []
    for _ in range(int(repeats)):
        rows, wall, cpu = _proxy_pass(sectors, windows, labels, c, fast.evaluate_point)
        warm.append((rows, wall, cpu))
    rows = warm[-1][0] if warm else rows_cold
    walls = [w for _, w, _ in warm] or [wall_cold]
    cpus = [x for _, _, x in warm] or [cpu_cold]

    # chi^2 at the fixed baseline c: D3 matching per irrep (one root per level here, so the pairing is the identity)
    E_lat = np.asarray(ls.Ecm, dtype=float)
    E_qc = np.full(len(labels), np.nan)
    for irrep in IRREPS:
        idx = [i for i, l in enumerate(labels) if l[0] == irrep]
        if not idx:
            continue
        lat = E_lat[idx]
        qc_roots = np.array([next(r["E_root"] for r in rows if (r["irrep"], r["level_index"]) == labels[i]) for i in idx])
        for a, b in lattice.match_levels(lat, qc_roots):
            E_qc[idx[a]] = qc_roots[b]
    resid = E_qc - E_lat
    sigma = np.sqrt(np.diag(ls.cov))
    chi2 = lattice.chi2(resid, ls.cov)
    by_label = {(r["irrep"], r["level_index"]): r for r in rows}
    levels = []
    for i, l in enumerate(labels):
        r = dict(by_label[l])
        r.update({"lattice_Ecm": float(E_lat[i]), "lattice_sigma": float(sigma[i]),
                  "residual_qc_minus_lattice": float(resid[i]), "residual_over_sigma": float(resid[i] / sigma[i])})
        levels.append(r)
    ev_total = sum(r["n_evaluations"] for r in rows)
    wall_med = float(np.median(walls))
    return {
        "description": "time to refine the 12 labelled L20 roots (kind='zero', fast backend, 1 core) + chi^2 at FIXED baseline c",
        "is_a_fit": False,
        "note": "chi2 is evaluated at the label-generation vector C_BASELINE; no parameter is varied or minimised. "
                "It is the cost of ONE chi2 evaluation, not a fitted chi2.",
        "k2_pwave_form": k2_pwave_form, "c": list(c), "Lbyas": Lbyas, "Ecm_cutoff": cutoff,
        "n_roots": len(rows), "n_zero": sum(1 for r in rows if r["kind"] == "zero"),
        "all_kind_zero": all(r["kind"] == "zero" for r in rows),
        "evaluations_total": ev_total, "evaluations_per_root": [r["n_evaluations"] for r in rows],
        "evaluations_per_root_median": float(np.median([r["n_evaluations"] for r in rows])),
        "cold": {"wall_s": wall_cold, "cpu_s": cpu_cold},
        "warm": {"wall_s": _stats(walls), "cpu_s": _stats(cpus), "cpu_over_wall_median": float(np.median(cpus) / wall_med),
                 "per_evaluation_ms_median": 1e3 * wall_med / ev_total},
        "reference_times_s": dict(REFERENCE_TIMES_S),
        "ratio_vs_reference": {"old_direct_fcn_over_proxy_warm": REFERENCE_TIMES_S["old_direct_fcn"] / wall_med,
                               "hot_window_fcn_over_proxy_warm": REFERENCE_TIMES_S["hot_window_fcn"] / wall_med},
        "chi2_at_baseline_c": chi2, "n_levels": len(labels), "chi2_covariance": "D1 (jackknife of per-configuration values)",
        "Ncfg": int(ls.Ncfg), "levels": levels,
    }


# ---------------------------------------------------------------- the benchmark

def _points(irreps, energies, include_anchors):
    pts = [{"tag": "grid", "irrep": ir, "Ecm": float(E)} for ir in irreps for E in energies]
    if include_anchors:
        if not ANCHOR_FILE.exists():
            raise FileNotFoundError(f"anchor file {ANCHOR_FILE} not found (use include_anchors=False to skip the 12 anchors)")
        with open(ANCHOR_FILE) as fh:
            for p in json.load(fh)["points"]:
                pts.append({"tag": "anchor", "irrep": p["irrep"], "Ecm": float(p["Ecm"])})
    return pts


def run(irreps=IRREPS, energies=(0.2700, 0.3100, 0.3599), repeats=5, warmup=1, threads=1, backends=("oracle", "fast"),
        c=C_BASELINE, out_json=None, k2_pwave_form=K2_PWAVE_FORM_DEFAULT, include_anchors=True, include_proxy=True,
        oracle_repeats=1, oracle_warmup=0, proxy_repeats=3, cold=True, core=None, verbose=True) -> dict:
    """Run the benchmark of the module docstring; returns the result dict and writes it to `out_json` if given.
    The process is pinned to `threads` core(s) for the duration of the call; the previous affinity is restored afterwards
    (the thread environment variables are left set)."""
    affinity_before = sorted(os.sched_getaffinity(0)) if hasattr(os, "sched_getaffinity") else None
    try:
        return _run(irreps, energies, repeats, warmup, threads, backends, c, out_json, k2_pwave_form, include_anchors,
                    include_proxy, oracle_repeats, oracle_warmup, proxy_repeats, cold, core, verbose)
    finally:
        if affinity_before is not None:
            os.sched_setaffinity(0, set(affinity_before))


def _run(irreps, energies, repeats, warmup, threads, backends, c, out_json, k2_pwave_form, include_anchors, include_proxy,
         oracle_repeats, oracle_warmup, proxy_repeats, cold, core, verbose) -> dict:
    t_start = time.perf_counter()
    numpy_before = "numpy" in sys.modules
    env = set_thread_env(threads)
    affinity = pin_cores(threads, core)
    import numpy as np

    from . import fast, oracle
    from .sector import Sector

    c = tuple(float(x) for x in c)
    unknown = set(backends) - {"oracle", "fast"}
    if unknown:
        raise ValueError(f"unknown backends {sorted(unknown)}")
    evals = {"fast": fast.evaluate_point, "oracle": oracle.evaluate_point}
    pts = _points(irreps, energies, include_anchors)
    log = (lambda *a: print(*a, flush=True)) if verbose else (lambda *a: None)

    results = [dict(p) for p in pts]
    peak_rss = {}
    for backend in ("fast", "oracle"):          # fast first so that its peak RSS is not hidden by the oracle's
        if backend not in backends:
            continue
        reps, wup = (repeats, warmup) if backend == "fast" else (oracle_repeats, oracle_warmup)
        for rec in results:
            sec = Sector.production(rec["irrep"], LBYAS, k2_pwave_form=k2_pwave_form)
            ev = evals[backend]
            entry = {}
            if backend == "fast" and cold:
                fast.clear_caches()
                t0 = time.perf_counter()
                r = ev(sec, rec["Ecm"], c, want=("Q",))
                entry["cold"] = {"wall": time.perf_counter() - t0, "stages": {s: float(r.timings.get(s, 0.0)) for s in STAGES}}
            entry.update(_time_point(backend, ev, sec, rec["Ecm"], c, reps, wup))
            rec[backend] = entry
            rec.setdefault("N", entry["N"]); rec.setdefault("Nproj", entry["Nproj"])
            log(f"[{backend}] {rec['tag']:6s} {rec['irrep']} E={rec['Ecm']:.6f} N={entry['N']:3d} Nproj={entry['Nproj']:3d} "
                f"median wall {1e3 * entry['wall']['median']:9.2f} ms")
        peak_rss[backend] = _peak_rss_mb()

    for rec in results:
        if "fast" in rec and "oracle" in rec:
            fo, ff = rec["oracle"], rec["fast"]
            rec["speedup_total"] = fo["stages"]["total"]["median"] / ff["stages"]["total"]["median"]
            rec["speedup_wall"] = fo["wall"]["median"] / ff["wall"]["median"]
            proj = fo["stages"]["projector"]["median"]
            rec["speedup_total_without_oracle_projector"] = (fo["stages"]["total"]["median"] - proj) / ff["stages"]["total"]["median"]
            rec["sign_lnabsdet_agree"] = {"sign": fo["sign"] == ff["sign"], "dln": abs(fo["lnabsdet"] - ff["lnabsdet"])}

    summary = {}
    if all("speedup_total" in r for r in results) and results and "oracle" in backends and "fast" in backends:
        for key in ("speedup_total", "speedup_wall", "speedup_total_without_oracle_projector"):
            summary[key] = {"min": float(min(r[key] for r in results)), "median": float(np.median([r[key] for r in results])),
                            "max": float(max(r[key] for r in results))}
    if "fast" in backends and results:
        summary["fast_warm_total_ms"] = {"min": 1e3 * min(r["fast"]["stages"]["total"]["median"] for r in results),
                                         "max": 1e3 * max(r["fast"]["stages"]["total"]["median"] for r in results)}
        if cold:
            summary["fast_cold_minus_warm_ms_max"] = 1e3 * max(r["fast"]["cold"]["wall"] - r["fast"]["wall"]["median"] for r in results)

    proxy = None
    if include_proxy:
        log("[proxy] refining the 12 labelled roots (coded_v2, fast, kind='zero') ...")
        proxy = chi2_proxy(c, repeats=proxy_repeats)
        peak_rss["after_proxy"] = _peak_rss_mb()
        log(f"[proxy] n_zero={proxy['n_zero']}/{proxy['n_roots']} evals={proxy['evaluations_total']} "
            f"cold {proxy['cold']['wall_s']:.3f} s, warm median {proxy['warm']['wall_s']['median']:.3f} s "
            f"(cpu/wall {proxy['warm']['cpu_over_wall_median']:.2f}); chi2(fixed baseline c, NOT a fit) = {proxy['chi2_at_baseline_c']:.4f}")

    out = {
        "kind": "qc3.bench", "protocol": "prompts/001 section 7 (+ prompts/002 section 2.5, sub-step 7.8)",
        "provenance": provenance(), "machine": machine_info(),
        "env": {**env, "threads_requested": int(threads), "numpy_imported_before_env": numpy_before,
                "env_applies_to_blas": not numpy_before, "affinity": affinity,
                "fast_threads_option": 1, "taskset_equivalent": "os.sched_setaffinity"},
        "constants": {"production": {k: (list(v) if isinstance(v, tuple) else v) for k, v in PRODUCTION.items()},
                      "c": list(c), "Lbyas": LBYAS, "grid": {"E0": 0.26310, "dE": 0.0969 / 19999, "n": 20000}},
        "k2_pwave_form": k2_pwave_form, "k2_pwave_form_proxy": PROXY_K2_FORM,
        "settings": {"irreps": list(irreps), "energies": [float(e) for e in energies], "include_anchors": bool(include_anchors),
                     "backends": list(backends), "fast_repeats": int(repeats), "fast_warmup": int(warmup),
                     "oracle_repeats": int(oracle_repeats), "oracle_warmup": int(oracle_warmup), "cold": bool(cold),
                     "proxy_repeats": int(proxy_repeats), "include_proxy": bool(include_proxy)},
        "seeds": SEEDS, "peak_rss_mb": peak_rss, "summary": summary, "points": results, "chi2_proxy": proxy,
        "total_seconds": time.perf_counter() - t_start,
    }
    if out_json is not None:
        p = pathlib.Path(out_json)
        p.parent.mkdir(parents=True, exist_ok=True)
        p.write_text(json.dumps(out, indent=2))
        log(f"wrote {p}")
    return out


def build_parser() -> argparse.ArgumentParser:
    ap = argparse.ArgumentParser(prog="python -m qc3.bench", description=__doc__.split("\n\n")[0])
    ap.add_argument("--out", default=None, help="JSON output path, e.g. results/bench_<tag>.json")
    ap.add_argument("--no-oracle", action="store_true", help="skip the (slow) oracle backend")
    ap.add_argument("--threads", type=int, default=1, help="OMP/BLAS threads and pinned cores (1 = the metric)")
    ap.add_argument("--k2-pwave-form", default=K2_PWAVE_FORM_DEFAULT, choices=("coded_v2", "q0_consistent"),
                    help="p-wave K2^-1 form of the timing sectors (the chi2 proxy always pins coded_v2)")
    ap.add_argument("--repeats", type=int, default=5, help="timed repeats per point for the fast backend")
    ap.add_argument("--warmup", type=int, default=1, help="warm-up calls per point for the fast backend")
    ap.add_argument("--oracle-repeats", type=int, default=1)
    ap.add_argument("--proxy-repeats", type=int, default=3)
    ap.add_argument("--core", type=int, default=None, help="first core to pin to (default: first allowed core)")
    ap.add_argument("--no-anchors", action="store_true", help="skip the 12 anchor points")
    ap.add_argument("--no-proxy", action="store_true", help="skip the chi2 proxy")
    ap.add_argument("--no-cold", action="store_true", help="skip the cold fast call")
    return ap


def main(argv=None):
    args = build_parser().parse_args(argv)
    # must happen before numpy is imported (BLAS reads the variables at import)
    return run(threads=args.threads, backends=("fast",) if args.no_oracle else ("oracle", "fast"), out_json=args.out,
               k2_pwave_form=args.k2_pwave_form, repeats=args.repeats, warmup=args.warmup, oracle_repeats=args.oracle_repeats,
               proxy_repeats=args.proxy_repeats, core=args.core, include_anchors=not args.no_anchors,
               include_proxy=not args.no_proxy, cold=not args.no_cold)


if __name__ == "__main__":
    main()

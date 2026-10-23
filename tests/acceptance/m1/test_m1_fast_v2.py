"""M1 step 7 (prompts/002): fast path v2 — fast-owned projector cache with a complete key, orbit/character projector,
threshold basis builder and dimension jumps, C++ tail, single-core behaviour, speed.

Tolerances and their justification: prompts/002 section 4.  Nothing here touches the locked step-6 tests.
Every production call goes through qc3.fast (and qc3.oracle for the oracle subspace); the reference basis rule comes
from the planner's qc3ref package.  All energies in a_t^{-1}.
"""
import time

import numpy as np
import pytest

from m1_helpers import C_BASELINE, C_ZERO, IRREPS, L_BOX, NNP, PRODUCTION, anchor_points, as_c, oracle_point, sector

TOL_SUBSPACE = 1e-12     # ||V V^dag - V' V'^dag||_F / ||V' V'^dag||_F ; entries of P are multiples of 1/|LG|, measured 6.5e-16
TOL_ORTH = 1e-12         # ||V^dag V - 1||_F ; measured 1.2e-15
TOL_LNDET = 1e-6         # same LU round-off argument as test_m1_fast_vs_oracle (floor 1.4e-9, measured 1.9e-9 worst)
TOL_JUMP = 1e-12         # relative, on jump energies: same float64 formula, different operation order (~1e-15)
JUMP_EPS = 2e-9          # energy offset used to see N change across a jump (grid step is 4.845e-6)
COLD_PROJECTOR_MAX = 0.02    # s; orbit build is O(|LG| N) + a 103x3 SVD (microseconds); old value 5.4 s
COLD_EXTRA_MAX = 0.1         # s; cold - warm total: threshold tables (~21 ms once per sector) + projector build
SPEEDUP_MIN = 150.0          # 111_A2 0.35103 (N=160): expected ~400x, margin 2.7x
E_MIN, E_MAX = 0.26310, 0.36000
NMAX = 20
Z_MIN = 0.02


def _fast(irrep, Ecm, c=C_ZERO, want=("F", "G", "K2inv", "F3", "F3inv", "Vsel", "K3", "Q"), **options):
    """Fast-path call; the projector method is always passed explicitly (default "orbit", the production default)."""
    from qc3 import fast
    opts = {"projector": "orbit"}
    opts.update(options)
    return fast.evaluate_point(sector(irrep), float(Ecm), tuple(c), want=want, options=opts)


def _proj(V):
    V = as_c(V)
    return V @ V.conj().T


def _rel(a, b):
    nb = np.linalg.norm(b)
    return float(np.linalg.norm(a - b) / nb) if nb > 0 else float(np.linalg.norm(a - b))


def _thresholds_rule(irrep):
    """E_lab thresholds E_min^+(n) of docs/METHODS.md section 2 for both channels, loop order, from the rule (not from qc3)."""
    mK, mpi = PRODUCTION["atmK"], PRODUCTION["atmpi"]
    P = 2 * np.pi / L_BOX * np.array(NNP[irrep], dtype=float)
    g = np.arange(-NMAX, NMAX + 1)
    nx, ny, nz = np.meshgrid(g, g, g, indexing="ij")
    k = 2 * np.pi / L_BOX * np.stack([nx.ravel(), ny.ravel(), nz.ravel()], axis=1).astype(float)
    out = []
    for (mi, mj, mk) in ((mK, mK, mpi), (mpi, mK, mK)):   # channel 1 (K spectator), channel 2 (pi spectator)
        sig_min = abs(mj * mj - mk * mk) + 2.0 * Z_MIN * (mj + mk) * min(mj, mk)
        omg = np.sqrt((k * k).sum(axis=1) + mi * mi)
        out.append(omg + np.sqrt(sig_min + ((P - k) ** 2).sum(axis=1)))
    return np.concatenate(out), P


def _jumps_rule(irrep, E_min=E_MIN, E_max=E_MAX):
    """Sorted distinct E_cm (strictly inside (E_min, E_max)) at which the basis changes, merged within 1e-12 relative."""
    t, P = _thresholds_rule(irrep)
    Psq = float(P @ P)
    lo, hi = np.sqrt(E_min ** 2 + Psq), np.sqrt(E_max ** 2 + Psq)
    inside = np.sort(t[(t > lo) & (t < hi)])
    merged = []
    for v in inside:
        if not merged or abs(v - merged[-1]) > 1e-12 * v:
            merged.append(float(v))
    return np.sqrt(np.array(merged) ** 2 - Psq)


# ----------------------------------------------------------------------------------------------------------------------
# 7.0 — complete cache key (fast-owned projector cache)
# ----------------------------------------------------------------------------------------------------------------------

def test_projector_cache_key_complete():
    """Report 006 note 2: the legacy key omitted eig_tol, so eig_tol=1e-30 after a default call returned the default's 17
    columns instead of the 5 a cold build gives.  The fast-owned cache must make the result order-independent."""
    from qc3 import fast
    irrep, Ecm = "100_A2", 0.2915561
    o = oracle_point(irrep, Ecm, C_ZERO)

    def nproj(eig_tol):
        return fast.evaluate_point(sector(irrep, eig_tol=eig_tol), Ecm, C_ZERO, want=("Vsel",),
                                   options={"projector": "eigen"}).Nproj

    fast.clear_caches()
    a1 = nproj(0.05)
    b1 = nproj(1e-30)
    info1 = fast.cache_info()
    fast.clear_caches()
    b2 = nproj(1e-30)
    a2 = nproj(0.05)
    info2 = fast.cache_info()
    assert a1 == a2 == o.Nproj, (a1, a2, o.Nproj)
    assert b1 == b2, (b1, b2)
    assert b1 != a1, "eig_tol=1e-30 must not reproduce the default selection (report 006: 5 vs 17)"
    assert info1["projector_entries"] == 2 and info2["projector_entries"] == 2, (info1, info2)
    assert info2["projector_misses"] >= 2


# ----------------------------------------------------------------------------------------------------------------------
# 7.2 — orbit projector
# ----------------------------------------------------------------------------------------------------------------------

def test_projector_orbit_matches_eigen_and_oracle_subspace():
    pts = anchor_points(small_only=True)
    assert len(pts) >= 3
    for p in pts:
        irrep, Ecm = p["irrep"], p["Ecm"]
        fo = _fast(irrep, Ecm, C_BASELINE, projector="orbit")
        fe = _fast(irrep, Ecm, C_BASELINE, projector="eigen")
        o = oracle_point(irrep, Ecm, C_BASELINE)
        tag = (irrep, Ecm)
        assert fo.metadata["projector_method"] == "orbit" and fe.metadata["projector_method"] == "eigen", tag
        assert fo.Nproj == fe.Nproj == o.Nproj, (tag, fo.Nproj, fe.Nproj, o.Nproj)
        assert fo.Vsel.shape == (fo.N, fo.Nproj), tag
        Vo = as_c(fo.Vsel)
        assert np.linalg.norm(Vo.conj().T @ Vo - np.eye(fo.Nproj)) <= TOL_ORTH, tag
        Pe, Po, Por = _proj(fe.Vsel), _proj(fo.Vsel), _proj(o.Vsel)
        assert _rel(Po, Pe) <= TOL_SUBSPACE, (tag, "orbit vs eigen", _rel(Po, Pe))
        assert _rel(Po, Por) <= TOL_SUBSPACE, (tag, "orbit vs oracle", _rel(Po, Por))
        assert int(fo.sign) == int(fe.sign) == int(o.sign), tag
        assert abs(fo.lnabsdet - fe.lnabsdet) <= TOL_LNDET, (tag, fo.lnabsdet, fe.lnabsdet)
        assert abs(fo.lnabsdet - o.lnabsdet) <= TOL_LNDET, (tag, fo.lnabsdet, o.lnabsdet)


def test_projector_cold_build_is_fast():
    """The 5.4 s case of report 006 (000_A1m 0.34538, N=103, |LG|=48).  Cold projector build and cold-minus-warm total."""
    from qc3 import fast
    irrep, Ecm = "000_A1m", 0.34538
    fast.clear_caches()
    cold = _fast(irrep, Ecm, C_BASELINE, want=("Q",), projector="orbit")
    assert cold.metadata["projector_cache_hit"] is False
    assert cold.N == 103 and cold.Nproj == 9, (cold.N, cold.Nproj)
    assert cold.timings["projector"] <= COLD_PROJECTOR_MAX, cold.timings
    warm = [_fast(irrep, Ecm, C_BASELINE, want=("Q",), projector="orbit") for _ in range(3)]
    assert all(w.metadata["projector_cache_hit"] is True for w in warm)
    t_warm = float(np.median([w.timings["total"] for w in warm]))
    assert cold.timings["total"] - t_warm <= COLD_EXTRA_MAX, (cold.timings["total"], t_warm)
    assert abs(cold.lnabsdet - warm[0].lnabsdet) <= TOL_LNDET and cold.sign == warm[0].sign


def test_fast_options_projector_values():
    from qc3 import fast
    s = sector("100_A2")
    with pytest.raises(ValueError):
        fast.evaluate_point(s, 0.2915561, C_ZERO, want=("Q",), options={"projector": "bogus"})
    r = fast.evaluate_point(s, 0.2915561, C_ZERO, want=("Q",))
    for key in ("projector_method", "projector_cache_hit", "tail", "trunc_eps", "use_classes", "threads"):
        assert key in r.metadata, key
    assert r.metadata["projector_method"] in ("orbit", "eigen")


# ----------------------------------------------------------------------------------------------------------------------
# 7.3 — threshold basis builder and dimension jumps
# ----------------------------------------------------------------------------------------------------------------------

def test_dimension_jumps_match_rule_and_N():
    from qc3 import fast
    total = 0
    for irrep in IRREPS:
        got = np.asarray(fast.dimension_jumps(sector(irrep), E_MIN, E_MAX), dtype=float)
        exp = _jumps_rule(irrep)
        assert got.ndim == 1 and np.all(np.diff(got) > 0), (irrep, got)
        assert np.all((got > E_MIN) & (got < E_MAX)), irrep
        assert got.shape == exp.shape, (irrep, got.shape, exp.shape)
        assert np.max(np.abs(got - exp) / exp) <= TOL_JUMP, (irrep, np.max(np.abs(got - exp) / exp))
        for Ej in got:
            n_lo = fast.evaluate_point(sector(irrep), float(Ej - JUMP_EPS), C_ZERO, want=("Q",)).N
            n_hi = fast.evaluate_point(sector(irrep), float(Ej + JUMP_EPS), C_ZERO, want=("Q",)).N
            assert n_hi > n_lo, (irrep, Ej, n_lo, n_hi)
        total += len(got)
    assert total >= 40, total   # research 008: 4+11+15+12+11 = 53 jumps (58 bases) over the 5 irreps


@pytest.mark.slow
def test_fast_basis_all_signatures_vs_reference():
    """Every basis the production range can produce: the fast path's basis equals the reference rule at the midpoint
    of each constant-N segment (exact integers)."""
    from qc3 import fast
    from qc3ref.basis import expected_basis
    n_bases = 0
    for irrep in IRREPS:
        jumps = np.asarray(fast.dimension_jumps(sector(irrep), E_MIN, E_MAX), dtype=float)
        edges = np.concatenate([[E_MIN], jumps, [E_MAX]])
        mids = 0.5 * (edges[1:] + edges[:-1])
        P = 2 * np.pi / L_BOX * np.array(NNP[irrep], dtype=float)
        last_N = -1
        for Ecm in mids:
            r = fast.evaluate_point(sector(irrep), float(Ecm), C_ZERO, want=("Q",))
            E_lab = float(np.sqrt(Ecm ** 2 + P @ P))
            exp = expected_basis(E_lab, P, L_BOX, PRODUCTION["atmK"], PRODUCTION["atmpi"])
            assert len(r.basis) == len(exp) == r.N, (irrep, Ecm, len(r.basis), len(exp))
            for name in exp.dtype.names:
                assert np.array_equal(np.asarray(r.basis[name]), exp[name]), (irrep, Ecm, name)
            assert r.N > last_N, (irrep, Ecm, r.N, last_N)   # N strictly increases from segment to segment
            last_N = r.N
            n_bases += 1
    assert n_bases >= 45, n_bases


# ----------------------------------------------------------------------------------------------------------------------
# 7.6 — C++ tail, want subsets, single core, speed
# ----------------------------------------------------------------------------------------------------------------------

def test_fast_tail_consistent_with_numpy_and_want_subset():
    from qc3.qc import slogdet
    for p in anchor_points(small_only=True):
        irrep, Ecm = p["irrep"], p["Ecm"]
        full = _fast(irrep, Ecm, C_BASELINE)
        sub = _fast(irrep, Ecm, C_BASELINE, want=("Q",))
        tag = (irrep, Ecm)
        assert full.metadata["tail"] == "cpp" and sub.metadata["tail"] == "cpp", tag
        assert full.F3inv is not None and sub.F3inv is None and sub.Q is not None, tag
        s_np, l_np = slogdet(as_c(full.Q))
        assert int(full.sign) == s_np and abs(full.lnabsdet - l_np) <= TOL_LNDET, (tag, full.lnabsdet, l_np)
        assert int(sub.sign) == int(full.sign) and abs(sub.lnabsdet - full.lnabsdet) <= TOL_LNDET, (tag, sub.lnabsdet, full.lnabsdet)
        assert sub.Q.shape == (sub.Nproj, sub.Nproj) and _rel(as_c(sub.Q), as_c(full.Q)) <= 1e-9, tag


CPU_WALL_MAX = 1.3            # unchanged bar (prompts/002 section 4); measured 0.97-1.00 in steady state (reports/tests/010)
CPU_WORK_WINDOW = 0.6         # s of fast-path work per point; >= 2x the 0.1-0.3 s OpenBLAS spin-down, so a residual spin cannot
                              # push a one-core loop over the bar by itself (1 extra core for 0.3 s of a 0.6 s window = 1.5 only
                              # if it were all inside the window, which the settle step below excludes)
CPU_SETTLE_MAX = 3.0          # s; upper limit on the idle wait for the process to go quiet (measured need: <= 0.3 s)
CPU_POINTS = (("111_A2", 0.35103), ("000_A1m", 0.34538))   # N = 160 and N = 103: where the pre-7.6 NumPy tail DID thread
                                                            # (unpinned CPU/wall 3.06 and 2.29, reports/tests/009)


def _settle_idle(max_wait=CPU_SETTLE_MAX, step=0.1, quiet=0.05):
    """Wait until the process is idle: an idle `step` window uses < `quiet` of a core.  OpenBLAS worker threads spin for
    0.1-0.3 s after `import numpy` (and after a BLAS call); this makes the measurement independent of test order."""
    waited = 0.0
    while waited < max_wait:
        c0, w0 = time.process_time(), time.perf_counter()
        time.sleep(step)
        dc, dw = time.process_time() - c0, time.perf_counter() - w0
        waited += dw
        if dc / dw < quiet:
            return waited
    return waited


def _cpu_wall_ratio(irrep, Ecm, min_wall=CPU_WORK_WINDOW, min_calls=20):
    """Process CPU time / wall time of warm `want=("Q",)` calls only, over at least `min_wall` s and `min_calls` calls."""
    _fast(irrep, Ecm, C_BASELINE, want=("Q",))   # warm-up (projector, threshold tables, F class tables)
    settled = _settle_idle()
    n = 0
    t0c, t0w = time.process_time(), time.perf_counter()
    while True:
        _fast(irrep, Ecm, C_BASELINE, want=("Q",))
        n += 1
        wall = time.perf_counter() - t0w
        if wall >= min_wall and n >= min_calls:
            break
    cpu = time.process_time() - t0c
    return cpu, wall, n, settled


def test_fast_single_thread_cpu_time():
    """Report 006 note 1: the NumPy tail used ~4 cores (CPU/wall 3.9 at N = 160).  With the tail in C++ the fast path
    must stay on one core without any environment pinning (no OPENBLAS_NUM_THREADS / OMP_NUM_THREADS in the test).

    Amendment acceptance-v5 (audit 001 ISSUE C): the earlier version measured a ~40 ms window at N = 35, which (i) fell
    inside OpenBLAS's post-import thread spin when the test ran alone (CPU/wall 3.3-4.0, a false failure) and (ii) could
    not detect real BLAS threading, which only sets in at N >~ 100.  Now: settle until the process is idle, then measure
    process CPU / wall over >= 0.6 s of warm calls at the two large points where the NumPy tail did thread.  The bar
    (1.3) is unchanged.  Verified to FAIL (ratio > 2) against a variant whose tail calls NumPy/BLAS matmul + slogdet."""
    for irrep, Ecm in CPU_POINTS:
        cpu, wall, n, settled = _cpu_wall_ratio(irrep, Ecm)
        assert wall >= CPU_WORK_WINDOW and n >= 20, (irrep, Ecm, cpu, wall, n)
        assert cpu / wall <= CPU_WALL_MAX, (irrep, Ecm, "cpu", cpu, "wall", wall, "calls", n, "ratio", cpu / wall,
                                            "settle_s", settled)


@pytest.mark.slow
def test_fast_speedup_single_thread_largest_point_ge_150x():
    from qc3 import fast, oracle
    sec = sector("111_A2")
    Ecm = 0.35103   # N = 160, the largest anchor
    opts = {"projector": "orbit"}
    fast.evaluate_point(sec, Ecm, C_BASELINE, want=("Q",), options=opts)   # warm-up (projector, threshold tables)
    t_f = []
    for _ in range(3):
        t0 = time.perf_counter(); fast.evaluate_point(sec, Ecm, C_BASELINE, want=("Q",), options=opts); t_f.append(time.perf_counter() - t0)
    t0 = time.perf_counter(); o = oracle.evaluate_point(sec, Ecm, C_BASELINE); t_o = time.perf_counter() - t0
    t_o_excl = t_o - float(o.timings.get("projector", 0.0))   # the oracle's own cold projector build is not the comparison
    ratio = t_o_excl / float(np.median(t_f))
    assert ratio >= SPEEDUP_MIN, (t_o, t_o_excl, np.median(t_f), ratio)

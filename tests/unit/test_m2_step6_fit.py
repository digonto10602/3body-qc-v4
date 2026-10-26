"""Unit tests for qc3.fit.fit_cutoff / verify_fit / save_results / the CLI helpers (M2 step 6, prompts/003 section 2.7).
Synthetic trackers whose roots are affine in c, so the chi^2 is exactly quadratic and the minimum, the Hesse errors and the
Gauss-Newton errors are known in closed form; no randomness, no C++ core except for the provenance block of the JSON."""
import json
from types import SimpleNamespace

import numpy as np
import pytest

from qc3 import fit
from qc3.point import LevelSet, MissingRootError

C0 = np.array([1.0, -2.0, 0.5, 0.25])
A = {"000_A1m": [0.280, 0.300, 0.330], "100_A2": [0.290, 0.310, 0.335]}
JJ = {"000_A1m": np.array([[-1.0e-6, 2e-7, 0, 1e-7], [-2e-6, 0, 3e-7, 0], [-3e-6, 1e-7, 1e-7, 5e-7]]),
      "100_A2": np.array([[-1.5e-6, 1e-7, 2e-7, 0], [-2.5e-6, 3e-7, 0, 2e-7], [-1e-6, -1e-7, 4e-7, 3e-7]])}


class Tracker:
    def __init__(self, irrep, n_keep=3, form="q0_consistent"):
        self.sector = SimpleNamespace(irrep=irrep, Lbyas=20, k2_pwave_form=form)
        self.a, self.J, self.n_keep = np.asarray(A[irrep]), JJ[irrep], n_keep
        self.c, self.E_max = C0.copy(), 0.34
        self.n_evaluations = self.n_local_rescans = self.n_full_rescans = 0

    def update(self, c):
        self.c = np.asarray(c, float).copy()
        self.n_evaluations += 3
        return self.energies()

    def energies(self):
        return (self.a + self.J @ (self.c - C0))[: self.n_keep]

    def jacobian(self):
        return self.J[: self.n_keep].copy()

    def wrong_sign_roots(self):
        return []


def make_levels():
    labels = [("000_A1m", 0), ("000_A1m", 1), ("000_A1m", 2), ("100_A2", 0), ("100_A2", 1), ("100_A2", 2)]
    E = np.array([0.2803, 0.3004, 0.3297, 0.2899, 0.3103, 0.3352])
    s = 1e-3 * np.ones(6)
    cov = np.outer(s, s) * (0.3 + 0.7 * np.eye(6))
    return LevelSet(labels=labels, irreps=["000_A1m", "100_A2"], nP=[(0, 0, 0)] * 6, files=[""] * 6, Ncfg=10, Ecm=E, cov=cov,
                    Ecm_cutoff=0.335, Lbyas=20)


def setup(**kw):
    return make_levels(), {ir: Tracker(ir, **kw) for ir in A}


def analytic(lv, free):
    """Weighted linear least squares of E_lat = a + J (c - C0) over the free parameters (others at C0)."""
    J = np.vstack([JJ["000_A1m"], JJ["100_A2"]])
    a = np.concatenate([A["000_A1m"], A["100_A2"]])
    Ci = np.linalg.inv(lv.cov)
    Jf = J[:, free]
    N = Jf.T @ Ci @ Jf
    dc = np.linalg.solve(N, Jf.T @ Ci @ (lv.Ecm - a))
    c = C0.copy(); c[free] += dc
    r = lv.Ecm - (a + J @ (c - C0))
    return c, float(r @ Ci @ r), np.sqrt(np.diag(np.linalg.inv(N)))


def test_four_parameter_fit_finds_the_exact_minimum_and_errors():
    lv, tr = setup()
    r = fit.fit_cutoff(lv, tr, C0 + np.array([2e3, -1e3, 5e2, 3e3]) * 1.0, verify=False)
    c_ref, chi2_ref, err_ref = analytic(lv, [0, 1, 2, 3])
    assert r.valid and r.edm < 2e-3 and r.n_free == 4 and r.ndof == 2 and r.fixed == ()
    assert np.allclose(r.c, c_ref, rtol=0, atol=1e-3 * err_ref.max()) and r.chi2 == pytest.approx(chi2_ref, abs=1e-6)
    assert np.allclose(r.errors, err_ref, rtol=1e-2) and np.allclose(r.errors_gauss_newton, err_ref, rtol=1e-8)   # linear model: GN == Hesse
    assert r.chi2 <= r.chi2_lm + 1e-6 and r.chi2 <= r.chi2_start and r.n_penalised_calls_final_stage == 0
    assert np.allclose(r.errors_fd_hessian, err_ref, rtol=1e-6)                       # FD of the analytic gradient is exact here
    assert r.chi2_per_dof * r.ndof == pytest.approx(r.chi2, abs=1e-12)
    assert np.allclose(np.diag(r.correlation), 1.0) and np.allclose(np.sqrt(np.diag(r.covariance)), r.errors)
    assert r.unconstrained == (len(r.unconstrained_reasons) > 0)
    assert r.n_wrong_sign_zeros == 0 and r.wrong_sign_zeros == [] and r.verification is None


def test_fixed_parameters_stay_exactly_and_have_zero_rows():
    lv, tr = setup()
    c_start = (1.0 + 2e3, -2.0 - 1e3, 0.0, 0.0)
    r = fit.fit_cutoff(lv, tr, c_start, fixed=("KB", "KE"), verify=False)
    assert r.n_free == 2 and r.ndof == 4 and tuple(r.fixed) == ("KB", "KE")
    assert r.c[2] == 0.0 and r.c[3] == 0.0                                   # fixed at their START values (0), not at C0
    assert np.all(r.errors[2:] == 0.0) and np.all(r.errors_gauss_newton[2:] == 0.0) and np.all(r.errors_fd_hessian[2:] == 0.0)
    assert np.all(r.covariance[2:, :] == 0.0) and np.all(r.covariance[:, 2:] == 0.0) and np.all(r.correlation[2:, :] == 0.0)
    assert np.allclose(np.diag(r.correlation)[:2], 1.0) and np.all(np.diag(r.correlation)[2:] == 0.0)
    # closed form of the 2-parameter weighted linear least squares: E_lat = a + J (c - C0), c2 = c3 = 0
    J = np.vstack([JJ["000_A1m"], JJ["100_A2"]]); a = np.concatenate([A["000_A1m"], A["100_A2"]])
    Ci = np.linalg.inv(lv.cov); Jf = J[:, :2]
    N = Jf.T @ Ci @ Jf
    target = lv.Ecm - a - J[:, 2:] @ (np.zeros(2) - C0[2:]) + Jf @ C0[:2]    # E_lat - a - J_fixed (c_fixed - C0_fixed) = Jf c_free - Jf C0_free
    sol = np.linalg.solve(N, Jf.T @ Ci @ target)
    assert np.allclose(r.c[:2], sol, rtol=0, atol=1e-3 * r.errors[:2].max())
    assert np.allclose(r.errors[:2], np.sqrt(np.diag(np.linalg.inv(N))), rtol=1e-2)


def test_unknown_fixed_name_and_all_fixed_raise():
    lv, tr = setup()
    with pytest.raises(ValueError, match="unknown parameter"):
        fit.fit_cutoff(lv, tr, tuple(C0), fixed=("KX",), verify=False)
    with pytest.raises(ValueError, match="all four"):
        fit.fit_cutoff(lv, tr, tuple(C0), fixed=fit.PARAM_NAMES, verify=False)
    with pytest.raises(ValueError, match="minimizer"):
        fit.fit_cutoff(lv, tr, tuple(C0), minimizer="migrad", verify=False)


def test_start_vector_with_a_missing_level_raises_loudly():
    lv, tr = setup(n_keep=2)                                                 # 2 roots per sector for 3 levels
    with pytest.raises(MissingRootError, match="fewer QC zeros"):
        fit.fit_cutoff(lv, tr, tuple(C0), verify=False)


def test_unconstrained_reasons_flag_exactly_the_documented_conditions():
    c = np.array([1.0, 10.0, -5.0, 2.0])
    err = np.array([0.5, 12.0, 1.0, 1.0])                                    # Kiso1 and KE: error >= |value|? 12 >= 10 yes; 1 < 2 no
    corr = np.eye(4); corr[0, 2] = corr[2, 0] = -0.995; corr[1, 3] = corr[3, 1] = 0.98
    reasons = fit._unconstrained([0, 1, 2, 3], c, err, corr, accurate=False)
    assert reasons == ["hesse_not_accurate", "error_exceeds_value:Kiso1", "correlation_ge_0.99:Kiso0,KB"]
    assert fit._unconstrained([0, 1], c, np.array([0.5, 1.0, 0, 0]), corr, accurate=True) == []
    # fixed parameters are never reported
    assert fit._unconstrained([0, 1], c, np.array([50.0, 1.0, 9e9, 9e9]), corr, accurate=True) == ["error_exceeds_value:Kiso0"]


def test_verify_fit_refuses_trackers_that_are_not_at_the_fit_parameters():
    lv, tr = setup()
    r = fit.fit_cutoff(lv, tr, tuple(C0), verify=False)
    tr["000_A1m"].c = r.c + 1.0
    with pytest.raises(ValueError, match="not at the fit"):
        fit.verify_fit(r, tr, lv)


def test_chi2_from_roots_uses_d3_and_the_penalty():
    lv, tr = setup()
    f = fit.QCChi2(lv, tr)
    E, _ = f.model(C0)
    roots = {ir: tr[ir].energies() for ir in A}
    chi2, E_m = f.chi2_from_roots(roots)
    assert np.array_equal(E_m, E) and chi2 == pytest.approx(f(C0), rel=1e-14)
    roots_short = {"000_A1m": roots["000_A1m"][:2], "100_A2": np.append(roots["100_A2"], 0.36)}     # n<m and n>m
    chi2_s, E_s = f.chi2_from_roots(roots_short)
    assert E_s[2] == pytest.approx(lv.Ecm[2] - 100.0 * f.sigma[2]) and np.array_equal(E_s[3:], E[3:])
    assert chi2_s > chi2


def test_save_results_writes_provenance_constants_and_roundtrips(tmp_path):
    lv, tr = setup()
    r = fit.fit_cutoff(lv, tr, tuple(C0), fixed=("KE",), verify=False)
    r.verification = {"counts_equal": True, "max_root_diff": float("inf"), "chi2_fresh": r.chi2}
    path = tmp_path / "sub" / "fit.json"
    fit.save_results([r], path, extra={"note": {"x": np.float64(1.5)}})
    d = json.loads(path.read_text())
    assert d["kind"] == "qc3.fit" and d["is_a_fit"] is True and d["seeds"]["stochastic_inputs"] is None
    for k in ("git_commit", "git_tree_dirty", "core_built_at_commit", "core_version", "numpy", "scipy", "iminuit", "date", "machine",
              "command", "threads_env"):
        assert k in d["provenance"], k
    assert d["constants"]["production"]["xi"] == 3.444 and d["constants"]["k2_pwave_form"] == "q0_consistent"
    assert d["constants"]["covariance"].startswith("D1") and d["constants"]["matching"].startswith("D3")
    assert d["constants"]["exclusion"].startswith("D8") and d["constants"]["residue_sign"].startswith("D9")
    f0 = d["fits"][0]
    assert f0["chi2"] == r.chi2 and f0["fixed"] == ["KE"] and f0["n_free"] == 3 and f0["c"] == r.c.tolist()
    assert f0["verification"]["max_root_diff"] == "inf" and d["note"] == {"x": 1.5}


def test_cli_start_vector_presets_and_fixed_defaults_to_zero():
    assert fit._start_vector(["baseline"], ()) == fit.C_BASELINE
    c = fit._start_vector(["baseline"], ("KB", "KE"))
    assert c == (fit.C_BASELINE[0], fit.C_BASELINE[1], 0.0, 0.0)
    assert fit._start_vector(["v32w"], ()) == fit.C_V32W and fit._start_vector(["legacy0335"], ()) == fit.C_LEGACY_0335
    assert fit._start_vector(["1", "2", "3", "4"], ("KB",)) == (1.0, 2.0, 3.0, 4.0)
    with pytest.raises(SystemExit):
        fit._start_vector(["nonsense"], ())
    args = fit.build_parser().parse_args(["--fix", "KB", "KE", "--cutoffs", "0.335"])
    assert args.fix == ["KB", "KE"] and args.cutoffs == [0.335] and args.form == "q0_consistent"


class FlakyTracker(Tracker):
    """Raises AmbiguousRootError once, on its first move away from the start vector (an LM trial point on a zero next to a pole)."""

    def __init__(self, *a, **kw):
        super().__init__(*a, **kw)
        self.fail_once = True

    def update(self, c):
        from qc3.roots import AmbiguousRootError
        if self.fail_once and np.any(np.asarray(c, float) != self.c):
            self.fail_once = False
            raise AmbiguousRootError("synthetic ambiguous candidate at a trial point")
        return super().update(c)


def test_lm_trial_point_with_an_ambiguous_candidate_is_rejected_counted_and_reported():
    from qc3.roots import AmbiguousRootError
    lv = make_levels()
    start = tuple(C0)
    tr = {"000_A1m": FlakyTracker("000_A1m"), "100_A2": Tracker("100_A2")}
    r = fit.fit_cutoff(lv, tr, start, verify=False)
    c_ref, chi2_ref, _ = analytic(lv, [0, 1, 2, 3])
    assert len(r.ambiguous_trials_lm) == 1 and "synthetic" in r.ambiguous_trials_lm[0]["message"]
    assert r.valid and r.chi2 == pytest.approx(chi2_ref, abs=1e-6) and r.n_penalised_calls_final_stage == 0
    assert "WARNING: LM trial point rejected" in fit.format_fit(r)
    tr2 = {"000_A1m": FlakyTracker("000_A1m"), "100_A2": Tracker("100_A2")}
    with pytest.raises(AmbiguousRootError):                                    # the switch restores the loud behaviour
        fit.fit_cutoff(lv, tr2, start, verify=False, reject_ambiguous_trials=False)


def test_error_gap_flag_synthetic():
    err = np.array([1.0, 2.0, 0.0, 4.0])
    ok = np.array([1.1, 1.8, 0.0, 4.4])                       # ratios 1.1, 0.9, (fixed), 1.1: inside [1/1.2, 1.2]
    gap, dis, names = fit.error_gap_flags(err, ok, (0, 1, 3))
    assert np.allclose(gap[[0, 1, 3]], [1.1, 0.9, 1.1]) and np.isnan(gap[2]) and not dis and names == []
    bad = np.array([1.3, 2.0, 0.0, 3.0])                      # 1.3 > 1.2 ; 0.75 < 1/1.2 ; middle one fine
    gap, dis, names = fit.error_gap_flags(err, bad, (0, 1, 3))
    assert dis and names == [fit.PARAM_NAMES[0], fit.PARAM_NAMES[3]]
    gap, dis, names = fit.error_gap_flags(err, np.full(4, np.nan), (0, 1, 3))      # FD Hessian failed: cannot confirm -> flagged
    assert dis and len(names) == 3
    gap, dis, names = fit.error_gap_flags(err, bad, (1,))      # only the free parameter counts
    assert not dis and names == []

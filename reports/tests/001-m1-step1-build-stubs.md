# Test report 001: M1 step 1 "Build + stubs"

- **Date:** 2026-10-21
- **Git commit:** `5239739` for the checks; HEAD moved to `cb90300` during the run (that commit touches only `CLAUDE.md` and `graphify-out/`, no code)
- **Prompt in force:** `prompts/001-m1-cache-free-single-energy-qc.md` (§3 row 1, §5)
- **Mode:** step
- **Related:** `reports/research/007-single-energy-qc-methods-and-reference.md`
- **Raw numbers:** `results/test_m1_step1_build_stubs.json`
- **Machine:** Linux 6.18, g++ 13.3.0, Python 3.13.16, NumPy 2.5.3, SciPy 1.18.1, pybind11 3.1.0. No stochastic tests, so no seeds.

## Commands run (repo root, `PYTHONPATH=src/python`)
1. `python3 -c "import qc3, qc3.oracle, qc3.fast, qc3.lattice, qc3.qc, qc3.roots, qc3.bench, qc3.point, qc3.sector, qc3.constants"`
2. `rm -f src/python/qc3/_core*.so; ./build.sh core; python3 -c "import qc3._core as c; print(c.__file__, c.version(), c.git_commit())"`
3. `python3 -m pytest --collect-only -q tests/acceptance -p no:cacheprovider`
4. `python3 -m pytest -q tests/acceptance/m1 -p no:cacheprovider --tb=line --junitxml=...` (failure types parsed from the XML)
5. Programmatic comparison of `qc3.constants` and `Sector` defaults to plan §5.1/§5.2 and invariant 3; grep of imports in all stubs.
6. `./build.sh all` (5 m 32 s), `ls bin`, `git check-ignore -v`.
7. `.claude/bin/acceptance-check`
8. `python3 -m pytest -q tests --ignore=tests/acceptance -p no:cacheprovider`

## Results
| # | check | expected | got | tolerance | verdict |
|---|---|---|---|---|---|
| 1 | import of 10 qc3 modules | works | `import OK`, rc 0 | exact | PASS |
| 2 | `./build.sh core` from clean (`.so` deleted first) | rc 0 | rc 0, 7.4 s | exact | PASS |
| 2 | `qc3._core` import | file, version, commit | `.../src/python/qc3/_core.cpython-313-x86_64-linux-gnu.so 0.1.0-step1 5239739` (commit equals HEAD at build time) | exact | PASS |
| 3 | collect acceptance | 33 collected, 0 errors | 33 collected, 0 errors (0.28 s) | exact | PASS |
| 4 | `pytest tests/acceptance/m1` | 33 fail, all NotImplementedError | 33 failed, 0 passed; types: NotImplementedError 33, any other 0 | exact | PASS |
| 5a | `PRODUCTION`, `C_BASELINE`, `IRREPS`, `NNP`, `GRID_*` vs §5.1 / invariant 3 | identical | `PRODUCTION == expected dict` True; `C_BASELINE` True; `GRID_DE == 0.0969/19999` True; `GRID_E0=0.2631`, `GRID_N=20000`; IRREPS and NNP equal §5.1 | exact | PASS |
| 5b | `Sector` fields/defaults vs §5.2 and `PRODUCTION` | identical | all 12 numeric defaults equal `PRODUCTION`; `a0_Kpi, a1_Kpi, a0_KK = 4.04, -43.2, 4.12`; `Lbyas=20`; `k2_pwave_form="coded_v2"`; `frozen=True` | exact | PASS |
| 5c | no stub imports `_core` | none | grep: only the docstring in `__init__.py` mentions it; no import | exact | PASS |
| 5d | no physics in stubs | only constants/dataclasses/dtype | `oracle`, `fast`, `qc`, `roots`, `lattice`, `bench`, `Sector` methods all only `raise NotImplementedError`; `point.py` has dtype, dataclasses, `MissingRootError`, `ALL` | n/a | PASS |
| 6 | `./build.sh all` | rc 0, fitter skipped with a message | rc 0; message `[build] Minuit2 not found: skipping fitter`; `bin/` has 4 executables (hotwindow, runtime_k3basis_v33g, det_scan_v33h, generate_cpu_component_cache); no fitter binary; `_core` rebuilt | exact | PASS |
| 6 | gitignore | `bin/`, `*.so` ignored | `.gitignore:1:/bin/`, `.gitignore:4:*.so` | exact | PASS |
| 7 | `acceptance-check` | PASS | `PASS: tests/acceptance matches acceptance-v1` | exact | PASS |
| 8 | fast check (no acceptance) | n/a | `no tests ran`, pytest rc 5 (no tests collected); `tests/unit` empty, `tests/reference` has scripts only | n/a | NOTE |

## Failure breakdown of check 4 (33 tests, all expected at this step)
Stub that raised: `Sector.production` 19, `qc.project_Q` 4, `lattice.load_levels` 4, `lattice.default_data_dir` 2, `lattice.read_jack` 1, `lattice.jackknife_resample` 1, `lattice.match_levels` 1, `lattice.chi2` 1.

- `test_m1_portability.py::test_no_absolute_home_paths_in_m1_code` fails with `NotImplementedError` (`qc3.lattice.default_data_dir: prompts/001 step 5`). This is the known, accepted plan slip. Its sibling `test_data_root_env_var` also fails, as expected.
- Pending (later-step NotImplementedError): all 33. None is a regression, since no milestone is complete.

## Observations (not failures)
- The stub messages read `qc3.qc.project_Q: prompts/001 step 4`, not exactly the template `"... not implemented (prompts/001 step N)"` in plan §3. They name the step, which is all that matters here.
- `qc3.bench.main()` ignores the CLI arguments (`--out`, `--no-oracle`, `--threads`) and just calls `run()`. It is a stub; step 8 must implement the arguments.
- The step-1 row says "`pytest tests/acceptance/m1 -x -q` fails only with NotImplementedError". Confirmed (33/33).
- Fast check has no unit tests yet, as agreed, so check 8 has no content. `affected` coverage: nothing beyond these stubs exists yet.
- `./build.sh all` takes 5 m 32 s, dominated by the four C++ tools. `core` alone takes 7.4 s.

## Verdict: PASS
33/33 expected failures, 0 unexpected, 0 build errors, `acceptance-check` PASS.

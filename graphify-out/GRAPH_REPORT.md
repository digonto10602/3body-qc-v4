# Graph Report - 3body-qc-v4  (2026-10-22)

## Corpus Check
- 292 files · ~460,637 words
- Verdict: corpus is large enough that graph structure adds value.

## Summary
- 6236 nodes · 11389 edges · 356 communities (294 shown, 60 thin omitted)
- Extraction: 96% EXTRACTED · 4% INFERRED · 0% AMBIGUOUS · INFERRED: 438 edges (avg confidence: 0.85)
- Token cost: 0 input · 0 output

## Graph Freshness
- Built from commit: `8f8e9d32`
- Run `git rev-parse HEAD` and compare to check if the graph is stale.
- Run `graphify update .` after code changes (no API cost).

## Community Hubs (Navigation)
- Code and function index
- Package `corrected_reader`
- Package `production_fitter`
- Package `cache_generator`
- f
- projections_v1.hpp
- generate_swave_cpu_cache_qc.cpp
- G_gpu_safe_builder.cuh
- Options
- generate_cpu_component_cache.cpp
- FitSettings
- FitSettings
- run_v33n_cutoff_scan_fits.py
- nonint3body_general_doublecover.py
- EigZeroRecordV30s
- lattice_data_covariance_cpp.hpp
- fv_projector_cartesian_l1_v30m.hpp
- KKπ Three-Body QC Code Version History
- MultiLFCN
- MultiConfig
- qc_fitter_norm_refine_v2_multiL.cpp
- compare_gpu_cache_vs_cpu_openmp_det_oldscale.cpp
- C. Fitter lineage and fast-fit caching design (K K-bar pi, RFT three-body QC)
- RowOut
- V31aEigPoint
- The physics loop (active when `notes/STATUS.md` exists in the project)
- DeviceBuffer
- v32w_single_L_fitter/K3df_minuit_fit_v32f_fullF3inv_QCfull_cached_classifier.hpp
- run
- projections_gpu_safe.cuh
- v33g_build_runtime_k3basis_cache.cpp
- run_block
- v33h_patched_gpu_cache_oldscale_det_scan.cpp
- EvalData
- plot_classifier_algorithm_sweep_v33l.py
- F2_functions_v2.h
- CpuGpuChunk
- functions.h
- AGENTS.md
- PipelineGpuCache
- F2_gpu_safe_builder.cuh
- string
- RuntimeCacheMeta
- Test report 008: M1 step 7.0 (projector cache key), 7.1 (pole/zero discrimination), 7.8 (D7 default q0_consistent)
- ZeroRecord
- F3_gpu_omp_cublas_pipeline_v9_inward_shape_refinement.cu
- Full F3 GPU pipeline versions
- v33g_runtime_k3basis_cache.hpp
- KKpi Codebase Map for Codex
- E — `threebody_fullmatrix_hpc_cuda_bench_v3_codex`: what it is and whether it belongs in 3body-qc-v4
- K3_functions_2plus1.hpp
- MomentumIrrepSpec
- fitter/K3df_minuit_fit_v32f_fullF3inv_QCfull_cached_classifier.hpp
- DeviceResidentResult
- F3_gpu_cachegen_v32zu.cu
- PhysicsParams
- vector
- PhysicsParams
- Report B — F3inv/Vsel cache generation (CPU + GPU) and root-finding lineage for the projected QC determinant
- K2_functions_gpu_safe.cuh
- SvdMinRecordV30r
- vector
- v32zt_generate_one_irrep
- Version history summary
- KKpi v33p Cache-Anchored Classifier/Fitter Checkpoint
- slogdet
- run_groups_concurrent
- QCPointValue
- FitResult
- ProjectedQCCacheEntry
- inspect_v32_f3inv_cache_v32f_fixed.py
- load_multil_jack_targets
- vector
- make_qc_state_lines_from_scattered_spectrum.py
- BlockSummary
- Cand
- parse_label
- `corrected_reader`
- F3_device_resident_true_batch_v21_side_abs_trend_refine.cu
- ProjectedQCCacheEntry
- F3_device_resident_true_batch_v16_custom_batched_getrs.cu
- F3_cpu_openmp_v25_K3QC_cached_core.hpp
- QCRefinedFCN
- gpu_cache_reader_convention.hpp
- EnergyGpuResult
- QCRefinedFCN
- bench.py
- V32ztAllObjRowPayload
- CachedProjectorV30q
- Tester report 004: M1 step 4, QC helpers + lattice reader (D1, D3)
- acceptance-check
- FinalBracket
- Cand
- L20_spectrum_plotter_model_true_zeros_yaxis_fixed.py
- QC_functions_v2.h
- fitter/qc_fitter_norm_refine_v2.cpp
- AcceptedZeroWindow
- inspect_v32_f3inv_cache_mmap_progress.py
- 3. Results
- F3MatrixSaveOptions
- SignFlipCandidate
- FitResult
- PhysicsParams
- DeviceConfig
- find_QC_zeros_refined
- Codex Report: Full Cache Generator Package
- ClassifyResult
- InwardShapeResult
- BuiltRow
- ScanArgs
- ClassifyResult
- roots.py
- METHODS — numerical methods, equations, conventions, sources
- 3. Findings
- debug_v30r_svd_validated_cartesian_l1_100_A2_zero_scan.cpp
- Faddeeva.cc
- v33d CPU/GPU Component Cache Generator Package
- acceptance-lock
- EvalPoint
- v32w_single_L_fitter/qc_fitter_norm_refine_v2.cpp
- ConfigMaker4Params
- run_v33m_classifier_speed_search.py
- RowResult
- process_same_shape_group_true_batch_v15
- ofstream
- next-num
- F3iProjectedCacheEntry
- GroupStream
- Full L20 fit and validation report
- KKpi v33o Production Fitter
- F3_cpu_openmp_v24_K3_QC_multi_signflip_side_abs_refinement.cpp
- EnergyPrepassInfo
- PhysicsParams
- generate_K3iso.cpp
- 5. Key physics and numerics lessons (Question 2)
- CombinedRawGpuCacheRowHeader
- V32ztCacheRowPayload
- V32ztTiming
- map
- ThreadSafeQueue
- OPEN ISSUES — read before trusting any number
- .upload
- K3iso_maker.py
- scan_projF3inv_digonto_classifier_v2.cpp
- IrrepResult
- TargetLevel
- string
- V32zuRunOptions
- process_same_shape_group_true_batch_v15
- 002 — M1 step 7 re-planned: cold projector from orbits, threshold basis builder, pair kernels for $\mathcal K_{\rm df,3}$ and $G$, C++ tail, and pole/zero discrimination in the root refiner
- EvalResult
- ZeroPole
- functions_gpu_config_maker4.cuh
- gpu_config_maker_4_many_entries
- 001 — M1: a cache-free, CPU-only, single-energy three-body QC with an oracle and an independent reference
- EvalResult
- QCPointValue
- ZeroPole
- PIPELINE — from physics parameters to a $\mathcal K_{\rm df,3}$ fit
- KKpi L20 user-validated hot-window fitter
- KKpi v33o production fitter
- ScanPoint
- make_refined_mesh_v32s
- JumpRow
- Input
- Tester report 003: M1 step 3, reference checks and the `q0_consistent` p-wave form
- AssignmentRow
- KKpi v32w + v32y reproducibility package
- Key numbers
- InwardShapeRefineOptions
- 3.3 Question 2 — equations implemented
- L20 full cutoff-dependence final report
- run_v33p_classifier_sweep_L20_100_A2.py
- ScanPoint
- SideClass
- load_binary_f3inv_cache_v32f
- save_binary_f3inv_cache_v32f
- find_QC_zeros_refined
- Tester report 002: M1 step 2, oracle wrapper
- L20 K3B/K3E zero tests and chi-square audit
- comp
- SignFlip
- save_binary_f3inv_cache_v32f
- DetInfo
- config_gpu_component_cachegen_sep_f2k2.sh
- Final shutdown handoff: L20 user-validated hot-window fitter
- Classifier Restart Notes
- Function Map
- Core function inventory with descriptions
- DetInfo
- __device__
- Output
- STATUS
- main
- Kin
- F_channel
- ComponentMatrixCacheHeader
- EnergyPrepassInfo
- K2inv_EREord2_2plus1_mat
- ClassifierDispatchInfo
- DeviceResidentResult
- build.sh
- Version summary and working status
- GPU Cache Reader Fix
- DevicePtr
- F3ScopedTimer
- F3ScanRunFiles
- WorkBracket
- L20 user-validated pass-only fit final report
- audit_l20_chi2_reconstruction.py
- DevicePtr
- fast.py
- SingleThread
- WindowGridRow
- compare_K3_cpp_vs_python.py
- CLAUDE.md — project rules for 3body-qc-v4
- LocalScopedTimer
- RefineTask
- EnergyShapeKey
- Independent L20 chi-square reconstruction
- L20 pass-only completed checkpoint
- GitHub checkpoint push failure
- L20_LATTICE_VS_MODEL_ALL_CUTOFFS.md
- L20_PARAMETER_COVARIANCE_COMPARISON.md
- main
- build_l20_user_validated_windows.py
- ConfigEntry
- ScopedTimer
- EnergyShapeKey
- V32zuF2K2ValidationStats
- 000 — Consolidation and roadmap (starter plan from chat, 2026-10-21)
- RefineTask
- Tools quick guide
- Hot-window FCN checkpoint
- Where the scatter parameters are used
- GPU determinant backend audit
- run_v32w_QC_fitter_norm_refine_v2_fixedcoarse.sh
- run_gpu_component_cachegen.sh
- plot_v33p_classifier_det_family_interactive.py
- sector
- module.cpp
- Options
- ROOT FINDING for the three-body QC
- v33p_gpu_det_lu.cu
- fgk.py
- Test report 001: M1 step 1 "Build + stubs"
- load_labels
- make_spectrum_config_from_fit.py
- main
- materialize_v33p_hot_windows.py
- main
- SvdMinCandidateV30r
- FCNBase
- vector
- QCSearchTiming
- L20 full Ecm-cutoff fit suite
- L20 user-validated prefit coverage audit
- plot_fit_levels.py
- load_runtime_cache
- Options
- test_qc_lattice_roots.py
- run_nonint3body_general_doublecover.sh
- build_gpu.sh
- .val500_100A2_gpu_component_cachegen.sh
- L20_000_A1m/PROJECTED_BASIS_HOTWINDOW_CACHE_REPORT.md
- L20_100_A2/PROJECTED_BASIS_HOTWINDOW_CACHE_REPORT.md
- L20_110_A2/PROJECTED_BASIS_HOTWINDOW_CACHE_REPORT.md
- L20_111_A2/PROJECTED_BASIS_HOTWINDOW_CACHE_REPORT.md
- L20_200_A2/PROJECTED_BASIS_HOTWINDOW_CACHE_REPORT.md
- data/README.md
- L20_L24_NEXT_STEPS.md
- L20_LARGE_FILES_AND_CACHE_POLICY.md
- v31_v32_ARCHIVE_VERSION_GUIDE.md
- v33p_GPU_VS_CPU_DET_VALIDATION_STATUS.md
- PROVENANCE.md
- legacy/README.md
- CHI2_COMPARISON_TO_OLD_ISOTROPIC.md
- PREFIT_AUDIT_ZERO_TESTS.md
- ZERO_TEST_FIT_SUMMARY_COMPARISON.md
- ZERO_TEST_ROOT_SHIFT_COMPARISON.md
- Ecm_cut_0p315/FIT_SUMMARY.md
- Ecm_cut_0p315/LATTICE_INPUT_PROVENANCE.md
- Ecm_cut_0p315/lattice_vs_model_levels.md
- Ecm_cut_0p315/PARAMETER_IMPACT_SUMMARY.md
- Ecm_cut_0p325/FIT_SUMMARY.md
- Ecm_cut_0p325/LATTICE_INPUT_PROVENANCE.md
- Ecm_cut_0p325/lattice_vs_model_levels.md
- Ecm_cut_0p325/PARAMETER_IMPACT_SUMMARY.md
- Ecm_cut_0p335/FIT_SUMMARY.md
- Ecm_cut_0p335/LATTICE_INPUT_PROVENANCE.md
- Ecm_cut_0p335/lattice_vs_model_levels.md
- Ecm_cut_0p335/PARAMETER_IMPACT_SUMMARY.md
- Ecm_cut_0p345/FIT_SUMMARY.md
- Ecm_cut_0p345/LATTICE_INPUT_PROVENANCE.md
- Ecm_cut_0p345/lattice_vs_model_levels.md
- Ecm_cut_0p345/PARAMETER_IMPACT_SUMMARY.md
- L20_K3DF_PARAMETER_DRIFT.md
- L20_PASSONLY_FIT_SUMMARY.md
- L20_PASSONLY_LATTICE_VS_MODEL_ALL_CUTOFFS.md
- L20_PASSONLY_PARAMETER_DRIFT.md
- L20_PASSONLY_TRUE_ZEROS_ALL_WINDOWS_ALL_CUTOFFS.md
- PREFIT_COVERAGE_AUDIT.md
- README_v32w_FIXED_COARSE.md
- results/README.md
- compile_cpu_component_cachegen.sh
- run_cpu_component_cachegen.sh
- compile_gpu_component_cachegen.sh
- mutex
- DeviceConfig
- src/README.md
- real_wigner_d.hpp
- test_m1_fast_v2.py
- NeumaierComplexSum
- nonint3body_general_doublecover/README.md
- F — Development history v31 → v32y: lessons, bugs fixed, dead ends
- anchor_points
- dig_tools.hpp
- 4. Q2 - Exactly how the fast FCN works (L20 code)
- RepDiagnostics
- ZeroCandidate
- Sector
- /kickoff — start the physics loop here
- InwardShapeResult
- test_m1_reference.py
- Options
- k2_pwave_shift.py
- fixer.md
- tester.md
- guard.py
- Tester report 005: M1 step 5, portability + three no-silent-failure fixes
- string
- Research 008 — M1 step 7 re-plan: cold projector, basis builder, $\mathcal K_3$/$G$ kernels, tail, and pole/zero discrimination
- test_m1_lattice.py
- auditor.md
- builder.md
- EnergyMatrixPack
- unit/conftest.py
- session_start.py
- BasisRow
- chi2_covariance_convention_check.py
- test_K3mat_2plus1_cpp_dump.cpp
- Report 007: M1 step 9, p-wave $\mathcal K_2^{-1}$ level-shift table (D7)
- DeviceConfigView
- Tester report 006: M1 step 6, fast path v1 (rungs 1-4)
- default_data_dir
- 7. Q4 - Fit results achieved
- 9. Q6 - Recommended minimal fitter file set for the consolidated repo
- conftest.py
- BenchmarkTiming
- string
- DetRow
- load_binary_f3inv_cache_v32f
- F2GpuOptions
- check_cuda

## God Nodes (most connected - your core abstractions)
1. `Code and function index` - 225 edges
2. `Package `corrected_reader`` - 170 edges
3. `Package `production_fitter`` - 169 edges
4. `Package `cache_generator`` - 158 edges
5. `Options` - 94 edges
6. `FitSettings` - 75 edges
7. `FitSettings` - 75 edges
8. `MultiConfig` - 75 edges
9. `map` - 72 edges
10. `EvalData` - 55 edges

## Surprising Connections (you probably didn't know these)
- `test_find_sign_changes_and_zeros_synthetic()` --uses--> `MissingRootError`  [INFERRED]
  tests/unit/test_qc_lattice_roots.py → src/python/qc3/point.py
- `test_covariance_symmetric_psd_and_correlated()` --calls--> `load_levels()`  [INFERRED]
  tests/acceptance/m1/test_m1_lattice.py → src/python/qc3/lattice.py
- `test_real_file_000_A1m_n0_vs_independent_numpy()` --calls--> `load_levels()`  [INFERRED]
  tests/acceptance/m1/test_m1_lattice.py → src/python/qc3/lattice.py
- `find_root()` --uses--> `MissingRootError`  [INFERRED]
  scripts/m1/k2_pwave_shift.py → src/python/qc3/point.py
- `test_match_levels_edge_cases()` --uses--> `MissingRootError`  [INFERRED]
  tests/unit/test_qc_lattice_roots.py → src/python/qc3/point.py

## Import Cycles
- None detected.

## Communities (356 total, 60 thin omitted)

### Community 0 - "Code and function index"
Cohesion: 0.01
Nodes (225): Code and function index, `L20_spectrum_plotter_model_true_zeros_yaxis_fixed.py`, `packaging/generate_function_index.py`, `scripts/assemble_l20_passonly_reports.py`, `scripts/assemble_l20_zero_test_reports.py`, `scripts/audit_l20_chi2_reconstruction.py`, `scripts/audit_l20_user_validated_prefit.py`, `scripts/build_l20_user_validated_windows.py` (+217 more)

### Community 1 - "Package `corrected_reader`"
Cohesion: 0.01
Nodes (170): Package `corrected_reader`, `scripts/check_v33f_truezero_counts.py`, `scripts/make_spectrum_config_from_fit.py`, `scripts/make_v33e_cache_file_list.py`, `scripts/plot_classifier_algorithm_sweep_v33i.py`, `scripts/plot_classifier_algorithm_sweep_v33j.py`, `scripts/plot_classifier_algorithm_sweep_v33l.py`, `scripts/plot_F3inv_iso_nscale_with_v3_zeros.py` (+162 more)

### Community 2 - "Package `production_fitter`"
Cohesion: 0.01
Nodes (169): Package `production_fitter`, `scripts/check_v33f_truezero_counts.py`, `scripts/make_spectrum_config_from_fit.py`, `scripts/make_v33e_cache_file_list.py`, `scripts/plot_classifier_algorithm_sweep_v33i.py`, `scripts/plot_classifier_algorithm_sweep_v33j.py`, `scripts/plot_classifier_algorithm_sweep_v33l.py`, `scripts/plot_F3inv_iso_nscale_with_v3_zeros.py` (+161 more)

### Community 3 - "Package `cache_generator`"
Cohesion: 0.01
Nodes (158): `cpu_component_cachegen/source/generate_cpu_component_cache.cpp`, `gpu_component_cachegen/scripts/run_gpu_component_cachegen.sh`, `gpu_component_cachegen/source/F2_gpu_safe_builder.cuh`, `gpu_component_cachegen/source/F3_device_resident_true_batch_v16_custom_batched_getrs.cu`, `gpu_component_cachegen/source/F3_gpu_cachegen_v32zu.cu`, `gpu_component_cachegen/source/functions_gpu_config_maker4.cuh`, `gpu_component_cachegen/source/functions.h`, `gpu_component_cachegen/source/G_gpu_safe_builder.cuh` (+150 more)

### Community 4 - "f"
Cohesion: 0.06
Nodes (67): lattice_counts(), main(), covariance_status(), d(), main(), Path, read_csv(), write_csv() (+59 more)

### Community 5 - "projections_v1.hpp"
Cohesion: 0.13
Nodes (44): Derived, EigenBase, Matrix, Matrix4d, MatrixBase, all_perms_of_3(), analyze_projector_eigensystem(), blockDiag() (+36 more)

### Community 6 - "generate_swave_cpu_cache_qc.cpp"
Cohesion: 0.05
Nodes (71): compute_diagnostics(), compute_F3iso_inv_from_F3inv_full(), Config, cache_file, coarseN, combined_output_file, const_norm_power, Ecm_max (+63 more)

### Community 7 - "G_gpu_safe_builder.cuh"
Cohesion: 0.07
Nodes (72): boost_gpu(), cabs(), cabs2(), cadd(), cdiv(), cdiv_d(), cexp_c(), cmul() (+64 more)

### Community 8 - "Options"
Cohesion: 0.03
Nodes (71): Options, accept_only_likely_zero, alpha, atmK, atmpi, coarse_parallel, debug, dimension_jump_guard (+63 more)

### Community 9 - "generate_cpu_component_cache.cpp"
Cohesion: 0.06
Nodes (65): build_component_entry(), ComponentEntry, dim1, dim2, Ecm, En, error, F2 (+57 more)

### Community 10 - "FitSettings"
Cohesion: 0.03
Nodes (61): FitSettings, alpha, atmK, atmpi, binary_f3inv_cache_file, build_cache_once, chi_square_mode, classifier_peak_ratio_min (+53 more)

### Community 11 - "FitSettings"
Cohesion: 0.03
Nodes (61): FitSettings, alpha, atmK, atmpi, binary_f3inv_cache_file, build_cache_once, chi_square_mode, classifier_peak_ratio_min (+53 more)

### Community 12 - "run_v33n_cutoff_scan_fits.py"
Cohesion: 0.11
Nodes (57): Axes, all_candidates_file(), append_report_section(), base_config_text(), BenchSummary, bestfit_qc_file(), build_block_index(), build_block_index_by_file_irrep() (+49 more)

### Community 13 - "nonint3body_general_doublecover.py"
Cohesion: 0.08
Nodes (48): build_double_group(), class_algebra_irreps(), class_C2(), class_C2v(), class_C3v(), class_C4v(), class_Oh(), Config (+40 more)

### Community 14 - "EigZeroRecordV30s"
Cohesion: 0.05
Nodes (56): bisect_eigen_min_candidate_v30s(), closest_zero_eigenpair_from_matrix_v30s(), collect_eigen_min_signflips_v30s(), MatrixXcd, PhysicsParams, string, vector, VectorXcd (+48 more)

### Community 15 - "lattice_data_covariance_cpp.hpp"
Cohesion: 0.06
Nodes (52): CovarianceTuple, tuple, MomentumList, convert_jackknife_energy_to_ecm(), covariance_between_states_szscl21_based(), covariance_between_states_szscl21_based_result(), CovarianceResult, correlation_mat (+44 more)

### Community 16 - "fv_projector_cartesian_l1_v30m.hpp"
Cohesion: 0.15
Nodes (38): all_conventions(), best_projector(), BestProjectorResult, conv, diag, P, cart_axis_from_m_v30m(), cartesian_signed_perm_element_v30m() (+30 more)

### Community 17 - "KKπ Three-Body QC Code Version History"
Cohesion: 0.04
Nodes (47): C++ K3mat_2plus1 Header, C++ Lattice Covariance Header, Cache key, Cached object, Current Recommended Workflow, Features, First small test, `fit_levels_v25.dat` (+39 more)

### Community 18 - "MultiLFCN"
Cohesion: 0.07
Nodes (24): array, atomic, Minuit2::FCNBase, ofstream, time_point, MultiLFCN, accepted_windows, audit_normalized (+16 more)

### Community 19 - "MultiConfig"
Cohesion: 0.04
Nodes (48): MultiConfig, accepted_zeros_files, auto_build_missing_coarse, base, benchmark_repeat, benchmark_warmup, build_v33g_runtime_if_missing, cache_backend (+40 more)

### Community 20 - "qc_fitter_norm_refine_v2_multiL.cpp"
Cohesion: 0.18
Nodes (28): CacheKey, L, label, coarse_path_for(), contains_label(), convert_samples_to_ecm(), string, default_gpu_cache_path() (+20 more)

### Community 21 - "compare_gpu_cache_vs_cpu_openmp_det_oldscale.cpp"
Cohesion: 0.08
Nodes (51): optional, Args, gpu_cache_dir, irreps, k3, k3df_zero, Lbyas, outdir (+43 more)

### Community 22 - "C. Fitter lineage and fast-fit caching design (K K-bar pi, RFT three-body QC)"
Cohesion: 0.09
Nodes (22): 0. Path abbreviations used in this report, 10. What the results mean (summary), 11. Problems and assumptions, 1. What was asked, 2. Concrete actions taken, 3.1 Hash/size table of the multi-L fitter `qc_fitter_norm_refine_v2_multiL.cpp`, 3.2 What changed V33P (1373) $\to$ MF (2021): the hot-window FCN (+648 / -~60 lines), 3.3 What changed MF (2021) $\to$ L20 (2107): chi-square audit trace only (88 diff lines) (+14 more)

### Community 23 - "RowOut"
Cohesion: 0.05
Nodes (44): RowOut, const_norm_scale, critical_score, det, det_abs, det_imag, det_imag_abs, det_imag_rel (+36 more)

### Community 24 - "V31aEigPoint"
Cohesion: 0.06
Nodes (41): collect_v31a_flips(), MatrixXcd, string, vector, eigen_summary_from_projF3inv_v31a(), hermitize_v31a(), v31a_point_from_eval(), V31aEigPoint (+33 more)

### Community 25 - "The physics loop (active when `notes/STATUS.md` exists in the project)"
Cohesion: 0.13
Nodes (14): Build loop, Code understanding: graphify first, Coordinator rules (main session), Escalation — objective triggers, counted in STATUS, Global rules (Digonto) — loaded in every Claude Code session, graphify, Hard rules for correctness, How to report to me (every session, every project) (+6 more)

### Community 26 - "DeviceBuffer"
Cohesion: 0.10
Nodes (26): compare_F2_cpu_gpu(), complex, MatrixXcd, size_t, T, vector, DeviceBuffer, count (+18 more)

### Community 27 - "v32w_single_L_fitter/K3df_minuit_fit_v32f_fullF3inv_QCfull_cached_classifier.hpp"
Cohesion: 0.16
Nodes (36): build_cache_entry(), build_F3inv_only_cache(), chi_square_v32f(), clean_label(), covariance_to_correlation_v32f(), det_real_projected_F3inv_v32f(), find_all_zero_poles(), find_zero_poles_from_cache() (+28 more)

### Community 28 - "run"
Cohesion: 0.08
Nodes (41): Args, coarseN, Emax, Emin, git_commit, irrep, L, outdir (+33 more)

### Community 29 - "projections_gpu_safe.cuh"
Cohesion: 0.15
Nodes (28): Vsel, all_perms_host(), build_projector_from_eigenvectors_near_one_gpu_real(), chi_host(), conj_class_host(), array, MatrixXcd, MatrixXd (+20 more)

### Community 30 - "v33g_build_runtime_k3basis_cache.cpp"
Cohesion: 0.12
Nodes (39): BuildConfig, base, coarse_cache_root, explicit_coarse, force_rebuild, gpu_coarseN, gpu_Ecm_max, gpu_Ecm_min (+31 more)

### Community 31 - "run_block"
Cohesion: 0.09
Nodes (41): set, boolean_value(), BridgeArgs, accepted_truezeros_csv, cpu_fallback, gpu_batch_rows, progress_every, reference_csv (+33 more)

### Community 32 - "v33h_patched_gpu_cache_oldscale_det_scan.cpp"
Cohesion: 0.13
Nodes (38): assemble_scaled_projected_QC_for_row(), compute_row(), comp, FitSettings, GpuComplexReadConvention, MatrixXcd, ofstream, pair (+30 more)

### Community 33 - "EvalData"
Cohesion: 0.05
Nodes (40): EvalData, cond_F3, cond_projF3inv, detF3, detF3inv, detProjF3inv, Ecm, En (+32 more)

### Community 34 - "plot_classifier_algorithm_sweep_v33l.py"
Cohesion: 0.08
Nodes (67): RuntimeError, analyze_algorithm(), BlockResult, Candidate, cluster_candidates(), collect_matrix_diagnostics(), load_nonint_levels(), load_runtime_matrices() (+59 more)

### Community 35 - "F2_functions_v2.h"
Cohesion: 0.13
Nodes (35): debug_F2_2plus1_component_for_lm_steps(), debug_F2_ang_mom_steps(), debug_I_int_ang_mom_steps(), debug_I_sum_ang_mom_steps(), F2_2plus1_mat(), F2_ang_mom(), F2_i_ang_mom_mat(), comp (+27 more)

### Community 36 - "CpuGpuChunk"
Cohesion: 0.16
Nodes (20): append_failed_results_for_chunk(), bytes_to_mib(), CpuGpuChunk, chunk_id, Ecm_end, Ecm_start, estimated_gpu_bytes, i_end (+12 more)

### Community 37 - "functions.h"
Cohesion: 0.20
Nodes (35): config_maker(), config_maker_1(), config_maker_2(), config_maker_3(), config_maker_4(), config_maker_4_momentum_first(), config_maker_5(), cutoff_function_1() (+27 more)

### Community 38 - "AGENTS.md"
Cohesion: 0.06
Nodes (34): 10. code-reviewer, 11. report-writer, 1. repo-explorer, 2. implementation-worker, 3. validation-scientist, 4. gpu-performance-engineer, 5. cpu-openmp-reference-keeper, 6. cache-format-guardian (+26 more)

### Community 39 - "PipelineGpuCache"
Cohesion: 0.06
Nodes (31): PipelineGpuCache, cap, dDet, dF2, dF2_array, dF2X1, dF2X1_array, dF3 (+23 more)

### Community 40 - "F2_gpu_safe_builder.cuh"
Cohesion: 0.28
Nodes (33): build_F2_2plus1_kernel(), cuDoubleComplex erfi_series_d(), cuDoubleComplex F2_ang_mom_d(), cuDoubleComplex I0F_d(), cuDoubleComplex I1F_d(), cuDoubleComplex I2F_d(), cuDoubleComplex I_int_ang_mom_d(), cuDoubleComplex I_sum_ang_mom_d() (+25 more)

### Community 41 - "string"
Cohesion: 0.20
Nodes (36): bisection_eigenbranch_minabs_zero(), bisection_eigenbranch_reference_zero(), bisection_f3iso_zero(), bisection_method_zero(), classify_side_abs(), collect_det_eigenbranch_candidates_from_adjacent_pairs(), collect_zero_candidates_for_method(), pair (+28 more)

### Community 42 - "RuntimeCacheMeta"
Cohesion: 0.06
Nodes (31): size_t, RuntimeCacheMeta, alpha, atmK, atmpi, coarseN, created_by, Ecm_max (+23 more)

### Community 43 - "Test report 008: M1 step 7.0 (projector cache key), 7.1 (pole/zero discrimination), 7.8 (D7 default q0_consistent)"
Cohesion: 0.17
Nodes (11): 7.0 details (cache key), 7.1 physics review, 7.8 details (D7), Commands run (repo root, `OMP_NUM_THREADS=1`), Failures, Observations for the coordinator (not failures; guesses labelled), Pending (7 failures, all `NotImplementedError`), Raw numbers (+3 more)

### Community 44 - "ZeroRecord"
Cohesion: 0.07
Nodes (32): classify_zeros(), dedup_eigenbranch_candidate_records(), nearest_index_within(), same_coarse_interval(), unique_sorted_zeros(), zero_record_residual_for_merge(), ZeroRecord, abs_y (+24 more)

### Community 45 - "F3_gpu_omp_cublas_pipeline_v9_inward_shape_refinement.cu"
Cohesion: 0.14
Nodes (30): add_unique_zero(), build_one_energy_pack(), classify_refined_file_by_inward_shape(), classify_sign_flips_zero_vs_pole(), comp, cublasStatus_t, cudaError_t, string (+22 more)

### Community 46 - "Full F3 GPU pipeline versions"
Cohesion: 0.07
Nodes (28): Current bottleneck, Current status after v16, F3 GPU Pipeline Version History, Full F3 GPU pipeline versions, GPU-safe `config_maker_4`, GPU-safe `F2`, GPU-safe `G`, GPU-safe ingredient builders (+20 more)

### Community 47 - "v33g_runtime_k3basis_cache.hpp"
Cohesion: 0.24
Nodes (17): gpu_complex_read_convention_name(), GpuComplexReadConvention, ofstream, string, vector, meta_double(), meta_int(), meta_string() (+9 more)

### Community 48 - "KKpi Codebase Map for Codex"
Cohesion: 0.07
Nodes (27): `cache_generator` file tree, Canonical build and run commands, Change-routing table for Codex, Codex start-here prompt, Configuration option map, `corrected_reader` file tree, CPU component cache generator, CPU component cache generator (+19 more)

### Community 49 - "E — `threebody_fullmatrix_hpc_cuda_bench_v3_codex`: what it is and whether it belongs in 3body-qc-v4"
Cohesion: 0.07
Nodes (26): External data locations (source machine), QC3_release (reference Python implementation), RELATED CODES (not part of the build), Three-body full-matrix HPC/CUDA benchmark (2026-08-21), 0. Summary, 1.1 Ingredients (all in `PKG/include/threebody_physics.hpp`), 1.2 The equation and how it is turned into a matrix problem, 1.3 Parameters and energy window (+18 more)

### Community 50 - "K3_functions_2plus1.hpp"
Cohesion: 0.25
Nodes (27): boost_cpp(), dot3(), ell_from_config(), flavor_mass(), from_stdvec(), comp, MatrixXcd, vector (+19 more)

### Community 51 - "MomentumIrrepSpec"
Cohesion: 0.07
Nodes (29): append_zero_records(), array, ofstream, MomentumIrrepSpec, irrep, irrep_tag, label, nnP (+21 more)

### Community 52 - "fitter/K3df_minuit_fit_v32f_fullF3inv_QCfull_cached_classifier.hpp"
Cohesion: 0.17
Nodes (24): assemble_QC(), build_cache_entry(), build_F3inv_only_cache(), closest_zero_eigenvalue(), det_real_projected_F3inv_v32f(), find_zero_poles_from_cache(), hermitize(), comp (+16 more)

### Community 53 - "DeviceResidentResult"
Cohesion: 0.18
Nodes (11): DeviceResidentResult, build_sec, det, dim1, dim2, Ecm, En, i (+3 more)

### Community 54 - "F3_gpu_cachegen_v32zu.cu"
Cohesion: 0.17
Nodes (30): all_perms_of_3(), build_projector_from_eigenvectors_near_one(), cart_axis_from_m_v30m(), chi(), conj_class(), Convention, comp, cublasStatus_t (+22 more)

### Community 55 - "PhysicsParams"
Cohesion: 0.08
Nodes (25): PhysicsParams, alpha, atmK, atmpi, eig_tol, epsilon_h, eta_1, eta_2 (+17 more)

### Community 56 - "vector"
Cohesion: 0.20
Nodes (21): BlockInfo, file_irrep_hint, internal_irrep, L, target_indices, complex, pair, vector (+13 more)

### Community 57 - "PhysicsParams"
Cohesion: 0.08
Nodes (25): PhysicsParams, alpha, atmK, atmpi, eig_tol, epsilon_h, eta_1, eta_2 (+17 more)

### Community 58 - "Report B — F3inv/Vsel cache generation (CPU + GPU) and root-finding lineage for the projected QC determinant"
Cohesion: 0.08
Nodes (26): 1. What was asked, 2. Actions taken, 3.1 Physics objects stored per energy row, 3.2 Energy grid, 3.3 Combined raw cache `v32zu_Lbyas<L>_xi<xi>_irrep<label>_coarse<N>_F3inv_Vsel_gpu.bin`, 3.4 Combined-cache sidecar `<bin>.meta.json`, 3.5 Component caches (format `KKPI_COMPONENT_CACHE_V001`), CPU and GPU, 3.6 The `variant_04_real_imag_swapped` convention: where it comes from (+18 more)

### Community 59 - "K2_functions_gpu_safe.cuh"
Cohesion: 0.23
Nodes (23): abs_cx(), build_K2inv_2plus1_kernel(), __global__, __host__, cutoff_function_1_gpu(), Cx, im, re (+15 more)

### Community 60 - "SvdMinRecordV30r"
Cohesion: 0.14
Nodes (14): SvdMinRecordV30r, coarse_i, E_min, E_min_double, E_window_left, E_window_right, near_dimension_jump, nearest_dimension_jump (+6 more)

### Community 61 - "vector"
Cohesion: 0.10
Nodes (37): chi_square_v32f(), clean_label(), covariance_to_correlation_v32f(), diagonal_chi2(), find_all_zero_poles(), fit_K3df_parameters_v32f(), get_or_build_F3inv_cache_v32f(), FunctionMinimum (+29 more)

### Community 62 - "v32zt_generate_one_irrep"
Cohesion: 0.13
Nodes (27): cudaDeviceProp, array, PhysicsParams, string, Entry, object, os, path (+19 more)

### Community 63 - "Version history summary"
Cohesion: 0.08
Nodes (24): Conversation Context Reconstruction, `F2iso` and `Giso`, Group-theory degeneracy, Important decisions made, Key definitions, Known issues and open questions, Momentum/irrep label convention, Non-interacting three-body energy (+16 more)

### Community 64 - "KKpi v33p Cache-Anchored Classifier/Fitter Checkpoint"
Cohesion: 0.08
Nodes (24): Build, Build and smoke test, Cache format notes, Checkpoint status, Critical checkpoint, Current classifier state, Exact old scaling, External cache requirement (+16 more)

### Community 65 - "slogdet"
Cohesion: 0.23
Nodes (13): project_Q(), QC helpers (prompts/001 step 4): projection of the inverse and (sign, ln|det|)…, (L^-6) V^dag (F3inv + K3) V ; L = xi*Lbyas (so L**6 == (Lbyas*xi)**6). The…, (sign, ln|det Q|) by partial-pivot LU (LAPACK via NumPy). sign is an int in…, slogdet(), M1 step 4 — qc3.qc implements Q = (L xi)^-6 V^dag (F3inv + K3) V and (sign,…, _synthetic(), test_lndet_invariant_under_unitary_mixing_of_Vsel() (+5 more)

### Community 66 - "run_groups_concurrent"
Cohesion: 0.14
Nodes (24): EnergyPrepassInfo, EnergyShapeKey, size_t, estimate_group_vram_bytes(), get_free_vram_bytes(), run_groups_concurrent(), v32zt_auto_batch_cap(), v32zt_baseline_rows_per_second() (+16 more)

### Community 67 - "QCPointValue"
Cohesion: 0.10
Nodes (22): assemble_QC(), closest_zero_eigenvalue(), hermitize(), MatrixXcd, K3dfParameters, K3B, K3E, K3iso0 (+14 more)

### Community 68 - "FitResult"
Cohesion: 0.12
Nodes (18): FitResult, best, chi2, chi2_dof, data_correlation, data_covariance, errors, model_levels (+10 more)

### Community 69 - "ProjectedQCCacheEntry"
Cohesion: 0.08
Nodes (27): canonical_shell_momentum(), array, infer_row_spec_index_from_covariance_order_v32f(), ProjectedQCCacheEntry, Ecm, En, En_c, error (+19 more)

### Community 70 - "inspect_v32_f3inv_cache_v32f_fixed.py"
Cohesion: 0.23
Nodes (22): CacheEOF, CacheFormatError, expected_coarse_grid(), file_size(), has_near(), main(), parse_cache(), Exception (+14 more)

### Community 71 - "load_multil_jack_targets"
Cohesion: 0.12
Nodes (32): Cand, assign_block_model(), build_assignment_rows(), CandidateWithBlock, c, file_irrep, internal_irrep, K3dfParameters (+24 more)

### Community 72 - "vector"
Cohesion: 0.12
Nodes (21): classify_sign_flips_zero_vs_pole(), string, vector, median_of_vector(), SignFlipCandidate, abs_gap, classification, Ecm1 (+13 more)

### Community 73 - "make_qc_state_lines_from_scattered_spectrum.py"
Cohesion: 0.22
Nodes (18): Namespace, collect_dat_files(), deduplicate_points(), group_by_label_L(), interp_track(), main(), make_plot(), match_tracks_to_points() (+10 more)

### Community 74 - "BlockSummary"
Cohesion: 0.09
Nodes (22): BlockSummary, avg_abs_det_cpu, avg_abs_det_gpu, avg_raw_det_cpu_real, avg_raw_det_gpu_real, cpu_fail_rows, duplicate_rows_removed, gpu_fail_rows (+14 more)

### Community 75 - "Cand"
Cohesion: 0.10
Nodes (22): Cand, BL, BR, depth, E, FL, FR, init (+14 more)

### Community 76 - "parse_label"
Cohesion: 0.12
Nodes (25): F3inv, parse_label(), block_cnorm(), build_runtime_block(), L, FitSettings, IrrepCache, PhysicsParams (+17 more)

### Community 77 - "`corrected_reader`"
Cohesion: 0.10
Nodes (21): `cache_generator`, `configs/`, `configs/`, `corrected_reader`, `cpu_component_cachegen/scripts/`, `cpu_component_cachegen/source/`, `diagnostics/`, `diagnostics/` (+13 more)

### Community 78 - "F3_device_resident_true_batch_v21_side_abs_trend_refine.cu"
Cohesion: 0.22
Nodes (18): build_F3_from_F2_and_F2X1_kernel(), c_make(), combine_H_colmajor_kernel(), ConfigEntry, cuDoubleComplex, __global__, gather_real_evecs_to_complex_Vsel_kernel(), v14_build_F2_batch_same_shape_kernel() (+10 more)

### Community 79 - "ProjectedQCCacheEntry"
Cohesion: 0.10
Nodes (21): ProjectedQCCacheEntry, Ecm, En, En_c, error, F3_equiv, F3_leak, F3inv_full (+13 more)

### Community 80 - "F3_device_resident_true_batch_v16_custom_batched_getrs.cu"
Cohesion: 0.22
Nodes (18): build_F3_from_F2_and_F2X1_kernel(), c_make(), combine_H_colmajor_kernel(), ConfigEntry, cuDoubleComplex, __global__, gather_real_evecs_to_complex_Vsel_kernel(), v14_build_F2_batch_same_shape_kernel() (+10 more)

### Community 81 - "F3_cpu_openmp_v25_K3QC_cached_core.hpp"
Cohesion: 0.28
Nodes (16): add_unique_zero(), classify_one_subflip_side_abs(), evaluate_one_energy_cpu_eigen_v24_QC(), find_likely_zeros_v24_K3QC_silent(), find_likely_zeros_v25_K3QC_cached(), find_sign_flips(), pair, vector (+8 more)

### Community 82 - "QCRefinedFCN"
Cohesion: 0.10
Nodes (22): cov_to_corr(), atomic, MatrixD, Minuit2::FCNBase, shared_ptr, TargetLevel, QCRefinedFCN, cn (+14 more)

### Community 83 - "gpu_cache_reader_convention.hpp"
Cohesion: 0.27
Nodes (19): GpuCacheFileFormat, format_from_header(), comp, GpuComplexReadConvention, istream, string, uint64_t, name() (+11 more)

### Community 84 - "EnergyGpuResult"
Cohesion: 0.10
Nodes (20): EnergyGpuResult, det_F3inv_projected, Ecm, eig_F3inv_projected, En, F2, F3, F3inv_projected (+12 more)

### Community 85 - "QCRefinedFCN"
Cohesion: 0.11
Nodes (17): atomic, Minuit2::FCNBase, shared_ptr, TargetLevel, QCRefinedFCN, cn, coarse_cache, corr (+9 more)

### Community 86 - "bench.py"
Cohesion: 0.67
Nodes (3): main(), Benchmark driver (prompts/001 step 8)., run()

### Community 87 - "V32ztAllObjRowPayload"
Cohesion: 0.10
Nodes (20): V32ztAllObjRowPayload, Ecm, F2, F3, F3inv, G, grid_i, K2inv (+12 more)

### Community 88 - "CachedProjectorV30q"
Cohesion: 0.07
Nodes (34): basis_signature_v30q(), build_cached_projector_v30q(), CachedProjectorV30q, fv_best_convention, fv_rep_closure_best, fv_rep_unitarity, P_I_herm_res, P_I_idem_res (+26 more)

### Community 89 - "Tester report 004: M1 step 4, QC helpers + lattice reader (D1, D3)"
Cohesion: 0.25
Nodes (7): Commands run, Failures, Independent physics / statistics checks, Notes (non-blocking; for the builder or auditor, no fixes done by me), Test suites, Tester report 004: M1 step 4, QC helpers + lattice reader (D1, D3), Verdict: PASS (with 3 non-blocking notes)

### Community 91 - "FinalBracket"
Cohesion: 0.11
Nodes (19): FinalBracket, BL, BR, depth, Ezero, flipL, flipR, id (+11 more)

### Community 92 - "Cand"
Cohesion: 0.11
Nodes (19): Cand, BL, BR, depth, E, FL, FR, init (+11 more)

### Community 93 - "L20_spectrum_plotter_model_true_zeros_yaxis_fixed.py"
Cohesion: 0.20
Nodes (18): find_nonint_script(), generate_recent_nonint_levels(), main(), _parse_bool(), _parse_float(), parse_nonint_output_file(), plot_spectrum(), Parse v32w bestfit_QC_spectrum.dat. Expected columns: label index Ecm kind… (+10 more)

### Community 94 - "QC_functions_v2.h"
Cohesion: 0.19
Nodes (18): boost(), G_2plus1_mat(), G_ij_lm(), get_parity_block_from_plm_config(), comp, MatrixXcd, vector, print_Gij_boosts() (+10 more)

### Community 95 - "fitter/qc_fitter_norm_refine_v2.cpp"
Cohesion: 0.17
Nodes (20): classify(), Eval, det, E, proj_dim, success, y, find_QC_zeros_v3_from_coarse_grid() (+12 more)

### Community 96 - "AcceptedZeroWindow"
Cohesion: 0.11
Nodes (19): AcceptedZeroWindow, bracket_id, center_row, E_left_bracket, E_right_bracket, has_previous_model_zero, inside_Ecm_cutoff, irrep (+11 more)

### Community 97 - "inspect_v32_f3inv_cache_mmap_progress.py"
Cohesion: 0.29
Nodes (18): expected_grid(), f64(), has_near(), i32(), inspect_cache(), main(), need(), parse_one_entry() (+10 more)

### Community 98 - "3. Results"
Cohesion: 0.14
Nodes (13): 007 — Single-energy QC: independent reference, speed-up candidates, p-wave $\mathcal K_2$ convention, 1. What was asked, 2. Actions taken, 3.1 The oracle: existing CPU path, 1 core, 1 energy, 3.2 Independent reference for $F$, $G$, $\mathcal K_2^{-1}$ (NumPy, from the paper), 3.3 Independent reference for $\mathcal K_{\rm df,3}$ (numerical angular projection), 3.4 Sensitivity of $\ln|\det Q|$ → tolerances, 3.5 Cost of the $F$ sum and the truncation bound (+5 more)

### Community 99 - "F3MatrixSaveOptions"
Cohesion: 0.11
Nodes (18): F3Hdf5OpenMode, F3MatrixSaveOptions, debug, Ecm_max, Ecm_min, enabled, filename, gzip_level (+10 more)

### Community 100 - "SignFlipCandidate"
Cohesion: 0.11
Nodes (18): SignFlipCandidate, classification, Ecm_left, Ecm_mid, Ecm_right, Ecm_window_left, Ecm_window_right, Ecm_zero_linear (+10 more)

### Community 101 - "FitResult"
Cohesion: 0.12
Nodes (18): FitResult, best, chi2, chi2_dof, data_correlation, data_covariance, errors, model_levels (+10 more)

### Community 102 - "PhysicsParams"
Cohesion: 0.11
Nodes (18): PhysicsParams, alpha, atmK, atmpi, epsilon_h, eta_1, eta_2, K3B_par (+10 more)

### Community 103 - "DeviceConfig"
Cohesion: 0.12
Nodes (14): ConfigView, ell, m, n, px, py, pz, DeviceConfig (+6 more)

### Community 104 - "find_QC_zeros_refined"
Cohesion: 0.25
Nodes (17): build_entry(), FitSettings, IrrepCache, K3dfParameters, pair, PhysicsParams, ProjectedQCCacheEntry, vector (+9 more)

### Community 105 - "Codex Report: Full Cache Generator Package"
Cohesion: 0.12
Nodes (16): Benchmark validation, Build validation, Changes made earlier in this session, Codex Report: Full Cache Generator Package, Current status, Decision tree used, Files changed for the package, Files included in the package (+8 more)

### Community 106 - "ClassifyResult"
Cohesion: 0.12
Nodes (17): ClassifyResult, classification, E_linear, edge_abs_scale, flip_index, j, jp1, left_dec_frac (+9 more)

### Community 107 - "InwardShapeResult"
Cohesion: 0.12
Nodes (17): InwardShapeResult, classification, decreasing_fraction, decreasing_steps, Ecm_left, Ecm_right, Ecm_zero_linear, flat_steps (+9 more)

### Community 108 - "BuiltRow"
Cohesion: 0.12
Nodes (17): BuiltRow, B0, B1, BB, BE, det_basis, det_diff, det_direct (+9 more)

### Community 109 - "ScanArgs"
Cohesion: 0.11
Nodes (18): K3dfParameters, vector, ScanArgs, all_default_blocks, coarseN, complex_read_convention, Emax, Emin (+10 more)

### Community 110 - "ClassifyResult"
Cohesion: 0.12
Nodes (17): ClassifyResult, classification, E_linear, edge_abs_scale, flip_index, j, jp1, left_dec_frac (+9 more)

### Community 111 - "roots.py"
Cohesion: 0.12
Nodes (27): MissingRootError, Raised when a root is not bracketed / fewer QC levels than lattice levels., AmbiguousRootError, _brent(), _check_kind(), _classify(), classify_sign_change(), _Counter (+19 more)

### Community 112 - "METHODS — numerical methods, equations, conventions, sources"
Cohesion: 0.11
Nodes (18): 0. System and formalism, 10. Energy grid, 11. Data, covariance, χ² — `src/fitter/qc_fitter_norm_refine_v2_multiL.cpp` (`load_multil_jack_targets`), `src/fitter/lattice_data_covariance_cpp.hpp`, `chi_square_v32f`, 12. Minimiser, 13. Root finding, 14. Non-interacting levels — `tools/nonint3body_general_doublecover/`, 15. Cache-free single-energy evaluation: oracle, fast path, tolerances — `src/python/qc3/`, `src/qc3_core/` (plan `prompts/001`), 16. Independent reference implementation — `tests/acceptance/m1/qc3ref/` (planner-owned) (+10 more)

### Community 113 - "3. Findings"
Cohesion: 0.12
Nodes (15): 1. What was asked, 2. Actions taken, 3.1 Lattice spectrum data: what exists, 3.2 Critical: the production fitter inflates the data errors by a factor of 602, 3.3 User true-zero label CSVs and accepted windows, 3.4 Config files and keys, 3.5 Two-body input ($K_2^{-1}$) and where the scattering parameters come from, 3.6 $K_{\rm df,3}$ parametrization (exact, from `PKG/source/K3_functions_2plus1.hpp`) (+7 more)

### Community 114 - "debug_v30r_svd_validated_cartesian_l1_100_A2_zero_scan.cpp"
Cohesion: 0.09
Nodes (45): eigmin_point_from_eval_v30s(), main(), main(), BranchPoint, eigvecs, success, coarse_eigensystems_for_method(), collect_svd_local_min_candidates_v30r() (+37 more)

### Community 115 - "Faddeeva.cc"
Cohesion: 0.23
Nodes (11): cmplx, cpolar(), erfcx_y100(), FADDEEVA(), FADDEEVA_RE(), main(), relerr(), sinc() (+3 more)

### Community 116 - "v33d CPU/GPU Component Cache Generator Package"
Cohesion: 0.13
Nodes (14): Build and run, CPU, CPU path, Current status at packaging time, Decision tree used to get to the current generator, Execution tree, GPU, GPU path (+6 more)

### Community 118 - "EvalPoint"
Cohesion: 0.19
Nodes (15): classify_final_v32s(), pair, vector, EvalPoint, det, E, from_refinement, proj_dim (+7 more)

### Community 119 - "v32w_single_L_fitter/qc_fitter_norm_refine_v2.cpp"
Cohesion: 0.13
Nodes (21): classify(), MatrixXcd, det_info(), DetInfo, det_re, logabs, ok, sign (+13 more)

### Community 120 - "ConfigMaker4Params"
Cohesion: 0.17
Nodes (12): ConfigMaker4Params, epsilon_h, L, max_num_shell, mi, mj, mk, Px (+4 more)

### Community 121 - "run_v33m_classifier_speed_search.py"
Cohesion: 0.29
Nodes (13): assignment_rows(), assignment_signature(), BenchRow, compare_assignments(), main(), parse_summary(), Path, run() (+5 more)

### Community 122 - "RowResult"
Cohesion: 0.13
Nodes (15): RowResult, abs_diff, cache_index, det_cpu, det_gpu, Ecm, irrep, Lbyas (+7 more)

### Community 123 - "process_same_shape_group_true_batch_v15"
Cohesion: 0.18
Nodes (26): adaptive_final_zeros_device_resident(), build_projection_Vsel_device(), classify_inward_shape(), compute_one_energy_device_resident(), array, comp, cublasHandle_t, string (+18 more)

### Community 124 - "ofstream"
Cohesion: 0.20
Nodes (12): ofstream, T, time_point, v32zt_write_cache_record(), v32zt_write_exact(), v32zt_write_matrix_record(), v32zt_write_scalar(), v32zt_write_staged_allobj_rows() (+4 more)

### Community 126 - "F3iProjectedCacheEntry"
Cohesion: 0.11
Nodes (18): build_F3i_projected_cache_entry_v25(), determinant_via_partial_piv_lu(), F3iProjectedCacheEntry, dim1, dim2, Ecm, En, error (+10 more)

### Community 127 - "GroupStream"
Cohesion: 0.15
Nodes (14): cudaStream_t, cublasHandle_t, GroupStream, batch_size, completed, dim1, dim2, estimated_vram_bytes (+6 more)

### Community 128 - "Full L20 fit and validation report"
Cohesion: 0.14
Nodes (13): Accepted-window package, Chi-square audit, Classifier validation, Errors and correlations, Exact reproducibility inputs and outputs, Full L20 fit and validation report, K3B/K3E zero test, L20 pass-only fits (+5 more)

### Community 129 - "KKpi v33o Production Fitter"
Cohesion: 0.14
Nodes (13): Build and smoke test, Cache format notes, Critical checkpoint, Current classifier state, Exact old scaling, KKpi v33o Production Fitter, Known limits, Main workflow commands (+5 more)

### Community 130 - "F3_cpu_openmp_v24_K3_QC_multi_signflip_side_abs_refinement.cpp"
Cohesion: 0.38
Nodes (13): add_unique_zero(), classify_one_subflip_side_abs(), pair, string, vector, evaluate_one_energy_cpu_eigen_v24_QC(), find_sign_flips(), L_tag() (+5 more)

### Community 131 - "EnergyPrepassInfo"
Cohesion: 0.25
Nodes (8): EnergyPrepassInfo, dim1, dim2, Ecm, En, i, total_dim, vdim

### Community 132 - "PhysicsParams"
Cohesion: 0.14
Nodes (14): PhysicsParams, alpha, atmK, atmpi, epsilon_h, eta_1, eta_2, Lbyas (+6 more)

### Community 133 - "generate_K3iso.cpp"
Cohesion: 0.21
Nodes (12): comp, string, generate_K3iso_L20_from_F3(), generate_K3iso_L20_from_F3inv(), generate_K3iso_using_F3_vs_En_KKpi_omp_single_irrep_centralvalue(), generate_K3iso_using_F3_vs_En_KKpi_omp_single_irrep_with_bounds(), main(), det_F3mat_poles_secant_method() (+4 more)

### Community 134 - "5. Key physics and numerics lessons (Question 2)"
Cohesion: 0.14
Nodes (14): 5.10 Covariance fix (v32l), cache-only fix (v32k), namespace fixes (v31z, v32i, v32r, v32y), 5.11 Lattice input conventions (v32p → v32s), 5.12 Memory (v32x → v32y), 5.13 Scattering-parameter plumbing, 5.1 Use slogdet, never the raw determinant, 5.2 Projection audit: the projector was not the problem, 5.3 Pole vs zero: a sign flip is not a level, 5.4 Inverse convention: project the QC, not $F_3$ (+6 more)

### Community 135 - "CombinedRawGpuCacheRowHeader"
Cohesion: 0.14
Nodes (14): CombinedRawGpuCacheRowHeader, dim1, dim2, Ecm, En, grid_i, n, vdim (+6 more)

### Community 136 - "V32ztCacheRowPayload"
Cohesion: 0.14
Nodes (15): cuDoubleComplex, __global__, v32zt_identity_batch_kernel(), V32ztCacheRowPayload, dim1, dim2, Ecm, En (+7 more)

### Community 137 - "V32ztTiming"
Cohesion: 0.14
Nodes (14): V32ztTiming, basis_sec, bytes_written, chunks, f2k2_sec, gf3_build_sec, max_total_dim, max_vdim (+6 more)

### Community 138 - "map"
Cohesion: 0.33
Nodes (12): cnorm(), string, gd(), gi(), gs(), is_v3_like_mode(), is_v4_like_mode(), main() (+4 more)

### Community 139 - "ThreadSafeQueue"
Cohesion: 0.18
Nodes (13): condition_variable, bytes_complex_matrix(), mutex, T, estimate_gpu_bytes_one_energy(), flush_cpu_buffer_to_gpu_queue(), gpu_consumer_thread_func(), ThreadSafeQueue (+5 more)

### Community 140 - "OPEN ISSUES — read before trusting any number"
Cohesion: 0.14
Nodes (14): 10. LOW — minor physics-core bugs (not on the production path), 11. LOW — portability, 12. Historical warnings (do not reuse these results), 13. LOW — legacy projector cache key omits the selection tolerances (`basis_signature_v30q`), 1. BLOCKER (CONFIRMED by D1) — lattice covariance too large by $(N-1)^2=602^2\approx3.6\times10^5$, 2. HIGH — p-wave $\mathcal K_2^{-1}$ differs from the reference implementation, 3. HIGH — complex-number read convention (`variant_04_real_imag_swapped`), 4. HIGH — inputs needed to re-run fits are not in this repo (+6 more)

### Community 141 - ".upload"
Cohesion: 0.30
Nodes (8): host_comp, MatrixXcd, vector, DeviceScatterParams, d_im, d_re, K2inv_EREord2_2plus1_many_mat_gpu_safe(), K2inv_EREord2_2plus1_mat_gpu_safe()

### Community 142 - "K3iso_maker.py"
Cohesion: 0.29
Nodes (10): E_to_Ecm(), Ecmsq_to_Esq(), Esq_to_Ecmsq(), jackknifeavg_centralvalue_lattice_data(), jackknifeavg_lattice_data(), P000(), P100(), P110() (+2 more)

### Community 143 - "scan_projF3inv_digonto_classifier_v2.cpp"
Cohesion: 0.41
Nodes (11): const_norm_scale_v32s(), string, TargetLevel, gd(), gi(), gs(), main(), read_kv_v32s() (+3 more)

### Community 144 - "IrrepResult"
Cohesion: 0.15
Nodes (13): InitialFlip, El, Er, id, il, ir, yl, yr (+5 more)

### Community 145 - "TargetLevel"
Cohesion: 0.13
Nodes (17): canonical_shell_momentum(), diagonal_chi2(), array, infer_row_spec_index_from_covariance_order_v32f(), same_momentum_shell_v32f(), TargetLevel, atP, E_read (+9 more)

### Community 146 - "string"
Cohesion: 0.29
Nodes (13): cnorm(), cov_to_corr(), MatrixD, string, gd(), gi(), gs(), main() (+5 more)

### Community 147 - "V32zuRunOptions"
Cohesion: 0.15
Nodes (13): V32zuRunOptions, concurrent_streams_safety_fraction, device_count, gpu_memory_safety_fraction, irrep_filter_csv, max_batch_energies, use_concurrent_streams, use_fused_f2k2 (+5 more)

### Community 148 - "process_same_shape_group_true_batch_v15"
Cohesion: 0.45
Nodes (11): build_projection_Vsel_device(), compute_one_energy_device_resident(), comp, cublasHandle_t, vector, cublas_lu_solve_inplace_batched_same_n(), determinant_colmajor_device_matrix(), device_lu_solve_inplace() (+3 more)

### Community 149 - "002 — M1 step 7 re-planned: cold projector from orbits, threshold basis builder, pair kernels for $\mathcal K_{\rm df,3}$ and $G$, C++ tail, and pole/zero discrimination in the root refiner"
Cohesion: 0.09
Nodes (22): 002 — M1 step 7 re-planned: cold projector from orbits, threshold basis builder, pair kernels for $\mathcal K_{\rm df,3}$ and $G$, C++ tail, and pole/zero discrimination in the root refiner, 0. Diagnosis — why step 7 is re-planned, 1. Goal of step 7 (revised), 2.1 Projector from orbits (rung 7.2) — `src/qc3_core/fast_projector.hpp`, 2.2 Threshold basis builder (rung 7.3) — `src/qc3_core/fast_basis.hpp`, 2.3 Pair kernels for $\mathcal K_{\rm df,3}$ and $G$ (rungs 7.4, 7.5) — `src/qc3_core/fast_k3.hpp`, `fast_g.hpp`, 2.4 C++ tail, LU solve, lazy conversion (rung 7.6), 2.5 Pole/zero discrimination (rung 7.1) — `src/python/qc3/roots.py` (+14 more)

### Community 150 - "EvalResult"
Cohesion: 0.17
Nodes (12): EvalResult, det, dim1, dim2, Ecm, eig_min_proj, En, error (+4 more)

### Community 151 - "ZeroPole"
Cohesion: 0.17
Nodes (12): ZeroPole, E_left, E_mid, E_right, eig_index_left, eig_index_right, i_left, i_right (+4 more)

### Community 152 - "functions_gpu_config_maker4.cuh"
Cohesion: 0.47
Nodes (8): __global__, __host__, en_min_plus_real(), f4_pi(), mark_config_maker_4_candidates_kernel(), norm3(), omega_real(), selected_by_config_maker_4_condition()

### Community 153 - "gpu_config_maker_4_many_entries"
Cohesion: 0.42
Nodes (12): build_lm_list_from_waves(), candidate_count_per_energy(), compact_flags_to_entries_cpu_order(), config_maker_4_gpu_safe(), complex, size_t, vector, entries_to_cpu_vectors() (+4 more)

### Community 154 - "001 — M1: a cache-free, CPU-only, single-energy three-body QC with an oracle and an independent reference"
Cohesion: 0.07
Nodes (27): 001 — M1: a cache-free, CPU-only, single-energy three-body QC with an oracle and an independent reference, 0. Diagnosis (why the old pipeline is slow and what is wrong), 10. Expected cost, 1. Goal of M1, 2.1 Equations and conventions (unchanged physics; see `docs/METHODS.md` §§1–8 for the full set), 2.2 p-wave $\mathcal K_2^{-1}$: two named forms, production default unchanged until D6, 2.3 Speed-up candidates, ranked (expected gain, how verified), 2.4 Lattice data (D1, D3) (+19 more)

### Community 155 - "EvalResult"
Cohesion: 0.17
Nodes (12): EvalResult, det, dim1, dim2, Ecm, eig_min_proj, En, error (+4 more)

### Community 156 - "QCPointValue"
Cohesion: 0.17
Nodes (12): QCPointValue, det_re, det_sign, Ecm, eig_index, herm_rel, logabsdet, min_abs_eig (+4 more)

### Community 157 - "ZeroPole"
Cohesion: 0.17
Nodes (12): ZeroPole, E_left, E_mid, E_right, eig_index_left, eig_index_right, i_left, i_right (+4 more)

### Community 158 - "PIPELINE — from physics parameters to a $\mathcal K_{\rm df,3}$ fit"
Cohesion: 0.11
Nodes (14): (A) Cache generation, (B) Determinant scan, (C) Labelling (the current gate for "which roots exist"), Combined raw cache format (`*_F3inv_Vsel_gpu.bin`), (D) Accepted windows, (E) Hot-window projected-basis cache, (F) Fit, (G) Audits and plots (+6 more)

### Community 159 - "KKpi L20 user-validated hot-window fitter"
Cohesion: 0.18
Nodes (10): Build, Current status, Included and excluded, KKpi L20 user-validated hot-window fitter, Known blockers and next action, Layout, Physics and production defaults, Plotting (+2 more)

### Community 160 - "KKpi v33o production fitter"
Cohesion: 0.18
Nodes (10): Build and smoke test, Canonical workflow, Generated outputs, How to modify the workflow, Important parameters, KKpi v33o production fitter, Known limits, Logic flow (+2 more)

### Community 161 - "ScanPoint"
Cohesion: 0.18
Nodes (11): ScanPoint, abs_norm, det, Ecm, En, i, success, total_dim (+3 more)

### Community 162 - "make_refined_mesh_v32s"
Cohesion: 0.36
Nodes (11): almost_same_E(), build_refined_entry_v32s(), FitSettings, IrrepCache, PhysicsParams, ProjectedQCCacheEntry, eval_from_entry(), linear_zero() (+3 more)

### Community 163 - "JumpRow"
Cohesion: 0.18
Nodes (11): JumpRow, delta_Nfull, delta_Nproj, Ecm, jump_id, Nfull, Nproj, prev_Nfull (+3 more)

### Community 164 - "Input"
Cohesion: 0.09
Nodes (22): Input, a0_KK, a0_Kpi, a1_Kpi, alpha, atmK, atmpi, c (+14 more)

### Community 165 - "Tester report 003: M1 step 3, reference checks and the `q0_consistent` p-wave form"
Cohesion: 0.18
Nodes (10): 4b. Term-by-term check of `K2_inv_ERE_ang_mom_q0_consistent` against prompts/001 §2.2, 4d. Independent implementation, 5. Max relative Frobenius deviation oracle vs `qc3ref` ($\lVert a-b\rVert_F/\lVert b\rVert_F$, reference in the denominator, $\vec c=0$), 6. Fast tally (`pytest tests/acceptance/m1 -m "not slow" -rA --tb=line`), Commands run, Failures, Notes (not failures), Results (+2 more)

### Community 166 - "AssignmentRow"
Cohesion: 0.18
Nodes (11): AssignmentRow, assignment_status, classifier_candidate_id, irrep, Lbyas, mode, model_root_Ecm, residual (+3 more)

### Community 167 - "KKpi v32w + v32y reproducibility package"
Cohesion: 0.20
Nodes (9): Cache inspection, Continuous QC state-line interpolation, Environment, KKpi v32w + v32y reproducibility package, Main trusted fitter: v32w, Memory-bounded spectrum-vs-L: v32y, Plotting helpers, Version guidance (+1 more)

### Community 168 - "Key numbers"
Cohesion: 0.20
Nodes (9): CPU subset validation, Current checkpoint, Gate 7 plot, Identity and cache reader, Key numbers, Passed validations, Production benchmark, Residual risk (+1 more)

### Community 169 - "InwardShapeRefineOptions"
Cohesion: 0.20
Nodes (10): InwardShapeRefineOptions, coarse_window_left, coarse_window_right, debug, duplicate_Ecm_tolerance, initial_refine_points, max_refine_rounds, refine_growth_factor (+2 more)

### Community 170 - "3.3 Question 2 — equations implemented"
Cohesion: 0.12
Nodes (17): 3.1 Formalism and references (as stated in the code), 3.2 Question 1 — the canonical, newest physics core, 3.3.1 Cutoff function $H$ (`L20/functions.h:97-141`), 3.3.2 Basis: spectator momentum × $(\ell,m)$, shell cutoff (`L20/functions.h:976-1113`), 3.3.3 $F_2$ — sum minus integral (`L20/F2_functions_v2.h`), 3.3.4 $G$ — switch-state (exchange) propagator (`L20/G_functions_v2.h`), 3.3.5 $\mathcal K_2^{-1}$ (`L20/K2_functions_v2.h`), 3.3.6 $F_3$ and $F_3^{-1}$ (+9 more)

### Community 171 - "L20 full cutoff-dependence final report"
Cohesion: 0.20
Nodes (9): Ecm_cutoff=0.315, Ecm_cutoff=0.325, Ecm_cutoff=0.335, Ecm_cutoff=0.345, Interpretation and caveats, L20 full cutoff-dependence final report, Per-cutoff results, Prefit coverage (+1 more)

### Community 172 - "run_v33p_classifier_sweep_L20_100_A2.py"
Cohesion: 0.51
Nodes (9): f(), linear_zero(), main(), previous_matrix_row(), Path, sign(), sign_flips(), v3_classify() (+1 more)

### Community 173 - "ScanPoint"
Cohesion: 0.18
Nodes (11): ScanPoint, abs_norm, det, Ecm, En, i, success, total_dim (+3 more)

### Community 174 - "SideClass"
Cohesion: 0.20
Nodes (10): SideClass, classifier, left_decreasing, left_increasing, left_score, likely_pole, likely_zero, right_decreasing (+2 more)

### Community 175 - "load_binary_f3inv_cache_v32f"
Cohesion: 0.64
Nodes (9): comp, istream, load_binary_f3inv_cache_v32f(), read_binary_raw_v32f(), read_comp_v32f(), read_matrix_v32f(), read_nested_vec_comp_v32f(), read_string_v32f() (+1 more)

### Community 176 - "save_binary_f3inv_cache_v32f"
Cohesion: 0.58
Nodes (9): ostream, T, save_binary_f3inv_cache_v32f(), write_binary_raw_v32f(), write_comp_v32f(), write_matrix_v32f(), write_nested_vec_comp_v32f(), write_string_v32f() (+1 more)

### Community 177 - "find_QC_zeros_refined"
Cohesion: 0.56
Nodes (9): build_entry(), FitSettings, K3dfParameters, PhysicsParams, ProjectedQCCacheEntry, eval_entry_QC(), find_QC_zeros_refined(), mesh() (+1 more)

### Community 178 - "Tester report 002: M1 step 2, oracle wrapper"
Cohesion: 0.20
Nodes (9): Check 5: review of `oracle_point.hpp` against `build_component_entry` (cachegen L223-288), Check 6 table (16 cases; full data in the JSON), Check 7: single-thread timings (`PointResult.timings`, $\vec c=\vec c_{\rm base}$, seconds), Commands run, Failures, Notes (not failures; for the auditor / planner), Results, Tester report 002: M1 step 2, oracle wrapper (+1 more)

### Community 179 - "L20 K3B/K3E zero tests and chi-square audit"
Cohesion: 0.22
Nodes (8): Chi-square normalization and FCN audit, Exact source audit points, L20 K3B/K3E zero tests and chi-square audit, Required output files, Table A, Table B, Table C, Table D

### Community 180 - "comp"
Cohesion: 0.50
Nodes (3): comp, MatrixXcd, determinant_via_partial_piv_lu()

### Community 181 - "SignFlip"
Cohesion: 0.50
Nodes (4): SignFlip, E_linear, left, right

### Community 182 - "save_binary_f3inv_cache_v32f"
Cohesion: 0.58
Nodes (9): ostream, T, save_binary_f3inv_cache_v32f(), write_binary_raw_v32f(), write_comp_v32f(), write_matrix_v32f(), write_nested_vec_comp_v32f(), write_string_v32f() (+1 more)

### Community 183 - "DetInfo"
Cohesion: 0.22
Nodes (9): MatrixXcd, det_info(), DetInfo, det_im, det_re, logabs, ok, sign (+1 more)

### Community 184 - "config_gpu_component_cachegen_sep_f2k2.sh"
Cohesion: 0.25
Nodes (4): config_gpu_component_cachegen_sep_f2k2_coarse20000_L20_000_A1m.sh script, config_gpu_component_cachegen_sep_f2k2_coarse20000.sh script, config_gpu_component_cachegen_sep_f2k2.sh script, config_gpu_component_cachegen.sh script

### Community 185 - "Final shutdown handoff: L20 user-validated hot-window fitter"
Cohesion: 0.25
Nodes (7): Checkpoint, Completed artifacts, Exact next action, Final shutdown handoff: L20 user-validated hot-window fitter, Physics and defaults, Restart commands, Restore

### Community 186 - "Classifier Restart Notes"
Cohesion: 0.25
Nodes (7): Candidate brackets, Classifier Restart Notes, Current state, Dimension jump, Files to inspect next, Input used, Rule to keep

### Community 187 - "Function Map"
Cohesion: 0.25
Nodes (7): Classifier core, Combined raw cache reader, Dispatch and speed checks, Function Map, Isotropic diagnostics, Production fitter orchestration, Runtime cache writer / reader

### Community 188 - "Core function inventory with descriptions"
Cohesion: 0.25
Nodes (8): `cache_generator/cpu_component_cachegen/source/generate_cpu_component_cache.cpp`, `cache_generator/gpu_component_cachegen/source/F3_gpu_cachegen_v32zu.cu`, Core function inventory with descriptions, `production_fitter/source/digonto_classifier_v3.hpp`, `production_fitter/source/extract_F3inv_isotropic_from_gpu_cache_v33a.cpp`, `production_fitter/source/qc_fitter_norm_refine_v2_multiL.cpp`, `production_fitter/source/v33g_build_runtime_k3basis_cache.cpp`, `production_fitter/source/v33g_runtime_k3basis_cache.hpp`

### Community 189 - "DetInfo"
Cohesion: 0.25
Nodes (8): MatrixXcd, det_info(), DetInfo, det_re, logabsdet, ok, sign, signed_logabsdet

### Community 190 - "__device__"
Cohesion: 0.15
Nodes (24): __device__, build_PI_kernel(), cubic_transform_dev(), __global__, D_lookup_dev(), orbit_key_dev(), sign_dev(), sort2_int() (+16 more)

### Community 191 - "Output"
Cohesion: 0.09
Nodes (22): string, Output, basis, E_lab, Ecm, F3inv, G, has_k3_basis (+14 more)

### Community 192 - "STATUS"
Cohesion: 0.25
Nodes (8): Decisions answered at the cloud session (2026-10-22), Decisions from Digonto (answered 2026-10-21), Direction for the cloud session (D4, Digonto's words, paraphrased), Known-unverified, Next action, Notes for the M1 auditor (from tester 002), Open failures and attempt counters, STATUS

### Community 193 - "main"
Cohesion: 0.61
Nodes (7): active_files(), all_config_files(), equal(), main(), parse(), Path, write_audit()

### Community 194 - "Kin"
Cohesion: 0.09
Nodes (27): Clock, t0_, compute_k3(), evaluate(), array, comp, Mat, PhysicsParams (+19 more)

### Community 195 - "F_channel"
Cohesion: 0.15
Nodes (19): MatR, evaluate(), F_channel(), comp, Mat, string, vector, max_abs_imag() (+11 more)

### Community 196 - "ComponentMatrixCacheHeader"
Cohesion: 0.25
Nodes (8): ComponentMatrixCacheHeader, coarseN, created_by, irrep, Lbyas, object_name, row_count, xi

### Community 197 - "EnergyPrepassInfo"
Cohesion: 0.25
Nodes (8): EnergyPrepassInfo, dim1, dim2, Ecm, En, i, total_dim, vdim

### Community 198 - "K2inv_EREord2_2plus1_mat"
Cohesion: 0.25
Nodes (16): comp, MatrixXcd, vector, K2_inv_00(), K2_inv_ERE_ang_mom(), K2inv_EREord2_2plus1_mat(), K2inv_EREord2_i_mat(), comp (+8 more)

### Community 199 - "ClassifierDispatchInfo"
Cohesion: 0.25
Nodes (8): classifier_dispatch_info(), ClassifierDispatchInfo, alias, dispatch_unique, implementation_tag, requested_mode, resolved_mode, version_tag

### Community 200 - "DeviceResidentResult"
Cohesion: 0.18
Nodes (11): DeviceResidentResult, build_sec, det, dim1, dim2, Ecm, En, i (+3 more)

### Community 201 - "build.sh"
Cohesion: 0.42
Nodes (8): build_cachegen_cpu(), build_core(), build_fitter(), build_hotwindow(), build_runtime(), build_scan(), have_minuit2(), build.sh script

### Community 202 - "Version summary and working status"
Cohesion: 0.29
Nodes (6): Available older packages, Best current versions, v32w known good fit, Version summary and working status, Why v32w is trusted, Why v32y exists

### Community 203 - "GPU Cache Reader Fix"
Cohesion: 0.29
Nodes (6): Correct format, Correct scaling, Evidence, GPU Cache Reader Fix, Problem, Relevant code

### Community 204 - "DevicePtr"
Cohesion: 0.43
Nodes (4): T, DevicePtr, n, p

### Community 205 - "F3ScopedTimer"
Cohesion: 0.29
Nodes (6): count, time_point, F3ScopedTimer, debug, name, t0

### Community 206 - "F3ScanRunFiles"
Cohesion: 0.29
Nodes (7): F3ScanRunFiles, Ecm_final, Ecm_initial, Ecm_points, normalized_output_filename, raw_output_filename, scan_tag

### Community 207 - "WorkBracket"
Cohesion: 0.29
Nodes (7): WorkBracket, BL, BR, depth, init_id, initL, initR

### Community 208 - "L20 user-validated pass-only fit final report"
Cohesion: 0.29
Nodes (6): All fitted-model true zeros, Blocked next cutoff, L20 user-validated pass-only fit final report, Required artifacts, Results, Scope

### Community 209 - "audit_l20_chi2_reconstruction.py"
Cohesion: 0.52
Nodes (6): fit_dir(), main(), parse_matrix(), ndarray, Path, read_csv()

### Community 210 - "DevicePtr"
Cohesion: 0.43
Nodes (4): T, DevicePtr, n, p

### Community 211 - "fast.py"
Cohesion: 0.16
Nodes (15): cache_info(), clear_caches(), dimension_jumps(), evaluate_point(), Fast backend, cache-free, single core (prompts/001 step 6: rungs 1-4).…, {"projector_entries", "projector_hits", "projector_misses",…, Sorted E_cm strictly inside (E_min, E_max) where the basis dimension N changes…, Drop the fast path's in-memory caches (projector V_sel; threshold tables from… (+7 more)

### Community 213 - "WindowGridRow"
Cohesion: 0.29
Nodes (7): WindowGridRow, det_imag, det_real, Ecm, Nfull, Nproj, row

### Community 214 - "compare_K3_cpp_vs_python.py"
Cohesion: 0.52
Nodes (6): add_qc3_paths(), expected_element_from_python(), local_index_from_ell_m(), main(), read_cpp_configs(), read_cpp_matrix()

### Community 215 - "CLAUDE.md — project rules for 3body-qc-v4"
Cohesion: 0.22
Nodes (8): CLAUDE.md — project rules for 3body-qc-v4, Decisions already made by Digonto (2026-10-21/22; details in `notes/STATUS.md`), Orientation, Physics invariants (never change silently), Pushing to GitHub (pre-authorised by Digonto, 2026-10-21), Setting up a fresh (cloud) checkout, Stop and ask Digonto, Things that are known to be wrong or unverified (see `docs/OPEN_ISSUES.md`)

### Community 216 - "LocalScopedTimer"
Cohesion: 0.33
Nodes (5): time_point, LocalScopedTimer, debug, name, t0

### Community 217 - "RefineTask"
Cohesion: 0.33
Nodes (6): RefineTask, branch, candidate, E_left, E_right, round

### Community 218 - "EnergyShapeKey"
Cohesion: 0.33
Nodes (5): EnergyShapeKey, dim1, dim2, total_dim, vdim

### Community 219 - "Independent L20 chi-square reconstruction"
Cohesion: 0.33
Nodes (5): fixed_K3B0_K3E0, cutoff 0.315, fixed_K3B0_K3E0, cutoff 0.325, float_from_K3B0_K3E0, cutoff 0.315, float_from_K3B0_K3E0, cutoff 0.325, Independent L20 chi-square reconstruction

### Community 220 - "L20 pass-only completed checkpoint"
Cohesion: 0.33
Nodes (5): Fit results, Key artifacts, L20 pass-only completed checkpoint, Next action, Scope

### Community 221 - "GitHub checkpoint push failure"
Cohesion: 0.33
Nodes (5): Failure, GitHub checkpoint push failure, Local checkpoint, Recommended manual resolution, Safety result

### Community 222 - "L20_LATTICE_VS_MODEL_ALL_CUTOFFS.md"
Cohesion: 0.33
Nodes (5): L20 lattice versus model across cutoffs, Lattice versus model, Ecm_cutoff=0.315, Lattice versus model, Ecm_cutoff=0.325, Lattice versus model, Ecm_cutoff=0.335, Lattice versus model, Ecm_cutoff=0.345

### Community 223 - "L20_PARAMETER_COVARIANCE_COMPARISON.md"
Cohesion: 0.33
Nodes (5): L20 covariance and error comparison, L20 fit summary, Ecm_cutoff=0.315, L20 fit summary, Ecm_cutoff=0.325, L20 fit summary, Ecm_cutoff=0.335, L20 fit summary, Ecm_cutoff=0.345

### Community 224 - "main"
Cohesion: 0.67
Nodes (5): ltag(), main(), Path, read_summary(), suffix_base()

### Community 225 - "build_l20_user_validated_windows.py"
Cohesion: 0.60
Nodes (5): lattice_levels(), main(), Path, read_label_rows(), sha256()

### Community 226 - "ConfigEntry"
Cohesion: 0.22
Nodes (9): ConfigEntry, ell, m, nx, ny, nz, px, py (+1 more)

### Community 227 - "ScopedTimer"
Cohesion: 0.33
Nodes (5): time_point, ScopedTimer, debug, name, t0

### Community 228 - "EnergyShapeKey"
Cohesion: 0.33
Nodes (5): EnergyShapeKey, dim1, dim2, total_dim, vdim

### Community 229 - "V32zuF2K2ValidationStats"
Cohesion: 0.33
Nodes (6): V32zuF2K2ValidationStats, chunks, max_abs_diff_F2, max_abs_diff_K2inv, max_rel_diff_F2, max_rel_diff_K2inv

### Community 230 - "000 — Consolidation and roadmap (starter plan from chat, 2026-10-21)"
Cohesion: 0.11
Nodes (16): 000 — Consolidation and roadmap (starter plan from chat, 2026-10-21), Goal, Milestone 0 — consolidation (DONE, 2026-10-21), Milestone 1 — trustworthy baseline (decisions answered 2026-10-21; step 6 replaced, see the update above), Milestone 2 — root finding without human labels ("gradient change"), Milestone 3 — L24 and combined fit, Milestone 4 — speed and robustness, Milestone 5 — physics output (+8 more)

### Community 231 - "RefineTask"
Cohesion: 0.33
Nodes (6): RefineTask, branch, candidate, E_left, E_right, round

### Community 232 - "Tools quick guide"
Cohesion: 0.33
Nodes (5): Cache inspectors, Continuous QC state-line maker, Plotting scripts, Tools quick guide, v32x all-irrep cache readiness checker

### Community 233 - "Hot-window FCN checkpoint"
Cohesion: 0.40
Nodes (4): Correctness fix, FCN safety, Hot-window FCN checkpoint, Next step

### Community 234 - "Where the scatter parameters are used"
Cohesion: 0.40
Nodes (4): 1. Config parsing, 2. Stored in PhysicsParams, 3. Used in K2 inverse inside F3/F3inv construction, Where the scatter parameters are used

### Community 235 - "GPU determinant backend audit"
Cohesion: 0.40
Nodes (4): GPU determinant backend audit, Required next patch, Result, Safety decision

### Community 236 - "run_v32w_QC_fitter_norm_refine_v2_fixedcoarse.sh"
Cohesion: 0.40
Nodes (4): EIGEN_INC, MINUIT2_ROOT, OMP_NUM_THREADS, run_v32w_QC_fitter_norm_refine_v2_fixedcoarse.sh script

### Community 237 - "run_gpu_component_cachegen.sh"
Cohesion: 0.60
Nodes (3): run_gpu_component_cachegen.sh script, tagL(), usage()

### Community 238 - "plot_v33p_classifier_det_family_interactive.py"
Cohesion: 0.60
Nodes (4): energy(), main(), Path, rows()

### Community 239 - "sector"
Cohesion: 0.19
Nodes (18): accepted_windows(), Rows of data/windows/L20/L20_<irrep>/accepted_windows.csv as dicts (strings)., Production sector for the tests, with the p-wave K2^{-1} form PINNED to…, sector(), slow, Every basis the production range can produce: the fast path's basis equals the…, test_fast_basis_all_signatures_vs_reference(), test_fast_speedup_single_thread_largest_point_ge_150x() (+10 more)

### Community 240 - "module.cpp"
Cohesion: 0.22
Nodes (14): array_t, _core, dict, m, b, complex, MatrixXcd, string (+6 more)

### Community 241 - "Options"
Cohesion: 0.40
Nodes (5): Options, chop_tol, debug, eig_tol, threads_per_block

### Community 242 - "ROOT FINDING for the three-body QC"
Cohesion: 0.12
Nodes (12): HISTORY — where this code came from, Lessons that must not be relearned, Phases (condensed), Source packages, 1. Why it is hard (measured), 2. What was tried (oldest → newest), 3.1 How to compute it cheaply with the existing caches, 3.2 Smooth root position for Minuit (proposal) (+4 more)

### Community 243 - "v33p_gpu_det_lu.cu"
Cohesion: 0.33
Nodes (12): cublasStatus_t, cudaError_t, size_t, string, uint64_t, cublas_error(), cublas_ok(), cuda_error() (+4 more)

### Community 244 - "fgk.py"
Cohesion: 0.08
Nodes (44): kept_momenta(), Basis rule (docs/METHODS.md section 2; prompts/001 section 2.1): spectator…, Integer vectors kept for the given channel, in loop order., cube_points(), F_block(), F_matrix(), G_elem(), G_matrix() (+36 more)

### Community 245 - "Test report 001: M1 step 1 "Build + stubs""
Cohesion: 0.29
Nodes (6): Commands run (repo root, `PYTHONPATH=src/python`), Failure breakdown of check 4 (33 tests, all expected at this step), Observations (not failures), Results, Test report 001: M1 step 1 "Build + stubs", Verdict: PASS

### Community 246 - "load_labels"
Cohesion: 0.83
Nodes (3): load_labels(), main(), Path

### Community 247 - "make_spectrum_config_from_fit.py"
Cohesion: 0.83
Nodes (3): main(), read_summary(), replace_key()

### Community 248 - "main"
Cohesion: 0.83
Nodes (3): main(), Path, write_report()

### Community 249 - "materialize_v33p_hot_windows.py"
Cohesion: 1.00
Nodes (3): main(), read_rows(), read_targets()

### Community 250 - "main"
Cohesion: 0.83
Nodes (3): main(), Path, records()

### Community 251 - "SvdMinCandidateV30r"
Cohesion: 0.22
Nodes (9): SvdMinCandidateV30r, coarse_i, E_left, E_mid, E_right, proj_dim, sigma_left, sigma_mid (+1 more)

### Community 253 - "vector"
Cohesion: 0.36
Nodes (6): IrrepCache, pair, vector, print_cache_summary(), print_split_cache_summary(), split_fixed_coarse_and_refined()

### Community 254 - "QCSearchTiming"
Cohesion: 0.50
Nodes (4): QCSearchTiming, candidate_generation_sec, classifier_sec, determinant_scan_sec

### Community 258 - "load_runtime_cache"
Cohesion: 0.18
Nodes (16): comp, ifstream, int32_t, IrrepCache, istream, MatrixXcd, ProjectedQCCacheEntry, uint64_t (+8 more)

### Community 259 - "Options"
Cohesion: 0.67
Nodes (3): Options, debug, threads_per_block

### Community 260 - "test_qc_lattice_roots.py"
Cohesion: 0.09
Nodes (36): chi2(), elab_to_ecm(), jackknife_mean_cov(), jackknife_resample(), _lattice_dir(), load_levels(), Lattice .jack reader and level matching (prompts/001 step 4, 5). Decision D1: a…, chi^2 = r^T C^+ r with the eigenvalue pseudo-inverse dropping eigenvalues <=… (+28 more)

### Community 307 - "mutex"
Cohesion: 0.19
Nodes (14): fast_projector_cached(), FastProjector, method, Nproj, Vsel, Mat, mutex, string (+6 more)

### Community 308 - "DeviceConfig"
Cohesion: 0.12
Nodes (13): DeviceConfig, ell, m, n, nx, ny, nz, DeviceProjectionTables (+5 more)

### Community 310 - "real_wigner_d.hpp"
Cohesion: 0.20
Nodes (28): wigner_d_tests(), approx_equal(), approx_equal_matrix(), check_m_range(), D_real_element(), D_real_matrix(), D_real_matrix_proper(), factorial_int() (+20 more)

### Community 311 - "test_m1_fast_v2.py"
Cohesion: 0.15
Nodes (19): _fast(), _jumps_rule(), _proj(), M1 step 7 (prompts/002): fast path v2 — fast-owned projector cache with a…, The 5.4 s case of report 006 (000_A1m 0.34538, N=103, |LG|=48). Cold projector…, Report 006 note 1: the NumPy tail used ~4 cores (CPU/wall 3.9). With the tail…, Fast-path call; the projector method is always passed explicitly (default…, E_lab thresholds E_min^+(n) of docs/METHODS.md section 2 for both channels,… (+11 more)

### Community 312 - "NeumaierComplexSum"
Cohesion: 0.33
Nodes (5): NeumaierComplexSum, c_im, c_re, sum_im, sum_re

### Community 316 - "F — Development history v31 → v32y: lessons, bugs fixed, dead ends"
Cohesion: 0.10
Nodes (20): 0. Glossary (read first), 1. What was asked, 2. What I did, 3. Pre-v31 background (from the v31u context-transfer package), 4. Chronological table v31v → v32y (Question 1), 6. Gradient / slope / derivative use (Question 3), 7.1 Recorded in the archive at the time, 7.2 Noticed during this review (**[my observations]**, need confirmation by running) (+12 more)

### Community 317 - "anchor_points"
Cohesion: 0.14
Nodes (30): anchor_points(), E_lab_of(), fast_point(), grid_energy(), _key(), load_anchors(), oracle_point(), Shared helpers for the M1 acceptance tests (planner-owned, locked). Nothing… (+22 more)

### Community 318 - "dig_tools.hpp"
Cohesion: 0.14
Nodes (11): Func, string, T, time_point, ManualTimer, t0, print_function_output_to_file(), printer() (+3 more)

### Community 319 - "4. Q2 - Exactly how the fast FCN works (L20 code)"
Cohesion: 0.20
Nodes (10): 4.1 The quantization condition as implemented, 4.2 What is precomputed and cached (three cache levels), 4.3 Startup (outside the FCN), 4.4 Accepted windows (how each lattice level gets its search window), 4.5 Inside one FCN call, 4.6 Chi-square and covariance, 4.7 Minuit2 setup and parameter mask, 4.8 Run modes of the binary `bin/v33f_k3df_fitter_multiL_v33e <config> <mode>` (+2 more)

### Community 320 - "RepDiagnostics"
Cohesion: 0.20
Nodes (10): RepDiagnostics, best_convention_index, closure_AB_max, closure_BA_max, closure_best_max, group_size, projector_diff_old_rel, projector_herm_res (+2 more)

### Community 321 - "ZeroCandidate"
Cohesion: 0.33
Nodes (6): ZeroCandidate, a, branch_index, coarse_i, kind, method

### Community 322 - "Sector"
Cohesion: 0.11
Nodes (26): parametrize, Production constants (CLAUDE.md invariant 3; prompts/001 section 5.1)., ndarray, Sector definition (prompts/001 section 5.2)., Production sector (CLAUDE.md invariant 3); `overrides` replace individual…, Sector, test_production_default_k2_pwave_form_is_q0_consistent(), Unit tests for the fast path v1 (prompts/001 step 6). Small points only. (+18 more)

### Community 323 - "/kickoff — start the physics loop here"
Cohesion: 0.22
Nodes (8): 1. Scaffold (idempotent — never overwrite existing files), 2. Starter plan → prompts/000, 3. Code graph, 4. First commit, 5. Planner (kickoff mode), 6. Stubs and lock, 7. STATUS and go, /kickoff — start the physics loop here

### Community 324 - "InwardShapeResult"
Cohesion: 0.25
Nodes (8): InwardShapeResult, classification, decreasing_fraction, ezero, flip_left_index, flip_right_index, increasing_fraction, trend_score

### Community 325 - "test_m1_reference.py"
Cohesion: 0.21
Nodes (16): as_c(), basis_list(), ||a-b||_F / ||b||_F (b is the reference)., PointResult.basis (structured) -> [(flavour, pvec, ell, m)] for qc3ref., rel_frob(), basis_to_list(), Structured basis array -> list of (flavour, pvec(3,), ell, m) with pvec = 2 pi…, test_oracle_reality_hermiticity_isometry() (+8 more)

### Community 326 - "Options"
Cohesion: 0.40
Nodes (5): Options, cuda_threads_per_block, debug, max_energies_per_chunk, vram_safety_fraction

### Community 327 - "k2_pwave_shift.py"
Cohesion: 0.25
Nodes (9): _ev(), find_root(), main(), M1 step 9: p-wave K2^{-1} form quantification for decision D7 (prompts/001 sec…, Plan rule: refine_root_info (default kind="zero": poles raise PoleError) in…, Supplementary (not the plan rule): qc3.roots.find_sign_changes on npts points;…, read_windows(), scan_roots() (+1 more)

### Community 328 - "fixer.md"
Cohesion: 0.29
Nodes (6): If you cannot fix it, Method, Not allowed, Report (always), Return to the caller (short), Start

### Community 329 - "tester.md"
Cohesion: 0.29
Nodes (6): Hard rules, Modes (the caller says which; default = step), Report (always), Return to the caller (short), Start, What "verified" means (check every item that applies)

### Community 330 - "guard.py"
Cohesion: 0.57
Nodes (6): check_bash(), check_write(), deny(), is_acc(), main(), segments()

### Community 331 - "Tester report 005: M1 step 5, portability + three no-silent-failure fixes"
Cohesion: 0.25
Nodes (7): Commands run, Failures, Informational (not a failure): legacy files outside the M1 path still containing the literals, Notes (non-blocking), Results, Tester report 005: M1 step 5, portability + three no-silent-failure fixes, Verdict: PASS (2 non-blocking notes)

### Community 333 - "string"
Cohesion: 0.15
Nodes (14): cache_key_v25(), evaluate_one_energy_cpu_eigen_v25_QC_cached(), F3iProjectedCache, hits, misses, table, size_t, string (+6 more)

### Community 334 - "Research 008 — M1 step 7 re-plan: cold projector, basis builder, $\mathcal K_3$/$G$ kernels, tail, and pole/zero discrimination"
Cohesion: 0.20
Nodes (9): 1. Where the time goes now (report 006, 1 core), 2. Projector from orbits (candidate (a)) — method and evidence, 3. Basis builder (candidate (b)), 4. $\mathcal K_{\rm df,3}$ and $G$ kernels (candidates (c), (d)), 5. Tail (candidate (e)), 6. $F$ (optional), 7. Pole/zero discrimination (coordinator addendum, step 9 findings), 8. Ranking (gain on realistic use = a scan or a fit evaluation over a window) (+1 more)

### Community 335 - "test_m1_lattice.py"
Cohesion: 0.16
Nodes (14): match_levels(), D3: sort both lists ascending; pair the lowest m QC levels with the m lattice…, _P_of(), M1 step 4 — lattice data per decision D1 (per-configuration values: resample,…, Legacy (trusted v32w) order: resample E_lab, convert each resample to E*, then…, Counts recomputed from the files; must also equal the table of…, test_covariance_symmetric_psd_and_correlated(), test_level_counts_per_cutoff_L20() (+6 more)

### Community 336 - "auditor.md"
Cohesion: 0.40
Nodes (4): Checklist (state each item: OK / ISSUE / N/A, with file:line), Report (always), Return to the caller (short), Start

### Community 337 - "builder.md"
Cohesion: 0.40
Nodes (4): Build, Not allowed (a guard hook enforces most of this), Return to the coordinator (short), Start

### Community 338 - "EnergyMatrixPack"
Cohesion: 0.14
Nodes (15): cuDoubleComplex, MatrixXcd, EnergyMatrixPack, Ecm, En, F2, G, H (+7 more)

### Community 340 - "session_start.py"
Cohesion: 0.83
Nodes (3): latest(), main(), sh()

### Community 341 - "BasisRow"
Cohesion: 0.14
Nodes (15): BasisRow, B0, B1, BB, BE, Ecm, Fproj, proj_dim (+7 more)

### Community 344 - "test_K3mat_2plus1_cpp_dump.cpp"
Cohesion: 0.39
Nodes (7): comp, MatrixXcd, string, vector, main(), write_config_file(), write_matrix_file()

### Community 345 - "Report 007: M1 step 9, p-wave $\mathcal K_2^{-1}$ level-shift table (D7)"
Cohesion: 0.29
Nodes (6): Commands, Findings (not step-9 failures), Report 007: M1 step 9, p-wave $\mathcal K_2^{-1}$ level-shift table (D7), Results, Sign argument (check 5), Verdict

### Community 346 - "DeviceConfigView"
Cohesion: 0.29
Nodes (7): DeviceConfigView, ell, m, n, nx, ny, nz

### Community 347 - "Tester report 006: M1 step 6, fast path v1 (rungs 1-4)"
Cohesion: 0.20
Nodes (9): Check 3: step-6 test ids (pinned to one core, `OMP_NUM_THREADS=1`), Check 6: 12 anchors, $\vec c$=baseline, 1 thread, pinned (`results/m1_step6_timings_1thread.json`), Commands run, Failures, Notes (non-blocking; none changes a physics number), Results, Tester report 006: M1 step 6, fast path v1 (rungs 1-4), Verdict line (+1 more)

### Community 349 - "default_data_dir"
Cohesion: 0.38
Nodes (6): default_data_dir(), Data root: $QC_DATA_ROOT if set (not checked for existence), else <repo>/data.…, _files(), M1 step 5 — no machine-specific absolute paths on the M1 production path;…, test_data_root_env_var(), test_no_absolute_home_paths_in_m1_code()

### Community 350 - "7. Q4 - Fit results achieved"
Cohesion: 0.33
Nodes (6): 7.1 v32w (trusted single-L baseline, 15 Jun), 7.2 v33n cutoff scan (L20+L24, 9 Jul; computed with `raw_sign_only`, v33g runtime cache, full scan), 7.3 L20 pass-only projected-basis hot-window fits (15 Jul), 7.4 K3B/K3E zero tests (`L20/output/l20_k3b_k3e_zero_tests/L20_K3B_K3E_ZERO_TEST_FINAL_REPORT.md:7-61`), 7.5 Blockers, 7. Q4 - Fit results achieved

### Community 351 - "9. Q6 - Recommended minimal fitter file set for the consolidated repo"
Cohesion: 0.33
Nodes (6): 9.1 C++ sources (fitter executable `v33f_k3df_fitter_multiL_v33e`), 9.2 Python scripts (under `L20/scripts/`), 9.3 Configs and data, 9.4 Not needed in the consolidated fitter, 9.5 Suggested consolidation design for a fast FCN (based on what works), 9. Q6 - Recommended minimal fitter file set for the consolidated repo

### Community 353 - "BenchmarkTiming"
Cohesion: 0.17
Nodes (12): BenchmarkTiming, cache_load_sec, chi2, chisq_sec, fallback_occurred, full_scan_occurred, model_found, precompute_sec (+4 more)

### Community 354 - "string"
Cohesion: 0.23
Nodes (13): adaptive_final_zeros_device_resident(), choose_projection_parity_v20(), array, string, time_point, get_energy_shape_prepass(), main(), process_same_shape_group_stream_parallel() (+5 more)

### Community 356 - "DetRow"
Cohesion: 0.25
Nodes (9): classify_inward_shape(), classify_inward_shape_result(), DetRow, det_norm, det_raw, Ecm, i, linear_zero() (+1 more)

### Community 357 - "load_binary_f3inv_cache_v32f"
Cohesion: 0.71
Nodes (8): istream, load_binary_f3inv_cache_v32f(), read_binary_raw_v32f(), read_comp_v32f(), read_matrix_v32f(), read_nested_vec_comp_v32f(), read_string_v32f(), read_vec_comp_v32f()

### Community 359 - "F2GpuOptions"
Cohesion: 0.40
Nodes (5): F2GpuOptions, debug, erfi_max_terms, erfi_tol, threads_per_block

## Knowledge Gaps
- **3350 isolated node(s):** `Direction for the cloud session (D4, Digonto's words, paraphrased)`, `Decisions answered at the cloud session (2026-10-22)`, `Open failures and attempt counters`, `Notes for the M1 auditor (from tester 002)`, `Known-unverified` (+3345 more)
  These have ≤1 connection - possible missing edges or undocumented components. (Counts symbols only; 3723 node(s) total have ≤1 connection when file, concept and rationale nodes are included.)
- **60 thin communities (<3 nodes) omitted from report** — run `graphify query` to explore isolated nodes.

## Suggested Questions
_Questions this graph is uniquely positioned to answer:_

- **Why does `map` connect `map` to `F3_cpu_openmp_v24_K3_QC_multi_signflip_side_abs_refinement.cpp`, `generate_swave_cpu_cache_qc.cpp`, `generate_cpu_component_cache.cpp`, `scan_projF3inv_digonto_classifier_v2.cpp`, `string`, `MultiConfig`, `qc_fitter_norm_refine_v2_multiL.cpp`, `v32w_single_L_fitter/K3df_minuit_fit_v32f_fullF3inv_QCfull_cached_classifier.hpp`, `run`, `v33g_build_runtime_k3basis_cache.cpp`, `run_block`, `F3_gpu_omp_cublas_pipeline_v9_inward_shape_refinement.cu`, `v33g_runtime_k3basis_cache.hpp`, `mutex`, `fitter/K3df_minuit_fit_v32f_fullF3inv_QCfull_cached_classifier.hpp`, `vector`, `vector`, `Output`, `run_groups_concurrent`, `Kin`, `FitResult`, `parse_label`, `string`, `F3_device_resident_true_batch_v21_side_abs_trend_refine.cu`, `F3_device_resident_true_batch_v16_custom_batched_getrs.cu`, `F3_cpu_openmp_v25_K3QC_cached_core.hpp`, `QCRefinedFCN`, `fitter/qc_fitter_norm_refine_v2.cpp`, `FitResult`, `debug_v30r_svd_validated_cartesian_l1_100_A2_zero_scan.cpp`, `v32w_single_L_fitter/qc_fitter_norm_refine_v2.cpp`, `vector`?**
  _High betweenness centrality (0.084) - this node is a cross-community bridge._
- **Why does `Options` connect `Options` to `string`, `debug_v30r_svd_validated_cartesian_l1_100_A2_zero_scan.cpp`, `ZeroRecord`?**
  _High betweenness centrality (0.018) - this node is a cross-community bridge._
- **Why does `PipelineGpuCache` connect `PipelineGpuCache` to `EnergyMatrixPack`, `CpuGpuChunk`, `F3_gpu_omp_cublas_pipeline_v9_inward_shape_refinement.cu`?**
  _High betweenness centrality (0.017) - this node is a cross-community bridge._
- **What connects `Direction for the cloud session (D4, Digonto's words, paraphrased)`, `Decisions answered at the cloud session (2026-10-22)`, `Open failures and attempt counters` to the rest of the system?**
  _3350 weakly-connected nodes found - possible documentation gaps or missing edges._
- **Should `Code and function index` be split into smaller, more focused modules?**
  _Cohesion score 0.008849557522123894 - nodes in this community are weakly interconnected._
- **Should `Package `corrected_reader`` be split into smaller, more focused modules?**
  _Cohesion score 0.011764705882352941 - nodes in this community are weakly interconnected._
- **Should `Package `production_fitter`` be split into smaller, more focused modules?**
  _Cohesion score 0.011834319526627219 - nodes in this community are weakly interconnected._
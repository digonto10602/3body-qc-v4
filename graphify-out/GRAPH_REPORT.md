# Graph Report - 3body-qc-v4  (2026-10-27)

## Corpus Check
- 432 files · ~661,961 words
- Verdict: corpus is large enough that graph structure adds value.

## Summary
- 7188 nodes · 13085 edges · 400 communities (332 shown, 64 thin omitted)
- Extraction: 95% EXTRACTED · 5% INFERRED · 0% AMBIGUOUS · INFERRED: 638 edges (avg confidence: 0.86)
- Token cost: 0 input · 0 output

## Graph Freshness
- Built from commit: `c0e7a91c`
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
- process_same_shape_group_true_batch_v15
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
- F3_device_resident_true_batch_v16_custom_batched_getrs.cu
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
- K3Mom
- F3_gpu_cachegen_v32zu.cu
- PhysicsParams
- vector
- Input
- Report B — F3inv/Vsel cache generation (CPU + GPU) and root-finding lineage for the projected QC determinant
- K2_functions_gpu_safe.cuh
- SvdMinRecordV30r
- TargetLevel
- v32zt_generate_one_irrep
- Version history summary
- KKpi v33p Cache-Anchored Classifier/Fitter Checkpoint
- vector
- run_groups_concurrent
- QCPointValue
- FitResult
- ProjectedQCCacheEntry
- inspect_v32_f3inv_cache_v32f_fixed.py
- MultiTarget
- vector
- make_qc_state_lines_from_scattered_spectrum.py
- BlockSummary
- Cand
- RuntimeBlock
- `corrected_reader`
- PhysicsParams
- ProjectedQCCacheEntry
- functions_gpu_config_maker4.cuh
- F3_cpu_openmp_v25_K3QC_cached_core.hpp
- QCRefinedFCN
- gpu_cache_reader_convention.hpp
- EnergyGpuResult
- sector
- m1_helpers.py
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
- debug_v30r_svd_validated_cartesian_l1_100_A2_zero_scan.cpp
- EvalResult
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
- string
- Faddeeva.cc
- v33d CPU/GPU Component Cache Generator Package
- acceptance-lock
- EvalPoint
- fast_point_debug_scratch.hpp
- ZeroCandidate
- run_v33m_classifier_speed_search.py
- RowResult
- process_same_shape_group_true_batch_v15
- 3. Findings
- next-num
- F3iProjectedCacheEntry
- GroupStream
- Full L20 fit and validation report
- KKpi v33o Production Fitter
- PhysicsParams
- DeviceBuffer
- PhysicsParams
- slogdet
- 5. Key physics and numerics lessons (Question 2)
- CombinedRawGpuCacheRowHeader
- V32ztCacheRowPayload
- V32ztTiming
- string
- ThreadSafeQueue
- OPEN ISSUES — read before trusting any number
- ValueError
- K3iso_maker.py
- scan_projF3inv_digonto_classifier_v2.cpp
- IrrepResult
- TargetLevel
- 003 — M2: label-free root finding with D8, root tracking for the fit, and the L20 $\mathcal K_{\rm df,3}$ fit
- V32zuRunOptions
- bench.py
- 002 — M1 step 7 re-planned: cold projector from orbits, threshold basis builder, pair kernels for $\mathcal K_{\rm df,3}$ and $G$, C++ tail, and pole/zero discrimination in the root refiner
- test_tracking.py
- ZeroPole
- tracking.py
- test_m1_fast_v2.py
- 001 — M1: a cache-free, CPU-only, single-energy three-body QC with an oracle and an independent reference
- EvalResult
- QCPointValue
- ZeroPole
- PIPELINE — from physics parameters to a $\mathcal K_{\rm df,3}$ fit
- KKpi L20 user-validated hot-window fitter
- KKpi v33o production fitter
- indep.py
- make_refined_mesh_v32s
- JumpRow
- Output
- Tester report 003: M1 step 3, reference checks and the `q0_consistent` p-wave form
- evaluate_point
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
- 009 — M2 research: label-free root finding, free-level (D8) detection, root tracking for the fit, and the minimiser
- Tester report 002: M1 step 2, oracle wrapper
- L20 K3B/K3E zero tests and chi-square audit
- M2 step 5 tester report: qc3.fit.QCChi2, make_trackers, lattice.cov_eig_pinv, tracker fallback
- Test report 015: M2 step 2, `qc3.rootscan.scan_sector`
- save_binary_f3inv_cache_v32f
- DetInfo
- config_gpu_component_cachegen_sep_f2k2.sh
- Final shutdown handoff: L20 user-validated hot-window fitter
- Classifier Restart Notes
- Function Map
- Core function inventory with descriptions
- DetInfo
- DeviceConfig
- projections_gpu_v32y_exact_device.cuh
- STATUS
- main
- Kin
- fast_point.hpp
- sector
- ComponentMatrixCacheHeader
- K2inv_2plus1_by_form
- ClassifierDispatchInfo
- test_m2_step6_fit.py
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
- Tester report 018: M2 step 6 (fit_cutoff / verify_fit / cutoff_scan / save_results / CLI / code equivalence)
- ofstream
- mutex
- fit.py
- CLAUDE.md — project rules for 3body-qc-v4
- DeviceConfig
- RefineTask
- EnergyShapeKey
- Independent L20 chi-square reconstruction
- L20 pass-only completed checkpoint
- GitHub checkpoint push failure
- L20_LATTICE_VS_MODEL_ALL_CUTOFFS.md
- L20_PARAMETER_COVARIANCE_COMPARISON.md
- main
- build_l20_user_validated_windows.py
- fast_fclass.hpp
- read_gpu_cache_complex
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
- 014 M2 step 1: interface stubs + `qc3.free_levels` (tester, mode step)
- module.cpp
- Options
- HISTORY.md
- v33p_gpu_det_lu.cu
- fgk.py
- Test report 001: M1 step 1 "Build + stubs"
- load_labels
- make_spectrum_config_from_fit.py
- main
- materialize_v33p_hot_windows.py
- main
- Mom
- FCNBase
- 3. Results
- QCSearchTiming
- L20 full Ecm-cutoff fit suite
- L20 user-validated prefit coverage audit
- plot_fit_levels.py
- Test report 009: M1 step 7.2 (orbit/character projector) and 7.3 (threshold basis builder, `dimension_jumps`)
- gpu_complex_read_convention_name
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
- fast_projector.hpp
- test_m2_step5_chi2.py
- src/README.md
- real_wigner_d.hpp
- DeviceResidentResult
- WindowGridRow
- nonint3body_general_doublecover/README.md
- F — Development history v31 → v32y: lessons, bugs fixed, dead ends
- anchor_points
- dig_tools.hpp
- test_rootscan.py
- RepDiagnostics
- MatrixCompareStats
- Sector
- /kickoff — start the physics loop here
- AssignmentRow
- QCRefinedFCN
- Test report 010: M1 step 7.4 ($\mathcal K_{\rm df,3}$ pair kernel), 7.5 ($G$ pair kernel), 7.6 (C++ tail)
- Tester report 012: M1 sub-step 7.7, class moment tables + shared (canonical $\mathbf n_P$, $|\mathbf p|^2$) sums for the $F$ lattice sum
- fixer.md
- tester.md
- guard.py
- Tester report 005: M1 step 5, portability + three no-silent-failure fixes
- 3body-qc-v4 — fitting the three-particle K matrix to lattice spectra
- F3iProjectedCache
- Research 008 — M1 step 7 re-plan: cold projector, basis builder, $\mathcal K_3$/$G$ kernels, tail, and pole/zero discrimination
- F3_device_resident_true_batch_v21_side_abs_trend_refine.cu
- auditor.md
- builder.md
- generate_K3iso.cpp
- unit/conftest.py
- session_start.py
- LocalScopedTimer
- chi2_covariance_convention_check.py
- point.py
- test_K3mat_2plus1_cpp_dump.cpp
- Report 007: M1 step 9, p-wave $\mathcal K_2^{-1}$ level-shift table (D7)
- 3. Q1 - Which fitter source is newest, and what changed
- Tester report 006: M1 step 6, fast path v1 (rungs 1-4)
- V3Shoulder
- 7. Q4 - Fit results achieved
- 6. Problems found during this review (read before trusting any L20/v33n chi-square)
- Audit 002: M2 milestone (label-free roots, tracking, L20 fit)
- conftest.py
- SignFlip
- EnergyPrepassInfo
- Tester report 016: M2 step 4 (qc3.tracking)
- Benchmark report 001: M1 single-energy QC, fast path vs oracle (milestone gate)
- CacheKey
- test_m1_fast_vs_oracle.py
- 5. The M2 method (plan `prompts/003`; research `reports/research/009`)
- EnergyMatrixPack
- v32w_single_L_fitter/qc_fitter_norm_refine_v2.cpp
- F3_cpu_openmp_v24_K3_QC_multi_signflip_side_abs_refinement.cpp
- scanN.py
- EnergyPrepassInfo
- string
- Tester report 011: M1 step 8, benchmark harness `qc3.bench` and the $\chi^2$ proxy
- acc.py
- Tester report 013: M1 milestone gate, mode `full`
- score_labels.py
- .upload
- slow
- InwardShapeResult
- blasprobe.py
- find_QC_zeros_refined
- rootscan.py
- Tester report 019: M2 gate, mode full
- vector
- DetInfo
- ScanPoint
- MatrixCompareStats
- default_data_dir
- 4. Q2 - Exactly how the fast FCN works (L20 code)
- DetRow
- string
- ScopedTimer
- 4. Chronological table v31v → v32y (Question 1)
- 3. Cache formats and stored objects (Question 1)
- run.sh
- LocalScopedTimer
- 9. Q6 - Recommended minimal fitter file set for the consolidated repo
- F2GpuOptions
- test_bench.py
- comp
- SignFlip
- projcheck.py
- Options
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
9. `Sector` - 63 edges
10. `EvalData` - 55 edges

## Surprising Connections (you probably didn't know these)
- `test_ambiguous_candidate_is_a_loud_event()` --uses--> `AmbiguousRootError`  [INFERRED]
  tests/unit/test_tracking.py → src/python/qc3/roots.py
- `test_nonconstant_dimension_inside_a_segment_is_loud()` --uses--> `DimensionJumpError`  [INFERRED]
  tests/unit/test_rootscan.py → src/python/qc3/roots.py
- `FakeQ` --uses--> `SingularQError`  [INFERRED]
  tests/unit/test_tracking.py → src/python/qc3/point.py
- `test_tracker_errors_propagate_and_verify_roots_calls_rescan()` --uses--> `AmbiguousRootError`  [INFERRED]
  tests/unit/test_m2_step5_chi2.py → src/python/qc3/roots.py
- `FlakyTracker` --uses--> `AmbiguousRootError`  [INFERRED]
  tests/unit/test_m2_step6_fit.py → src/python/qc3/roots.py

## Import Cycles
- None detected.

## Communities (400 total, 64 thin omitted)

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
Nodes (72): compute_diagnostics(), compute_F3iso_inv_from_F3inv_full(), Config, cache_file, coarseN, combined_output_file, const_norm_power, Ecm_max (+64 more)

### Community 7 - "G_gpu_safe_builder.cuh"
Cohesion: 0.14
Nodes (48): G, boost_gpu(), cabs(), cabs2(), cadd(), cdiv(), cdiv_d(), cexp_c() (+40 more)

### Community 8 - "Options"
Cohesion: 0.03
Nodes (71): Options, accept_only_likely_zero, alpha, atmK, atmpi, coarse_parallel, debug, dimension_jump_guard (+63 more)

### Community 9 - "generate_cpu_component_cache.cpp"
Cohesion: 0.06
Nodes (66): build_component_entry(), ComponentEntry, dim1, dim2, Ecm, En, error, F2 (+58 more)

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
Cohesion: 0.05
Nodes (41): BenchmarkTiming, cache_load_sec, chi2, chisq_sec, fallback_occurred, full_scan_occurred, model_found, precompute_sec (+33 more)

### Community 19 - "MultiConfig"
Cohesion: 0.04
Nodes (48): MultiConfig, accepted_zeros_files, auto_build_missing_coarse, base, benchmark_repeat, benchmark_warmup, build_v33g_runtime_if_missing, cache_backend (+40 more)

### Community 20 - "qc_fitter_norm_refine_v2_multiL.cpp"
Cohesion: 0.15
Nodes (44): F3inv, parse_label(), block_cnorm(), build_runtime_block(), L, coarse_path_for(), contains_label(), convert_samples_to_ecm() (+36 more)

### Community 21 - "compare_gpu_cache_vs_cpu_openmp_det_oldscale.cpp"
Cohesion: 0.08
Nodes (51): optional, Args, gpu_cache_dir, irreps, k3, k3df_zero, Lbyas, outdir (+43 more)

### Community 22 - "C. Fitter lineage and fast-fit caching design (K K-bar pi, RFT three-body QC)"
Cohesion: 0.18
Nodes (11): 0. Path abbreviations used in this report, 10. What the results mean (summary), 11. Problems and assumptions, 1. What was asked, 2. Concrete actions taken, 5.1 Classifier generations, 5.2 Why human labels were needed, 5.3 Workflow (from `L20/README.md:64-75` and scripts) (+3 more)

### Community 23 - "RowOut"
Cohesion: 0.05
Nodes (44): RowOut, const_norm_scale, critical_score, det, det_abs, det_imag, det_imag_abs, det_imag_rel (+36 more)

### Community 24 - "V31aEigPoint"
Cohesion: 0.06
Nodes (42): collect_v31a_flips(), MatrixXcd, string, vector, eigen_summary_from_projF3inv_v31a(), hermitize_v31a(), main(), v31a_point_from_eval() (+34 more)

### Community 25 - "The physics loop (active when `notes/STATUS.md` exists in the project)"
Cohesion: 0.13
Nodes (14): Build loop, Code understanding: graphify first, Coordinator rules (main session), Escalation — objective triggers, counted in STATUS, Global rules (Digonto) — loaded in every Claude Code session, graphify, Hard rules for correctness, How to report to me (every session, every project) (+6 more)

### Community 26 - "process_same_shape_group_true_batch_v15"
Cohesion: 0.49
Nodes (10): compute_one_energy_device_resident(), comp, cublasHandle_t, vector, cublas_lu_solve_inplace_batched_same_n(), determinant_colmajor_device_matrix(), device_lu_solve_inplace(), make_scatter_params_1() (+2 more)

### Community 27 - "v32w_single_L_fitter/K3df_minuit_fit_v32f_fullF3inv_QCfull_cached_classifier.hpp"
Cohesion: 0.14
Nodes (40): build_cache_entry(), build_F3inv_only_cache(), canonical_shell_momentum(), chi_square_v32f(), clean_label(), covariance_to_correlation_v32f(), det_real_projected_F3inv_v32f(), find_all_zero_poles() (+32 more)

### Community 28 - "run"
Cohesion: 0.08
Nodes (41): Args, coarseN, Emax, Emin, git_commit, irrep, L, outdir (+33 more)

### Community 29 - "projections_gpu_safe.cuh"
Cohesion: 0.11
Nodes (39): build_projection_Vsel_device(), build_projection_Vsel_device(), Vsel, __device__, all_perms_host(), build_PI_kernel(), build_projector_from_eigenvectors_near_one_gpu_real(), chi_host() (+31 more)

### Community 30 - "v33g_build_runtime_k3basis_cache.cpp"
Cohesion: 0.11
Nodes (40): BuildConfig, base, coarse_cache_root, explicit_coarse, force_rebuild, gpu_coarseN, gpu_Ecm_max, gpu_Ecm_min (+32 more)

### Community 31 - "run_block"
Cohesion: 0.10
Nodes (39): set, boolean_value(), BridgeArgs, accepted_truezeros_csv, cpu_fallback, gpu_batch_rows, progress_every, reference_csv (+31 more)

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
Cohesion: 0.11
Nodes (40): debug_F2_2plus1_component_for_lm_steps(), debug_F2_ang_mom_steps(), debug_I_int_ang_mom_steps(), debug_I_sum_ang_mom_steps(), F2_2plus1_mat(), F2_ang_mom(), F2_i_ang_mom_mat(), comp (+32 more)

### Community 36 - "CpuGpuChunk"
Cohesion: 0.16
Nodes (20): append_failed_results_for_chunk(), bytes_to_mib(), CpuGpuChunk, chunk_id, Ecm_end, Ecm_start, estimated_gpu_bytes, i_end (+12 more)

### Community 37 - "functions.h"
Cohesion: 0.19
Nodes (36): boost(), config_maker(), config_maker_1(), config_maker_2(), config_maker_3(), config_maker_4(), config_maker_4_momentum_first(), config_maker_5() (+28 more)

### Community 38 - "AGENTS.md"
Cohesion: 0.06
Nodes (34): 10. code-reviewer, 11. report-writer, 1. repo-explorer, 2. implementation-worker, 3. validation-scientist, 4. gpu-performance-engineer, 5. cpu-openmp-reference-keeper, 6. cache-format-guardian (+26 more)

### Community 39 - "PipelineGpuCache"
Cohesion: 0.06
Nodes (31): PipelineGpuCache, cap, dDet, dF2, dF2_array, dF2X1, dF2X1_array, dF3 (+23 more)

### Community 40 - "F2_gpu_safe_builder.cuh"
Cohesion: 0.28
Nodes (33): build_F2_2plus1_kernel(), cuDoubleComplex erfi_series_d(), cuDoubleComplex F2_ang_mom_d(), cuDoubleComplex I0F_d(), cuDoubleComplex I1F_d(), cuDoubleComplex I2F_d(), cuDoubleComplex I_int_ang_mom_d(), cuDoubleComplex I_sum_ang_mom_d() (+25 more)

### Community 41 - "F3_device_resident_true_batch_v16_custom_batched_getrs.cu"
Cohesion: 0.22
Nodes (18): build_F3_from_F2_and_F2X1_kernel(), c_make(), combine_H_colmajor_kernel(), ConfigEntry, cuDoubleComplex, __global__, gather_real_evecs_to_complex_Vsel_kernel(), v14_build_F2_batch_same_shape_kernel() (+10 more)

### Community 42 - "RuntimeCacheMeta"
Cohesion: 0.06
Nodes (32): size_t, vector, RuntimeCacheMeta, alpha, atmK, atmpi, coarseN, created_by (+24 more)

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
Cohesion: 0.17
Nodes (28): ifstream, int32_t, IrrepCache, map, MatrixXcd, ofstream, ProjectedQCCacheEntry, string (+20 more)

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
Nodes (28): append_zero_records(), array, ofstream, MomentumIrrepSpec, irrep, irrep_tag, label, nnP (+20 more)

### Community 52 - "fitter/K3df_minuit_fit_v32f_fullF3inv_QCfull_cached_classifier.hpp"
Cohesion: 0.20
Nodes (26): assemble_QC(), build_cache_entry(), build_F3inv_only_cache(), closest_zero_eigenvalue(), det_real_projected_F3inv_v32f(), find_zero_poles_from_cache(), hermitize(), comp (+18 more)

### Community 53 - "K3Mom"
Cohesion: 0.04
Nodes (66): g_matrix(), g_momentum(), g_pair_block(), GMom, bp, cut, E2, flav (+58 more)

### Community 54 - "F3_gpu_cachegen_v32zu.cu"
Cohesion: 0.17
Nodes (32): all_perms_of_3(), build_projector_from_eigenvectors_near_one(), cart_axis_from_m_v30m(), chi(), conj_class(), Convention, comp, cublasStatus_t (+24 more)

### Community 55 - "PhysicsParams"
Cohesion: 0.07
Nodes (25): PhysicsParams, alpha, atmK, atmpi, eig_tol, epsilon_h, eta_1, eta_2 (+17 more)

### Community 56 - "vector"
Cohesion: 0.16
Nodes (26): assign_block_model(), BlockInfo, file_irrep_hint, internal_irrep, L, target_indices, complex, ifstream (+18 more)

### Community 57 - "Input"
Cohesion: 0.06
Nodes (39): MatrixXcd, main(), rel(), compute_thresholds(), dimension_jumps(), fast_config(), fast_kinematics(), comp (+31 more)

### Community 58 - "Report B — F3inv/Vsel cache generation (CPU + GPU) and root-finding lineage for the projected QC determinant"
Cohesion: 0.11
Nodes (19): 1. What was asked, 2. Actions taken, 4.1 Pipeline per irrep (`v32zt_generate_one_irrep`, `F3_gpu_cachegen_v32zu.cu:2308-2519`), 4.2 Streams and the host-pressure guard, 4.3 Fused vs separate F2/K2, 4.4 Performance numbers found (wall-clock seconds, GPU and host unspecified unless stated), 4. GPU generator internals and performance (Question 2), 5.1 Table of methods (oldest → newest) (+11 more)

### Community 59 - "K2_functions_gpu_safe.cuh"
Cohesion: 0.23
Nodes (23): abs_cx(), build_K2inv_2plus1_kernel(), __global__, __host__, cutoff_function_1_gpu(), Cx, im, re (+15 more)

### Community 60 - "SvdMinRecordV30r"
Cohesion: 0.08
Nodes (26): collect_svd_local_min_candidates_v30r(), near_any_dimension_jump_v30r(), refine_svd_candidates_parallel_v30r(), SvdMinCandidateV30r, coarse_i, E_left, E_mid, E_right (+18 more)

### Community 61 - "TargetLevel"
Cohesion: 0.10
Nodes (26): canonical_shell_momentum(), chi_square_v32f(), covariance_to_correlation_v32f(), diagonal_chi2(), array, FunctionMinimum, MatrixD, tuple (+18 more)

### Community 62 - "v32zt_generate_one_irrep"
Cohesion: 0.19
Nodes (20): cudaDeviceProp, array, PhysicsParams, string, v32zt_energy_shape_prepass_momentum_first(), v32zt_finalize_allobj_sidecars(), v32zt_generate_one_irrep(), v32zt_hostname_string() (+12 more)

### Community 63 - "Version history summary"
Cohesion: 0.08
Nodes (24): Conversation Context Reconstruction, `F2iso` and `Giso`, Group-theory degeneracy, Important decisions made, Key definitions, Known issues and open questions, Momentum/irrep label convention, Non-interacting three-body energy (+16 more)

### Community 64 - "KKpi v33p Cache-Anchored Classifier/Fitter Checkpoint"
Cohesion: 0.08
Nodes (24): Build, Build and smoke test, Cache format notes, Checkpoint status, Critical checkpoint, Current classifier state, Exact old scaling, External cache requirement (+16 more)

### Community 65 - "vector"
Cohesion: 0.24
Nodes (17): clean_label(), find_all_zero_poles(), fit_K3df_parameters_v32f(), get_or_build_F3inv_cache_v32f(), map, string, vector, IrrepCache (+9 more)

### Community 66 - "run_groups_concurrent"
Cohesion: 0.13
Nodes (24): EnergyPrepassInfo, EnergyShapeKey, map, size_t, estimate_group_vram_bytes(), get_free_vram_bytes(), run_groups_concurrent(), v32zt_auto_batch_cap() (+16 more)

### Community 67 - "QCPointValue"
Cohesion: 0.09
Nodes (23): assemble_QC(), closest_zero_eigenvalue(), hermitize(), MatrixXcd, K3dfParameters, K3B, K3E, K3iso0 (+15 more)

### Community 68 - "FitResult"
Cohesion: 0.09
Nodes (23): FitResult, best, chi2, chi2_dof, data_correlation, data_covariance, errors, model_levels (+15 more)

### Community 69 - "ProjectedQCCacheEntry"
Cohesion: 0.09
Nodes (23): ProjectedQCCacheEntry, Ecm, En, En_c, error, F3_equiv, F3_leak, F3inv_full (+15 more)

### Community 70 - "inspect_v32_f3inv_cache_v32f_fixed.py"
Cohesion: 0.23
Nodes (22): CacheEOF, CacheFormatError, expected_coarse_grid(), file_size(), has_near(), main(), parse_cache(), Exception (+14 more)

### Community 71 - "MultiTarget"
Cohesion: 0.13
Nodes (22): Cand, build_assignment_rows(), CandidateWithBlock, c, file_irrep, internal_irrep, MatrixD, TargetLevel (+14 more)

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

### Community 76 - "RuntimeBlock"
Cohesion: 0.09
Nodes (22): BasisRow, B0, B1, BB, BE, Ecm, Fproj, proj_dim (+14 more)

### Community 77 - "`corrected_reader`"
Cohesion: 0.10
Nodes (21): `cache_generator`, `configs/`, `configs/`, `corrected_reader`, `cpu_component_cachegen/scripts/`, `cpu_component_cachegen/source/`, `diagnostics/`, `diagnostics/` (+13 more)

### Community 78 - "PhysicsParams"
Cohesion: 0.11
Nodes (18): PhysicsParams, alpha, atmK, atmpi, epsilon_h, eta_1, eta_2, K3B_par (+10 more)

### Community 79 - "ProjectedQCCacheEntry"
Cohesion: 0.10
Nodes (21): ProjectedQCCacheEntry, Ecm, En, En_c, error, F3_equiv, F3_leak, F3inv_full (+13 more)

### Community 80 - "functions_gpu_config_maker4.cuh"
Cohesion: 0.08
Nodes (46): build_lm_list_from_waves(), candidate_count_per_energy(), compact_flags_to_entries_cpu_order(), config_maker_4_gpu_safe(), ConfigEntry, ell, m, nx (+38 more)

### Community 81 - "F3_cpu_openmp_v25_K3QC_cached_core.hpp"
Cohesion: 0.31
Nodes (19): add_unique_zero(), build_F3i_projected_cache_entry_v25(), cache_key_v25(), classify_one_subflip_side_abs(), determinant_via_partial_piv_lu(), evaluate_one_energy_cpu_eigen_v24_QC(), evaluate_one_energy_cpu_eigen_v25_QC_cached(), find_likely_zeros_v24_K3QC_silent() (+11 more)

### Community 82 - "QCRefinedFCN"
Cohesion: 0.10
Nodes (21): cov_to_corr(), atomic, MatrixD, Minuit2::FCNBase, shared_ptr, TargetLevel, QCRefinedFCN, cn (+13 more)

### Community 83 - "gpu_cache_reader_convention.hpp"
Cohesion: 0.27
Nodes (19): GpuCacheFileFormat, format_from_header(), comp, GpuComplexReadConvention, istream, string, uint64_t, name() (+11 more)

### Community 84 - "EnergyGpuResult"
Cohesion: 0.10
Nodes (20): EnergyGpuResult, det_F3inv_projected, Ecm, eig_F3inv_projected, En, F2, F3, F3inv_projected (+12 more)

### Community 85 - "sector"
Cohesion: 0.06
Nodes (59): brute_force_free_levels(), chi2_bound_from_root_error(), label_path(), lambda_min_cov(), Shared helpers for the M2 acceptance tests (planner-owned, locked with the…, Non-interacting K K pi energies with total momentum d = NNP[irrep], any irrep,…, Rigorous bound on |delta chi2| when every model energy moves by at most dE…, [(irrep, level_index, lattice_Ecm, model_Ecm)] from the legacy fit-levels table… (+51 more)

### Community 86 - "m1_helpers.py"
Cohesion: 0.06
Nodes (15): ln|det Q| orbit vs eigen vs oracle (and an independent eigenvalue-based value…, Tester 010: K3 basis, G, tail vs oracle at 12 anchors + random energies (seed…, The locked test body, run in a fresh process, optionally after sleeping `delay`…, CPU/wall of 20 warm calls (the locked-test protocol) for a given point, FIRST…, Which G fill mode is closer to the exact ln|det Q| (60-digit) of ITS OWN…, 111_A2 at 0.3382011860823352 (kappa(Q)=1.5e13): who is right, fast or oracle?…, Non-production parameters (eps_H != 0, other xi, masses, eta, q0_consistent):…, fast K3 basis and G vs the planner's independent NumPy reference qc3ref (not… (+7 more)

### Community 87 - "V32ztAllObjRowPayload"
Cohesion: 0.11
Nodes (19): V32ztAllObjRowPayload, Ecm, F2, F3, F3inv, grid_i, K2inv, klm_config (+11 more)

### Community 88 - "CachedProjectorV30q"
Cohesion: 0.10
Nodes (20): CachedProjectorV30q, fv_best_convention, fv_rep_closure_best, fv_rep_unitarity, P_I_herm_res, P_I_idem_res, PI_gap_selected_to_rejected, PI_largest_rejected_eval (+12 more)

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
Cohesion: 0.15
Nodes (24): G_2plus1_mat(), G_ij_lm(), get_parity_block_from_plm_config(), comp, MatrixXcd, vector, print_Gij_boosts(), comp (+16 more)

### Community 95 - "fitter/qc_fitter_norm_refine_v2.cpp"
Cohesion: 0.17
Nodes (19): classify(), cnorm(), Eval, det, E, proj_dim, success, y (+11 more)

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

### Community 102 - "debug_v30r_svd_validated_cartesian_l1_100_A2_zero_scan.cpp"
Cohesion: 0.10
Nodes (41): eigmin_point_from_eval_v30s(), basis_signature_v30q(), BranchPoint, eigvecs, success, build_cached_projector_v30q(), coarse_eigensystems_for_method(), comp (+33 more)

### Community 103 - "EvalResult"
Cohesion: 0.17
Nodes (12): EvalResult, det, dim1, dim2, Ecm, eig_min_proj, En, error (+4 more)

### Community 104 - "find_QC_zeros_refined"
Cohesion: 0.28
Nodes (16): build_entry(), FitSettings, IrrepCache, K3dfParameters, pair, PhysicsParams, ProjectedQCCacheEntry, vector (+8 more)

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
Nodes (27): _brent(), _check_kind(), _check_xtol(), _classify(), classify_sign_change(), _Counter, DimensionJumpError, find_sign_changes() (+19 more)

### Community 112 - "METHODS — numerical methods, equations, conventions, sources"
Cohesion: 0.11
Nodes (19): 0. System and formalism, 10. Energy grid, 11. Data, covariance, χ² — `src/fitter/qc_fitter_norm_refine_v2_multiL.cpp` (`load_multil_jack_targets`), `src/fitter/lattice_data_covariance_cpp.hpp`, `chi_square_v32f`, 12. Minimiser, 13. Root finding, 14. Non-interacting levels — `tools/nonint3body_general_doublecover/`, `qc3.free_levels` (M2), 15. Cache-free single-energy evaluation: oracle, fast path, tolerances — `src/python/qc3/`, `src/qc3_core/` (plan `prompts/001`), 16. Independent reference implementation — `tests/acceptance/m1/qc3ref/` (planner-owned) (+11 more)

### Community 113 - "3. Findings"
Cohesion: 0.12
Nodes (15): 1. What was asked, 2. Actions taken, 3.1 Lattice spectrum data: what exists, 3.2 Critical: the production fitter inflates the data errors by a factor of 602, 3.3 User true-zero label CSVs and accepted windows, 3.4 Config files and keys, 3.5 Two-body input ($K_2^{-1}$) and where the scattering parameters come from, 3.6 $K_{\rm df,3}$ parametrization (exact, from `PKG/source/K3_functions_2plus1.hpp`) (+7 more)

### Community 114 - "string"
Cohesion: 0.20
Nodes (34): bisection_eigenbranch_minabs_zero(), bisection_eigenbranch_reference_zero(), bisection_f3iso_zero(), bisection_method_zero(), classify_side_abs(), collect_det_eigenbranch_candidates_from_adjacent_pairs(), collect_zero_candidates_for_method(), PhysicsParams (+26 more)

### Community 115 - "Faddeeva.cc"
Cohesion: 0.23
Nodes (12): cmplx, cpolar(), erfcx_y100(), FADDEEVA(), FADDEEVA_RE(), main(), relerr(), sinc() (+4 more)

### Community 116 - "v33d CPU/GPU Component Cache Generator Package"
Cohesion: 0.13
Nodes (14): Build and run, CPU, CPU path, Current status at packaging time, Decision tree used to get to the current generator, Execution tree, GPU, GPU path (+6 more)

### Community 118 - "EvalPoint"
Cohesion: 0.19
Nodes (15): classify_final_v32s(), pair, vector, EvalPoint, det, E, from_refinement, proj_dim (+7 more)

### Community 119 - "fast_point_debug_scratch.hpp"
Cohesion: 0.08
Nodes (33): main(), evaluate(), F_channel(), comp, M, Mat, MatR, string (+25 more)

### Community 120 - "ZeroCandidate"
Cohesion: 0.29
Nodes (7): ZeroCandidate, a, b, branch_index, coarse_i, kind, method

### Community 121 - "run_v33m_classifier_speed_search.py"
Cohesion: 0.29
Nodes (13): assignment_rows(), assignment_signature(), BenchRow, compare_assignments(), main(), parse_summary(), Path, run() (+5 more)

### Community 122 - "RowResult"
Cohesion: 0.13
Nodes (15): RowResult, abs_diff, cache_index, det_cpu, det_gpu, Ecm, irrep, Lbyas (+7 more)

### Community 123 - "process_same_shape_group_true_batch_v15"
Cohesion: 0.11
Nodes (36): adaptive_final_zeros_device_resident(), classify_inward_shape(), compute_one_energy_device_resident(), array, comp, cublasHandle_t, string, vector (+28 more)

### Community 124 - "3. Findings"
Cohesion: 0.11
Nodes (17): 001 — M1 milestone audit: cache-free single-energy QC, D1 covariance, D3 matching, D7 p-wave form, 1. What I did (commands), 2. Checklist, 3.1 Fast path vs oracle (priority 1), 3.2 p-wave $\mathcal K_2^{-1}$ (D7) (priority 2), 3.3 Statistics (priority 3), 3.4 Roots (priority 4), 3.5 Tests and tolerances (priority 5) (+9 more)

### Community 126 - "F3iProjectedCacheEntry"
Cohesion: 0.13
Nodes (15): F3iProjectedCacheEntry, dim1, dim2, Ecm, En, error, F3inv_projected, klm_config (+7 more)

### Community 127 - "GroupStream"
Cohesion: 0.15
Nodes (14): cudaStream_t, cublasHandle_t, GroupStream, batch_size, completed, dim1, dim2, estimated_vram_bytes (+6 more)

### Community 128 - "Full L20 fit and validation report"
Cohesion: 0.14
Nodes (13): Accepted-window package, Chi-square audit, Classifier validation, Errors and correlations, Exact reproducibility inputs and outputs, Full L20 fit and validation report, K3B/K3E zero test, L20 pass-only fits (+5 more)

### Community 129 - "KKpi v33o Production Fitter"
Cohesion: 0.14
Nodes (13): Build and smoke test, Cache format notes, Critical checkpoint, Current classifier state, Exact old scaling, KKpi v33o Production Fitter, Known limits, Main workflow commands (+5 more)

### Community 130 - "PhysicsParams"
Cohesion: 0.08
Nodes (25): PhysicsParams, alpha, atmK, atmpi, eig_tol, epsilon_h, eta_1, eta_2 (+17 more)

### Community 131 - "DeviceBuffer"
Cohesion: 0.10
Nodes (26): compare_F2_cpu_gpu(), complex, MatrixXcd, size_t, T, vector, DeviceBuffer, count (+18 more)

### Community 132 - "PhysicsParams"
Cohesion: 0.14
Nodes (14): PhysicsParams, alpha, atmK, atmpi, epsilon_h, eta_1, eta_2, Lbyas (+6 more)

### Community 133 - "slogdet"
Cohesion: 0.25
Nodes (12): project_Q(), QC helpers (prompts/001 step 4): projection of the inverse and (sign, ln|det|)…, (L^-6) V^dag (F3inv + K3) V ; L = xi*Lbyas (so L**6 == (Lbyas*xi)**6). The…, (sign, ln|det Q|) by partial-pivot LU (LAPACK via NumPy). sign is an int in…, slogdet(), M1 step 4 — qc3.qc implements Q = (L xi)^-6 V^dag (F3inv + K3) V and (sign,…, _synthetic(), test_lndet_invariant_under_unitary_mixing_of_Vsel() (+4 more)

### Community 134 - "5. Key physics and numerics lessons (Question 2)"
Cohesion: 0.14
Nodes (14): 5.10 Covariance fix (v32l), cache-only fix (v32k), namespace fixes (v31z, v32i, v32r, v32y), 5.11 Lattice input conventions (v32p → v32s), 5.12 Memory (v32x → v32y), 5.13 Scattering-parameter plumbing, 5.1 Use slogdet, never the raw determinant, 5.2 Projection audit: the projector was not the problem, 5.3 Pole vs zero: a sign flip is not a level, 5.4 Inverse convention: project the QC, not $F_3$ (+6 more)

### Community 135 - "CombinedRawGpuCacheRowHeader"
Cohesion: 0.14
Nodes (14): CombinedRawGpuCacheRowHeader, dim1, dim2, Ecm, En, grid_i, n, vdim (+6 more)

### Community 136 - "V32ztCacheRowPayload"
Cohesion: 0.15
Nodes (14): cuDoubleComplex, __global__, v32zt_identity_batch_kernel(), V32ztCacheRowPayload, dim1, dim2, Ecm, En (+6 more)

### Community 137 - "V32ztTiming"
Cohesion: 0.14
Nodes (14): V32ztTiming, basis_sec, bytes_written, chunks, f2k2_sec, gf3_build_sec, max_total_dim, max_vdim (+6 more)

### Community 138 - "string"
Cohesion: 0.51
Nodes (9): map, string, gd(), gi(), gs(), main(), read_kv(), settings_from_config() (+1 more)

### Community 139 - "ThreadSafeQueue"
Cohesion: 0.18
Nodes (13): condition_variable, bytes_complex_matrix(), mutex, T, estimate_gpu_bytes_one_energy(), flush_cpu_buffer_to_gpu_queue(), gpu_consumer_thread_func(), ThreadSafeQueue (+5 more)

### Community 140 - "OPEN ISSUES — read before trusting any number"
Cohesion: 0.13
Nodes (15): 10. LOW — minor physics-core bugs (not on the production path), 11. LOW — portability, 12. Historical warnings (do not reuse these results), 13. LOW — legacy projector cache key omits the selection tolerances (`basis_signature_v30q`), 14. HIGH (physics, open) — M2 L20 fit quality, error reliability and L24 (added 2026-10-26; `reports/tests/018`, `notes/STATUS.md`), 1. BLOCKER (CONFIRMED by D1) — lattice covariance too large by $(N-1)^2=602^2\approx3.6\times10^5$, 2. HIGH — p-wave $\mathcal K_2^{-1}$ differs from the reference implementation, 3. HIGH — complex-number read convention (`variant_04_real_imag_swapped`) (+7 more)

### Community 141 - "ValueError"
Cohesion: 0.21
Nodes (12): main(), Path, records(), {"tp","fp","fn","tn": int, "mismatches": list[dict], "unlabelled_zeros":…, score_labels(), add_qc3_paths(), expected_element_from_python(), local_index_from_ell_m() (+4 more)

### Community 142 - "K3iso_maker.py"
Cohesion: 0.29
Nodes (10): E_to_Ecm(), Ecmsq_to_Esq(), Esq_to_Ecmsq(), jackknifeavg_centralvalue_lattice_data(), jackknifeavg_lattice_data(), P000(), P100(), P110() (+2 more)

### Community 143 - "scan_projF3inv_digonto_classifier_v2.cpp"
Cohesion: 0.41
Nodes (12): const_norm_scale_v32s(), map, string, TargetLevel, gd(), gi(), gs(), main() (+4 more)

### Community 144 - "IrrepResult"
Cohesion: 0.15
Nodes (13): InitialFlip, El, Er, id, il, ir, yl, yr (+5 more)

### Community 145 - "TargetLevel"
Cohesion: 0.15
Nodes (13): diagonal_chi2(), TargetLevel, atP, E_read, Ecm, err, err_read, label (+5 more)

### Community 146 - "003 — M2: label-free root finding with D8, root tracking for the fit, and the L20 $\mathcal K_{\rm df,3}$ fit"
Cohesion: 0.07
Nodes (27): 003 — M2: label-free root finding with D8, root tracking for the fit, and the L20 $\mathcal K_{\rm df,3}$ fit, 0. Diagnosis — what M1 left open and why the old approach cannot continue, 10. Expected cost, 1. Goal of M2, 2.1 Free levels (D8 input) — `qc3.free_levels`, 2.2 Label-free root finder — `qc3.rootscan.scan_sector`, 2.3 Exact free-level solutions (D8) — what the literature says, 2.4 Pole/zero discrimination and conventions (unchanged from M1) (+19 more)

### Community 147 - "V32zuRunOptions"
Cohesion: 0.15
Nodes (13): V32zuRunOptions, concurrent_streams_safety_fraction, device_count, gpu_memory_safety_fraction, irrep_filter_csv, max_batch_energies, use_concurrent_streams, use_fused_f2k2 (+5 more)

### Community 148 - "bench.py"
Cohesion: 0.12
Nodes (27): build_parser(), chi2_proxy(), _cpu_model(), _git(), machine_info(), main(), _peak_rss_mb(), pin_cores() (+19 more)

### Community 149 - "002 — M1 step 7 re-planned: cold projector from orbits, threshold basis builder, pair kernels for $\mathcal K_{\rm df,3}$ and $G$, C++ tail, and pole/zero discrimination in the root refiner"
Cohesion: 0.08
Nodes (23): 002 — M1 step 7 re-planned: cold projector from orbits, threshold basis builder, pair kernels for $\mathcal K_{\rm df,3}$ and $G$, C++ tail, and pole/zero discrimination in the root refiner, 0. Diagnosis — why step 7 is re-planned, 1. Goal of step 7 (revised), 2.1 Projector from orbits (rung 7.2) — `src/qc3_core/fast_projector.hpp`, 2.2 Threshold basis builder (rung 7.3) — `src/qc3_core/fast_basis.hpp`, 2.3 Pair kernels for $\mathcal K_{\rm df,3}$ and $G$ (rungs 7.4, 7.5) — `src/qc3_core/fast_k3.hpp`, `fast_g.hpp`, 2.4 C++ tail, LU solve, lazy conversion (rung 7.6), 2.5 Pole/zero discrimination (rung 7.1) — `src/python/qc3/roots.py` (+15 more)

### Community 150 - "test_tracking.py"
Cohesion: 0.17
Nodes (16): brent(), FakeQ, _fresh(), Unit tests for qc3.tracking (M2 step 4, prompts/003 section 2.5): a fake two-…, Q = [[sg*50 (E - 0.2800), 0.05], [0.05, 30 (E - 0.2900) + 1e-4/(E - 0.2850)]] +…, test_ambiguous_candidate_is_a_loud_event(), test_corrector_ambiguous_recheck_falls_back_and_local_rescan_still_raises_for_a_genuinely_ambiguous_candidate(), test_rescan_raises_when_the_tracked_set_is_stale() (+8 more)

### Community 151 - "ZeroPole"
Cohesion: 0.17
Nodes (12): ZeroPole, E_left, E_mid, E_right, eig_index_left, eig_index_right, i_left, i_right (+4 more)

### Community 152 - "tracking.py"
Cohesion: 0.06
Nodes (31): dimension_jumps(), Sorted E_cm strictly inside (E_min, E_max) where the basis dimension N changes…, The fast path's LU of Q met an exactly zero pivot (C++ "fast: Q is singular"):…, SingularQError, ndarray, E of kind == "zero", sorted., ScanResult, _CountingEvaluate (+23 more)

### Community 153 - "test_m1_fast_v2.py"
Cohesion: 0.14
Nodes (21): _cpu_wall_ratio(), _fast(), _jumps_rule(), _proj(), M1 step 7 (prompts/002): fast path v2 — fast-owned projector cache with a…, The 5.4 s case of report 006 (000_A1m 0.34538, N=103, |LG|=48). Cold projector…, Wait until the process is idle: an idle `step` window uses < `quiet` of a core.…, Process CPU time / wall time of warm `want=("Q",)` calls only, over at least… (+13 more)

### Community 154 - "001 — M1: a cache-free, CPU-only, single-energy three-body QC with an oracle and an independent reference"
Cohesion: 0.07
Nodes (28): 001 — M1: a cache-free, CPU-only, single-energy three-body QC with an oracle and an independent reference, 0. Diagnosis (why the old pipeline is slow and what is wrong), 10. Expected cost, 1. Goal of M1, 2.1 Equations and conventions (unchanged physics; see `docs/METHODS.md` §§1–8 for the full set), 2.2 p-wave $\mathcal K_2^{-1}$: two named forms, production default unchanged until D6, 2.3 Speed-up candidates, ranked (expected gain, how verified), 2.4 Lattice data (D1, D3) (+20 more)

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
Cohesion: 0.22
Nodes (9): (A) Cache generation, (B) Determinant scan, (C) Labelling (the current gate for "which roots exist"), Combined raw cache format (`*_F3inv_Vsel_gpu.bin`), (D) Accepted windows, (E) Hot-window projected-basis cache, (F) Fit, (G) Audits and plots (+1 more)

### Community 159 - "KKpi L20 user-validated hot-window fitter"
Cohesion: 0.18
Nodes (10): Build, Current status, Included and excluded, KKpi L20 user-validated hot-window fitter, Known blockers and next action, Layout, Physics and production defaults, Plotting (+2 more)

### Community 160 - "KKpi v33o production fitter"
Cohesion: 0.18
Nodes (10): Build and smoke test, Canonical workflow, Generated outputs, How to modify the workflow, Important parameters, KKpi v33o production fitter, Known limits, Logic flow (+2 more)

### Community 161 - "indep.py"
Cohesion: 0.27
Nodes (9): go(), chi2_at(), D1_levels(), make_chi(), labels list of (irrep,n). Own numpy implementation of the D1 recipe., D3: lowest m QC zeros (kind zero only) matched to sorted lattice levels per…, zeros_fresh(), job() (+1 more)

### Community 162 - "make_refined_mesh_v32s"
Cohesion: 0.36
Nodes (11): almost_same_E(), build_refined_entry_v32s(), FitSettings, IrrepCache, PhysicsParams, ProjectedQCCacheEntry, eval_from_entry(), linear_zero() (+3 more)

### Community 163 - "JumpRow"
Cohesion: 0.18
Nodes (11): JumpRow, delta_Nfull, delta_Nproj, Ecm, jump_id, Nfull, Nproj, prev_Nfull (+3 more)

### Community 164 - "Output"
Cohesion: 0.09
Nodes (23): map, Output, basis, E_lab, Ecm, F, F3inv, G (+15 more)

### Community 165 - "Tester report 003: M1 step 3, reference checks and the `q0_consistent` p-wave form"
Cohesion: 0.18
Nodes (10): 4b. Term-by-term check of `K2_inv_ERE_ang_mom_q0_consistent` against prompts/001 §2.2, 4d. Independent implementation, 5. Max relative Frobenius deviation oracle vs `qc3ref` ($\lVert a-b\rVert_F/\lVert b\rVert_F$, reference in the denominator, $\vec c=0$), 6. Fast tally (`pytest tests/acceptance/m1 -m "not slow" -rA --tb=line`), Commands run, Failures, Notes (not failures), Results (+2 more)

### Community 166 - "evaluate_point"
Cohesion: 0.18
Nodes (21): ev(), ev(), evaluate_point(), _params(), _oracle(), parametrize, Unit tests for prompts/002 sub-steps 7.4 (K3 pair kernel), 7.5 (G pair kernel)…, prompts/002 7.6: qc3.fast must not call into qc3.qc (NumPy/BLAS) any more. (+13 more)

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
Cohesion: 0.20
Nodes (10): 3.3.1 Cutoff function $H$ (`L20/functions.h:97-141`), 3.3.2 Basis: spectator momentum × $(\ell,m)$, shell cutoff (`L20/functions.h:976-1113`), 3.3.3 $F_2$ — sum minus integral (`L20/F2_functions_v2.h`), 3.3.4 $G$ — switch-state (exchange) propagator (`L20/G_functions_v2.h`), 3.3.5 $\mathcal K_2^{-1}$ (`L20/K2_functions_v2.h`), 3.3.6 $F_3$ and $F_3^{-1}$, 3.3.7 $\mathcal K_{\rm df,3}$ parametrisation (`L20/K3_functions_2plus1.hpp`), 3.3.8 Units and constants (+2 more)

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

### Community 177 - "009 — M2 research: label-free root finding, free-level (D8) detection, root tracking for the fit, and the minimiser"
Cohesion: 0.13
Nodes (14): 009 — M2 research: label-free root finding, free-level (D8) detection, root tracking for the fit, and the minimiser, 1. What was asked, 2.1 Finding every zero of $\det Q(E)$ on $[E_{\min},E_{\max}]$ (label-free), 2.2 Exact free-level solutions (D8), 2.3 Fewer evaluations per root inside a fit, 2.4 Minimiser, 2. Candidate methods compared, 3.1 Label gate (coded_v2, $\vec c_{\rm base}$, $[0.26310,0.36]$, coarse step 4 rows $=1.938\times10^{-5}$) (+6 more)

### Community 178 - "Tester report 002: M1 step 2, oracle wrapper"
Cohesion: 0.20
Nodes (9): Check 5: review of `oracle_point.hpp` against `build_component_entry` (cachegen L223-288), Check 6 table (16 cases; full data in the JSON), Check 7: single-thread timings (`PointResult.timings`, $\vec c=\vec c_{\rm base}$, seconds), Commands run, Failures, Notes (not failures; for the auditor / planner), Results, Tester report 002: M1 step 2, oracle wrapper (+1 more)

### Community 179 - "L20 K3B/K3E zero tests and chi-square audit"
Cohesion: 0.22
Nodes (8): Chi-square normalization and FCN audit, Exact source audit points, L20 K3B/K3E zero tests and chi-square audit, Required output files, Table A, Table B, Table C, Table D

### Community 180 - "M2 step 5 tester report: qc3.fit.QCChi2, make_trackers, lattice.cov_eig_pinv, tracker fallback"
Cohesion: 0.40
Nodes (4): Commands, M2 step 5 tester report: qc3.fit.QCChi2, make_trackers, lattice.cov_eig_pinv, tracker fallback, Notes, Results

### Community 181 - "Test report 015: M2 step 2, `qc3.rootscan.scan_sector`"
Cohesion: 0.33
Nodes (5): Cost, Judgement of the builder's choices, Problems, Test report 015: M2 step 2, `qc3.rootscan.scan_sector`, Verdict: PASS

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

### Community 190 - "DeviceConfig"
Cohesion: 0.12
Nodes (14): ConfigView, ell, m, n, px, py, pz, DeviceConfig (+6 more)

### Community 191 - "projections_gpu_v32y_exact_device.cuh"
Cohesion: 0.24
Nodes (14): build_basis_arrays_exact_v32y(), build_projection_Vsel_device_exact_v32y(), build_projector_v32y_exact_kernel(), cart_axis_from_m_dev(), copy_int_vector_to_device_exact_v32y(), cuDoubleComplex, __global__, string (+6 more)

### Community 192 - "STATUS"
Cohesion: 0.20
Nodes (10): Decisions answered at the cloud session (2026-10-22), Decisions from Digonto (answered 2026-10-21), Digonto's answers for M2 (2026-10-24), Direction for the cloud session (D4, Digonto's words, paraphrased), Known-unverified, Next action, Notes for the M1 auditor (from tester 002), Open failures and attempt counters (+2 more)

### Community 193 - "main"
Cohesion: 0.61
Nodes (7): active_files(), all_config_files(), equal(), main(), parse(), Path, write_audit()

### Community 194 - "Kin"
Cohesion: 0.08
Nodes (29): Clock, t0_, compute_k3(), evaluate(), array, comp, Mat, PhysicsParams (+21 more)

### Community 195 - "fast_point.hpp"
Cohesion: 0.08
Nodes (32): evaluate(), F_channel(), comp, M, Mat, MatR, string, vector (+24 more)

### Community 196 - "sector"
Cohesion: 0.17
Nodes (19): accepted_windows(), Rows of data/windows/L20/L20_<irrep>/accepted_windows.csv as dicts (strings)., Production sector for the tests, with the p-wave K2^{-1} form PINNED to…, sector(), slow, Every basis the production range can produce: the fast path's basis equals the…, Report 006 note 2: the legacy key omitted eig_tol, so eig_tol=1e-30 after a…, test_fast_basis_all_signatures_vs_reference() (+11 more)

### Community 197 - "ComponentMatrixCacheHeader"
Cohesion: 0.25
Nodes (8): ComponentMatrixCacheHeader, coarseN, created_by, irrep, Lbyas, object_name, row_count, xi

### Community 198 - "K2inv_2plus1_by_form"
Cohesion: 0.44
Nodes (9): comp, MatrixXcd, string, vector, K2_inv_ERE_ang_mom_q0_consistent(), K2inv_2plus1_by_form(), K2inv_EREord2_2plus1_mat_q0_consistent(), K2inv_EREord2_i_mat_q0_consistent() (+1 more)

### Community 199 - "ClassifierDispatchInfo"
Cohesion: 0.25
Nodes (8): classifier_dispatch_info(), ClassifierDispatchInfo, alias, dispatch_unique, implementation_tag, requested_mode, resolved_mode, version_tag

### Community 200 - "test_m2_step6_fit.py"
Cohesion: 0.12
Nodes (16): analytic(), FlakyTracker, make_levels(), Unit tests for qc3.fit.fit_cutoff / verify_fit / save_results / the CLI helpers…, Raises AmbiguousRootError once, on its first move away from the start vector…, Weighted linear least squares of E_lat = a + J (c - C0) over the free…, setup(), test_chi2_from_roots_uses_d3_and_the_penalty() (+8 more)

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

### Community 211 - "Tester report 018: M2 step 6 (fit_cutoff / verify_fit / cutoff_scan / save_results / CLI / code equivalence)"
Cohesion: 0.18
Nodes (10): Builder deviations, Commands, Findings for the coordinator and Digonto, Fresh-scan verification (all 9 fits in the 3 JSONs), Independent chi^2 (own fresh scans, own NumPy D1 covariance from the `.jack` files, own D3; no tracker), Is it a minimum? Gradient and $\Delta\chi^2$, Provenance (JSONs), Test suites (+2 more)

### Community 212 - "ofstream"
Cohesion: 0.13
Nodes (19): ofstream, T, time_point, Entry, object, os, path, v32zt_close_allobj_writers() (+11 more)

### Community 213 - "mutex"
Cohesion: 0.11
Nodes (23): main(), build_orbit_V(), clear_fast_caches(), Entry, tgt, w, fast_projector_cached(), FastProjector (+15 more)

### Community 214 - "fit.py"
Cohesion: 0.06
Nodes (48): build_parser(), code_equivalence(), cutoff_scan(), _document(), error_gap_flags(), _fd_hessian_errors(), fit_cutoff(), FitResult (+40 more)

### Community 215 - "CLAUDE.md — project rules for 3body-qc-v4"
Cohesion: 0.22
Nodes (8): CLAUDE.md — project rules for 3body-qc-v4, Decisions already made by Digonto (2026-10-21/22; details in `notes/STATUS.md`), Orientation, Physics invariants (never change silently), Pushing to GitHub (pre-authorised by Digonto, 2026-10-21), Setting up a fresh (cloud) checkout, Stop and ask Digonto, Things that are known to be wrong or unverified (see `docs/OPEN_ISSUES.md`)

### Community 216 - "DeviceConfig"
Cohesion: 0.08
Nodes (21): DeviceConfig, ell, m, n, nx, ny, nz, DeviceConfigView (+13 more)

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

### Community 226 - "fast_fclass.hpp"
Cohesion: 0.06
Nodes (38): main(), build_ftable(), class_moments(), FClass, cnt, ep, nb, perp2 (+30 more)

### Community 227 - "read_gpu_cache_complex"
Cohesion: 0.67
Nodes (3): comp, istream, read_gpu_cache_complex()

### Community 228 - "EnergyShapeKey"
Cohesion: 0.33
Nodes (5): EnergyShapeKey, dim1, dim2, total_dim, vdim

### Community 229 - "V32zuF2K2ValidationStats"
Cohesion: 0.33
Nodes (6): V32zuF2K2ValidationStats, chunks, max_abs_diff_F2, max_abs_diff_K2inv, max_rel_diff_F2, max_rel_diff_K2inv

### Community 230 - "000 — Consolidation and roadmap (starter plan from chat, 2026-10-21)"
Cohesion: 0.22
Nodes (9): 000 — Consolidation and roadmap (starter plan from chat, 2026-10-21), Goal, Milestone 0 — consolidation (DONE, 2026-10-21), Milestone 1 — trustworthy baseline (decisions answered 2026-10-21; step 6 replaced, see the update above), Milestone 2 — root finding without human labels ("gradient change"), Milestone 3 — L24 and combined fit, Milestone 4 — speed and robustness, Milestone 5 — physics output (+1 more)

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

### Community 239 - "014 M2 step 1: interface stubs + `qc3.free_levels` (tester, mode step)"
Cohesion: 0.40
Nodes (4): 014 M2 step 1: interface stubs + `qc3.free_levels` (tester, mode step), Commands, Problems, Results

### Community 240 - "module.cpp"
Cohesion: 0.20
Nodes (17): array_t, _core, dict, list, compare_kinematics(), complex, m, MatrixXcd (+9 more)

### Community 241 - "Options"
Cohesion: 0.40
Nodes (5): Options, chop_tol, debug, eig_tol, threads_per_block

### Community 242 - "HISTORY.md"
Cohesion: 0.25
Nodes (4): HISTORY — where this code came from, Lessons that must not be relearned, Phases (condensed), Source packages

### Community 243 - "v33p_gpu_det_lu.cu"
Cohesion: 0.25
Nodes (14): cublasStatus_t, cudaError_t, size_t, string, uint64_t, cublas_error(), cublas_ok(), cuda_error() (+6 more)

### Community 244 - "fgk.py"
Cohesion: 0.08
Nodes (46): expected_basis(), kept_momenta(), Basis rule (docs/METHODS.md section 2; prompts/001 section 2.1): spectator…, Integer vectors kept for the given channel, in loop order., Structured array (BASIS_DTYPE) in the production order., cube_points(), F_block(), F_matrix() (+38 more)

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
Cohesion: 0.19
Nodes (18): main(), map, evaluate_full_grid_parallel(), get_double(), get_int(), get_int_list(), get_string(), main() (+10 more)

### Community 251 - "Mom"
Cohesion: 0.21
Nodes (12): LD, brute(), direct_double(), main(), Mom, A0, A1, A2 (+4 more)

### Community 253 - "3. Results"
Cohesion: 0.17
Nodes (12): 1. What was asked, 2. Actions taken, 3.1 Formalism and references (as stated in the code), 3.2 Question 1 — the canonical, newest physics core, 3.4 Question 3 — variants and dead ends, 3.5 Question 4a — fixes visible by diffing old (V3) vs newest, 3.6 Question 4b — problems I found by reading (not noted in the code), 3.7 Question 5 — recommended minimal physics-core file set (+4 more)

### Community 254 - "QCSearchTiming"
Cohesion: 0.50
Nodes (4): QCSearchTiming, candidate_generation_sec, classifier_sec, determinant_scan_sec

### Community 258 - "Test report 009: M1 step 7.2 (orbit/character projector) and 7.3 (threshold basis builder, `dimension_jumps`)"
Cohesion: 0.15
Nodes (12): 7.2 physics review, 7.3 review, Classification of the three tests the caller flagged, Commands run (repo root; scripts are copied to `results/tester_009_scripts/`), Failures (both pending on later sub-steps), Item 6: timings (taskset -c 0, OMP=OPENBLAS=1, $\vec c_{\rm base}$, `want=("Q",)`), $\ln\lvert\det Q\rvert$, $\vec c=\vec c_{\rm base}$, `coded_v2`, 2 energies per irrep, Notes (non-blocking) (+4 more)

### Community 260 - "test_qc_lattice_roots.py"
Cohesion: 0.05
Nodes (58): chi2(), cov_eig_pinv(), elab_to_ecm(), jackknife_mean_cov(), jackknife_resample(), _lattice_dir(), load_levels(), match_levels() (+50 more)

### Community 307 - "fast_projector.hpp"
Cohesion: 0.11
Nodes (22): s, build_orbit_V(), Entry, tgt, w, fast_projector_cached(), FastProjector, method (+14 more)

### Community 308 - "test_m2_step5_chi2.py"
Cohesion: 0.15
Nodes (16): FakeTracker, make_levels(), Unit tests for qc3.fit.QCChi2 / make_trackers (M2 step 5, prompts/003 section…, Affine roots E_k(c) = a_k + J_k . (c - C0); `slopes` < 0 marks wrong-sign…, 110_A2, coded_v2, baseline, window [0.3190, 0.3200]: one genuine zero…, setup(), test_chi2_equals_lattice_chi2_of_the_matched_roots_and_d3_uses_lowest(), test_constructor_validation() (+8 more)

### Community 310 - "real_wigner_d.hpp"
Cohesion: 0.20
Nodes (28): wigner_d_tests(), approx_equal(), approx_equal_matrix(), check_m_range(), D_real_element(), D_real_matrix(), D_real_matrix_proper(), factorial_int() (+20 more)

### Community 311 - "DeviceResidentResult"
Cohesion: 0.18
Nodes (11): DeviceResidentResult, build_sec, det, dim1, dim2, Ecm, En, i (+3 more)

### Community 312 - "WindowGridRow"
Cohesion: 0.29
Nodes (7): WindowGridRow, det_imag, det_real, Ecm, Nfull, Nproj, row

### Community 316 - "F — Development history v31 → v32y: lessons, bugs fixed, dead ends"
Cohesion: 0.18
Nodes (11): 0. Glossary (read first), 1. What was asked, 2. What I did, 3. Pre-v31 background (from the v31u context-transfer package), 6. Gradient / slope / derivative use (Question 3), 7.1 Recorded in the archive at the time, 7.2 Noticed during this review (**[my observations]**, need confirmation by running), 7. Open questions (Question 4) (+3 more)

### Community 317 - "anchor_points"
Cohesion: 0.11
Nodes (34): slow, anchor_points(), as_c(), basis_list(), E_lab_of(), load_anchors(), oracle_point(), ||a-b||_F / ||b||_F (b is the reference). (+26 more)

### Community 318 - "dig_tools.hpp"
Cohesion: 0.14
Nodes (11): Func, string, T, time_point, ManualTimer, t0, print_function_output_to_file(), printer() (+3 more)

### Community 319 - "test_rootscan.py"
Cohesion: 0.12
Nodes (30): _cand(), Fake, parametrize, Unit tests for qc3.rootscan (M2 step 2, prompts/003 section 2.2): synthetic det…, det Q(E) = s0 prod (E-z)^m / prod (E-p) * exp(bg(E)) with Q = diag(det, 1) (so…, scan(), test_dimension_jump_is_a_segment_boundary_not_a_root(), test_even_order_zero_without_sign_change_is_found_by_the_minimum_rule() (+22 more)

### Community 320 - "RepDiagnostics"
Cohesion: 0.20
Nodes (10): RepDiagnostics, best_convention_index, closure_AB_max, closure_BA_max, closure_best_max, group_size, projector_diff_old_rel, projector_herm_res (+2 more)

### Community 321 - "MatrixCompareStats"
Cohesion: 0.15
Nodes (14): compare_matrices(), MatrixCompareStats, frob_abs, frob_rel, max_abs, max_i, max_j, max_rel (+6 more)

### Community 322 - "Sector"
Cohesion: 0.07
Nodes (56): Scan and score the labelled sectors; returns {"L20/000_A1m": record, ...}.…, run_gate(), ndarray, Production sector (CLAUDE.md invariant 3); `overrides` replace individual…, Sector, test_production_default_k2_pwave_form_is_q0_consistent(), _proj(), parametrize (+48 more)

### Community 323 - "/kickoff — start the physics loop here"
Cohesion: 0.22
Nodes (8): 1. Scaffold (idempotent — never overwrite existing files), 2. Starter plan → prompts/000, 3. Code graph, 4. First commit, 5. Planner (kickoff mode), 6. Stubs and lock, 7. STATUS and go, /kickoff — start the physics loop here

### Community 324 - "AssignmentRow"
Cohesion: 0.18
Nodes (11): AssignmentRow, assignment_status, classifier_candidate_id, irrep, Lbyas, mode, model_root_Ecm, residual (+3 more)

### Community 325 - "QCRefinedFCN"
Cohesion: 0.11
Nodes (17): atomic, Minuit2::FCNBase, shared_ptr, TargetLevel, QCRefinedFCN, cn, coarse_cache, corr (+9 more)

### Community 326 - "Test report 010: M1 step 7.4 ($\mathcal K_{\rm df,3}$ pair kernel), 7.5 ($G$ pair kernel), 7.6 (C++ tail)"
Cohesion: 0.15
Nodes (12): Commands run (repo root; own scripts are in `results/tester_010_scripts/`, logs in `results/tester_010_logs/`), Failing output, Items not done / limits, N1. `test_fast_single_thread_cpu_time` alone: cause is OpenBLAS start-up, not the fast path (verified), N2. At $\kappa(Q)\gtrsim10^{13}$ the oracle is less accurate than the fast path (verified with 60 digits), N3. $G$ fill mode: recommendation to keep the default "every ordered pair", N4. Minor, Notes (+4 more)

### Community 327 - "Tester report 012: M1 sub-step 7.7, class moment tables + shared (canonical $\mathbf n_P$, $|\mathbf p|^2$) sums for the $F$ lattice sum"
Cohesion: 0.20
Nodes (9): 1. Gate checks, 2. Math review, 3. Independent numerics (coded_v2, $\vec c_{\rm base}$, seed 20261022), 4. Memory and file I/O, 5. Timings (`taskset -c 0`, OMP and OPENBLAS = 1, median of 21 warm / 9 cold, `want=("Q",)`, coded_v2), Commands run, Notes (none blocks), Return summary (+1 more)

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

### Community 332 - "3body-qc-v4 — fitting the three-particle K matrix to lattice spectra"
Cohesion: 0.29
Nodes (7): 3body-qc-v4 — fitting the three-particle K matrix to lattice spectra, Build, Layout, Running a fit needs external inputs, Status (2026-10-21), The hard part: finding the right roots, The idea that makes fits fast

### Community 333 - "F3iProjectedCache"
Cohesion: 0.29
Nodes (6): F3iProjectedCache, hits, misses, table, map, size_t

### Community 334 - "Research 008 — M1 step 7 re-plan: cold projector, basis builder, $\mathcal K_3$/$G$ kernels, tail, and pole/zero discrimination"
Cohesion: 0.20
Nodes (9): 1. Where the time goes now (report 006, 1 core), 2. Projector from orbits (candidate (a)) — method and evidence, 3. Basis builder (candidate (b)), 4. $\mathcal K_{\rm df,3}$ and $G$ kernels (candidates (c), (d)), 5. Tail (candidate (e)), 6. $F$ (optional), 7. Pole/zero discrimination (coordinator addendum, step 9 findings), 8. Ranking (gain on realistic use = a scan or a fit evaluation over a window) (+1 more)

### Community 335 - "F3_device_resident_true_batch_v21_side_abs_trend_refine.cu"
Cohesion: 0.22
Nodes (18): build_F3_from_F2_and_F2X1_kernel(), c_make(), combine_H_colmajor_kernel(), ConfigEntry, cuDoubleComplex, __global__, gather_real_evecs_to_complex_Vsel_kernel(), v14_build_F2_batch_same_shape_kernel() (+10 more)

### Community 336 - "auditor.md"
Cohesion: 0.40
Nodes (4): Checklist (state each item: OK / ISSUE / N/A, with file:line), Report (always), Return to the caller (short), Start

### Community 337 - "builder.md"
Cohesion: 0.40
Nodes (4): Build, Not allowed (a guard hook enforces most of this), Return to the coordinator (short), Start

### Community 338 - "generate_K3iso.cpp"
Cohesion: 0.21
Nodes (12): comp, string, generate_K3iso_L20_from_F3(), generate_K3iso_L20_from_F3inv(), generate_K3iso_using_F3_vs_En_KKpi_omp_single_irrep_centralvalue(), generate_K3iso_using_F3_vs_En_KKpi_omp_single_irrep_with_bounds(), main(), det_F3mat_poles_secant_method() (+4 more)

### Community 340 - "session_start.py"
Cohesion: 0.83
Nodes (3): latest(), main(), sh()

### Community 341 - "LocalScopedTimer"
Cohesion: 0.33
Nodes (5): time_point, LocalScopedTimer, debug, name, t0

### Community 343 - "point.py"
Cohesion: 0.13
Nodes (16): _ev(), find_root(), main(), M1 step 9: p-wave K2^{-1} form quantification for decision D7 (prompts/001 sec…, Plan rule: refine_root_info (default kind="zero": poles raise PoleError) in…, Supplementary (not the plan rule): qc3.roots.find_sign_changes on npts points;…, read_windows(), scan_roots() (+8 more)

### Community 344 - "test_K3mat_2plus1_cpp_dump.cpp"
Cohesion: 0.39
Nodes (7): comp, MatrixXcd, string, vector, main(), write_config_file(), write_matrix_file()

### Community 345 - "Report 007: M1 step 9, p-wave $\mathcal K_2^{-1}$ level-shift table (D7)"
Cohesion: 0.29
Nodes (6): Commands, Findings (not step-9 failures), Report 007: M1 step 9, p-wave $\mathcal K_2^{-1}$ level-shift table (D7), Results, Sign argument (check 5), Verdict

### Community 346 - "3. Q1 - Which fitter source is newest, and what changed"
Cohesion: 0.33
Nodes (6): 3.1 Hash/size table of the multi-L fitter `qc_fitter_norm_refine_v2_multiL.cpp`, 3.2 What changed V33P (1373) $\to$ MF (2021): the hot-window FCN (+648 / -~60 lines), 3.3 What changed MF (2021) $\to$ L20 (2107): chi-square audit trace only (88 diff lines), 3.4 Other L20-only refinements outside the fitter, 3.5 Earlier evolution (from docs), 3. Q1 - Which fitter source is newest, and what changed

### Community 347 - "Tester report 006: M1 step 6, fast path v1 (rungs 1-4)"
Cohesion: 0.20
Nodes (9): Check 3: step-6 test ids (pinned to one core, `OMP_NUM_THREADS=1`), Check 6: 12 anchors, $\vec c$=baseline, 1 thread, pinned (`results/m1_step6_timings_1thread.json`), Commands run, Failures, Notes (non-blocking; none changes a physics number), Results, Tester report 006: M1 step 6, fast path v1 (rungs 1-4), Verdict line (+1 more)

### Community 348 - "V3Shoulder"
Cohesion: 0.40
Nodes (5): V3Shoulder, kind, n, reason, score

### Community 349 - "7. Q4 - Fit results achieved"
Cohesion: 0.33
Nodes (6): 7.1 v32w (trusted single-L baseline, 15 Jun), 7.2 v33n cutoff scan (L20+L24, 9 Jul; computed with `raw_sign_only`, v33g runtime cache, full scan), 7.3 L20 pass-only projected-basis hot-window fits (15 Jul), 7.4 K3B/K3E zero tests (`L20/output/l20_k3b_k3e_zero_tests/L20_K3B_K3E_ZERO_TEST_FINAL_REPORT.md:7-61`), 7.5 Blockers, 7. Q4 - Fit results achieved

### Community 350 - "6. Problems found during this review (read before trusting any L20/v33n chi-square)"
Cohesion: 0.40
Nodes (5): 6.1 The multi-L covariance is probably inflated by a factor $(N-1)^2 \approx 3.6\times10^{5}$, 6.2 Fit-quality caveats independent of 6.1, 6.3 Complex read convention: two swaps that cancel, 6.4 Portability of the release tree, 6. Problems found during this review (read before trusting any L20/v33n chi-square)

### Community 351 - "Audit 002: M2 milestone (label-free roots, tracking, L20 fit)"
Cohesion: 0.33
Nodes (5): Answers to the focus questions, Audit 002: M2 milestone (label-free roots, tracking, L20 fit), Checklist, Notes (all NON-BLOCKING), VERDICT: PASS WITH NOTES (no blocking issue)

### Community 353 - "SignFlip"
Cohesion: 0.50
Nodes (4): SignFlip, E_linear, left, right

### Community 354 - "EnergyPrepassInfo"
Cohesion: 0.25
Nodes (8): EnergyPrepassInfo, dim1, dim2, Ecm, En, i, total_dim, vdim

### Community 355 - "Tester report 016: M2 step 4 (qc3.tracking)"
Cohesion: 0.40
Nodes (4): Code reading (robustness and loudness), Results, Tester report 016: M2 step 4 (qc3.tracking), Verdict

### Community 356 - "Benchmark report 001: M1 single-energy QC, fast path vs oracle (milestone gate)"
Cohesion: 0.20
Nodes (9): Benchmark report 001: M1 single-energy QC, fast path vs oracle (milestone gate), Caveats, $\chi^2$ proxy (12 labelled L20 roots, cutoff 0.335, `coded_v2`, fast backend, `refine_root_info(kind="zero")`, 1 core, no caches on disk), Commands, Comparison with the history (minima where noise matters), Machine and build, Memory, Per-point table (official run; core 0, 1 thread) (+1 more)

### Community 357 - "CacheKey"
Cohesion: 0.50
Nodes (3): CacheKey, L, label

### Community 358 - "test_m1_fast_vs_oracle.py"
Cohesion: 0.24
Nodes (14): fast_point(), grid_energy(), _key(), _compare(), slow, M1 steps 6-7 — the cache-free fast path reproduces the oracle at every energy…, Generalisation check: seeded random energies over the whole production range,…, test_fast_basis_identical_to_oracle() (+6 more)

### Community 360 - "5. The M2 method (plan `prompts/003`; research `reports/research/009`)"
Cohesion: 0.13
Nodes (15): 1. Why it is hard (measured), 2. What was tried (oldest → newest), 3.1 How to compute it cheaply with the existing caches, 3.2 Smooth root position for Minuit (proposal), 3.3 Acceptance gate for any new method, 3. The "gradient change" idea, stated precisely, 4. Rules that apply to every method, 5.1 Label-free scan — `qc3.rootscan.scan_sector` (+7 more)

### Community 361 - "EnergyMatrixPack"
Cohesion: 0.14
Nodes (15): cuDoubleComplex, MatrixXcd, EnergyMatrixPack, Ecm, En, F2, G, H (+7 more)

### Community 362 - "v32w_single_L_fitter/qc_fitter_norm_refine_v2.cpp"
Cohesion: 0.21
Nodes (14): classify(), cnorm(), Eval, det, E, proj_dim, success, y (+6 more)

### Community 363 - "F3_cpu_openmp_v24_K3_QC_multi_signflip_side_abs_refinement.cpp"
Cohesion: 0.38
Nodes (13): add_unique_zero(), classify_one_subflip_side_abs(), pair, string, vector, evaluate_one_energy_cpu_eigen_v24_QC(), find_sign_flips(), L_tag() (+5 more)

### Community 364 - "scanN.py"
Cohesion: 0.21
Nodes (8): ev(), ev(), grid(), N(), Ng(), Independent scan of the basis dimension N(E) with the ORACLE's config_maker…, rec(), ev()

### Community 365 - "EnergyPrepassInfo"
Cohesion: 0.25
Nodes (8): EnergyPrepassInfo, dim1, dim2, Ecm, En, i, total_dim, vdim

### Community 366 - "string"
Cohesion: 0.36
Nodes (13): cov_to_corr(), map, MatrixD, string, gd(), gi(), gs(), main() (+5 more)

### Community 367 - "Tester report 011: M1 step 8, benchmark harness `qc3.bench` and the $\chi^2$ proxy"
Cohesion: 0.25
Nodes (7): Commands run, Counts, How the independent $\chi^2$ was computed, Notes (none blocking), Results, Tester report 011: M1 step 8, benchmark harness `qc3.bench` and the $\chi^2$ proxy, Verdict: PASS (with minor notes)

### Community 368 - "acc.py"
Cohesion: 0.50
Nodes (3): maxblock(), Tester 012: F (classes) vs oracle vs use_classes=False; per spectator block;…, rel()

### Community 369 - "Tester report 013: M1 milestone gate, mode `full`"
Cohesion: 0.40
Nodes (4): Commands run, Notes, Results, Tester report 013: M1 milestone gate, mode `full`

### Community 370 - "score_labels.py"
Cohesion: 0.24
Nodes (12): candidate_record(), exit_code(), gate_failed(), _git(), _jsonable(), main(), metadata(), M2 step 3 (prompts/003 section 3): the label gate. Runs the label-free scan… (+4 more)

### Community 371 - ".upload"
Cohesion: 0.30
Nodes (8): host_comp, MatrixXcd, vector, DeviceScatterParams, d_im, d_re, K2inv_EREord2_2plus1_many_mat_gpu_safe(), K2inv_EREord2_2plus1_mat_gpu_safe()

### Community 373 - "InwardShapeResult"
Cohesion: 0.25
Nodes (8): InwardShapeResult, classification, decreasing_fraction, ezero, flip_left_index, flip_right_index, increasing_fraction, trend_score

### Community 375 - "find_QC_zeros_refined"
Cohesion: 0.56
Nodes (9): build_entry(), FitSettings, K3dfParameters, PhysicsParams, ProjectedQCCacheEntry, eval_entry_QC(), find_QC_zeros_refined(), mesh() (+1 more)

### Community 376 - "rootscan.py"
Cohesion: 0.08
Nodes (32): free_level_states(), free_levels(), FreeLevel, _kaon_momenta(), ndarray, Non-interacting K K pi levels per sector (prompts/003 section 2.1; D8 input;…, All integer vectors n (lexicographic order) that a kaon can carry below…, [(E*, (n1, n2, n3))] for every state with E_min - margin <= E* <= E_max +… (+24 more)

### Community 377 - "Tester report 019: M2 gate, mode full"
Cohesion: 0.29
Nodes (6): Commands and results, Hesse reproducibility (Tables A and B), Return notes, Tester report 019: M2 gate, mode full, The one failure, VERDICT: FAIL (78 of 79 acceptance tests pass; 1 fails, an environment-dependent anchor tolerance, not a code change)

### Community 378 - "vector"
Cohesion: 0.36
Nodes (6): IrrepCache, pair, vector, print_cache_summary(), print_split_cache_summary(), split_fixed_coarse_and_refined()

### Community 379 - "DetInfo"
Cohesion: 0.25
Nodes (8): MatrixXcd, det_info(), DetInfo, det_re, logabs, ok, sign, slog

### Community 380 - "ScanPoint"
Cohesion: 0.18
Nodes (11): ScanPoint, abs_norm, det, Ecm, En, i, success, total_dim (+3 more)

### Community 381 - "MatrixCompareStats"
Cohesion: 0.18
Nodes (11): MatrixCompareStats, frob_abs, frob_ref, frob_rel, max_abs, max_i, max_j, max_rel (+3 more)

### Community 382 - "default_data_dir"
Cohesion: 0.22
Nodes (10): default_data_dir(), Data root: $QC_DATA_ROOT if set (not checked for existence), else <repo>/data.…, Rows {"bracket_id": int, "E_left": float, "E_right": float, "label":…, read_labels(), _files(), M1 step 5 — no machine-specific absolute paths on the M1 production path;…, test_data_root_env_var(), test_no_absolute_home_paths_in_m1_code() (+2 more)

### Community 383 - "4. Q2 - Exactly how the fast FCN works (L20 code)"
Cohesion: 0.20
Nodes (10): 4.1 The quantization condition as implemented, 4.2 What is precomputed and cached (three cache levels), 4.3 Startup (outside the FCN), 4.4 Accepted windows (how each lattice level gets its search window), 4.5 Inside one FCN call, 4.6 Chi-square and covariance, 4.7 Minuit2 setup and parameter mask, 4.8 Run modes of the binary `bin/v33f_k3df_fitter_multiL_v33e <config> <mode>` (+2 more)

### Community 384 - "DetRow"
Cohesion: 0.25
Nodes (9): classify_inward_shape(), classify_inward_shape_result(), DetRow, det_norm, det_raw, Ecm, i, linear_zero() (+1 more)

### Community 385 - "string"
Cohesion: 0.23
Nodes (13): adaptive_final_zeros_device_resident(), choose_projection_parity_v20(), array, string, time_point, get_energy_shape_prepass(), main(), process_same_shape_group_stream_parallel() (+5 more)

### Community 386 - "ScopedTimer"
Cohesion: 0.33
Nodes (5): time_point, ScopedTimer, debug, name, t0

### Community 387 - "4. Chronological table v31v → v32y (Question 1)"
Cohesion: 0.22
Nodes (9): 4. Chronological table v31v → v32y (Question 1), Phase A — $F_2$/$G$/$FG$ diagnostics, projection audit, blind pole classifier (default `111_A2` → `110_A2`, $L=20$), Phase B — projected $F_3^{-1}$ zero finder (`110_A2`, $L=24$, waves $\{0,1\}$), Phase C — toy-mass sanity tests and scattering-parameter mapping, Phase D — classifier v1, s vs s+p waves, p-wave strength, $v_{\rm iso}$, Phase E — iterative refinement and lattice masses, Phase F — K3df fitter, caching, diagnostics, Phase G — classifier v2, normalization, true refinement, fitter baseline (+1 more)

### Community 389 - "3. Cache formats and stored objects (Question 1)"
Cohesion: 0.29
Nodes (7): 3.1 Physics objects stored per energy row, 3.2 Energy grid, 3.3 Combined raw cache `v32zu_Lbyas<L>_xi<xi>_irrep<label>_coarse<N>_F3inv_Vsel_gpu.bin`, 3.4 Combined-cache sidecar `<bin>.meta.json`, 3.5 Component caches (format `KKPI_COMPONENT_CACHE_V001`), CPU and GPU, 3.6 The `variant_04_real_imag_swapped` convention: where it comes from, 3. Cache formats and stored objects (Question 1)

### Community 391 - "LocalScopedTimer"
Cohesion: 0.33
Nodes (5): time_point, LocalScopedTimer, debug, name, t0

### Community 392 - "9. Q6 - Recommended minimal fitter file set for the consolidated repo"
Cohesion: 0.33
Nodes (6): 9.1 C++ sources (fitter executable `v33f_k3df_fitter_multiL_v33e`), 9.2 Python scripts (under `L20/scripts/`), 9.3 Configs and data, 9.4 Not needed in the consolidated fitter, 9.5 Suggested consolidation design for a fast FCN (based on what works), 9. Q6 - Recommended minimal fitter file set for the consolidated repo

### Community 393 - "F2GpuOptions"
Cohesion: 0.40
Nodes (5): F2GpuOptions, debug, erfi_max_terms, erfi_tol, threads_per_block

### Community 394 - "test_bench.py"
Cohesion: 0.25
Nodes (3): fixture, Unit tests for qc3.bench (prompts/001 step 8): CLI parsing, JSON schema, chi^2…, restore_env()

### Community 395 - "comp"
Cohesion: 0.50
Nodes (3): comp, MatrixXcd, determinant_via_partial_piv_lu()

### Community 396 - "SignFlip"
Cohesion: 0.50
Nodes (4): SignFlip, E_linear, left, right

### Community 398 - "Options"
Cohesion: 0.67
Nodes (3): Options, debug, threads_per_block

## Knowledge Gaps
- **3594 isolated node(s):** `0. Diagnosis (why the old pipeline is slow and what is wrong)`, `1. Goal of M1`, `2.1 Equations and conventions (unchanged physics; see `docs/METHODS.md` §§1–8 for the full set)`, `2.2 p-wave $\mathcal K_2^{-1}$: two named forms, production default unchanged until D6`, `2.3 Speed-up candidates, ranked (expected gain, how verified)` (+3589 more)
  These have ≤1 connection - possible missing edges or undocumented components. (Counts symbols only; 4204 node(s) total have ≤1 connection when file, concept and rationale nodes are included.)
- **64 thin communities (<3 nodes) omitted from report** — run `graphify query` to explore isolated nodes.

## Suggested Questions
_Questions this graph is uniquely positioned to answer:_

- **Why does `FitSettings` connect `FitSettings` to `vector`, `fitter/K3df_minuit_fit_v32f_fullF3inv_QCfull_cached_classifier.hpp`, `TargetLevel`, `FitResult`?**
  _High betweenness centrality (0.016) - this node is a cross-community bridge._
- **Why does `FitSettings` connect `FitSettings` to `QCPointValue`, `v32w_single_L_fitter/K3df_minuit_fit_v32f_fullF3inv_QCfull_cached_classifier.hpp`, `FitResult`?**
  _High betweenness centrality (0.016) - this node is a cross-community bridge._
- **Why does `EnergyMatrixPack` connect `EnergyMatrixPack` to `ThreadSafeQueue`, `CpuGpuChunk`, `F3_gpu_omp_cublas_pipeline_v9_inward_shape_refinement.cu`?**
  _High betweenness centrality (0.012) - this node is a cross-community bridge._
- **What connects `0. Diagnosis (why the old pipeline is slow and what is wrong)`, `1. Goal of M1`, `2.1 Equations and conventions (unchanged physics; see `docs/METHODS.md` §§1–8 for the full set)` to the rest of the system?**
  _3594 weakly-connected nodes found - possible documentation gaps or missing edges._
- **Should `Code and function index` be split into smaller, more focused modules?**
  _Cohesion score 0.008849557522123894 - nodes in this community are weakly interconnected._
- **Should `Package `corrected_reader`` be split into smaller, more focused modules?**
  _Cohesion score 0.011764705882352941 - nodes in this community are weakly interconnected._
- **Should `Package `production_fitter`` be split into smaller, more focused modules?**
  _Cohesion score 0.011834319526627219 - nodes in this community are weakly interconnected._
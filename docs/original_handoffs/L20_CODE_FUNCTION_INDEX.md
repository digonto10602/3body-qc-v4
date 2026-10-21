# Code and function index

Generated from the packaged C++/header/Python sources.

## `source/F2_functions.h`

Purpose: C++/header implementation.

Important functions/classes:
- `F2_i1`
- `F2_i1_1`
- `F2_i_mat`
- `F2_i_mat_1`
- `I00_sum_F`
- `I00_sum_F_1`
- `I0F`
- `if`

Inputs/outputs: inspect the file CLI and adjacent reports; this index is intentionally mechanical.

## `source/F2_functions_v2.h`

Purpose: C++/header implementation.

Important functions/classes:
- `F2_2plus1_mat`
- `F2_ang_mom`
- `F2_i_ang_mom_mat`
- `I0F`
- `I1F`
- `I2F`
- `I_int_ang_mom`
- `I_sum_ang_mom`
- `add`
- `comp`
- `debug_F2_2plus1_component_for_lm_steps`
- `debug_F2_ang_mom_steps`
- `debug_I_int_ang_mom_steps`
- `debug_I_sum_ang_mom_steps`
- `for`
- `if`
- `precise_vector_sum`
- `precise_vector_sum_debug`
- `print_Ylm1_Ylm2_for_ijk`
- `print_and_test_vector_sum`
- `print_comp_vec3_debug`
- `result`
- `smallest_eigenvalue`
- `sorted_neumaier_sum`
- `value`

Inputs/outputs: inspect the file CLI and adjacent reports; this index is intentionally mechanical.

## `source/F2_gpu_safe_builder.cuh`

Purpose: C++/header implementation.

Important functions/classes:
- `DeviceBuffer`
- `DeviceFlatConfig`
- `F2_2plus1_many_mat_gpu_safe`
- `F2_2plus1_mat_gpu_safe`
- `F2_ang_mom_d`
- `I0F_d`
- `I1F_d`
- `I2F_d`
- `I_int_ang_mom_d`
- `I_sum_ang_mom_d`
- `Jfunc_d`
- `allocate`
- `build_F2_2plus1_kernel`
- `check_cuda`
- `compare_F2_cpu_gpu`
- `copy_from_host`
- `copy_to_host`
- `cutoff_function_1_d`
- `erfi_series_d`
- `flatten_config`
- `if`
- `kallentriangle_d`
- `omega_func_d`
- `q2psq_star_d`
- `release`
- `sigma_pvec_based_d`
- `spherical_harmonics_d`
- `zabs`
- `zadd`
- `zdiv`
- `zdiv_real`
- `zexp`
- `zim`
- `zmake`
- `zmul`
- `zmul_real`
- `zneg`
- `zpow_int`
- `zpow_real_positive_base`
- `zre`
- `zsqrt`
- `zsub`

Inputs/outputs: inspect the file CLI and adjacent reports; this index is intentionally mechanical.

## `source/F3_cpu_openmp_v21_side_abs_refinement.cu`

Purpose: C++/header implementation.

Important functions/classes:
- `F3ScopedTimer`
- `add_unique_zero`
- `add_unique_zero_cpu_v21`
- `append_failed_results_for_chunk`
- `build_one_energy_pack`
- `bytes_complex_matrix`
- `bytes_to_mib`
- `catch`
- `classify_refined_file_by_side_abs_trend`
- `classify_sign_flips_zero_vs_pole`
- `close`
- `ensure`
- `estimate_gpu_bytes_for_items`
- `estimate_gpu_bytes_one_energy`
- `f3_cudaSyncCheck`
- `f3_exception_looks_like_cuda_oom`
- `f3_has_sign_flip`
- `f3_make_auto_hdf5_matrix_filename`
- `f3_packEigenToHostCublas`
- `f3_should_save_matrix_point`
- `f3_sign_double`
- `f3_throwOnCublas`
- `f3_throwOnCuda`
- `f3_unpackHostCublasToEigen`
- `find_all_sign_flips_in_rows`
- `find_sign_flip_nearest_middle`
- `flush_cpu_buffer_to_gpu_queue`
- `for`
- `gpu_consumer_thread_func`
- `gpu_process_chunk_oom_safe_recursive`
- `gpu_process_same_dim_chunk_cublas`
- `gpu_process_same_dim_same_vdim_group_cublas`
- `if`
- `linear_zero_estimate`
- `main`
- `make_scan_tag`
- `make_subchunk_copy`
- `pop`
- `push`
- `read_normalized_det_file_for_sign_flips`
- `release`
- `run_F3_energy_scan_gpu_omp_cublas_pipeline_v4`
- `test_F3_with_pwave_all_energy_cpu_openmp_v21_side_abs_refinement`
- `test_F3_with_pwave_all_energy_gpu_omp_cublas_pipeline_v3`
- `test_F3_with_pwave_all_energy_gpu_omp_cublas_pipeline_v8_save_matrices`
- `while`

Inputs/outputs: inspect the file CLI and adjacent reports; this index is intentionally mechanical.

## `source/F3_cpu_openmp_v22_multi_signflip_side_abs_refinement.cu`

Purpose: C++/header implementation.

Important functions/classes:
- `F3ScopedTimer`
- `add_unique_zero`
- `add_unique_zero_cpu_v21`
- `append_failed_results_for_chunk`
- `build_one_energy_pack`
- `bytes_complex_matrix`
- `bytes_to_mib`
- `catch`
- `classify_one_signflip_side_abs_trend_v22`
- `classify_refined_file_by_multi_side_abs_trend_v22`
- `classify_sign_flips_zero_vs_pole`
- `close`
- `ensure`
- `estimate_gpu_bytes_for_items`
- `estimate_gpu_bytes_one_energy`
- `f3_cudaSyncCheck`
- `f3_exception_looks_like_cuda_oom`
- `f3_has_sign_flip`
- `f3_make_auto_hdf5_matrix_filename`
- `f3_packEigenToHostCublas`
- `f3_should_save_matrix_point`
- `f3_sign_double`
- `f3_throwOnCublas`
- `f3_throwOnCuda`
- `f3_unpackHostCublasToEigen`
- `find_all_sign_flips_in_rows`
- `find_sign_flip_nearest_middle`
- `flush_cpu_buffer_to_gpu_queue`
- `for`
- `format_Lbyas_for_filename`
- `gpu_consumer_thread_func`
- `gpu_process_chunk_oom_safe_recursive`
- `gpu_process_same_dim_chunk_cublas`
- `gpu_process_same_dim_same_vdim_group_cublas`
- `if`
- `linear_zero_estimate`
- `main`
- `make_scan_tag`
- `make_subchunk_copy`
- `pop`
- `push`
- `read_normalized_det_file_for_sign_flips`
- `release`
- `run_F3_energy_scan_gpu_omp_cublas_pipeline_v4`
- `test_F3_with_pwave_all_energy_cpu_openmp_v22_multi_signflip_side_abs_refinement`
- `test_F3_with_pwave_all_energy_gpu_omp_cublas_pipeline_v3`
- `test_F3_with_pwave_all_energy_gpu_omp_cublas_pipeline_v8_save_matrices`
- `while`

Inputs/outputs: inspect the file CLI and adjacent reports; this index is intentionally mechanical.

## `source/F3_cpu_openmp_v22_multi_signflip_side_abs_refinement_Lbyas_arg.cu`

Purpose: C++/header implementation.

Important functions/classes:
- `F3ScopedTimer`
- `add_unique_zero`
- `add_unique_zero_cpu_v21`
- `append_failed_results_for_chunk`
- `build_one_energy_pack`
- `bytes_complex_matrix`
- `bytes_to_mib`
- `catch`
- `classify_one_signflip_side_abs_trend_v22`
- `classify_refined_file_by_multi_side_abs_trend_v22`
- `classify_sign_flips_zero_vs_pole`
- `close`
- `ensure`
- `estimate_gpu_bytes_for_items`
- `estimate_gpu_bytes_one_energy`
- `f3_cudaSyncCheck`
- `f3_exception_looks_like_cuda_oom`
- `f3_has_sign_flip`
- `f3_make_auto_hdf5_matrix_filename`
- `f3_packEigenToHostCublas`
- `f3_should_save_matrix_point`
- `f3_sign_double`
- `f3_throwOnCublas`
- `f3_throwOnCuda`
- `f3_unpackHostCublasToEigen`
- `find_all_sign_flips_in_rows`
- `find_sign_flip_nearest_middle`
- `flush_cpu_buffer_to_gpu_queue`
- `for`
- `format_Lbyas_for_filename`
- `gpu_consumer_thread_func`
- `gpu_process_chunk_oom_safe_recursive`
- `gpu_process_same_dim_chunk_cublas`
- `gpu_process_same_dim_same_vdim_group_cublas`
- `if`
- `linear_zero_estimate`
- `main`
- `make_scan_tag`
- `make_subchunk_copy`
- `pop`
- `push`
- `read_normalized_det_file_for_sign_flips`
- `release`
- `run_F3_energy_scan_gpu_omp_cublas_pipeline_v4`
- `test_F3_with_pwave_all_energy_cpu_openmp_v22_multi_signflip_side_abs_refinement`
- `test_F3_with_pwave_all_energy_gpu_omp_cublas_pipeline_v3`
- `test_F3_with_pwave_all_energy_gpu_omp_cublas_pipeline_v8_save_matrices`
- `while`

Inputs/outputs: inspect the file CLI and adjacent reports; this index is intentionally mechanical.

## `source/F3_cpu_openmp_v22_multi_signflip_side_abs_refinement_Lbyas_arg_compilefix2.cu`

Purpose: C++/header implementation.

Important functions/classes:
- `F3ScopedTimer`
- `add_unique_zero`
- `add_unique_zero_cpu_v21`
- `append_failed_results_for_chunk`
- `build_one_energy_pack`
- `bytes_complex_matrix`
- `bytes_to_mib`
- `catch`
- `classify_one_signflip_side_abs_trend_v22`
- `classify_refined_file_by_multi_side_abs_trend_v22`
- `classify_sign_flips_zero_vs_pole`
- `close`
- `ensure`
- `estimate_gpu_bytes_for_items`
- `estimate_gpu_bytes_one_energy`
- `f3_cudaSyncCheck`
- `f3_exception_looks_like_cuda_oom`
- `f3_has_sign_flip`
- `f3_make_auto_hdf5_matrix_filename`
- `f3_packEigenToHostCublas`
- `f3_should_save_matrix_point`
- `f3_sign_double`
- `f3_throwOnCublas`
- `f3_throwOnCuda`
- `f3_unpackHostCublasToEigen`
- `find_all_sign_flips_in_rows`
- `find_sign_flip_nearest_middle`
- `flush_cpu_buffer_to_gpu_queue`
- `for`
- `format_Lbyas_for_filename`
- `gpu_consumer_thread_func`
- `gpu_process_chunk_oom_safe_recursive`
- `gpu_process_same_dim_chunk_cublas`
- `gpu_process_same_dim_same_vdim_group_cublas`
- `if`
- `linear_zero_estimate`
- `main`
- `make_scan_tag`
- `make_subchunk_copy`
- `pop`
- `push`
- `read_normalized_det_file_for_sign_flips`
- `release`
- `run_F3_energy_scan_gpu_omp_cublas_pipeline_v4`
- `test_F3_with_pwave_all_energy_cpu_openmp_v22_multi_signflip_side_abs_refinement`
- `test_F3_with_pwave_all_energy_gpu_omp_cublas_pipeline_v3`
- `test_F3_with_pwave_all_energy_gpu_omp_cublas_pipeline_v8_save_matrices`
- `while`

Inputs/outputs: inspect the file CLI and adjacent reports; this index is intentionally mechanical.

## `source/F3_cpu_openmp_v23_pure_cpu_multi_signflip_side_abs_refinement.cpp`

Purpose: C++/header implementation.

Important functions/classes:
- `L`
- `L_tag`
- `LocalScopedTimer`
- `PhysicsParams`
- `add_unique_zero`
- `catch`
- `classify_one_subflip_side_abs`
- `determinant_via_partial_piv_lu`
- `evaluate_one_energy_cpu_eigen`
- `find_sign_flips`
- `for`
- `if`
- `main`
- `nnP_tag_string`
- `run_adaptive_cpu_only_v23`
- `run_scan_cpu_openmp`
- `trend_fraction`
- `while`

Inputs/outputs: inspect the file CLI and adjacent reports; this index is intentionally mechanical.

## `source/F3_cpu_openmp_v23_pure_cpu_multi_signflip_side_abs_refinement_compilefix1.cpp`

Purpose: C++/header implementation.

Important functions/classes:
- `L`
- `L_tag`
- `LocalScopedTimer`
- `PhysicsParams`
- `add_unique_zero`
- `catch`
- `classify_one_subflip_side_abs`
- `determinant_via_partial_piv_lu`
- `evaluate_one_energy_cpu_eigen`
- `find_sign_flips`
- `for`
- `if`
- `main`
- `nnP_tag_string`
- `run_adaptive_cpu_only_v23`
- `run_scan_cpu_openmp`
- `trend_fraction`
- `while`

Inputs/outputs: inspect the file CLI and adjacent reports; this index is intentionally mechanical.

## `source/F3_cpu_openmp_v24_K3QC_core.hpp`

Purpose: C++/header implementation.

Important functions/classes:
- `L`
- `L_tag`
- `LocalScopedTimer`
- `PhysicsParams`
- `add_unique_zero`
- `catch`
- `classify_one_subflip_side_abs`
- `determinant_via_partial_piv_lu`
- `evaluate_one_energy_cpu_eigen_v24_QC`
- `find_likely_zeros_v24_K3QC_silent`
- `find_sign_flips`
- `for`
- `if`
- `nnP_tag_string`
- `run_scan_cpu_openmp_silent`
- `trend_fraction`
- `while`

Inputs/outputs: inspect the file CLI and adjacent reports; this index is intentionally mechanical.

## `source/F3_cpu_openmp_v24_K3_QC_multi_signflip_side_abs_refinement.cpp`

Purpose: C++/header implementation.

Important functions/classes:
- `L`
- `L_tag`
- `LocalScopedTimer`
- `PhysicsParams`
- `add_unique_zero`
- `catch`
- `classify_one_subflip_side_abs`
- `determinant_via_partial_piv_lu`
- `evaluate_one_energy_cpu_eigen_v24_QC`
- `find_sign_flips`
- `for`
- `if`
- `main`
- `nnP_tag_string`
- `run_adaptive_cpu_only_v24_K3QC`
- `run_scan_cpu_openmp`
- `trend_fraction`
- `while`

Inputs/outputs: inspect the file CLI and adjacent reports; this index is intentionally mechanical.

## `source/F3_cpu_openmp_v25_K3QC_cached_core.hpp`

Purpose: C++/header implementation.

Important functions/classes:
- `L`
- `L_tag`
- `LocalScopedTimer`
- `PhysicsParams`
- `add_unique_zero`
- `build_F3i_projected_cache_entry_v25`
- `cache_key_v25`
- `catch`
- `classify_one_subflip_side_abs`
- `clear`
- `determinant_via_partial_piv_lu`
- `evaluate_one_energy_cpu_eigen_v24_QC`
- `evaluate_one_energy_cpu_eigen_v25_QC_cached`
- `find_likely_zeros_v24_K3QC_silent`
- `find_likely_zeros_v25_K3QC_cached`
- `find_sign_flips`
- `for`
- `if`
- `nnP_tag_string`
- `run_scan_cpu_openmp_cached`
- `run_scan_cpu_openmp_silent`
- `size`
- `trend_fraction`
- `while`

Inputs/outputs: inspect the file CLI and adjacent reports; this index is intentionally mechanical.

## `source/F3_device_resident_grouped_batched_lu_v13.cu`

Purpose: C++/header implementation.

Important functions/classes:
- `DevicePtr`
- `ScopedTimer`
- `adaptive_final_zeros_device_resident`
- `alloc`
- `build_F3_from_F2_and_F2X1_kernel`
- `build_projection_Vsel_device`
- `c_make`
- `classify_inward_shape`
- `combine_H_colmajor_kernel`
- `compute_one_energy_device_resident`
- `cublas_lu_solve_inplace_batched_same_n`
- `determinant_colmajor_device_matrix`
- `device_lu_solve_inplace`
- `for`
- `free`
- `gather_real_evecs_to_complex_Vsel_kernel`
- `get_energy_shape_prepass`
- `if`
- `linear_zero`
- `main`
- `make_scatter_params_1`
- `make_scatter_params_2`
- `process_same_shape_group_stream_parallel`
- `run_scan_device_resident`
- `sign_flip`
- `v11_set_stage`
- `v11_sigsegv_handler`
- `while`

Inputs/outputs: inspect the file CLI and adjacent reports; this index is intentionally mechanical.

## `source/F3_device_resident_grouped_stream_v12.cu`

Purpose: C++/header implementation.

Important functions/classes:
- `DevicePtr`
- `ScopedTimer`
- `adaptive_final_zeros_device_resident`
- `alloc`
- `build_F3_from_F2_and_F2X1_kernel`
- `build_projection_Vsel_device`
- `c_make`
- `classify_inward_shape`
- `combine_H_colmajor_kernel`
- `compute_one_energy_device_resident`
- `determinant_colmajor_device_matrix`
- `device_lu_solve_inplace`
- `for`
- `free`
- `gather_real_evecs_to_complex_Vsel_kernel`
- `get_energy_shape_prepass`
- `if`
- `linear_zero`
- `main`
- `make_scatter_params_1`
- `make_scatter_params_2`
- `process_same_shape_group_stream_parallel`
- `run_scan_device_resident`
- `sign_flip`
- `v11_set_stage`
- `v11_sigsegv_handler`
- `while`

Inputs/outputs: inspect the file CLI and adjacent reports; this index is intentionally mechanical.

## `source/F3_device_resident_speed_step_v11.cu`

Purpose: C++/header implementation.

Important functions/classes:
- `DevicePtr`
- `ScopedTimer`
- `adaptive_final_zeros_device_resident`
- `alloc`
- `build_F3_from_F2_and_F2X1_kernel`
- `build_projection_Vsel_device`
- `c_make`
- `classify_inward_shape`
- `combine_H_colmajor_kernel`
- `compute_one_energy_device_resident`
- `determinant_colmajor_device_matrix`
- `device_lu_solve_inplace`
- `for`
- `free`
- `gather_real_evecs_to_complex_Vsel_kernel`
- `if`
- `linear_zero`
- `main`
- `make_scatter_params_1`
- `make_scatter_params_2`
- `run_scan_device_resident`
- `sign_flip`
- `v11_set_stage`
- `v11_sigsegv_handler`

Inputs/outputs: inspect the file CLI and adjacent reports; this index is intentionally mechanical.

## `source/F3_device_resident_speed_step_v11_v2_compilefix.cu`

Purpose: C++/header implementation.

Important functions/classes:
- `DevicePtr`
- `ScopedTimer`
- `adaptive_final_zeros_device_resident`
- `alloc`
- `build_F3_from_F2_and_F2X1_kernel`
- `build_projection_Vsel_device`
- `c_make`
- `classify_inward_shape`
- `combine_H_colmajor_kernel`
- `compute_one_energy_device_resident`
- `determinant_colmajor_device_matrix`
- `device_lu_solve_inplace`
- `for`
- `free`
- `gather_real_evecs_to_complex_Vsel_kernel`
- `if`
- `linear_zero`
- `main`
- `make_scatter_params_1`
- `make_scatter_params_2`
- `run_scan_device_resident`
- `sign_flip`

Inputs/outputs: inspect the file CLI and adjacent reports; this index is intentionally mechanical.

## `source/F3_device_resident_speed_step_v11_v3_guarded.cu`

Purpose: C++/header implementation.

Important functions/classes:
- `DevicePtr`
- `ScopedTimer`
- `adaptive_final_zeros_device_resident`
- `alloc`
- `build_F3_from_F2_and_F2X1_kernel`
- `build_projection_Vsel_device`
- `c_make`
- `classify_inward_shape`
- `combine_H_colmajor_kernel`
- `compute_one_energy_device_resident`
- `determinant_colmajor_device_matrix`
- `device_lu_solve_inplace`
- `for`
- `free`
- `gather_real_evecs_to_complex_Vsel_kernel`
- `if`
- `linear_zero`
- `main`
- `make_scatter_params_1`
- `make_scatter_params_2`
- `run_scan_device_resident`
- `sign_flip`
- `v11_set_stage`
- `v11_sigsegv_handler`

Inputs/outputs: inspect the file CLI and adjacent reports; this index is intentionally mechanical.

## `source/F3_device_resident_speed_step_v11_v4_cusolver_lu.cu`

Purpose: C++/header implementation.

Important functions/classes:
- `DevicePtr`
- `ScopedTimer`
- `adaptive_final_zeros_device_resident`
- `alloc`
- `build_F3_from_F2_and_F2X1_kernel`
- `build_projection_Vsel_device`
- `c_make`
- `classify_inward_shape`
- `combine_H_colmajor_kernel`
- `compute_one_energy_device_resident`
- `determinant_colmajor_device_matrix`
- `device_lu_solve_inplace`
- `for`
- `free`
- `gather_real_evecs_to_complex_Vsel_kernel`
- `if`
- `linear_zero`
- `main`
- `make_scatter_params_1`
- `make_scatter_params_2`
- `run_scan_device_resident`
- `sign_flip`
- `v11_set_stage`
- `v11_sigsegv_handler`

Inputs/outputs: inspect the file CLI and adjacent reports; this index is intentionally mechanical.

## `source/F3_device_resident_true_batch_v14.cu`

Purpose: C++/header implementation.

Important functions/classes:
- `DevicePtr`
- `ScopedTimer`
- `adaptive_final_zeros_device_resident`
- `alloc`
- `build_F3_from_F2_and_F2X1_kernel`
- `build_projection_Vsel_device`
- `c_make`
- `classify_inward_shape`
- `combine_H_colmajor_kernel`
- `compute_one_energy_device_resident`
- `cublas_lu_solve_inplace_batched_same_n`
- `determinant_colmajor_device_matrix`
- `device_lu_solve_inplace`
- `for`
- `free`
- `gather_real_evecs_to_complex_Vsel_kernel`
- `get_energy_shape_prepass`
- `if`
- `linear_zero`
- `main`
- `make_scatter_params_1`
- `make_scatter_params_2`
- `process_same_shape_group_stream_parallel`
- `run_scan_device_resident`
- `sign_flip`
- `v11_set_stage`
- `v11_sigsegv_handler`
- `v14_build_F2_batch_same_shape_kernel`
- `v14_build_F3_batch_kernel`
- `v14_build_G_batch_same_shape_kernel`
- `v14_build_K2_batch_same_shape_kernel`
- `v14_combine_H_batch_kernel`
- `v14_make_ptrs_kernel`
- `while`

Inputs/outputs: inspect the file CLI and adjacent reports; this index is intentionally mechanical.

## `source/F3_device_resident_true_batch_v15_force_batch_diagnostics.cu`

Purpose: C++/header implementation.

Important functions/classes:
- `DevicePtr`
- `ScopedTimer`
- `adaptive_final_zeros_device_resident`
- `alloc`
- `build_F3_from_F2_and_F2X1_kernel`
- `build_projection_Vsel_device`
- `c_make`
- `catch`
- `classify_inward_shape`
- `combine_H_colmajor_kernel`
- `compute_one_energy_device_resident`
- `cublas_lu_solve_inplace_batched_same_n`
- `determinant_colmajor_device_matrix`
- `device_lu_solve_inplace`
- `for`
- `free`
- `gather_real_evecs_to_complex_Vsel_kernel`
- `get_energy_shape_prepass`
- `if`
- `linear_zero`
- `main`
- `make_scatter_params_1`
- `make_scatter_params_2`
- `process_same_shape_group_stream_parallel`
- `process_same_shape_group_true_batch_v15`
- `run_scan_device_resident`
- `sign_flip`
- `v11_set_stage`
- `v11_sigsegv_handler`
- `v14_build_F2_batch_same_shape_kernel`
- `v14_build_F3_batch_kernel`
- `v14_build_G_batch_same_shape_kernel`
- `v14_build_K2_batch_same_shape_kernel`
- `v14_combine_H_batch_kernel`
- `v14_make_ptrs_kernel`
- `while`

Inputs/outputs: inspect the file CLI and adjacent reports; this index is intentionally mechanical.

## `source/F3_device_resident_true_batch_v16_custom_batched_getrs.cu`

Purpose: C++/header implementation.

Important functions/classes:
- `DevicePtr`
- `ScopedTimer`
- `adaptive_final_zeros_device_resident`
- `alloc`
- `build_F3_from_F2_and_F2X1_kernel`
- `build_projection_Vsel_device`
- `c_make`
- `catch`
- `classify_inward_shape`
- `combine_H_colmajor_kernel`
- `compute_one_energy_device_resident`
- `cublas_lu_solve_inplace_batched_same_n`
- `determinant_colmajor_device_matrix`
- `device_lu_solve_inplace`
- `for`
- `free`
- `gather_real_evecs_to_complex_Vsel_kernel`
- `get_energy_shape_prepass`
- `if`
- `linear_zero`
- `main`
- `make_scatter_params_1`
- `make_scatter_params_2`
- `process_same_shape_group_stream_parallel`
- `process_same_shape_group_true_batch_v15`
- `run_scan_device_resident`
- `sign_flip`
- `v11_set_stage`
- `v11_sigsegv_handler`
- `v14_build_F2_batch_same_shape_kernel`
- `v14_build_F3_batch_kernel`
- `v14_build_G_batch_same_shape_kernel`
- `v14_build_K2_batch_same_shape_kernel`
- `v14_combine_H_batch_kernel`
- `v14_make_ptrs_kernel`
- `v16_getrs_from_getrf_batched_kernel`
- `v16_zadd`
- `v16_zdiv`
- `v16_zmul`
- `v16_zsub`
- `while`

Inputs/outputs: inspect the file CLI and adjacent reports; this index is intentionally mechanical.

## `source/F3_device_resident_true_batch_v17_vsel_reuse_chunk512_refine200.cu`

Purpose: C++/header implementation.

Important functions/classes:
- `DevicePtr`
- `ScopedTimer`
- `adaptive_final_zeros_device_resident`
- `alloc`
- `build_F3_from_F2_and_F2X1_kernel`
- `build_projection_Vsel_device`
- `c_make`
- `catch`
- `classify_inward_shape`
- `combine_H_colmajor_kernel`
- `compute_one_energy_device_resident`
- `cublas_lu_solve_inplace_batched_same_n`
- `determinant_colmajor_device_matrix`
- `device_lu_solve_inplace`
- `for`
- `free`
- `gather_real_evecs_to_complex_Vsel_kernel`
- `get_energy_shape_prepass`
- `if`
- `linear_zero`
- `main`
- `make_scatter_params_1`
- `make_scatter_params_2`
- `process_same_shape_group_stream_parallel`
- `process_same_shape_group_true_batch_v15`
- `run_scan_device_resident`
- `sign_flip`
- `v11_set_stage`
- `v11_sigsegv_handler`
- `v14_build_F2_batch_same_shape_kernel`
- `v14_build_F3_batch_kernel`
- `v14_build_G_batch_same_shape_kernel`
- `v14_build_K2_batch_same_shape_kernel`
- `v14_combine_H_batch_kernel`
- `v14_make_ptrs_kernel`
- `v16_getrs_from_getrf_batched_kernel`
- `v16_zadd`
- `v16_zdiv`
- `v16_zmul`
- `v16_zsub`
- `while`

Inputs/outputs: inspect the file CLI and adjacent reports; this index is intentionally mechanical.

## `source/F3_device_resident_true_batch_v18_zoom_refine_constantN.cu`

Purpose: C++/header implementation.

Important functions/classes:
- `DevicePtr`
- `ScopedTimer`
- `adaptive_final_zeros_device_resident`
- `alloc`
- `build_F3_from_F2_and_F2X1_kernel`
- `build_projection_Vsel_device`
- `c_make`
- `catch`
- `classify_inward_shape`
- `classify_inward_shape_result`
- `combine_H_colmajor_kernel`
- `compute_one_energy_device_resident`
- `cublas_lu_solve_inplace_batched_same_n`
- `determinant_colmajor_device_matrix`
- `device_lu_solve_inplace`
- `for`
- `free`
- `gather_real_evecs_to_complex_Vsel_kernel`
- `get_energy_shape_prepass`
- `if`
- `linear_zero`
- `main`
- `make_scatter_params_1`
- `make_scatter_params_2`
- `process_same_shape_group_stream_parallel`
- `process_same_shape_group_true_batch_v15`
- `run_scan_device_resident`
- `sign_flip`
- `v11_set_stage`
- `v11_sigsegv_handler`
- `v14_build_F2_batch_same_shape_kernel`
- `v14_build_F3_batch_kernel`
- `v14_build_G_batch_same_shape_kernel`
- `v14_build_K2_batch_same_shape_kernel`
- `v14_combine_H_batch_kernel`
- `v14_make_ptrs_kernel`
- `v16_getrs_from_getrf_batched_kernel`
- `v16_zadd`
- `v16_zdiv`
- `v16_zmul`
- `v16_zsub`
- `while`

Inputs/outputs: inspect the file CLI and adjacent reports; this index is intentionally mechanical.

## `source/F3_device_resident_true_batch_v19_parityfix_v9match.cu`

Purpose: C++/header implementation.

Important functions/classes:
- `DevicePtr`
- `ScopedTimer`
- `adaptive_final_zeros_device_resident`
- `alloc`
- `build_F3_from_F2_and_F2X1_kernel`
- `build_projection_Vsel_device`
- `c_make`
- `catch`
- `choose_projection_parity_v19`
- `classify_inward_shape`
- `classify_inward_shape_result`
- `combine_H_colmajor_kernel`
- `compute_one_energy_device_resident`
- `cublas_lu_solve_inplace_batched_same_n`
- `determinant_colmajor_device_matrix`
- `device_lu_solve_inplace`
- `for`
- `free`
- `gather_real_evecs_to_complex_Vsel_kernel`
- `get_energy_shape_prepass`
- `if`
- `linear_zero`
- `main`
- `make_scatter_params_1`
- `make_scatter_params_2`
- `process_same_shape_group_stream_parallel`
- `process_same_shape_group_true_batch_v15`
- `run_scan_device_resident`
- `sign_flip`
- `v11_set_stage`
- `v11_sigsegv_handler`
- `v14_build_F2_batch_same_shape_kernel`
- `v14_build_F3_batch_kernel`
- `v14_build_G_batch_same_shape_kernel`
- `v14_build_K2_batch_same_shape_kernel`
- `v14_combine_H_batch_kernel`
- `v14_make_ptrs_kernel`
- `v16_getrs_from_getrf_batched_kernel`
- `v16_zadd`
- `v16_zdiv`
- `v16_zmul`
- `v16_zsub`
- `while`

Inputs/outputs: inspect the file CLI and adjacent reports; this index is intentionally mechanical.

## `source/F3_device_resident_true_batch_v20_exact_signflip_refine.cu`

Purpose: C++/header implementation.

Important functions/classes:
- `DevicePtr`
- `ScopedTimer`
- `adaptive_final_zeros_device_resident`
- `alloc`
- `build_F3_from_F2_and_F2X1_kernel`
- `build_projection_Vsel_device`
- `c_make`
- `catch`
- `choose_projection_parity_v20`
- `classify_inward_shape`
- `classify_inward_shape_result`
- `combine_H_colmajor_kernel`
- `compute_one_energy_device_resident`
- `cublas_lu_solve_inplace_batched_same_n`
- `determinant_colmajor_device_matrix`
- `device_lu_solve_inplace`
- `for`
- `free`
- `gather_real_evecs_to_complex_Vsel_kernel`
- `get_energy_shape_prepass`
- `if`
- `linear_zero`
- `main`
- `make_scatter_params_1`
- `make_scatter_params_2`
- `process_same_shape_group_stream_parallel`
- `process_same_shape_group_true_batch_v15`
- `run_scan_device_resident`
- `sign_flip`
- `v11_set_stage`
- `v11_sigsegv_handler`
- `v14_build_F2_batch_same_shape_kernel`
- `v14_build_F3_batch_kernel`
- `v14_build_G_batch_same_shape_kernel`
- `v14_build_K2_batch_same_shape_kernel`
- `v14_combine_H_batch_kernel`
- `v14_make_ptrs_kernel`
- `v16_getrs_from_getrf_batched_kernel`
- `v16_zadd`
- `v16_zdiv`
- `v16_zmul`
- `v16_zsub`
- `while`

Inputs/outputs: inspect the file CLI and adjacent reports; this index is intentionally mechanical.

## `source/F3_device_resident_true_batch_v21_side_abs_trend_refine.cu`

Purpose: C++/header implementation.

Important functions/classes:
- `DevicePtr`
- `ScopedTimer`
- `adaptive_final_zeros_device_resident`
- `alloc`
- `build_F3_from_F2_and_F2X1_kernel`
- `build_projection_Vsel_device`
- `c_make`
- `catch`
- `choose_projection_parity_v20`
- `classify_inward_shape`
- `classify_inward_shape_result`
- `combine_H_colmajor_kernel`
- `compute_one_energy_device_resident`
- `cublas_lu_solve_inplace_batched_same_n`
- `determinant_colmajor_device_matrix`
- `device_lu_solve_inplace`
- `for`
- `free`
- `gather_real_evecs_to_complex_Vsel_kernel`
- `get_energy_shape_prepass`
- `if`
- `linear_zero`
- `main`
- `make_scatter_params_1`
- `make_scatter_params_2`
- `process_same_shape_group_stream_parallel`
- `process_same_shape_group_true_batch_v15`
- `run_scan_device_resident`
- `sign_flip`
- `v11_set_stage`
- `v11_sigsegv_handler`
- `v14_build_F2_batch_same_shape_kernel`
- `v14_build_F3_batch_kernel`
- `v14_build_G_batch_same_shape_kernel`
- `v14_build_K2_batch_same_shape_kernel`
- `v14_combine_H_batch_kernel`
- `v14_make_ptrs_kernel`
- `v16_getrs_from_getrf_batched_kernel`
- `v16_zadd`
- `v16_zdiv`
- `v16_zmul`
- `v16_zsub`
- `while`

Inputs/outputs: inspect the file CLI and adjacent reports; this index is intentionally mechanical.

## `source/F3_gpu_full_matrix_builder_v10_adaptive_zeros.cu`

Purpose: C++/header implementation.

Important functions/classes:
- `F3ScopedTimer`
- `add_unique_zero`
- `append_failed_results_for_chunk`
- `build_one_energy_pack`
- `bytes_complex_matrix`
- `bytes_to_mib`
- `catch`
- `classify_refined_file_by_inward_shape`
- `classify_sign_flips_zero_vs_pole`
- `close`
- `ensure`
- `estimate_gpu_bytes_for_items`
- `estimate_gpu_bytes_one_energy`
- `f3_cudaSyncCheck`
- `f3_exception_looks_like_cuda_oom`
- `f3_has_sign_flip`
- `f3_make_auto_hdf5_matrix_filename`
- `f3_packEigenToHostCublas`
- `f3_should_save_matrix_point`
- `f3_sign_double`
- `f3_throwOnCublas`
- `f3_throwOnCuda`
- `f3_unpackHostCublasToEigen`
- `find_all_sign_flips_in_rows`
- `find_sign_flip_nearest_middle`
- `flush_cpu_buffer_to_gpu_queue`
- `for`
- `gpu_consumer_thread_func`
- `gpu_process_chunk_oom_safe_recursive`
- `gpu_process_same_dim_chunk_cublas`
- `gpu_process_same_dim_same_vdim_group_cublas`
- `if`
- `linear_zero_estimate`
- `main`
- `make_scan_tag`
- `make_subchunk_copy`
- `normalize_det_vector_by_max`
- `pop`
- `push`
- `read_normalized_det_file_for_sign_flips`
- `release`
- `run_F3_energy_scan_gpu_omp_cublas_pipeline_v4`
- `smallest_eigenvalue`
- `test_F3_with_pwave_all_energy_gpu_omp_cublas_adaptive_zero_search_v1`
- `test_F3_with_pwave_all_energy_gpu_omp_cublas_pipeline_v3`
- `test_F3_with_pwave_all_energy_gpu_omp_cublas_pipeline_v8_save_matrices`
- `while`

Inputs/outputs: inspect the file CLI and adjacent reports; this index is intentionally mechanical.

## `source/F3_gpu_full_matrix_builder_v10_adaptive_zeros_v2_compilefix.cu`

Purpose: C++/header implementation.

Important functions/classes:
- `F3ScopedTimer`
- `add_unique_zero`
- `append_failed_results_for_chunk`
- `build_one_energy_pack`
- `bytes_complex_matrix`
- `bytes_to_mib`
- `catch`
- `classify_refined_file_by_inward_shape`
- `classify_sign_flips_zero_vs_pole`
- `close`
- `ensure`
- `estimate_gpu_bytes_for_items`
- `estimate_gpu_bytes_one_energy`
- `f3_cudaSyncCheck`
- `f3_exception_looks_like_cuda_oom`
- `f3_has_sign_flip`
- `f3_make_auto_hdf5_matrix_filename`
- `f3_packEigenToHostCublas`
- `f3_should_save_matrix_point`
- `f3_sign_double`
- `f3_throwOnCublas`
- `f3_throwOnCuda`
- `f3_unpackHostCublasToEigen`
- `find_all_sign_flips_in_rows`
- `find_sign_flip_nearest_middle`
- `flush_cpu_buffer_to_gpu_queue`
- `for`
- `gpu_consumer_thread_func`
- `gpu_process_chunk_oom_safe_recursive`
- `gpu_process_same_dim_chunk_cublas`
- `gpu_process_same_dim_same_vdim_group_cublas`
- `if`
- `linear_zero_estimate`
- `main`
- `make_scan_tag`
- `make_subchunk_copy`
- `normalize_det_vector_by_max`
- `pop`
- `push`
- `read_normalized_det_file_for_sign_flips`
- `release`
- `run_F3_energy_scan_gpu_omp_cublas_pipeline_v4`
- `smallest_eigenvalue`
- `test_F3_with_pwave_all_energy_gpu_omp_cublas_adaptive_zero_search_v1`
- `test_F3_with_pwave_all_energy_gpu_omp_cublas_pipeline_v3`
- `test_F3_with_pwave_all_energy_gpu_omp_cublas_pipeline_v8_save_matrices`
- `while`

Inputs/outputs: inspect the file CLI and adjacent reports; this index is intentionally mechanical.

## `source/F3_gpu_omp_cublas_pipeline_v3.cu`

Purpose: C++/header implementation.

Important functions/classes:
- `ScopedTimer`
- `build_f3_kernel`
- `build_one_energy_pack`
- `build_pointer_array_kernel`
- `bytes_complex_matrix`
- `catch`
- `close`
- `compute_det_from_lu_kernel`
- `cudaSyncCheck`
- `ensure`
- `estimate_gpu_bytes_one_energy`
- `flush_cpu_buffer_to_gpu_queue`
- `for`
- `gpu_consumer_thread_func`
- `gpu_process_same_dim_chunk_cublas`
- `gpu_process_same_dim_same_vdim_group_cublas`
- `if`
- `main`
- `packEigenToHostCublas`
- `pop`
- `push`
- `release`
- `test_F3_with_pwave_all_energy_gpu_omp_cublas_pipeline_v3`
- `throwOnCublas`
- `throwOnCuda`
- `unpackHostCublasToEigen`
- `while`

Inputs/outputs: inspect the file CLI and adjacent reports; this index is intentionally mechanical.

## `source/F3_gpu_omp_cublas_pipeline_v4_include_safe.cu`

Purpose: C++/header implementation.

Important functions/classes:
- `F3ScopedTimer`
- `build_one_energy_pack`
- `bytes_complex_matrix`
- `catch`
- `close`
- `ensure`
- `estimate_gpu_bytes_one_energy`
- `f3_cudaSyncCheck`
- `f3_packEigenToHostCublas`
- `f3_throwOnCublas`
- `f3_throwOnCuda`
- `f3_unpackHostCublasToEigen`
- `flush_cpu_buffer_to_gpu_queue`
- `for`
- `gpu_consumer_thread_func`
- `gpu_process_same_dim_chunk_cublas`
- `gpu_process_same_dim_same_vdim_group_cublas`
- `if`
- `main`
- `pop`
- `push`
- `release`
- `test_F3_with_pwave_all_energy_gpu_omp_cublas_pipeline_v3`
- `while`

Inputs/outputs: inspect the file CLI and adjacent reports; this index is intentionally mechanical.

## `source/F3_gpu_omp_cublas_pipeline_v6_oom_safe_chunks.cu`

Purpose: C++/header implementation.

Important functions/classes:
- `F3ScopedTimer`
- `append_failed_results_for_chunk`
- `build_one_energy_pack`
- `bytes_complex_matrix`
- `bytes_to_mib`
- `catch`
- `close`
- `ensure`
- `estimate_gpu_bytes_for_items`
- `estimate_gpu_bytes_one_energy`
- `f3_cudaSyncCheck`
- `f3_exception_looks_like_cuda_oom`
- `f3_packEigenToHostCublas`
- `f3_throwOnCublas`
- `f3_throwOnCuda`
- `f3_unpackHostCublasToEigen`
- `flush_cpu_buffer_to_gpu_queue`
- `for`
- `gpu_consumer_thread_func`
- `gpu_process_chunk_oom_safe_recursive`
- `gpu_process_same_dim_chunk_cublas`
- `gpu_process_same_dim_same_vdim_group_cublas`
- `if`
- `main`
- `make_subchunk_copy`
- `pop`
- `push`
- `release`
- `test_F3_with_pwave_all_energy_gpu_omp_cublas_pipeline_v3`
- `while`

Inputs/outputs: inspect the file CLI and adjacent reports; this index is intentionally mechanical.

## `source/F3_gpu_omp_cublas_pipeline_v7_adaptive_zero_search.cu`

Purpose: C++/header implementation.

Important functions/classes:
- `F3ScopedTimer`
- `add_unique_zero`
- `append_failed_results_for_chunk`
- `build_one_energy_pack`
- `bytes_complex_matrix`
- `bytes_to_mib`
- `catch`
- `classify_sign_flips_zero_vs_pole`
- `close`
- `ensure`
- `estimate_gpu_bytes_for_items`
- `estimate_gpu_bytes_one_energy`
- `f3_cudaSyncCheck`
- `f3_exception_looks_like_cuda_oom`
- `f3_packEigenToHostCublas`
- `f3_sign_double`
- `f3_throwOnCublas`
- `f3_throwOnCuda`
- `f3_unpackHostCublasToEigen`
- `flush_cpu_buffer_to_gpu_queue`
- `for`
- `gpu_consumer_thread_func`
- `gpu_process_chunk_oom_safe_recursive`
- `gpu_process_same_dim_chunk_cublas`
- `gpu_process_same_dim_same_vdim_group_cublas`
- `if`
- `linear_zero_estimate`
- `main`
- `make_scan_tag`
- `make_subchunk_copy`
- `pop`
- `push`
- `read_normalized_det_file_for_sign_flips`
- `release`
- `run_F3_energy_scan_gpu_omp_cublas_pipeline_v4`
- `test_F3_with_pwave_all_energy_gpu_omp_cublas_adaptive_zero_search_v1`
- `test_F3_with_pwave_all_energy_gpu_omp_cublas_pipeline_v3`
- `while`

Inputs/outputs: inspect the file CLI and adjacent reports; this index is intentionally mechanical.

## `source/F3_gpu_omp_cublas_pipeline_v8_adaptive_matrix_saving.cu`

Purpose: C++/header implementation.

Important functions/classes:
- `F3ScopedTimer`
- `add_unique_zero`
- `append_failed_results_for_chunk`
- `build_one_energy_pack`
- `bytes_complex_matrix`
- `bytes_to_mib`
- `catch`
- `classify_sign_flips_zero_vs_pole`
- `close`
- `ensure`
- `estimate_gpu_bytes_for_items`
- `estimate_gpu_bytes_one_energy`
- `f3_cudaSyncCheck`
- `f3_exception_looks_like_cuda_oom`
- `f3_make_auto_hdf5_matrix_filename`
- `f3_packEigenToHostCublas`
- `f3_should_save_matrix_point`
- `f3_sign_double`
- `f3_throwOnCublas`
- `f3_throwOnCuda`
- `f3_unpackHostCublasToEigen`
- `flush_cpu_buffer_to_gpu_queue`
- `for`
- `gpu_consumer_thread_func`
- `gpu_process_chunk_oom_safe_recursive`
- `gpu_process_same_dim_chunk_cublas`
- `gpu_process_same_dim_same_vdim_group_cublas`
- `if`
- `linear_zero_estimate`
- `main`
- `make_scan_tag`
- `make_subchunk_copy`
- `pop`
- `push`
- `read_normalized_det_file_for_sign_flips`
- `release`
- `run_F3_energy_scan_gpu_omp_cublas_pipeline_v4`
- `test_F3_with_pwave_all_energy_gpu_omp_cublas_adaptive_zero_search_v1`
- `test_F3_with_pwave_all_energy_gpu_omp_cublas_pipeline_v3`
- `test_F3_with_pwave_all_energy_gpu_omp_cublas_pipeline_v8_save_matrices`
- `while`

Inputs/outputs: inspect the file CLI and adjacent reports; this index is intentionally mechanical.

## `source/F3_gpu_omp_cublas_pipeline_v8_matrix_saving.cu`

Purpose: C++/header implementation.

Important functions/classes:
- `F3ScopedTimer`
- `add_unique_zero`
- `append_failed_results_for_chunk`
- `build_one_energy_pack`
- `bytes_complex_matrix`
- `bytes_to_mib`
- `catch`
- `classify_sign_flips_zero_vs_pole`
- `close`
- `ensure`
- `estimate_gpu_bytes_for_items`
- `estimate_gpu_bytes_one_energy`
- `f3_cudaSyncCheck`
- `f3_exception_looks_like_cuda_oom`
- `f3_make_auto_hdf5_matrix_filename`
- `f3_packEigenToHostCublas`
- `f3_should_save_matrix_point`
- `f3_sign_double`
- `f3_throwOnCublas`
- `f3_throwOnCuda`
- `f3_unpackHostCublasToEigen`
- `flush_cpu_buffer_to_gpu_queue`
- `for`
- `gpu_consumer_thread_func`
- `gpu_process_chunk_oom_safe_recursive`
- `gpu_process_same_dim_chunk_cublas`
- `gpu_process_same_dim_same_vdim_group_cublas`
- `if`
- `linear_zero_estimate`
- `main`
- `make_scan_tag`
- `make_subchunk_copy`
- `pop`
- `push`
- `read_normalized_det_file_for_sign_flips`
- `release`
- `run_F3_energy_scan_gpu_omp_cublas_pipeline_v4`
- `test_F3_with_pwave_all_energy_gpu_omp_cublas_adaptive_zero_search_v1`
- `test_F3_with_pwave_all_energy_gpu_omp_cublas_pipeline_v3`
- `test_F3_with_pwave_all_energy_gpu_omp_cublas_pipeline_v8_save_matrices`
- `while`

Inputs/outputs: inspect the file CLI and adjacent reports; this index is intentionally mechanical.

## `source/F3_gpu_omp_cublas_pipeline_v9_inward_shape_refinement.cu`

Purpose: C++/header implementation.

Important functions/classes:
- `F3ScopedTimer`
- `add_unique_zero`
- `append_failed_results_for_chunk`
- `build_one_energy_pack`
- `bytes_complex_matrix`
- `bytes_to_mib`
- `catch`
- `classify_refined_file_by_inward_shape`
- `classify_sign_flips_zero_vs_pole`
- `close`
- `ensure`
- `estimate_gpu_bytes_for_items`
- `estimate_gpu_bytes_one_energy`
- `f3_cudaSyncCheck`
- `f3_exception_looks_like_cuda_oom`
- `f3_has_sign_flip`
- `f3_make_auto_hdf5_matrix_filename`
- `f3_packEigenToHostCublas`
- `f3_should_save_matrix_point`
- `f3_sign_double`
- `f3_throwOnCublas`
- `f3_throwOnCuda`
- `f3_unpackHostCublasToEigen`
- `find_all_sign_flips_in_rows`
- `find_sign_flip_nearest_middle`
- `flush_cpu_buffer_to_gpu_queue`
- `for`
- `gpu_consumer_thread_func`
- `gpu_process_chunk_oom_safe_recursive`
- `gpu_process_same_dim_chunk_cublas`
- `gpu_process_same_dim_same_vdim_group_cublas`
- `if`
- `linear_zero_estimate`
- `main`
- `make_scan_tag`
- `make_subchunk_copy`
- `pop`
- `push`
- `read_normalized_det_file_for_sign_flips`
- `release`
- `run_F3_energy_scan_gpu_omp_cublas_pipeline_v4`
- `test_F3_with_pwave_all_energy_gpu_omp_cublas_adaptive_zero_search_v1`
- `test_F3_with_pwave_all_energy_gpu_omp_cublas_pipeline_v3`
- `test_F3_with_pwave_all_energy_gpu_omp_cublas_pipeline_v8_save_matrices`
- `while`

Inputs/outputs: inspect the file CLI and adjacent reports; this index is intentionally mechanical.

## `source/F3_gpu_omp_pipeline_v2_eigenbased.cpp`

Purpose: C++/header implementation.

Important functions/classes:
- `ScopedTimer`
- `build_one_energy_pack`
- `bytes_complex_matrix`
- `catch`
- `close`
- `estimate_gpu_bytes_one_energy`
- `flush_cpu_buffer_to_gpu_queue`
- `for`
- `gpu_consumer_thread_func`
- `gpu_process_same_dim_chunk_cusolver`
- `if`
- `main`
- `pop`
- `push`
- `test_F3_with_pwave_all_energy_gpu_omp_pipeline_v2`
- `while`

Inputs/outputs: inspect the file CLI and adjacent reports; this index is intentionally mechanical.

## `source/F3_matrix_hdf5_saver.hpp`

Purpose: C++/header implementation.

Important functions/classes:
- `F3MatrixHdf5Saver`
- `catch`
- `create_or_open_group`
- `delete_link_if_exists`
- `epoint_group_name`
- `filename`
- `has_energy_point`
- `hdf5_link_exists`
- `if`
- `load_energy_point`
- `make_matrix_dump_filename`
- `open_file`
- `read_complex_matrix_dataset`
- `read_double_pair_dataset`
- `read_int_vector_dataset`
- `read_metadata_from_group`
- `read_scalar_dataset`
- `read_string_dataset`
- `safe_chunk_dim`
- `save_energy_point`
- `save_energy_point_minimal`
- `save_single_energy_matrices_to_hdf5`
- `variable_string_type`
- `write_complex_matrix_dataset`
- `write_double_pair_dataset`
- `write_global_metadata_if_missing`
- `write_int_vector_dataset`
- `write_metadata_to_group`
- `write_scalar_dataset`
- `write_string_dataset`

Inputs/outputs: inspect the file CLI and adjacent reports; this index is intentionally mechanical.

## `source/F3iso_maker.py`

Purpose: Python workflow/utility.

Important functions/classes:
- `E_to_Ecm`
- `Ecmsq_to_Esq`
- `Esq_to_Ecmsq`
- `P000`
- `P100`
- `P110`
- `P111`
- `P200`
- `jackknifeavg_centralvalue_lattice_data`
- `jackknifeavg_lattice_data`

Inputs/outputs: inspect the file CLI and adjacent reports; this index is intentionally mechanical.

## `source/Faddeeva.cc`

Purpose: C++/header implementation.

Important functions/classes:
- `C`
- `FADDEEVA`
- `FADDEEVA_RE`
- `complex`
- `dawson`
- `erf`
- `erfcx`
- `exp`
- `faster`
- `files`
- `if`
- `main`
- `my_copysign`
- `my_isinf`
- `my_isnan`
- `relerr`
- `return`
- `sinc`
- `sinh_taylor`
- `sqr`
- `sqrt`
- `while`

Inputs/outputs: inspect the file CLI and adjacent reports; this index is intentionally mechanical.

## `source/Faddeeva.h`

Purpose: C++/header implementation.

Important functions/classes:
- No function definitions detected by the lightweight indexer.

Inputs/outputs: inspect the file CLI and adjacent reports; this index is intentionally mechanical.

## `source/G_functions.h`

Purpose: C++/header implementation.

Important functions/classes:
- `G_ij`
- `G_ij_mat`
- `if`

Inputs/outputs: inspect the file CLI and adjacent reports; this index is intentionally mechanical.

## `source/G_functions_v2.h`

Purpose: C++/header implementation.

Important functions/classes:
- `G_2plus1_mat`
- `G_ij_lm`
- `get_parity_block_from_plm_config`
- `if`
- `print_Gij_boosts`

Inputs/outputs: inspect the file CLI and adjacent reports; this index is intentionally mechanical.

## `source/G_gpu_safe_builder.cuh`

Purpose: C++/header implementation.

Important functions/classes:
- `G_2plus1_kernel`
- `G_2plus1_many_mat_gpu_safe`
- `G_2plus1_mat_gpu_safe`
- `G_ij_lm_gpu`
- `Jfunc_gpu`
- `boost_gpu`
- `cabs`
- `cabs2`
- `cadd`
- `cdiv`
- `cdiv_d`
- `cexp_c`
- `cmul`
- `cmul_d`
- `cneg`
- `compare_matrices`
- `cpow_int`
- `csqrt_c`
- `csub`
- `cutoff_gpu`
- `flatten_config_host`
- `if`
- `kallen_gpu`
- `make_c`
- `omega_gpu`
- `q2psq_star_gpu`
- `sigma_pvec_gpu`
- `spherical_harmonics_gpu`
- `while`

Inputs/outputs: inspect the file CLI and adjacent reports; this index is intentionally mechanical.

## `source/K2_functions.h`

Purpose: C++/header implementation.

Important functions/classes:
- `K2_inv_00`
- `K2_inv_00_test_FRL`
- `K2inv_i_mat`
- `if`
- `tilde_K2_00`

Inputs/outputs: inspect the file CLI and adjacent reports; this index is intentionally mechanical.

## `source/K2_functions_gpu_safe.cuh`

Purpose: C++/header implementation.

Important functions/classes:
- `Cx`
- `Jfunc_gpu`
- `K2_inv_ERE_ang_mom_gpu`
- `K2inv_EREord2_2plus1_many_mat_gpu_safe`
- `K2inv_EREord2_2plus1_mat_gpu_safe`
- `K2inv_diag_element_gpu`
- `abs_cx`
- `build_K2inv_2plus1_kernel`
- `compare_matrices`
- `cutoff_function_1_gpu`
- `exp_cx`
- `get`
- `if`
- `kallen_gpu`
- `make_cx`
- `move_from`
- `norm2`
- `omega_func_gpu`
- `pow_int_cx`
- `q2psq_star_gpu`
- `release`
- `sqrt_cx`
- `upload`
- `view`

Inputs/outputs: inspect the file CLI and adjacent reports; this index is intentionally mechanical.

## `source/K2_functions_v2.h`

Purpose: C++/header implementation.

Important functions/classes:
- `K2_inv_00`
- `K2_inv_ERE_ang_mom`
- `K2inv_EREord2_2plus1_mat`
- `K2inv_EREord2_i_mat`
- `if`

Inputs/outputs: inspect the file CLI and adjacent reports; this index is intentionally mechanical.

## `source/K2_gpu_safe_builder.cuh`

Purpose: C++/header implementation.

Important functions/classes:
- `CudaComplex`
- `Jfunc_dev`
- `K2_inv_ERE_ang_mom_dev`
- `K2inv_EREord2_2plus1_many_mat_gpu_safe`
- `K2inv_EREord2_2plus1_mat_gpu_safe`
- `abs_c`
- `allocate_and_copy`
- `compare_matrices`
- `cutoff_function_1_dev`
- `exp_c`
- `flatten_config`
- `flatten_scatter_params`
- `if`
- `k2_fill_full_kernel`
- `kallen_dev`
- `make_c`
- `omega_func_dev`
- `pow_int_c`
- `q2psq_star_dev`
- `release`
- `sigma_pvec_based_dev`
- `sqrt_c`

Inputs/outputs: inspect the file CLI and adjacent reports; this index is intentionally mechanical.

## `source/K3_functions_2plus1.hpp`

Purpose: C++/header implementation.

Important functions/classes:
- `K3B_element`
- `K3E_element`
- `K3_element_raw`
- `K3mat_2plus1`
- `boost_cpp`
- `dot3`
- `ell_from_config`
- `flavor_mass`
- `from_stdvec`
- `if`
- `is_p`
- `is_s`
- `m_from_config`
- `make_vec3`
- `norm2`
- `norm3`
- `omega_vec`
- `qst2_i_cpp`
- `sigma_i_cpp`
- `third_mass_for_E_piece`
- `to_stdvec`
- `two_plus_one_block_factor`
- `vec_from_config`
- `y1real_cpp`

Inputs/outputs: inspect the file CLI and adjacent reports; this index is intentionally mechanical.

## `source/K3df_minuit_fit_v24.hpp`

Purpose: C++/header implementation.

Important functions/classes:
- `K3dfFCN`
- `canonical_shell_momentum`
- `catch`
- `chi_square_from_model`
- `covariance_to_correlation`
- `finite_difference_energy_errors_and_covariance`
- `fit_K3df_parameters_minuit_v24`
- `for`
- `if`
- `infer_row_spec_index_from_covariance_order`
- `make_physics_params_from_settings`
- `minuit_covariance_to_eigen`
- `parse_momentum_irrep_label`
- `print_fit_result_summary`
- `same_momentum`
- `same_momentum_shell`
- `solve_model_energies_for_rows`
- `spec_key`
- `symmetric_pseudoinverse`
- `while`

Inputs/outputs: inspect the file CLI and adjacent reports; this index is intentionally mechanical.

## `source/K3df_minuit_fit_v24_NEW_FIXED.hpp`

Purpose: C++/header implementation.

Important functions/classes:
- `K3dfFCN`
- `canonical_shell_momentum`
- `catch`
- `chi_square_from_model`
- `covariance_to_correlation`
- `finite_difference_energy_errors_and_covariance`
- `fit_K3df_parameters_minuit_v24`
- `for`
- `if`
- `infer_row_spec_index_from_covariance_order`
- `make_physics_params_from_settings`
- `minuit_covariance_to_eigen`
- `parse_momentum_irrep_label`
- `print_fit_result_summary`
- `same_momentum`
- `same_momentum_shell`
- `solve_model_energies_for_rows`
- `spec_key`
- `symmetric_pseudoinverse`
- `while`

Inputs/outputs: inspect the file CLI and adjacent reports; this index is intentionally mechanical.

## `source/K3df_minuit_fit_v25_cached.hpp`

Purpose: C++/header implementation.

Important functions/classes:
- `K3dfFCN`
- `canonical_shell_momentum`
- `catch`
- `chi_square_from_model`
- `covariance_to_correlation`
- `finite_difference_energy_errors_and_covariance`
- `fit_K3df_parameters_minuit_v24`
- `for`
- `if`
- `infer_row_spec_index_from_covariance_order`
- `make_physics_params_from_settings`
- `minuit_covariance_to_eigen`
- `parse_momentum_irrep_label`
- `print_fit_result_summary`
- `same_momentum`
- `same_momentum_shell`
- `solve_model_energies_for_rows`
- `solve_model_energies_for_rows_cached`
- `spec_key`
- `symmetric_pseudoinverse`
- `while`

Inputs/outputs: inspect the file CLI and adjacent reports; this index is intentionally mechanical.

## `source/K3df_minuit_fit_v26_reuse_cache.hpp`

Purpose: C++/header implementation.

Important functions/classes:
- `K3dfFCN`
- `canonical_shell_momentum`
- `catch`
- `chi_square_from_model`
- `covariance_to_correlation`
- `finite_difference_energy_errors_and_covariance`
- `fit_K3df_parameters_minuit_v24`
- `fit_K3df_parameters_minuit_v26`
- `for`
- `get_shared_cache`
- `if`
- `infer_row_spec_index_from_covariance_order`
- `make_output_path`
- `make_physics_params_from_settings`
- `minuit_covariance_to_eigen`
- `parse_momentum_irrep_label`
- `print_fit_result_summary`
- `same_momentum`
- `same_momentum_shell`
- `solve_model_energies_for_rows`
- `solve_model_energies_for_rows_cached`
- `spec_key`
- `symmetric_pseudoinverse`
- `while`
- `write_fit_outputs`
- `write_matrix_dat`

Inputs/outputs: inspect the file CLI and adjacent reports; this index is intentionally mechanical.

## `source/K3df_minuit_fit_v31l_lattice_covariance.hpp`

Purpose: C++/header implementation.

Important functions/classes:
- `K3dfFCN_v31l`
- `assemble_QC`
- `build_F3inv_only_cache`
- `build_cache_entry`
- `canonical_shell_momentum`
- `chi_square_v31l`
- `clean_label`
- `closest_zero_eigenvalue`
- `covariance_to_correlation_v31l`
- `diagonal_chi2`
- `find_all_zero_poles`
- `find_zero_poles_from_cache`
- `fit_K3df_parameters_v31l`
- `for`
- `get_or_build_F3inv_cache_v31l`
- `hermitize`
- `if`
- `infer_row_spec_index_from_covariance_order_v31l`
- `load_binary_f3inv_cache_v31l`
- `load_targets_and_covariance_v31l`
- `make_K3_projected`
- `make_base_physics`
- `minuit_covariance_to_eigen_v31l`
- `model_levels_for_targets`
- `operator`
- `point_value`
- `print_fit_result_summary`
- `read_binary_raw_v31l`
- `read_comp_v31l`
- `read_matrix_v31l`
- `read_nested_vec_comp_v31l`
- `read_string_v31l`
- `read_target_levels`
- `read_vec_comp_v31l`
- `same_momentum_shell_v31l`
- `save_binary_f3inv_cache_v31l`
- `symmetric_pseudoinverse_v31l`
- `while`
- `write_binary_raw_v31l`
- `write_cache_grid_files`
- `write_comp_v31l`
- `write_fit_summary`
- `write_matrix_v31l`
- `write_nested_vec_comp_v31l`
- `write_qc_eig_grid_files`
- `write_string_v31l`
- `write_vec_comp_v31l`

Inputs/outputs: inspect the file CLI and adjacent reports; this index is intentionally mechanical.

## `source/K3df_minuit_fit_v32f_fullF3inv_QCfull_cached_classifier.hpp`

Purpose: C++/header implementation.

Important functions/classes:
- `K3dfFCN_v32f`
- `assemble_QC`
- `build_F3inv_only_cache`
- `build_cache_entry`
- `canonical_shell_momentum`
- `chi_square_v32f`
- `clean_label`
- `closest_zero_eigenvalue`
- `covariance_to_correlation_v32f`
- `det_real_projected_F3inv_v32f`
- `diagonal_chi2`
- `find_all_zero_poles`
- `find_zero_poles_from_cache`
- `fit_K3df_parameters_v32f`
- `for`
- `get_or_build_F3inv_cache_v32f`
- `hermitize`
- `if`
- `infer_row_spec_index_from_covariance_order_v32f`
- `linear_root_v32f`
- `load_binary_f3inv_cache_v32f`
- `load_targets_and_covariance_v32f`
- `make_K3_projected`
- `make_base_physics`
- `minuit_covariance_to_eigen_v32f`
- `model_levels_for_targets`
- `operator`
- `point_value`
- `precompute_projected_k3_basis`
- `print_fit_result_summary`
- `read_binary_raw_v32f`
- `read_comp_v32f`
- `read_matrix_v32f`
- `read_nested_vec_comp_v32f`
- `read_string_v32f`
- `read_target_levels`
- `read_vec_comp_v32f`
- `refine_F3inv_det_signflips_v32f`
- `resolve_existing_binary_f3inv_cache_v32k`
- `same_momentum_shell_v32f`
- `save_binary_f3inv_cache_v32f`
- `signflip_v32f`
- `sort_unique_cache_grid_v32f`
- `symmetric_pseudoinverse_v32f`
- `while`
- `write_binary_raw_v32f`
- `write_cache_grid_files`
- `write_comp_v32f`
- `write_fit_summary`
- `write_matrix_v32f`
- `write_nested_vec_comp_v32f`
- `write_qc_eig_grid_files`
- `write_string_v32f`
- `write_vec_comp_v32f`

Inputs/outputs: inspect the file CLI and adjacent reports; this index is intentionally mechanical.

## `source/K3iso_maker.py`

Purpose: Python workflow/utility.

Important functions/classes:
- `E_to_Ecm`
- `Ecmsq_to_Esq`
- `Esq_to_Ecmsq`
- `P000`
- `P100`
- `P110`
- `P111`
- `P200`
- `jackknifeavg_centralvalue_lattice_data`
- `jackknifeavg_centralvalue_lattice_data_L20_L24`
- `jackknifeavg_lattice_data`
- `jackknifeavg_lattice_data_L20_L24`
- `print`

Inputs/outputs: inspect the file CLI and adjacent reports; this index is intentionally mechanical.

## `source/QC_fitter_1.py`

Purpose: Python workflow/utility.

Important functions/classes:
- `E_to_Ecm`
- `Ecmsq_to_Esq`
- `Esq_to_Ecmsq`
- `K3iso_fitting_function`
- `K3iso_fitting_function_all_moms_one_parameter`
- `K3iso_fitting_function_all_moms_two_parameter`
- `K3iso_fitting_function_all_moms_two_parameter_interp1d_based`
- `K3iso_fitting_function_all_moms_two_parameter_secant`
- `K3iso_fitting_function_multiL_oneparameter_interp1d_based`
- `K3iso_fitting_function_multiL_twoparameter_interp1d_based`
- `K3iso_fitting_function_multiL_twoparameter_interp1d_based_K3iso0_fixed`
- `QC3`
- `QC3_bissection_eigen_based`
- `QC3_bissection_interp1d_based`
- `QC3_bissection_interp1d_based_multiL`
- `QC3_bissection_spline_based`
- `QC3_secant_eigen_based`
- `QC_spectrum_one_parameter`
- `QC_spectrum_one_parameter_multiLs_multiK3df`
- `QC_spectrum_two_parameter`
- `QC_spectrum_two_parameter_multiLs`
- `QC_spectrum_two_parameter_multiLs_multiK3df`
- `QC_states_jackknife_resampler_oneparam`
- `QC_states_jackknife_resampler_twoparams`
- `if`
- `main`
- `one_param_fitting`
- `one_param_fitting_state0`
- `one_param_fitting_state0_000_A1m`
- `print`
- `sign_func`
- `spectrum_checker_for_QC`
- `spectrum_checker_for_splines`
- `test`
- `test1`
- `test1_K3df_fitting_and_QC_state_generator_one_param`
- `test1_K3df_fitting_and_QC_state_generator_two_params`
- `test1_K3df_fitting_twoLs_one_param`
- `test1_K3df_fitting_twoLs_one_param_state0_000_A1m_only`
- `test1_K3df_fitting_twoLs_one_param_state0_only`
- `test1_K3df_fitting_twoLs_two_params`
- `test1_K3df_fitting_twoLs_two_params_K3iso0_fixed`
- `test1_QC_spectroscopy_twoLs_one_param`
- `test1_QC_spectroscopy_twoLs_two_params`
- `test1_two_params`
- `test2`
- `test3`
- `test5`
- `two_params_fitting`
- `two_params_fitting_K3iso0_fixed`

Inputs/outputs: inspect the file CLI and adjacent reports; this index is intentionally mechanical.

## `source/QC_functions.h`

Purpose: C++/header implementation.

Important functions/classes:
- `F2_mat_builder`
- `F3_ID`
- `F3_ID_mat`
- `F3_ND_2plus1_mat`
- `LinearSolver_3`
- `LinearSolver_4`
- `det_F3_ND_2plus1_mat`
- `function_F3_ND_2plus1_mat`
- `function_F3inv_ND_2plus1_mat`
- `function_for_pole_F3inv_ND_2plus1_mat`
- `if`
- `testF3_additionalpoles_1`
- `testF3_additionalpoles_2`
- `testF3_additionalpoles_3`
- `test_F3_ID_zeroK2`
- `test_F3_ND_2plus1_mat`
- `test_F3iso_ND_2plus1_mat`
- `test_F3iso_ND_2plus1_mat_with_normalization`
- `test_F3iso_ND_2plus1_mat_with_normalization_twobody_var_strength_alpha`

Inputs/outputs: inspect the file CLI and adjacent reports; this index is intentionally mechanical.

## `source/QC_functions_v2.h`

Purpose: C++/header implementation.

Important functions/classes:
- `LinearSolver_3`
- `LinearSolver_4`
- `check_invertibility`
- `if`
- `normalize_det_vector_by_max`
- `solve_Hmat_X_eq_F2mat`
- `test_F3iso_ND_2plus1_mat_with_normalization_single_En`

Inputs/outputs: inspect the file CLI and adjacent reports; this index is intentionally mechanical.

## `source/calc_v29g_F3iso_inputfile.cpp`

Purpose: C++/header implementation.

Important functions/classes:
- `catch`
- `evaluate_F3iso_one_energy`
- `get_double`
- `get_int`
- `get_int_list`
- `get_string`
- `if`
- `main`
- `parse_label`
- `read_kv_file`
- `read_options`
- `split_csv`
- `strip_inline_comment`
- `trim`
- `usage`
- `while`

Inputs/outputs: inspect the file CLI and adjacent reports; this index is intentionally mechanical.

## `source/compare_K3_cpp_vs_python.py`

Purpose: Python workflow/utility.

Important functions/classes:
- `add_qc3_paths`
- `expected_element_from_python`
- `local_index_from_ell_m`
- `main`
- `print`
- `read_cpp_configs`
- `read_cpp_matrix`

Inputs/outputs: inspect the file CLI and adjacent reports; this index is intentionally mechanical.

## `source/compare_gpu_cache_vs_cpu_openmp_det_oldscale.cpp`

Purpose: C++/header implementation.

Important functions/classes:
- `block_prefix`
- `build_param_block`
- `det_complex`
- `die_usage`
- `for`
- `gpu_cache_filename`
- `if`
- `load_gpu_row`
- `ltag`
- `main`
- `make_physics`
- `make_settings`
- `parse_args`
- `print_block_summary`
- `resolve_gpu_cache`
- `run_plot_script`
- `scan_gpu_cache`
- `split_csv`
- `summarize`
- `trim`
- `while`
- `write_csv_header`

Inputs/outputs: inspect the file CLI and adjacent reports; this index is intentionally mechanical.

## `source/debug_L20_100_A2_oldscale_gpu_vs_cpu.cpp`

Purpose: C++/header implementation.

Important functions/classes:
- `block_key`
- `build_curve_row_from_entry`
- `build_matched_row`
- `csv_value`
- `det_complex`
- `die_usage`
- `dimension_jump_indices`
- `failure_message`
- `for`
- `format_value`
- `frob_norm`
- `gpu_cache_filename`
- `if`
- `join_strings`
- `load_gpu_row`
- `ltag`
- `main`
- `make_csv_block`
- `make_inline_table_header`
- `make_physics`
- `make_settings`
- `orthonormality_error`
- `pad3`
- `parse_args`
- `parse_bool`
- `projector_error`
- `quote_csv`
- `read_csv_text`
- `read_file`
- `resolve_gpu_cache`
- `scan_gpu_cache`
- `sign_real`
- `sign_str`
- `split_csv`
- `stable_slogdet`
- `trim`
- `while`
- `window_tag`
- `write_block_summary`
- `write_curve_header`
- `write_curve_row`
- `write_det_summary`
- `write_jump_table`
- `write_matched_header`
- `write_matched_row`
- `write_matched_summary`
- `write_report`
- `write_text_file`
- `write_worst_table`

Inputs/outputs: inspect the file CLI and adjacent reports; this index is intentionally mechanical.

## `source/debug_covariance_input_v24.cpp`

Purpose: C++/header implementation.

Important functions/classes:
- `E_to_Ecm_local`
- `catch`
- `for`
- `if`
- `jackknife_average_local`
- `jackknife_resampling_local`
- `main`
- `mom_label_to_nP_python_covariance_convention`
- `read_second_column_skip_header`
- `while`

Inputs/outputs: inspect the file CLI and adjacent reports; this index is intentionally mechanical.

## `source/debug_fresh_gpu_cpu_qc_elementwise_L20_100_A2.cpp`

Purpose: C++/header implementation.

Important functions/classes:
- `augment_gpu_row`
- `basis_signature`
- `build_k3_full`
- `cache_fname`
- `compare_row`
- `compute_projected_spectral`
- `determinant_complex`
- `dump_elementwise_matrix`
- `for`
- `format_complex`
- `frob_norm`
- `if`
- `main`
- `matrix_diff_stats`
- `max_abs_imag`
- `old_cache_reason_from_rows`
- `parse_bool_fresh`
- `parse_fresh_args`
- `pick_top_rows`
- `projector_error`
- `projector_idempotence_error`
- `rel_hermiticity`
- `resolve_cache_file`
- `row_stats_line`
- `run_fresh_debug`
- `select_debug_rows`
- `select_old_cache_rows`
- `selected_phase1_rows`
- `sort_unique_rows`
- `spectral_stats`
- `stable_slogdet`
- `vsel_alignment_stats`
- `write_element_dump_header`
- `write_manifest`
- `write_report_file`
- `write_row_summary_header`
- `write_row_summary_row`
- `write_text_file`
- `xi_tag`

Inputs/outputs: inspect the file CLI and adjacent reports; this index is intentionally mechanical.

## `source/debug_matrix_read_conventions_L20_100_A2.cpp`

Purpose: C++/header implementation.

Important functions/classes:
- `build_cpu_entry`
- `build_k3_full`
- `classify_full_row`
- `collect_manifest_files`
- `compare_matrices`
- `dense_from_pairs`
- `det_complex_scaled`
- `ensure_parent`
- `evaluate_exact_variant`
- `exact_targets_small`
- `for`
- `format_complex`
- `fresh_cache_name`
- `hermitian_from_triangle`
- `hermiticity_rel`
- `if`
- `join_strings`
- `lookup_variant`
- `main`
- `matrix_classification`
- `max_diag_imag`
- `nearest_indices_for_targets`
- `old_targets_selected`
- `parse_bool_mr`
- `parse_options`
- `point_classification`
- `projection_stats`
- `read_raw_gpu_row`
- `reconstruct_matrix`
- `reconstruct_variant_pair`
- `resolve_cache_with_coarse`
- `run_full_phase`
- `run_small_exact_phase`
- `scan_cache_rows_sorted`
- `select_exact_rows`
- `singular_stats`
- `swap_real_imag`
- `tag_ecm`
- `tag_xi`
- `top_indices_by_metric`
- `unique_rows_by_ecm`
- `validate_csv_file`
- `variant_score`
- `variant_specs`
- `while`
- `worst_hermiticity_pairs`
- `write_final_report`
- `write_manifest`
- `write_matrix_csv`
- `write_matrix_diff_csv`
- `write_matrix_dump_and_diff`
- `write_selected_element_dumps`
- `write_text_file`

Inputs/outputs: inspect the file CLI and adjacent reports; this index is intentionally mechanical.

## `source/debug_v29f_eigenbranch_tracker_sigma_sort.cpp`

Purpose: C++/header implementation.

Important functions/classes:
- `apply_final_filters`
- `catch`
- `cluster_zero_events`
- `evaluate_qc_eigen`
- `for`
- `greedy_assignment`
- `if`
- `k3_parameters_are_zero`
- `linroot`
- `main`
- `make_params`
- `max_overlap_assignment`
- `parse_double`
- `parse_int`
- `parse_label`
- `usage`
- `while`
- `write_cluster_and_final_files`
- `write_headers`

Inputs/outputs: inspect the file CLI and adjacent reports; this index is intentionally mechanical.

## `source/debug_v29h_F3inv_zero_compare_inputfile.cpp`

Purpose: C++/header implementation.

Important functions/classes:
- `append_zero_records`
- `bisection_eigenbranch_minabs_zero`
- `bisection_eigenbranch_reference_zero`
- `bisection_f3iso_zero`
- `bisection_method_zero`
- `catch`
- `classify_side_abs`
- `classify_zeros`
- `coarse_eigensystems_for_method`
- `collect_det_eigenbranch_candidates_from_adjacent_pairs`
- `collect_zero_candidates_for_method`
- `dedup_eigenbranch_candidate_records`
- `determinant_scaled_partial_piv_lu`
- `dump_matrix`
- `eigensystem_for_method`
- `eval_minabs_eigenvalue_real`
- `eval_reference_tracked_eigenvalue_real`
- `evaluate_data`
- `evaluate_full`
- `evaluate_full_grid_parallel`
- `evaluate_grid_parallel`
- `evaluate_refine_grid_data_parallel`
- `evaluate_refine_grid_full_parallel`
- `extract_eval_data_grid`
- `find_candidate_records_for_det_method`
- `find_candidate_records_for_method`
- `find_zeros_for_method`
- `finite_comp_v29k`
- `finite_complex`
- `finite_double_v29k`
- `finite_matrix_v29k`
- `for`
- `get_double`
- `get_int`
- `get_int_list`
- `get_string`
- `if`
- `local_zero_like_crossing`
- `main`
- `mark_failure_v29k`
- `method_label`
- `method_name`
- `method_value`
- `nearest_index_within`
- `parse_label`
- `read_kv_file`
- `read_options`
- `refine_det_candidate_eigenbranch`
- `refine_f3iso_candidate_sideclass`
- `refine_zero_candidates_parallel`
- `same_coarse_interval`
- `sanitize_error`
- `sign_flip`
- `split_csv`
- `strip_inline_comment`
- `tracked_eigenvalue_branches`
- `trim`
- `unique_sorted_zeros`
- `validate_eval_v29k`
- `while`
- `write_grid_file`
- `write_match_summary`
- `write_refine_log_branch`
- `write_refine_log_data`
- `write_zero_classification`
- `zero_record_residual_for_merge`

Inputs/outputs: inspect the file CLI and adjacent reports; this index is intentionally mechanical.

## `source/debug_v29q_detProjF3inv_spwave_111_A2.cpp`

Purpose: C++/header implementation.

Important functions/classes:
- `append_zero_records`
- `bisection_eigenbranch_minabs_zero`
- `bisection_eigenbranch_reference_zero`
- `bisection_f3iso_zero`
- `bisection_method_zero`
- `catch`
- `classify_side_abs`
- `classify_zeros`
- `coarse_eigensystems_for_method`
- `collect_det_eigenbranch_candidates_from_adjacent_pairs`
- `collect_zero_candidates_for_method`
- `dedup_eigenbranch_candidate_records`
- `determinant_scaled_partial_piv_lu`
- `dump_matrix`
- `eigensystem_for_method`
- `eval_minabs_eigenvalue_real`
- `eval_reference_tracked_eigenvalue_real`
- `evaluate_data`
- `evaluate_full`
- `evaluate_full_grid_parallel`
- `evaluate_grid_parallel`
- `evaluate_refine_grid_data_parallel`
- `evaluate_refine_grid_full_parallel`
- `extract_eval_data_grid`
- `find_candidate_records_for_det_method`
- `find_candidate_records_for_method`
- `find_zeros_for_method`
- `finite_comp_v29k`
- `finite_complex`
- `finite_double_v29k`
- `finite_matrix_v29k`
- `for`
- `get_double`
- `get_int`
- `get_int_list`
- `get_string`
- `if`
- `local_zero_like_crossing`
- `main`
- `mark_failure_v29k`
- `method_label`
- `method_name`
- `method_value`
- `nearest_index_within`
- `parse_label`
- `read_kv_file`
- `read_options`
- `refine_det_candidate_eigenbranch`
- `refine_f3iso_candidate_sideclass`
- `refine_zero_candidates_parallel`
- `same_coarse_interval`
- `sanitize_error`
- `sign_flip`
- `split_csv`
- `strip_inline_comment`
- `tracked_eigenvalue_branches`
- `trim`
- `unique_sorted_zeros`
- `validate_eval_v29k`
- `while`
- `write_grid_file`
- `write_match_summary`
- `write_refine_log_branch`
- `write_refine_log_data`
- `write_zero_classification`
- `zero_record_residual_for_merge`

Inputs/outputs: inspect the file CLI and adjacent reports; this index is intentionally mechanical.

## `source/debug_v30d_detProjF3inv_spwave_111_A2_whole_region.cpp`

Purpose: C++/header implementation.

Important functions/classes:
- `V30dScopedTimer`
- `append_zero_records`
- `bisection_eigenbranch_minabs_zero`
- `bisection_eigenbranch_reference_zero`
- `bisection_f3iso_zero`
- `bisection_method_zero`
- `catch`
- `classify_side_abs`
- `classify_zeros`
- `coarse_eigensystems_for_method`
- `collect_det_eigenbranch_candidates_from_adjacent_pairs`
- `collect_zero_candidates_for_method`
- `dedup_eigenbranch_candidate_records`
- `determinant_scaled_partial_piv_lu`
- `dump_matrix`
- `eigensystem_for_method`
- `eval_minabs_eigenvalue_real`
- `eval_reference_tracked_eigenvalue_real`
- `evaluate_data`
- `evaluate_full`
- `evaluate_full_grid_parallel`
- `evaluate_grid_parallel`
- `evaluate_refine_grid_data_parallel`
- `evaluate_refine_grid_full_parallel`
- `extract_eval_data_grid`
- `find_candidate_records_for_det_method`
- `find_candidate_records_for_method`
- `find_zeros_for_method`
- `finite_comp_v29k`
- `finite_complex`
- `finite_double_v29k`
- `finite_matrix_v29k`
- `for`
- `get_double`
- `get_int`
- `get_int_list`
- `get_string`
- `gram_orth_res`
- `hermitian_rel_res`
- `idempotent_rel_res`
- `if`
- `local_zero_like_crossing`
- `main`
- `mark_failure_v29k`
- `method_label`
- `method_name`
- `method_value`
- `modified_gram_schmidt_complex`
- `nearest_index_within`
- `parse_label`
- `read_kv_file`
- `read_options`
- `refine_det_candidate_eigenbranch`
- `refine_f3iso_candidate_sideclass`
- `refine_zero_candidates_parallel`
- `same_coarse_interval`
- `sanitize_error`
- `sign_flip`
- `split_csv`
- `stage_log`
- `strip_inline_comment`
- `tracked_eigenvalue_branches`
- `trim`
- `unique_sorted_zeros`
- `validate_eval_v29k`
- `wall_seconds_now`
- `while`
- `write_grid_file`
- `write_match_summary`
- `write_refine_log_branch`
- `write_refine_log_data`
- `write_zero_classification`
- `zero_record_residual_for_merge`

Inputs/outputs: inspect the file CLI and adjacent reports; this index is intentionally mechanical.

## `source/debug_v30r_svd_validated_cartesian_l1_100_A2_zero_scan.cpp`

Purpose: C++/header implementation.

Important functions/classes:
- `V30eScopedTimer`
- `append_zero_records`
- `basis_signature_v30q`
- `bisection_eigenbranch_minabs_zero`
- `bisection_eigenbranch_reference_zero`
- `bisection_f3iso_zero`
- `bisection_method_zero`
- `build_cached_projector_v30q`
- `catch`
- `classify_side_abs`
- `classify_zeros`
- `coarse_eigensystems_for_method`
- `collect_det_eigenbranch_candidates_from_adjacent_pairs`
- `collect_svd_local_min_candidates_v30r`
- `collect_zero_candidates_for_method`
- `dedup_eigenbranch_candidate_records`
- `determinant_scaled_partial_piv_lu`
- `dimension_jump_midpoints_v30r`
- `dimension_jump_mids_v30q`
- `dump_matrix`
- `eigensystem_for_method`
- `eval_minabs_eigenvalue_real`
- `eval_reference_tracked_eigenvalue_real`
- `evaluate_data`
- `evaluate_full`
- `evaluate_full_grid_parallel`
- `evaluate_grid_parallel`
- `evaluate_refine_grid_data_parallel`
- `evaluate_refine_grid_full_parallel`
- `extract_eval_data_grid`
- `find_candidate_records_for_det_method`
- `find_candidate_records_for_method`
- `find_zeros_for_method`
- `finite_comp_v29k`
- `finite_complex`
- `finite_double_v29k`
- `finite_matrix_v29k`
- `for`
- `get_double`
- `get_int`
- `get_int_list`
- `get_projector_cached_v30q`
- `get_string`
- `gram_orth_res`
- `hermitian_rel_res`
- `idempotent_rel_res`
- `if`
- `interval_hits_reference_window_v30q`
- `interval_near_dimension_jump_v30q`
- `local_zero_like_crossing`
- `main`
- `mark_failure_v29k`
- `method_label`
- `method_name`
- `method_value`
- `modified_gram_schmidt_complex`
- `near_any_dimension_jump_v30r`
- `nearest_index_within`
- `parse_label`
- `progress_percent_log`
- `read_kv_file`
- `read_options`
- `reference_swave_100_A2_defaults_v30q`
- `refine_det_candidate_eigenbranch`
- `refine_f3iso_candidate_sideclass`
- `refine_svd_candidates_parallel_v30r`
- `refine_svd_min_window_v30r`
- `refine_zero_candidates_parallel`
- `same_coarse_interval`
- `sanitize_error`
- `serialize_int_config_v30q`
- `serialize_lm_config_v30q`
- `sign_flip`
- `split_csv`
- `stage_log`
- `strip_inline_comment`
- `tracked_eigenvalue_branches`
- `trim`
- `unique_sorted_zeros`
- `validate_eval_v29k`

Inputs/outputs: inspect the file CLI and adjacent reports; this index is intentionally mechanical.

## `source/debug_v30s_eigen_min_bisect_cartesian_l1_100_A2_zero_scan.cpp`

Purpose: C++/header implementation.

Important functions/classes:
- `bisect_eigen_min_candidate_v30s`
- `closest_zero_eigenpair_from_matrix_v30s`
- `collect_eigen_min_signflips_v30s`
- `eigmin_point_from_eval_v30s`
- `eval_closest_eig_tracked_v30s`
- `for`
- `hermitize_v30s`
- `if`
- `main`
- `refine_candidates_parallel_v30s`
- `sign_flip_v30s`
- `write_v30s_outputs`

Inputs/outputs: inspect the file CLI and adjacent reports; this index is intentionally mechanical.

## `source/debug_v30t_eigen_min_bisect_residual_checked_100_A2_zero_scan.cpp`

Purpose: C++/header implementation.

Important functions/classes:
- `bisect_eigen_min_candidate_v30t`
- `closest_zero_eigenpair_from_matrix_v30t`
- `collect_eigen_min_signflips_v30t`
- `eigmin_point_from_eval_v30t`
- `eval_closest_eig_tracked_v30t`
- `eval_sigma_min_projected_v30t`
- `for`
- `hermitize_v30t`
- `if`
- `main`
- `refine_candidates_parallel_v30t`
- `sign_flip_v30t`
- `write_v30t_outputs`

Inputs/outputs: inspect the file CLI and adjacent reports; this index is intentionally mechanical.

## `source/debug_v30z_eigen_min_localmin_pole_reject_A2_zero_scan.cpp`

Purpose: C++/header implementation.

Important functions/classes:
- `closest_zero_eigenpair_from_matrix_v30v`
- `collect_eigen_min_signflips_v30v`
- `eigmin_point_from_eval_v30v`
- `eval_closest_eig_tracked_v30v`
- `eval_sigma_min_projected_v30v`
- `for`
- `hermitize_v30v`
- `if`
- `localmin_eigen_min_candidate_v30v`
- `main`
- `refine_candidates_parallel_v30v`
- `scan_local_min_abs_closest_eig_v30v`
- `sign_flip_v30v`
- `write_v30v_outputs`

Inputs/outputs: inspect the file CLI and adjacent reports; this index is intentionally mechanical.

## `source/debug_v31b_closest_eig_midpoint_A2_zero_pole_scan.cpp`

Purpose: C++/header implementation.

Important functions/classes:
- `collect_v31a_flips`
- `eigen_summary_from_projF3inv_v31a`
- `for`
- `hermitize_v31a`
- `if`
- `main`
- `v31a_point_from_eval`
- `write_v31a_outputs`

Inputs/outputs: inspect the file CLI and adjacent reports; this index is intentionally mechanical.

## `source/diagnose_QC_mixed_sp_checks_v31p.cpp`

Purpose: C++/header implementation.

Important functions/classes:
- `assemble_QC_fullinv_projected`
- `block_norm`
- `closest_eig_v31p`
- `flattened_ell`
- `for`
- `gd`
- `gi`
- `gs`
- `list_int_str_v31p`
- `main`
- `offdiag_ell_norm`
- `parse_int_list_v31p`
- `parse_wave_sets_v31p`
- `read_kv_v31p`
- `relnorm`
- `safe_comp`
- `safe_int`
- `settings_from_config_v31p`
- `sp_norm`
- `split_ws_v31p`
- `trim_v31p`
- `wave_tag_v31p`
- `write_flat_compare`

Inputs/outputs: inspect the file CLI and adjacent reports; this index is intentionally mechanical.

## `source/dig_tools.hpp`

Purpose: C++/header implementation.

Important functions/classes:
- `ScopedTimer`
- `catch`
- `if`
- `ms`
- `print_function_output_to_file`
- `printer`
- `start`

Inputs/outputs: inspect the file CLI and adjacent reports; this index is intentionally mechanical.

## `source/digonto_classifier_v3.hpp`

Purpose: C++/header implementation.

Important functions/classes:
- `classify_one_flip`
- `classify_series`
- `classify_shoulder`
- `finite`
- `if`
- `linear_zero`
- `same_nonzero_sign`
- `sign`
- `while`

Inputs/outputs: inspect the file CLI and adjacent reports; this index is intentionally mechanical.

## `source/eigenvalue_tracker.hpp`

Purpose: C++/header implementation.

Important functions/classes:
- `compute_eigs_at_E`
- `for`
- `if`
- `match_levels`
- `overlap_matrix`
- `print_eigenvalue_tracks`
- `print_spectrum`
- `smallest_eigenvalue`
- `track_eigenvalues`

Inputs/outputs: inspect the file CLI and adjacent reports; this index is intentionally mechanical.

## `source/extract_F3inv_isotropic_from_gpu_cache_v33a.cpp`

Purpose: C++/header implementation.

Important functions/classes:
- `if`
- `isotropic_projected_vector`
- `main`
- `read_scalar`
- `scaled_signed_logdet`
- `while`

Inputs/outputs: inspect the file CLI and adjacent reports; this index is intentionally mechanical.

## `source/file_processor.hpp`

Purpose: C++/header implementation.

Important functions/classes:
- `classify_sign_flips_zero_vs_pole`
- `for`
- `if`
- `median_of_vector`
- `while`

Inputs/outputs: inspect the file CLI and adjacent reports; this index is intentionally mechanical.

## `source/filter_v29f_zeros_by_sigma.cpp`

Purpose: C++/header implementation.

Important functions/classes:
- `catch`
- `if`
- `main`
- `parse_double`
- `parse_int`
- `usage`
- `while`

Inputs/outputs: inspect the file CLI and adjacent reports; this index is intentionally mechanical.

## `source/functions.h`

Purpose: C++/header implementation.

Important functions/classes:
- `ERFI_func`
- `E_to_Ecm`
- `Ecm_to_E`
- `En_min_minus_for_config`
- `En_min_plus_for_config`
- `Jfunc`
- `boost`
- `config_maker`
- `config_maker_1`
- `config_maker_2`
- `config_maker_3`
- `config_maker_4`
- `config_maker_4_momentum_first`
- `config_maker_5`
- `constexpr`
- `cutoff_function_1`
- `dawson_func`
- `if`
- `kallentriangle`
- `kmax_for_P0`
- `mysqrt`
- `omega_func`
- `particle_energy`
- `pmom`
- `print_config_5col`
- `q2psq_star`
- `sigma`
- `sigma_pvec_based`
- `threebody_Gpoles`
- `threebody_non_int_energy_lab`
- `threebody_non_int_spectrum`
- `threebody_non_int_spectrum_with_multiplicity`

Inputs/outputs: inspect the file CLI and adjacent reports; this index is intentionally mechanical.

## `source/functions_gpu_config_maker4.cuh`

Purpose: C++/header implementation.

Important functions/classes:
- `build_lm_list_from_waves`
- `candidate_count_per_energy`
- `catch`
- `compact_flags_to_entries_cpu_order`
- `config_maker_4_gpu_safe`
- `en_min_plus_real`
- `entries_to_cpu_vectors`
- `f4_pi`
- `for`
- `gpu_config_maker_4_many_Ecm_to_cpu_vectors`
- `gpu_config_maker_4_many_entries`
- `gpu_config_maker_4_many_to_cpu_vectors`
- `gpu_config_maker_4_single_to_cpu_vectors`
- `if`
- `mark_config_maker_4_candidates_kernel`
- `norm3`
- `omega_real`
- `selected_by_config_maker_4_condition`

Inputs/outputs: inspect the file CLI and adjacent reports; this index is intentionally mechanical.

## `source/functions_gpu_config_maker4_v2_fixed.cuh`

Purpose: C++/header implementation.

Important functions/classes:
- `build_lm_list_from_waves`
- `candidate_count_per_energy`
- `catch`
- `compact_flags_to_entries_cpu_order`
- `config_maker_4_gpu_safe`
- `en_min_plus_real`
- `entries_to_cpu_vectors`
- `f4_pi`
- `for`
- `gpu_config_maker_4_many_Ecm_to_cpu_vectors`
- `gpu_config_maker_4_many_entries`
- `gpu_config_maker_4_many_to_cpu_vectors`
- `gpu_config_maker_4_single_to_cpu_vectors`
- `if`
- `mark_config_maker_4_candidates_kernel`
- `norm3`
- `omega_real`
- `selected_by_config_maker_4_condition`

Inputs/outputs: inspect the file CLI and adjacent reports; this index is intentionally mechanical.

## `source/fv_projector_cartesian_l1_v30m.hpp`

Purpose: C++/header implementation.

Important functions/classes:
- `D_real_element_v30m`
- `P_irrep_projection_2plus1_best`
- `S_single_flavor`
- `U_2plus1`
- `U_single_flavor`
- `W_single_flavor`
- `all_conventions`
- `best_projector`
- `cart_axis_from_m_v30m`
- `cartesian_signed_perm_element_v30m`
- `compose_AB`
- `diagnose_convention`
- `equivariance_rel_res`
- `find_group_index`
- `for`
- `hermitian_rel_res`
- `idempotent_rel_res`
- `if`
- `inverse_signed_perm`
- `irrep_dimension`
- `max_equivariance_best_over_little_group`
- `max_equivariance_over_little_group`
- `name`
- `parity_factor_for_l`
- `projector_for_convention`
- `signed_perm_equal`

Inputs/outputs: inspect the file CLI and adjacent reports; this index is intentionally mechanical.

## `source/generate_F3iso.cpp`

Purpose: C++/header implementation.

Important functions/classes:
- `generate_F3iso_L20_from_F3`
- `main`
- `real`

Inputs/outputs: inspect the file CLI and adjacent reports; this index is intentionally mechanical.

## `source/generate_K3iso.cpp`

Purpose: C++/header implementation.

Important functions/classes:
- `generate_K3iso_L20_from_F3`
- `generate_K3iso_using_F3_vs_En_KKpi_omp_single_irrep_centralvalue`
- `generate_K3iso_using_F3_vs_En_KKpi_omp_single_irrep_with_bounds`
- `if`
- `main`
- `real`

Inputs/outputs: inspect the file CLI and adjacent reports; this index is intentionally mechanical.

## `source/generate_eigen_based_F3inv.cpp`

Purpose: C++/header implementation.

Important functions/classes:
- `main`
- `real`

Inputs/outputs: inspect the file CLI and adjacent reports; this index is intentionally mechanical.

## `source/generate_likely_zeros_K3df_v26.cpp`

Purpose: C++/header implementation.

Important functions/classes:
- `catch`
- `for`
- `if`
- `main`
- `parse_double`
- `parse_int`
- `print_usage`

Inputs/outputs: inspect the file CLI and adjacent reports; this index is intentionally mechanical.

## `source/generate_pole.cpp`

Purpose: C++/header implementation.

Important functions/classes:
- `F3_inv_mat_poles_secant_method`
- `generate_pole_L20`
- `main`

Inputs/outputs: inspect the file CLI and adjacent reports; this index is intentionally mechanical.

## `source/generate_spline_based_F3inv.cpp`

Purpose: C++/header implementation.

Important functions/classes:
- `if`
- `main`
- `real`
- `while`

Inputs/outputs: inspect the file CLI and adjacent reports; this index is intentionally mechanical.

## `source/gpu_cache_reader_convention.hpp`

Purpose: C++/header implementation.

Important functions/classes:
- `format_from_header`
- `if`
- `name`
- `read_combined_raw_cache_row_header`
- `read_complex`
- `read_component_cache_header`
- `read_component_matrix_record`
- `read_d`
- `read_header_line`
- `read_i32`
- `read_matrix_col_major`
- `read_string`
- `read_u64`
- `require_gpu_cache_format_header`

Inputs/outputs: inspect the file CLI and adjacent reports; this index is intentionally mechanical.

## `source/gpu_check/cusolver_batched_varsize_solve_AXeqI.cu`

Purpose: C++/header implementation.

Important functions/classes:
- `CHECK_CUDA`
- `CHECK_SOLVER`
- `bytes_per_matrix_estimate`
- `cabs2`
- `cadd`
- `cim`
- `cmul`
- `cre`
- `csub`
- `fill_fn`
- `for`
- `if`
- `main`
- `residual_AX_minus_I`
- `set_identity_batched`
- `size_fn`
- `solve_group_same_n_cusolver`
- `while`

Inputs/outputs: inspect the file CLI and adjacent reports; this index is intentionally mechanical.

## `source/gpu_solvers_batched_streams.cpp`

Purpose: C++/header implementation.

Important functions/classes:
- `checkCublasError`
- `checkCudaError`
- `checkCusolverError`
- `cusolverBatchedQR_pipeline`
- `cusolverBatchedQR_withStreams`
- `cusolverComplex_mat`
- `if`

Inputs/outputs: inspect the file CLI and adjacent reports; this index is intentionally mechanical.

## `source/gpu_solvers_batched_streams_v2.cpp`

Purpose: C++/header implementation.

Important functions/classes:
- `checkCublasError`
- `checkCudaError`
- `checkCusolverError`
- `cublasStatusToString`
- `cusolverBatchedQR_pipeline_fixed_m`
- `cusolverBatchedQR_pipeline_variable_inverse`
- `cusolverBatchedQR_pipeline_variable_m`
- `cusolverStatusToString`
- `for`
- `from_cu`
- `if`
- `to_cu`

Inputs/outputs: inspect the file CLI and adjacent reports; this index is intentionally mechanical.

## `source/gpu_varsize_batched_inverse.cu`

Purpose: C++/header implementation.

Important functions/classes:
- `build_and_invert_energy_sweep_varsize_batched_gpu`
- `build_and_invert_energy_sweep_varsize_batched_gpu_v2`
- `build_and_invert_energy_sweep_varsize_batched_gpu_v3`
- `bytesForOneMatWorstCase`
- `clear`
- `cudaSyncCheck`
- `ensureCapacity`
- `for`
- `get`
- `get_avail_ram_bytes_linux`
- `get_free_vram_bytes`
- `if`
- `invertSameSizeGroup_cublasBatched`
- `invertSameSizeGroup_cublasBatched_cached`
- `invert_and_scatter_chunk`
- `main`
- `makeIdentityBatchHost`
- `packEigenToHostCublas`
- `ram_bytes_per_energy`
- `release`
- `size_t`
- `suggest_chunk_size_from_En_final`
- `throwOnCublas`
- `throwOnCuda`
- `unpackHostCublasToEigen`
- `vram_bytes_per_mat`
- `while`

Inputs/outputs: inspect the file CLI and adjacent reports; this index is intentionally mechanical.

## `source/lattice_data_covariance_cpp.hpp`

Purpose: C++/header implementation.

Important functions/classes:
- `E_to_Ecm`
- `Ecmsq_to_Esq`
- `Esq_to_Ecmsq`
- `convert_jackknife_energy_to_ecm`
- `correlation_matrix`
- `covariance_between_states_szscl21_based`
- `covariance_between_states_szscl21_based_result`
- `default_szscl21_mass_path`
- `default_threebody_path`
- `for`
- `if`
- `jackknife_average`
- `jackknife_error`
- `jackknife_resampling`
- `momentum_magnitude_from_nP`
- `nP_from_mom_label`
- `normalize_lattice_energy_type`
- `print_covariance_summary`
- `read_second_column_skip_header`
- `same_nP`
- `while`

Inputs/outputs: inspect the file CLI and adjacent reports; this index is intentionally mechanical.

## `source/non_int_spectrum_maker_3body_trivial_irrep.py`

Purpose: Python workflow/utility.

Important functions/classes:
- `canonical_mom_maker`
- `config_maker`
- `config_maker_positive_only`
- `energy`
- `full_nonint_spectrum_maker_final`
- `irrep_energy_list_maker`
- `irrep_list_maker`

Inputs/outputs: inspect the file CLI and adjacent reports; this index is intentionally mechanical.

## `source/nonint_degeneracy_F2Giso_v31t.cpp`

Purpose: C++/header implementation.

Important functions/classes:
- `canonical_key`
- `char_irrep`
- `det3`
- `det_lu_v31t`
- `eval_grid`
- `for`
- `gd`
- `gi`
- `gs`
- `if`
- `int_list_str_v31t`
- `irrep_dim_v31t_local`
- `iso_vec`
- `little_group`
- `main`
- `mat_vec`
- `mom2`
- `nonint_shells`
- `oneE`
- `parse_int_list_v31t`
- `parse_wave_sets_v31t`
- `read_kv_v31t`
- `safe_comp`
- `settings_from_config_v31t`
- `signed_perm_mats`
- `split_ws_v31t`
- `three_Ecm`
- `trace3`
- `transform_key`
- `trim_v31t`
- `wave_tag_v31t`

Inputs/outputs: inspect the file CLI and adjacent reports; this index is intentionally mechanical.

## `source/nonint_degeneracy_F2Giso_v31u.cpp`

Purpose: C++/header implementation.

Important functions/classes:
- `canonical_key`
- `char_irrep`
- `det3`
- `eval_grid`
- `for`
- `gd`
- `gi`
- `gs`
- `if`
- `int_list_str_v31u`
- `irrep_dim_v31u_local`
- `iso_vec`
- `little_group`
- `main`
- `mat_vec`
- `mom2`
- `nonint_shells`
- `oneE`
- `parse_int_list_v31u`
- `parse_wave_sets_v31u`
- `read_kv_v31u`
- `safe_comp`
- `settings_from_config_v31u`
- `signed_perm_mats`
- `slogdet_lu_v31u`
- `split_ws_v31u`
- `three_Ecm`
- `transform_key`
- `trim_v31u`
- `wave_tag_v31u`

Inputs/outputs: inspect the file CLI and adjacent reports; this index is intentionally mechanical.

## `source/pole_maker.py`

Purpose: Python workflow/utility.

Important functions/classes:
- `E_to_Ecm`
- `Ecmsq_to_Esq`
- `Esq_to_Ecmsq`
- `P100`
- `P110`
- `P111`
- `P200`
- `Run_polefinder`
- `jackknifeavg_centralvalue_lattice_data`
- `jackknifeavg_lattice_data`
- `pole_P000`
- `pole_finding_by_reading_data_file_F3_P000`
- `pole_finding_by_reading_data_file_F3_P100`
- `pole_finding_by_reading_data_file_F3_P110`
- `pole_finding_by_reading_data_file_F3_P111`
- `pole_finding_by_reading_data_file_F3_P200`
- `pole_finding_by_reading_data_file_F3_all_boost`
- `pole_finding_by_reading_data_file_F3inv_P000`
- `pole_finding_by_reading_data_file_F3inv_P100`
- `pole_finding_by_reading_data_file_F3inv_P110`
- `pole_finding_by_reading_data_file_F3inv_P111`
- `pole_finding_by_reading_data_file_F3inv_P200`
- `pole_finding_by_reading_data_file_F3inv_all_boost`
- `pole_region_finding_by_reading_data_file_F3inv_all_boost`

Inputs/outputs: inspect the file CLI and adjacent reports; this index is intentionally mechanical.

## `source/pole_searching.h`

Purpose: C++/header implementation.

Important functions/classes:
- `F3_inv_mat_poles_secant_method`
- `det_F3mat_poles_secant_method`
- `function_for_pole_F3inv_ND_2plus1_mat`
- `if`
- `test_F3inv_pole_searching_vs_L`

Inputs/outputs: inspect the file CLI and adjacent reports; this index is intentionally mechanical.

## `source/printer_F3_acc.cpp`

Purpose: C++/header implementation.

Important functions/classes:
- `if`
- `main`
- `test_detF3inv_vs_En_KKpi_acc`

Inputs/outputs: inspect the file CLI and adjacent reports; this index is intentionally mechanical.

## `source/printer_F3_omp.cpp`

Purpose: C++/header implementation.

Important functions/classes:
- `F3_fixing_L24`
- `if`
- `main`
- `test_3body_non_int`
- `test_3body_non_int_with_multiplicity`
- `test_F2_for_missing_poles`
- `test_F3_pole_datagenerator_for_residue_vs_En_KKpi_omp`
- `test_F3_vs_En_KKpi_6_diff_ma_omp`
- `test_F3_vs_En_KKpi_omp`
- `test_F3_vs_En_KKpi_variable_2body_strength_omp`
- `test_F3_vs_En_L24_KKpi_omp`
- `test_F3_vs_En_L24_KKpi_omp_single_irrep_with_bounds`
- `test_F3inv_with_splines`
- `test_F3tilde_vs_En_KKpi_omp`
- `while`

Inputs/outputs: inspect the file CLI and adjacent reports; this index is intentionally mechanical.

## `source/printer_function.cpp`

Purpose: C++/header implementation.

Important functions/classes:
- `K2printer`
- `activated_shell`
- `if`
- `main`
- `p_in_lattice_units`
- `poles_of_G_in_Ecm`
- `pvec_by_hand`
- `test_3body_non_int`
- `test_F2_i1_mombased_vs_En`
- `test_F2_i_mat`
- `test_F2_sum_func`
- `test_F2_vs_sigp`
- `test_F3_ID_printer`
- `test_F3_mat`
- `test_F3_mat_vs_En`
- `test_F3_nd_2plus1`
- `test_G_ij_mat`
- `test_Gmat_vs_sigp`
- `test_K2_i_mat`
- `test_additionalpoles_in_F3_vs_En_KKpi`
- `test_config_maker`
- `test_detF2inv_vs_En_KKpi`
- `test_detF3_vs_En`
- `test_detF3inv_vs_En`
- `test_detF3inv_vs_En_KKpi`
- `test_detF3inv_vs_En_KKpi_test_nonintpoles`
- `test_functions_with_FRL_codebase_2plus1`
- `test_functions_with_FRL_codebase_ID`
- `test_individual_functions`
- `test_individual_functions_KKpi`
- `test_mass_dependences_F3_2plus1_vs_En`
- `test_uneven_matrix`

Inputs/outputs: inspect the file CLI and adjacent reports; this index is intentionally mechanical.

## `source/printer_function_temp.cpp`

Purpose: C++/header implementation.

Important functions/classes:
- `I00_sum_F_test`
- `K2printer`
- `main`
- `test_F2_i1_mombased`
- `test_F2_i1_mombased_vs_En`
- `test_QC3_vs_En`
- `test_cutoff_function_1`

Inputs/outputs: inspect the file CLI and adjacent reports; this index is intentionally mechanical.

## `source/projected_F3inv_zero_finder_v31z.cpp`

Purpose: C++/header implementation.

Important functions/classes:
- `add_iterative_signflip_mesh_v32`
- `add_refined_points_for_bracket_v31zm`
- `add_refined_signflip_mesh_v31zm`
- `apply_scatter_params_from_config_v31zl`
- `basis_label`
- `canonical_key`
- `char_irrep`
- `classify_zero`
- `comp_get_int`
- `count_near_zero`
- `det3`
- `detproj_re_for_refine_v32`
- `detproj_sign_for_refine_v32`
- `eval_detail`
- `eval_grid`
- `for`
- `gd`
- `gi`
- `gs`
- `if`
- `int_list_str`
- `interp_detproj_candidate_v32`
- `interp_zero`
- `inverse_scaled_slogdet_from_matrix`
- `irrep_dim_local`
- `ivec_get`
- `little_group`
- `main`
- `mat_vec`
- `max_abs_eig`
- `min_abs_eig`
- `min_abs_eig_index`
- `mom2`
- `nearest_eig_index`
- `nonint_shells`
- `oneE`
- `parse_i3`
- `parse_int_list`
- `print_active_scatter_params_v31zl`
- `read_kv`
- `refine_det_root`
- `rel_matrix_res`
- `settings_from_config`
- `sign_flip_v31zm`
- `sign_nonzero`
- `signed_perm_mats`
- `slogdet_lu`
- `split_ws`
- `three_Ecm`
- `transform_key`
- `trim_v31z`
- `valid_detproj_bracket_v32`
- `wave_tag`
- `while`
- `write_projection_audit`

Inputs/outputs: inspect the file CLI and adjacent reports; this index is intentionally mechanical.

## `source/projections.hpp`

Purpose: C++/header implementation.

Important functions/classes:
- `Dmat11`
- `Dmat11_real`
- `Oh_list`
- `P_irrep_subspace_flavor`
- `block_diag`
- `block_diag_sq`
- `build_full_projector`
- `chi`
- `chop`
- `chop_real`
- `conj_class`
- `cubic_transf`
- `for`
- `get_lm_size`
- `if`
- `irrep_dim`
- `irrep_list`
- `irrep_proj_2plus1`
- `is_proper_rotation`
- `little_group`
- `rotations_list`
- `subspace_from_projector`

Inputs/outputs: inspect the file CLI and adjacent reports; this index is intentionally mechanical.

## `source/projections_from_config.hpp`

Purpose: C++/header implementation.

Important functions/classes:
- `build_reorder_perm`
- `for`
- `if`
- `irrep_proj_2plus1_from_config`
- `parse_config`
- `permute_matrix`
- `print_irrep_eigenvalues`
- `unique_lm_list`
- `unique_nnk_list`

Inputs/outputs: inspect the file CLI and adjacent reports; this index is intentionally mechanical.

## `source/projections_gpu_safe.cuh`

Purpose: C++/header implementation.

Important functions/classes:
- `D_lookup_dev`
- `Oh_list_host`
- `P_irrep_projection_2plus1_and_Vsel_gpu_safe`
- `P_irrep_projection_2plus1_gpu_real`
- `all_perms_host`
- `build_PI_kernel`
- `build_projector_from_eigenvectors_near_one_gpu_real`
- `chi_host`
- `conj_class_host`
- `cubic_transform_dev`
- `flatten_n_config`
- `if`
- `irrep_dim_host`
- `is_rotation_host`
- `little_group_host`
- `orbit_key_dev`
- `release`
- `rotations_list_host`
- `sign_dev`
- `sign_int`
- `sort2_int`
- `sort3_abs_order`
- `sort3_int`
- `upload_config`
- `upload_projection_tables`
- `view`

Inputs/outputs: inspect the file CLI and adjacent reports; this index is intentionally mechanical.

## `source/projections_gpu_safe_v2_fixed_move.cuh`

Purpose: C++/header implementation.

Important functions/classes:
- `D_lookup_dev`
- `Oh_list_host`
- `P_irrep_projection_2plus1_and_Vsel_gpu_safe`
- `P_irrep_projection_2plus1_gpu_real`
- `all_perms_host`
- `build_PI_kernel`
- `build_projector_from_eigenvectors_near_one_gpu_real`
- `chi_host`
- `conj_class_host`
- `cubic_transform_dev`
- `flatten_n_config`
- `if`
- `irrep_dim_host`
- `is_rotation_host`
- `little_group_host`
- `orbit_key_dev`
- `release`
- `rotations_list_host`
- `sign_dev`
- `sign_int`
- `sort2_int`
- `sort3_abs_order`
- `sort3_int`
- `upload_config`
- `upload_projection_tables`
- `view`

Inputs/outputs: inspect the file CLI and adjacent reports; this index is intentionally mechanical.

## `source/projections_v1.hpp`

Purpose: C++/header implementation.

Important functions/classes:
- `Dmat`
- `Dmat11`
- `Dmat22`
- `Oh_list`
- `P_irrep_projection_2plus1`
- `P_irrep_projection_single_flavor`
- `RealScalar`
- `SignType`
- `all_perms_of_3`
- `analyze_projector_eigensystem`
- `blockDiag`
- `build_projector_from_eigenvectors_near_one`
- `chi`
- `chop`
- `chop_scalar`
- `complex_to_string`
- `conj_class`
- `constexpr`
- `cubic_transf`
- `for`
- `get_orbit`
- `if`
- `irrep_dim`
- `irrep_list`
- `is_in_rotations_list`
- `is_rotation`
- `little_group`
- `main`
- `make_permutation_matrix_18`
- `orbit_maker`
- `print_bad_entries`
- `print_nonzero_elements_and_counts`
- `projector_subspace_basis`
- `rotations_list`
- `slogdet`
- `test_P_I_v1`
- `test_P_irrep_maker_v1`
- `test_projector_idempotency`
- `unique_sort_orbits`
- `wigner_d_tests`

Inputs/outputs: inspect the file CLI and adjacent reports; this index is intentionally mechanical.

## `source/qc_fitter_norm_refine_v2.cpp`

Purpose: C++/header implementation.

Important functions/classes:
- `build_entry`
- `cand_zero_proxy`
- `classifier_dispatch_info`
- `classify`
- `cnorm`
- `cov_to_corr`
- `det_info`
- `eval_entry_QC`
- `find_QC_zeros_refined`
- `find_QC_zeros_v3_from_coarse_grid`
- `flip`
- `flips`
- `for`
- `gd`
- `gs`
- `if`
- `is_v3_like_mode`
- `is_v4_like_mode`
- `linzero`
- `main`
- `merge_QC_zeros_v4`
- `mesh`
- `nearest`
- `print_cache_summary`
- `print_split_cache_summary`
- `progress_bar_line`
- `read_kv`
- `settings_from_config`
- `sgn`
- `split_fixed_coarse_and_refined`
- `v3_classify_shoulder`
- `v3_finite_eval`
- `v3_same_sign`
- `while`
- `write_matrix`
- `write_outputs`

Inputs/outputs: inspect the file CLI and adjacent reports; this index is intentionally mechanical.

## `source/qc_fitter_norm_refine_v2_multiL.cpp`

Purpose: C++/header implementation.

Important functions/classes:
- `MultiLFCN`
- `assign_block_model`
- `benchmark_fcn`
- `block_cnorm`
- `build_assignment_rows`
- `build_runtime_block`
- `cache_load_sec`
- `catch`
- `coarse_path_for`
- `contains_label`
- `convert_samples_to_ecm`
- `csv_at2`
- `csv_double2`
- `csv_index2`
- `csv_int2`
- `csv_true2`
- `default_gpu_cache_path`
- `dump_fitter_targets_v33k`
- `ensure_gpu_coarse_cache`
- `evaluate_accepted_window`
- `file_exists`
- `for`
- `if`
- `internal_alias`
- `json_escape`
- `key_string`
- `load_accepted_windows`
- `load_basis_metadata_grid`
- `load_gpu_coarse_cache_one`
- `load_multil_jack_targets`
- `load_projected_basis_hotwindows`
- `load_refined_cache_or_empty`
- `load_v33g_runtime_cache_one`
- `load_window_brackets`
- `load_window_grid`
- `ltag`
- `main`
- `mean_vec`
- `model_for`
- `multiconfig_from_config`
- `nearest_candidate_id`
- `nearest_window_row`
- `parse_cache_block`
- `parse_doubles`
- `parse_jack_filename`
- `parse_words`
- `precompute_sec`
- `print_parameter_mask`
- `read_basis_matrix`
- `read_exact`
- `read_jack_values`
- `read_metadata_row_indices`
- `read_raw_lines`
- `read_scalar`
- `refined_path_for`
- `run_l20_hot_window_audit`
- `run_shell`
- `runtime_path_for`
- `settings_for_block`
- `sha256_file2`
- `split_csv2`
- `strip_comment2`
- `swap_real_imag_for_trusted_v33h`
- `trim2`
- `while`
- `write_assignment_files`
- `write_chi2_audit`
- `write_matrix2`
- `write_outputs_multi`
- `write_parameter_sensitivity`
- `write_refined_meta`
- `xi_tag`

Inputs/outputs: inspect the file CLI and adjacent reports; this index is intentionally mechanical.

## `source/real_wigner_d.hpp`

Purpose: C++/header implementation.

Important functions/classes:
- `D_real_element`
- `D_real_matrix`
- `D_real_matrix_proper`
- `approx_equal`
- `approx_equal_matrix`
- `check_m_range`
- `factorial_int`
- `if`
- `midx`
- `print_matrix`
- `real_to_complex_U`
- `rotation_matrix_to_zyz`
- `run_basic_tests`
- `signed_perm_to_matrix`
- `test_D_real_element_vs_matrix`
- `test_D_real_identity`
- `test_D_real_orthogonality`
- `test_D_real_reflection_parity`
- `test_D_real_swap_xy`
- `test_signed_perm_to_matrix`
- `wigner_D_complex_element_from_euler`
- `wigner_D_complex_matrix_from_rotation`
- `wigner_small_d`

Inputs/outputs: inspect the file CLI and adjacent reports; this index is intentionally mechanical.

## `source/scan_K3df_levels_v31l_fixed_params.cpp`

Purpose: C++/header implementation.

Important functions/classes:
- `for`
- `gd`
- `gi`
- `gs`
- `main`
- `read_simple_kv_scan`
- `settings_from_config`
- `split_ws`
- `while`

Inputs/outputs: inspect the file CLI and adjacent reports; this index is intentionally mechanical.

## `source/scan_K3df_levels_v32f_fixed_params.cpp`

Purpose: C++/header implementation.

Important functions/classes:
- `for`
- `gd`
- `gi`
- `gs`
- `main`
- `read_simple_kv_scan`
- `settings_from_config`
- `split_ws`
- `while`

Inputs/outputs: inspect the file CLI and adjacent reports; this index is intentionally mechanical.

## `source/scan_QC_all_eigenvalues_v31l_fixed_params.cpp`

Purpose: C++/header implementation.

Important functions/classes:
- `eval_all_eigenvalues`
- `for`
- `gd`
- `gi`
- `gs`
- `if`
- `main`
- `read_simple_kv_scan_all`
- `settings_from_config`
- `split_ws`
- `while`
- `write_all_eigenvalue_outputs`

Inputs/outputs: inspect the file CLI and adjacent reports; this index is intentionally mechanical.

## `source/scan_QC_fullF3inv_projected_tracked_eigenbranches_v31l.cpp`

Purpose: C++/header implementation.

Important functions/classes:
- `assemble_QC_fullF3inv_projected_v31l`
- `eval_qc_eigendecomp`
- `for`
- `gd`
- `gi`
- `gs`
- `hungarian_minimize_square_v31l`
- `if`
- `main`
- `read_simple_kv_tr`
- `settings_from_config`
- `split_ws`
- `track_by_eigenvector_overlap`
- `while`
- `write_tracked_outputs_for_label`

Inputs/outputs: inspect the file CLI and adjacent reports; this index is intentionally mechanical.

## `source/scan_defaultK3df_QC_v32n_cache_plot.cpp`

Purpose: C++/header implementation.

Important functions/classes:
- `for`
- `gd`
- `gi`
- `gs`
- `main`
- `read_kv_v32n`
- `settings_from_config_v32n`
- `split_ws`
- `while`
- `write_lattice_targets_v32n`

Inputs/outputs: inspect the file CLI and adjacent reports; this index is intentionally mechanical.

## `source/scan_defaultK3df_QC_v32q_4panel.cpp`

Purpose: C++/header implementation.

Important functions/classes:
- `for`
- `gd`
- `gi`
- `gs`
- `if`
- `main`
- `read_kv_v32q`
- `settings_from_config_v32q`
- `slogdet_matrix_v32q`
- `split_ws`
- `while`
- `write_4panel_grid_v32q`
- `write_lattice_file_audit_v32q`
- `write_lattice_targets_v32q`

Inputs/outputs: inspect the file CLI and adjacent reports; this index is intentionally mechanical.

## `source/scan_projF3inv_digonto_classifier_v2.cpp`

Purpose: C++/header implementation.

Important functions/classes:
- `almost_same_E`
- `build_refined_entry_v32s`
- `classify_final_v32s`
- `const_norm_scale_v32s`
- `det_info`
- `eval_from_entry`
- `find_flips_in_pts`
- `find_initial_flips`
- `for`
- `gd`
- `gi`
- `gs`
- `if`
- `linear_zero`
- `main`
- `make_refined_mesh_v32s`
- `nearest_point`
- `read_kv_v32s`
- `refine_irrep_v32s`
- `settings_from_config_v32s`
- `sign_of`
- `signflip_eval`
- `split_ws`
- `while`
- `write_lattice_targets_v32s`
- `write_outputs_for_irrep`

Inputs/outputs: inspect the file CLI and adjacent reports; this index is intentionally mechanical.

## `source/spherical_functions.h`

Purpose: C++/header implementation.

Important functions/classes:
- `ell_m_vector`
- `if`
- `spherical_harmonics`

Inputs/outputs: inspect the file CLI and adjacent reports; this index is intentionally mechanical.

## `source/splines.h`

Purpose: C++/header implementation.

Important functions/classes:
- `Aji`
- `Bji`
- `Cji`
- `Sij`
- `Sij_builder`
- `if`
- `mu_j`
- `pj`
- `print_vec_1D`
- `print_vec_2D`
- `qj`
- `random_test_function_1`
- `test_spline`
- `test_spline_function_1`
- `xj`

Inputs/outputs: inspect the file CLI and adjacent reports; this index is intentionally mechanical.

## `source/splines_test.cpp`

Purpose: C++/header implementation.

Important functions/classes:
- `Aji`
- `Bji`
- `Cji`
- `Sij`
- `Sij_builder`
- `if`
- `main`
- `mu_j`
- `pj`
- `print_vec_1D`
- `print_vec_2D`
- `qj`
- `random_test_function_1`
- `test_spline`
- `test_spline_function_1`
- `xj`

Inputs/outputs: inspect the file CLI and adjacent reports; this index is intentionally mechanical.

## `source/temp_solvers/varsize_grouped_batched_inverse.cu`

Purpose: C++/header implementation.

Important functions/classes:
- `checkCublas`
- `checkCuda`
- `for`
- `if`
- `main`
- `packEigenColMajor`
- `unpackEigenColMajor`

Inputs/outputs: inspect the file CLI and adjacent reports; this index is intentionally mechanical.

## `source/temp_solvers/varsize_grouped_batched_inverse_lib.cu`

Purpose: C++/header implementation.

Important functions/classes:
- `checkCublas`
- `checkCuda`
- `for`
- `if`
- `invert_varsize_mats_batched_gpu`
- `packEigenColMajor`
- `unpackEigenColMajor`

Inputs/outputs: inspect the file CLI and adjacent reports; this index is intentionally mechanical.

## `source/test.cpp`

Purpose: C++/header implementation.

Important functions/classes:
- `main`
- `test_F2_2plus1_mat`
- `test_F2_ang_mat`
- `test_F2_ang_mom_function`
- `test_F3iso_2plus1_mat`
- `test_K2inv_2plus1_mat`
- `test_config_maker_3`
- `test_spherical_functions`

Inputs/outputs: inspect the file CLI and adjacent reports; this index is intentionally mechanical.

## `source/test_F2_gpu_safe_builder.cu`

Purpose: C++/header implementation.

Important functions/classes:
- `catch`
- `compare_matrices`
- `if`
- `main`
- `make_total_P`
- `print_usage`
- `wall_seconds_since`

Inputs/outputs: inspect the file CLI and adjacent reports; this index is intentionally mechanical.

## `source/test_F3_v10_gpu_matrix_adaptive_zeros.cu`

Purpose: C++/header implementation.

Important functions/classes:
- `catch`
- `if`
- `main`

Inputs/outputs: inspect the file CLI and adjacent reports; this index is intentionally mechanical.

## `source/test_F3_v9_cpu_builder_adaptive_zeros.cu`

Purpose: C++/header implementation.

Important functions/classes:
- `catch`
- `if`
- `main`

Inputs/outputs: inspect the file CLI and adjacent reports; this index is intentionally mechanical.

## `source/test_G_gpu_safe_builder.cu`

Purpose: C++/header implementation.

Important functions/classes:
- `if`
- `main`
- `status_from_stats`
- `write_mismatch_entries`

Inputs/outputs: inspect the file CLI and adjacent reports; this index is intentionally mechanical.

## `source/test_K2_functions_gpu_safe.cu`

Purpose: C++/header implementation.

Important functions/classes:
- `if`
- `main`
- `pass_fail`
- `write_mismatches`

Inputs/outputs: inspect the file CLI and adjacent reports; this index is intentionally mechanical.

## `source/test_K2_gpu_safe_builder.cu`

Purpose: C++/header implementation.

Important functions/classes:
- `if`
- `main`
- `status_from_compare`
- `wall_time_sec`
- `write_mismatch_details`

Inputs/outputs: inspect the file CLI and adjacent reports; this index is intentionally mechanical.

## `source/test_K3df_minuit_fit_v24.cpp`

Purpose: C++/header implementation.

Important functions/classes:
- `main`

Inputs/outputs: inspect the file CLI and adjacent reports; this index is intentionally mechanical.

## `source/test_K3df_minuit_fit_v24_DEBUG_SMALL.cpp`

Purpose: C++/header implementation.

Important functions/classes:
- `catch`
- `for`
- `main`

Inputs/outputs: inspect the file CLI and adjacent reports; this index is intentionally mechanical.

## `source/test_K3df_minuit_fit_v25_cached.cpp`

Purpose: C++/header implementation.

Important functions/classes:
- `catch`
- `main`

Inputs/outputs: inspect the file CLI and adjacent reports; this index is intentionally mechanical.

## `source/test_K3df_minuit_fit_v26_reuse_cache.cpp`

Purpose: C++/header implementation.

Important functions/classes:
- `catch`
- `main`

Inputs/outputs: inspect the file CLI and adjacent reports; this index is intentionally mechanical.

## `source/test_K3df_minuit_fit_v31l_lattice_covariance.cpp`

Purpose: C++/header implementation.

Important functions/classes:
- `gd`
- `gi`
- `gs`
- `main`
- `read_simple_kv_main`
- `split_ws`
- `while`

Inputs/outputs: inspect the file CLI and adjacent reports; this index is intentionally mechanical.

## `source/test_K3df_minuit_fit_v32f_fullF3inv_QCfull_cached_classifier.cpp`

Purpose: C++/header implementation.

Important functions/classes:
- `gd`
- `gi`
- `gs`
- `main`
- `read_simple_kv_main`
- `split_ws`
- `while`

Inputs/outputs: inspect the file CLI and adjacent reports; this index is intentionally mechanical.

## `source/test_K3mat_2plus1_cpp_dump.cpp`

Purpose: C++/header implementation.

Important functions/classes:
- `comp`
- `main`
- `write_config_file`
- `write_matrix_file`

Inputs/outputs: inspect the file CLI and adjacent reports; this index is intentionally mechanical.

## `source/test_gpu.cpp`

Purpose: C++/header implementation.

Important functions/classes:
- `build_and_invert_energy_sweep_varsize_batched_gpu_v3`
- `comp`
- `evals`
- `for`
- `if`
- `klm_vec`
- `main`
- `matrix_size_generator`
- `nconfig_check`
- `plm_vec`
- `test_F2_2plus1_mat`
- `test_F2_ang_mat`
- `test_F2_ang_mom_function`
- `test_F3iso_2plus1_mat`
- `test_F3iso_gpu`
- `test_F3iso_gpu_2`
- `test_F3iso_gpu_3`
- `test_K2inv_2plus1_mat`
- `test_config_maker_3`
- `test_projections`
- `test_projections_1`
- `test_projections_cpu`
- `test_projections_gpu`
- `test_projections_gpu_v2`
- `test_spherical_functions`

Inputs/outputs: inspect the file CLI and adjacent reports; this index is intentionally mechanical.

## `source/test_gpu_config_maker4_two_flavor.cu`

Purpose: C++/header implementation.

Important functions/classes:
- `cabs_diff`
- `catch`
- `compare_one_config`
- `if`
- `main`
- `print_usage`
- `test_gpu_config_maker4_two_flavor_over_energy_region`
- `vec3_to_string`

Inputs/outputs: inspect the file CLI and adjacent reports; this index is intentionally mechanical.

## `source/test_gpu_config_maker4_two_flavor_v2_fixed.cu`

Purpose: C++/header implementation.

Important functions/classes:
- `cabs_diff`
- `catch`
- `compare_one_config`
- `if`
- `main`
- `print_usage`
- `test_gpu_config_maker4_two_flavor_over_energy_region`
- `vec3_to_string`

Inputs/outputs: inspect the file CLI and adjacent reports; this index is intentionally mechanical.

## `source/test_gpu_projections.cpp`

Purpose: C++/header implementation.

Important functions/classes:
- `catch`
- `debug_Gij_and_Ylm_for_ijk`
- `det_proj_F2_all`
- `det_proj_F3i_all`
- `det_proj_G_all`
- `eig_proj_F3i_all`
- `find_QC_points_from_normalized_file`
- `for`
- `if`
- `klm_vec`
- `main`
- `make_QC_output_filename`
- `nk_vec`
- `np_vec`
- `plm_vec`
- `print_F2_2plus1_component_for_lm`
- `print_function_output_to_file`
- `print_nan_inf_entries`
- `print_total_dim_vs_energy`
- `test_F3_with_pwave_all_energy_gpu_omp_normalized_v1`
- `test_F3_with_pwave_all_energy_v1`
- `test_F3_with_pwave_single_energy_v1`
- `test_projections_cpu_v3`
- `test_projections_gpu_v2`
- `test_projections_gpu_v3`
- `while`

Inputs/outputs: inspect the file CLI and adjacent reports; this index is intentionally mechanical.

## `source/test_lattice_data_covariance_cpp.cpp`

Purpose: C++/header implementation.

Important functions/classes:
- `main`

Inputs/outputs: inspect the file CLI and adjacent reports; this index is intentionally mechanical.

## `source/test_projection_v9_vs_v18_gpu_vsel_issue.cu`

Purpose: C++/header implementation.

Important functions/classes:
- `catch`
- `compare_matrix`
- `first_eigs_string`
- `if`
- `last_eigs_string`
- `main`
- `print_config_sample`
- `summarize_projector`

Inputs/outputs: inspect the file CLI and adjacent reports; this index is intentionally mechanical.

## `source/test_projections_gpu_safe.cu`

Purpose: C++/header implementation.

Important functions/classes:
- `compare_matrix`
- `if`
- `main`

Inputs/outputs: inspect the file CLI and adjacent reports; this index is intentionally mechanical.

## `source/v33f_F3inv_iso_visual_test.cpp`

Purpose: C++/header implementation.

Important functions/classes:
- `determinant_real`
- `diagnose_row`
- `find_cache`
- `for`
- `gd`
- `gi`
- `gs`
- `if`
- `ints`
- `isotropic_projected_vector`
- `ltag`
- `main`
- `read_kv`
- `read_merged_cache`
- `read_scalar`
- `signed_logabsdet_real`
- `trim`
- `while`
- `words`
- `xi_tag`

Inputs/outputs: inspect the file CLI and adjacent reports; this index is intentionally mechanical.

## `source/v33g_build_runtime_k3basis_cache.cpp`

Purpose: C++/header implementation.

Important functions/classes:
- `coarse_path_for`
- `file_exists`
- `for`
- `if`
- `internal_alias`
- `key_string`
- `load_gpu_coarse_cache_one`
- `ltag`
- `main`
- `make_meta`
- `parse_cache_block`
- `parse_doubles`
- `parse_words`
- `read_build_config`
- `read_exact`
- `read_scalar`
- `runtime_path_for`
- `settings_for_block`
- `strip_comment`
- `trim`
- `while`
- `xi_tag`

Inputs/outputs: inspect the file CLI and adjacent reports; this index is intentionally mechanical.

## `source/v33g_runtime_k3basis_cache.hpp`

Purpose: C++/header implementation.

Important functions/classes:
- `for`
- `gpu_complex_read_convention_name`
- `if`
- `load_runtime_cache`
- `meta_double`
- `meta_int`
- `meta_string`
- `read_compact_matrix`
- `read_gpu_cache_complex`
- `read_gpu_cache_matrix_col_major`
- `read_raw_d`
- `read_raw_i32`
- `read_raw_u64`
- `read_runtime_meta_kv`
- `runtime_meta_path`
- `to_projected_entry_v33g`
- `while`
- `write_compact_matrix`
- `write_json_vec`
- `write_raw_d`
- `write_raw_i32`
- `write_raw_u64`
- `write_runtime_cache`
- `write_runtime_meta_json`

Inputs/outputs: inspect the file CLI and adjacent reports; this index is intentionally mechanical.

## `source/v33g_validate_runtime_vs_v33f.cpp`

Purpose: C++/header implementation.

Important functions/classes:
- `block_tag`
- `build_blocks`
- `for`
- `if`
- `ltag`
- `main`
- `matrix_diff`
- `matrix_scale`
- `model_from_zero_lists`
- `write_text_file`

Inputs/outputs: inspect the file CLI and adjacent reports; this index is intentionally mechanical.

## `source/v33h_patched_gpu_cache_oldscale_det_scan.cpp`

Purpose: C++/header implementation.

Important functions/classes:
- `assemble_scaled_projected_QC_for_row`
- `compute_row`
- `for`
- `format_double`
- `hermiticity_rel`
- `if`
- `load_gpu_row_direct_variant04`
- `main`
- `max_diag_imag`
- `orthonormality_error_local`
- `parse_bool`
- `parse_complex_read_convention`
- `parse_scan_args`
- `process_block`
- `projector_error_p2`
- `quote_csv_local`
- `rank_from_svd`
- `scan_prefix`
- `select_window_rows`
- `stable_slogdet_local`
- `swap_real_imag`
- `tag_window`
- `validate_metadata`
- `write_jump_csv`
- `write_jump_csv_header`
- `write_row_csv`
- `write_row_csv_header`

Inputs/outputs: inspect the file CLI and adjacent reports; this index is intentionally mechanical.

## `source/v33i_single_row_gpu_reader_path_identity_test.cpp`

Purpose: C++/header implementation.

Important functions/classes:
- `bytes_to_string`
- `dense_from_pairs`
- `fmt`
- `for`
- `gpu_cache_filename`
- `hermitian_from_triangle`
- `hermiticity_rel`
- `if`
- `ltag`
- `main`
- `make_physics`
- `make_settings`
- `max_abs_diff`
- `max_diag_imag`
- `orthonormality_error`
- `projector_error`
- `read_raw_row`
- `reconstruct_variant`
- `rel_frob`
- `resolve_gpu_cache`
- `scan_rows`
- `seek_to_matrix_start`
- `select_rows`
- `split_csv`
- `swap_real_imag`
- `top_antiherm_pairs`
- `trim`
- `while`
- `write_csv_header`
- `write_diag_preview`
- `write_flat_preview`
- `write_matrix_preview`

Inputs/outputs: inspect the file CLI and adjacent reports; this index is intentionally mechanical.

## `source/v33p_build_projected_basis_hotwindow_cache.cpp`

Purpose: C++/header implementation.

Important functions/classes:
- `for`
- `if`
- `json_array`
- `json_array_d`
- `json_array_i`
- `k3_matrix`
- `main`
- `now_utc`
- `parse`
- `read_windows`
- `run`
- `split`
- `while`
- `write_matrix`

Inputs/outputs: inspect the file CLI and adjacent reports; this index is intentionally mechanical.

## `source/v33p_cpuassembled_gpu_det_scan.cpp`

Purpose: C++/header implementation.

Important functions/classes:
- `boolean_value`
- `for`
- `if`
- `main`
- `pack_matrices`
- `parse_bridge_args`
- `read_reference`
- `read_zero_centers`
- `run_block`
- `split_line`
- `subset_ids`
- `while`
- `write_report`
- `write_subset_debug`

Inputs/outputs: inspect the file CLI and adjacent reports; this index is intentionally mechanical.

## `source/v33p_gpu_det_lu.cu`

Purpose: C++/header implementation.

Important functions/classes:
- `cublas_error`
- `cublas_ok`
- `cuda_error`
- `cuda_ok`
- `if`
- `set_error`
- `v33p_gpu_available`
- `v33p_gpu_batched_determinants`

Inputs/outputs: inspect the file CLI and adjacent reports; this index is intentionally mechanical.

## `source/v33p_gpu_det_lu.hpp`

Purpose: C++/header implementation.

Important functions/classes:
- No function definitions detected by the lightweight indexer.

Inputs/outputs: inspect the file CLI and adjacent reports; this index is intentionally mechanical.

## `source/v33p_projected_basis_chunk_rawdet_test.cpp`

Purpose: C++/header implementation.

Important functions/classes:
- `basis_parse_args`
- `catch`
- `for`
- `if`
- `k3_matrix`
- `main`
- `read_selected_ids`
- `run`
- `split_csv_line`
- `while`

Inputs/outputs: inspect the file CLI and adjacent reports; this index is intentionally mechanical.

## `scripts/assemble_l20_passonly_reports.py`

Purpose: Python workflow/utility.

Important functions/classes:
- `SystemExit`
- `covariance_status`
- `d`
- `main`
- `read_csv`
- `with`
- `write_csv`

Inputs/outputs: inspect the file CLI and adjacent reports; this index is intentionally mechanical.

## `scripts/assemble_l20_zero_test_reports.py`

Purpose: Python workflow/utility.

Important functions/classes:
- `d`
- `in`
- `main`
- `params`
- `read`
- `write`

Inputs/outputs: inspect the file CLI and adjacent reports; this index is intentionally mechanical.

## `scripts/audit_l20_chi2_reconstruction.py`

Purpose: Python workflow/utility.

Important functions/classes:
- `SystemExit`
- `fit_dir`
- `main`
- `parse_matrix`
- `read_csv`

Inputs/outputs: inspect the file CLI and adjacent reports; this index is intentionally mechanical.

## `scripts/audit_l20_user_validated_prefit.py`

Purpose: Python workflow/utility.

Important functions/classes:
- `lattice_counts`
- `main`
- `with`

Inputs/outputs: inspect the file CLI and adjacent reports; this index is intentionally mechanical.

## `scripts/build_l20_user_validated_windows.py`

Purpose: Python workflow/utility.

Important functions/classes:
- `in`
- `lattice_levels`
- `main`
- `read_label_rows`
- `sha256`
- `with`

Inputs/outputs: inspect the file CLI and adjacent reports; this index is intentionally mechanical.

## `scripts/check_v33f_truezero_counts.py`

Purpose: Python workflow/utility.

Important functions/classes:
- `ltag`
- `main`
- `read_summary`
- `suffix_base`

Inputs/outputs: inspect the file CLI and adjacent reports; this index is intentionally mechanical.

## `scripts/check_v33p_production_scatter_params.py`

Purpose: Python workflow/utility.

Important functions/classes:
- `active_files`
- `all_config_files`
- `equal`
- `main`
- `parse`
- `print`
- `write_audit`

Inputs/outputs: inspect the file CLI and adjacent reports; this index is intentionally mechanical.

## `scripts/diagnose_highEcm_cache_vs_cpu_F3inv.py`

Purpose: Python workflow/utility.

Important functions/classes:
- `as_float`
- `choose_param_values`
- `classify_overall`
- `compile_helper`
- `configure_matplotlib_backend`
- `enumerate`
- `load_meta`
- `main`
- `make_plot`
- `make_raw_det_plots`
- `plot_scaled`
- `read_rows`
- `resolve_cache_file`
- `run`
- `run_helper`
- `stat`
- `write_report`

Inputs/outputs: inspect the file CLI and adjacent reports; this index is intentionally mechanical.

## `scripts/diagnose_highEcm_cache_vs_cpu_F3inv_helper.cpp`

Purpose: C++/header implementation.

Important functions/classes:
- `classify`
- `diagnose_row`
- `for`
- `hermiticity`
- `if`
- `join_complex_values`
- `join_double_values`
- `max_absdiff`
- `notes_for`
- `orthonormality_error`
- `projector`
- `projector_reldiff_norm`
- `relnorm`
- `select_rows`
- `singular_values_desc`
- `sorted_eigenvalues`
- `sv_metrics`
- `write_csv`

Inputs/outputs: inspect the file CLI and adjacent reports; this index is intentionally mechanical.

## `scripts/finalize_l20_k3b_k3e_zero_tests.py`

Purpose: Python workflow/utility.

Important functions/classes:
- `abs`
- `and`
- `d`
- `empty_artifacts`
- `f`
- `main`
- `normalize_traces`
- `not`
- `params`
- `read`
- `repair_blocked_metrics`
- `summarize`
- `write_final`

Inputs/outputs: inspect the file CLI and adjacent reports; this index is intentionally mechanical.

## `scripts/interactive_classifier_compare.py`

Purpose: Python workflow/utility.

Important functions/classes:
- `build_block_cache_path`
- `cluster_candidates`
- `collect_matrix_diagnostics`
- `compare_energy_sets`
- `configure_matplotlib_backend`
- `diag`
- `discover_fit_params`
- `ensure_review_artifacts`
- `evaluate_method`
- `filter_range`
- `load_nonint_levels`
- `load_runtime_cache_with_params`
- `load_runtime_matrices`
- `load_target_dump`
- `main`
- `nearest_energy`
- `parse_float_list`
- `parse_summary_file`
- `plot_block`
- `print`
- `print_block_summary`
- `read_cpx_matrix`
- `read_f64`
- `read_i32`
- `read_u64`
- `runtime_meta`
- `select_candidates`
- `selfcheck`
- `sign`
- `sign_flip_candidates`
- `split_words`
- `summarize_values`
- `to_json`
- `to_row`
- `unique_sorted`
- `write_review_templates`

Inputs/outputs: inspect the file CLI and adjacent reports; this index is intentionally mechanical.

## `scripts/make_spectrum_config_from_fit.py`

Purpose: Python workflow/utility.

Important functions/classes:
- `main`
- `read_summary`
- `replace_key`

Inputs/outputs: inspect the file CLI and adjacent reports; this index is intentionally mechanical.

## `scripts/make_v33e_cache_file_list.py`

Purpose: Python workflow/utility.

Important functions/classes:
- `find_one`
- `main`
- `xi_tag`

Inputs/outputs: inspect the file CLI and adjacent reports; this index is intentionally mechanical.

## `scripts/materialize_v33p_hot_windows.py`

Purpose: Python workflow/utility.

Important functions/classes:
- `main`
- `read_rows`
- `read_targets`

Inputs/outputs: inspect the file CLI and adjacent reports; this index is intentionally mechanical.

## `scripts/merge_projected_basis_rows.py`

Purpose: Python workflow/utility.

Important functions/classes:
- `main`
- `print`
- `records`
- `val`

Inputs/outputs: inspect the file CLI and adjacent reports; this index is intentionally mechanical.

## `scripts/plot_F3inv_iso_nscale_with_v3_zeros.py`

Purpose: Python workflow/utility.

Important functions/classes:
- `load_true_zeros`
- `loadtxt_safe`
- `main`
- `read_rows`

Inputs/outputs: inspect the file CLI and adjacent reports; this index is intentionally mechanical.

## `scripts/plot_F3inv_isotropic_nscale_v33a.py`

Purpose: Python workflow/utility.

Important functions/classes:
- `F3inv_isotropic_rayleigh`
- `main`
- `map`

Inputs/outputs: inspect the file CLI and adjacent reports; this index is intentionally mechanical.

## `scripts/plot_L20_100_A2_oldscale_gpu_vs_cpu_debug.py`

Purpose: Python workflow/utility.

Important functions/classes:
- `add_dim_lines`
- `as_float`
- `as_int`
- `configure_backend`
- `dim_change_rows`
- `infer_stem_parts`
- `main`
- `plot_matched_diff`
- `plot_overlay`
- `read_rows`

Inputs/outputs: inspect the file CLI and adjacent reports; this index is intentionally mechanical.

## `scripts/plot_classifier_algorithm_sweep_v33i.py`

Purpose: Python workflow/utility.

Important functions/classes:
- `analyze_algorithm`
- `assign_branches`
- `branch_windows`
- `build_candidates`
- `classify_block`
- `cluster_candidates`
- `convert_samples_to_ecm`
- `f1`
- `f2`
- `f3`
- `jackknife_average`
- `jackknife_error`
- `jackknife_resampling`
- `load_lattice_levels`
- `load_nonint_candidates`
- `load_nonint_levels`
- `load_rows`
- `load_runtime_cache`
- `local_window_indices`
- `main`
- `median_abs_deviation`
- `nP_from_irrep`
- `nearest_energy`
- `parse_float_list`
- `parse_jack_path`
- `parse_summary`
- `plot_block`
- `read_cpx_matrix`
- `read_f64`
- `read_i32`
- `read_jack_values`
- `read_kv`
- `read_u64`
- `runtime_meta`
- `selfcheck`
- `sign`
- `slope`
- `split_words`
- `to_json`
- `to_row`
- `unique_sorted`
- `write_csv`

Inputs/outputs: inspect the file CLI and adjacent reports; this index is intentionally mechanical.

## `scripts/plot_classifier_algorithm_sweep_v33j.py`

Purpose: Python workflow/utility.

Important functions/classes:
- `adjacent_spacing`
- `analyze_algorithm`
- `assign_branches`
- `branch_assignment_cost`
- `branch_windows`
- `build_branch_anchors`
- `build_candidates`
- `classify_block`
- `classify_block_v6`
- `cluster_candidates`
- `convert_samples_to_ecm`
- `jackknife_average`
- `jackknife_error`
- `jackknife_resampling`
- `load_lattice_levels`
- `load_nonint_candidates`
- `load_nonint_levels`
- `load_rows`
- `load_runtime_cache`
- `local_discontinuity_score`
- `local_noise_score`
- `local_window_indices`
- `main`
- `median_abs_deviation`
- `mk_levels`
- `nP_from_irrep`
- `nearest_energy`
- `parse_float_list`
- `parse_jack_path`
- `parse_summary`
- `plot_block`
- `read_cpx_matrix`
- `read_f64`
- `read_i32`
- `read_jack_values`
- `read_kv`
- `read_u64`
- `run_case`
- `runtime_meta`
- `select_for_upper`
- `selfcheck`
- `sign`
- `slope`
- `split_words`
- `to_json`
- `to_row`
- `unique_sorted`
- `update_window_flags`
- `window_bounds_from_anchors`
- `write_csv`

Inputs/outputs: inspect the file CLI and adjacent reports; this index is intentionally mechanical.

## `scripts/plot_classifier_algorithm_sweep_v33l.py`

Purpose: Python workflow/utility.

Important functions/classes:
- `analyze_algorithm`
- `cluster_candidates`
- `collect_matrix_diagnostics`
- `diag`
- `load_nonint_levels`
- `load_runtime_matrices`
- `load_target_dump`
- `main`
- `median_abs_deviation`
- `nearest_energy`
- `parse_float_list`
- `parse_summary`
- `plot_block`
- `q_matrix_from_components`
- `read_cpx_matrix`
- `read_f64`
- `read_i32`
- `read_kv`
- `read_u64`
- `run_selfcheck`
- `runtime_meta`
- `scan_runtime_cache`
- `select_candidates`
- `sign`
- `sign_flip_candidates`
- `split_words`
- `to_json`
- `to_row`
- `unique_sorted`
- `write_csv`

Inputs/outputs: inspect the file CLI and adjacent reports; this index is intentionally mechanical.

## `scripts/plot_debug_fresh_gpu_cpu_qc_elementwise_L20_100_A2.py`

Purpose: Python workflow/utility.

Important functions/classes:
- `col`
- `configure_backend`
- `load_summary`
- `main`
- `save_fig`

Inputs/outputs: inspect the file CLI and adjacent reports; this index is intentionally mechanical.

## `scripts/plot_fit_levels.py`

Purpose: Python workflow/utility.

Important functions/classes:
- `load`
- `main`

Inputs/outputs: inspect the file CLI and adjacent reports; this index is intentionally mechanical.

## `scripts/plot_matrix_read_convention_debug.py`

Purpose: Python workflow/utility.

Important functions/classes:
- `col`
- `configure_backend`
- `load_rows`
- `main`
- `save_fig`

Inputs/outputs: inspect the file CLI and adjacent reports; this index is intentionally mechanical.

## `scripts/plot_oldscale_det_coarse20000_interactive.py`

Purpose: Python workflow/utility.

Important functions/classes:
- `load_csv`
- `main`
- `parse_args`
- `parse_ylim`
- `select_backend`

Inputs/outputs: inspect the file CLI and adjacent reports; this index is intentionally mechanical.

## `scripts/plot_oldscale_det_compare.py`

Purpose: Python workflow/utility.

Important functions/classes:
- `configure_backend`
- `main`
- `plot_block`
- `read_rows`
- `to_float`

Inputs/outputs: inspect the file CLI and adjacent reports; this index is intentionally mechanical.

## `scripts/plot_v31z_projected_F3inv_zero.py`

Purpose: Python workflow/utility.

Important functions/classes:
- `SystemExit`

Inputs/outputs: inspect the file CLI and adjacent reports; this index is intentionally mechanical.

## `scripts/plot_v31zd_F3inv_5panel.py`

Purpose: Python workflow/utility.

Important functions/classes:
- `col`
- `enumerate`
- `load_grid`
- `load_nonint`
- `load_raw_eigs`
- `main`
- `parse_columns`
- `plot_scaled`
- `scatter_nonint`
- `track_eigenbranches`
- `zip`

Inputs/outputs: inspect the file CLI and adjacent reports; this index is intentionally mechanical.

## `scripts/plot_v31zh_compare_two_scatter_5panel.py`

Purpose: Python workflow/utility.

Important functions/classes:
- `add_case`
- `col`
- `data_paths`
- `load_grid`
- `load_nonint`
- `load_raw_eigs`
- `main`
- `parse_columns`
- `plot_scaled`
- `scatter_nonint`
- `track_eigenbranches`

Inputs/outputs: inspect the file CLI and adjacent reports; this index is intentionally mechanical.

## `scripts/plot_v31zi_F3inv_6panel.py`

Purpose: Python workflow/utility.

Important functions/classes:
- `col`
- `digonto_classifier_v1`
- `draw_digonto_classifier_v1_lines`
- `enumerate`
- `len`
- `load_grid`
- `load_nonint`
- `load_raw_eigs`
- `local_peak_ratio_v31zn`
- `main`
- `open`
- `orientation_ok_v31zn`
- `parse_columns`
- `plot_scaled`
- `print`
- `range`
- `scatter_nonint`
- `track_eigenbranches`
- `write_digonto_classifier_v1_report`

Inputs/outputs: inspect the file CLI and adjacent reports; this index is intentionally mechanical.

## `scripts/plot_v31zi_compare_two_scatter_6panel.py`

Purpose: Python workflow/utility.

Important functions/classes:
- `add_case`
- `col`
- `data_paths`
- `digonto_classifier_v1`
- `draw_digonto_classifier_v1_lines`
- `load_grid`
- `load_nonint`
- `load_raw_eigs`
- `local_peak_ratio_v31zn`
- `main`
- `orientation_ok_v31zn`
- `parse_columns`
- `plot_scaled`
- `range`
- `scatter_nonint`
- `track_eigenbranches`
- `write_digonto_classifier_v1_report`

Inputs/outputs: inspect the file CLI and adjacent reports; this index is intentionally mechanical.

## `scripts/plot_v32_digonto_classifier_v1_6panel_sorted.py`

Purpose: Python workflow/utility.

Important functions/classes:
- `add_case`
- `digonto_classifier_v1`
- `find_eigs`
- `find_grid`
- `finite_masked_y`
- `interp_zero`
- `load_eigen_branches`
- `load_table`
- `main`
- `parse_columns`
- `peak_ratio`
- `plot_scaled`
- `safe_sign`
- `signflip_indices`
- `write_candidate_report`

Inputs/outputs: inspect the file CLI and adjacent reports; this index is intentionally mechanical.

## `scripts/plot_v32m_bestfit_QC_3panel.py`

Purpose: Python workflow/utility.

Important functions/classes:
- `finite_mask`
- `label_from_grid`
- `load_table`
- `main`
- `nonint_2Kpi_levels`
- `omega`
- `parse_cols`
- `parse_total_momentum`
- `read_fit_levels`

Inputs/outputs: inspect the file CLI and adjacent reports; this index is intentionally mechanical.

## `scripts/plot_v32n_defaultK3df_QC_3panel.py`

Purpose: Python workflow/utility.

Important functions/classes:
- `abs`
- `clean_label_from_grid`
- `finite_xy`
- `load_lattice_targets`
- `load_predictions`
- `load_table`
- `main`
- `nonint_2Kpi_levels`
- `omega`
- `parse_cols`
- `parse_total_momentum`

Inputs/outputs: inspect the file CLI and adjacent reports; this index is intentionally mechanical.

## `scripts/plot_v32s_projF3inv_norm_refine_v2_irrep.py`

Purpose: Python workflow/utility.

Important functions/classes:
- `gf`
- `infer_nonint_file`
- `len`
- `load_grid`
- `load_lattice`
- `load_nonint`
- `load_true_zero_candidates`
- `main`
- `parse_args`
- `parse_label`
- `print`
- `read_rows`

Inputs/outputs: inspect the file CLI and adjacent reports; this index is intentionally mechanical.

## `scripts/plot_v32t_projF3inv_norm_refine_v2_irrep_orange.py`

Purpose: Python workflow/utility.

Important functions/classes:
- `gf`
- `infer_nonint_file`
- `len`
- `load_grid`
- `load_lattice`
- `load_nonint`
- `load_true_zero_candidates`
- `main`
- `parse_args`
- `parse_label`
- `print`
- `read_rows`

Inputs/outputs: inspect the file CLI and adjacent reports; this index is intentionally mechanical.

## `scripts/plot_v33g_runtime_qc_det_compare.py`

Purpose: Python workflow/utility.

Important functions/classes:
- `load_runtime_cache`
- `load_runtime_cache_with_params`
- `load_runtime_meta`
- `main`
- `parse_fit_levels`
- `parse_fit_summary`
- `parse_zero_file`
- `read_cpx_matrix`
- `read_f64`
- `read_i32`
- `read_rows`
- `read_u64`
- `same_sign`
- `sign`
- `v3_classify_shoulder`
- `v3_true_zeros`
- `v4_merge_true_zeros`

Inputs/outputs: inspect the file CLI and adjacent reports; this index is intentionally mechanical.

## `scripts/plot_v33p_classifier_det_family_interactive.py`

Purpose: Python workflow/utility.

Important functions/classes:
- `energy`
- `main`
- `rows`

Inputs/outputs: inspect the file CLI and adjacent reports; this index is intentionally mechanical.

## `scripts/prepare_v33p_L20_200_from_existing_raw.py`

Purpose: Python workflow/utility.

Important functions/classes:
- `main`
- `with`

Inputs/outputs: inspect the file CLI and adjacent reports; this index is intentionally mechanical.

## `scripts/prepare_v33p_all_sector_validation.py`

Purpose: Python workflow/utility.

Important functions/classes:
- `SystemExit`
- `candidate_rows`
- `f`
- `main`
- `meta_for`
- `process`
- `sign`

Inputs/outputs: inspect the file CLI and adjacent reports; this index is intentionally mechanical.

## `scripts/run_l20_cutoff_dependence.py`

Purpose: Python workflow/utility.

Important functions/classes:
- `SystemExit`
- `last_float`
- `last_int`
- `main`
- `make_config`
- `parse_summary`
- `run_logged`
- `set_key`

Inputs/outputs: inspect the file CLI and adjacent reports; this index is intentionally mechanical.

## `scripts/run_l20_full_cutoff_fits_v2.py`

Purpose: Python workflow/utility.

Important functions/classes:
- `SystemExit`
- `accepted_source_paths`
- `alias`
- `all`
- `cfg_value`
- `create_per_level`
- `enumerate`
- `filtered_zero_files`
- `in`
- `lattice_records`
- `main`
- `make_config`
- `matrix_dat`
- `nearest_row`
- `parse_fcn`
- `parse_summary`
- `read_jack`
- `run`
- `set_key`
- `sha256`
- `with`
- `write_lattice_files`
- `write_matrix_csv`
- `zip`

Inputs/outputs: inspect the file CLI and adjacent reports; this index is intentionally mechanical.

## `scripts/run_l20_k3b_k3e_zero_tests.py`

Purpose: Python workflow/utility.

Important functions/classes:
- `except`
- `main`
- `make_config`
- `prefit_audit`
- `run_test`
- `set_key`
- `with`
- `zero_sources`

Inputs/outputs: inspect the file CLI and adjacent reports; this index is intentionally mechanical.

## `scripts/run_l20_user_validated_passonly.py`

Purpose: Python workflow/utility.

Important functions/classes:
- `main`
- `write_fitter_zero_sources`

Inputs/outputs: inspect the file CLI and adjacent reports; this index is intentionally mechanical.

## `scripts/run_v33m_classifier_speed_search.py`

Purpose: Python workflow/utility.

Important functions/classes:
- `assignment_rows`
- `assignment_signature`
- `compare_assignments`
- `main`
- `parse_summary`
- `run`
- `run_benchmark`
- `run_time`
- `write_bench_reports`
- `write_correctness_report`
- `write_dispatch_audit`
- `write_temp_cfg`

Inputs/outputs: inspect the file CLI and adjacent reports; this index is intentionally mechanical.

## `scripts/run_v33m_cpu_subset_validation.py`

Purpose: Python workflow/utility.

Important functions/classes:
- `add_selection`
- `f`
- `i`
- `int`
- `load_scan_rows`
- `main`
- `parse_args`
- `rmax`
- `rmin`
- `run_helper`
- `select_rows`

Inputs/outputs: inspect the file CLI and adjacent reports; this index is intentionally mechanical.

## `scripts/run_v33n_cutoff_scan_fits.py`

Purpose: Python workflow/utility.

Important functions/classes:
- `all_candidates_file`
- `append_report_section`
- `base_config_text`
- `bestfit_qc_file`
- `build_block_index`
- `build_block_index_by_file_irrep`
- `build_qc_block_index`
- `build_qc_outputs`
- `choose_seed_order`
- `clean_warnings`
- `copy_target_dump`
- `ctag`
- `fit_levels_file`
- `fit_summary_file`
- `load_csv_rows`
- `load_nonint_rows`
- `load_runtime_cache_with_params`
- `main`
- `parse_all_candidates`
- `parse_benchmark_output`
- `parse_bestfit_qc`
- `parse_fit_levels`
- `parse_fit_summary`
- `parse_keyvals`
- `parse_nonint_file`
- `parse_target_dump`
- `parse_whitespace_table`
- `plot_det_block`
- `plot_spectrum`
- `plot_spectrum_panel`
- `read_cpx_matrix`
- `read_f64`
- `read_i32`
- `read_runtime_meta`
- `read_u64`
- `replace_line`
- `run_cmd`
- `run_fit_attempt`
- `run_spectrum_only_check`
- `run_target_dump`
- `target_groups`
- `write_benchmark_report`
- `write_config`
- `write_csv`
- `write_dispatch_report`
- `write_json`
- `write_parameter_mask_report`
- `write_summary_csv`
- `write_summary_json`
- `write_target_parity_report`
- `write_target_summary_report`

Inputs/outputs: inspect the file CLI and adjacent reports; this index is intentionally mechanical.

## `scripts/run_v33p_classifier_sweep_L20_100_A2.py`

Purpose: Python workflow/utility.

Important functions/classes:
- `f`
- `linear_zero`
- `main`
- `previous_matrix_row`
- `print`
- `sign`
- `sign_flips`
- `v3_classify`
- `v3_shoulder`

Inputs/outputs: inspect the file CLI and adjacent reports; this index is intentionally mechanical.

## `scripts/run_v33p_projected_basis_hotwindow_parallel.py`

Purpose: Python workflow/utility.

Important functions/classes:
- `main`
- `merge`
- `run_chunk`
- `split_windows`
- `write_chunk`

Inputs/outputs: inspect the file CLI and adjacent reports; this index is intentionally mechanical.

## `scripts/save_l20_passonly_root_artifacts.py`

Purpose: Python workflow/utility.

Important functions/classes:
- `candidate_rows`
- `fit_levels`
- `lattice_sources`
- `load_basis`
- `main`
- `params`
- `root_for_window`
- `windows`
- `with`

Inputs/outputs: inspect the file CLI and adjacent reports; this index is intentionally mechanical.

## `scripts/score_v33p_classifier_against_user_labels.py`

Purpose: Python workflow/utility.

Important functions/classes:
- `SystemExit`
- `load_labels`
- `main`

Inputs/outputs: inspect the file CLI and adjacent reports; this index is intentionally mechanical.

## `scripts/v33p_nscale50_classifier_test_plot.py`

Purpose: Python workflow/utility.

Important functions/classes:
- `as_float`
- `configure_matplotlib_backend`
- `main`
- `read_rows`
- `resolve_cache_file`

Inputs/outputs: inspect the file CLI and adjacent reports; this index is intentionally mechanical.

## `scripts/v33p_plot_full_cache_det_scaled_projected_QC.py`

Purpose: Python workflow/utility.

Important functions/classes:
- `configure_matplotlib_backend`
- `discover_fit_params`
- `load_runtime_cache_with_params`
- `load_runtime_meta`
- `main`
- `parse_summary_file`
- `read_cpx_matrix`
- `read_f64`
- `read_i32`
- `read_rows`
- `read_u64`

Inputs/outputs: inspect the file CLI and adjacent reports; this index is intentionally mechanical.

## `scripts/v33p_plot_full_cache_det_scaled_projected_QC_autok.py`

Purpose: Python workflow/utility.

Important functions/classes:
- `choose_k`
- `configure_matplotlib_backend`
- `discover_fit_params`
- `format_cmd`
- `load_runtime_cache_with_params`
- `load_runtime_meta`
- `main`
- `nscale_cmd`
- `parse_summary_file`
- `percentile`
- `raw_det_stats`
- `read_cpx_matrix`
- `read_f64`
- `read_i32`
- `read_rows`
- `read_u64`
- `resolve_runtime_cache`
- `write_report`

Inputs/outputs: inspect the file CLI and adjacent reports; this index is intentionally mechanical.

## `scripts/v33p_plot_full_cache_det_scaled_projected_QC_helper.cpp`

Purpose: C++/header implementation.

Important functions/classes:
- `for`
- `if`
- `in_window`
- `join_doubles`
- `write_csv`

Inputs/outputs: inspect the file CLI and adjacent reports; this index is intentionally mechanical.

## `scripts/v33p_plot_full_cache_det_scaled_projected_QC_multik.py`

Purpose: Python workflow/utility.

Important functions/classes:
- `configure_matplotlib_backend`
- `load_runtime_cache_with_params`
- `load_runtime_meta`
- `main`
- `parse_summary_file`
- `read_cpx_matrix`
- `read_f64`
- `read_i32`
- `read_rows`
- `read_u64`

Inputs/outputs: inspect the file CLI and adjacent reports; this index is intentionally mechanical.

## `scripts/v33p_plot_fullcache_raw_projectedQC_global_kshift_family.py`

Purpose: Python workflow/utility.

Important functions/classes:
- `build_cmd`
- `configure_matplotlib_backend`
- `discover_fit_params`
- `finite_stats`
- `load_rows`
- `load_runtime_meta`
- `main`
- `parse_summary_file`
- `read_cpx_matrix`
- `read_f64`
- `read_i32`
- `read_rows`
- `read_u64`
- `resolve_runtime_cache`
- `sign_change_count`
- `write_report`

Inputs/outputs: inspect the file CLI and adjacent reports; this index is intentionally mechanical.

## `scripts/v33p_plot_fullcache_scaleQC_Lxi1_multik.py`

Purpose: Python workflow/utility.

Important functions/classes:
- `build_rerun_cmd`
- `configure_matplotlib_backend`
- `discover_fit_params`
- `finite_stats`
- `load_rows`
- `load_runtime_meta`
- `main`
- `parse_summary_file`
- `read_cpx_matrix`
- `read_f64`
- `read_i32`
- `read_rows`
- `read_u64`
- `resolve_runtime_cache`
- `summarize_matrix`
- `write_report`

Inputs/outputs: inspect the file CLI and adjacent reports; this index is intentionally mechanical.

## `scripts/v33p_plot_fullcache_scaled_projected_QC_multik.py`

Purpose: Python workflow/utility.

Important functions/classes:
- `SystemExit`

Inputs/outputs: inspect the file CLI and adjacent reports; this index is intentionally mechanical.

## `scripts/v33p_plot_fullcache_slogdet_c_scaled_projectedQC.py`

Purpose: Python workflow/utility.

Important functions/classes:
- `build_cmd`
- `configure_matplotlib_backend`
- `count_slogdet_signs`
- `discover_fit_params`
- `finite_stats`
- `load_rows`
- `load_runtime_meta`
- `main`
- `parse_summary_file`
- `read_cpx_matrix`
- `read_f64`
- `read_i32`
- `read_rows`
- `read_u64`
- `resolve_runtime_cache`
- `sign_change_count`
- `write_report`

Inputs/outputs: inspect the file CLI and adjacent reports; this index is intentionally mechanical.

## `scripts/validate_fitter_target_counts_v33k.py`

Purpose: Python workflow/utility.

Important functions/classes:
- `main`
- `write_report`

Inputs/outputs: inspect the file CLI and adjacent reports; this index is intentionally mechanical.

## `scripts/verify_classifier_dispatch_v33m.py`

Purpose: Python workflow/utility.

Important functions/classes:
- `main`
- `make_temp_config`
- `run_mode`
- `write_outputs`

Inputs/outputs: inspect the file CLI and adjacent reports; this index is intentionally mechanical.

## `L20_spectrum_plotter_model_true_zeros_yaxis_fixed.py`

Purpose: Python workflow/utility.

Important functions/classes:
- `_parse_bool`
- `_parse_float`
- `find_nonint_script`
- `generate_recent_nonint_levels`
- `main`
- `parse_nonint_output_file`
- `plot_spectrum`
- `print`
- `read_lattice_v20`
- `read_model_true_zeros_csv`
- `read_old_qc_file`
- `read_recent_bestfit_spectrum`
- `resolve_file`
- `resolve_recent_fit_file`
- `safe_usetex`
- `write_nonint_config`

Inputs/outputs: inspect the file CLI and adjacent reports; this index is intentionally mechanical.

## `packaging/generate_function_index.py`

Purpose: Python workflow/utility.

Important functions/classes:
- `functions`
- `sorted`

Inputs/outputs: inspect the file CLI and adjacent reports; this index is intentionally mechanical.

## `tests/run_package_tests.py`

Purpose: Python workflow/utility.

Important functions/classes:
- `main`
- `print`
- `run`

Inputs/outputs: inspect the file CLI and adjacent reports; this index is intentionally mechanical.

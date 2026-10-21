# v33d GPU component cache generator config
# Compile/run through scripts/run_gpu_component_cachegen.sh to save plm_config, klm_config, K2inv, F2, G, F3, F3inv, Vsel, projected_F3inv.
# The run script sources this file. Command-line Lbyas values override any Lbyas setting here.

# Physics and grid
Ecm_min=0.26310
Ecm_max=0.36000
coarseN=10000   # set to 50000 for production
xi=3.444

# Safe default two-body scatter parameters used by the compiled CUDA code:
# scatter1_00 = 4.04
# scatter1_10 = -43.2
# scatter2_00 = 4.12
# waves_vec_1 = {0,1}
# waves_vec_2 = {0}

# Default cache saving place requested by Digonto
output_root=/media/digonto/Data/F3inv_cache

# GPU execution controls
use_fused_f2k2=true
validate_fused_f2k2=false
use_kinematic_precompute=true
use_stream_pipeline=true
# max_batch_energies=0 means auto from available VRAM.
max_batch_energies=0
gpu_memory_safety_fraction=0.85
use_concurrent_streams=true
concurrent_streams_safety_fraction=0.85
use_multi_gpu=false
use_mixed_precision_f2k2=false
irreps=000_A1m,100_A2,110_A2,111_A2,200_A2

default_step_size=0.2

debug=n

# v32zu optimization note: safety fraction raised to 0.85; guarded OPT-1..8 flags keep validated fallback paths available.

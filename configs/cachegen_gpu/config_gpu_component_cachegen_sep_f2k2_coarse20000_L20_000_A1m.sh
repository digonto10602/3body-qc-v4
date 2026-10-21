# v33d GPU component cache generator config
# Narrow run for the validated Lbyas=20 / 000_A1m coarse20000 check.

CONFIG_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=/dev/null
source "$CONFIG_DIR/config_gpu_component_cachegen_sep_f2k2_coarse20000.sh"

irreps=000_A1m
output_root=/media/digonto/Data/F3inv_cache/v33e_coarse20000_sep_f2k2_L20_000_A1m

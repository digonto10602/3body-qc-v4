# v33d GPU component cache generator config
# Dedicated coarse20000 production config under the default F3inv_cache mount.

CONFIG_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=/dev/null
source "$CONFIG_DIR/config_gpu_component_cachegen_sep_f2k2.sh"

coarseN=20000
output_root="${V33E_OUTPUT_ROOT:-/media/digonto/Data/F3inv_cache/v33e_coarse20000_sep_f2k2}"
irreps="${V33E_IRREPS:-000_A1m,100_A2,110_A2,111_A2,200_A2}"

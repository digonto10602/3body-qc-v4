# v33d GPU component cache generator config
# This opt-in variant keeps the default physics/settings but disables the fused F2/K2 kernel.
# It is the currently validated faster wrapper behavior for the L20 / 000_A1m smoke benchmark.

CONFIG_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=/dev/null
source "$CONFIG_DIR/config_gpu_component_cachegen.sh"

use_fused_f2k2=false


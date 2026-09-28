#!/usr/bin/env bash
# Optional wrapper. The documented primary interface is Rscript.
set -euo pipefail
SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
REPO_DIR="$(cd -- "$SCRIPT_DIR/../.." && pwd)"
action="${1:---help}"
if [[ $# -gt 0 ]]; then shift; fi
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
export EIVGP_CORES="${EIVGP_CORES:-4}"
RSCRIPT_BIN="${RSCRIPT_BIN:-Rscript}"
run_r() {
  if [[ "$(uname -s)" == Darwin ]] && command -v caffeinate >/dev/null 2>&1; then
    exec caffeinate -i "$RSCRIPT_BIN" "$@"
  else exec "$RSCRIPT_BIN" "$@"; fi
}
case "$action" in
  adni|ocean|data) run_r "$SCRIPT_DIR/run_real_data.R" "$action" "$@" ;;
  smoke) export EIVGP_SMOKE_TEST=1; run_r "$SCRIPT_DIR/run_real_data.R" adni "$@" ;;
  report) run_r "$SCRIPT_DIR/ADNI_report.R" "$@" ;;
  check) run_r "$SCRIPT_DIR/check_real_data.R" "$@" ;;
  setup)
    mkdir -p "$REPO_DIR/real_data_application_outputs/setup"
    cd "$REPO_DIR/real_data_application_outputs/setup"
    run_r "$SCRIPT_DIR/setup_real_data.R" "$@" ;;
  --help|-h) echo 'Optional: bash codes/applications/run_real_data_mac.sh adni|ocean|data|smoke [--cores N] [--plan]'
    echo 'Also: check, setup, report. ADNI is one fixed five-fold analysis.' ;;
  *) echo "Unknown action: $action" >&2; exit 2 ;;
esac

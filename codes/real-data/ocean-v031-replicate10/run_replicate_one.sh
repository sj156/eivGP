#!/usr/bin/env bash
# Run one frozen Ocean redraw from this bundle.
# Usage: ./run_replicate_one.sh REP_ID N_CALIB [smoke|full] [n_cores] [target_diagnostics 0|1]
set -euo pipefail
bundle="$(cd "$(dirname "$0")" && pwd)"
rep_id="${1:?REP_ID 1-10}"
n_calib="${2:?N_CALIB 10, 25, or 50}"
mode="${3:-full}"
n_cores="${4:-4}"
target_diagnostics="${5:-0}"

ocean_data="${OCEAN_DATA_DIR:-$bundle/../ocean_data}"
if [[ ! -f "$ocean_data/study1_data.csv" ]]; then
  echo "Missing $ocean_data/study1_data.csv. Keep this folder next to ocean_data/ or set OCEAN_DATA_DIR." >&2
  exit 1
fi

exec Rscript "$bundle/01_run_replicate_one.R" \
  "$ocean_data" \
  "$bundle/calib_designs" \
  "$bundle/outputs/replicate10" \
  "$rep_id" "$n_calib" "$mode" "$target_diagnostics" "$n_cores"

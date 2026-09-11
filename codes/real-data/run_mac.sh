#!/usr/bin/env bash
# macOS/Unix launcher: set thread budgets before starting R.
set -euo pipefail
SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
REPO_DIR="$(cd -- "$SCRIPT_DIR/../.." && pwd)"
action="${1:---help}"
if [[ $# -gt 0 ]]; then shift; fi
core_budget="${EIVGP_CORES:-2}"
repeat_spec="${ADNI_REPEATS:-${ADNI_REPEAT_ID:-2}}"
while [[ $# -gt 0 ]]; do
  case "$1" in
    --cores)
      if [[ $# -lt 2 ]]; then echo '--cores requires a positive integer.' >&2; exit 2; fi
      core_budget="$2"; shift 2 ;;
    --cores=*) core_budget="${1#--cores=}"; shift ;;
    --repeats)
      if [[ $# -lt 2 ]]; then echo '--repeats requires 2, 3, or 2,3.' >&2; exit 2; fi
      repeat_spec="$2"; shift 2 ;;
    --repeats=*) repeat_spec="${1#--repeats=}"; shift ;;
    *) echo "Unknown option: $1" >&2; exit 2 ;;
  esac
done
case "$action" in
  --help|-h) cat <<'HELP'
Usage: bash codes/real-data/run_mac.sh {check|setup|smoke|adni|ocean|report} [--cores N] [--repeats 2,3]
  check   Validate prepared data and print runtime/package availability; no fitting.
  setup   Install dependencies and tested eivGP 0.3.1 (requires internet).
  smoke   Short ADNI pipeline test, isolated from production output.
  adni    Production ADNI run; default repeat 2, four chains, two workers.
  ocean   Ocean validation; default quick mode. Set OCEAN_QUICK=false for production.
  report  Combine completed ADNI repeats 2 and 3; no fitting.
Overrides: ADNI_DATA_DIR, OCEAN_DATA_DIR, ADNI_OUTPUT_DIR, OCEAN_OUTPUT_DIR,
           ADNI_REPEAT_ID, EIVGP_CORES (default budget: 2), RSCRIPT_BIN.
--cores N overrides EIVGP_CORES and controls both analysis worker settings.
See MAC_MINI.md for foreground/background commands and checkpoint resumption.
HELP
    exit 0 ;;
  check|setup|smoke|adni|ocean|report) ;;
  *) echo "Unknown action: $action" >&2; exit 2 ;;
esac
RSCRIPT_BIN="${RSCRIPT_BIN:-Rscript}"
if ! command -v "$RSCRIPT_BIN" >/dev/null 2>&1; then
  echo 'Rscript was not found. Install R or set RSCRIPT_BIN to its full path.' >&2
  exit 1
fi
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
if [[ ! "$core_budget" =~ ^[1-9][0-9]*$ ]] || [[ ${#core_budget} -gt 6 ]]; then
  echo '--cores must be a positive integer no larger than 999999.' >&2
  exit 2
fi
case "$repeat_spec" in 2|3|2,3|3,2) ;; *) echo '--repeats must be 2, 3, or 2,3.' >&2; exit 2 ;; esac
export ADNI_REPEATS="$repeat_spec"
export EIVGP_CORES="$core_budget"
# The notebook divides this total budget across folds and chains.
# Keep legacy chain settings as a fallback for older code.
chain_workers="$core_budget"
if [[ "$chain_workers" -gt 4 ]]; then chain_workers=4; fi
export EIVGP_N_CORES="$chain_workers"
export MIXEDGP_CORES="$core_budget"
if [[ -z "${ADNI_DATA_DIR:-}" ]]; then
  if [[ -f "$REPO_DIR/real-data/adni/toledo_adni_cohort_n495.csv" ]]; then
    export ADNI_DATA_DIR="$REPO_DIR/real-data/adni"
  elif [[ -f "$REPO_DIR/real-data/toledo_adni_cohort_n495.csv" ]]; then
    export ADNI_DATA_DIR="$REPO_DIR/real-data"
  else
    export ADNI_DATA_DIR="$SCRIPT_DIR/ADNI-toledo/data"
  fi
fi
export OCEAN_DATA_DIR="${OCEAN_DATA_DIR:-$REPO_DIR/real-data/ocean}"
export ADNI_OUTPUT_DIR="${ADNI_OUTPUT_DIR:-$REPO_DIR/results/adni}"
export OCEAN_OUTPUT_DIR="${OCEAN_OUTPUT_DIR:-$REPO_DIR/results/ocean}"
run_r() {
  if [[ "$(uname -s)" == Darwin ]] && command -v caffeinate >/dev/null 2>&1; then
    exec caffeinate -i "$RSCRIPT_BIN" "$@"
  else
    exec "$RSCRIPT_BIN" "$@"
  fi
}
case "$action" in
  setup)
    mkdir -p "$REPO_DIR/results/setup"
    cd "$REPO_DIR/results/setup"
    run_r "$SCRIPT_DIR/setup.R" ;;
  check) run_r "$SCRIPT_DIR/check_runtime.R" ;;
  smoke)
    export EIVGP_SMOKE_TEST=1
    run_r "$SCRIPT_DIR/run_application.R" adni --cores "$core_budget" ;;
  adni)
    export EIVGP_SMOKE_TEST=0
    run_r "$SCRIPT_DIR/run_application.R" adni --cores "$core_budget" ;;
  ocean) run_r "$SCRIPT_DIR/run_application.R" ocean --cores "$core_budget" ;;
  report) run_r "$SCRIPT_DIR/ADNI-toledo/report_adni.R" ;;
esac

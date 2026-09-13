#!/usr/bin/env bash
# macOS/Unix launcher: set thread budgets before starting R.
set -euo pipefail
SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
REPO_DIR="$(cd -- "$SCRIPT_DIR/../.." && pwd)"
action="${1:---help}"
if [[ $# -gt 0 ]]; then shift; fi
core_budget="${EIVGP_CORES:-2}"
validation_default="1"
if [[ "$action" == "report" ]]; then validation_default="1-20"; fi
validation_spec="${ADNI_VALIDATIONS:-${ADNI_VALIDATION_ID:-$validation_default}}"

while [[ $# -gt 0 ]]; do
  case "$1" in
    --cores)
      if [[ $# -lt 2 ]]; then echo '--cores requires a positive integer.' >&2; exit 2; fi
      core_budget="$2"; shift 2 ;;
    --cores=*) core_budget="${1#--cores=}"; shift ;;
    --validations)
      if [[ $# -lt 2 ]]; then echo '--validations requires an ID selection.' >&2; exit 2; fi
      validation_spec="$2"; shift 2 ;;
    --validations=*) validation_spec="${1#--validations=}"; shift ;;
    *) echo "Unknown option: $1" >&2; exit 2 ;;
  esac
done

case "$action" in
  --help|-h) cat <<'HELP'
Usage: bash codes/real-data/run_mac.sh {check|setup|smoke|adni|ocean|report} [--cores N] [--validations 1-20]
  check   Validate prepared data and print runtime/package availability; no fitting.
  setup   Install dependencies and tested eivGP 0.3.1 (requires internet).
  smoke   Short ADNI pipeline test, isolated from production output.
  adni    Production ADNI run for one or more validation IDs.
  ocean   Ocean validation; default quick mode. Set OCEAN_QUICK=false for production.
  report  Combine completed ADNI validations; no fitting.
Examples:
  bash codes/real-data/run_mac.sh smoke --cores 2 --validations 10
  bash codes/real-data/run_mac.sh adni --cores 4 --validations 10-12
  bash codes/real-data/run_mac.sh report --validations 1-20
Overrides: ADNI_DATA_DIR, OCEAN_DATA_DIR, ADNI_OUTPUT_DIR, OCEAN_OUTPUT_DIR,
           ADNI_VALIDATIONS, ADNI_VALIDATION_ID, EIVGP_CORES, RSCRIPT_BIN.
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
if [[ -z "$validation_spec" ]]; then echo '--validations cannot be empty.' >&2; exit 2; fi
export ADNI_VALIDATIONS="$validation_spec"
unset ADNI_REPEATS ADNI_REPEAT_ID || true
export EIVGP_CORES="$core_budget"
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
    run_r "$SCRIPT_DIR/run_application.R" adni --cores "$core_budget" --validations "$validation_spec" ;;
  adni)
    export EIVGP_SMOKE_TEST=0
    run_r "$SCRIPT_DIR/run_application.R" adni --cores "$core_budget" --validations "$validation_spec" ;;
  ocean)
    unset ADNI_VALIDATIONS ADNI_VALIDATION_ID || true
    run_r "$SCRIPT_DIR/run_application.R" ocean --cores "$core_budget" ;;
  report) run_r "$SCRIPT_DIR/ADNI-toledo/report_adni.R" --validations "$validation_spec" ;;
esac

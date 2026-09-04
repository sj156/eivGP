############################################################
## Master script: Study II
##
## Safe default:
##   Rscript run_study2_simulation.R
##
## Workstation publication run:
##   MIXEDGP_RUN_MODE=publication MIXEDGP_WORKERS=8 \
##     Rscript run_study2_simulation.R
############################################################

master_args <- commandArgs(FALSE)
master_file_arg <- grep("^--file=", master_args, value = TRUE)
MASTER_CODE_DIR <- if (length(master_file_arg) > 0L) {
  dirname(normalizePath(
    sub("^--file=", "", master_file_arg[[1L]]),
    winslash = "/", mustWork = TRUE
  ))
} else {
  normalizePath(getwd(), winslash = "/", mustWork = TRUE)
}
source(file.path(MASTER_CODE_DIR, "reproduction_workflows.R"))

RUN_MODE <- Sys.getenv("MIXEDGP_RUN_MODE", unset = "dry_run")
WORKERS <- suppressWarnings(as.integer(Sys.getenv(
  "MIXEDGP_WORKERS", unset = as.character(mixedgp_available_workers())
)))
OUTPUT_ROOT <- Sys.getenv("MIXEDGP_OUTPUT_ROOT", unset = "")
DATA_ROOT <- Sys.getenv("MIXEDGP_DATA_ROOT", unset = "")

## Single user-editable configuration block. The publication constructor
## freezes the q progression, calibration curve, focused controls, methods,
## MCMC controls, seeds, and reportable run sizes.
CONFIG <- study2_simulation_config(
  mode = RUN_MODE,
  code_dir = MASTER_CODE_DIR,
  workers = WORKERS,
  output_root = if (nzchar(OUTPUT_ROOT)) OUTPUT_ROOT else NULL,
  data_root = if (nzchar(DATA_ROOT)) DATA_ROOT else NULL
)

stage_override <- Sys.getenv("MIXEDGP_STAGES", unset = "")
if (nzchar(stage_override)) {
  if (RUN_MODE %in% c("dry_run", "data")) {
    stop(
      "MIXEDGP_STAGES may be overridden only with ",
      "MIXEDGP_RUN_MODE=smoke or publication."
    )
  }
  CONFIG$stages <- trimws(strsplit(stage_override, ",", fixed = TRUE)[[1L]])
  CONFIG <- validate_simulation_config(CONFIG)
}

STUDY2_RUN <- run_study2_simulation(CONFIG)

############################################################
## run_ocean_all.R
##
## Master script for the ocean real-data case.
## Same pattern as new-codes/run_study1_all.R.
############################################################

rm(list = ls())

file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
ocean_runner_file <- if (length(file_arg) > 0L) {
  normalizePath(sub("^--file=", "", file_arg[1L]), mustWork = TRUE)
} else {
  normalizePath(sys.frame(1)$ofile, mustWork = TRUE)
}
OCEAN_REALDATA_DIR <- dirname(ocean_runner_file)

## Set TRUE for a short test run.
## Set FALSE for paper-quality runs (4 x 20k, 10 redraws).
OCEAN_QUICK <- TRUE

## Use cached .rds results if they exist.
OCEAN_USE_CACHE <- TRUE

## Use the same covariance family for EIV-GP and every GP comparator.
OCEAN_KERNEL <- "se"              # "se" or "matern"
OCEAN_MATERN_NU <- 2.5
OCEAN_CACHE_VERSION <- "v2_same-target"

## Output prefix. Scripts live in new-codes/real_data/, so ".."
## writes figures/tables/results next to new-codes/.
OCEAN_OUT_PREFIX <- normalizePath(
  file.path(OCEAN_REALDATA_DIR, "..", ".."), mustWork = FALSE
)

source(file.path(OCEAN_REALDATA_DIR, "..", "00_study1_functions.R"))
source(file.path(OCEAN_REALDATA_DIR, "ocean_data_helpers.R"))

cat("\nRunning representative ocean figures...\n")
source(file.path(OCEAN_REALDATA_DIR, "01_ocean_representative_figures.R"))

cat("\nRunning ocean nested-calibration redraws...\n")
source(file.path(OCEAN_REALDATA_DIR, "02_ocean_replicates.R"))

cat("\nDone.\n")

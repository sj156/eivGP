############################################################
## run_ocean_all.R
##
## Master script for the ocean real-data case.
## Same pattern as new-codes/run_study1_all.R.
############################################################



if (!exists("OCEAN_REALDATA_DIR")) {
file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
ocean_runner_file <- if (length(file_arg) > 0L) {
  normalizePath(sub("^--file=", "", file_arg[1L]), mustWork = TRUE)
} else {
  normalizePath(sys.frame(1)$ofile, mustWork = TRUE)
}
OCEAN_REALDATA_DIR <- dirname(ocean_runner_file)

}

## Set TRUE for a short test run.
## Set FALSE for paper-quality runs (4 x 20k, 10 redraws).
env_flag <- function(name, default) {
  value <- tolower(Sys.getenv(name, default))
  if (!value %in% c("true", "false", "1", "0")) stop(name, " must be true/false or 1/0")
  value %in% c("true", "1")
}
OCEAN_QUICK <- env_flag("OCEAN_QUICK", "true")

## Use cached .rds results if they exist.
OCEAN_USE_CACHE <- env_flag("OCEAN_USE_CACHE", "true")

## Use the same covariance family for EIV-GP and every GP comparator.
OCEAN_KERNEL <- Sys.getenv("OCEAN_KERNEL", "se")              # "se" or "matern"
OCEAN_MATERN_NU <- as.numeric(Sys.getenv("OCEAN_MATERN_NU", "2.5"))
stopifnot(OCEAN_KERNEL %in% c("se", "matern"), OCEAN_MATERN_NU %in% c(0.5, 1.5, 2.5))
OCEAN_CACHE_VERSION <- "v3_eivGP031"

OCEAN_OUT_PREFIX <- Sys.getenv("OCEAN_OUTPUT_DIR", file.path(OCEAN_REALDATA_DIR, "outputs"))

source(file.path(OCEAN_REALDATA_DIR, "shared", "require_eivgp.R"))
source(file.path(OCEAN_REALDATA_DIR, "shared", "00_study1_functions.R"))
source(file.path(OCEAN_REALDATA_DIR, "ocean_data_helpers.R"))

cat("\nRunning representative ocean figures...\n")
source(file.path(OCEAN_REALDATA_DIR, "01_ocean_representative_figures.R"))

cat("\nRunning ocean nested-calibration redraws...\n")
source(file.path(OCEAN_REALDATA_DIR, "02_ocean_replicates.R"))

cat("\nDone.\n")

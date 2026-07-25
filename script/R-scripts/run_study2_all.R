############################################################
## run_study2_all.R
##
## Master script for revised Study II.
##
## Fully Bayesian EIV-GP with unknown ordinal-probit
## measurement model.
############################################################

rm(list = ls())

############################################################
## Main knobs
############################################################

## Choose one of:
##   "quick"    : short run for code checking
##   "balanced" : moderate run for development
##   "thorough" : paper-quality run
STUDY2_CONFIG <- "quick"

## Use cached .rds results if they exist.
STUDY2_USE_CACHE <- FALSE

## Output prefix.
## If scripts are in a folder such as "code/", this writes to "../figures".
## Change to "." if you want output in the current directory.
STUDY2_OUT_PREFIX <- ".."

source("00_study2_functions.R")

cat("\nRunning representative Study II figures...\n")
source("01_study2_representative_figures.R")

cat("\nRunning Study II Monte Carlo comparison...\n")
source("02_study2_monte_carlo.R")

cat("\nDone.\n")
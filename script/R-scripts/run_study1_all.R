############################################################
## run_study1_all.R
##
## Master script for revised Study I.
############################################################

rm(list = ls())

## Set TRUE for a short test run.
## Set FALSE for paper-quality runs.
STUDY1_QUICK <- FALSE

## Use cached .rds results if they exist.
STUDY1_USE_CACHE <- FALSE

## Output prefix.
## If scripts are in a folder such as "code/", this writes to "../figures".
## Change to "." if you want figures/tables/results in the current directory.
STUDY1_OUT_PREFIX <- ".."

source("00_study1_functions.R")

cat("\nRunning representative Study I figures...\n")
source("01_study1_representative_figures.R")

cat("\nRunning Study I Monte Carlo comparison...\n")
source("02_study1_monte_carlo.R")

cat("\nDone.\n")
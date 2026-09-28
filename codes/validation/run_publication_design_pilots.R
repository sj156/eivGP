############################################################
## Run both redesigned numerical-study pilots.
## Run from revision/codes/. These outputs are diagnostics, not final results.
############################################################

rm(list = ls())
source("setup/00_project_setup.R")

## Optional overrides can be defined here before sourcing the pilot scripts.
## PILOT_N_REP <- 2L
## PILOT_N_TEST <- 250L
## PILOT_N_DRAW <- 250L
## PILOT_N_ITER <- 900L
## PILOT_BURN <- 350L
## PILOT_N_CHAINS <- 2L

cat("Running redesigned Study I pilot.\n")
source("validation/10_study1_design_pilot.R")

cat("Running redesigned Study II pilot.\n")
source("validation/11_study2_design_pilot.R")

cat("Both publication-design pilots completed.\n")

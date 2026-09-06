############################################################
## 09_experiment_design_validation.R
##
## Fast invariance checks for the targeted Study I controls and paired
## Study II scenarios. This script checks the experimental design only; it
## does not fit a model or replace the Monte Carlo publication run.
############################################################

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
codes_dir <- if (length(script_arg) > 0L) {
  dirname(normalizePath(sub("^--file=", "", script_arg[[1L]])))
} else {
  getwd()
}

source(file.path(codes_dir, "00_parallel_utils.R"))
source(file.path(codes_dir, "00_study1_functions.R"))
source(file.path(codes_dir, "00_study2_functions.R"))
source(file.path(codes_dir, "00_synthetic_data.R"))

############################################################
## Study I: the within-level heterogeneity continuum is exactly paired
############################################################

study1_primary <- simulate_1d_data(
  n = 100L,
  n_test = 200L,
  scenario = "heterogeneity_continuum",
  heterogeneity_eta = 1,
  threshold_design = "balanced",
  seed = 100001L
)
study1_control <- simulate_1d_data(
  n = 100L,
  n_test = 200L,
  scenario = "heterogeneity_continuum",
  heterogeneity_eta = 0,
  threshold_design = "balanced",
  seed = 100001L
)
study1_mid <- simulate_1d_data(
  n = 100L,
  n_test = 200L,
  scenario = "heterogeneity_continuum",
  heterogeneity_eta = 0.5,
  threshold_design = "balanced",
  seed = 100001L
)

stopifnot(
  identical(study1_primary$train$x, study1_control$train$x),
  identical(study1_primary$train$u, study1_control$train$u),
  identical(study1_primary$train$c, study1_control$train$c),
  identical(study1_primary$test$x, study1_control$test$x),
  identical(study1_primary$test$u, study1_control$test$u),
  identical(study1_primary$test$c, study1_control$test$c),
  isTRUE(all.equal(
    study1_primary$train$y - study1_primary$train$f,
    study1_control$train$y - study1_control$train$f
  )),
  isTRUE(all.equal(
    study1_primary$test$y - study1_primary$test$f,
    study1_control$test$y - study1_control$test$f
  )),
  isTRUE(all.equal(
    study1_mid$train$f,
    0.5 * (study1_primary$train$f + study1_control$train$f)
  )),
  isTRUE(all.equal(
    study1_mid$test$f,
    0.5 * (study1_primary$test$f + study1_control$test$f)
  ))
)

invalid_calibration <- try(
  make_nested_calibration_sets(30L, c(0L, 31L), seed = 100002L),
  silent = TRUE
)
stopifnot(inherits(invalid_calibration, "try-error"))

study1_balanced <- simulate_1d_data(
  n = 100L,
  n_test = 200L,
  scenario = "active",
  threshold_design = "balanced",
  seed = 100001L
)
stopifnot(
  identical(study1_balanced$min_class_count, 0L),
  isTRUE(all.equal(
    study1_balanced$tau_true,
    qnorm((1:5) / 6)
  ))
)

############################################################
## Study II: scenarios and nested proxy dimensions share primitives
############################################################

study2_seed <- 1010001L
study2_primary <- simulate_study2_data(
  n = 40L,
  n_test = 60L,
  scenario = "primary",
  seed = study2_seed
)
study2_additive <- simulate_study2_data(
  n = 40L,
  n_test = 60L,
  scenario = "latent_additive_control",
  seed = study2_seed
)
study2_uncertain <- simulate_study2_data(
  n = 40L,
  n_test = 60L,
  scenario = "high_uncertainty",
  seed = study2_seed
)
study2_q2 <- simulate_study2_data(
  n = 40L, n_test = 60L, scenario = "primary", q = 2L,
  seed = study2_seed
)
study2_q6 <- simulate_study2_data(
  n = 40L, n_test = 60L, scenario = "primary", q = 6L,
  seed = study2_seed
)

stopifnot(
  identical(study2_primary$train$X, study2_additive$train$X),
  identical(study2_primary$train$U, study2_additive$train$U),
  identical(study2_primary$train$C, study2_additive$train$C),
  identical(study2_primary$test$X, study2_additive$test$X),
  identical(study2_primary$test$U, study2_additive$test$U),
  identical(study2_primary$test$C, study2_additive$test$C),
  isTRUE(all.equal(
    study2_primary$train$y - study2_primary$train$f,
    study2_additive$train$y - study2_additive$train$f
  )),
  isTRUE(all.equal(
    study2_primary$test$y - study2_primary$test$f,
    study2_additive$test$y - study2_additive$test$f
  )),
  identical(study2_primary$train$X, study2_uncertain$train$X),
  identical(study2_primary$train$U, study2_uncertain$train$U),
  identical(study2_primary$test$X, study2_uncertain$test$X),
  identical(study2_primary$test$U, study2_uncertain$test$U),
  isTRUE(all.equal(
    study2_primary$train$y - study2_primary$train$f,
    study2_uncertain$train$y - study2_uncertain$train$f
  )),
  isTRUE(all.equal(
    study2_primary$test$y - study2_primary$test$f,
    study2_uncertain$test$y - study2_uncertain$test$f
  )),
  identical(study2_q2$train$X, study2_primary$train$X),
  identical(study2_q2$train$U, study2_primary$train$U),
  identical(study2_q2$train$y, study2_primary$train$y),
  identical(study2_q2$train$C, study2_primary$train$C[, 1:2, drop = FALSE]),
  identical(study2_primary$train$C, study2_q6$train$C[, 1:4, drop = FALSE]),
  identical(study2_q2$test$X, study2_q6$test$X),
  identical(study2_q2$test$U, study2_q6$test$U),
  identical(study2_q2$test$y, study2_q6$test$y)
)

cat("Experiment-design invariance checks passed.\n")

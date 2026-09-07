source("codes/simulation_helpers.R")
code_dir <- normalizePath("codes")
engine <- mixedgp_simulation_engine(code_dir)
expected_ids <- list(
  study1 = c("eta0_balanced", "eta1_balanced"),
  study2 = c("primary_q2", "primary_q4_calibration", "logistic_q4"))
fit_counts <- c(study1 = 300L, study2 = 300L)
for (study in names(expected_ids)) {
  ctor <- get(paste0(study, "_simulation_config"))
  publication <- ctor("publication", code_dir = code_dir)
  development <- ctor("development", code_dir = code_dir)
  stopifnot(identical(vapply(publication$cells, `[[`, "", "id"), expected_ids[[study]]))
  expected_grids <- if (study == "study1") rep(list(c(0L, 20L, 50L)), 2L) else
    list(50L, c(0L, 20L, 50L, 80L), 50L)
  for (i in seq_along(publication$cells)) {
    p <- publication$cells[[i]]
    d <- development$cells[[i]]
    stopifnot(p$n == 100L, p$n_test == 100L, p$n_rep == 50L, d$n_rep == 3L,
              identical(p$calibration_grid, expected_grids[[i]]))
    p$n_rep <- d$n_rep <- NULL
    stopifnot(identical(p, d))
  }
  stopifnot(sum(vapply(publication$cells, function(x)
    x$n_rep * length(x$calibration_grid), integer(1))) == fit_counts[[study]])
  ## Data-only test: first replication for every setting in a temporary root.
  ## No MCMC or existing frozen files are touched.
  cfg <- development
  cfg$data_root <- tempfile(paste0(study, "-design-test-"))
  dir.create(cfg$data_root)
  cfg$cells <- lapply(cfg$cells, function(x) {x$n_rep <- 1L; x})
  cfg$parallel$workers <- 1L
  for (cell in cfg$cells) {
    mixedgp_generate_cell_data(cfg, cell, engine)
    mixedgp_verify_cell_data(cfg, cell, engine)
    frozen <- mixedgp_read_cell_replication(cfg, cell, 1L)
    sets <- frozen$calibration_sets
    stopifnot(identical(as.integer(names(sets)), cell$calibration_grid))
    if (length(sets) > 1L) for (j in seq_len(length(sets) - 1L)) {
      stopifnot(all(sets[[j]] %in% sets[[j + 1L]]))
    }
  }
  checks <- mixedgp_validate_common_random_numbers(cfg, engine)
  stopifnot(nrow(checks) > 0L, all(checks$pass))
}
message("Five settings, 100/100 sizes, 50 versus 3 replications, nested grids and paired data validated.")

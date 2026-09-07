## Run from the repository root: Rscript codes/tests/test_core_budget_drivers.R
source("codes/simulation_helpers.R")
source("codes/00_parallel_utils.R")
for (study in c("study1", "study2")) {
  constructor <- get(paste0(study, "_simulation_config"))
  for (cores in c(1L, 3L, 12L, 16L)) {
    cfg <- constructor("development", code_dir = "codes", core_budget = cores)
    stopifnot(cfg$mcmc$n_iter == 1750L, cfg$mcmc$burn == 500L,
              cfg$mcmc$thin == 1L, cfg$mcmc$n_chains == 4L,
              cfg$parallel$workers * cfg$parallel$chain_workers <= cores)
    controls <- get(paste0("mixedgp_cell_controls_", study))(cfg, cfg$cells[[1]], tempdir())
    prefix <- toupper(study)
    stopifnot(controls[[paste0(prefix, "_CHAIN_WORKERS")]] == min(cores, 4L),
              controls[[paste0(prefix, "_PARALLEL_LEVEL")]] == "hybrid")
  }
  publication <- constructor("publication", code_dir = "codes", core_budget = 12L)
  stopifnot(!publication$strict_competitors, !publication$fail_closed,
            !publication$mcmc$require_gate)
}
if (.Platform$OS.type != "windows") {
  parent <- Sys.getpid()
  result <- mixedgp_parallel_lapply(1:2, function(i) {
    dataset_pid <- Sys.getpid()
    children <- mixedgp_parallel_lapply(1:2, function(j) Sys.getpid(),
                                        n_cores = 2L, seeds = 10L + 1:2)
    c(dataset_pid, unlist(children))
  }, n_cores = 2L, seeds = 1:2)
  stopifnot(all(vapply(result, function(ids) {
    length(unique(ids)) == 3L && !parent %in% ids
  }, logical(1))))
}
message("Both driver allocations and actual nested fork execution passed.")

## Execute the actual diagnostic-warning branches without costly simulation.
find_if <- function(x, condition) {
  if (is.call(x) && identical(x[[1L]], as.name("if")) &&
      identical(paste(deparse(x[[2L]]), collapse = " "), condition)) return(list(x))
  if (!is.call(x) && !is.expression(x) && !is.pairlist(x)) return(list())
  unlist(lapply(as.list(x), find_if, condition = condition), recursive = FALSE)
}
warning_eval <- function(expr, env) {
  notices <- character()
  withCallingHandlers(eval(expr, env), warning = function(w) {
    notices <<- c(notices, conditionMessage(w))
    invokeRestart("muffleWarning")
  })
  stopifnot(length(notices) > 0L, any(grepl("Inspect diagnostic", notices)))
}
s1 <- parse("codes/02_study1_monte_carlo.R")
s2 <- parse("codes/02_study2_monte_carlo.R")
env <- new.env(parent = globalenv())
env$gate_pass <- FALSE
env$rep_id <- 1L
env$n_calib <- 10L
env$RES_DIR <- tempdir()
env$fit_eiv <- list(retained = 1:5)
env$parameter_diag <- data.frame(rhat = 1.3)
env$panel <- 1L
env$STUDY1_CACHE_SPEC <- list()
env$max_rhat <- 1.3
env$min_ess <- 2
warning_eval(find_if(s1, "!isTRUE(gate_pass)")[[1L]], env)
stopifnot(identical(readRDS(env$failure_file)$fit$retained, 1:5))

env$diag_row <- list(mcmc_pass = FALSE)
env$fit <- list(retained = 1:5)
env$scenario <- "primary"
warning_eval(find_if(s2, "!isTRUE(diag_row$mcmc_pass)")[[1L]], env)
env$measurement_warning <- TRUE
env$measurement_fit <- env$fit
env$measurement_diag <- list(convergence_pass = FALSE)
env$measurement_advice <- mixedgp_simulation_diagnostic_advice()
measurement_blocks <- Filter(function(x) any(grepl("warning\\(", deparse(x))),
                             find_if(s2, "measurement_warning"))
warning_eval(measurement_blocks[[1L]], env)
stopifnot(!any(grepl("failed_convergence", readLines("codes/02_study2_monte_carlo.R"))))

for (study in c("study1", "study2")) {
  ctor <- get(paste0(study, "_simulation_config"))
  for (mode in c("development", "publication")) {
    cfg <- ctor(mode, code_dir = "codes", core_budget = 12L)
    stopifnot(!cfg$strict_competitors, !cfg$fail_closed, !cfg$mcmc$require_gate)
    stopifnot(all(vapply(cfg$cells, function(cell) {
      cell$n_test == 100L && cell$n == 100L
    }, logical(1))))
  }
  bad <- ctor("publication", code_dir = "codes")
  bad$mcmc$burn <- bad$mcmc$n_iter
  stopifnot(inherits(try(validate_simulation_config(bad), silent = TRUE), "try-error"))
}
message("Actual diagnostic branches warn and preserve fits; both policies match; invalid inputs still stop.")

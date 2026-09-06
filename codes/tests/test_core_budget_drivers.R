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
  stopifnot(publication$strict_competitors, publication$fail_closed,
            publication$mcmc$require_gate)
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

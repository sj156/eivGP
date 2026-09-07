## Configuration-only regression test; does not generate data or fit models.
## Run from the repository root: Rscript experiments/test_publication_budget.R
source(file.path("codes", "simulation_helpers.R"))
for (constructor in list(study1_simulation_config, study2_simulation_config)) {
  config <- constructor(code_dir = "codes", mode = "publication", workers = 1L)
  stopifnot(config$mcmc$n_iter == 20000L,
            config$mcmc$burn == 5000L,
            config$mcmc$n_chains == 4L,
            config$mcmc$thin == 1L,
            (config$mcmc$n_iter - config$mcmc$burn) * config$mcmc$n_chains == 60000L,
            !config$mcmc$require_gate, !config$fail_closed)
  smoke <- constructor(code_dir = "codes", mode = "smoke", workers = 1L)
  stopifnot(smoke$mcmc$n_iter == 120L, smoke$mcmc$burn == 40L)
}
cat("Publication budgets, no thinning, warning-only diagnostics, and smoke controls passed.\n")

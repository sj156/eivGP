source("codes/00_diagnostics.R")
source("codes/simulation_helpers.R")
source("codes/experiment_reporting.R")
root <- tempfile("reporting-test-"); dir.create(root)
config <- study1_simulation_config("development", code_dir = normalizePath("codes"),
  core_budget = 1L, output_root = root, data_root = file.path(root, "data"))
config$cells <- config$cells[1L]
config$cells[[1]]$n_rep <- 1L
config$cells[[1]]$n <- 30L
config$cells[[1]]$n_test <- 8L
config$cells[[1]]$calibration_grid <- 0L
config$cells[[1]]$run_ablations <- FALSE
config$cells[[1]]$evaluate_f <- config$cells[[1]]$evaluate_u <- FALSE
config$mcmc$n_iter <- 12L; config$mcmc$burn <- 4L; config$mcmc$n_chains <- 2L
config$evaluation <- list(n_pred_draw=4L, n_m_eval=2L, n_m_draw=4L, n_m_latent=4L)
engine <- mixedgp_simulation_engine(config$code_dir)
## Avoid competitor optimization in this narrow reporting test.
engine$run_study1_published_competitors <- function(...) list(draws=list(), latent_means=list(),
  status=data.frame(method="EzGP", status="unavailable", message="test fixture", elapsed_seconds=0))
mixedgp_generate_cell_data(config, config$cells[[1]], engine)
dir.create(file.path(root,"config")); saveRDS(config,file.path(root,"config","resolved_config.rds"))
mixedgp_run_study1_cell(config, config$cells[[1]], engine, root)
stopifnot(!length(list.files(root, pattern="[.](pdf|png)$", recursive=TRUE)))
checkpoint <- file.path(root,"cells",config$cells[[1]]$id,"report_inputs.rds")
stopifnot(file.exists(checkpoint))
hash <- tools::md5sum(checkpoint)
## Trap fitting entry points in every reporting engine: any call fails the test.
original_engine <- mixedgp_simulation_engine
mixedgp_simulation_engine <- function(code_dir) {
  e <- original_engine(code_dir)
  for (name in ls(e)) if (grepl("^(fit_|run_one_|mixedgp_v030_fit)",name))
    assign(name,function(...) stop("Reporting attempted to fit!"),e)
  e
}
for (action in c("summarize", "plot")) {
  status <- mixedgp_report_run(root, action, config$code_dir)
  stopifnot(all(status$status == "success"), identical(tools::md5sum(checkpoint),hash))
  if (action == "summarize") stopifnot(!length(list.files(file.path(root,"reporting"),pattern="[.](pdf|png)$",recursive=TRUE)))
}
stopifnot(length(list.files(file.path(root,"reporting"),pattern="[.]pdf$",recursive=TRUE)) > 0L)
message("Fit-only wrote no figures; summaries and plots replayed saved data without fitting.")

############################################################
## Development loader
##
## Experiment scripts source this file. A future package build copies the
## same canonical modules into R/ and does not need this loader.
############################################################

mixedgp_code_dir <- function() {
  candidates <- c(".", "codes", file.path("revision", "codes"))
  hit <- candidates[file.exists(file.path(candidates, "00_parallel_utils.R"))]
  if (length(hit) == 0L) stop("Cannot locate the mixed-input GP code directory.")
  normalizePath(hit[1L], winslash = "/")
}

mixedgp_source_core <- function(code_dir = mixedgp_code_dir(),
                                include_competitors = TRUE,
                                include_synthetic = TRUE) {
  modules <- c(
    "00_parallel_utils.R",
    "00_study1_functions.R",
    "00_study2_functions.R"
  )
  if (isTRUE(include_competitors)) {
    modules <- c(modules, "03_study2_published_competitors.R")
  }
  if (isTRUE(include_synthetic)) {
    modules <- c(modules, "00_synthetic_data.R")
  }
  modules <- c(modules, "00_public_api.R", "00_diagnostics.R", "00_mcmc_workflow.R")
  modules <- c(modules, "00_experiment_runner.R")
  for (module in modules) {
    sys.source(file.path(code_dir, module), envir = .GlobalEnv)
  }
  invisible(modules)
}

mixedgp_source_core()

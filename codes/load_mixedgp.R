############################################################
## Development loader
##
## Experiment scripts source this file. The eivGP package build copies the
## same canonical model modules into R/ and does not need this loader.
############################################################

mixedgp_code_dir <- function() {
  candidates <- c(".", "codes", file.path("revision", "codes"))
  hit <- candidates[file.exists(file.path(candidates, "core/00_parallel_utils.R"))]
  if (length(hit) == 0L) stop("Cannot locate the mixed-input GP code directory.")
  normalizePath(hit[1L], winslash = "/")
}

mixedgp_source_core <- function(code_dir = mixedgp_code_dir(),
                                include_competitors = TRUE,
                                include_synthetic = TRUE) {
  modules <- c(
    "core/00_parallel_utils.R",
    "core/00_study1_functions.R",
    "core/00_study2_functions.R",
    "core/00_sampler_v030.R",
    "core/00_sampler_v030_api.R"
  )
  if (isTRUE(include_competitors)) {
    modules <- c(modules, "core/03_study2_published_competitors.R", "simulations/competitor_cache.R")
    assign(".mixedgp_competitor_code_dir", normalizePath(code_dir), envir = .GlobalEnv)
    assign(".mixedgp_competitor_cache_root", Sys.getenv("EIVGP_COMPETITOR_CACHE",
      file.path(dirname(normalizePath(code_dir)), "reproduction", "competitor-cache")), envir = .GlobalEnv)
  }
  if (isTRUE(include_synthetic)) {
    modules <- c(modules, "simulations/00_synthetic_data.R")
  }
  modules <- c(modules, "core/00_public_api.R", "core/00_diagnostics.R", "core/00_mcmc_workflow.R")
  modules <- c(modules, "simulations/00_experiment_runner.R")
  for (module in modules) {
    sys.source(file.path(code_dir, module), envir = .GlobalEnv)
  }
  invisible(modules)
}

mixedgp_source_core()

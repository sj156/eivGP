############################################################
## simulation_helpers.R
##
## The single publication-facing helper bundle for both numerical studies.
## Scientific posterior engines remain in their package-ready source modules;
## all simulation design, validation, data freezing, orchestration, manifests,
## task eligibility, aggregation, and fail-closed checks live here.
############################################################

MIXEDGP_SIMULATION_SCHEMA <- "2.1.1"
MIXEDGP_PUBLISHED_METHODS <- c("UC-GP", "LVGP", "EzGP")

## Shared with the installed package; keep a single diagnostic implementation.
if (!exists("mixedgp_summarize_diagnostic_series", mode = "function")) {
  diagnostic_source_file <- tryCatch(sys.frame(1)$ofile, error = function(e) NULL)
  diagnostic_candidates <- c(
    if (!is.null(diagnostic_source_file)) file.path(dirname(diagnostic_source_file), "core/00_diagnostics.R"),
    "core/00_diagnostics.R", file.path("codes", "core/00_diagnostics.R"),
    file.path("revision", "codes", "core/00_diagnostics.R")
  )
  diagnostic_hit <- diagnostic_candidates[file.exists(diagnostic_candidates)]
  if (!length(diagnostic_hit)) stop("Source 00_diagnostics.R before simulation_helpers.R.")
  sys.source(diagnostic_hit[1L], envir = environment())
}

mixedgp_simulation_code_dir <- function(code_dir = NULL) {
  marker <- "simulation_helpers.R"
  if (!is.null(code_dir)) {
    code_dir <- normalizePath(code_dir, winslash = "/", mustWork = TRUE)
    if (!file.exists(file.path(code_dir, marker))) {
      stop("code_dir does not contain ", marker, ": ", code_dir)
    }
    return(code_dir)
  }
  command_file <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  command_dir <- if (length(command_file) > 0L) {
    dirname(sub("^--file=", "", command_file[[1L]]))
  } else {
    character(0)
  }
  candidates <- unique(c(
    command_dir, getwd(), file.path(getwd(), "codes"),
    file.path(getwd(), "revision", "codes")
  ))
  hit <- candidates[file.exists(file.path(candidates, marker))]
  if (length(hit) == 0L) {
    stop("Cannot locate revision/codes; supply code_dir explicitly.")
  }
  normalizePath(hit[[1L]], winslash = "/", mustWork = TRUE)
}

mixedgp_simulation_modules <- function() {
  c(
    "core/00_parallel_utils.R",
    "core/00_study1_functions.R",
    "core/00_study2_functions.R",
    "core/00_sampler_v030.R",
    "core/00_sampler_v030_api.R",
    "core/00_public_api.R",
    "core/00_diagnostics.R",
    "core/00_mcmc_workflow.R",
    "core/03_study2_published_competitors.R",
    "simulations/competitor_cache.R",
    "simulations/00_synthetic_data.R",
    "simulations/04_study1_ablations.R",
    "simulations/04_study2_ablations.R"
  )
}

mixedgp_simulation_engine <- function(code_dir) {
  code_dir <- mixedgp_simulation_code_dir(code_dir)
  old_wd <- setwd(code_dir)
  on.exit(setwd(old_wd), add = TRUE)
  configured_library <- Sys.getenv("MIXEDGP_R_LIBRARY", unset = "")
  if (nzchar(configured_library) && !dir.exists(configured_library)) {
    stop("MIXEDGP_R_LIBRARY does not exist: ", configured_library)
  }
  library_candidates <- c(
    if (nzchar(configured_library)) configured_library else character(0),
    file.path(dirname(code_dir), "R-library")
  )
  library_candidates <- library_candidates[dir.exists(library_candidates)]
  if (length(library_candidates) > 0L) {
    .libPaths(unique(c(
      normalizePath(library_candidates, winslash = "/", mustWork = TRUE),
      .libPaths()
    )))
  }
  engine <- new.env(parent = .GlobalEnv)
  engine$.mixedgp_competitor_code_dir <- code_dir
  engine$.mixedgp_competitor_cache_root <- Sys.getenv("EIVGP_COMPETITOR_CACHE",
    file.path(dirname(code_dir), "reproduction", "competitor-cache"))
  for (module in mixedgp_simulation_modules()) {
    path <- file.path(code_dir, module)
    if (!file.exists(path)) stop("Missing simulation module: ", path)
    sys.source(path, envir = engine, chdir = FALSE)
  }
  engine
}

# Keep this entry point stable; implementation is grouped by responsibility.
.mixedgp_helper_source <- unlist(lapply(sys.frames(), function(frame) {
  get0("ofile", envir = frame, inherits = FALSE, ifnotfound = character(0))
}), use.names = FALSE)
.mixedgp_helper_candidates <- c(
  if (!is.null(.mixedgp_helper_source)) dirname(.mixedgp_helper_source),
  getwd(), file.path(getwd(), "codes"), file.path(getwd(), ".."))
.mixedgp_helper_candidates <- .mixedgp_helper_candidates[
  file.exists(file.path(.mixedgp_helper_candidates, "simulations", "helpers", "design.R"))]
if (!length(.mixedgp_helper_candidates)) stop("Cannot locate simulation helper modules.")
.mixedgp_helper_dir <- normalizePath(.mixedgp_helper_candidates[1L])
for (.mixedgp_helper in c("design.R", "data.R", "execution.R", "gates.R", "summaries.R", "workflow.R")) {
  sys.source(file.path(.mixedgp_helper_dir, "simulations", "helpers", .mixedgp_helper),
             envir = environment())
}
rm(.mixedgp_helper_source, .mixedgp_helper_candidates, .mixedgp_helper_dir, .mixedgp_helper)

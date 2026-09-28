############################################################
## Project-local R setup
##
## Source from revision/codes/. The code never installs packages. When the
## archived project library is present, make it visible to a clean R session;
## sessionInfo() and the competitor preflight still record exact versions.
############################################################

setup_sources <- unlist(lapply(sys.frames(), function(frame) {
  get0("ofile", envir = frame, inherits = FALSE, ifnotfound = character(0))
}), use.names = FALSE)
setup_sources <- setup_sources[file.exists(setup_sources)]
setup_repo <- if (length(setup_sources)) {
  dirname(dirname(dirname(normalizePath(tail(setup_sources, 1L)))))
} else {
  setup_candidates <- c(".", "..", "../..")
  setup_candidates <- setup_candidates[file.exists(file.path(setup_candidates, "codes", "load_mixedgp.R"))]
  if (!length(setup_candidates)) stop("Cannot locate the repository for project setup.")
  normalizePath(setup_candidates[1L])
}
project_library <- normalizePath(
  file.path(setup_repo, "R-library"),
  winslash = "/",
  mustWork = FALSE
)
if (dir.exists(project_library)) {
  .libPaths(unique(c(project_library, .libPaths())))
}

Sys.setenv(
  OMP_NUM_THREADS = "1",
  OPENBLAS_NUM_THREADS = "1",
  MKL_NUM_THREADS = "1",
  VECLIB_MAXIMUM_THREADS = "1",
  NUMEXPR_NUM_THREADS = "1"
)

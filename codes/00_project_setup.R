############################################################
## Project-local R setup
##
## Source from revision/codes/. The code never installs packages. When the
## archived project library is present, make it visible to a clean R session;
## sessionInfo() and the competitor preflight still record exact versions.
############################################################

project_library <- normalizePath(
  file.path("..", "R-library"),
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

## Explicit experiment setup, not a package load/install hook.
## Run from the repository root:
##   Rscript --vanilla experiments/install_eivgp_dependencies.R
## Simulation/reporting commands only check dependencies; they never install.

## Match the experiment engine's library precedence. For an explicit custom
## library, create it here and retain MIXEDGP_R_LIBRARY when running experiments.
configured_library <- Sys.getenv("MIXEDGP_R_LIBRARY", unset = "")
file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
repo <- if (length(file_arg) == 1L) {
  normalizePath(file.path(dirname(sub("^--file=", "", file_arg)), ".."), mustWork = TRUE)
} else getwd()
local_library <- file.path(repo, "R-library")
if (nzchar(configured_library)) {
  dir.create(configured_library, recursive = TRUE, showWarnings = FALSE)
  if (!dir.exists(configured_library)) stop("Cannot create MIXEDGP_R_LIBRARY: ", configured_library)
  .libPaths(unique(c(normalizePath(configured_library),
    if (dir.exists(local_library)) local_library, .libPaths())))
} else if (dir.exists(local_library)) {
  .libPaths(unique(c(local_library, .libPaths())))
}
cat("Experiment R library search path:\n", paste(.libPaths(), collapse = "\n"), "\n")

required <- c(
  "posterior", "TruncatedNormal", # eivGP Imports
  "rmarkdown", "knitr",           # document rendering
  "ggplot2", "dplyr", "tidyr", "patchwork", # tables and figures
  "kergp", "LVGP", "EzGP" # UC-GP, LVGP and EzGP comparisons
)

missing <- required[!vapply(required, requireNamespace, logical(1),
                             quietly = TRUE)]
if (!length(missing)) {
  cat("All required eivGP reproduction dependencies are already installed.\n")
} else {
  cat("Installing required eivGP reproduction dependencies:\n  ",
      paste(missing, collapse = ", "), "\n", sep = "")
  utils::install.packages(missing, lib = .libPaths()[1L], repos = "https://cloud.r-project.org")

  still_missing <- required[!vapply(required, requireNamespace, logical(1),
                                     quietly = TRUE)]
  if (length(still_missing)) {
    stop(
      "Installation did not complete for: ", paste(still_missing, collapse = ", "),
      ". Check the R console output and your network/proxy configuration."
    )
  }
  cat("All required eivGP reproduction dependencies are installed.\n")
}

print(data.frame(package = required,
  version = vapply(required, function(p) as.character(utils::packageVersion(p)), ""),
  library = vapply(required, function(p) dirname(find.package(p)), "")), row.names = FALSE)
cat("\nNext, inspect the standalone competitor plan (no fitting):\n",
    "  Rscript --vanilla experiments/run_competitors.R study1 plan publication\n",
    "Use study2 for Study II. See README for frozen-data and Overleaf setup.\n", sep = "")

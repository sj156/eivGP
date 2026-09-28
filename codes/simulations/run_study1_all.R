## Backward-compatible alias for the audited Study I master.
alias_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
alias_dir <- if (length(alias_arg) > 0L) {
  dirname(normalizePath(sub("^--file=", "", alias_arg[[1L]]), mustWork = TRUE))
}
if (!file.exists(file.path(alias_dir, "load_mixedgp.R"))) alias_dir <- dirname(alias_dir) else {
  getwd()
}
source(file.path(alias_dir, "simulations/run_study1_simulation.R"))

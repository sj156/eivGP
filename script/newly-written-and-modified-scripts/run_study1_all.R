## Backward-compatible alias for the audited Study I master.
alias_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
alias_dir <- if (length(alias_arg) > 0L) {
  dirname(normalizePath(sub("^--file=", "", alias_arg[[1L]]), mustWork = TRUE))
} else {
  getwd()
}
source(file.path(alias_dir, "run_study1_simulation.R"))

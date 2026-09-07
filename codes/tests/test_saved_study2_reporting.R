## Optional integration test: supply an existing Study II run with primary_q2 caches.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) stop("Supply the saved Study II run directory.")
source("codes/00_diagnostics.R")
source("codes/simulation_helpers.R")
source("codes/experiment_reporting.R")
original <- normalizePath(args[1L], mustWork = TRUE)
config <- readRDS(file.path(original,"config","resolved_config.rds"))
stopifnot(config$study == "study2")
config$code_dir <- normalizePath("codes")
config$cells <- Filter(function(x) x$id == "primary_q2", config$cells)
stopifnot(length(config$cells) == 1L)
root <- tempfile("study2-report-test-"); dir.create(root)
source_cell <- file.path(original,"cells","primary_q2")
target_cell <- file.path(root,"cells","primary_q2")
dir.create(target_cell, recursive=TRUE)
stopifnot(file.copy(file.path(source_cell,"results"), target_cell, recursive=TRUE))
dir.create(file.path(root,"config"))
saveRDS(config,file.path(root,"config","resolved_config.rds"))
input_files <- list.files(source_cell, recursive=TRUE, full.names=TRUE)
hash <- tools::md5sum(input_files)
original_engine <- mixedgp_simulation_engine
mixedgp_simulation_engine <- function(code_dir) {
  e <- original_engine(code_dir)
  for (name in ls(e)) if (grepl("^(fit_|run_one_|mixedgp_v030_fit|run_study[12]_published)",name))
    assign(name,function(...) stop("Unexpected fitting during saved-result replay!"),e)
  e$ggsave <- function(...) stop("Unexpected plotting in fitting driver!")
  e
}
engine <- mixedgp_simulation_engine(config$code_dir)
saved_input <- file.path(source_cell, "report_inputs.rds")
if (file.exists(saved_input)) {
  stopifnot(file.copy(saved_input, target_cell))
} else {
  invisible(mixedgp_run_study2_cell(config,config$cells[[1]],engine,root))
}
stopifnot(file.exists(file.path(target_cell,"report_inputs.rds")))
status <- mixedgp_report_run(root, "report", config$code_dir)
stopifnot(all(status$status == "success"), identical(tools::md5sum(input_files),hash))
for (suffix in c("diagnostics", "parameter_diagnostics", "target_diagnostics", "diagnostic_details"))
  stopifnot(file.exists(file.path(root, "reporting", "primary_q2", "tables",
    paste0("study2_mcmc_", suffix, ".csv"))))
## Exercise an empty primary plot with logistic-only saved summary inputs.
path <- file.path(target_cell,"report_inputs.rds")
snapshot <- readRDS(path)
for (name in names(snapshot$context)) {
  x <- snapshot$context[[name]]
  if (is.data.frame(x) && "scenario" %in% names(x)) {
    x$scenario <- rep("logistic_misspec", nrow(x)); snapshot$context[[name]] <- x
  }
}
snapshot$context$STUDY2_SCENARIOS <- "logistic_misspec"
saveRDS(snapshot,path)
status <- mixedgp_report_run(root, "plot", config$code_dir, file.path(root,"logistic-report"))
stopifnot(all(status$status == "success"))
plots <- read.csv(file.path(root,"logistic-report","primary_q2","plot_status.csv"))
stopifnot(any(plots$status == "skipped_empty"))
message("Saved Study II caches and report checkpoints replayed without fitting; empty facets skipped.")

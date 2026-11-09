#!/usr/bin/env Rscript
# Recover deterministic CC predictions/scores from existing artifacts only.
# The source experiment is read-only; use a new output directory.
args <- commandArgs(trailingOnly = TRUE)
if (any(!grepl("^--[^=]+=.+$", args))) stop("Use --name=value arguments")
a <- setNames(as.list(sub("^--[^=]+=", "", args)), sub("^--([^=]+)=.*$", "\\1", args))
allowed <- c("action", "source-root", "data-root", "output", "jobs", "library")
if (anyDuplicated(names(a)) || length(setdiff(names(a), allowed))) stop("Duplicate or unknown option")
get <- function(k, default = NULL) if (is.null(a[[k]])) default else a[[k]]
action <- get("action", "plan")
if (!action %in% c("plan", "run")) stop("action must be plan or run")
for (k in c("source-root", "data-root", "output")) if (is.null(get(k))) stop("Required: --", k)
if (!is.null(get("library"))) .libPaths(c(normalizePath(get("library"), mustWork = TRUE), .libPaths()))
suppressPackageStartupMessages(library(eivGP))
api <- getFromNamespace("eivgp_numerical_api", "eivGP")(system.file(package = "eivGP"))
if (grepl("mixedgp_stage_lapply", paste(deparse(body(api$score_result)), collapse = "\n"), fixed = TRUE))
  stop("Install the updated serial-scoring package before recovery")
source_root <- normalizePath(get("source-root"), mustWork = TRUE)
data_root <- normalizePath(get("data-root"), mustWork = TRUE)
# Resolve through a real parent to reject output symlinks into the source tree.
out_arg <- path.expand(get("output"))
if (!dir.exists(dirname(out_arg))) stop("Output parent must already exist")
output <- if (dir.exists(out_arg)) normalizePath(out_arg, mustWork = TRUE) else
  file.path(normalizePath(dirname(out_arg), mustWork = TRUE), basename(out_arg))
inside <- function(x, parent) identical(x, parent) || startsWith(x, paste0(parent, "/"))
if (inside(output, source_root) || inside(source_root, output) || inside(output, data_root))
  stop("Output must be separate from the source experiment and frozen data")
jobs <- readRDS(file.path(source_root, "job_registry.rds"))
if (anyDuplicated(jobs$job) || any(!grepl("^rep[0-9]+-cal[0-9]+-(UC-GP|LVGP|EzGP)-CC$", jobs$job)))
  stop("Recovery supports calibrated competitor jobs only")
if (!is.null(get("jobs"))) {
  selected <- strsplit(get("jobs"), ",", fixed = TRUE)[[1]]
  if (length(setdiff(selected, jobs$job))) stop("Unknown job selection")
  jobs <- jobs[match(selected, jobs$job), , drop = FALSE]
}
jobs$has_scores <- file.exists(file.path(source_root, jobs$job, "scores.rds"))
jobs$has_predictions <- file.exists(file.path(source_root, jobs$job, "predictions.rds"))
jobs$has_fit <- file.exists(file.path(source_root, jobs$job, "fit.rds"))
shared <- file.path(source_root, sprintf("shared-rep%03d-cal%02d", jobs$rep, jobs$calib))
jobs$recoverable <- jobs$has_predictions | (jobs$has_fit &
  file.exists(file.path(shared, "prospective.rds")))
jobs$recoverable <- jobs$recoverable & file.exists(file.path(shared, "truth.rds"))
print(jobs, row.names = FALSE)
if (action == "plan") quit(status = 0L)
if (is.null(get("jobs"))) stop("Run requires an explicit --jobs selection")
if (any(!jobs$recoverable)) stop("Selected job lacks saved predictions/model, prospective draws, or truth; no fitting is allowed")
protocol_path <- file.path(source_root, "protocol.rds")
protocol <- readRDS(protocol_path)
versions <- vapply(names(protocol$versions), function(p) as.character(packageVersion(p)), character(1))
if (!identical(unname(versions), unname(protocol$versions)) || as.character(getRversion()) != protocol$R)
  stop("R/package versions differ from the source protocol")
dir.create(output, showWarnings = FALSE)
# This recovery is deliberately serial, including native prediction calls.
api$mixedgp_stage_lapply <- function(X, FUN, ...) lapply(X, FUN, ...)
api$stage_time <- function(folder, stage, expr) {
  if (stage == "optimizer") stop("Recovery cannot fit or optimize")
  force(expr)
}
api$atomic_rds <- function(...) stop("Unexpected write from prediction recovery")
hash <- function(p) setNames(unname(tools::md5sum(p)), normalizePath(p, mustWork = TRUE))
for (j in seq_len(nrow(jobs))) {
  job <- jobs[j, ]; src <- file.path(source_root, job$job)
  dest <- file.path(output, job$job)
  if (dir.exists(dest)) stop("Recovery job output already exists: ", dest)
  data_record <- file.path(source_root, sprintf("data-rep%03d.rds", job$rep))
  record <- readRDS(data_record)
  if (basename(record$file) != record$file) stop("Invalid frozen data filename")
  data_path <- file.path(data_root, record$file)
  if (!identical(unname(tools::md5sum(data_path)), record$md5)) stop("Frozen data hash differs")
  artifact <- readRDS(data_path); dat <- artifact$data
  stopifnot(artifact$study == "study2", artifact$scenario == "primary", artifact$design$q == 4)
  cal <- artifact$calibration_sets[[as.character(job$calib)]]
  stopifnot(length(cal) == job$calib)
  shared <- file.path(source_root, sprintf("shared-rep%03d-cal%02d", job$rep, job$calib))
  truth_path <- file.path(shared, "truth.rds")
  pred_path <- file.path(src, "predictions.rds")
  paths <- c(protocol_path, data_record, data_path, truth_path,
    if (file.exists(pred_path)) pred_path else c(file.path(src, "fit.rds"), file.path(shared, "prospective.rds")))
  diagnostics <- file.path(src, "diagnostics.rds")
  old_scores <- file.path(src, "scores.rds")
  paths <- c(paths, c(diagnostics, old_scores)[file.exists(c(diagnostics, old_scores))])
  before <- hash(paths)
  dir.create(dest)
  message("Recovering ", job$job, " from saved artifacts; no fitting")
  started <- proc.time()[["elapsed"]]
  if (file.exists(pred_path)) pred <- readRDS(pred_path) else {
    # Read the model directly from its existing path. cc_evaluate is forbidden
    # to optimize or write anything, even if the cache disappears mid-run.
    pred <- api$cc_evaluate(job$method, dat$train, dat$test, cal,
      readRDS(file.path(shared, "prospective.rds")), protocol$settings,
      protocol$settings$seed + job$rep * 1000L + job$calib + 500000L,
      file.path(src, "fit.rds"))
  }
  saveRDS(pred, file.path(dest, "predictions.rds"))
  message("Predictions saved; scoring serially")
  scores <- api$score_result(pred, dat, cal, readRDS(truth_path), job$method)
  after <- hash(paths)
  if (!identical(before, after)) stop("A source artifact changed during recovery")
  comparison <- if (file.exists(old_scores)) all.equal(scores, readRDS(old_scores), tolerance = 1e-10) else NULL
  saveRDS(scores, file.path(dest, "scores.rds"))
  if (file.exists(diagnostics)) saveRDS(readRDS(diagnostics), file.path(dest, "source-diagnostics.rds"))
  saveRDS(list(source_root = source_root, source_md5 = before, source_protocol = protocol,
    recovery_package_manifest = getFromNamespace("eivgp_module_identity", "eivGP")(),
    comparison = comparison, no_fitting = TRUE, serial = TRUE, time = Sys.time(),
    elapsed_seconds = proc.time()[["elapsed"]] - started), file.path(dest, "recovery-provenance.rds"))
  if (!is.null(comparison) && !isTRUE(comparison)) stop("Recovered scores differ from existing scores")
  message("Recovered ", job$job, " in ", round(proc.time()[["elapsed"]] - started, 2), " seconds")
}

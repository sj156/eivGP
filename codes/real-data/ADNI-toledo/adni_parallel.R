# Fold/chain scheduling. Statistical settings and seed formulas are unchanged.
adni_core_plan <- function(total_cores, n_folds = 3L, n_chains = 4L,
                           max_fold_workers = 3L, os_type = .Platform$OS.type) {
  valid <- function(x) is.numeric(x) && length(x) == 1L && is.finite(x) && x >= 1 && x == floor(x)
  if (!all(vapply(list(total_cores, n_folds, n_chains, max_fold_workers), valid, logical(1))))
    stop("Core, fold, and chain counts must be positive integers.")
  candidates <- seq_len(min(n_folds, max_fold_workers))
  chains <- pmin(n_chains, floor(total_cores / candidates))
  active <- candidates * chains
  # Maximize concurrent chain workers; ties prefer fewer in-memory folds.
  best <- which.max(active)
  f <- candidates[best]; c <- chains[best]
  if (os_type == "windows") { f <- 1L; c <- 1L }
  list(total_cores = as.integer(total_cores), fold_workers = as.integer(f),
       chain_workers = as.integer(c), active_chain_limit = as.integer(f * c),
       n_folds = as.integer(n_folds), n_chains = as.integer(n_chains))
}

adni_print_core_plan <- function(plan) {
  cat("Total core budget:", plan$total_cores, "| concurrent folds:", plan$fold_workers,
      "| chain workers per fold:", plan$chain_workers,
      "| maximum concurrent chains:", plan$active_chain_limit, "\n")
}

adni_run_jobs <- function(jobs, FUN, plan, status_path) {
  worker <- function(i) {
    fold <- jobs$fold[i]; repeat_id <- jobs$repeat_id[i]
    dir <- file.path(jobs$output_root[i], paste0("fold_", fold))
    dir.create(dir, recursive = TRUE, showWarnings = FALSE)
    lock <- file.path(dir, ".fold-lock")
    if (!dir.create(lock, showWarnings = FALSE))
      return(list(ok = FALSE, fold = fold, message = paste("Fold is locked:", lock)))
    on.exit(unlink(lock, recursive = TRUE), add = TRUE)
    writeLines(paste(Sys.info()[["nodename"]], Sys.getpid()), file.path(lock, "owner"))
    log <- file(file.path(dir, "worker.log"), open = "at")
    sink(log); sink(log, type = "message")
    on.exit({ sink(type = "message"); sink(); close(log) }, add = TRUE)
    cat("\nWorker start:", as.character(Sys.time()), "PID", Sys.getpid(), "fold", fold, "\n")
    unlink(file.path(dir, "WORKER_FAILURE.txt"))
    tryCatch({
      value <- FUN(repeat_id, fold)
      cat("Worker finished:", as.character(Sys.time()), "\n")
      list(ok = TRUE, fold = fold, value = value)
    }, error = function(e) {
      text <- paste("Fold", fold, "failed:", conditionMessage(e))
      writeLines(text, file.path(dir, "WORKER_FAILURE.txt"))
      cat(text, "\n")
      list(ok = FALSE, fold = fold, message = text)
    })
  }
  if (plan$fold_workers > 1L) {
    results <- parallel::mclapply(seq_len(nrow(jobs)), worker, mc.cores = plan$fold_workers,
      mc.preschedule = FALSE, mc.set.seed = FALSE, mc.allow.recursive = TRUE)
  } else results <- lapply(seq_len(nrow(jobs)), worker)
  good <- vapply(results, function(x) is.list(x) && isTRUE(x$ok), logical(1))
  status <- data.frame(repeat_id = jobs$repeat_id, fold = jobs$fold, success = good,
    message = vapply(seq_along(results), function(i) {
      x <- results[[i]]
      if (good[i]) "finished" else if (is.list(x) && !is.null(x$message)) x$message else
        "Worker exited without a result; inspect worker.log and any stale fold lock."
    }, character(1)))
  write.csv(status, status_path, row.names = FALSE)
  if (any(!good)) stop("One or more folds failed. Completed fold checkpoints are retained. See ",
                       status_path, ". No combined report generated.")
  lapply(results, `[[`, "value")
}

# Direct Rmd rendering keeps a single-repeat entry point.
adni_run_folds <- function(folds, FUN, plan, output_root) {
  jobs <- data.frame(repeat_id = 1L, fold = folds, output_root = output_root)
  values <- adni_run_jobs(jobs, function(repeat_id, fold) FUN(fold), plan,
                         file.path(output_root, "fold_worker_status.csv"))
  setNames(values, as.character(folds))
}

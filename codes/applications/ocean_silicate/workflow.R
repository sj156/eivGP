# GLODAP silicate/oxygen workflow: public frozen inputs, local reproducible outputs.
ocean_read_inputs <- function(data_dir) {
  cohort <- read.csv(file.path(data_dir, "cohort_manifest.csv"))
  stopifnot(nrow(cohort) == 200L, !anyDuplicated(cohort$profile))
  for (fold in 1:5) {
    train <- read.csv(file.path(data_dir, sprintf("fold%02d_train.csv", fold)))
    test <- read.csv(file.path(data_dir, sprintf("fold%02d_test.csv", fold)))
    stopifnot(nrow(train) == 160L, nrow(test) == 40L,
      !length(intersect(train$profile, test$profile)),
      sum(train$calib_o20) == 32L, sum(train$calib_o50) == 80L,
      all(train$calib_o20 <= train$calib_o50),
      all(test$U < 180), all(train$U < 180))
  }
  cohort$oxygen_class <- factor(findInterval(cohort$U, c(50, 100, 150)) + 1L,
    levels = 1:4, labels = c("<50", "50–100", "100–150", "150–180"))
  cohort
}

ocean_run_analysis <- function(repository_dir, data_dir, output_dir, profile = c("quick", "paper")) {
  profile <- match.arg(profile)
  source_dir <- file.path(repository_dir, "codes", "applications", "ocean_silicate")
  files <- list.files(data_dir, "[.]csv$", full.names = TRUE)
  hashes <- tools::md5sum(files); names(hashes) <- basename(files)
  source_hashes <- tools::md5sum(file.path(source_dir, c("fit.R", "continue.R", "workflow.R")))
  names(source_hashes) <- basename(names(source_hashes))
  notebook_cache_identity(output_dir, list(profile = profile, data = hashes, code = source_hashes,
    eivGP = as.character(packageVersion("eivGP"))))
  dir.create(file.path(output_dir, "data"), showWarnings = FALSE)
  for (file in files) {
    target <- file.path(output_dir, "data", basename(file))
    if (!file.exists(target)) stopifnot(file.copy(file, target))
    stopifnot(identical(unname(tools::md5sum(file)), unname(tools::md5sum(target))))
  }
  env <- new.env(parent = environment()); env$OCEAN_CASE_ROOT <- output_dir
  run_dir <- file.path(output_dir, if (profile == "quick") "smoke_class4" else "full5_class4")
  notebook_with_env(c(SI_O2_CLASS4_SMOKE = if (profile == "quick") "1" else "0"), {
    if (!file.exists(file.path(run_dir, "complete.flag"))) sys.source(file.path(source_dir, "fit.R"), env)
  })
  if (!file.exists(file.path(run_dir, "complete.flag"))) stop("Ocean fitting is incomplete; inspect task_status.csv.")
  if (profile == "paper" && !file.exists(file.path(output_dir, "continuation_plus20k", "complete.flag")))
    sys.source(file.path(source_dir, "continue.R"), env)
  invisible(run_dir)
}

ocean_read_results <- function(output_dir, profile) {
  initial <- file.path(output_dir, if (profile == "quick") "smoke_class4" else "full5_class4")
  metrics <- predictions <- diagnostics <- list()
  for (fold in if (profile == "quick") 1L else 1:5) {
    for (cell in if (profile == "quick") c("O0", "O50") else c("O0", "O20", "O50")) {
      continued <- file.path(output_dir, "continuation_plus20k", sprintf("fold%02d", fold), cell)
      path <- if (profile == "paper" && file.exists(file.path(continued, "predictive_metrics.csv"))) continued else
        file.path(initial, sprintf("fold%02d", fold), cell)
      key <- paste(fold, cell)
      metrics[[key]] <- read.csv(file.path(path, "predictive_metrics.csv"))
      z <- read.csv(file.path(path, "predictions.csv")); z$fold <- fold; z$cell <- cell
      predictions[[key]] <- z
      z <- read.csv(file.path(path, "convergence_summary.csv")); z$fold <- fold; z$cell <- cell
      diagnostics[[key]] <- z
    }
  }
  list(metrics = do.call(rbind, metrics), predictions = do.call(rbind, predictions),
       diagnostics = do.call(rbind, diagnostics))
}

# Reader workflows use the existing scientific implementations.
adni_run_analysis <- function(repository_dir, profile = c("quick", "paper"), cores = 4L,
    data_dir = file.path(repository_dir, "real-data"),
    output_dir = file.path(repository_dir, "real_data_application_outputs", "adni")) {
  profile <- match.arg(profile)
  application_dir <- file.path(repository_dir, "codes", "applications")
  data_hash <- unname(tools::md5sum(file.path(data_dir, "toledo_adni_cohort_n495.csv")))
  code_files <- list.files(application_dir, "^ADNI.*[.]R$", full.names = TRUE)
  code_hash <- tools::md5sum(code_files); names(code_hash) <- basename(code_files)
  notebook_cache_identity(file.path(output_dir, if (profile == "quick") "smoke" else "five_fold"),
    list(profile = profile, data = data_hash, code = code_hash, eivGP = as.character(packageVersion("eivGP"))))
  result_dir <- file.path(output_dir, if (profile == "quick") "smoke" else "five_fold")
  result_file <- file.path(result_dir, "main", "table1_prediction.csv")
  folds <- if (profile == "quick") 1L else 1:5
  completed <- file.path(result_dir, paste0("fold_", folds), "CASE_STUDY_COMPLETED.txt")
  if (file.exists(result_file) && all(file.exists(completed)))
    return(list(output_dir = result_dir, folds = folds, metrics = read.csv(result_file)))
  env <- new.env(parent = environment())
  on.exit({
    if (isTRUE(get0("LOCK_OWNED", env, inherits = FALSE, ifnotfound = FALSE)))
      unlink(env$RUN_LOCK, recursive = TRUE)
  }, add = TRUE)
  notebook_with_env(c(ADNI_PROJECT_DIR = application_dir, ADNI_DATA_DIR = data_dir,
    ADNI_OUTPUT_DIR = output_dir, EIVGP_SMOKE_TEST = if (profile == "quick") "1" else "0",
    ADNI_SMOKE_ALL_FOLDS = "0", EIVGP_CORES = as.character(cores)), {
    sys.source(file.path(application_dir, "ADNI_analysis.R"), envir = env)
  })
  list(output_dir = env$OUTPUT_ROOT, folds = env$FOLDS_TO_RUN,
    metrics = read.csv(file.path(env$OUTPUT_ROOT, "main", "table1_prediction.csv")))
}

prepare_adni_inputs <- function(repository_dir, raw_dir, output_dir) {
  env <- new.env(parent = environment())
  env$repository_dir <- repository_dir
  notebook_with_env(c(ADNI_RAW_DIR = raw_dir, ADNI_TOLEDO_DATA_DIR = output_dir), {
    sys.source(file.path(repository_dir, "codes", "applications", "preprocessing", "ADNI_preprocess.R"), env)
  })
  invisible(output_dir)
}

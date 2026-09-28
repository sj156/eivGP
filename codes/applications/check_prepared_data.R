check_prepared_data <- function(application_dir) {
  source(file.path(application_dir, "real_data_paths.R"), local = TRUE)
  paths <- application_paths(application_dir)
  adni_file <- file.path(paths$adni_data, "toledo_adni_cohort_n495.csv")
  if (file.exists(adni_file)) {
    source(file.path(application_dir, "ADNI_cv.R"), local = TRUE)
    dat <- adni_read_cohort(paths$adni_data)
    folds <- adni_make_folds(dat)
    cat("ADNI: verified", nrow(dat), "participants and", length(unique(folds$fold)), "folds.\n")
  } else cat("ADNI inputs unavailable; approved-access data are required.\n")
  source(file.path(application_dir, "ocean_silicate", "workflow.R"), local = TRUE)
  dat <- ocean_read_inputs(paths$ocean_data)
  cat("Ocean: verified", nrow(dat), "profiles and five frozen train/test splits.\n")
  invisible(TRUE)
}

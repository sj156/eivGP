# All paths are anchored to codes/, never the caller's working directory.
application_paths <- function(code_dir) {
  repo <- normalizePath(file.path(code_dir, ".."), mustWork = TRUE)
  list(repo = repo,
    adni_data = Sys.getenv("ADNI_DATA_DIR", file.path(repo, "real-data")),
    ocean_data = Sys.getenv("OCEAN_DATA_DIR", file.path(repo, "real-data", "ocean")),
    adni_output = Sys.getenv("ADNI_OUTPUT_DIR", file.path(repo, "real_data_application_outputs", "adni")),
    ocean_output = Sys.getenv("OCEAN_OUTPUT_DIR", file.path(repo, "real_data_application_outputs", "ocean")))
}
resolve_adni_data_dir <- function(project_dir) {
  path <- application_paths(project_dir)$adni_data
  if (!file.exists(file.path(path, "toledo_adni_cohort_n495.csv")))
    stop("Missing cleaned ADNI cohort: ", file.path(path, "toledo_adni_cohort_n495.csv"),
         ". Put the collaborator's cleaned cohort in real-data/ or set ADNI_DATA_DIR.")
  normalizePath(path, mustWork = TRUE)
}

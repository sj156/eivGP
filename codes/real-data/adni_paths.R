# Resolve private cleaned data identically for Rscript and direct Rmd rendering.
resolve_adni_data_dir <- function(project_dir) {
  required <- c("toledo_adni_cohort_n495.csv", "toledo_adni_balanced_repeated_3fold.csv")
  explicit <- Sys.getenv("ADNI_DATA_DIR", "")
  if (nzchar(explicit)) {
    missing <- required[!file.exists(file.path(explicit, required))]
    if (length(missing)) stop("ADNI_DATA_DIR=", explicit, " is missing: ",
      paste(missing, collapse = ", "), ". Correct this path or supply both frozen CSVs.")
    return(normalizePath(explicit, mustWork = TRUE))
  }
  repo <- normalizePath(file.path(project_dir, "..", "..", ".."), mustWork = TRUE)
  candidates <- c(file.path(repo, "real-data", "adni"), file.path(repo, "real-data"),
                  file.path(project_dir, "data"))
  for (path in candidates) {
    if (all(file.exists(file.path(path, required)))) return(normalizePath(path, mustWork = TRUE))
  }
  stop("ADNI cleaned data were not found. Put both ", paste(required, collapse = " and "),
       " in ", file.path(repo, "real-data"), " or set ADNI_DATA_DIR to their directory. Searched: ",
       paste(candidates, collapse = "; "))
}

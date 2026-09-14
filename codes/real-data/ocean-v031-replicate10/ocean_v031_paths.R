## Shared path helpers for the ocean 0.3.1 replicate-10 bundle.
## Cohort files live in codes/real-data/ocean_data/ (sibling of this folder).
## Override with OCEAN_DATA_DIR if needed.

ocean_v031_bundle_dir <- function() {
  file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(file_arg) >= 1L) {
    return(dirname(normalizePath(sub("^--file=", "", file_arg[[1L]]), mustWork = TRUE)))
  }
  normalizePath(getwd(), mustWork = TRUE)
}

ocean_v031_data_dir <- function(bundle_dir = ocean_v031_bundle_dir(), explicit = NULL) {
  env <- Sys.getenv("OCEAN_DATA_DIR", unset = "")
  candidates <- c(
    explicit,
    if (nzchar(env)) env,
    file.path(bundle_dir, "..", "ocean_data")
  )
  candidates <- unique(candidates[nzchar(candidates)])
  for (p in candidates) {
    if (file.exists(file.path(p, "study1_data.csv")) &&
        file.exists(file.path(p, "meta.json"))) {
      return(normalizePath(p, mustWork = TRUE))
    }
  }
  stop(
    "Could not find ocean_data/study1_data.csv. ",
    "Keep this folder next to codes/real-data/ocean_data, or set OCEAN_DATA_DIR."
  )
}

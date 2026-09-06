############################################################
## Generate and freeze publication synthetic datasets
##
## Run from revision/codes. Existing files are verified and retained unless
## SYNTHETIC_OVERWRITE is explicitly set to TRUE.
############################################################

source("00_project_setup.R")
source("load_mixedgp.R")

if (!exists("SYNTHETIC_CONFIG")) SYNTHETIC_CONFIG <- "thorough"
if (!SYNTHETIC_CONFIG %in% c("quick", "balanced", "thorough")) {
  stop("SYNTHETIC_CONFIG must be quick, balanced, or thorough.")
}
if (!exists("SYNTHETIC_OVERWRITE")) SYNTHETIC_OVERWRITE <- FALSE
if (!exists("SYNTHETIC_N_CORES")) SYNTHETIC_N_CORES <- NULL
if (!exists("SYNTHETIC_OUT_DIR")) {
  SYNTHETIC_OUT_DIR <- file.path("..", "data-synthetic")
}

config <- switch(
  SYNTHETIC_CONFIG,
  quick = list(n_rep = 2L, n_test_study1 = 150L, n_test_study2 = 150L),
  balanced = list(n_rep = 10L, n_test_study1 = 300L, n_test_study2 = 400L),
  thorough = list(n_rep = 50L, n_test_study1 = 500L, n_test_study2 = 500L)
)

study1_manifest <- generate_study1_synthetic_datasets(
  n_rep = config$n_rep,
  directory = file.path(SYNTHETIC_OUT_DIR, "study1"),
  n_cores = SYNTHETIC_N_CORES,
  overwrite = SYNTHETIC_OVERWRITE,
  n = 100L,
  n_test = config$n_test_study1,
  m = 6L,
  scenario = "active",
  threshold_design = "imbalanced",
  calib_grid = c(0L, 5L, 10L, 20L, 50L)
)

study2_manifest <- generate_study2_synthetic_datasets(
  n_rep = config$n_rep,
  scenarios = c(
    "primary", "latent_additive_control", "high_uncertainty",
    "logistic_misspec"
  ),
  directory = file.path(SYNTHETIC_OUT_DIR, "study2"),
  n_cores = SYNTHETIC_N_CORES,
  overwrite = SYNTHETIC_OVERWRITE,
  n = 120L,
  n_test = config$n_test_study2,
  q = 4L,
  m = 4L,
  calib_grid = c(0L, 10L, 25L, 50L, 80L)
)

combined_manifest <- rbind(study1_manifest, study2_manifest)
dir.create(SYNTHETIC_OUT_DIR, recursive = TRUE, showWarnings = FALSE)
combined_csv <- file.path(SYNTHETIC_OUT_DIR, "manifest.csv")
combined_rds <- file.path(SYNTHETIC_OUT_DIR, "manifest.rds")
combined_csv_tmp <- tempfile(
  "manifest-", tmpdir = SYNTHETIC_OUT_DIR, fileext = ".csv"
)
combined_rds_tmp <- tempfile(
  "manifest-", tmpdir = SYNTHETIC_OUT_DIR, fileext = ".rds"
)
on.exit({
  if (file.exists(combined_csv_tmp)) unlink(combined_csv_tmp)
  if (file.exists(combined_rds_tmp)) unlink(combined_rds_tmp)
}, add = TRUE)
utils::write.csv(
  combined_manifest,
  combined_csv_tmp,
  row.names = FALSE
)
saveRDS(
  list(
    schema_version = MIXEDGP_DATA_SCHEMA_VERSION,
    config = SYNTHETIC_CONFIG,
    manifest = combined_manifest,
    session = utils::sessionInfo()
  ),
  combined_rds_tmp,
  version = 3L
)
mixedgp_atomic_replace(combined_rds_tmp, combined_rds)
mixedgp_atomic_replace(combined_csv_tmp, combined_csv)

print(aggregate(file ~ study + scenario, combined_manifest, length))

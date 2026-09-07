## Configuration/persistence tests only; no publication fits.
source("codes/simulation_helpers.R")
engine <- mixedgp_simulation_engine(normalizePath("codes"))
for (cores in c(1L, 2L)) {
  root <- tempfile("dataset-restart-"); dir.create(root)
  run <- function(fail) suppressWarnings(mixedgp_run_replications(
    as.list(1:3), function(i) {
      path <- file.path(root, paste0(i, ".rds"))
      if (file.exists(path)) return(readRDS(path))
      if (fail && i == 2L) stop("interrupted dataset")
      result <- list(id = i, value = runif(1))
      mixedgp_save_replication(result, path)
      result
    }, n_cores = cores, seeds = 101:103,
    status_path = file.path(root, "status.csv"), study = "restart test",
    parallel_map = engine$mixedgp_parallel_lapply))
  first <- run(TRUE)
  saved <- tools::md5sum(file.path(root, c("1.rds", "3.rds")))
  statuses <- list.files(file.path(root, "status.csv.tasks"), full.names = TRUE)
  stopifnot(length(first) == 2L, length(statuses) == 3L,
    identical(unname(vapply(statuses, function(p) read.csv(p)$status, "")),
              c("success", "failed", "success")))
  second <- run(FALSE)
  stopifnot(length(second) == 3L,
    identical(saved, tools::md5sum(names(saved))),
    all(read.csv(file.path(root, "status.csv"))$status == "success"))
}
prior <- engine$mixedgp_v030_priors(p = 1L, d = 1L, m_vec = 3L)
stopifnot(identical(prior$signal_shape, c(13, 3)),
  abs(diff(pbeta(c(.65, .95), 13, 3)) - .90206567) < 1e-7)
cat("Dataset restart, live bookkeeping, cache preservation, and Beta(13,3) passed.\n")

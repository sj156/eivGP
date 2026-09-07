#!/usr/bin/env Rscript
## Repository-only case study; not part of the reusable package API.
## Rscript experiments/case_study_mixing.R [additional_iterations=2500] [cores=4] [output_dir]
## Use --smoke for a cheap checkpoint-continuation test, not a scientific run.

case_continue <- function(fit, additional, cores) {
  cp <- fit$checkpoint
  stopifnot(identical(fit$sampler_version, "0.3.1"), identical(cp$version, 2L),
    identical(cp$sampler_version, "0.3.1"), identical(cp$control, fit$control),
    cp$thin == 1L, length(cp$states) == length(fit$mcmc$samples_by_chain$u))
  old <- fit$mcmc$samples_by_chain
  expected <- cp$iteration - cp$burn
  stopifnot(all(vapply(old$u, nrow, integer(1)) == expected))
  cp$old_chains <- lapply(seq_along(cp$states), function(i) {
    z <- lapply(old, function(chains) chains[[i]])
    names(z) <- paste0("samples_", names(z))
    z$stats <- fit$mcmc$chain_stats[i, , drop = FALSE]
    if (!is.null(fit$mcmc$chain_initial)) z$initial_state <- fit$mcmc$chain_initial[[i]]
    z
  })
  priors <- fit$priors
  priors$dictionary_type <- NULL # derived descriptor, not a prior input
  d <- fit$data
  ans <- fit_eivgp_1d_v030(x_raw = d$x_raw, y_raw = d$y_raw, c_ord = d$c_ord,
    u_obs = d$u_obs, calib_idx = d$calib_idx, m = d$m,
    u_true = d$u_true, tau_true = d$tau_true,
    n_iter = cp$iteration + additional, burn = cp$burn, thin = 1L,
    n_chains = length(cp$states), parallel_chains = cores > 1L, n_cores = cores,
    kernel = fit$kernel$name, matern_nu = fit$kernel$matern_nu,
    priors = priors, sampler_control = fit$control, .resume = cp)
  stopifnot(identical(ans$priors, fit$priors), identical(ans$data, fit$data))
  for (nm in names(old)) for (i in seq_along(old[[nm]])) {
    x <- ans$mcmc$samples_by_chain[[nm]][[i]]
    prefix <- if (is.null(dim(x))) x[seq_len(expected)] else x[seq_len(expected), , drop = FALSE]
    stopifnot(NROW(x) == expected + additional,
      isTRUE(all.equal(prefix, old[[nm]][[i]],
        check.attributes = FALSE)))
  }
  ans
}

case_series <- function(fit, panel, n_latent) {
  raw <- mixedgp_v030_raw_series(fit)
  targets <- mixedgp_study1_target_series(fit, X = panel$X, C = panel$C,
    U = panel$U, max_draws_per_chain = Inf, n_latent = n_latent, seed = 481517L)
  J <- fit$mcmc$samples_by_chain$J
  k <- length(fit$priors$u_dictionary[[1L]])
  for (j in seq_len(k)) raw[[paste0("dictionary_state_", j)]] <-
    lapply(J, function(x) as.numeric(x[, 1L] == j))
  list(raw = raw, target = targets)
}

case_report <- function(fit, panel, out, label, n_latent) {
  series <- case_series(fit, panel, n_latent)
  detail <- do.call(rbind, lapply(names(series), function(role) {
    x <- mixedgp_summarize_diagnostic_series(series[[role]])
    x$diagnostic_role <- role
    x
  }))
  mixedgp_write_diagnostic_tables(fit$diagnostics$summary, detail, out, label)
  all <- c(series$raw, series$target)
  chain_summary <- do.call(rbind, lapply(names(all), function(nm) {
    do.call(rbind, lapply(seq_along(all[[nm]]), function(i) {
      z <- all[[nm]][[i]]
      data.frame(parameter = nm, chain = i, draws = length(z), mean = mean(z),
        sd = sd(z), q025 = unname(quantile(z, .025)), q975 = unname(quantile(z, .975)))
    }))
  }))
  write.csv(chain_summary, file.path(out, paste0(label, "_chain_summary.csv")), row.names = FALSE)
  ## Chain means of dictionary indicators are occupancy probabilities; retain
  ## unvisited states and undefined diagnostics rather than silently dropping them.
  write.csv(chain_summary[grepl("^dictionary_state_", chain_summary$parameter), ],
    file.path(out, paste0(label, "_dictionary_occupancy.csv")), row.names = FALSE)
  selected <- names(all)[!grepl("^[uU]\\[|^dictionary_state_", names(all))]
  grDevices::pdf(file.path(out, paste0(label, "_traces.pdf")), width = 10, height = 6)
  tryCatch(for (nm in selected) {
    matplot(do.call(cbind, all[[nm]]), type = "l", lty = 1, col = seq_along(all[[nm]]),
      xlab = "Retained iteration (no thinning)", ylab = nm, main = paste(label, nm))
    legend("topright", legend = paste("Chain", seq_along(all[[nm]])),
      col = seq_along(all[[nm]]), lty = 1, bty = "n")
  }, finally = grDevices::dev.off())
  invisible(detail)
}

case_main <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  smoke <- identical(args, "--smoke")
  positive <- function(x) {
    y <- suppressWarnings(as.numeric(x))
    if (length(y) != 1L || !is.finite(y) || y < 1 || y != floor(y)) stop("Expected positive integer.")
    as.integer(y)
  }
  additional <- if (smoke) 2L else if (length(args)) positive(args[1L]) else 2500L
  cores <- if (smoke) 1L else if (length(args) >= 2L) positive(args[2L]) else 4L
  cli <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  script <- normalizePath(sub("^--file=", "", cli[1L]))
  repo <- Sys.getenv("EIVGP_REPO", unset = dirname(dirname(script)))
  if (!file.exists(file.path(repo, "codes/load_mixedgp.R"))) stop("Set EIVGP_REPO to your repository root.")
  oldwd <- getwd(); on.exit(setwd(oldwd), add = TRUE)
  setwd(repo); source("codes/load_mixedgp.R")
  input <- file.path(repo, "reproduction/development/results/study1/study1-development-ef651d075f78/cells/eta1_balanced/results/study1_publication/mcmc_failure_rep001_cal020_103801fa52f9.rds")
  hash <- tools::md5sum(input)
  fit <- readRDS(input)$fit
  stopifnot(length(fit$data$calib_idx) == 20L, fit$sampler_version == "0.3.1")
  if (smoke) {
    two <- case_continue(fit, 2L, 1L)
    split <- case_continue(case_continue(fit, 1L, 1L), 1L, 1L)
    stopifnot(identical(two$mcmc$samples_by_chain, split$mcmc$samples_by_chain),
      identical(tools::md5sum(input), hash))
    message("Smoke passed: two steps equal one plus one; original draws and input unchanged.")
    return(invisible(NULL))
  }
  out <- if (length(args) >= 3L) args[3L] else file.path(repo, "reproduction/case-studies",
    paste0("eta1-rep1-cal20-", format(Sys.time(), "%Y%m%d-%H%M%S")))
  if (dir.exists(out) || file.exists(out)) stop("Use a new output directory; existing results are never overwritten.")
  dir.create(out, recursive = TRUE)
  ## Fixed, calibration-anchored panel: training locations, not held-out accuracy.
  ids <- fit$data$calib_idx
  ids <- ids[unique(round(seq(1, length(ids), length.out = min(5L, length(ids)))))]
  panel <- list(X = fit$data$x_raw[ids, , drop = FALSE], C = fit$data$c_ord[ids],
    U = matrix(fit$data$u_obs[ids], ncol = 1L), training_rows = ids)
  n_latent <- 32L
  saveRDS(list(input = input, input_md5 = hash, panel = panel, n_latent = n_latent,
    integration_seed = 481517L, additional = additional, cores = cores,
    session = sessionInfo(), source_md5 = tools::md5sum(list.files("codes", pattern = "[.]R$", full.names = TRUE))),
    file.path(out, "provenance.rds"))
  message("Computing baseline diagnostics...")
  before <- case_report(fit, panel, out, "before", n_latent)
  elapsed <- system.time(updated <- case_continue(fit, additional, cores))
  saveRDS(updated, file.path(out, "extended_fit.rds")) # save before costly reporting
  saveRDS(elapsed, file.path(out, "extension_time.rds"))
  message("Computing diagnostics on all old and new draws...")
  after <- case_report(updated, panel, out, "after", n_latent)
  comparison <- merge(before, after, by = c("parameter", "diagnostic_role"),
    suffixes = c("_before", "_after"), all = TRUE)
  write.csv(comparison, file.path(out, "before_after.csv"), row.names = FALSE)
  stopifnot(identical(tools::md5sum(input), hash))
  message("Finished bounded case study: ", out,
    "\nWarnings are retained. No automatic continuation or convergence claim.")
}

if (sys.nframe() == 0L) case_main()

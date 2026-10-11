# 04: U imputation, category means, continuous functions, and prediction diagnostics.
# Loading defines functions only. Direct Rscript execution uses a parameter RDS.

# Final paper-profile diagnostics. Source this file after 01_setup.R; sourcing
# only defines functions. It never fits models, changes masks, or installs packages.

# 1. Frozen inputs and final-fit provenance ------------------------------------
ocean_diag_require <- function(ok, message) {
  if (!isTRUE(ok)) stop(message, call. = FALSE)
}

ocean_diag_csv <- function(path, columns) {
  ocean_diag_require(file.exists(path), paste("Missing input:", path))
  x <- utils::read.csv(path, stringsAsFactors = FALSE)
  ocean_diag_require(all(columns %in% names(x)), paste("Unexpected columns:", path))
  x
}

ocean_diag_flag <- function(x) {
  ocean_diag_require(!anyNA(x) && all(x %in% c(0, 1, FALSE, TRUE)),
    "Observation masks must contain only 0/1 or FALSE/TRUE.")
  as.logical(x)
}

ocean_diag_inputs <- function(data_dir) {
  columns <- c("cohort_id", "profile", "G2salinity", "G2pressure", "log_pressure", "U", "Y")
  cohort <- ocean_diag_csv(file.path(data_dir, "cohort_manifest.csv"), columns)
  ocean_diag_require(nrow(cohort) == 200L && !anyDuplicated(cohort$cohort_id) &&
    !anyDuplicated(cohort$profile), "Expected the frozen 200-profile cohort.")
  check_values <- function(d) {
    ocean_diag_require(all(is.finite(as.matrix(d[, columns[3:7]]))) &&
      all(d$U > 0 & d$U < 180) && all(d$G2pressure >= 300 & d$G2pressure <= 1500) &&
      max(abs(d$log_pressure - log(d$G2pressure))) < 1e-10,
      "Inputs must use physical units and natural-log pressure in the frozen domain.")
    d$C <- findInterval(d$U, c(50, 100, 150)) + 1L # Derive the four oxygen classes from measured oxygen.
    d
  }
  cohort <- check_values(cohort)
  ocean_diag_require(identical(as.integer(table(factor(cohort$C, levels = 1:4))), rep(50L, 4)),
    "Expected 50 reference observations in each recorded oxygen band.")
  train <- test <- audit <- vector("list", 5L)
  for (fold in 1:5) {
    tr <- check_values(ocean_diag_csv(file.path(data_dir, sprintf("fold%02d_train.csv", fold)),
      c(columns, "calib_o20", "calib_o50")))
    te <- check_values(ocean_diag_csv(file.path(data_dir, sprintf("fold%02d_test.csv", fold)), columns))
    tr$calib_o20 <- ocean_diag_flag(tr$calib_o20)
    tr$calib_o50 <- ocean_diag_flag(tr$calib_o50)
    ocean_diag_require(nrow(tr) == 160L && nrow(te) == 40L &&
      !anyDuplicated(c(tr$cohort_id, te$cohort_id)) &&
      setequal(c(tr$cohort_id, te$cohort_id), cohort$cohort_id) &&
      !length(intersect(tr$profile, te$profile)), "Train/test separation is inconsistent.")
    for (d in list(tr, te)) {
      ref <- cohort[match(d$cohort_id, cohort$cohort_id), columns]
      ocean_diag_require(isTRUE(all.equal(d[, columns], ref, check.attributes = FALSE)),
        "Fold data differ from the frozen cohort.")
    }
    ocean_diag_require(all(tr$calib_o20 <= tr$calib_o50) &&
      all(table(factor(tr$C[tr$calib_o20], levels = 1:4)) == 8L) &&
      all(table(factor(tr$C[tr$calib_o50], levels = 1:4)) == 20L) &&
      all(table(factor(tr$C, levels = 1:4)) == 40L) &&
      all(table(factor(te$C, levels = 1:4)) == 10L), "Calibration or class counts changed.")
    train[[fold]] <- tr; test[[fold]] <- te
    audit[[fold]] <- data.frame(fold = fold, n_train = nrow(tr), n_test = nrow(te),
      n_calib_O20 = sum(tr$calib_o20), n_calib_O50 = sum(tr$calib_o50),
      n_common_missing = sum(!tr$calib_o50))
  }
  ocean_diag_require(!anyDuplicated(unlist(lapply(test, `[[`, "cohort_id"))),
    "Each cohort observation must be held out exactly once across five folds.")
  list(cohort = cohort, train = train, test = test, audit = do.call(rbind, audit))
}

ocean_diag_sources <- function(output_dir) {
  selected<-ocean_final_eiv_sources(output_dir,"paper")
  ocean_diag_require(all(selected$available),"Final fits selected from current initial diagnostics are incomplete.")
  sources<-selected[,c("fold","cell","source","continued","fit_file")]
  for(i in seq_len(nrow(sources))) {
    s<-sources[i,];path<-file.path(output_dir,s$source)
    for(file in c(s$fit_file,"predictive_metrics.csv","prediction_distributions.rds",
      "U_posterior_intervals.csv","convergence_summary.csv"))
      ocean_diag_require(file.exists(file.path(path,file)),paste("Incomplete paper output:",path,file))
  }
  manifest <- file.path(output_dir, "continuation_plus20k", "final_review_20260927", "fit_sources.csv")
  if (file.exists(manifest)) {
    recorded <- ocean_diag_csv(manifest, c("fold", "cell", "source", "continued"))
    i <- match(paste(sources$fold, sources$cell), paste(recorded$fold, recorded$cell))
    ocean_diag_require(nrow(recorded) == 15L && !anyNA(i) &&
      identical(sources$source, recorded$source[i]), "Final-review source manifest disagrees with the paper run.")
  }
  sources
}

# 2. Training U imputation: five metrics on the same 80 hidden rows ------------
ocean_diag_five_metrics <- function(d) {
  error <- d$posterior_mean - d$reference_U
  data.frame(n = nrow(d), Bias = mean(error), RMSE = sqrt(mean(error^2)),
    MAE = mean(abs(error)), Coverage95 = mean(d$reference_U >= d$q025 & d$reference_U <= d$q975),
    Width95 = mean(d$q975 - d$q025))
}

ocean_diag_fold_summary <- function(d, groups, metrics) {
  # Fold SD is descriptive variability across overlapping CV training sets.
  # It is neither posterior SD nor an independent-replicate Monte Carlo SE.
  key <- do.call(paste, c(d[groups], sep = "/"))
  do.call(rbind, lapply(split(seq_len(nrow(d)), key), function(i) {
    z <- d[i, , drop = FALSE]
    ocean_diag_require(nrow(z) == 5L && setequal(z$fold, 1:5), "Expected one score per fold.")
    do.call(rbind, lapply(metrics, function(metric) {
      ocean_diag_require(all(is.finite(z[[metric]])), paste("Non-finite score:", metric))
      data.frame(z[1, groups, drop = FALSE], metric = metric, mean = mean(z[[metric]]),
        sd = stats::sd(z[[metric]]), n_folds = nrow(z), row.names = NULL)
    }))
  }))
}

# Display all 160 training rows in each panel. This is deliberately separate
# from the common-hidden scoring set below: calibration values are NOT imputations.
ocean_diag_u_display <- function(output_dir, inputs, sources) {
  parts <- scales <- counts <- list()
  for (i in seq_len(nrow(sources))) {
    s <- sources[i, ]; tr <- inputs$train[[s$fold]]
    observed <- switch(s$cell, O0 = rep(FALSE, nrow(tr)),
      O20 = tr$calib_o20, O50 = tr$calib_o50)
    z <- ocean_diag_csv(file.path(output_dir, s$source, "U_posterior_intervals.csv"),
      c("cohort_id", "C", "posterior_mean", "q025", "q975", "scale"))
    expected_scale <- if (s$cell == "O0") "model" else "raw"
    ocean_diag_require(!anyDuplicated(z$cohort_id) &&
      setequal(z$cohort_id, tr$cohort_id[!observed]) && all(z$scale == expected_scale),
      "U display rows or scales disagree with the fitted observation mask.")
    # The full TRAINING reference mean/SD are an external plotting convention.
    # Hidden U is used here only for evaluation/axes, never to fit/align O0.
    center <- mean(tr$U); scale <- stats::sd(tr$U)
    ocean_diag_require(is.finite(scale) && scale > 0, "Invalid training reference SD.")
    d <- data.frame(fold = s$fold, cell = s$cell, cohort_id = tr$cohort_id,
      C = tr$C, observed = observed, common_hidden = !tr$calib_o50,
      reference_U = tr$U, posterior_mean = NA_real_, q025 = NA_real_, q975 = NA_real_,
      posterior_scale = expected_scale, reference_center = center, reference_sd = scale,
      plot_reference = (tr$U - center) / scale,
      plot_estimate = NA_real_, plot_q025 = NA_real_, plot_q975 = NA_real_)
    j <- match(z$cohort_id, tr$cohort_id)
    ocean_diag_require(!anyNA(j) && all(z$C == tr$C[j]) &&
      all(is.finite(as.matrix(z[, c("posterior_mean", "q025", "q975")]))) &&
      all(z$q975 > z$q025), "Invalid saved U posterior display intervals.")
    d$posterior_mean[j] <- z$posterior_mean
    d$q025[j] <- z$q025; d$q975[j] <- z$q975
    # O0 is left in its ORIGINAL model coordinate. Standardizing its posterior
    # or fitting a truth-based affine alignment would conceal scale ambiguity.
    shift <- if (s$cell == "O0") 0 else center
    divisor <- if (s$cell == "O0") 1 else scale
    d$plot_estimate[j] <- (z$posterior_mean - shift) / divisor
    d$plot_q025[j] <- (z$q025 - shift) / divisor
    d$plot_q975[j] <- (z$q975 - shift) / divisor
    d$plot_estimate[observed] <- d$plot_reference[observed]
    # NA posterior summaries on observed rows prevent them being mistaken for
    # inferred values or entering imputation scores. Only display values coincide.
    parts[[i]] <- d
    scales[[i]] <- data.frame(fold = s$fold, cell = s$cell,
      reference_center = center, reference_sd = scale, posterior_scale = expected_scale,
      interpretation = if (s$cell == "O0") "Unanchored model coordinate; illustrative only" else
        "Raw posterior transformed to training-reference SD units")
    counts[[i]] <- data.frame(fold = s$fold, cell = s$cell, n_train = nrow(tr),
      n_masked = sum(!observed), n_observed = sum(observed), n_common_hidden_scored = if (s$cell == "O0") 0L else 80L)
  }
  list(points = do.call(rbind, parts), scales = do.call(rbind, scales), counts = do.call(rbind, counts))
}

ocean_diag_u <- function(output_dir, inputs, sources) {
  points <- list()
  for (i in which(sources$cell != "O0")) {
    s <- sources[i, ]; tr <- inputs$train[[s$fold]]
    z <- ocean_diag_csv(file.path(output_dir, s$source, "U_posterior_intervals.csv"),
      c("cohort_id", "C", "posterior_mean", "q025", "q975", "scale"))
    mask <- tr[[if (s$cell == "O20") "calib_o20" else "calib_o50"]]
    ocean_diag_require(!anyDuplicated(z$cohort_id) && setequal(z$cohort_id, tr$cohort_id[!mask]) &&
      all(z$scale == "raw"), "Saved U intervals do not match the training observation mask.")
    # Nested masks imply that !calib_o50 is hidden under BOTH fitted models.
    ref <- tr[!tr$calib_o50, ]; z <- z[match(ref$cohort_id, z$cohort_id), ]
    ocean_diag_require(nrow(z) == 80L && all(z$C == ref$C) &&
      all(is.finite(as.matrix(z[, c("posterior_mean", "q025", "q975")]))) &&
      all(z$q975 > z$q025), "Invalid common-missing U intervals.")
    points[[length(points) + 1L]] <- data.frame(fold = s$fold, cell = s$cell,
      cohort_id = ref$cohort_id, C = ref$C, reference_U = ref$U,
      posterior_mean = z$posterior_mean, q025 = z$q025, q975 = z$q975)
  }
  points <- do.call(rbind, points)
  scores <- list()
  for (fold in 1:5) for (cell in c("O20", "O50")) for (C in 0:4) {
    z <- points[points$fold == fold & points$cell == cell & (C == 0 | points$C == C), ]
    scores[[length(scores) + 1L]] <- data.frame(fold = fold, cell = cell,
      group = if (C == 0) "overall" else paste0("class", C), ocean_diag_five_metrics(z))
  }
  scores <- do.call(rbind, scores)
  display <- ocean_diag_u_display(output_dir, inputs, sources)
  list(points = points, fold_metrics = scores, summary = ocean_diag_fold_summary(scores,
    c("cell", "group"), c("Bias", "RMSE", "MAE", "Coverage95", "Width95")),
    display_points = display$points, display_scales = display$scales, display_counts = display$counts)
}

# 3. Y-blind display domains for m(x,c) and f(x,u) ------------------------------
ocean_diag_grids <- function(cohort) {
  # Fixed S=34.20. Within each recorded band, select the 25
  # profiles closest in salinity; use their central 80% log-pressure range.
  # Ties are resolved by frozen cohort ID, without looking at Y.
  parts <- neighbors <- list()
  for (C in 1:4) {
    d <- cohort[cohort$C == C, ]
    near <- d[order(abs(d$G2salinity - 34.20), d$cohort_id)[1:25], ]
    bounds <- stats::quantile(near$log_pressure, c(.1, .9), type = 7)
    lp <- seq(bounds[1], bounds[2], length.out = 9L)
    parts[[C]] <- data.frame(design = "candidate_A", C = C,
      G2salinity = 34.20, log_pressure = lp, pressure = exp(lp))
    neighbors[[C]] <- data.frame(design = "candidate_A", class = C, cohort_id = near$cohort_id)
  }
  m <- do.call(rbind, parts)
  # Selected fixed X: S=34.25, P=650. Reference U is used ONLY to set an
  # offline illustration domain and rug; it is never supplied as calibration.
  center <- c(34.25, log(650))
  x <- as.matrix(cohort[, c("G2salinity", "log_pressure")])
  distance <- sqrt(rowSums(sweep(sweep(x, 2, center, "-"), 2, apply(x, 2, stats::sd), "/")^2))
  near <- cohort[order(distance, cohort$cohort_id)[1:30], ]
  bounds <- stats::quantile(near$U, c(.1, .9), type = 7)
  f <- data.frame(G2salinity = center[1], log_pressure = center[2], pressure = 650,
    U = seq(bounds[1], bounds[2], length.out = 41L))
  same_x <- data.frame(design = "same_x", C = 1:4, G2salinity = center[1],
    log_pressure = center[2], pressure = 650)
  m <- rbind(m, same_x); m$grid_id <- sprintf("m_%02d", seq_len(nrow(m)))
  f$grid_id <- sprintf("f_%02d", seq_len(nrow(f)))
  list(m = m, f = f, near = near, m_neighbors = do.call(rbind, neighbors),
    support = data.frame(S = 34.25, pressure = 650, n_neighbors = 30L,
      U_low = bounds[1], U_high = bounds[2],
      C1 = sum(near$C == 1L), C2 = sum(near$C == 2L), C3 = sum(near$C == 3L), C4 = sum(near$C == 4L)))
}

ocean_diag_balanced_ids <- function(fit, n_per_chain) {
  lens <- vapply(fit$mcmc$samples_by_chain$u, nrow, integer(1))
  ocean_diag_require(length(lens) == 4L && all(lens >= n_per_chain),
    "Final diagnostics require four retained paper chains.")
  offsets <- c(0L, head(cumsum(lens), -1L))
  unlist(Map(function(n, off) off + unique(round(seq(1, n, length.out = n_per_chain))), lens, offsets),
    use.names = FALSE)
}

ocean_diag_point_summary <- function(draws, means, grid) {
  ocean_diag_require(is.matrix(draws) && is.matrix(means) && ncol(draws) == nrow(grid) &&
    ncol(means) == nrow(grid) && all(is.finite(draws)) && all(is.finite(means)),
    "Package returned invalid function draws.")
  q <- apply(draws, 2L, stats::quantile, probs = c(.025, .975), type = 7)
  # Average conditional means for the point estimate; do not add GP simulation
  # noise to this estimate. Quantiles include conditional GP uncertainty.
  data.frame(grid, posterior_mean = colMeans(means), q025 = q[1, ], q975 = q[2, ],
    Width95 = q[2, ] - q[1, ])
}

# 4. Package-only functional posterior calculations, cached separately per fit -
ocean_diag_functionals <- function(i, sources, output_dir, inputs, grids, settings,
  cache_dir, overwrite) {
  s <- sources[i, ]; path <- file.path(output_dir, s$source)
  fit_path <- file.path(path, s$fit_file)
  prediction_path <- file.path(path, "prediction_distributions.rds")
  tr <- inputs$train[[s$fold]]
  mg <- grids$m[grids$m$design == "candidate_A" | s$cell != "O0", ]
  identity <- list(version = "ocean-diagnostics-1", package = as.character(utils::packageVersion("eivGP")),
    source = s, hashes = unname(tools::md5sum(c(fit_path, prediction_path))),
    train = tr[, c("cohort_id", "G2salinity", "log_pressure", "C", "Y", "U", "calib_o20", "calib_o50")],
    m_grid = mg, f_grid = grids$f, settings = settings,
    code = list(body(ocean_diag_functionals), body(ocean_diag_balanced_ids), body(ocean_diag_point_summary)))
  cache_path <- file.path(cache_dir, sprintf("fold%02d_%s.rds", s$fold, s$cell))
  if (file.exists(cache_path) && !overwrite) {
    saved <- readRDS(cache_path)
    ocean_diag_require(identical(saved$identity, identity),
      paste("Diagnostic cache inputs changed. Choose a new report_dir or set overwrite=TRUE:", cache_path))
    return(saved$result)
  }
  message("Function diagnostics: fold ", s$fold, "/", s$cell)
  fit <- readRDS(fit_path)
  mask <- switch(s$cell, O0 = rep(FALSE, nrow(tr)), O20 = tr$calib_o20, O50 = tr$calib_o50)
  xmat <- function(d) as.matrix(d[, c("G2salinity", "log_pressure"), drop = FALSE])
  ocean_diag_require(inherits(fit, "eivgp_fit") && identical(fit$sampler_version, "0.3.1") &&
    isTRUE(all.equal(unname(fit$data$x_raw), unname(xmat(tr)), tolerance = 1e-10)) &&
    isTRUE(all.equal(as.numeric(fit$data$y_raw), tr$Y, tolerance = 1e-10)) &&
    identical(as.integer(fit$data$c_ord), tr$C) &&
    setequal(fit$data$calib_idx, which(mask)) && setequal(fit$data$miss_idx, which(!mask)),
    "Fit object does not match the frozen training rows, model version, or calibration mask.")
  if (any(mask)) ocean_diag_require(
    isTRUE(all.equal(as.numeric(fit$data$u_obs_raw[mask]), tr$U[mask], tolerance = 1e-10)),
    "Exact oxygen measurements in the fitted object disagree with the calibration rows.")
  ids <- ocean_diag_balanced_ids(fit, settings$mean_states_per_chain)
  seed <- settings$seed + 100L * s$fold + match(s$cell, c("O0", "O20", "O50"))
  md <- eivGP::predict_eivgp(fit, new_X = xmat(mg), new_C = matrix(mg$C, ncol = 1),
    target = "mean", draw_ids = ids, n_latent = settings$n_latent,
    n_per_draw = settings$mean_gp_draws, include_process_uncertainty = TRUE,
    joint = FALSE, return_components = TRUE, seed = seed)
  conditional_means <- attr(md, "conditional_means", exact = TRUE)
  ms <- ocean_diag_point_summary(md, conditional_means, mg)
  ms$fold <- s$fold; ms$cell <- s$cell
  # Small numerical integration sensitivity check, not a recovery score. The
  # same 120 parameter states are evaluated with 256 and 512 latent nodes.
  check_columns <- match(1:4, mg$C)
  check_states <- unique(round(seq(1, length(ids), length.out = 120L)))
  check_grid <- mg[check_columns, ]
  alt <- eivGP::predict_eivgp(fit, new_X = xmat(check_grid),
    new_C = matrix(check_grid$C, ncol = 1), target = "mean", draw_ids = ids[check_states],
    n_latent = 2L * settings$n_latent, n_per_draw = 1L,
    include_process_uncertainty = FALSE, return_components = TRUE, seed = seed + 20000L)
  delta <- attr(alt, "conditional_means", exact = TRUE) -
    conditional_means[check_states, check_columns, drop = FALSE]
  check <- data.frame(fold = s$fold, cell = s$cell, check_grid,
    delta_mean_double_nodes = colMeans(delta), paired_state_mcse = apply(delta, 2, stats::sd) / sqrt(nrow(delta)))
  fs <- NULL; fids <- integer()
  if (s$cell != "O0") {
    fids <- readRDS(prediction_path)$draw_ids
    ocean_diag_require(length(fids) >= 1000L && !anyDuplicated(fids) &&
      all(fids >= 1 & fids <= nrow(fit$mcmc$samples_tau)), "Unexpected saved posterior state IDs.")
    fm <- eivGP::predict_eivgp(fit, new_X = xmat(grids$f), new_U = matrix(grids$f$U, ncol = 1),
      new_U_scale = "raw", target = "surface", draw_ids = fids, include_process_uncertainty = FALSE)
    fd <- eivGP::predict_eivgp(fit, new_X = xmat(grids$f), new_U = matrix(grids$f$U, ncol = 1),
      new_U_scale = "raw", target = "surface", draw_ids = fids, n_per_draw = settings$surface_gp_draws,
      include_process_uncertainty = TRUE, joint = FALSE, seed = seed + 10000L)
    fs <- ocean_diag_point_summary(fd, fm, grids$f)
    fs$fold <- s$fold; fs$cell <- s$cell
  }
  result <- list(m = ms, f = fs, integration_check = check,
    draw_ids = list(mean = ids, surface = fids), seed = seed)
  saveRDS(list(identity = identity, result = result), cache_path)
  result
}

# 5. Original held-out Y tasks; reuse saved scores without new predictions -----
ocean_diag_prediction <- function(output_dir, inputs, sources) {
  metrics <- c("RMSE", "MAE", "CRPS", "NLPD", "Coverage95", "Width95", "IntervalScore95")
  rows <- masks <- list()
  for (fold in 1:5) {
    path <- file.path(output_dir, "full5_class4", sprintf("fold%02d", fold))
    z <- ocean_diag_csv(file.path(path, "baseline_metrics.csv"), c("fold", "cell", "method", "n_test", metrics))
    ocean_diag_require(nrow(z) == 3L && setequal(z$method, c("UC-GP", "LVGP", "EzGP")) &&
      all(z$fold == fold & z$n_test == 40L), "Unexpected baseline scores.")
    z$task <- "class_only"; z$test_exact_fraction <- 0
    rows[[length(rows) + 1L]] <- z[, c("fold", "cell", "method", "task", "test_exact_fraction", "n_test", metrics)]
    mask <- ocean_diag_csv(file.path(path, "test_U_masks.csv"),
      c("cohort_id", "C", "observed_U_O20", "observed_U_O50"))
    te <- inputs$test[[fold]]
    ocean_diag_require(nrow(mask) == 40L && !anyDuplicated(mask$cohort_id) &&
      setequal(mask$cohort_id, te$cohort_id), "Unexpected held-out oxygen mask rows.")
    mask <- mask[match(te$cohort_id, mask$cohort_id), ]
    m20 <- ocean_diag_flag(mask$observed_U_O20); m50 <- ocean_diag_flag(mask$observed_U_O50)
    ocean_diag_require(all(mask$C == te$C) && all(m20 <= m50) &&
      all(table(factor(te$C[m20], levels = 1:4)) == 2L) &&
      all(table(factor(te$C[m50], levels = 1:4)) == 5L), "Held-out masks changed.")
    masks[[fold]] <- data.frame(fold = fold, mask)
    for (i in which(sources$fold == fold)) {
      s <- sources[i, ]
      z <- ocean_diag_csv(file.path(output_dir, s$source, "predictive_metrics.csv"),
        c("fold", "cell", "method", "n_test", metrics))
      ocean_diag_require(nrow(z) == 2L && setequal(z$method, c("EIV-GP xc", "EIV-GP mixed")) &&
        all(z$fold == fold & z$cell == s$cell & z$n_test == 40L), "Unexpected EIV prediction scores.")
      saved <- readRDS(file.path(output_dir, s$source, "prediction_distributions.rds"))
      known <- switch(s$cell, O0 = rep(FALSE, 40L), O20 = m20, O50 = m50)
      ocean_diag_require(identical(as.logical(saved$known), known) &&
        ncol(saved$xc) == 40L && ncol(saved$mixed) == 40L &&
        identical(saved$xc[, !known, drop = FALSE], saved$mixed[, !known, drop = FALSE]),
        "Saved mixed predictions disagree with the frozen exact-oxygen masks.")
      if (s$cell == "O0") z <- z[z$method == "EIV-GP xc", ] # O0 mixed is identical.
      z$task <- ifelse(z$method == "EIV-GP xc", "class_only", "matched_mixed")
      z$test_exact_fraction <- ifelse(z$task == "class_only", 0, if (s$cell == "O20") .2 else .5)
      z$method <- "EIV-GP"
      rows[[length(rows) + 1L]] <- z[, c("fold", "cell", "method", "task", "test_exact_fraction", "n_test", metrics)]
    }
  }
  raw <- do.call(rbind, rows)
  normalized <- raw
  s_y <- stats::sd(inputs$cohort$Y) # Scale by the SD of Y in the frozen cohort.
  for (metric in c("RMSE", "MAE", "CRPS", "Width95", "IntervalScore95"))
    normalized[[metric]] <- normalized[[metric]] / s_y
  normalized$NLPD <- normalized$NLPD - log(s_y) # Density transformation, not division.
  list(fold_metrics_raw = raw, fold_metrics_SD_Y = normalized,
    summary_raw = ocean_diag_fold_summary(raw, c("task", "method", "cell"), metrics),
    summary_SD_Y = ocean_diag_fold_summary(normalized, c("task", "method", "cell"), metrics),
    test_masks = do.call(rbind, masks), SD_Y = s_y)
}

# 6. Figures: one fixed fold in the main text, all folds in the appendix --------
ocean_diag_plots <- function(u, m, f, grids, representative_fold, baseline_same_x = NULL) {
  colors <- c(O0 = "#777777", O20 = "#2365a3", O50 = "#c76326")
  labels <- c(O0 = "0% exact training U", O20 = "20% exact training U", O50 = "50% exact training U")
  wrap <- function(x, width = 105L) paste(strwrap(x, width = width), collapse = "\n")
  style <- function(p, fill = FALSE) {
    p <- p + ggplot2::theme_minimal(base_size = 11) +
      ggplot2::theme(legend.position = "bottom", panel.grid.minor = ggplot2::element_blank(),
        plot.caption = ggplot2::element_text(hjust = 0)) +
      ggplot2::scale_color_manual(values = colors, labels = labels, name = NULL)
    if (fill) p <- p + ggplot2::scale_fill_manual(values = colors, labels = labels, name = NULL)
    p
  }
  interval_note <- "Pointwise 95% credible intervals for the mean; response noise excluded."
  u_plot <- function(all_folds) {
    d <- u$display_points
    if (!all_folds) d <- d[d$fold == representative_fold, ]
    d$class <- factor(d$C, levels = 1:4)
    d$role <- factor(ifelse(d$observed, "Observed calibration (not imputed)", "Masked training U"),
      levels = c("Masked training U", "Observed calibration (not imputed)"))
    panel_labels <- c(O0 = "0% calibration\n160 masked / 0 observed\nUnanchored model coordinate",
      O20 = "20% calibration\n128 masked / 32 observed\nTraining-reference SD units",
      O50 = "50% calibration\n80 masked / 80 observed\nTraining-reference SD units")
    d$calibration <- factor(d$cell, levels = names(panel_labels), labels = unname(panel_labels))
    d$fold_label <- factor(paste("Fold", d$fold), levels = paste("Fold", 1:5))
    p <- ggplot2::ggplot(d, ggplot2::aes(plot_reference, plot_estimate, color = class)) +
      ggplot2::geom_abline(slope = 1, intercept = 0, linetype = 2, color = "grey50") +
      ggplot2::geom_linerange(data = d[!d$observed, ],
        ggplot2::aes(ymin = plot_q025, ymax = plot_q975),
        alpha = .20, linewidth = .35, show.legend = FALSE) +
      ggplot2::geom_point(ggplot2::aes(shape = role), size = 1.5, alpha = .85, stroke = .65) +
      ggplot2::scale_color_manual(values = c("1" = "#1b9e77", "2" = "#d95f02",
        "3" = "#7570b3", "4" = "#0072B2"), name = "Oxygen class", drop = FALSE) +
      ggplot2::scale_shape_manual(values = c("Masked training U" = 16,
        "Observed calibration (not imputed)" = 1), name = NULL) +
      ggplot2::guides(shape = ggplot2::guide_legend(order = 1, override.aes = list(color = "grey25", alpha = 1)),
        color = ggplot2::guide_legend(order = 2, override.aes = list(shape = 16, alpha = 1))) +
      ggplot2::coord_cartesian(xlim = range(d$plot_reference) + c(-.1, .1),
        ylim = range(c(d$plot_q025, d$plot_q975, d$plot_estimate), na.rm = TRUE) + c(-.1, .1)) +
      ggplot2::labs(x = "Reference oxygen (training-reference SD units)", y = "U estimate (display units)",
        title = "Training oxygen: masked estimates and observed calibration inputs",
        subtitle = if (all_folds) "All five folds; 160 training profiles in every panel" else
          paste("Fixed fold", representative_fold, "| 160 training profiles in every panel"),
        caption = wrap(paste("Solid: masked training U with 95% posterior intervals. Hollow: exact calibration, not imputed.",
          "Reference is centered/scaled by each fold's full training reference U; O20/O50 posteriors use that same transformation.",
          "O0 stays on its unanchored model scale: identity-line agreement is illustrative, not physical-scale recovery.",
          "Training Y participates. Scores remain restricted to the common 80 hidden rows at O20/O50."))) +
      ggplot2::theme_minimal(base_size = 11) +
      ggplot2::theme(legend.position = "bottom", panel.grid.minor = ggplot2::element_blank(),
        plot.caption = ggplot2::element_text(hjust = 0))
    if (all_folds) p + ggplot2::facet_grid(fold_label ~ calibration) else
      p + ggplot2::facet_wrap(~calibration, nrow = 1)
  }
  m_plot <- function(all_folds) {
    d <- m[m$design == "candidate_A", ]
    if (!all_folds) d <- d[d$fold == representative_fold, ]
    d$class <- factor(d$C, levels = 1:4, labels = paste("Class", 1:4))
    p <- ggplot2::ggplot(d, ggplot2::aes(pressure, posterior_mean, color = cell, fill = cell)) +
      ggplot2::geom_ribbon(ggplot2::aes(ymin = q025, ymax = q975), alpha = .13, color = NA) +
      ggplot2::geom_line(linewidth = .65) +
      ggplot2::labs(x = "Pressure (dbar)", y = "Mean silicate m(x,c) (micromol/kg)",
        title = "Category-averaged response at salinity 34.20",
        subtitle = if (all_folds) "All five folds; each band belongs to its own fitted model" else
          paste("Fixed fold", representative_fold, "| candidate A"),
        caption = wrap(paste(interval_note, "Pressure ranges follow the prespecified local display rule.")))
    p <- if (all_folds) p + ggplot2::facet_wrap(~fold + class, ncol = 4, scales = "free") else
      p + ggplot2::facet_wrap(~class, ncol = 2, scales = "free")
    style(p, fill = TRUE)
  }
  fm_plot <- function(all_folds) {
    df <- f; dm <- m[m$design == "same_x", ]; db <- baseline_same_x
    if (!all_folds) {
      df <- df[df$fold == representative_fold, ]; dm <- dm[dm$fold == representative_fold, ]
      if(!is.null(db)) db <- db[db$fold == representative_fold, ]
    }
    df$model <- paste0("EIV-GP ",df$cell); dm$model <- paste0("EIV-GP ",dm$cell)
    means <- ocean_bind(list(dm,db))
    palette <- c("EIV-GP O20"="#2365a3","EIV-GP O50"="#c76326",
      "UC-GP (x,c)"="#7b3294","LVGP (x,c)"="#008837","EzGP (x,c)"="#a6761d")
    models <- names(palette)[names(palette) %in% means$model]
    offsets <- setNames(seq(-.26,.26,length.out=length(models)),models)
    means$class_position <- means$C + offsets[means$model]
    means$type <- ifelse(grepl("^EIV",means$model),"EIV posterior mean","Baseline fitted mean")
    # Match the physical response scale across panels; give the discrete panel
    # less width, and use larger points/capped intervals to reduce empty space.
    ylim <- range(c(df$q025,df$q975,dm$q025,dm$q975,means$posterior_mean),na.rm=TRUE)
    rug <- grids$near[grids$near$U >= min(df$U) & grids$near$U <= max(df$U), ]
    pf <- ggplot2::ggplot(df, ggplot2::aes(U, posterior_mean, color = model, fill = model)) +
      ggplot2::geom_vline(xintercept = c(50, 100, 150), linetype = 3, color = "grey65", linewidth = .4) +
      ggplot2::geom_ribbon(ggplot2::aes(ymin = q025, ymax = q975), alpha = .14, color = NA) +
      ggplot2::geom_line(linewidth = .9) +
      ggplot2::geom_rug(data = rug, ggplot2::aes(x = U), inherit.aes = FALSE,
        sides = "b", color = "grey45", alpha = .65) +
      ggplot2::scale_x_continuous(expand=ggplot2::expansion(mult=.025)) +
      ggplot2::coord_cartesian(xlim = range(df$U), ylim = ylim) +
      ggplot2::labs(x = "Oxygen (micromol/kg)", y = "Mean silicate (micromol/kg)", title = "(a) Continuous response f(x,u)")
    pm <- ggplot2::ggplot(means, ggplot2::aes(class_position,posterior_mean,color=model)) +
      ggplot2::geom_errorbar(data=means[grepl("^EIV",means$model),],
        ggplot2::aes(ymin=q025,ymax=q975),width=.10,linewidth=.85,show.legend=FALSE) +
      ggplot2::geom_point(ggplot2::aes(shape=type),size=3,stroke=.8) +
      ggplot2::scale_shape_manual(values=c("EIV posterior mean"=16,"Baseline fitted mean"=17),guide="none") +
      ggplot2::scale_x_continuous(breaks=1:4,limits=c(.5,4.5),expand=c(0,0)) +
      ggplot2::coord_cartesian(ylim = ylim) +
      ggplot2::labs(x = "Recorded oxygen class", y = NULL, title = "(b) Category mean m(x,c)")
    theme_fm <- function(p) p + ggplot2::theme_minimal(base_size=12) +
      ggplot2::scale_color_manual(values=palette,limits=models,breaks=models,drop=FALSE,name=NULL) +
      ggplot2::scale_fill_manual(values=palette,guide="none") +
      ggplot2::theme(panel.grid.minor=ggplot2::element_blank(),legend.position="bottom",
        plot.title=ggplot2::element_text(size=12),legend.key.width=grid::unit(1.1,"cm")) +
      ggplot2::guides(color=ggplot2::guide_legend(nrow=1,override.aes=list(shape=16,linewidth=.8)))
    if (all_folds) { pf <- pf + ggplot2::facet_grid(fold ~ .); pm <- pm + ggplot2::facet_grid(fold ~ .) }
    patchwork::wrap_plots(theme_fm(pf)+ggplot2::guides(color="none"),
      theme_fm(pm)+ggplot2::guides(color=ggplot2::guide_legend(nrow=1,
        override.aes=list(shape=ifelse(grepl("^EIV",models),16,17),size=3))),
      widths=c(1.65,1),guides="collect") +
      patchwork::plot_annotation(title = "Salinity 34.25, pressure 650 dbar",
        subtitle = if (all_folds) "All five folds; separate fitted posteriors" else paste("Fixed fold", representative_fold),
        caption = paste("EIV bars/shading: pointwise 95% mean credible intervals. Triangles: original (x,c) baseline fitted means.",
          "\nm integrates the full fitted oxygen distribution. Dotted lines: recorded band boundaries; rug: nearby reference oxygen.")) &
      ggplot2::theme(legend.position="bottom",plot.caption=ggplot2::element_text(hjust=0,size=10))
  }
  list(u_main = u_plot(FALSE), m_main = m_plot(FALSE), fm_main = fm_plot(FALSE),
    u_all_folds = u_plot(TRUE), m_all_folds = m_plot(TRUE), fm_all_folds = fm_plot(TRUE))
}

# 7. Rmd entry point: postprocess and export; no fitting ------------------------
ocean_build_diagnostics <- function(output_dir, data_dir = file.path(output_dir, "data"),
  report_dir = file.path(output_dir, "diagnostics"), representative_fold = 1L,
  workers = 1L, overwrite = FALSE) {
  for (package in c("eivGP", "ggplot2", "patchwork"))
    ocean_diag_require(requireNamespace(package, quietly = TRUE), paste("Required package:", package))
  ocean_diag_require(as.character(utils::packageVersion("eivGP")) == "0.3.1",
    "This paper component was validated with eivGP 0.3.1; use the paper package environment.")
  ocean_diag_require(length(representative_fold) == 1L && representative_fold %in% 1:5,
    "representative_fold must be one of 1:5; the default main display uses fold 1.")
  ocean_diag_require(length(workers) == 1L && is.finite(workers) && workers >= 1 && workers == as.integer(workers),
    "workers must be a positive integer.")
  output_dir <- normalizePath(output_dir, mustWork = TRUE)
  inputs <- ocean_diag_inputs(data_dir)
  sources <- ocean_diag_sources(output_dir)
  grids <- ocean_diag_grids(inputs$cohort)
  u <- ocean_diag_u(output_dir, inputs, sources)
  prediction <- ocean_diag_prediction(output_dir, inputs, sources)
  settings <- list(mean_states_per_chain = 60L, n_latent = 256L, mean_gp_draws = 20L,
    surface_gp_draws = 4L, seed = 202610080L)
  for (subdir in c("tables", "figures", "internal", "cache"))
    dir.create(file.path(report_dir, subdir), recursive = TRUE, showWarnings = FALSE)
  run_one <- function(i) ocean_diag_functionals(i, sources, output_dir, inputs, grids, settings,
    cache_dir = file.path(report_dir, "cache"), overwrite = overwrite)
  if (workers > 1L && .Platform$OS.type != "windows") {
    parts <- parallel::mclapply(seq_len(nrow(sources)), run_one, mc.cores = min(workers, nrow(sources)),
      mc.preschedule = FALSE, mc.set.seed = FALSE) # Explicit seed for every package call.
    ocean_diag_require(!any(vapply(parts, inherits, logical(1), "try-error")),
      "At least one function calculation failed. Inspect messages; completed cell caches can be reused.")
  } else parts <- lapply(seq_len(nrow(sources)), run_one)
  m <- do.call(rbind, lapply(parts, `[[`, "m"))
  f <- do.call(rbind, lapply(parts, `[[`, "f"))
  checks <- do.call(rbind, lapply(parts, `[[`, "integration_check"))
  ocean_diag_require(nrow(m) == 15L * 36L + 10L * 4L && nrow(f) == 10L * 41L &&
    all(m$Width95 > 0) && all(f$Width95 > 0), "Incomplete functional diagnostics.")
  # Internal width averages are descriptive precision checks, not true m/f
  # recovery metrics or new main-text tables. Equal weights over the fixed grid.
  widths <- stats::aggregate(Width95 ~ fold + cell + C, m[m$design == "candidate_A", ], mean)
  width_summary <- ocean_diag_fold_summary(widths, c("cell", "C"), "Width95")
  diag <- lapply(seq_len(nrow(sources)), function(i) {
    s <- sources[i, ]; d <- utils::read.csv(file.path(output_dir, s$source, "convergence_summary.csv"))
    data.frame(fold = s$fold, cell = s$cell, d)
  })
  baseline_paths <- file.path(output_dir,"full5_class4",sprintf("fold%02d",1:5),"baseline_package_outputs.rds")
  baseline_same_x <- if(all(file.exists(baseline_paths)))
    ocean_baseline_mean_grid(output_dir,grids$m[grids$m$design=="same_x",]) else NULL
  plots <- ocean_diag_plots(u,m,f,grids,representative_fold,baseline_same_x)
  write_csv <- function(x, subdir, name) utils::write.csv(x, file.path(report_dir, subdir, name), row.names = FALSE)
  write_csv(u$points, "tables", "u_common_missing_points.csv")
  write_csv(u$display_points, "tables", "u_display_points.csv")
  write_csv(u$display_scales, "internal", "u_display_scales.csv")
  write_csv(u$display_counts, "internal", "u_display_counts.csv")
  write_csv(u$fold_metrics, "tables", "u_fold_metrics.csv")
  write_csv(u$summary, "tables", "u_summary.csv")
  write_csv(m, "tables", "m_points.csv"); write_csv(f, "tables", "f_points.csv")
  if(!is.null(baseline_same_x)) write_csv(baseline_same_x,"tables","baseline_m_same_x.csv")
  paper_u <- ocean_paper_u_tables(u$summary)
  ocean_export_paper_table(paper_u$main,"paper_u_imputation_main","Oxygen imputation: common 80 masked training profiles",file.path(report_dir,"tables"),paper_u$note)
  ocean_export_paper_table(paper_u$class,"paper_u_imputation_class","Class-specific oxygen imputation",file.path(report_dir,"tables"),paper_u$note)
  write_csv(prediction$fold_metrics_raw, "tables", "y_prediction_fold_raw.csv")
  write_csv(prediction$fold_metrics_SD_Y, "tables", "y_prediction_fold_SD_Y.csv")
  write_csv(prediction$summary_raw, "tables", "y_prediction_summary_raw.csv")
  write_csv(prediction$summary_SD_Y, "tables", "y_prediction_summary_SD_Y.csv")
  write_csv(inputs$audit, "internal", "input_audit.csv")
  write_csv(sources, "internal", "fit_sources.csv")
  write_csv(grids$m, "internal", "m_grid.csv"); write_csv(grids$f, "internal", "f_grid.csv")
  write_csv(grids$m_neighbors, "internal", "m_display_neighbors.csv")
  write_csv(grids$near, "internal", "f_display_neighbors.csv")
  write_csv(grids$support, "internal", "f_display_support.csv")
  write_csv(checks, "internal", "m_integration_sensitivity.csv")
  write_csv(widths, "internal", "m_width_fold.csv")
  write_csv(width_summary, "internal", "m_width_summary.csv")
  write_csv(prediction$test_masks, "internal", "test_U_masks.csv")
  write_csv(do.call(rbind, diag), "internal", "mcmc_existing_summary.csv")
  dimensions <- list(u_main = c(13.2, 5.2), m_main = c(10, 7), fm_main = c(11, 4.9),
    u_all_folds = c(13.2, 14.5), m_all_folds = c(13, 13), fm_all_folds = c(11, 14))
  for (name in names(plots)) for (extension in c("png", "pdf"))
    ggplot2::ggsave(file.path(report_dir, "figures", paste0(name, ".", extension)), plots[[name]],
      width = dimensions[[name]][1], height = dimensions[[name]][2], dpi = 180, bg = "white")
  provenance <- list(settings = settings, representative_fold = representative_fold,
    package = as.character(utils::packageVersion("eivGP")), upstream_commit = "671ce98339807cd3a13347e8468f36acf097d50a",
    input_hashes = tools::md5sum(file.path(data_dir, c("cohort_manifest.csv",
      sprintf("fold%02d_train.csv", 1:5), sprintf("fold%02d_test.csv", 1:5)))),
    source_hashes = lapply(seq_len(nrow(sources)), function(i) {
      s <- sources[i, ]
      tools::md5sum(file.path(output_dir, s$source, c(s$fit_file,
        "prediction_distributions.rds", "predictive_metrics.csv", "U_posterior_intervals.csv")))
    }),
    draw_ids = lapply(parts, `[[`, "draw_ids"), SD_Y = prediction$SD_Y,
    main_targets = c("training U imputation", "m(x,c)", "f(x,u) and same-X m(x,c)", "held-out Y prediction"))
  saveRDS(provenance, file.path(report_dir, "internal", "provenance.rds"))
  capture.output(utils::sessionInfo(), file = file.path(report_dir, "internal", "sessionInfo.txt"))
  invisible(list(u = u, m = m, f = f, prediction = prediction, grids = grids, plots = plots, baseline_same_x=baseline_same_x,
    integration_check = checks, sources = sources, report_dir = report_dir, provenance = provenance))
}

# 8. Compact HTML/Rmd presentation: the four diagnostic blocks ----------
ocean_diag_display_table <- function(summary, groups, metrics, digits = 3L) {
  rows <- unique(summary[, groups, drop = FALSE])
  key <- function(d) do.call(paste, c(d[groups], sep = "/"))
  for (metric in metrics) {
    z <- summary[summary$metric == metric, ]
    i <- match(key(rows), key(z))
    rows[[metric]] <- sprintf(paste0("%.", digits, "f +/- %.", digits, "f"), z$mean[i], z$sd[i])
  }
  rownames(rows) <- NULL
  rows
}

ocean_show_diagnostics <- function(report, appendix = FALSE, show_prediction = TRUE) {
  ocean_diag_require(requireNamespace("knitr", quietly = TRUE), "knitr is required for notebook presentation.")
  table <- function(x) cat(as.character(knitr::kable(x, format = "pipe", row.names = FALSE)), sep = "\n")
  cat("\n\n### Training oxygen imputation\n\n",
    "Both calibration levels are scored on the same 80 hidden training observations per fold. ",
    "Reference oxygen is the masked source measurement. Entries are the mean +/- sample SD of five fold scores; ",
    "training sets overlap, and these SDs are descriptive. Training Y contributes to imputation.\n\n", sep = "")
  table(ocean_paper_u_tables(report$u$summary)$main)
  cat("\n\nErrors and widths use micromol/kg; coverage is percent. O0 has no anchored raw oxygen scale.\n\n")
  cat("\n\nThe figure displays all 160 training rows per calibration level (O0/O20/O50). ",
    "Solid points with intervals are masked inputs; hollow points are observed calibration values, not imputations. ",
    "Reference U is centered/scaled by the full training-reference mean/SD within each fold. ",
    "O20/O50 posteriors are transformed to those same units; O0 remains in its original unanchored model coordinate. ",
    "The O0 identity-line comparison is illustrative only and supplies no physical-scale recovery score. ",
    "Hidden references enter only this offline display, never fitting or posterior alignment. ",
    "Displayed masked-row counts are 160/128/80; the table still scores the same common 80 rows at O20/O50.\n\n", sep = "")
  print(report$plots$u_main)
  cat("\n\n### Category-averaged silicate response m(x,c)\n\n",
    "The pressure slice fixes salinity at 34.20 and varies pressure within a class-specific local display range. ",
    "The package integrates f against the fitted conditional oxygen distribution. Bands describe uncertainty ",
    "in this mean and exclude response noise. True m is unknown, so no mean-recovery error is reported.\n\n", sep = "")
  print(report$plots$m_main)
  cat("\n\n### Continuous response f(x,u) and the category average at the same X\n\n",
    "Salinity is 34.25 and pressure is 650 dbar. The curve varies oxygen continuously; the four category means ",
    "average over their full fitted conditional oxygen distributions. The displayed oxygen range is the ",
    "10th-90th percentile among 30 nearest reference profiles in standardized (salinity, log pressure). ",
    "Reference oxygen defines only this offline display domain and rug, not additional calibration. ",
    "Nearby profiles indicate the available input support but do not guarantee dense support at every joint input. ",
    "The fitted end-class distributions may extend beyond 0-180 and beyond the plotted range. ",
    "True f is unknown; this is a descriptive fitted association, without a surface-recovery score or a causal interpretation.\n\n", sep = "")
  print(report$plots$fm_main)
  if (show_prediction) {
  cat("\n\n### Held-out silicate prediction\n\n",
    "The original prediction tasks are retained. Class-only uses no exact test oxygen. Matched mixed uses ",
    "20% or 50% exact test oxygen with the corresponding training calibration level. O0 mixed duplicates ",
    "O0 class-only and is shown once. The three categorical comparators use class-only information. ",
    "All entries are five-fold mean +/- SD; no new held-out predictions are generated here.\n\n", sep = "")
  metrics <- c("RMSE", "CRPS", "NLPD", "Coverage95", "Width95", "IntervalScore95")
  for (task in c("class_only", "matched_mixed")) {
    cat("\n\n**", if (task == "class_only") "Class-only" else "Matched mixed", "**\n\n", sep = "")
    s <- report$prediction$summary_SD_Y
    table(ocean_diag_display_table(s[s$task == task, ], c("method", "cell"), metrics))
  }
  cat("\n\nRMSE, CRPS, width and interval score are divided by the frozen cohort SD(Y); ",
    "NLPD is transformed by subtracting log(SD(Y)); coverage is unchanged. ",
    "Raw-unit tables and MAE are also exported. Mixed-versus-baseline differences combine model and input-information differences.\n\n", sep = "")
  }
  if (appendix) {
    cat("\n\n### Appendix: all five fitted folds\n\n",
      "Each interval belongs to one fitted fold. These plots do not pool the folds into a posterior.\n\n", sep = "")
    print(report$plots$u_all_folds); print(report$plots$m_all_folds); print(report$plots$fm_all_folds)
    cat("\n\n**Class-specific U imputation (same five metrics)**\n\n")
    table(ocean_diag_display_table(report$u$summary[report$u$summary$group != "overall", ], c("cell", "group"),
      c("Bias", "RMSE", "MAE", "Coverage95", "Width95")))
  }
  invisible(report)
}

# Native baseline mean predictors; no custom covariance or GP posterior algebra.
ocean_native_mean <- function(method,fit,grid) {
  X <- as.matrix(grid[,c("G2salinity","log_pressure")]); C <- matrix(grid$C,ncol=1L)
  if(method=="UC-GP") {
    d <- eivGP:::make_mixedgp_competitor_data(X,C,4L)
    as.numeric(stats::predict(fit,newdata=d,type="UK",seCompute=FALSE,
      covCompute=FALSE,lightReturn=TRUE)$mean)
  } else {
    d <- eivGP:::make_mixedgp_numeric_matrix(X,C,4L)
    if(method=="LVGP") as.numeric(LVGP::LVGP_predict(d,fit,MSE_on=0)$Y_hat)
    else as.numeric(EzGP::EzGP_predict(d,fit,MSE_on=0)$Y_hat)
  }
}

# The original three baselines predict a category response directly. Evaluate
# their native means on the fixed X grids; no R/W slice is relabeled as m.
ocean_baseline_mean_grid <- function(output_dir,grid,folds=1:5) {
  out <- list()
  for(fold in folds) {
    saved <- readRDS(file.path(output_dir,"full5_class4",sprintf("fold%02d",fold),"baseline_package_outputs.rds"))
    for(method in c("UC-GP","LVGP","EzGP")) {
      mean <- ocean_native_mean(method,saved$fits[[method]],grid)
      stopifnot(length(mean)==nrow(grid),all(is.finite(mean)))
      out[[length(out)+1L]] <- data.frame(grid,posterior_mean=mean,fold,model=paste0(method," (x,c)"))
    }
  }
  ocean_bind(out)
}

# Input-based display selection. Use the fitted grid and its local
# salinity pool; show at most 3 nearest-salinity observations per log-pressure
# bin, within the curve range. Y, residuals and held-out scores are never used.
ocean_m_display_points <- function(cohort,grids,n_bins=6L,per_bin=3L) {
  out <- list()
  for(C in 1:4) {
    g <- grids$m[grids$m$design=="candidate_A" & grids$m$C==C,]
    ids <- grids$m_neighbors$cohort_id[grids$m_neighbors$class==C]
    d <- cohort[cohort$cohort_id %in% ids & cohort$log_pressure>=min(g$log_pressure) &
      cohort$log_pressure<=max(g$log_pressure),]
    d$pressure_bin <- cut(d$log_pressure,breaks=seq(min(g$log_pressure),max(g$log_pressure),length.out=n_bins+1L),
      include.lowest=TRUE,labels=FALSE)
    d$salinity_offset <- d$G2salinity-34.20
    keep <- unlist(lapply(split(seq_len(nrow(d)),d$pressure_bin),function(i)
      head(i[order(abs(d$salinity_offset[i]),d$cohort_id[i])],per_bin)),use.names=FALSE)
    out[[C]] <- d[keep,,drop=FALSE]
  }
  ocean_bind(out)
}

ocean_paper_u_tables <- function(summary) {
  format_u <- function(groups) {
    z <- summary[summary$group %in% groups,]
    keys <- unique(z[c("cell","group")]); keys <- keys[order(keys$cell,keys$group),]
    d <- data.frame(Calibration=paste0(sub("O","",keys$cell),"%"),check.names=FALSE)
    if(length(groups)>1L) d$Class <- sub("class","",keys$group)
    ranked <- as.data.frame(matrix(NA_real_,nrow(keys),3L))
    names(ranked) <- c("Bias","RMSE","MAE")
    eligible <- rep(TRUE,nrow(keys))
    for(metric in c("Bias","RMSE","MAE","Coverage95","Width95")) {
      v <- z[z$metric==metric,]; i <- match(paste(keys$cell,keys$group),paste(v$cell,v$group))
      mult <- if(metric=="Coverage95") 100 else 1
      label <- if(metric=="Coverage95") "Coverage (%)" else if(metric=="Width95") "Width" else metric
      d[[label]] <- ocean_mean_sd(mult*v$mean[i],mult*v$sd[i],if(metric=="Coverage95") 1 else 2)
      eligible <- eligible & v$n_folds[i]==5L
      if(metric %in% names(ranked))
        ranked[[metric]] <- if(metric=="Bias") abs(v$mean[i]) else v$mean[i]
      if(any(v$n_folds[i]!=5L)) d[["n folds"]] <- v$n_folds[i]
    }
    attr(d,"best") <- ocean_best_mask(ranked,keys$group,eligible)
    d
  }
  list(main=format_u("overall"),class=format_u(paste0("class",1:4)),
    note=paste("Entries are five-fold means (sample SDs), computed on the same 80 masked training profiles per fold.",
      "Bias, RMSE, MAE and width are in micromol/kg; coverage is percent for individual-U 95% posterior intervals.",
      "Bold marks the smaller absolute mean bias, RMSE and MAE across calibration levels for the same target group, using unrounded values. Coverage and width are descriptive."))
}

# Descriptive data analysis: measured nearby responses provide data context.
ocean_surface_report <- function(output_dir,data_dir,diagnostics_dir=file.path(output_dir,"diagnostics"),
  representative_fold=1L,folds=1:5) {
  inputs <- ocean_diag_inputs(data_dir); grids <- ocean_diag_grids(inputs$cohort)
  m <- read.csv(file.path(diagnostics_dir,"tables/m_points.csv"))
  f <- read.csv(file.path(diagnostics_dir,"tables/f_points.csv"))
  baseline_all <- ocean_baseline_mean_grid(output_dir,grids$m,folds)
  baseline <- baseline_all[baseline_all$design=="candidate_A",]
  baseline_same <- baseline_all[baseline_all$design=="same_x",]
  selected_m <- ocean_m_display_points(inputs$cohort,grids)
  points_m <- points_f <- list()
  for(fold in folds) {
    training <- inputs$train[[fold]]$cohort_id
    dm <- selected_m; dm$role <- ifelse(dm$cohort_id %in% training,"Training","Held-out")
    dm$fold <- fold; points_m[[length(points_m)+1L]] <- dm
    df <- grids$near[grids$near$U>=min(grids$f$U) & grids$near$U<=max(grids$f$U),]
    df$role <- ifelse(df$cohort_id %in% training,"Training","Held-out")
    df$fold <- fold; points_f[[length(points_f)+1L]] <- df
  }
  points_m <- ocean_bind(points_m); points_f <- ocean_bind(points_f)
  em <- m[m$design=="candidate_A" & m$fold %in% folds,]; em$model <- paste0("EIV-GP ",em$cell)
  lines <- ocean_bind(list(em,baseline))
  palette <- c("EIV-GP O0"="#777777","EIV-GP O20"="#2365a3","EIV-GP O50"="#c76326",
    "UC-GP (x,c)"="#7b3294","LVGP (x,c)"="#008837","EzGP (x,c)"="#a6761d")
  plot_m <- function(all=FALSE) {
    dl <- lines; de <- em; dp <- points_m
    if(!all) {dl<-dl[dl$fold==representative_fold,];de<-de[de$fold==representative_fold,];dp<-dp[dp$fold==representative_fold,]}
    dl$class_label <- paste("Class",dl$C); de$class_label <- paste("Class",de$C); dp$class_label <- paste("Class",dp$C)
    p <- ggplot2::ggplot(dl,ggplot2::aes(exp(log_pressure),posterior_mean,color=model)) +
      ggplot2::geom_ribbon(data=de,ggplot2::aes(ymin=q025,ymax=q975,fill=model),alpha=.08,color=NA) +
      ggplot2::geom_point(data=dp,ggplot2::aes(G2pressure,Y,shape=role),inherit.aes=FALSE,
        color="grey45",size=2,alpha=.7,stroke=.65) +
      ggplot2::geom_line(linewidth=.8) +
      ggplot2::scale_color_manual(values=palette,name=NULL) +
      ggplot2::scale_fill_manual(values=palette,guide="none") +
      ggplot2::scale_shape_manual(values=c(Training=16,"Held-out"=1),name="Nearby profiles") +
      ggplot2::guides(color=ggplot2::guide_legend(nrow=2,order=2),shape=ggplot2::guide_legend(order=1)) +
      ggplot2::scale_x_continuous(expand=ggplot2::expansion(mult=.05)) +
      ggplot2::labs(x="Pressure (dbar)",y="Silicate (micromol/kg)",
        title="Category response: fitted means and nearby observations",
        subtitle="Fixed salinity 34.20; separate oxygen classes",
        caption="Points: up to 3 nearest-salinity profiles per pressure bin, within the displayed curve range.\nEIV shading: pointwise 95% mean credible intervals. Baseline lines: original (x,c) fitted means.") +
      ggplot2::theme_minimal(base_size=12) + ggplot2::theme(legend.position="bottom",
        panel.grid.minor=ggplot2::element_blank(),plot.caption=ggplot2::element_text(hjust=0,size=10))
    if(all) p+ggplot2::facet_wrap(~fold+class_label,ncol=4,scales="free") else p+ggplot2::facet_wrap(~class_label,ncol=2,scales="free")
  }
  plot_f <- function(all=FALSE) {
    df <- f[f$fold %in% folds,]; dp <- points_f
    if(!all) {df<-df[df$fold==representative_fold,];dp<-dp[dp$fold==representative_fold,]}
    p <- ggplot2::ggplot(df,ggplot2::aes(U,posterior_mean,color=cell,fill=cell)) +
      ggplot2::geom_vline(xintercept=c(50,100,150),linetype=3,color="grey60") +
      ggplot2::geom_ribbon(ggplot2::aes(ymin=q025,ymax=q975),alpha=.12,color=NA) +
      ggplot2::geom_line(linewidth=.7) +
      ggplot2::geom_point(data=dp,ggplot2::aes(U,Y,shape=role),inherit.aes=FALSE,
        color="grey35",size=1.8,alpha=.7) +
      ggplot2::scale_color_manual(values=c(O20="#2365a3",O50="#c76326"),name="Training calibration") +
      ggplot2::scale_fill_manual(values=c(O20="#2365a3",O50="#c76326"),guide="none") +
      ggplot2::scale_shape_manual(values=c(Training=16,"Held-out"=1),name="Nearby observed profiles") +
      ggplot2::labs(x="Oxygen (micromol/kg)",y="Silicate (micromol/kg)",
        title="Continuous response f(x,u): salinity 34.25, pressure 650 dbar",
        caption="Points: nearby profiles with measured oxygen and silicate; their salinity/pressure vary.\nShading: pointwise 95% credible intervals for the mean, excluding response noise.") +
      ggplot2::theme_minimal(base_size=11) + ggplot2::theme(legend.position="bottom")
    if(all) p+ggplot2::facet_wrap(~fold,ncol=2) else p
  }
  plots <- list(m_comparison=plot_m(),f_with_observations=plot_f(),
    m_comparison_all_folds=plot_m(TRUE),f_with_observations_all_folds=plot_f(TRUE))
  dest <- file.path(diagnostics_dir,"surface_comparison")
  ocean_atomic_csv(baseline,file.path(dest,"baseline_m_points.csv"))
  ocean_atomic_csv(baseline_same,file.path(dest,"baseline_m_same_x.csv"))
  ocean_atomic_csv(points_m,file.path(dest,"nearby_m_observations.csv"))
  ocean_atomic_csv(selected_m,file.path(dest,"m_selection_audit.csv"))
  ocean_atomic_csv(points_f,file.path(dest,"nearby_f_observations.csv"))
  for(n in names(plots)) for(ext in c("png","pdf")) ggplot2::ggsave(file.path(dest,paste0(n,".",ext)),plots[[n]],
    width=if(grepl("fm_",n)) 11 else if(grepl("m_comparison",n)) 11 else 10,
    height=if(grepl("all_folds",n)) 13 else if(grepl("fm_",n)) 4.9 else 7,dpi=180,bg="white")
  list(plots=plots,baseline=baseline,baseline_same_x=baseline_same,points_m=points_m,points_f=points_f)
}

# Unified postprocessing entry used by both Rmd and this script's CLI.
# No fitting routine is called. Optional copied caches retain their full identity
# checks inside ocean_build_diagnostics; changed fits/settings are rejected.
ocean_run_diagnostics <- function(output_dir,data_dir,
  report_dir=file.path(output_dir,"diagnostics"),augmented_output_dir="",
  representative_fold=1L,workers=1L,cache_from="",
  build_prediction=TRUE,build_surfaces=TRUE,prediction_report=NULL) {
  for(p in c("eivGP","ggplot2","patchwork","digest"))
    ocean_diag_require(requireNamespace(p,quietly=TRUE),paste("Required package:",p))
  output_dir <- normalizePath(output_dir,mustWork=TRUE)
  data_dir <- normalizePath(data_dir,mustWork=TRUE)
  dir.create(report_dir,recursive=TRUE,showWarnings=FALSE)
  report_dir <- normalizePath(report_dir,mustWork=TRUE)
  progress <- function(stage) ocean_atomic_lines(c(paste("stage:",stage),
    paste("updated:",format(Sys.time(),"%Y-%m-%d %H:%M:%S %z")),
    paste("pid:",Sys.getpid()),"Postprocessing saved fits only; no MCMC or baseline refit."),
    file.path(report_dir,"PROGRESS.txt"))
  progress("input and source audit")
  inputs <- ocean_prepare_inputs(data_dir)
  sources <- ocean_diag_sources(output_dir)
  copied <- character()
  if(nzchar(cache_from)) {
    from <- normalizePath(cache_from,mustWork=TRUE)
    dir.create(file.path(report_dir,"cache"),showWarnings=FALSE)
    for(n in sprintf("fold%02d_%s.rds",sources$fold,sources$cell)) {
      source <- file.path(from,n); target <- file.path(report_dir,"cache",n)
      if(file.exists(source) && !file.exists(target)) {
        stopifnot(file.copy(source,target)); copied <- c(copied,n)
      }
    }
  }
  progress("U imputation and cached package m/f calculations")
  diagnostics <- ocean_build_diagnostics(output_dir,data_dir,report_dir=report_dir,
    representative_fold=representative_fold,workers=workers)
  prediction <- prediction_report; surfaces <- NULL
  if(build_prediction) {
    progress("original and O20/O50 expanded baseline prediction summaries")
    prediction <- ocean_prediction_report(output_dir,inputs,augmented_output_dir=augmented_output_dir,
      report_dir=file.path(report_dir,"prediction"))
  }
  if(!build_prediction && !is.null(prediction))
    ocean_export_prediction_report(prediction,file.path(report_dir,"prediction"))
  if(build_surfaces) {
    progress("native baseline m comparison and nearby measured-profile overlays")
    surfaces <- ocean_surface_report(output_dir,data_dir,report_dir,
      representative_fold=representative_fold)
  }
  # Paper-facing outputs use the selected b=0 encoding only. Raw audit
  # exports retain scores for every available encoding.
  main <- c("U imputation (O0/O20/O50, hollow calibration)"="figures/u_main.png",
    "m: native baselines and pressure-balanced nearby profiles"="surface_comparison/m_comparison.png",
    "f: S=34.25, P=650, measured nearby profiles"="surface_comparison/f_with_observations.png")
  main <- main[file.exists(file.path(report_dir,main))]
  ut <- ocean_paper_u_tables(diagnostics$u$summary)
  yt <- if(!is.null(prediction)) ocean_paper_prediction_tables(prediction) else NULL
  unote <- ut$note
  ynote <- if(!is.null(prediction)) paste0("Entries are fold means (sample SDs). sY = ",sprintf("%.4f",prediction$SD_Y),
    " micromol/kg. Coverage is percent. * marks baselines receiving the same available exact oxygen measurements using (x,c,R,W) inputs and b=0.") else ""
  text <- c("# Ocean diagnostics review","","Saved final fits; matched-oxygen baselines are marked with *.","",
    "## Paper table: U imputation","",ocean_table_html(ut$main,"Training oxygen imputation",unote),"",
    "[CSV](tables/paper_u_imputation_main.csv) | [LaTeX](tables/paper_u_imputation_main.tex)","")
  if(!is.null(yt) && nrow(yt$main)) text <- c(text,"## Paper table: Y prediction","",
    ocean_prediction_panel_html(yt$main,"Silicate prediction under matched oxygen information",
      paste(ynote,"Bold marks the smallest RMSE and CRPS mean within each complete five-fold information block, using unrounded values. Coverage and width are descriptive.")),"",
    "[CSV](prediction/paper_y_prediction_main.csv) | [LaTeX](prediction/paper_y_prediction_main.tex)","")
  text <- c(text,"## Main displays","",unlist(Map(function(label,path) c(paste0("### ",label),"",
    paste0("![",label,"](",path,")"),""),names(main),unname(main))),
    "## Additional tables","",
    "- [Class-specific U imputation](tables/paper_u_imputation_class.csv)",
    "- [Full Y-prediction summaries, b=0](prediction/paper_y_prediction_full.csv)",
    "- [Unrounded primary Y summary](prediction/primary_summary_numeric.csv)",
    "- [Same-X native baseline means](tables/baseline_m_same_x.csv)",
    "- [m measured-point selection audit](surface_comparison/m_selection_audit.csv)","",
    "Nearby points are measured responses at varying X. Curves and fitted posteriors are unchanged by the display-only point selection.",
    "Baseline category means are original (x,c) fits. Expanded baselines appear in matched Y-prediction comparisons.")
  ocean_atomic_lines(text,file.path(report_dir,"INDEX.md"))
  escape <- function(x) {x<-gsub("&","&amp;",x,fixed=TRUE);x<-gsub("<","&lt;",x,fixed=TRUE);gsub('"',"&quot;",x,fixed=TRUE)}
  html <- c('<!doctype html><meta charset="utf-8"><title>Ocean diagnostics review</title>',
    '<style>body{max-width:1150px;margin:32px auto;padding:0 16px;font:16px/1.5 sans-serif;color:#222}img{max-width:100%}section{margin:32px 0}table{border-collapse:collapse;width:100%;font-size:14px}td,th{padding:7px 9px;border-bottom:1px solid #ddd;text-align:right}td:nth-child(-n+3),th:nth-child(-n+3){text-align:left}caption{text-align:left;font-weight:600;padding:10px 0}.note{font-size:13px;color:#555}.table-scroll{overflow-x:auto}a{color:#2365a3}</style>',
    '<h1>Ocean diagnostics review</h1><p>Saved final fits; matched-oxygen baselines are marked with *. Tables show fold means (sample SDs).</p>',
    '<section><h2>Paper table: oxygen imputation</h2><div class="table-scroll">',
    ocean_table_html(ut$main,"Common 80 masked training profiles per fold",unote),
    '</div><p><a href="tables/paper_u_imputation_main.csv">CSV</a> · <a href="tables/paper_u_imputation_main.tex">LaTeX</a> · <a href="tables/paper_u_imputation_class.html">Class-specific table</a></p></section>')
  if(!is.null(yt) && nrow(yt$main)) html <- c(html,
    '<section><h2>Paper table: silicate prediction</h2><div class="table-scroll">',
    ocean_table_html(yt$main,"Original O0 comparators and matched-oxygen O20/O50 baselines",ynote),
    '</div><p><a href="prediction/paper_y_prediction_main.csv">CSV</a> · <a href="prediction/paper_y_prediction_main.tex">LaTeX</a> · <a href="prediction/paper_y_prediction_full.html">Full metric table</a></p></section>')
  html <- c(html,unlist(Map(function(label,path) paste0('<section><h2>',escape(label),'</h2><img src="',escape(path),'"></section>'),names(main),unname(main))))
  ocean_atomic_lines(html,file.path(report_dir,"preview.html"))
  saveRDS(list(diagnostic_plots=diagnostics$plots,
    surface_plots=if(!is.null(surfaces)) surfaces$plots else NULL,
    prediction_plots=if(!is.null(prediction)) prediction$plots else NULL,
    u_tables=ut,y_tables=yt),file.path(report_dir,"internal/presentation_plots.rds"))
  saveRDS(list(source_fit_root=output_dir,augmented_output_root=augmented_output_dir,
    report_dir=report_dir,copied_cache_files=copied,representative_fold=representative_fold,
    time=as.character(Sys.time()),no_refit=TRUE),file.path(report_dir,"internal/review_run.rds"))
  progress("COMPLETE")
  invisible(list(diagnostics=diagnostics,prediction=prediction,surfaces=surfaces,report_dir=report_dir))
}

# Direct execution: Rscript .../04_diagnostics.R diagnostics_settings.rds
# Source this file from Rmd/R without executing the CLI below.
.ocean_file_arg <- grep("^--file=",commandArgs(FALSE),value=TRUE)
if(sys.nframe()==0L && length(.ocean_file_arg) &&
   basename(sub("^--file=","",.ocean_file_arg[[1]]))=="04_diagnostics.R") {
  entry <- normalizePath(sub("^--file=","",.ocean_file_arg[[1]]),mustWork=TRUE)
  repository_dir <- normalizePath(file.path(dirname(entry),"../../.."),mustWork=TRUE)
  source(file.path(dirname(entry),"01_setup.R"))
  ocean_source_modules(repository_dir,globalenv())
  args <- commandArgs(trailingOnly=TRUE)
  if(!length(args)) stop("Provide a parameter RDS with output_dir (saved fit root), and optional report_dir, augmented_output_dir, cache_from.")
  p <- readRDS(args[[1]])
  if(is.null(p$output_dir) || !nzchar(p$output_dir)) stop("output_dir must identify existing completed fits.")
  paths <- c("output_dir","data_dir","report_dir","augmented_output_dir","cache_from")
  for(n in intersect(paths,names(p))) if(nzchar(p[[n]])) p[[n]] <- ocean_resolve_path(p[[n]],repository_dir)
  if(is.null(p$data_dir)) p$data_dir <- file.path(repository_dir,"data/ocean_silicate")
  p <- p[intersect(names(p),names(formals(ocean_run_diagnostics)))]
  if(is.null(p$report_dir)) p$report_dir <- file.path(p$output_dir,"diagnostics")
  dir.create(p$report_dir,recursive=TRUE,showWarnings=FALSE)
  result <- ocean_logged(file.path(p$report_dir,"execution.log"),do.call(ocean_run_diagnostics,p))
  cat("Diagnostics ready:",result$report_dir,"\n")
}

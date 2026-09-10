# Reporting and comparison helpers. Inference remains in eivGP 0.3.1.
ADNI_CASE_SCHEMA <- "case-study-v3"
ADNI_METHODS <- c("EIV-GP", "UC-GP", "LVGP", "EzGP")

adni_write_csv <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  tmp <- tempfile(tmpdir = dirname(path))
  on.exit(unlink(tmp), add = TRUE)
  write.csv(x, tmp, row.names = FALSE)
  if (!file.rename(tmp, path)) stop("Could not write ", path)
}

adni_case_preflight <- function(smoke = FALSE) {
  pkgs <- c("ggplot2", "patchwork", "kergp", "LVGP", "EzGP")
  tab <- data.frame(package = pkgs,
    installed = vapply(pkgs, requireNamespace, logical(1), quietly = TRUE))
  if (any(!tab$installed)) stop("Install case-study dependencies via setup.R: ",
                               paste(tab$package[!tab$installed], collapse = ", "))
  tab$version <- vapply(pkgs, function(p) as.character(utils::packageVersion(p)), character(1))
  tab
}

adni_eda <- function(dat, root, smoke) {
  z <- dat
  z$CSF <- factor(ifelse(z$R == 1, "Observed CSF", "Missing CSF"),
                  levels = c("Observed CSF", "Missing CSF"))
  z$pattern <- interaction(z$c_ab42_ab40_3, z$c_gfap_3, sep = "/", lex.order = TRUE)
  p1 <- ggplot2::ggplot(z, ggplot2::aes(x = CSF, fill = CSF)) +
    ggplot2::geom_bar(width = .65) + ggplot2::labs(x = NULL, y = "Participants", title = "A  CSF availability") +
    ggplot2::theme(legend.position = "none")
  p2 <- ggplot2::ggplot(z[z$R == 1, ], ggplot2::aes(x = u_csf_abeta42, y = y_centiloid)) +
    ggplot2::geom_point(alpha = .55, size = 1.3, color = "#0072B2") +
    ggplot2::labs(x = "Observed CSF A-beta42 (source units)", y = "Amyloid PET (Centiloid)",
                  title = "B  Observed-CSF subset")
  p3 <- ggplot2::ggplot(z, ggplot2::aes(x = pattern, y = y_centiloid)) +
    ggplot2::geom_boxplot(outlier.size = .7, fill = "#B8D8E8") +
    ggplot2::labs(x = "Ordinal A-beta42/40 / GFAP pattern", y = "Amyloid PET (Centiloid)",
                  title = "C  Variation within proxy patterns")
  fig <- patchwork::wrap_plots(p1, p2, p3, nrow = 1) +
    patchwork::plot_annotation(title = if (smoke) "SMOKE RUN - descriptive cohort figure" else
      "ADNI: incomplete CSF measurements and ordinal plasma proxies",
      caption = "Exploratory associations, not causal effects. Proxy classes are constructed from continuous assays.")
  ggplot2::ggsave(file.path(root, "main", "fig1_exploration.pdf"), fig, width = 12, height = 3.8)
  ggplot2::ggsave(file.path(root, "main", "fig1_exploration.png"), fig, width = 12, height = 3.8, dpi = 150)
  counts <- as.data.frame(table(z$pattern, z$CSF))
  names(counts) <- c("pattern", "CSF_availability", "n")
  adni_write_csv(counts, file.path(root, "appendix", "pattern_counts.csv"))
  desc <- do.call(rbind, lapply(split(z, z$CSF), function(a) data.frame(
    group = as.character(a$CSF[1]), n = nrow(a), mean_age = mean(a$x_age_years),
    female_fraction = mean(a$x_female), mean_apoe4_dose = mean(a$x_apoe4_dose),
    mean_centiloid = mean(a$y_centiloid), sd_centiloid = sd(a$y_centiloid))))
  adni_write_csv(desc, file.path(root, "appendix", "cohort_by_CSF_availability.csv"))
}

# Empirical-distribution CRPS uses every draw, in O(S log S) per participant.
adni_draw_summary <- function(draws, truth, point_mean = NULL) {
  stopifnot(is.matrix(draws), ncol(draws) == length(truth), nrow(draws) > 1L,
            all(is.finite(draws)), all(is.finite(truth)))
  s <- nrow(draws)
  cm <- attr(draws, "conditional_means", exact = TRUE)
  cv <- attr(draws, "conditional_vars", exact = TRUE)
  if (is.null(point_mean) && is.matrix(cm) && ncol(cm) == length(truth) && all(is.finite(cm)))
    point_mean <- colMeans(cm)
  mu <- if (is.null(point_mean)) colMeans(draws) else as.numeric(point_mean)
  nlpd <- rep(NA_real_, length(truth))
  if (is.matrix(cm) && is.matrix(cv) && identical(dim(cm), dim(cv)) &&
      ncol(cm) == length(truth) && all(is.finite(cm)) && all(is.finite(cv)) && all(cv > 0)) {
    nlpd <- vapply(seq_along(truth), function(j) {
      ll <- dnorm(truth[j], cm[, j], sqrt(cv[, j]), log = TRUE)
      mx <- max(ll); -(mx + log(mean(exp(ll - mx))))
    }, numeric(1))
  }
  stopifnot(length(mu) == length(truth), all(is.finite(mu)))
  lo <- apply(draws, 2L, stats::quantile, probs = .025)
  hi <- apply(draws, 2L, stats::quantile, probs = .975)
  crps <- vapply(seq_along(truth), function(j) {
    a <- sort(draws[, j])
    mean(abs(a - truth[j])) - sum((2 * seq_len(s) - s - 1) * a) / s^2
  }, numeric(1))
  score <- hi - lo + 40 * pmax(lo - truth, 0) + 40 * pmax(truth - hi, 0)
  data.frame(y_true = truth, pred_mean = mu, pred_sd = apply(draws, 2L, stats::sd),
    pred_lo95 = lo, pred_hi95 = hi, squared_error = (truth - mu)^2,
    absolute_error = abs(truth - mu), CRPS = crps, NLPD = nlpd, IntervalScore95 = score,
    covered95 = truth >= lo & truth <= hi, Width95 = hi - lo,
    covered50 = truth >= apply(draws, 2L, quantile, .25) & truth <= apply(draws, 2L, quantile, .75),
    covered80 = truth >= apply(draws, 2L, quantile, .10) & truth <= apply(draws, 2L, quantile, .90))
}

adni_prediction_rows <- function(draws, test, method, scenario, repeat_id, fold,
                                 convergence = NA, point_mean = NULL) {
  z <- adni_draw_summary(draws, test$y_centiloid, point_mean)
  data.frame(repeat_id = repeat_id, fold = fold, RID = test$RID, R = test$R,
    diagnosis = test$diagnosis, method = method, scenario = scenario,
    convergence_passed = convergence, n_draws = nrow(draws), z)
}

adni_metrics <- function(z) {
  if (!nrow(z)) return(data.frame())
  do.call(rbind, lapply(c("overall", "R0", "R1"), function(g) {
    a <- if (g == "overall") z else z[z$R == as.integer(sub("R", "", g)), ]
    data.frame(group = g, n = nrow(a), RMSE = sqrt(mean(a$squared_error)),
      MAE = mean(a$absolute_error), CRPS = mean(a$CRPS),
      NLPD = if (all(is.finite(a$NLPD))) mean(a$NLPD) else NA_real_,
      Coverage95 = mean(a$covered95), Width95 = mean(a$Width95),
      IntervalScore95 = mean(a$IntervalScore95), Coverage50 = mean(a$covered50),
      Coverage80 = mean(a$covered80))
  }))
}

adni_competitors <- function(train, test, n_draw, seed, smoke, path) {
  xc <- c("x_age_years", "x_female", "x_apoe4_dose")
  cc <- c("c_ab42_ab40_3", "c_gfap_3")
  # All transformations use training rows only. Category codes are not scaled.
  center <- colMeans(train[, xc]); scale <- vapply(train[, xc], sd, numeric(1))
  scale[!is.finite(scale) | scale == 0] <- 1
  X <- sweep(sweep(as.matrix(train[, xc]), 2L, center), 2L, scale, "/")
  NX <- sweep(sweep(as.matrix(test[, xc]), 2L, center), 2L, scale, "/")
  yc <- mean(train$y_centiloid); ys <- sd(train$y_centiloid)
  stopifnot(is.finite(ys), ys > 0)
  key <- list(schema = ADNI_CASE_SCHEMA, train = train[, c("RID", xc, cc, "y_centiloid")],
    test = test[, c("RID", xc, cc)], n_draw = n_draw, seed = seed, smoke = smoke,
    eiv_code_md5 = unname(tools::md5sum(list.files(system.file("R", package = "eivGP"), full.names = TRUE))),
    packages = vapply(c("eivGP", "kergp", "LVGP", "EzGP"), function(p) as.character(packageVersion(p)), character(1)))
  if (file.exists(path)) {
    old <- readRDS(path)
    if (identical(old$key, key)) return(old$value)
    stop("Comparator cache configuration changed; remove or relocate: ", path)
  }
  out <- list(draws = list(), means = list(), status = list())
  controls <- if (smoke) list(
    `UC-GP` = list(n_starts = 1L),
    LVGP = list(n_starts = 1L, max_retries = 1L, max_iter_ini = 5L, max_iter_lat = 2L,
                rescue_iter_ini = 5L, rescue_iter_lat = 2L, max_elapsed_seconds = 20),
    EzGP = list(maxeval = 5L, tau_fractions = .01, cv_folds = 2L)
  ) else list()
  for (j in seq_along(ADNI_METHODS[-1])) {
    method <- ADNI_METHODS[-1][j]; started <- proc.time()[3]
    warnings <- character()
    cat("ADNI comparator:", method, "\n"); flush.console()
    ans <- tryCatch(withCallingHandlers(do.call(eivGP::fit_mixedgp_competitor,
      c(list(method = method, X = X, y = (train$y_centiloid - yc) / ys,
        C = as.matrix(train[, cc]), new_X = NX, new_C = as.matrix(test[, cc]),
        m_vec = c(3L, 3L), n_draw = as.integer(n_draw), seed = as.integer(seed + j)), controls[[method]])),
      warning = function(w) { warnings <<- c(warnings, conditionMessage(w)); invokeRestart("muffleWarning") }),
      error = function(e) e)
    good <- !inherits(ans, "error")
    if (good) {
      good <- is.matrix(ans$draws) && identical(dim(ans$draws), c(as.integer(n_draw), as.integer(nrow(test)))) &&
        all(is.finite(ans$draws)) && length(ans$predictive_mean) == nrow(test) && all(is.finite(ans$predictive_mean))
      if (!good) ans <- simpleError("Invalid public-package prediction dimensions or values.")
    }
    if (good) {
      d <- yc + ys * ans$draws
      cm <- attr(ans$draws, "conditional_means", exact = TRUE)
      cv <- attr(ans$draws, "conditional_vars", exact = TRUE)
      attr(d, "conditional_means") <- if (is.null(cm)) NULL else yc + ys * cm
      attr(d, "conditional_vars") <- if (is.null(cv)) NULL else ys^2 * cv
      out$draws[[method]] <- d
      out$means[[method]] <- yc + ys * ans$predictive_mean
    }
    out$status[[method]] <- data.frame(method = method, status = if (good) "success" else "failed",
      optimization_status = if (good && !is.null(ans$optimization_status)) ans$optimization_status else NA_character_,
      elapsed_seconds = unname(proc.time()[3] - started),
      message = if (good) "" else conditionMessage(ans), warnings = paste(unique(warnings), collapse = " | "))
  }
  out$status <- do.call(rbind, out$status)
  atomic_save_rds(list(key = key, value = out), path)
  out
}

adni_extended_fold <- function(fit, train, test, draw_ids, fold_dir, repeat_id,
                                fold, smoke, diagnostic_passed, available_draws) {
  out_dir <- file.path(fold_dir, "case-study")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  seed <- as.integer(202640000 + 100 * repeat_id + fold)
  cc <- c("c_ab42_ab40_3", "c_gfap_3")
  xc <- c("x_age_years", "x_female", "x_apoe4_dose")
  Cnew <- make_C(test)
  # Primary target: same test information X,C for every method.
  no_csf <- eivGP::predict_eivgp(fit, new_X = as.matrix(test[, xc]), new_C = Cnew,
    target = "response", draw_ids = draw_ids, seed = seed)
  # Check prediction sensitivity to the retained chain subset; these are not independent CV replications.
  chains <- fit$mcmc$mcmc_draw_info$chain[draw_ids]
  chain_scores <- do.call(rbind, lapply(sort(unique(chains)), function(k) {
    dd <- no_csf[chains == k, , drop = FALSE]
    if (nrow(dd) < 2L) return(NULL)
    cm <- attr(no_csf, "conditional_means", exact = TRUE)
    cv <- attr(no_csf, "conditional_vars", exact = TRUE)
    if (is.matrix(cm) && nrow(cm) == length(chains)) attr(dd, "conditional_means") <- cm[chains == k, , drop = FALSE]
    if (is.matrix(cv) && nrow(cv) == length(chains)) attr(dd, "conditional_vars") <- cv[chains == k, , drop = FALSE]
    zz <- adni_prediction_rows(dd, test, "EIV-GP", "no_test_CSF", repeat_id, fold, diagnostic_passed)
    mm <- adni_metrics(zz); mm$chain <- k; mm
  }))
  if (nrow(chain_scores)) adni_write_csv(chain_scores, file.path(out_dir, "prediction_scores_by_chain.csv"))
  rows <- list(adni_prediction_rows(no_csf, test, "EIV-GP", "no_test_CSF", repeat_id,
    fold, diagnostic_passed), adni_prediction_rows(available_draws, test, "EIV-GP", "available_CSF",
    repeat_id, fold, diagnostic_passed))
  comparison <- adni_competitors(train, test, nrow(no_csf), seed + 1000L, smoke,
                                file.path(out_dir, "competitors.rds"))
  for (method in names(comparison$draws)) for (scenario in c("no_test_CSF", "available_CSF")) {
    rows[[length(rows) + 1L]] <- adni_prediction_rows(comparison$draws[[method]], test, method,
      scenario, repeat_id, fold, point_mean = comparison$means[[method]])
  }
  adni_write_csv(do.call(rbind, rows), file.path(out_dir, "predictions.csv"))
  status <- rbind(data.frame(method = "EIV-GP", status = "success",
    optimization_status = if (diagnostic_passed) "mcmc_gate_passed" else "mcmc_gate_failed",
    elapsed_seconds = NA_real_, message = "Sampling time in PROGRESS.md", warnings = ""), comparison$status)
  status$repeat_id <- repeat_id; status$fold <- fold
  adni_write_csv(status, file.path(out_dir, "method_status.csv"))
  # Save only the bounded predictive draw subset, not another full MCMC fit.
  atomic_save_rds(list(no_test_CSF = no_csf, available_CSF = available_draws,
    draw_ids = draw_ids, RID = test$RID), file.path(out_dir, "eiv_predictive_draws.rds"))

  observed <- which(test$R == 1L)
  if (length(observed)) {
    # Prospective U|C; neither held-out U nor held-out Y is passed to imputation.
    ud <- eivGP::impute_eivgp(fit, new_C = Cnew[observed, , drop = FALSE],
      draw_ids = draw_ids, scale = "raw", seed = seed + 2000L)
    ud <- matrix(ud[, , 1L], nrow = length(draw_ids), ncol = length(observed))
    scores <- adni_draw_summary(ud, test$u_csf_abeta42[observed])
    scores$probability_nonpositive_CSF <- colMeans(ud <= 0)
    urows <- data.frame(repeat_id = repeat_id, fold = fold, RID = test$RID[observed],
      method = "EIV-GP: prospective U|C", convergence_passed = diagnostic_passed, scores)
    # Simple response-free benchmark illustrates that imputation is not exclusive to EIV-GP.
    tr <- train[train$R == 1L, ]
    for (v in cc) { tr[[v]] <- factor(tr[[v]], levels = 1:3) }
    te <- test[observed, ]
    for (v in cc) te[[v]] <- factor(te[[v]], levels = 1:3)
    ub <- tryCatch({
      model <- lm(u_csf_abeta42 ~ c_ab42_ab40_3 + c_gfap_3, data = tr)
      pr <- predict(model, te, se.fit = TRUE)
      set.seed(seed + 2001L)
      v <- matrix(rt(length(draw_ids) * nrow(te), df = model$df.residual), nrow = length(draw_ids))
      v <- sweep(sweep(v, 2, sqrt(pr$se.fit^2 + pr$residual.scale^2), "*"), 2, pr$fit, "+")
      result <- adni_draw_summary(v, test$u_csf_abeta42[observed], pr$fit)
      result$probability_nonpositive_CSF <- colMeans(v <= 0)
      result
    }, error = function(e) e)
    if (!inherits(ub, "error")) urows <- rbind(urows, data.frame(repeat_id = repeat_id, fold = fold,
      RID = test$RID[observed], method = "LM-CSF: response-free", convergence_passed = NA, ub))
    writeLines(if (inherits(ub, "error")) conditionMessage(ub) else "success",
               file.path(out_dir, "CSF_benchmark_status.txt"))
    adni_write_csv(urows, file.path(out_dir, "hidden_CSF_validation.csv"))
  }
  # Draws for missing training U are conditional on observed training Y; no truth score exists.
  missing <- which(train$R == 0L)
  ud <- eivGP::impute_eivgp(fit, rows = missing, draw_ids = draw_ids, scale = "raw")
  umat <- matrix(ud[, , 1L], nrow = length(draw_ids), ncol = length(missing))
  adni_write_csv(data.frame(RID = train$RID[missing],
    mean = colMeans(umat), probability_nonpositive_CSF = colMeans(umat <= 0), lo95 = apply(umat, 2, quantile, .025),
    hi95 = apply(umat, 2, quantile, .975)), file.path(out_dir, "training_missing_CSF_posterior.csv"))

  # Fixed, prespecified fold 1 for a physical-coordinate surface, never best-fold selection.
  if (fold == 1L) {
    ug <- seq(quantile(train$u_csf_abeta42[train$R == 1], .1),
              quantile(train$u_csf_abeta42[train$R == 1], .9), length.out = if (smoke) 8L else 40L)
    reference <- c(x_age_years = median(train$x_age_years), x_female = 0, x_apoe4_dose = 1)
    xx <- matrix(rep(reference, each = length(ug)), ncol = 3,
                 dimnames = list(NULL, names(reference)))
    fs <- eivGP::predict_eivgp(fit, new_X = xx, new_U = matrix(ug, ncol = 1),
      new_U_scale = "raw", target = "surface", draw_ids = draw_ids, seed = seed + 3000L)
    adni_write_csv(data.frame(U = ug, mean = colMeans(fs), lo95 = apply(fs, 2, quantile, .025),
      hi95 = apply(fs, 2, quantile, .975), age = unname(reference[1]), female = 0, apoe4_dose = 1,
      fold = fold, repeat_id = repeat_id, convergence_passed = diagnostic_passed),
      file.path(out_dir, "CSF_PET_surface.csv"))
  }
  # Include the dictionary states deliberately omitted from the lightweight stopping gate.
  adni_write_csv(fit$diagnostics$rhat_hyper, file.path(out_dir, "all_hyper_diagnostics.csv"))
  j <- fit$mcmc$samples_by_chain$J
  if (length(j)) {
    occupancy <- do.call(rbind, lapply(seq_along(j), function(k) do.call(rbind,
      lapply(seq_len(ncol(j[[k]])), function(d) {
        tb <- table(j[[k]][, d]); data.frame(chain = k, coordinate = d,
          index = as.integer(names(tb)), fraction = as.numeric(tb) / sum(tb))
      }))))
    adni_write_csv(occupancy, file.path(out_dir, "dictionary_occupancy.csv"))
  }
  invisible(NULL)
}

adni_case_report <- function(root, folds_expected = 1:3, smoke = FALSE) {
  case_dirs <- file.path(root, paste0("fold_", folds_expected), "case-study")
  files <- file.path(case_dirs, "predictions.csv")
  if (!all(file.exists(files))) stop("Missing fold prediction files; reporting requires every requested fold.")
  z <- do.call(rbind, lapply(files, read.csv, stringsAsFactors = FALSE))
  status <- do.call(rbind, lapply(file.path(case_dirs, "method_status.csv"), read.csv,
                                stringsAsFactors = FALSE))
  adni_write_csv(status, file.path(root, "appendix", "method_status.csv"))
  stopifnot(!anyDuplicated(z[, c("repeat_id", "RID", "method", "scenario")]))
  reference <- z[z$method == "EIV-GP" & z$scenario == "no_test_CSF", c("RID", "fold")]
  for (m in ADNI_METHODS[-1]) {
    a <- z[z$method == m & z$scenario == "no_test_CSF", c("RID", "fold")]
    for (f in unique(a$fold)) if (!setequal(a$RID[a$fold == f], reference$RID[reference$fold == f]))
      stop("Participant mismatch for ", m, " fold ", f)
  }
  coverage <- do.call(rbind, lapply(ADNI_METHODS, function(m) {
    a <- z[z$method == m & z$scenario == "no_test_CSF", ]
    data.frame(method = m, n_folds = length(unique(a$fold)), n_predictions = nrow(a),
      all_requested_folds = setequal(unique(a$fold), folds_expected),
      complete_repeat = !smoke && nrow(a) == 495L && setequal(unique(a$fold), 1:3),
      flagged_folds = sum(status$method == m & (status$optimization_status != "converged" &
        status$optimization_status != "mcmc_gate_passed"), na.rm = TRUE))
  }))
  adni_write_csv(coverage, file.path(root, "appendix", "completeness.csv"))
  # Main table requires all methods on the same fold set. Never compare each method's favorable subset.
  common_folds <- Reduce(intersect, lapply(ADNI_METHODS, function(m)
    unique(z$fold[z$method == m & z$scenario == "no_test_CSF"])))
  main <- list(); detailed <- list(); fold_metrics <- list()
  for (m in ADNI_METHODS) for (scenario in unique(z$scenario)) {
    a <- z[z$method == m & z$scenario == scenario, ]
    if (!nrow(a)) next
    metric <- adni_metrics(a); metric$method <- m; metric$scenario <- scenario
    metric$complete_repeat <- coverage$complete_repeat[coverage$method == m]
    detailed[[length(detailed) + 1L]] <- metric
    for (f in unique(a$fold)) {
      fm <- adni_metrics(a[a$fold == f, ]); fm$method <- m; fm$scenario <- scenario; fm$fold <- f
      fold_metrics[[length(fold_metrics) + 1L]] <- fm
    }
    if (scenario == "no_test_CSF") {
      matched <- a[a$fold %in% common_folds, ]
      if (nrow(matched)) {
        mm <- adni_metrics(matched)
        mm <- mm[mm$group %in% c("overall", "R0"), ]
        mm$method <- m; mm$matched_folds <- length(common_folds)
        mm$complete_comparison <- !smoke && all(coverage$complete_repeat) && setequal(common_folds, 1:3)
        mm$diagnostic_flag <- coverage$flagged_folds[coverage$method == m] > 0
        main[[length(main) + 1L]] <- mm
      }
    }
  }
  adni_write_csv(do.call(rbind, detailed), file.path(root, "appendix", "all_scenario_metrics.csv"))
  adni_write_csv(do.call(rbind, fold_metrics), file.path(root, "appendix", "fold_metrics.csv"))
  main_table <- if (length(main)) do.call(rbind, main) else data.frame(
    status = "No fold has valid predictions from all four methods; no comparative ranking available.")
  if ("method" %in% names(main_table)) main_table <- main_table[, c("method", "group", "n",
    "RMSE", "CRPS", "Coverage95", "Width95", "IntervalScore95", "complete_comparison", "diagnostic_flag")]
  adni_write_csv(main_table, file.path(root, "main", "table1_prediction.csv"))
  writeLines(c("# ADNI prediction comparison", "",
    if (smoke) "SMOKE OUTPUT - no inferential interpretation." else
      if (!all(coverage$complete_repeat) || !setequal(common_folds, 1:3)) "INCOMPLETE: descriptive common-fold subset only; no full-repeat superiority claim." else
        "All methods evaluated on identical out-of-fold participants.", "",
    "Primary task: PET response prediction with no test CSF. Overall scores and naturally missing-CSF subgroup are reported.",
    "EIV-GP additionally uses calibrated training CSF; competitors use training X,C,Y. This compares available workflows, not an isolated algorithmic effect.",
    "Diagnostic warnings are retained. Coverage and width are not ranked. See appendix for all failures and settings.", "",
    as.character(knitr::kable(main_table, format = "pipe", digits = 3))),
    file.path(root, "main", "table1_prediction.md"))

  # Per-participant paired contrasts. Bootstrap conditions on fitted folds, not independent refitted studies.
  paired <- list(); set.seed(202650000L)
  for (m in ADNI_METHODS[-1]) for (group in c("overall", "R0")) {
    e <- z[z$method == "EIV-GP" & z$scenario == "no_test_CSF", ]
    a <- z[z$method == m & z$scenario == "no_test_CSF", ]
    if (group == "R0") { e <- e[e$R == 0, ]; a <- a[a$R == 0, ] }
    merge_cols <- c("repeat_id", "fold", "RID")
    p <- merge(e[, c(merge_cols, "CRPS", "absolute_error")],
               a[, c(merge_cols, "CRPS", "absolute_error")], by = merge_cols, suffixes = c("_eiv", "_other"))
    if (!nrow(p)) next
    for (loss in c("CRPS", "absolute_error")) {
      difference <- p[[paste0(loss, "_other")]] - p[[paste0(loss, "_eiv")]]
      bs <- replicate(if (smoke) 50L else 1000L, mean(sample(difference, replace = TRUE)))
      paired[[length(paired) + 1L]] <- data.frame(method = m, group = group, loss = loss,
        n_pairs = length(difference), competitor_minus_EIV = mean(difference),
        lo_descriptive = quantile(bs, .025), hi_descriptive = quantile(bs, .975))
    }
  }
  if (length(paired)) adni_write_csv(do.call(rbind, paired), file.path(root, "appendix", "paired_losses.csv"))
  # Observable task shared by every method: held-out response predictions.
  zp <- z[z$scenario == "no_test_CSF", ]; zp$method <- factor(zp$method, levels = ADNI_METHODS)
  p <- ggplot2::ggplot(zp, ggplot2::aes(x = y_true, y = pred_mean)) +
    ggplot2::geom_point(alpha = .45, size = .8) + ggplot2::geom_abline(slope = 1, intercept = 0, linetype = 2) +
    ggplot2::facet_wrap(~method, nrow = 1) + ggplot2::theme_bw() +
    ggplot2::labs(x = "Observed Centiloid", y = "Out-of-fold prediction", title = "Common response prediction task")
  ggplot2::ggsave(file.path(root, "appendix", "prediction_scatter.pdf"), p, width = 11, height = 3.3)
  uf <- file.path(case_dirs, "hidden_CSF_validation.csv"); uf <- uf[file.exists(uf)]
  if (length(uf)) {
    uz <- do.call(rbind, lapply(uf, read.csv, stringsAsFactors = FALSE))
    adni_write_csv(uz, file.path(root, "appendix", "hidden_CSF_predictions.csv"))
    um <- do.call(rbind, lapply(split(uz, uz$method), function(a) data.frame(method = a$method[1],
      n = nrow(a), RMSE = sqrt(mean(a$squared_error)), CRPS = mean(a$CRPS),
      Coverage95 = mean(a$covered95), Width95 = mean(a$Width95),
      mean_probability_nonpositive_CSF = mean(a$probability_nonpositive_CSF))))
    adni_write_csv(um, file.path(root, "appendix", "hidden_CSF_metrics.csv"))
    e <- uz[uz$method == "EIV-GP: prospective U|C", ]
    pi <- ggplot2::ggplot(e, ggplot2::aes(x = y_true, y = pred_mean)) +
      ggplot2::geom_linerange(ggplot2::aes(ymin = pred_lo95, ymax = pred_hi95), alpha = .17) +
      ggplot2::geom_point(size = 1.2, color = "#0072B2") + ggplot2::geom_abline(slope = 1, intercept = 0, linetype = 2) +
      ggplot2::labs(x = "Held-out observed CSF A-beta42", y = "Prospective U|C mean and 95% interval",
                    title = "A  Masked-CSF validation (observed subset)") + ggplot2::theme_bw()
    sf <- file.path(root, "fold_1", "case-study", "CSF_PET_surface.csv")
    if (file.exists(sf)) {
      s <- read.csv(sf)
      ps <- ggplot2::ggplot(s, ggplot2::aes(x = U, y = mean)) +
        ggplot2::geom_ribbon(ggplot2::aes(ymin = lo95, ymax = hi95), fill = "#0072B2", alpha = .2) +
        ggplot2::geom_line(color = "#0072B2") + ggplot2::theme_bw() +
        ggplot2::labs(x = "CSF A-beta42 (source units)", y = "Latent PET response surface",
          title = "B  Physical-coordinate inference",
          subtitle = paste0("Prespecified fold 1; age ", round(s$age[1], 1), ", female=0, APOE4 dose=1"))
      combined <- patchwork::wrap_plots(pi, ps) + patchwork::plot_annotation(
        title = if (smoke) "SMOKE OUTPUT - not inferential results" else "EIV-GP: additional inference on the CSF scale",
        caption = paste("Pointwise posterior intervals; inspect convergence diagnostics.",
          "Surface is an association, not a causal effect. Masked-CSF validation does not verify naturally missing CSF."))
      ggplot2::ggsave(file.path(root, "main", "fig2_CSF_inference.pdf"), combined, width = 10, height = 4.5)
      ggplot2::ggsave(file.path(root, "main", "fig2_CSF_inference.png"), combined, width = 10, height = 4.5, dpi = 150)
    }
  }
  writeLines(c("# Reporting boundaries", "",
    "Do not claim superiority from smoke, incomplete, or inadequately mixed fits.",
    "Positive paired differences favor EIV-GP. Bootstrap intervals are descriptive conditional-on-fits intervals; they omit refitting and model-selection uncertainty.",
    "Repeat estimates share participants and are not independent experiments. Do not concatenate repeats as independent observations.",
    "The literature competitors do not estimate physical CSF posteriors or a CSF-coordinate surface. Other calibrated models can; these capabilities are not universally unique to EIV-GP.",
    "Background and task specifications are in the notebook. EDA is descriptive and must not drive post-hoc test-set tuning."),
    file.path(root, "appendix", "INTERPRETATION.md"))
  invisible(main_table)
}

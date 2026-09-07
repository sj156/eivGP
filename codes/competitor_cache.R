## Repository-only competitor protocol/cache. No eivGP MCMC or installation.
mixedgp_competitor_protocol <- function(study) {
  study <- match.arg(study, c("study1", "study2"))
  list(`UC-GP` = list(n_starts = 8L),
    LVGP = list(n_starts = 8L, max_retries = 3L, max_iter_ini = 100L,
      max_iter_lat = 20L, rescue_iter_ini = 300L, rescue_iter_lat = 100L,
      max_elapsed_seconds = if (study == "study1") 1800 else 3600, parallel = FALSE),
    EzGP = list(tau_fractions = c(1e-6, .0025, .01, .04, .16), cv_folds = 3L, maxeval = 100L))
}

mixedgp_competitor_hash <- function(x) {
  f <- tempfile(); on.exit(unlink(f), add = TRUE)
  saveRDS(x, f, version = 3L, compress = FALSE)
  unname(tools::md5sum(f))
}

mixedgp_cached_competitors <- function(X_train, y_train, C_train, X_test, C_test,
    n_draw, seed, m_vec = NULL, methods = c("UC-GP", "LVGP", "EzGP"),
    strict = FALSE, controls = list(), cache_root = .mixedgp_competitor_cache_root,
    allow_fit = FALSE, retry_failed = FALSE) {
  checked <- validate_mixedgp_inputs(X_train, C_train, m_vec)
  test <- validate_mixedgp_inputs(X_test, C_test, checked$m_vec)
  n_draw <- mixedgp_as_integer_strict(n_draw, "n_draw", 1L, 1L)
  registry <- mixedgp_competitor_registry()
  if (anyDuplicated(methods) || any(!methods %in% registry$method)) stop("Invalid competitor selection.")
  out <- list(draws=list(), means=list(), latent_means=list(), predictive_means=list(),
              predictive_variances=list(), fits=list(), status=NULL)
  statuses <- list()
  for (method in methods) {
    j <- match(method, registry$method)
    package <- registry$package[j]
    available <- requireNamespace(package, quietly = TRUE)
    identity <- list(schema = "competitor-cache-v1", method = method,
      X = unname(checked$X), y = as.numeric(y_train), C = unname(checked$C),
      X_test = unname(test$X), C_test = unname(test$C), m_vec = checked$m_vec,
      seed = as.integer(seed + j), controls = controls[[method]],
      package_version = if (available) as.character(utils::packageVersion(package)) else "missing",
      R_version = paste(R.version$major, R.version$minor, R.version$platform),
      adapter_md5 = unname(tools::md5sum(file.path(.mixedgp_competitor_code_dir,
        "03_study2_published_competitors.R"))))
    key <- mixedgp_competitor_hash(identity)
    path <- file.path(cache_root, method, paste0(key, ".rds"))
    cached <- if (file.exists(path)) readRDS(path) else NULL
    if (!is.null(cached) && !identical(cached$identity, identity)) stop("Competitor cache identity mismatch: ", path)
    hit <- !is.null(cached)
    need_fit <- is.null(cached) || (retry_failed && !all(cached$result$status$status == "success"))
    if (need_fit && allow_fit && available) {
      message("Competitor ", method, ": fitting ", substr(key,1L,12L))
      dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
      lock <- paste0(path, ".lock")
      if (!dir.create(lock, showWarnings = FALSE)) stop("Competitor cache is locked: ", lock,
        ". Another run may be active; do not remove its lock while it is running.")
      result <- tryCatch({
        ## The adapter returns a Gaussian predictor. Retain the model and exact
        ## moments; draws can subsequently be generated at any evaluation budget.
        run_published_mixedgp_competitors(checked$X, as.numeric(y_train), checked$C,
          test$X, test$C, n_draw = 1L, seed = as.integer(seed + j - 1L),
          m_vec = checked$m_vec, methods = method, strict = FALSE, controls = controls)
      }, error = function(e) list(status = mixedgp_competitor_status(method,
        status="unavailable_or_failed", optimization_status="error", message=conditionMessage(e))))
      tryCatch({
        if (file.exists(path)) {
          backup <- tempfile(paste0(key,"-previous-"), tmpdir=dirname(path), fileext=".rds")
          if (!file.copy(path, backup)) stop("Cannot preserve previous attempt.")
        }
        cached <- list(identity=identity, result=result, created=as.character(Sys.time()))
        tmp <- tempfile(tmpdir=dirname(path)); saveRDS(cached,tmp,version=3L)
        if (!file.rename(tmp,path)) stop("Cannot save competitor cache: ",path)
      }, finally = unlink(lock, recursive=TRUE))
      hit <- FALSE
      message("Competitor ", method, ": saved ", paste(result$status$status,collapse=","))
    }
    if (is.null(cached)) {
      status <- mixedgp_competitor_status(method, status="unavailable_or_failed",
        optimization_status=if (available) "missing_cache" else "missing_package",
        message=if (available) "Run experiments/run_competitors.R to prepare this method; simulation did not fit it."
          else paste("Install experiment dependencies; missing package", package))
    } else {
      result <- cached$result; status <- result$status
      if (all(status$status == "success")) {
        mu <- result$predictive_means[[method]]
        variance <- result$predictive_variances[[method]]
        draw <- sample_independent_predictive_marginals(mu, variance, n_draw, as.integer(seed+j+10000L))
        out$draws[[method]] <- attach_predictive_normal_components(draw, mu, variance)
        out$means[[method]] <- out$latent_means[[method]] <- result$latent_means[[method]]
        out$predictive_means[[method]] <- mu
        out$predictive_variances[[method]] <- variance
        out$fits[[method]] <- result$fits[[method]]
      }
    }
    status$cache_key <- key; status$cache_hit <- hit
    status$cache_file <- if (file.exists(path)) normalizePath(path) else path
    statuses[[method]] <- status
    if (strict && !all(status$status == "success")) stop(status$message[1L])
  }
  out$status <- if (length(statuses)) do.call(rbind,statuses) else
    mixedgp_competitor_status("UC-GP",status="success")[FALSE,,drop=FALSE]
  out
}

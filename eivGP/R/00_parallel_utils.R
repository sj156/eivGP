############################################################
## Reproducible parallel utilities
##
## Package-ready wrappers around parallel::mclapply() and
## parallel::mcmapply(). Forking is used on macOS/Linux and the same call
## falls back to serial evaluation on Windows or when one core is requested.
############################################################

mixedgp_as_integer_strict <- function(x,
                                      name,
                                      min_value = -.Machine$integer.max,
                                      length_expected = NULL) {
  if (!(is.numeric(x) || is.integer(x)) ||
      (!is.null(length_expected) && length(x) != length_expected) ||
      length(x) < 1L || anyNA(x) || any(!is.finite(x)) ||
      any(x != floor(x)) || any(x < min_value) ||
      any(abs(x) > .Machine$integer.max)) {
    length_text <- if (is.null(length_expected)) {
      "one or more"
    } else {
      as.character(length_expected)
    }
    stop(
      name, " must contain ", length_text,
      " integer-valued finite number(s) >= ", min_value, "."
    )
  }
  as.integer(x)
}

mixedgp_validate_flag <- function(x, name, allow_null = FALSE) {
  if (isTRUE(allow_null) && is.null(x)) return(NULL)
  if (!is.logical(x) || length(x) != 1L || is.na(x)) {
    stop(name, " must be TRUE or FALSE.")
  }
  x
}

mixedgp_validate_column_names <- function(x, name) {
  column_names <- colnames(x)
  if (is.null(column_names)) return(NULL)
  if (length(column_names) != ncol(x) || anyNA(column_names) ||
      any(!nzchar(column_names)) || anyDuplicated(column_names)) {
    stop(
      name,
      " column names must be complete and unique when any names are supplied."
    )
  }
  column_names
}

mixedgp_parse_integer_text <- function(x, name, min_value = 0L) {
  if (!is.character(x) || length(x) != 1L || is.na(x) ||
      !grepl("^[0-9]+$", x)) {
    stop(name, " must be a base-10 integer >= ", min_value, ".")
  }
  value <- suppressWarnings(as.numeric(x))
  mixedgp_as_integer_strict(
    value, name, min_value = min_value, length_expected = 1L
  )
}

#' Allocate a core budget across datasets and posterior chains
#'
#' @param core_budget Positive integer CPU budget supplied by the user.
#' @param pending_datasets Number of datasets available to process.
#' @param chains Number of independent chains per fitted model.
#' @param sampling_iterations Retained transitions per chain, excluding warmup.
#' @param warmup Warmup transitions per chain.
#' @return A list with dataset and chain worker counts and explicit fit arguments.
#'   Continuation is never automatic; use [continue_eivgp()] explicitly.
#' @details On macOS/Linux, dataset workers may fork chain workers within the
#'   supplied budget. On Windows the current fork backend runs serially.
#'   Numerical-library threads should be limited to one before starting R.
#' @export
eivgp_run_settings <- function(core_budget, pending_datasets = 1L, chains = 4L,
                               sampling_iterations = 1250L, warmup = 500L) {
  check <- function(x, name, minimum = 1L) {
    mixedgp_as_integer_strict(x, name, minimum, 1L)
  }
  core_budget <- check(core_budget, "core_budget")
  pending_datasets <- check(pending_datasets, "pending_datasets", 0L)
  chains <- check(chains, "chains")
  sampling_iterations <- check(sampling_iterations, "sampling_iterations")
  warmup <- check(warmup, "warmup", 0L)
  total <- check(as.double(sampling_iterations) + warmup, "total iterations")
  chain_workers <- min(chains, core_budget)
  dataset_workers <- min(pending_datasets, max(1L, core_budget %/% chain_workers))
  if (.Platform$OS.type == "windows") {
    chain_workers <- 1L
    dataset_workers <- min(pending_datasets, 1L)
  }
  list(core_budget = core_budget, workers = dataset_workers,
       chain_workers = chain_workers, level = "hybrid",
       active_chain_limit = dataset_workers * chain_workers,
       retained_draws = as.double(chains) * sampling_iterations,
       automatic_continuation = FALSE,
       fit_args = list(n_iter = total, burn = warmup, thin = 1L,
                       n_chains = chains, parallel = chain_workers > 1L,
                       n_cores = chain_workers))
}

mixedgp_resolve_cores <- function(n_cores = NULL, reserve = 2L) {
  validate <- function(x, label) {
    mixedgp_as_integer_strict(
      x, label, min_value = 1L, length_expected = 1L
    )
  }

  if (!is.null(n_cores)) return(validate(n_cores, "n_cores"))

  env_cores <- Sys.getenv("MIXEDGP_CORES", unset = "")
  if (nzchar(env_cores)) {
    return(mixedgp_parse_integer_text(
      env_cores, "MIXEDGP_CORES", min_value = 1L
    ))
  }

  option_cores <- getOption("mixedgp.cores", NULL)
  if (!is.null(option_cores)) {
    return(validate(option_cores, "getOption('mixedgp.cores')"))
  }

  detected <- suppressWarnings(parallel::detectCores(logical = TRUE))
  if (length(detected) != 1L || is.na(detected) || detected < 1L) {
    detected <- 1L
  }
  reserve <- mixedgp_as_integer_strict(
    reserve, "reserve", min_value = 0L, length_expected = 1L
  )
  max(1L, detected - reserve)
}

mixedgp_parallel_backend <- function(n_cores = NULL, reserve = 2L) {
  cores <- mixedgp_resolve_cores(n_cores, reserve)
  can_fork <- .Platform$OS.type != "windows" && cores > 1L
  list(
    backend = if (can_fork) "fork" else "serial",
    cores = if (can_fork) cores else 1L,
    requested_cores = cores,
    os_type = .Platform$OS.type
  )
}

mixedgp_validate_task_seeds <- function(seeds, n_tasks) {
  if (is.null(seeds)) return(NULL)
  mixedgp_as_integer_strict(
    seeds, "seeds", min_value = 0L, length_expected = n_tasks
  )
}

mixedgp_rng_state <- function() {
  exists_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  list(
    kind = RNGkind(),
    exists = exists_seed,
    value = if (exists_seed) get(".Random.seed", envir = .GlobalEnv) else NULL
  )
}

mixedgp_restore_rng_state <- function(state) {
  if (!is.null(state$kind) && !identical(RNGkind(), state$kind)) {
    do.call(RNGkind, as.list(state$kind))
  }
  if (isTRUE(state$exists)) {
    assign(".Random.seed", state$value, envir = .GlobalEnv)
  } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
    rm(".Random.seed", envir = .GlobalEnv)
  }
  invisible(NULL)
}

## Internal continuation support. Retained draws are not terminal states:
## the last transition need not fall on the thinning schedule.
mixedgp_checkpoint_hash <- function(object) {
  cp <- object$checkpoint
  path <- tempfile("eivgp-checkpoint-", fileext = ".rds")
  on.exit(unlink(path), add = TRUE)
  saveRDS(list(data = object$data, control = object$control, kernel = object$kernel,
    samples = object$mcmc$samples_by_chain, states = cp$states,
    arguments = cp$arguments, iteration = cp$iteration, burn = cp$burn,
    thin = cp$thin, checkpoint_control = cp$control), path, compress = FALSE,
    version = 2L)
  unname(tools::md5sum(path))
}

mixedgp_resume_offset <- function(resume, n_iter, burn, thin, control) {
  if (is.null(resume)) return(0L)
  if (!identical(resume$version, 1L) || resume$iteration <= burn ||
      resume$iteration >= n_iter || !identical(resume$burn, burn) ||
      !identical(resume$thin, thin) || !identical(resume$control, control)) {
    stop("Incompatible continuation checkpoint or changed sampler settings.")
  }
  resume$iteration
}

mixedgp_append_chain_samples <- function(chains, resume) {
  if (is.null(resume)) return(chains)
  for (i in seq_along(chains)) {
    old <- resume$old_chains[[i]]
    for (name in names(old)[startsWith(names(old), "samples_")]) {
      a <- old[[name]]
      b <- chains[[i]][[name]]
      if (is.null(a)) next
      if (length(dim(a)) == 3L) {
        joined <- array(NA_real_, c(dim(a)[1L] + dim(b)[1L], dim(a)[-1L]))
        joined[seq_len(dim(a)[1L]), , ] <- a
        joined[dim(a)[1L] + seq_len(dim(b)[1L]), , ] <- b
      } else if (is.matrix(a)) {
        joined <- rbind(a, b)
      } else {
        joined <- c(a, b)
      }
      chains[[i]][[name]] <- joined
    }
    if (!is.null(old$initial_state)) chains[[i]]$initial_state <- old$initial_state
    counters <- setdiff(names(old$stats), c("chain", "seed"))
    for (name in counters) {
      chains[[i]]$stats[[name]] <- chains[[i]]$stats[[name]] + old$stats[[name]]
    }
  }
  chains
}

mixedgp_latent_anchor_status <- function(U_obs, calib_idx, d = NULL) {
  U_obs <- as.matrix(U_obs)
  if (is.null(d)) d <- ncol(U_obs)
  d <- as.integer(d)
  calib_idx <- sort(unique(as.integer(calib_idx)))
  if (length(d) != 1L || is.na(d) || d < 1L || ncol(U_obs) != d ||
      anyNA(calib_idx) || any(calib_idx < 1L | calib_idx > nrow(U_obs))) {
    stop("Invalid inputs to mixedgp_latent_anchor_status().")
  }
  required_rank <- d + 1L
  affine_rank <- if (length(calib_idx) == 0L) {
    0L
  } else {
    U_calib <- U_obs[calib_idx, , drop = FALSE]
    if (any(!is.finite(U_calib))) {
      stop("Calibrated latent inputs must be finite to assess anchoring.")
    }
    as.integer(qr(cbind(1, U_calib))$rank)
  }
  list(
    anchored = affine_rank >= required_rank,
    affine_rank = affine_rank,
    required_rank = required_rank,
    n_calibrated = length(calib_idx)
  )
}

mixedgp_parallel_chains_enabled <- function(
    parallel_level = c("chains", "replications", "none"),
    n_cores = NULL) {
  parallel_level <- match.arg(parallel_level)
  if (parallel_level != "chains") return(FALSE)
  identical(mixedgp_parallel_backend(n_cores)$backend, "fork")
}

mixedgp_parallel_lapply <- function(X,
                                    FUN,
                                    ...,
                                    n_cores = NULL,
                                    seeds = NULL,
                                    mc.preschedule = FALSE) {
  if (!is.function(FUN)) stop("FUN must be a function.")
  n_tasks <- length(X)
  if (n_tasks == 0L) return(list())
  seeds <- mixedgp_validate_task_seeds(seeds, n_tasks)
  if (!is.null(seeds)) {
    caller_rng <- mixedgp_rng_state()
    on.exit(mixedgp_restore_rng_state(caller_rng), add = TRUE)
  }
  backend <- mixedgp_parallel_backend(n_cores)
  cores <- min(backend$cores, n_tasks)
  extra_args <- list(...)

  worker <- function(ii) {
    if (!is.null(seeds)) set.seed(seeds[ii])
    do.call(FUN, c(list(X[[ii]]), extra_args))
  }

  if (backend$backend == "fork" && cores > 1L) {
    return(parallel::mclapply(
      seq_len(n_tasks),
      worker,
      mc.cores = cores,
      mc.preschedule = isTRUE(mc.preschedule),
      mc.set.seed = is.null(seeds),
      mc.cleanup = TRUE
    ))
  }
  lapply(seq_len(n_tasks), worker)
}

mixedgp_parallel_mapply <- function(FUN,
                                    ...,
                                    MoreArgs = NULL,
                                    SIMPLIFY = FALSE,
                                    USE.NAMES = TRUE,
                                    n_cores = NULL,
                                    seeds = NULL,
                                    mc.preschedule = FALSE) {
  if (!is.function(FUN)) stop("FUN must be a function.")
  varying_args <- list(...)
  if (length(varying_args) == 0L) {
    stop("Provide at least one varying argument through ....")
  }
  lengths <- lengths(varying_args)
  n_tasks <- max(lengths)
  if (n_tasks == 0L) return(if (isTRUE(SIMPLIFY)) vector() else list())
  if (any(!lengths %in% c(1L, n_tasks))) {
    stop("Each varying argument must have length one or the common task count.")
  }
  seeds <- mixedgp_validate_task_seeds(seeds, n_tasks)
  if (!is.null(seeds)) {
    caller_rng <- mixedgp_rng_state()
    on.exit(mixedgp_restore_rng_state(caller_rng), add = TRUE)
  }
  seed_arg <- if (is.null(seeds)) rep(NA_integer_, n_tasks) else seeds
  backend <- mixedgp_parallel_backend(n_cores)
  cores <- min(backend$cores, n_tasks)
  if (is.null(MoreArgs)) MoreArgs <- list()

  worker <- function(..., .mixedgp_seed) {
    if (!is.na(.mixedgp_seed)) set.seed(.mixedgp_seed)
    do.call(FUN, c(list(...), MoreArgs))
  }
  call_args <- c(
    list(FUN = worker),
    varying_args,
    list(.mixedgp_seed = seed_arg)
  )

  if (backend$backend == "fork" && cores > 1L) {
    return(do.call(
      parallel::mcmapply,
      c(
        call_args,
        list(
          SIMPLIFY = SIMPLIFY,
          USE.NAMES = USE.NAMES,
          mc.preschedule = isTRUE(mc.preschedule),
          mc.set.seed = is.null(seeds),
          mc.cores = cores
        )
      )
    ))
  }
  do.call(
    base::mapply,
    c(call_args, list(SIMPLIFY = SIMPLIFY, USE.NAMES = USE.NAMES))
  )
}

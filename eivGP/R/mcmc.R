
# ---- core/00_parallel_utils.R ----
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


# Explicit requests are ceilings, never permission to oversubscribe the host.
mixedgp_available_cores <- function() {
  detected <- suppressWarnings(tryCatch(parallel::detectCores(logical=TRUE),
                                       error=function(e) NA_integer_))
  if(length(detected)!=1L || is.na(detected) || !is.finite(detected) || detected<1L)
    return(1L)
  as.integer(detected)
}

mixedgp_resolve_cores <- function(n_cores = NULL, reserve = 2L) {
  validate <- function(x, label) {
    mixedgp_as_integer_strict(
      x, label, min_value = 1L, length_expected = 1L
    )
  }

  available <- mixedgp_available_cores()
  if (!is.null(n_cores)) return(min(validate(n_cores, "n_cores"), available))

  env_cores <- Sys.getenv("MIXEDGP_CORES", unset = "")
  if (nzchar(env_cores)) {
    return(min(mixedgp_parse_integer_text(
      env_cores, "MIXEDGP_CORES", min_value = 1L
    ), available))
  }

  option_cores <- getOption("mixedgp.cores", NULL)
  if (!is.null(option_cores)) {
    return(min(validate(option_cores, "getOption('mixedgp.cores')"), available))
  }

  reserve <- mixedgp_as_integer_strict(
    reserve, "reserve", min_value = 0L, length_expected = 1L
  )
  max(1L, available - reserve)
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



# One process-wide budget inherited by forked workers. Waiting coordinators own
# no slot. Atomic directory creation arbitrates slots across nested maps.
mixedgp_budget <- function() getOption("eivGP.core_pool", NULL)
mixedgp_budget_cores <- function(default = 1L) {
  pool <- mixedgp_budget()
  if (is.null(pool)) default else pool$cores
}
mixedgp_budget_owned <- function() {
  lease <- getOption("eivGP.core_lease", NULL)
  !is.null(lease) && identical(lease$pid, Sys.getpid())
}
mixedgp_budget_event <- function(event, slot) {
  path <- getOption("eivGP.core_events", NULL)
  if (is.null(path)) return(invisible(NULL))
  dir.create(path, recursive=TRUE, showWarnings=FALSE)
  cat(sprintf("%.6f,%d,%s,%d\n", as.numeric(Sys.time()), Sys.getpid(), event, slot),
      file=file.path(path,paste0(Sys.getpid(),".csv")), append=TRUE)
}
mixedgp_budget_acquire <- function(wait=TRUE) {
  pool <- mixedgp_budget()
  if (is.null(pool) || mixedgp_budget_owned()) return(invisible(NULL))
  repeat {
    if (!dir.exists(pool$path)) stop("Shared core pool closed while a worker was waiting.")
    for (i in seq_len(pool$cores)) {
      path <- file.path(pool$path, as.character(i))
      if (dir.create(path, showWarnings=FALSE)) {
        nonce <- basename(tempfile("lease-"))
        writeLines(c(as.character(Sys.getpid()),nonce),file.path(path,"owner"))
        options(eivGP.core_lease=list(pid=Sys.getpid(),path=path,slot=i,nonce=nonce))
        mixedgp_budget_event("acquire",i)
        return(invisible(TRUE))
      }
      owner <- file.path(path,"owner")
      # A killed leaf must not permanently consume a slot. No live owner is
      # evicted; the directory is private to this invocation of the workflow.
      reap <- paste0(path,".reap")
      if (file.exists(owner) && dir.create(reap,showWarnings=FALSE)) {
        pid <- suppressWarnings(as.integer(tryCatch(readLines(owner,warn=FALSE)[1L],error=function(e)character())))
        if (length(pid)==1L && !is.na(pid) && !isTRUE(tools::pskill(pid,0L)))
          unlink(path,recursive=TRUE)
        unlink(reap,recursive=TRUE)
      }
    }
    if(!wait)return(FALSE)
    Sys.sleep(.01)
  }
}
mixedgp_budget_release <- function() {
  if (mixedgp_budget_owned()) {
    lease <- getOption("eivGP.core_lease")
    mixedgp_budget_retire(lease)
  }
  options(eivGP.core_lease=NULL)
  invisible(NULL)
}
mixedgp_budget_retire <- function(lease) {
  if(!dir.exists(dirname(lease$path)))return(invisible(NULL))
  reap <- paste0(lease$path,".reap")
  while(!dir.create(reap,showWarnings=FALSE))Sys.sleep(.001)
  on.exit(unlink(reap,recursive=TRUE),add=TRUE)
  owner <- suppressWarnings(tryCatch(readLines(file.path(lease$path,"owner"),warn=FALSE),error=function(e)character()))
  if(length(owner)==2L && identical(owner[2L],lease$nonce)) {
    mixedgp_budget_event("release",lease$slot)
    unlink(lease$path,recursive=TRUE)
  }
  invisible(NULL)
}

# Reserve a slot before forking, so nested stages do not create an unbounded
# population of waiting subprocesses. Only idle coordinators wait for slots.
mixedgp_budget_map <- function(n,worker,cores,preschedule,set_seed) {
  groups <- if(preschedule)split(seq_len(n),(seq_len(n)-1L)%%cores)else as.list(seq_len(n))
  mixedgp_budget_yield(function() {
    running <- list();next_job <- 1L;result <- vector("list",n)
    on.exit({
      for(item in running) {
        try(tools::pskill(item$job$pid),silent=TRUE)
        try(parallel::mccollect(item$job,wait=TRUE),silent=TRUE)
        mixedgp_budget_retire(item$lease)
      }
    },add=TRUE)
    while(next_job<=length(groups)||length(running)) {
      while(next_job<=length(groups)&&length(running)<cores) {
        if(!mixedgp_budget_acquire(wait=FALSE))break
        lease <- getOption("eivGP.core_lease");indices <- groups[[next_job]]
        child <- tryCatch(parallel::mcparallel({
          lease$pid <- Sys.getpid()
          writeLines(c(as.character(lease$pid),lease$nonce),file.path(lease$path,"owner"))
          options(eivGP.core_lease=lease)
          tryCatch(lapply(indices,worker),finally=mixedgp_budget_release())
        },mc.set.seed=set_seed,silent=FALSE),error=function(e) {
          mixedgp_budget_release();stop(e)
        })
        # Ownership has been transferred to this child, not released.
        options(eivGP.core_lease=NULL)
        running[[as.character(child$pid)]] <- list(job=child,lease=lease,indices=indices)
        next_job <- next_job+1L
      }
      done <- if(length(running))parallel::mccollect(lapply(running,`[[`,"job"),wait=FALSE)else NULL
      if(is.null(done)) {Sys.sleep(.01);next}
      for(pid in names(done)) {
        item <- running[[pid]];value <- done[[pid]]
        mixedgp_budget_retire(item$lease);running[[pid]] <- NULL
        if(!is.list(value)||inherits(value,"try-error")||length(value)!=length(item$indices))
          stop("A shared-budget worker failed; inspect the worker error and retained checkpoints.")
        result[item$indices] <- value
      }
    }
    result
  })
}

mixedgp_budget_enter <- function(cores) {
  if (!is.null(mixedgp_budget())) return(function() invisible(NULL))
  requested <- mixedgp_as_integer_strict(cores,"cores",1L,1L)
  cores <- mixedgp_resolve_cores(requested)
  if (.Platform$OS.type != "unix") cores <- 1L
  if(cores < requested) message("Shared CPU budget: requested ",requested,
    "; using ",cores," available worker(s). Chain counts are unchanged.")
  path <- tempfile("eivgp-core-pool-"); dir.create(path)
  old <- options(eivGP.core_pool=list(path=path,cores=cores),eivGP.core_lease=NULL)
  mixedgp_budget_acquire()
  function() {
    mixedgp_budget_release()
    unlink(path,recursive=TRUE)
    options(old)
  }
}
mixedgp_budget_worker <- function(fun) {
  owned <- mixedgp_budget_owned()
  if (!owned) {
    mixedgp_budget_acquire()
    on.exit(mixedgp_budget_release(),add=TRUE)
  }
  fun()
}
mixedgp_budget_yield <- function(fun) {
  owned <- mixedgp_budget_owned()
  if (owned) {
    mixedgp_budget_release()
    on.exit(mixedgp_budget_acquire(),add=TRUE)
  }
  fun()
}

# Deterministic stages can borrow the full pool. Seeds are deliberately omitted
# here: callers must precompute randomness, or explicitly set task-local seeds.
mixedgp_stage_lapply <- function(X, FUN, ..., workers=mixedgp_budget_cores()) {
  mixedgp_parallel_lapply(X,FUN,...,n_cores=workers,mc.preschedule=TRUE,
                         .preserve_rng=TRUE)
}

mixedgp_parallel_lapply <- function(X,
                                    FUN,
                                    ...,
                                    n_cores = NULL,
                                    seeds = NULL,
                                    mc.preschedule = FALSE,
                                    .preserve_rng = FALSE) {
  if (!is.function(FUN)) stop("FUN must be a function.")
  n_tasks <- length(X)
  if (n_tasks == 0L) return(list())
  seeds <- mixedgp_validate_task_seeds(seeds, n_tasks)
  if (!is.null(seeds)) {
    caller_rng <- mixedgp_rng_state()
    on.exit(mixedgp_restore_rng_state(caller_rng), add = TRUE)
  }
  backend <- mixedgp_parallel_backend(if(is.null(n_cores)) mixedgp_budget_cores(NULL) else n_cores)
  cores <- min(backend$cores, n_tasks, mixedgp_budget_cores(backend$cores))
  extra_args <- list(...)

  worker <- function(ii) {
    if (!is.null(seeds)) set.seed(seeds[ii])
    value <- mixedgp_budget_worker(function() do.call(FUN, c(list(X[[ii]]), extra_args)))
    if(is.null(mixedgp_budget()))value else list(ok=TRUE,value=value)
  }
  unpack <- function(result) {
    if(is.null(mixedgp_budget()))return(result)
    if(any(!vapply(result,function(x)is.list(x) && isTRUE(x$ok),logical(1))))
      stop("A shared-budget worker failed; inspect the worker error and retained checkpoints.")
    lapply(result,`[[`,"value")
  }

  if (backend$backend == "fork" && cores > 1L) {
    if(!is.null(mixedgp_budget()))return(unpack(mixedgp_budget_map(n_tasks,worker,cores,
      isTRUE(mc.preschedule),is.null(seeds)&&!isTRUE(.preserve_rng))))
    result <- mixedgp_budget_yield(function() parallel::mclapply(
      seq_len(n_tasks),
      worker,
      mc.cores = cores,
      mc.preschedule = isTRUE(mc.preschedule),
      mc.set.seed = is.null(seeds) && !isTRUE(.preserve_rng),
      mc.cleanup = TRUE
    ))
    return(unpack(result))
  }
  unpack(lapply(seq_len(n_tasks), worker))
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
  if(!is.null(mixedgp_budget())) {
    result <- mixedgp_parallel_lapply(seq_len(n_tasks),function(i) {
      args <- lapply(varying_args,function(x)x[[if(length(x)==1L)1L else i]])
      do.call(FUN,c(args,MoreArgs))
    },n_cores=n_cores,seeds=seeds,mc.preschedule=mc.preschedule)
    if(isTRUE(USE.NAMES)) {
      first <- varying_args[[1L]]
      if(!is.null(names(first)))names(result)<-rep(names(first),length.out=n_tasks)
      else if(is.character(first))names(result)<-rep(first,length.out=n_tasks)
    }
    return(if(!identical(SIMPLIFY,FALSE))simplify2array(result,higher=(SIMPLIFY=="array"))else result)
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


as_study1_numeric_matrix <- function(x,
                                     name,
                                     nrow_expected = NULL,
                                     ncol_expected = NULL,
                                     expected_names = NULL) {
  input_names <- if (is.null(dim(x))) NULL else colnames(x)

  if (is.null(dim(x))) {
    if (is.null(ncol_expected) || ncol_expected == 1L) {
      x <- matrix(x, ncol = 1L)
    } else if (length(x) == ncol_expected) {
      x <- matrix(x, nrow = 1L)
    } else {
      stop(
        name, " is a vector, but a matrix with ", ncol_expected,
        " columns is required."
      )
    }
  }

  if (is.data.frame(x)) {
    numeric_columns <- vapply(
      x,
      function(z) is.numeric(z) || is.integer(z) || is.logical(z),
      logical(1)
    )
    if (any(!numeric_columns)) {
      stop(name, " must contain only numeric columns; encode factors first.")
    }
  } else if (!(is.numeric(x) || is.integer(x) || is.logical(x))) {
    stop(name, " must be numeric.")
  }

  out <- as.matrix(x)
  storage.mode(out) <- "double"

  if (ncol(out) < 1L) stop(name, " must contain at least one column.")
  mixedgp_validate_column_names(out, name)
  if (!is.null(nrow_expected) && nrow(out) != nrow_expected) {
    stop(name, " must have ", nrow_expected, " rows.")
  }
  if (!is.null(ncol_expected) && ncol(out) != ncol_expected) {
    stop(name, " must have ", ncol_expected, " columns.")
  }
  if (any(!is.finite(out))) stop(name, " contains missing or non-finite values.")

  if (!is.null(expected_names)) {
    if (length(expected_names) != ncol(out)) {
      stop(name, " has the wrong number of columns.")
    }
    if (!is.null(input_names) && all(nzchar(input_names))) {
      if (!setequal(input_names, expected_names)) {
        stop(name, " column names do not match the training predictors.")
      }
      out <- out[, match(expected_names, colnames(out)), drop = FALSE]
    }
    colnames(out) <- expected_names
  }

  out
}


maximin_lhs_1d <- function(n, lower = -2, upper = 2) {
  z <- (seq_len(n) - runif(n)) / n
  z <- sample(z)
  lower + (upper - lower) * z
}


make_class <- function(u, tau) {
  as.integer(cut(u, breaks = c(-Inf, tau, Inf), labels = FALSE))
}


rtruncnorm_vec <- function(mean, sd, lower, upper) {
  n <- max(length(mean), length(lower), length(upper))
  mean <- rep(mean, length.out = n)
  sd <- rep(sd, length.out = n)
  lower <- rep(lower, length.out = n)
  upper <- rep(upper, length.out = n)

  if (any(!is.finite(mean)) || any(!is.finite(sd)) || any(sd <= 0) ||
      anyNA(lower) || anyNA(upper) || any(lower >= upper)) {
    stop("Invalid parameters in rtruncnorm_vec().")
  }

  a <- (lower - mean) / sd
  b <- (upper - mean) / sd
  z <- numeric(n)

  ## Draw uniformly between two probabilities without subtracting nearly
  ## equal CDF values.  Positive-tail intervals use survival probabilities;
  ## negative-tail intervals use lower-tail probabilities.  qnorm(log.p=TRUE)
  ## then remains accurate even when the ordinary probabilities underflow.
  log_uniform_between <- function(log_lo, log_hi, size) {
    ratio <- exp(log_lo - log_hi)
    log_hi + log(ratio + runif(size) * (1 - ratio))
  }

  right <- a >= 0
  if (any(right)) {
    log_lo <- pnorm(b[right], lower.tail = FALSE, log.p = TRUE)
    log_hi <- pnorm(a[right], lower.tail = FALSE, log.p = TRUE)
    log_prob <- log_uniform_between(log_lo, log_hi, sum(right))
    z[right] <- qnorm(log_prob, lower.tail = FALSE, log.p = TRUE)
  }

  left <- !right & b <= 0
  if (any(left)) {
    log_lo <- pnorm(a[left], log.p = TRUE)
    log_hi <- pnorm(b[left], log.p = TRUE)
    log_prob <- log_uniform_between(log_lo, log_hi, sum(left))
    z[left] <- qnorm(log_prob, log.p = TRUE)
  }

  middle <- !right & !left
  if (any(middle)) {
    p_lo <- pnorm(a[middle])
    p_hi <- pnorm(b[middle])
    prob <- p_lo + runif(sum(middle)) * (p_hi - p_lo)
    z[middle] <- qnorm(prob)
  }

  out <- mean + sd * z
  tolerance <- 64 * .Machine$double.eps * pmax(1, abs(lower), abs(upper))
  violates <- !is.finite(out) |
    out < lower - tolerance | out > upper + tolerance
  if (any(violates)) {
    stop("Numerical failure in tail-stable truncated-normal sampling.")
  }
  pmin(pmax(out, lower), upper)
}


pairwise_sqdist <- function(a, b = NULL) {
  a <- as.matrix(a)
  if (is.null(b)) b <- a
  b <- as.matrix(b)
  
  aa <- rowSums(a^2)
  bb <- rowSums(b^2)
  
  out <- outer(aa, bb, "+") - 2 * tcrossprod(a, b)
  pmax(out, 0)
}


solve_chol <- function(U, b) {
  backsolve(U, forwardsolve(t(U), b))
}


rmvnorm_psd <- function(n, mean, Sigma, tolerance = NULL) {
  mean <- as.numeric(mean)
  Sigma <- as.matrix(Sigma)
  if (!identical(dim(Sigma), c(length(mean), length(mean)))) {
    stop("Sigma has incompatible dimensions in rmvnorm_psd().")
  }
  Sigma <- 0.5 * (Sigma + t(Sigma))
  eig <- eigen(Sigma, symmetric = TRUE)
  scale <- max(1, max(abs(eig$values)))
  if (is.null(tolerance)) {
    tolerance <- 100 * length(mean) * .Machine$double.eps * scale
  }
  if (min(eig$values) < -tolerance) {
    stop("Predictive covariance is not positive semidefinite within tolerance.")
  }
  values <- pmax(eig$values, 0)
  factor <- diag(sqrt(values), nrow = length(values)) %*% t(eig$vectors)
  Z <- matrix(rnorm(as.integer(n) * length(mean)), as.integer(n), length(mean))
  out <- sweep(Z %*% factor, 2, mean, "+")
  attr(out, "psd_tolerance") <- tolerance
  attr(out, "clipped_eigenvalues") <- sum(eig$values < 0)
  out
}


bounded_slice_update_1d <- function(x0, logf, w = 1,
                                 lower = -Inf, upper = Inf,
                                 max_steps_out = 50,
                                 max_iter = 200,
                                 fail_on_limit = FALSE) {
  if (length(w) != 1L || !is.finite(w) || w <= 0 ||
      max_steps_out < 1L || max_iter < 1L) {
    stop("Slice width and search budgets must be positive.")
  }
  if (upper <= lower) {
    if (isTRUE(fail_on_limit)) stop("Slice interval is empty.")
    return(list(x = x0, n_eval = 0))
  }
  
  if (length(x0) != 1L || !is.finite(x0) || x0 < lower || x0 > upper) {
    stop("Current slice state is outside its declared interval.")
  }
  
  f0 <- logf(x0)
  n_eval <- 1L
  
  if (!is.finite(f0)) {
    if (isTRUE(fail_on_limit)) stop("Current slice state has non-finite density.")
    return(list(x = x0, n_eval = n_eval))
  }
  
  logy <- f0 + log(runif(1))
  
  L <- x0 - runif(1) * w
  R <- L + w
  
  if (is.finite(lower)) L <- max(L, lower)
  if (is.finite(upper)) R <- min(R, upper)
  
  J <- floor(runif(1) * max_steps_out)
  K <- max_steps_out - 1L - J
  
  while (J > 0 && (!is.finite(lower) || L > lower)) {
    fL <- logf(L)
    n_eval <- n_eval + 1L
    if (!is.finite(fL) || fL <= logy) break
    L <- L - w
    if (is.finite(lower)) L <- max(L, lower)
    J <- J - 1L
  }
  
  while (K > 0 && (!is.finite(upper) || R < upper)) {
    fR <- logf(R)
    n_eval <- n_eval + 1L
    if (!is.finite(fR) || fR <= logy) break
    R <- R + w
    if (is.finite(upper)) R <- min(R, upper)
    K <- K - 1L
  }
  
  for (iter in seq_len(max_iter)) {
    x1 <- runif(1, L, R)
    f1 <- logf(x1)
    n_eval <- n_eval + 1L
    
    if (is.finite(f1) && f1 >= logy) {
      return(list(x = x1, n_eval = n_eval))
    }
    
    if (x1 < x0) {
      L <- x1
    } else {
      R <- x1
    }
  }
  
  if (isTRUE(fail_on_limit)) {
    stop("Slice shrinkage exceeded max_iter.")
  }
  list(x = x0, n_eval = n_eval)
}


full_interval_slice_update <- function(x0, logf, lower, upper,
                                       max_iter = 200,
                                       fail_on_limit = FALSE) {
  if (!is.finite(lower) || !is.finite(upper)) {
    stop("full_interval_slice_update requires finite lower and upper bounds.")
  }
  
  if (upper <= lower) {
    if (isTRUE(fail_on_limit)) stop("Full slice interval is empty.")
    return(list(x = x0, n_eval = 0))
  }
  
  if (length(x0) != 1L || !is.finite(x0) || x0 < lower || x0 > upper) {
    stop("Current full-interval slice state is outside its declared interval.")
  }
  
  f0 <- logf(x0)
  n_eval <- 1L
  
  if (!is.finite(f0)) {
    if (isTRUE(fail_on_limit)) {
      stop("Current full-interval slice state has non-finite density.")
    }
    return(list(x = x0, n_eval = n_eval))
  }
  
  logy <- f0 + log(runif(1))
  L <- lower
  R <- upper
  
  for (iter in seq_len(max_iter)) {
    x1 <- runif(1, L, R)
    f1 <- logf(x1)
    n_eval <- n_eval + 1L
    
    if (is.finite(f1) && f1 >= logy) {
      return(list(x = x1, n_eval = n_eval))
    }
    
    if (x1 < x0) {
      L <- x1
    } else {
      R <- x1
    }
  }
  
  if (isTRUE(fail_on_limit)) {
    stop("Full-interval slice shrinkage exceeded max_iter.")
  }
  list(x = x0, n_eval = n_eval)
}


.study1_category_moment_cache <- new.env(parent = emptyenv())


study1_category_trig_moments <- function(tau,
                                         frequency = 0.9 * base::pi) {
  tau <- as.numeric(tau)
  if (length(tau) < 1L || any(!is.finite(tau)) ||
      is.unsorted(tau, strictly = TRUE)) {
    stop("tau must be a finite, strictly increasing vector.")
  }

  cache_key <- paste(
    formatC(c(frequency, tau), digits = 16L, format = "fg"),
    collapse = "|"
  )
  if (exists(cache_key, envir = .study1_category_moment_cache, inherits = FALSE)) {
    return(get(cache_key, envir = .study1_category_moment_cache, inherits = FALSE))
  }

  bounds <- c(-Inf, tau, Inf)
  m <- length(bounds) - 1L
  probability <- diff(pnorm(bounds))
  if (any(probability <= 0)) stop("tau creates an empty ordinal category.")

  integrate_moment <- function(fun, lower, upper) {
    stats::integrate(
      function(z) fun(frequency * z) * dnorm(z),
      lower = lower,
      upper = upper,
      rel.tol = 1e-11,
      subdivisions = 200L,
      stop.on.error = TRUE
    )$value
  }
  cosine <- sine <- numeric(m)
  for (r in seq_len(m)) {
    cosine[r] <- integrate_moment(cos, bounds[r], bounds[r + 1L]) /
      probability[r]
    sine[r] <- integrate_moment(sin, bounds[r], bounds[r + 1L]) /
      probability[r]
  }

  out <- list(cosine = cosine, sine = sine, probability = probability)
  assign(cache_key, out, envir = .study1_category_moment_cache)
  out
}


f0_1d <- function(x,
                  u,
                  scenario = c(
                    "active", "inactive", "category_sufficient",
                    "heterogeneity_continuum"
                  ),
                  c_ord = NULL,
                  tau = c(-5 / 3, -10 / 9, 0, 10 / 9, 5 / 3),
                  heterogeneity_eta = 1) {
  scenario <- match.arg(scenario)
  
  if (scenario == "inactive") {
    return(cos(0.9 * base::pi * u))
  }

  x_main <- 0.5 * sin(base::pi * x / 2)
  if (scenario %in% c("category_sufficient", "heterogeneity_continuum")) {
    if (is.null(c_ord)) c_ord <- make_class(u, tau)
    c_ord <- as.integer(c_ord)
    moments <- study1_category_trig_moments(tau)
    if (length(c_ord) != max(length(x), length(u)) ||
        anyNA(c_ord) || any(c_ord < 1L | c_ord > length(moments$cosine))) {
      stop("c_ord must contain one valid ordinal level per response input.")
    }
    conditional_mean <- moments$cosine[c_ord] +
      0.35 * (x / 2) * moments$sine[c_ord]
    if (scenario == "category_sufficient") return(x_main + conditional_mean)

    heterogeneity_eta <- as.numeric(heterogeneity_eta)
    if (length(heterogeneity_eta) != 1L || !is.finite(heterogeneity_eta) ||
        heterogeneity_eta < 0) {
      stop("heterogeneity_eta must be one finite nonnegative number.")
    }
    continuous_effect <- cos(0.9 * base::pi * u) +
      0.35 * (x / 2) * sin(0.9 * base::pi * u)
    return(
      x_main + conditional_mean +
        heterogeneity_eta * (continuous_effect - conditional_mean)
    )
  }
  
  x_main + cos(0.9 * base::pi * u) +
    0.35 * (x / 2) * sin(0.9 * base::pi * u)
}


m0_1d <- function(x,
                  c_ord,
                  tau = c(-5 / 3, -10 / 9, 0, 10 / 9, 5 / 3),
                  scenario = c(
                    "active", "inactive", "category_sufficient",
                    "heterogeneity_continuum"
                  )) {
  scenario <- match.arg(scenario)
  x <- as.numeric(x)
  c_ord <- as.integer(c_ord)
  n_out <- max(length(x), length(c_ord))
  x <- rep(x, length.out = n_out)
  c_ord <- rep(c_ord, length.out = n_out)
  moments <- study1_category_trig_moments(tau)
  if (anyNA(c_ord) || any(c_ord < 1L | c_ord > length(moments$cosine))) {
    stop("c_ord contains invalid ordinal levels in m0_1d().")
  }
  conditional_continuous <- moments$cosine[c_ord] +
    0.35 * (x / 2) * moments$sine[c_ord]
  if (scenario == "inactive") return(moments$cosine[c_ord])
  0.5 * sin(base::pi * x / 2) + conditional_continuous
}


simulate_1d_data <- function(n = 100,
                             n_test = 500,
                             m = 6,
                             scenario = c(
                               "active", "inactive", "category_sufficient",
                               "heterogeneity_continuum"
                             ),
                             sigma_eps = 0.1,
                             seed = NULL,
                             threshold_design = c("imbalanced", "balanced"),
                             min_class_count = NULL,
                             heterogeneity_eta = 1) {
  scenario <- match.arg(scenario)
  threshold_design <- match.arg(threshold_design)
  
  if (!is.null(seed)) {
    set.seed(seed)
  }

  if (threshold_design == "imbalanced") {
    if (m != 6L) stop("The imbalanced Study I design is defined for m = 6.")
    tau_true <- c(-5 / 3, -10 / 9, 0, 10 / 9, 5 / 3)
    if (is.null(min_class_count)) min_class_count <- 3L
  } else {
    tau_true <- qnorm(seq_len(m - 1L) / m)
    if (is.null(min_class_count)) min_class_count <- 0L
  }
  min_class_count <- as.integer(min_class_count)
  if (length(min_class_count) != 1L || is.na(min_class_count) ||
      min_class_count < 0L) {
    stop("min_class_count must be one nonnegative integer.")
  }
  
  repeat {
    x <- maximin_lhs_1d(n, lower = -2, upper = 2)
    u <- rnorm(n)
    c_ord <- make_class(u, tau_true)
    
    if (min_class_count == 0L ||
        all(tabulate(c_ord, nbins = m) >= min_class_count)) break
  }
  
  f <- f0_1d(
    x, u, scenario = scenario, c_ord = c_ord, tau = tau_true,
    heterogeneity_eta = heterogeneity_eta
  )
  y <- f + rnorm(n, mean = 0, sd = sigma_eps)
  
  x_test <- runif(n_test, -2, 2)
  u_test <- rnorm(n_test)
  c_test <- make_class(u_test, tau_true)
  f_test <- f0_1d(
    x_test, u_test, scenario = scenario, c_ord = c_test, tau = tau_true,
    heterogeneity_eta = heterogeneity_eta
  )
  y_test <- f_test + rnorm(n_test, mean = 0, sd = sigma_eps)
  
  list(
    train = data.frame(
      x = x,
      u = u,
      c = c_ord,
      y = y,
      f = f
    ),
    test = data.frame(
      x = x_test,
      u = u_test,
      c = c_test,
      y = y_test,
      f = f_test
    ),
    tau_true = tau_true,
    sigma_eps = sigma_eps,
    scenario = scenario,
    m = m,
    threshold_design = threshold_design,
    min_class_count = min_class_count,
    heterogeneity_eta = heterogeneity_eta
  )
}


make_nested_calibration_sets <- function(n, calib_grid, seed = NULL) {
  n <- as.integer(n)
  calib_grid <- as.integer(calib_grid)
  if (length(n) != 1L || is.na(n) || n < 1L ||
      length(calib_grid) < 1L || anyNA(calib_grid) ||
      any(calib_grid < 0L | calib_grid > n)) {
    stop("Calibration sizes must be integers between zero and n.")
  }
  if (!is.null(seed)) {
    set.seed(seed)
  }
  
  ord <- sample(seq_len(n))
  
  out <- lapply(calib_grid, function(k) {
    if (k == 0L) integer(0) else sort(ord[seq_len(k)])
  })
  
  names(out) <- as.character(calib_grid)
  out
}


make_theta_spec_multix <- function(p_x) {
  if (!is.finite(p_x) || p_x < 1L) stop("p_x must be positive")
  parameter <- c("log_rho", paste0("log_theta_x", seq_len(p_x)), "log_theta_u")
  bounds <- rbind(
    c(log(0.05), log(100)),
    matrix(rep(c(log(1e-4), log(100)), p_x + 1L), ncol = 2, byrow = TRUE)
  )
  rownames(bounds) <- parameter
  list(
    p_x = p_x,
    parameter = parameter,
    prior_mean = c(log(3), rep(log(0.5), p_x + 1L)),
    prior_sd = rep(1.5, p_x + 2L),
    bounds = bounds,
    x_index = 1L + seq_len(p_x),
    u_index = p_x + 2L
  )
}


new_gp_cache_1d <- function(use_blocks = TRUE) {
  cache <- new.env(parent = emptyenv())
  cache$entry <- NULL
  cache$use_blocks <- isTRUE(use_blocks)
  cache$full_factorizations <- 0L
  cache$cache_hits <- 0L
  cache$block_setups <- 0L
  cache$block_evaluations <- 0L
  cache$block_fallbacks <- 0L
  class(cache) <- "eivgp_cache_1d"
  cache
}


gp_cache_matches_1d <- function(cache, y, u, Dx_list, logtheta,
                                theta_spec, kernel, matern_nu) {
  if (is.null(cache)) return(FALSE)
  if (!is.environment(cache) || !inherits(cache, "eivgp_cache_1d")) {
    stop("gp_cache must be created by new_gp_cache_1d().")
  }
  entry <- cache$entry
  !is.null(entry) && identical(entry$y, y) && identical(entry$u, u) &&
    identical(entry$Dx_list, Dx_list) && identical(entry$logtheta, logtheta) &&
    identical(entry$theta_spec, theta_spec) && identical(entry$kernel, kernel) &&
    identical(entry$matern_nu, matern_nu)
}


gp_cache_store_1d <- function(cache, state, y, u, Dx_list, logtheta,
                              theta_spec, kernel, matern_nu) {
  if (!is.null(cache)) {
    cache$entry <- list(state = state, y = y, u = u, Dx_list = Dx_list,
                        logtheta = logtheta, theta_spec = theta_spec,
                        kernel = kernel, matern_nu = matern_nu)
  }
  invisible(state)
}


gp_update_loglik_1d <- function(state, n, collapse_noise = FALSE,
                                a_eps0 = 2, b_eps0 = 0.05) {
  if (!collapse_noise) return(state$loglik)
  ## A = I + rho^2 R and sigma2 ~ IG(a_eps0, b_eps0). Integrating
  ## sigma2 removes neither the latent-input prior nor the theta prior.
  ## The determinant and quadratic form are shared by dense/Schur caches.
  -0.5 * (n * log(2 * base::pi) + state$logdetA) +
    a_eps0 * log(b_eps0) - lgamma(a_eps0) + lgamma(a_eps0 + n / 2) -
    (a_eps0 + n / 2) * log(b_eps0 + 0.5 * state$quad)
}


gp_state_1d <- function(y, u, Dx_list, logtheta, sigma2_eps,
                        theta_spec = NULL,
                        kernel = "se",
                        matern_nu = 2.5,
                        gp_cache = NULL,
                        require_factor = TRUE) {
  if (!is.list(Dx_list)) Dx_list <- list(as.matrix(Dx_list))
  if (is.null(theta_spec)) {
    theta_spec <- make_theta_spec_multix(length(Dx_list))
  }
  n <- length(y)
  if (gp_cache_matches_1d(gp_cache, y, u, Dx_list, logtheta,
                          theta_spec, kernel, matern_nu) &&
      (!isTRUE(require_factor) || !is.null(gp_cache$entry$state$cholA))) {
    state <- gp_cache$entry$state
    gp_cache$cache_hits <- gp_cache$cache_hits + 1L
    state$loglik <- -0.5 * (n * log(2 * base::pi * sigma2_eps) +
                             state$logdetA + state$quad / sigma2_eps)
    return(state)
  }

  rho <- exp(logtheta[1])
  theta_x <- exp(logtheta[theta_spec$x_index])
  theta_u <- exp(logtheta[theta_spec$u_index])

  Du <- pairwise_sqdist(matrix(u, ncol = 1))
  D2 <- weighted_distance_multix(Dx_list, theta_x) + theta_u * Du
  R <- kernel_from_weighted_sqdist_1d(D2, kernel, matern_nu)
  
  A <- rho^2 * R + diag(n)
  
  U <- safe_chol_1d(A)
  if (!is.null(gp_cache)) {
    gp_cache$full_factorizations <- gp_cache$full_factorizations + 1L
  }
  Ainv_y <- solve_chol(U, y)
  
  logdetA <- 2 * sum(log(diag(U)))
  quad <- sum(y * Ainv_y)
  
  loglik <- -0.5 * (
    n * log(2 * base::pi * sigma2_eps) +
      logdetA +
      quad / sigma2_eps
  )
  
  state <- list(
    loglik = loglik,
    R = R,
    A = A,
    cholA = U,
    Ainv_y = Ainv_y,
    logdetA = logdetA,
    quad = quad
  )
  gp_cache_store_1d(gp_cache, state, y, u, Dx_list, logtheta,
                    theta_spec, kernel, matern_nu)
  state
}


gp_block_state_1d <- function(prepared, u_block, sigma2_eps, gp_cache) {
  a <- prepared
  if (length(u_block) != a$b || any(!is.finite(u_block))) {
    stop("A GP block proposal has incompatible or nonfinite latent values.")
  }
  u_prop <- a$u
  u_prop[a$block_idx] <- u_block
  if (gp_cache_matches_1d(gp_cache, a$y, u_prop, a$Dx_list, a$logtheta,
                          a$theta_spec, a$kernel, a$matern_nu)) {
    return(gp_state_1d(
      a$y, u_prop, a$Dx_list, a$logtheta, sigma2_eps, a$theta_spec,
      a$kernel, a$matern_nu, gp_cache = gp_cache, require_factor = FALSE
    ))
  }
  gp_cache$block_evaluations <- gp_cache$block_evaluations + 1L
  ## Only the proposed rows/columns change. Factor the fixed complement once,
  ## then evaluate the exact Gaussian conditional through its Schur complement.
  candidate <- tryCatch({
    cross <- a$rho2 * kernel_from_weighted_sqdist_1d(
      a$Dx_cross + a$theta_u * pairwise_sqdist(matrix(u_block, ncol = 1),
                                                matrix(a$u[a$rest], ncol = 1)),
      a$kernel, a$matern_nu
    )
    block <- diag(a$b) + a$rho2 * kernel_from_weighted_sqdist_1d(
      a$Dx_block + a$theta_u * pairwise_sqdist(matrix(u_block, ncol = 1)),
      a$kernel, a$matern_nu
    )
    V <- forwardsolve(t(a$chol_rest), t(cross))
    schur <- block - crossprod(V)
    schur <- 0.5 * (schur + t(schur))
    H <- chol(schur)
    residual <- a$y[a$block_idx] - as.numeric(cross %*% a$alpha_rest)
    quad <- a$quad_rest + sum(residual * solve_chol(H, residual))
    logdet <- a$logdet_rest + 2 * sum(log(diag(H)))
    if (!is.finite(quad) || quad < 0 || !is.finite(logdet)) {
      stop("Nonfinite Schur-complement likelihood calculation.")
    }
    A <- a$A
    A[a$block_idx, a$rest] <- cross
    A[a$rest, a$block_idx] <- t(cross)
    A[a$block_idx, a$block_idx] <- block
    list(A = A, cholA = NULL, R = NULL, Ainv_y = NULL,
         quad = quad, logdetA = logdet,
         loglik = -0.5 * (a$n * log(2 * base::pi * sigma2_eps) +
                           logdet + quad / sigma2_eps))
  }, error = function(e) NULL)
  if (is.null(candidate)) {
    ## Retry the same covariance with the dense factorization. This does not
    ## regularize the matrix or reject a proposal because the shortcut failed.
    gp_cache$block_fallbacks <- gp_cache$block_fallbacks + 1L
    return(gp_state_1d(
      a$y, u_prop, a$Dx_list, a$logtheta, sigma2_eps, a$theta_spec,
      a$kernel, a$matern_nu, gp_cache = gp_cache, require_factor = TRUE
    ))
  }
  gp_cache_store_1d(gp_cache, candidate, a$y, u_prop, a$Dx_list, a$logtheta,
                    a$theta_spec, a$kernel, a$matern_nu)
  candidate
}


make_default_control <- function(n, n_mis,
                                 preset = c("fast", "balanced", "thorough"),
                                 p_x = 1L) {
  preset <- match.arg(preset)
  
  if (preset == "fast") {
    local_frac <- 0.03
    theta_update_every <- 10L
    block_ess_every <- 2L
    n_blocks_per_iter <- 1L
    global_ess_every <- 50L
    full_local_every <- 100L
  }
  
  if (preset == "balanced") {
    local_frac <- 0.06
    theta_update_every <- 5L
    block_ess_every <- 1L
    n_blocks_per_iter <- 1L
    global_ess_every <- 25L
    full_local_every <- 50L
  }
  
  if (preset == "thorough") {
    local_frac <- 0.10
    theta_update_every <- 3L
    block_ess_every <- 1L
    n_blocks_per_iter <- 2L
    global_ess_every <- 10L
    full_local_every <- 25L
  }
  
  local_per_iter <- min(
    n_mis,
    max(3L, ceiling(local_frac * max(n_mis, 1L)))
  )
  
  ess_block_size <- min(
    n_mis,
    max(5L, ceiling(sqrt(max(n_mis, 1L))))
  )
  
  list(
    preset = preset,
    local_per_iter = local_per_iter,
    full_local_every = full_local_every,
    ess_block_size = ess_block_size,
    block_ess_every = block_ess_every,
    n_blocks_per_iter = n_blocks_per_iter,
    global_ess_every = global_ess_every,
    theta_update_every = theta_update_every,
    theta_slice_width_init = rep(1.0, p_x + 2L),
    adapt_theta_width = TRUE,
    adapt_every = 100L,
    adapt_window = 500L,
    theta_width_min = 0.20,
    theta_width_max = 2.50
  )
}


check_constraints_1d <- function(u, c_ord, tau) {
  lower <- c(-Inf, tau)[c_ord]
  upper <- c(tau, Inf)[c_ord]
  
  all(u > lower & u <= upper)
}


log_prior_logtheta <- function(logtheta, theta_spec = NULL,
                               mean = NULL, sd = NULL, bounds = NULL) {
  if (is.null(theta_spec)) {
    p_x <- length(logtheta) - 2L
    if (!is.finite(p_x) || p_x < 1L) {
      stop("logtheta must have length p_x + 2")
    }
    theta_spec <- make_theta_spec_multix(p_x)
  }
  if (!is.null(mean)) theta_spec$prior_mean <- mean
  if (!is.null(sd)) theta_spec$prior_sd <- sd
  if (!is.null(bounds)) theta_spec$bounds <- bounds
  if (any(logtheta < theta_spec$bounds[, 1]) ||
      any(logtheta > theta_spec$bounds[, 2])) {
    return(-Inf)
  }
  
  sum(dnorm(
    logtheta,
    mean = theta_spec$prior_mean,
    sd = theta_spec$prior_sd,
    log = TRUE
  ))
}


gp_loglik_with_constraints_1d <- function(y, u, Dx_list, c_ord, tau,
                                          logtheta, sigma2_eps,
                                          theta_spec = NULL,
                                          kernel = "se",
                                          matern_nu = 2.5,
                                          gp_cache = NULL,
                                          collapse_noise = FALSE) {
  if (!is.list(Dx_list)) Dx_list <- list(as.matrix(Dx_list))
  if (is.null(theta_spec)) {
    theta_spec <- make_theta_spec_multix(length(Dx_list))
  }
  if (!check_constraints_1d(u, c_ord, tau)) {
    return(-Inf)
  }
  
  state <- gp_state_1d(
    y, u, Dx_list, logtheta, if (collapse_noise) 1 else sigma2_eps, theta_spec,
    kernel = kernel, matern_nu = matern_nu,
    gp_cache = gp_cache, require_factor = FALSE
  )
  gp_update_loglik_1d(state, length(y), collapse_noise)
}


theta_logpost_1d <- function(y, u, Dx_list, logtheta, sigma2_eps,
                             theta_spec = NULL,
                             prior_mean = NULL, prior_sd = NULL, bounds = NULL,
                             kernel = "se", matern_nu = 2.5,
                             gp_cache = NULL,
                             collapse_noise = FALSE) {
  if (!is.list(Dx_list)) Dx_list <- list(as.matrix(Dx_list))
  if (is.null(theta_spec)) {
    theta_spec <- make_theta_spec_multix(length(Dx_list))
  }
  lp <- log_prior_logtheta(
    logtheta, theta_spec,
    mean = prior_mean, sd = prior_sd, bounds = bounds
  )
  
  if (!is.finite(lp)) {
    return(-Inf)
  }
  
  state <- gp_state_1d(
    y, u, Dx_list, logtheta, if (collapse_noise) 1 else sigma2_eps, theta_spec,
    kernel = kernel, matern_nu = matern_nu,
    gp_cache = gp_cache, require_factor = FALSE
  )
  gp_update_loglik_1d(state, length(y), collapse_noise) + lp
}


initialize_tau_1d <- function(c_ord, u_obs, calib_idx, m, tau_bound = 8) {
  n <- length(c_ord)
  
  counts <- tabulate(c_ord, nbins = m)
  probs <- cumsum(counts)[1:(m - 1)] / n
  probs <- pmin(pmax(probs, 0.03), 0.97)
  
  tau <- qnorm(probs)
  obs_c <- c_ord[calib_idx]
  obs_u <- u_obs[calib_idx]
  for (j in seq_len(m - 1)) {
    lower <- max(c(-tau_bound, if (j > 1L) tau[j - 1L],
                   obs_u[obs_c <= j]))
    upper <- min(c(tau_bound, obs_u[obs_c > j]))
    if (lower >= upper) {
      stop("Calibrated u values are incompatible with ordinal labels or threshold bounds.")
    }
    ## An absolute inward margin can exceed a valid calibration gap. A
    ## midpoint also leaves room for later thresholds sharing this gap.
    if (!is.finite(tau[j]) || tau[j] <= lower || tau[j] >= upper) {
      tau[j] <- lower + (upper - lower) / 2
    }
    if (tau[j] <= lower || tau[j] >= upper) {
      stop("No representable interior value exists in a threshold initialization interval.")
    }
  }
  
  tau
}


update_tau_1d <- function(tau, u, c_ord, m, tau_bound = 8) {
  tau_new <- tau
  
  for (j in seq_len(m - 1)) {
    L <- max(
      c(
        -tau_bound,
        if (j > 1) tau_new[j - 1] else -Inf,
        u[c_ord <= j]
      ),
      na.rm = TRUE
    )
    
    U <- min(
      c(
        tau_bound,
        if (j < m - 1) tau_new[j + 1] else Inf,
        u[c_ord > j]
      ),
      na.rm = TRUE
    )
    
    if (!is.finite(L) || !is.finite(U) || L >= U) {
      stop("Threshold Gibbs update encountered an empty conditional interval.")
    }
    tau_new[j] <- runif(1, L, U)
  }
  
  tau_new
}


assert_threshold_state_1d <- function(u,
                                      tau,
                                      c_ord,
                                      m,
                                      calib_idx = integer(0),
                                      u_obs = NULL,
                                      logtheta = NULL,
                                      sigma2_eps = NULL,
                                      tolerance = 1e-10) {
  if (length(u) != length(c_ord) || any(!is.finite(u)) ||
      length(tau) != m - 1L || any(!is.finite(tau)) ||
      is.unsorted(tau, strictly = TRUE)) {
    stop("The deterministic-threshold sampler state is invalid.")
  }
  lower <- c(-Inf, tau)[c_ord]
  upper <- c(tau, Inf)[c_ord]
  if (any(u <= lower | u > upper)) {
    stop("A latent input left the interval implied by its ordinal category.")
  }
  if (length(calib_idx) > 0L) {
    if (is.null(u_obs) || anyNA(u_obs[calib_idx]) ||
        any(u[calib_idx] != u_obs[calib_idx])) {
      stop("A calibrated latent input moved away from its observed value.")
    }
  }
  if (!is.null(logtheta) && any(!is.finite(logtheta))) {
    stop("A GP hyperparameter is nonfinite.")
  }
  if (!is.null(sigma2_eps) &&
      (length(sigma2_eps) != 1L || !is.finite(sigma2_eps) || sigma2_eps <= 0)) {
    stop("The observation-noise variance is not positive and finite.")
  }
  invisible(TRUE)
}


sample_sigma2_eps_1d <- function(y, u, Dx_list, logtheta, theta_spec = NULL,
                                 a_eps0 = 2, b_eps0 = 0.05,
                                 kernel = "se", matern_nu = 2.5,
                                 gp_cache = NULL) {
  if (!is.list(Dx_list)) Dx_list <- list(as.matrix(Dx_list))
  if (is.null(theta_spec)) {
    theta_spec <- make_theta_spec_multix(length(Dx_list))
  }
  n <- length(y)
  
  ## The factor, determinant and quadratic form do not depend on sigma2.
  quad <- gp_state_1d(
    y, u, Dx_list, logtheta, sigma2_eps = 1, theta_spec = theta_spec,
    kernel = kernel, matern_nu = matern_nu,
    gp_cache = gp_cache, require_factor = FALSE
  )$quad
  
  shape <- a_eps0 + n / 2
  rate <- b_eps0 + 0.5 * quad
  
  1 / rgamma(1, shape = shape, rate = rate)
}


update_u_ess_block_1d <- function(y, u, Dx_list, c_ord, tau, logtheta,
                                      sigma2_eps, theta_spec,
                                      block_idx, max_try = 300L,
                                      kernel = "se", matern_nu = 2.5,
                                      gp_cache = NULL,
                                      collapse_noise = FALSE) {
  if (length(block_idx) == 0) {
    return(list(u = u, n_eval = 0L, accepted = TRUE))
  }
  
  u_block <- u[block_idx]
  block_state <- prepare_gp_block_1d(
    y, u, Dx_list, logtheta, theta_spec, block_idx,
    kernel, matern_nu, gp_cache
  )
  
  loglik_fun <- function(u_block_prop) {
    u_prop <- u
    u_prop[block_idx] <- u_block_prop
    if (!is.null(block_state)) {
      if (!check_constraints_1d(u_prop, c_ord, tau)) return(-Inf)
      candidate <- gp_block_state_1d(
        block_state, u_block_prop, if (collapse_noise) 1 else sigma2_eps, gp_cache
      )
      return(gp_update_loglik_1d(candidate, length(y), collapse_noise))
    }
    
    gp_loglik_with_constraints_1d(
      y = y,
      u = u_prop,
      Dx_list = Dx_list,
      c_ord = c_ord,
      tau = tau,
      logtheta = logtheta,
      sigma2_eps = sigma2_eps,
      theta_spec = theta_spec,
      kernel = kernel,
      matern_nu = matern_nu,
      gp_cache = gp_cache,
      collapse_noise = collapse_noise
    )
  }
  
  loglik_cur <- loglik_fun(u_block)
  
  if (!is.finite(loglik_cur)) {
    stop("Current latent state violates constraints or has invalid likelihood.")
  }
  
  nu <- rnorm(length(u_block))
  
  logy <- loglik_cur + log(runif(1))
  
  angle <- runif(1, 0, 2 * base::pi)
  angle_min <- angle - 2 * base::pi
  angle_max <- angle
  
  n_eval <- 1L
  
  for (try_id in seq_len(max_try)) {
    u_block_prop <- u_block * cos(angle) + nu * sin(angle)
    
    loglik_prop <- loglik_fun(u_block_prop)
    n_eval <- n_eval + 1L
    
    if (is.finite(loglik_prop) && loglik_prop >= logy) {
      u_new <- u
      u_new[block_idx] <- u_block_prop
      return(list(u = u_new, n_eval = n_eval, accepted = TRUE))
    }
    
    if (angle < 0) {
      angle_min <- angle
    } else {
      angle_max <- angle
    }
    
    angle <- runif(1, angle_min, angle_max)
  }
  
  stop("Elliptical-slice update exceeded max_try.")
}


update_u_local_z_slice_1d <- function(y, u, Dx_list, c_ord, tau,
                                          logtheta, sigma2_eps, theta_spec,
                                          update_idx,
                                          kernel = "se", matern_nu = 2.5,
                                          gp_cache = NULL,
                                          collapse_noise = FALSE) {
  n_eval_total <- 0L
  update_idx <- update_idx[sample.int(length(update_idx))]
  
  for (i in update_idx) {
    lower_u <- c(-Inf, tau)[c_ord[i]]
    upper_u <- c(tau, Inf)[c_ord[i]]
    block_state <- prepare_gp_block_1d(
      y, u, Dx_list, logtheta, theta_spec, i,
      kernel, matern_nu, gp_cache
    )
    local_gp_loglik <- function(value) {
      if (!is.null(block_state)) {
        candidate <- gp_block_state_1d(
          block_state, value, if (collapse_noise) 1 else sigma2_eps, gp_cache
        )
        return(gp_update_loglik_1d(candidate, length(y), collapse_noise))
      }
      u_prop <- u
      u_prop[i] <- value
      candidate <- gp_state_1d(
        y, u_prop, Dx_list, logtheta, if (collapse_noise) 1 else sigma2_eps, theta_spec,
        kernel = kernel, matern_nu = matern_nu,
        gp_cache = gp_cache, require_factor = FALSE
      )
      gp_update_loglik_1d(candidate, length(y), collapse_noise)
    }
    
    z_lower <- pnorm(lower_u)
    z_upper <- pnorm(upper_u)
    
    ## Choose the coordinate system using only the fixed category bounds,
    ## not the current latent value. Infinite/tail/narrow intervals use an
    ## ordinary u-scale slice including the normal prior. This retains the
    ## whole interval instead of clipping its CDF tails or projecting u.
    probability_scale <- is.finite(lower_u) && is.finite(upper_u) &&
      z_lower > 1e-8 && z_upper < 1 - 1e-8 && z_upper - z_lower > 1e-8
    if (!probability_scale) {
      logf_u <- function(value) {
        if (!is.finite(value) || value <= lower_u || value > upper_u) return(-Inf)
        local_gp_loglik(value) + dnorm(value, log = TRUE)
      }
      width <- if (is.finite(lower_u) && is.finite(upper_u)) {
        min(1, (upper_u - lower_u) / 2)
      } else 1
      ans <- bounded_slice_update_1d(
        x0 = u[i], logf = logf_u, w = width,
        lower = lower_u, upper = upper_u, max_steps_out = 30L,
        max_iter = 200L, fail_on_limit = TRUE
      )
      n_eval_total <- n_eval_total + ans$n_eval
      u[i] <- ans$x
      next
    }
    z0 <- pnorm(u[i])
    
    logf_z <- function(z) {
      if (z <= z_lower || z > z_upper) {
        return(-Inf)
      }
      
      u_prop_i <- qnorm(z)
      
      if (!(u_prop_i > lower_u && u_prop_i <= upper_u)) {
        return(-Inf)
      }
      
      local_gp_loglik(u_prop_i)
    }
    
    ans <- full_interval_slice_update(
      x0 = z0,
      logf = logf_z,
      lower = z_lower,
      upper = z_upper,
      max_iter = 200L,
      fail_on_limit = TRUE
    )
    
    n_eval_total <- n_eval_total + ans$n_eval
    u[i] <- qnorm(ans$x)
  }
  
  list(u = u, n_eval = n_eval_total)
}


update_logtheta_slice_1d <- function(y, u, Dx_list, logtheta, sigma2_eps,
                                         theta_slice_width, theta_spec,
                                         kernel = "se", matern_nu = 2.5,
                                         gp_cache = NULL,
                                         collapse_noise = FALSE) {
  n_eval_total <- 0L
  
  for (j in seq_along(logtheta)) {
    logf_j <- function(val) {
      lt <- logtheta
      lt[j] <- val
      
      theta_logpost_1d(
        y, u, Dx_list, lt, sigma2_eps, theta_spec,
        kernel = kernel, matern_nu = matern_nu, gp_cache = gp_cache,
        collapse_noise = collapse_noise
      )
    }
    
    ans <- bounded_slice_update_1d(
      x0 = logtheta[j],
      logf = logf_j,
      w = theta_slice_width[j],
      lower = theta_spec$bounds[j, 1],
      upper = theta_spec$bounds[j, 2],
      max_steps_out = 30L,
      max_iter = 100L,
      fail_on_limit = TRUE
    )
    
    logtheta[j] <- ans$x
    n_eval_total <- n_eval_total + ans$n_eval
  }
  
  list(logtheta = logtheta, n_eval = n_eval_total)
}


split_rhat_1d <- function(chain_list) {
  chain_list <- lapply(chain_list, function(x) as.numeric(x[is.finite(x)]))
  
  lens <- vapply(chain_list, length, integer(1))
  n0 <- min(lens)
  
  if (length(chain_list) < 2 || n0 < 20) {
    return(NA_real_)
  }
  
  chain_list <- lapply(chain_list, function(x) tail(x, n0))
  
  n_half <- floor(n0 / 2)
  
  if (n_half < 10) {
    return(NA_real_)
  }
  
  split_mat <- do.call(
    rbind,
    lapply(chain_list, function(x) {
      rbind(
        x[seq_len(n_half)],
        x[(n0 - n_half + 1):n0]
      )
    })
  )
  
  n_split <- ncol(split_mat)
  
  chain_means <- rowMeans(split_mat)
  chain_vars <- apply(split_mat, 1, var)
  
  W <- mean(chain_vars)
  B <- n_split * var(chain_means)
  
  if (!is.finite(W) || W <= 0) {
    return(NA_real_)
  }
  
  var_hat <- ((n_split - 1) / n_split) * W + B / n_split
  
  sqrt(var_hat / W)
}


aligned_chain_matrix <- function(chain_list) {
  chain_list <- lapply(chain_list, as.numeric)
  lens <- vapply(chain_list, length, integer(1))
  n0 <- min(lens)
  if (length(chain_list) < 2L || n0 < 4L) return(NULL)
  do.call(cbind, lapply(chain_list, function(x) tail(x, n0)))
}


rank_rhat_1d <- function(chain_list) {
  mat <- aligned_chain_matrix(chain_list)
  if (is.null(mat)) return(NA_real_)
  if (requireNamespace("posterior", quietly = TRUE)) {
    return(as.numeric(posterior::rhat(mat)))
  }
  split_rhat_1d(chain_list)
}


bulk_ess_1d <- function(chain_list) {
  mat <- aligned_chain_matrix(chain_list)
  if (is.null(mat)) return(NA_real_)
  if (requireNamespace("posterior", quietly = TRUE)) {
    return(as.numeric(posterior::ess_bulk(mat)))
  }
  sum(vapply(chain_list, ess_ips_1d, numeric(1)))
}


tail_ess <- function(chain_list) {
  mat <- aligned_chain_matrix(chain_list)
  if (is.null(mat)) return(NA_real_)
  if (requireNamespace("posterior", quietly = TRUE)) {
    return(as.numeric(posterior::ess_tail(mat)))
  }
  NA_real_
}


ess_ips_1d <- function(x, max_lag = NULL) {
  x <- as.numeric(x)
  x <- x[is.finite(x)]
  
  n <- length(x)
  
  if (n < 5) {
    return(NA_real_)
  }
  
  sx <- sd(x)
  
  if (!is.finite(sx)) {
    return(NA_real_)
  }
  
  if (sx == 0) {
    return(n)
  }
  
  if (is.null(max_lag)) {
    max_lag <- min(n - 1L, 1000L)
  } else {
    max_lag <- min(n - 1L, as.integer(max_lag))
  }
  
  if (max_lag < 1L) {
    return(n)
  }
  
  ac <- tryCatch(
    as.numeric(stats::acf(
      x,
      lag.max = max_lag,
      plot = FALSE,
      demean = TRUE
    )$acf),
    error = function(e) NA_real_
  )
  
  if (length(ac) <= 1 || all(!is.finite(ac))) {
    return(NA_real_)
  }
  
  ac <- ac[-1]
  ac <- ac[is.finite(ac)]
  
  if (length(ac) == 0) {
    return(n)
  }
  
  mm <- floor(length(ac) / 2)
  
  if (mm >= 1) {
    pair_sums <- ac[2 * seq_len(mm) - 1] + ac[2 * seq_len(mm)]
    first_nonpos <- which(pair_sums <= 0)[1]
    
    if (is.na(first_nonpos)) {
      use_lag <- 2 * mm
    } else {
      use_lag <- 2 * (first_nonpos - 1)
    }
  } else {
    use_lag <- ifelse(ac[1] > 0, 1, 0)
  }
  
  if (use_lag > 0) {
    tau_int <- 1 + 2 * sum(ac[seq_len(use_lag)])
  } else {
    tau_int <- 1
  }
  
  tau_int <- max(tau_int, 1)
  ess <- n / tau_int
  min(max(ess, 1), n)
}


fit_eivgp_1d <- function(...) {
  if (!exists("fit_eivgp_1d_v030", mode = "function")) {
    stop("Install and load the current eivGP package before fitting.")
  }
  fit_eivgp_1d_v030(...)
}


gp_mle_fit_1d <- function(X, y, kernel = c("se", "matern"), matern_nu = 2.5,
                       n_starts = 5L, seed = 1L) {
  X <- as_study1_numeric_matrix(X, "X", nrow_expected = length(y))
  if (any(!is.finite(y)) || length(y) < 2L) {
    stop("y must contain at least two finite observations.")
  }
  y <- as.numeric(y)
  n <- nrow(X)
  d <- ncol(X)
  kernel_spec <- normalize_gp_kernel_1d(kernel, matern_nu)
  
  Dlist <- lapply(seq_len(d), function(j) {
    pairwise_sqdist(X[, j, drop = FALSE])
  })
  
  nll <- function(par) {
    log_sigma2 <- par[1]
    log_rho <- par[2]
    log_theta <- par[-c(1, 2)]
    
    sigma2 <- exp(log_sigma2)
    rho <- exp(log_rho)
    theta <- exp(log_theta)
    
    Rexp <- matrix(0, n, n)
    
    for (j in seq_len(d)) {
      Rexp <- Rexp + theta[j] * Dlist[[j]]
    }
    
    R <- kernel_from_weighted_sqdist_1d(
      Rexp, kernel_spec$name, kernel_spec$matern_nu
    )
    A <- rho^2 * R + diag(n)
    
    U <- try(safe_chol_1d(A), silent = TRUE)
    
    if (inherits(U, "try-error")) {
      return(1e20)
    }
    
    Ainv_y <- solve_chol(U, y)
    logdetA <- 2 * sum(log(diag(U)))
    quad <- sum(y * Ainv_y)
    
    0.5 * (n * log(2 * base::pi * sigma2) + logdetA + quad / sigma2)
  }
  
  lower <- c(log(1e-5), log(0.05), rep(log(1e-4), d))
  upper <- c(log(5), log(100), rep(log(100), d))
  n_starts <- as.integer(n_starts)
  if (length(n_starts) != 1L || is.na(n_starts) || n_starts < 1L) {
    stop("n_starts must be one positive integer.")
  }
  set.seed(as.integer(seed))
  initial <- matrix(NA_real_, nrow = n_starts, ncol = length(lower))
  initial[1L, ] <- c(log(0.05), log(3), rep(log(0.5), d))
  if (n_starts > 1L) {
    for (ss in 2:n_starts) {
      initial[ss, ] <- stats::runif(length(lower), lower, upper)
    }
  }

  attempts <- lapply(seq_len(n_starts), function(ss) {
    tryCatch(
      optim(
        par = initial[ss, ], fn = nll, method = "L-BFGS-B",
        lower = lower, upper = upper, control = list(maxit = 500)
      ),
      error = function(e) e
    )
  })
  report <- data.frame(
    start = seq_len(n_starts), convergence = NA_integer_, objective = Inf,
    message = "", stringsAsFactors = FALSE
  )
  for (ss in seq_along(attempts)) {
    if (inherits(attempts[[ss]], "error")) {
      report$message[ss] <- conditionMessage(attempts[[ss]])
    } else {
      report$convergence[ss] <- attempts[[ss]]$convergence
      report$objective[ss] <- attempts[[ss]]$value
      report$message[ss] <- if (is.null(attempts[[ss]]$message)) {
        ""
      } else {
        attempts[[ss]]$message
      }
    }
  }
  valid <- which(
    report$convergence == 0L & is.finite(report$objective)
  )
  if (length(valid) == 0L) {
    stop("Internal GP MLE did not converge from any prespecified start.")
  }
  best <- valid[which.min(report$objective[valid])]
  opt <- attempts[[best]]
  
  list(
    par = opt$par,
    value = opt$value,
    convergence = opt$convergence,
    optimizer_attempts = report,
    selected_start = best,
    X = X,
    y = y,
    kernel = kernel_spec
  )
}


gp_mle_fit_learned_embedding <- function(x, c_ord, y, m,
                                         n_starts = 8L,
                                         kernel = c("se", "matern"),
                                         matern_nu = 2.5) {
  n <- length(y)
  x <- as_study1_numeric_matrix(x, "x", nrow_expected = n)
  p_x <- ncol(x)
  kernel_spec <- normalize_gp_kernel_1d(kernel, matern_nu)
  Dx_list <- lapply(seq_len(p_x), function(j) {
    pairwise_sqdist(x[, j, drop = FALSE])
  })
  idx_x <- 2L + seq_len(p_x)
  idx_z <- p_x + 3L
  idx_a <- if (m > 2L) seq.int(p_x + 4L, p_x + m + 1L) else integer(0)
  
  nll <- function(par) {
    log_sigma2 <- par[1]
    log_rho <- par[2]
    log_theta_x <- par[idx_x]
    log_theta_z <- par[idx_z]

    a <- par[idx_a]
    
    sigma2 <- exp(log_sigma2)
    rho <- exp(log_rho)
    theta_x <- exp(log_theta_x)
    theta_z <- exp(log_theta_z)
    
    z_scores <- make_monotone_scores_1d(a, m)
    z <- z_scores[c_ord]
    
    Dz <- pairwise_sqdist(matrix(z, ncol = 1))
    
    D2 <- weighted_distance_multix(Dx_list, theta_x) + theta_z * Dz
    R <- kernel_from_weighted_sqdist_1d(
      D2, kernel_spec$name, kernel_spec$matern_nu
    )
    A <- rho^2 * R + diag(n)
    
    U <- try(safe_chol_1d(A), silent = TRUE)
    
    if (inherits(U, "try-error")) {
      return(1e20)
    }
    
    Ainv_y <- solve_chol(U, y)
    logdetA <- 2 * sum(log(diag(U)))
    quad <- sum(y * Ainv_y)
    
    0.5 * (n * log(2 * base::pi * sigma2) + logdetA + quad / sigma2)
  }
  
  lower <- c(
    log(1e-5),
    log(0.05),
    rep(log(1e-4), p_x),
    log(1e-4),
    rep(-6, max(m - 2, 0))
  )
  
  upper <- c(
    log(5),
    log(100),
    rep(log(100), p_x),
    log(100),
    rep(6, max(m - 2, 0))
  )
  
  make_init <- function() {
    c(
      log(0.05),
      log(3),
      rep(log(0.5), p_x),
      log(0.5),
      rnorm(max(m - 2, 0), 0, 0.5)
    )
  }
  
  starts <- replicate(n_starts, make_init(), simplify = FALSE)
  
  opts <- lapply(starts, function(init) {
    optim(
      par = init,
      fn = nll,
      method = "L-BFGS-B",
      lower = lower,
      upper = upper,
      control = list(maxit = 500)
    )
  })
  
  best_id <- which.min(vapply(opts, function(z) z$value, numeric(1)))
  opt <- opts[[best_id]]
  
  z_scores <- make_monotone_scores_1d(opt$par[idx_a], m)
  
  list(
    par = opt$par[seq_len(idx_z)],
    a = opt$par[idx_a],
    z_scores = z_scores,
    value = opt$value,
    convergence = opt$convergence,
    x = x,
    c_ord = c_ord,
    y = y,
    m = m,
    kernel = kernel_spec
  )
}


fit_embedding_baselines <- function(x_raw, y_raw, c_ord, m,
                                    n_starts_learned = 8L,
                                    kernel = c("se", "matern"),
                                    matern_nu = 2.5) {
  x_raw <- as_study1_numeric_matrix(
    x_raw, "x_raw", nrow_expected = length(y_raw)
  )
  kernel_spec <- normalize_gp_kernel_1d(kernel, matern_nu)
  x_center <- colMeans(x_raw)
  x_scale <- apply(x_raw, 2, sd)
  if (any(!is.finite(x_scale) | x_scale <= 0)) {
    stop("Every x_raw column must have positive finite variation")
  }
  x <- sweep(sweep(x_raw, 2, x_center, "-"), 2, x_scale, "/")
  
  y_center <- mean(y_raw)
  y_scale <- sd(y_raw)
  if (!is.finite(y_scale) || y_scale <= 0) {
    stop("y_raw must have positive finite variation.")
  }
  y <- as.numeric((y_raw - y_center) / y_scale)
  
  z_gauss <- qnorm(c_ord / (m + 1))
  
  fit_gauss <- gp_mle_fit_1d(
    X = cbind(x, z_gauss),
    y = y,
    kernel = kernel_spec$name,
    matern_nu = kernel_spec$matern_nu
  )
  
  z_cm_scores <- conditional_mean_scores(c_ord, m)
  z_cm <- z_cm_scores[c_ord]
  
  fit_cm <- gp_mle_fit_1d(
    X = cbind(x, z_cm),
    y = y,
    kernel = kernel_spec$name,
    matern_nu = kernel_spec$matern_nu
  )
  
  fit_learned <- gp_mle_fit_learned_embedding(
    x = x,
    c_ord = c_ord,
    y = y,
    m = m,
    n_starts = n_starts_learned,
    kernel = kernel_spec$name,
    matern_nu = kernel_spec$matern_nu
  )
  
  list(
    x_center = x_center,
    x_scale = x_scale,
    y_center = y_center,
    y_scale = y_scale,
    predictor_names = colnames(x_raw),
    kernel = kernel_spec,
    fit_gauss = fit_gauss,
    fit_cm = fit_cm,
    fit_learned = fit_learned,
    z_cm_scores = z_cm_scores
  )
}


crps_sample_one <- function(draws, y) {
  draws <- sort(as.numeric(draws))
  S <- length(draws)
  
  term1 <- mean(abs(draws - y))
  
  weights <- 2 * seq_len(S) - S - 1
  mean_abs_pair <- 2 * sum(weights * draws) / S^2
  
  term1 - 0.5 * mean_abs_pair
}


crps_sample_matrix <- function(draw_mat, y) {
  vapply(seq_along(y), function(j) {
    crps_sample_one(draw_mat[, j], y[j])
  }, numeric(1))
}


normal_mixture_nlpd <- function(y_true,
                                component_means,
                                component_vars) {
  y_true <- as.numeric(y_true)
  if (is.null(component_means) || is.null(component_vars)) {
    missing <- c(
      if (is.null(component_means)) "conditional_means" else NULL,
      if (is.null(component_vars)) "conditional_vars" else NULL
    )
    return(list(
      value = NA_real_,
      pointwise = rep(NA_real_, length(y_true)),
      reason = paste0(
        "predictive normal-mixture components unavailable: ",
        paste(missing, collapse = " and ")
      )
    ))
  }

  component_means <- as.matrix(component_means)
  component_vars <- as.matrix(component_vars)
  if (ncol(component_means) != length(y_true) ||
      ncol(component_vars) != length(y_true)) {
    stop("Predictive mixture components must have one column per test case.")
  }
  n_component <- max(nrow(component_means), nrow(component_vars))
  if (nrow(component_means) == 1L && n_component > 1L) {
    component_means <- component_means[
      rep(1L, n_component), , drop = FALSE
    ]
  }
  if (nrow(component_vars) == 1L && n_component > 1L) {
    component_vars <- component_vars[
      rep(1L, n_component), , drop = FALSE
    ]
  }
  if (nrow(component_means) != n_component ||
      nrow(component_vars) != n_component || n_component < 1L) {
    stop("Predictive component means and variances have incompatible rows.")
  }
  if (any(!is.finite(component_means)) ||
      any(!is.finite(component_vars)) || any(component_vars < 0)) {
    stop("Predictive normal-mixture components must have finite, nonnegative variances.")
  }
  if (any(component_vars == 0)) {
    return(list(
      value = NA_real_,
      pointwise = rep(NA_real_, length(y_true)),
      reason = paste0(
        "continuous NLPD undefined: predictive mixture contains ",
        "zero-variance components"
      )
    ))
  }

  pointwise <- vapply(seq_along(y_true), function(j) {
    log_component <- stats::dnorm(
      y_true[j],
      mean = component_means[, j],
      sd = sqrt(component_vars[, j]),
      log = TRUE
    )
    largest <- max(log_component)
    if (is.infinite(largest) && largest < 0) return(Inf)
    -(largest + log(mean(exp(log_component - largest))))
  }, numeric(1L))
  list(value = mean(pointwise), pointwise = pointwise, reason = "")
}


summarize_mean_recovery_1d <- function(draw_mat,
                                       m_true,
                                       method,
                                       rep_id,
                                       n_calib,
                                       scenario,
                                       valid_function_draws = FALSE) {
  draw_mat <- as.matrix(draw_mat)
  m_true <- as.numeric(m_true)
  if (ncol(draw_mat) != length(m_true) || any(!is.finite(m_true))) {
    stop("draw_mat and m_true are incompatible in summarize_mean_recovery_1d().")
  }
  conditional_means <- attr(draw_mat, "conditional_means")
  post_mean <- if (!is.null(conditional_means)) {
    colMeans(conditional_means)
  } else {
    colMeans(draw_mat)
  }
  lo <- hi <- rep(NA_real_, length(m_true))
  if (isTRUE(valid_function_draws)) {
    qs <- apply(
      draw_mat, 2L, stats::quantile,
      probs = c(0.025, 0.975), names = FALSE
    )
    lo <- qs[1L, ]
    hi <- qs[2L, ]
  }
  error <- post_mean - m_true
  data.frame(
    rep = rep_id,
    scenario = scenario,
    n_calib = n_calib,
    method = method,
    n_eval = length(m_true),
    RMSE = sqrt(mean(error^2)),
    MAE = mean(abs(error)),
    Bias = mean(error),
    Coverage95 = if (isTRUE(valid_function_draws)) {
      mean(m_true >= lo & m_true <= hi)
    } else {
      NA_real_
    },
    Width95 = if (isTRUE(valid_function_draws)) mean(hi - lo) else NA_real_,
    uncertainty_target = if (isTRUE(valid_function_draws)) "m(x,c)" else "none",
    stringsAsFactors = FALSE
  )
}


sample_oracle_y_star <- function(x_star, c_star, tau_true,
                                 scenario = "active",
                                 sigma_eps = 0.1,
                                 n_draw = 1000,
                                 heterogeneity_eta = 1) {
  lower <- c(-Inf, tau_true)[c_star]
  upper <- c(tau_true, Inf)[c_star]
  
  u_star <- rtruncnorm_vec(
    mean = rep(0, n_draw),
    sd = 1,
    lower = rep(lower, n_draw),
    upper = rep(upper, n_draw)
  )
  
  f_star <- f0_1d(
    x_star,
    u_star,
    scenario = scenario,
    c_ord = rep(c_star, n_draw),
    tau = tau_true,
    heterogeneity_eta = heterogeneity_eta
  )
  
  out <- f_star + rnorm(n_draw, 0, sigma_eps)
  ## Conditional on each Monte Carlo draw of U*, the DGM predictive law is
  ## exactly Gaussian.  Retain those components so that log predictive density
  ## can be evaluated without applying a KDE to the simulated Y* draws.
  attr(out, "conditional_means") <- matrix(f_star, ncol = 1L)
  attr(out, "conditional_vars") <- matrix(
    rep(sigma_eps^2, n_draw), ncol = 1L
  )
  out
}


sample_oracle_test_y <- function(x_test,
                                 c_test,
                                 tau_true,
                                 scenario,
                                 sigma_eps,
                                 n_draw = 1000,
                                 heterogeneity_eta = 1) {
  n_test <- length(c_test)
  out <- matrix(NA_real_, nrow = n_draw, ncol = n_test)
  conditional_means <- matrix(NA_real_, nrow = n_draw, ncol = n_test)
  conditional_vars <- matrix(NA_real_, nrow = n_draw, ncol = n_test)
  
  for (j in seq_len(n_test)) {
    draw_j <- sample_oracle_y_star(
      x_star = x_test[j],
      c_star = c_test[j],
      tau_true = tau_true,
      scenario = scenario,
      sigma_eps = sigma_eps,
      n_draw = n_draw,
      heterogeneity_eta = heterogeneity_eta
    )
    out[, j] <- draw_j
    conditional_means[, j] <- attr(draw_j, "conditional_means")[, 1L]
    conditional_vars[, j] <- attr(draw_j, "conditional_vars")[, 1L]
  }

  attr(out, "conditional_means") <- conditional_means
  attr(out, "conditional_vars") <- conditional_vars
  attr(out, "mixture_components") <- "DGM conditional on Monte Carlo U*"
  out
}


sample_eiv_test_y <- function(x_test_raw,
                              c_test,
                              fit_obj,
                              draw_ids,
                              n_per_draw = 1L,
                              u_test_obs = NULL,
                              u_input_scale = c("raw", "model"),
                              joint = FALSE) {
  if (!is.null(mixedgp_budget()) && !isTRUE(getOption("eivGP.prediction_block",FALSE))) {
    args <- lapply(as.list(match.call(expand.dots=TRUE))[-1L],eval,envir=parent.frame())
    return(mixedgp_predict_blocks(args,FUN=sample_eiv_test_y,object_name="fit_obj"))
  }

  u_input_scale <- match.arg(u_input_scale)
  samples_u <- fit_obj$mcmc$samples_u
  samples_tau <- fit_obj$mcmc$samples_tau
  samples_logtheta <- fit_obj$mcmc$samples_logtheta
  samples_sigma2 <- fit_obj$mcmc$samples_sigma2
  
  x_train <- fit_obj$data$x
  y_train <- fit_obj$data$y
  y_center <- fit_obj$data$y_center
  y_scale <- fit_obj$data$y_scale
  x_center <- fit_obj$data$x_center
  x_scale <- fit_obj$data$x_scale
  predictor_names <- fit_obj$data$predictor_names
  if (is.null(predictor_names)) predictor_names <- colnames(x_train)
  kernel_spec <- kernel_spec_from_fit_1d(fit_obj)

  x_test_raw <- as_study1_numeric_matrix(
    x_test_raw,
    "x_test_raw",
    ncol_expected = ncol(x_train),
    expected_names = predictor_names
  )
  x_test <- sweep(sweep(x_test_raw, 2, x_center, "-"), 2, x_scale, "/")

  n_test <- nrow(x_test)
  if (!(is.numeric(c_test) || is.integer(c_test)) ||
      length(c_test) != n_test || anyNA(c_test) ||
      any(c_test != as.integer(c_test)) ||
      any(c_test < 1L | c_test > fit_obj$data$m)) {
    stop("c_test must contain one ordinal code 1, ..., m per prediction row.")
  }
  c_test <- as.integer(c_test)
  observed_u <- integer(0)
  u_test_model <- rep(NA_real_, n_test)
  if (!is.null(u_test_obs)) {
    u_test_obs <- as.numeric(u_test_obs)
    if (length(u_test_obs) != n_test || any(is.infinite(u_test_obs))) {
      stop("u_test_obs must contain one finite value or NA per prediction row.")
    }
    observed_u <- which(is.finite(u_test_obs))
    u_center <- if (is.null(fit_obj$data$U_center)) 0 else fit_obj$data$U_center
    u_scale <- if (is.null(fit_obj$data$U_scale)) 1 else fit_obj$data$U_scale
    u_test_model[observed_u] <- if (u_input_scale == "raw") {
      (u_test_obs[observed_u] - u_center) / u_scale
    } else {
      u_test_obs[observed_u]
    }
  }
  draw_ids <- as.integer(draw_ids)
  if (length(draw_ids) < 1L || anyNA(draw_ids) ||
      any(draw_ids < 1L | draw_ids > nrow(samples_u))) {
    stop("draw_ids must index saved posterior draws.")
  }
  n_per_draw <- as.integer(n_per_draw)
  if (length(n_per_draw) != 1L || is.na(n_per_draw) || n_per_draw < 1L) {
    stop("n_per_draw must be one positive integer.")
  }
  n_draw <- length(draw_ids) * n_per_draw
  
  out <- matrix(NA_real_, nrow = n_draw, ncol = n_test)
  conditional_means <- matrix(NA_real_, nrow = n_draw, ncol = n_test)
  conditional_vars <- matrix(NA_real_, nrow = n_draw, ncol = n_test)
  
  row_id <- 0L
  
  for (s in draw_ids) {
    for (rr in seq_len(n_per_draw)) {
      row_id <- row_id + 1L
      
      tau_s <- samples_tau[s, ]
      
      lower <- c(-Inf, tau_s)[c_test]
      upper <- c(tau_s, Inf)[c_test]
      
      u_star <- rtruncnorm_vec(
        mean = rep(0, n_test),
        sd = 1,
        lower = lower,
        upper = upper
      )
      if (length(observed_u) > 0L) u_star[observed_u] <- u_test_model[observed_u]
      
      pred <- gp_predict_draw(
        x_train = x_train,
        u_train = samples_u[s, ],
        y_train = y_train,
        x_star = x_test,
        u_star = u_star,
        logtheta = samples_logtheta[s, ],
        sigma2_eps = samples_sigma2[s],
        noisy = TRUE,
        kernel = kernel_spec$name,
        matern_nu = kernel_spec$matern_nu,
        return_cov = isTRUE(joint)
      )

      y_std <- if (isTRUE(joint)) {
        as.numeric(rmvnorm_psd(1L, pred$mean, pred$cov))
      } else {
        pred$mean + sqrt(pred$var) * rnorm(n_test)
      }
      out[row_id, ] <- y_center + y_scale * y_std
      conditional_means[row_id, ] <- y_center + y_scale * pred$mean
      conditional_vars[row_id, ] <- y_scale^2 * pred$var
    }
  }

  attr(out, "joint") <- isTRUE(joint)
  attr(out, "latent_input_scale") <- u_input_scale
  attr(out, "conditional_means") <- conditional_means
  attr(out, "conditional_vars") <- conditional_vars
  attr(out, "mixture_components") <-
    "posterior and latent-input Monte Carlo Gaussian components"
  out
}


sample_eiv_f_given_xu_1d <- function(x_star_raw,
                                      u_star,
                                      fit_obj,
                                      draw_ids,
                                      n_per_draw = 1L,
                                      include_gp_uncertainty = TRUE,
                                      u_input_scale = c("raw", "model"),
                                      joint = FALSE) {
  if (!is.null(mixedgp_budget()) && !isTRUE(getOption("eivGP.prediction_block",FALSE))) {
    args <- lapply(as.list(match.call(expand.dots=TRUE))[-1L],eval,envir=parent.frame())
    return(mixedgp_predict_blocks(args,FUN=sample_eiv_f_given_xu_1d,object_name="fit_obj"))
  }

  u_input_scale <- match.arg(u_input_scale)
  samples_u <- fit_obj$mcmc$samples_u
  samples_logtheta <- fit_obj$mcmc$samples_logtheta
  samples_sigma2 <- fit_obj$mcmc$samples_sigma2
  x_train <- fit_obj$data$x
  predictor_names <- fit_obj$data$predictor_names
  if (is.null(predictor_names)) predictor_names <- colnames(x_train)
  kernel_spec <- kernel_spec_from_fit_1d(fit_obj)

  x_star_raw <- as_study1_numeric_matrix(
    x_star_raw,
    "x_star_raw",
    ncol_expected = ncol(x_train),
    expected_names = predictor_names
  )
  u_star <- as.numeric(u_star)
  if (length(u_star) != nrow(x_star_raw) || any(!is.finite(u_star))) {
    stop("u_star must contain one finite value per prediction row.")
  }
  u_center <- if (is.null(fit_obj$data$U_center)) 0 else fit_obj$data$U_center
  u_scale <- if (is.null(fit_obj$data$U_scale)) 1 else fit_obj$data$U_scale
  if (u_input_scale == "raw") {
    u_star <- (u_star - u_center) / u_scale
  }
  x_star <- sweep(
    sweep(x_star_raw, 2L, fit_obj$data$x_center, "-"),
    2L, fit_obj$data$x_scale, "/"
  )
  draw_ids <- as.integer(draw_ids)
  if (length(draw_ids) < 1L || anyNA(draw_ids) ||
      any(draw_ids < 1L | draw_ids > nrow(samples_u))) {
    stop("draw_ids must index saved posterior draws.")
  }
  n_per_draw <- as.integer(n_per_draw)
  if (length(n_per_draw) != 1L || is.na(n_per_draw) || n_per_draw < 1L) {
    stop("n_per_draw must be one positive integer.")
  }

  out <- matrix(
    NA_real_,
    nrow = length(draw_ids) * n_per_draw,
    ncol = nrow(x_star)
  )
  row_id <- 0L
  for (s in draw_ids) {
    pred <- gp_predict_draw(
      x_train = x_train,
      u_train = samples_u[s, ],
      y_train = fit_obj$data$y,
      x_star = x_star,
      u_star = u_star,
      logtheta = samples_logtheta[s, ],
      sigma2_eps = samples_sigma2[s],
      noisy = FALSE,
      kernel = kernel_spec$name,
      matern_nu = kernel_spec$matern_nu,
      return_cov = isTRUE(joint) && isTRUE(include_gp_uncertainty)
    )
    for (rr in seq_len(n_per_draw)) {
      row_id <- row_id + 1L
      f_std <- if (isTRUE(include_gp_uncertainty)) {
        if (isTRUE(joint)) {
          as.numeric(rmvnorm_psd(1L, pred$mean, pred$cov))
        } else {
          pred$mean + sqrt(pred$var) * rnorm(nrow(x_star))
        }
      } else {
        pred$mean
      }
      out[row_id, ] <- fit_obj$data$y_center + fit_obj$data$y_scale * f_std
    }
  }
  attr(out, "joint") <- isTRUE(joint)
  attr(out, "include_process_uncertainty") <- isTRUE(include_gp_uncertainty)
  attr(out, "latent_input_scale") <- u_input_scale
  out
}


gp_integrated_mean_state_1d <- function(x_train,
                                        u_train,
                                        y_train,
                                        x_star,
                                        u_mc,
                                        logtheta,
                                        sigma2_eps,
                                        kernel = "se",
                                        matern_nu = 2.5,
                                        return_cov = FALSE) {
  x_train <- as_study1_numeric_matrix(x_train, "x_train")
  x_star <- as_study1_numeric_matrix(
    x_star,
    "x_star",
    ncol_expected = ncol(x_train),
    expected_names = colnames(x_train)
  )
  u_mc <- as.matrix(u_mc)
  n <- nrow(x_train)
  N <- nrow(x_star)
  M <- nrow(u_mc)
  p_x <- ncol(x_train)
  if (length(u_train) != n || length(y_train) != n || ncol(u_mc) != N ||
      M < 1L || any(!is.finite(c(x_train, u_train, y_train, x_star, u_mc)))) {
    stop("Incompatible or nonfinite inputs in gp_integrated_mean_state_1d().")
  }

  kernel_spec <- normalize_gp_kernel_1d(kernel, matern_nu)
  theta_spec <- make_theta_spec_multix(p_x)
  rho <- exp(logtheta[1L])
  theta_x <- exp(logtheta[theta_spec$x_index])
  theta_u <- exp(logtheta[theta_spec$u_index])

  Dx_train <- lapply(seq_len(p_x), function(j) {
    pairwise_sqdist(x_train[, j, drop = FALSE])
  })
  D2_train <- weighted_distance_multix(Dx_train, theta_x) +
    theta_u * pairwise_sqdist(matrix(u_train, ncol = 1L))
  R_train <- kernel_from_weighted_sqdist_1d(
    D2_train, kernel_spec$name, kernel_spec$matern_nu
  )
  C_train <- sigma2_eps * (rho^2 * R_train + diag(n))
  chol_train <- safe_chol_1d(C_train)
  alpha <- solve_chol(chol_train, y_train)

  K_bar <- matrix(NA_real_, N, n)
  prior_diag <- numeric(N)
  for (i in seq_len(N)) {
    x_i <- x_star[rep(i, M), , drop = FALSE]
    u_i <- u_mc[, i]
    Dxi <- lapply(seq_len(p_x), function(j) {
      pairwise_sqdist(x_i[, j, drop = FALSE], x_train[, j, drop = FALSE])
    })
    D2_cross <- weighted_distance_multix(Dxi, theta_x) +
      theta_u * pairwise_sqdist(
        matrix(u_i, ncol = 1L), matrix(u_train, ncol = 1L)
      )
    K_cross <- rho^2 * sigma2_eps * kernel_from_weighted_sqdist_1d(
      D2_cross, kernel_spec$name, kernel_spec$matern_nu
    )
    K_bar[i, ] <- colMeans(K_cross)

    D2_ii <- theta_u * pairwise_sqdist(matrix(u_i, ncol = 1L))
    prior_diag[i] <- mean(
      rho^2 * sigma2_eps * kernel_from_weighted_sqdist_1d(
        D2_ii, kernel_spec$name, kernel_spec$matern_nu
      )
    )
  }

  conditional_mean <- as.numeric(K_bar %*% alpha)
  v <- forwardsolve(t(chol_train), t(K_bar))
  var_tolerance <- 100 * (n + N + M) * .Machine$double.eps *
    max(1, prior_diag)
  if (!isTRUE(return_cov)) {
    conditional_var <- prior_diag - colSums(v^2)
    if (any(conditional_var < -var_tolerance)) {
      stop("Integrated GP conditional variance is materially negative.")
    }
    conditional_var <- pmax(conditional_var, 0)
    return(list(mean = conditional_mean, var = conditional_var))
  }

  prior_cov <- matrix(NA_real_, N, N)
  diag(prior_cov) <- prior_diag
  if (N > 1L) {
    for (i in seq_len(N - 1L)) {
      u_i <- u_mc[, i]
      for (j in (i + 1L):N) {
        u_j <- u_mc[, j]
        D2_ij <- matrix(
          sum(theta_x * (x_star[i, ] - x_star[j, ])^2),
          M,
          M
        ) + theta_u * pairwise_sqdist(
          matrix(u_i, ncol = 1L), matrix(u_j, ncol = 1L)
        )
        value <- mean(
          rho^2 * sigma2_eps * kernel_from_weighted_sqdist_1d(
            D2_ij, kernel_spec$name, kernel_spec$matern_nu
          )
        )
        prior_cov[i, j] <- value
        prior_cov[j, i] <- value
      }
    }
  }
  conditional_cov <- prior_cov - crossprod(v)
  conditional_cov <- 0.5 * (conditional_cov + t(conditional_cov))
  if (any(diag(conditional_cov) < -var_tolerance)) {
    stop("Integrated GP conditional covariance has a materially negative diagonal.")
  }
  diag(conditional_cov) <- pmax(diag(conditional_cov), 0)
  list(
    mean = conditional_mean,
    var = diag(conditional_cov),
    cov = conditional_cov
  )
}


sample_eiv_m_given_xc_1d <- function(x_star_raw,
                                      c_star,
                                      fit_obj,
                                      draw_ids = NULL,
                                      n_per_draw = 1L,
                                      n_latent = 256L,
                                      include_process_uncertainty = TRUE,
                                      joint = FALSE,
                                      return_components = FALSE,
                                      seed = NULL) {
  if (!is.null(mixedgp_budget()) && !isTRUE(getOption("eivGP.prediction_block",FALSE))) {
    args <- lapply(as.list(match.call(expand.dots=TRUE))[-1L],eval,envir=parent.frame())
    return(mixedgp_predict_blocks(args,FUN=sample_eiv_m_given_xc_1d,object_name="fit_obj"))
  }

  if (!is.null(seed)) set.seed(as.integer(seed))
  x_train <- fit_obj$data$x
  predictor_names <- fit_obj$data$predictor_names
  if (is.null(predictor_names)) predictor_names <- colnames(x_train)
  x_star_raw <- as_study1_numeric_matrix(
    x_star_raw,
    "x_star_raw",
    ncol_expected = ncol(x_train),
    expected_names = predictor_names
  )
  N <- nrow(x_star_raw)
  if (!(is.numeric(c_star) || is.integer(c_star)) || length(c_star) != N ||
      anyNA(c_star) || any(c_star != as.integer(c_star)) ||
      any(c_star < 1L | c_star > fit_obj$data$m)) {
    stop("c_star must contain one ordinal code 1, ..., m per evaluation row.")
  }
  c_star <- as.integer(c_star)
  x_star <- sweep(
    sweep(x_star_raw, 2L, fit_obj$data$x_center, "-"),
    2L, fit_obj$data$x_scale, "/"
  )

  n_saved <- nrow(fit_obj$mcmc$samples_u)
  if (is.null(draw_ids)) draw_ids <- seq_len(n_saved)
  draw_ids <- as.integer(draw_ids)
  if (length(draw_ids) < 1L || anyNA(draw_ids) ||
      any(draw_ids < 1L | draw_ids > n_saved)) {
    stop("draw_ids must index saved posterior draws.")
  }
  n_per_draw <- as.integer(n_per_draw)
  n_latent <- as.integer(n_latent)
  if (length(n_per_draw) != 1L || is.na(n_per_draw) || n_per_draw < 1L ||
      length(n_latent) != 1L || is.na(n_latent) || n_latent < 1L) {
    stop("n_per_draw and n_latent must be positive integers.")
  }

  kernel_spec <- kernel_spec_from_fit_1d(fit_obj)
  out <- matrix(NA_real_, length(draw_ids) * n_per_draw, N)
  conditional_means <- matrix(NA_real_, length(draw_ids), N)
  conditional_vars <- matrix(NA_real_, length(draw_ids), N)
  row_id <- 0L
  for (ii in seq_along(draw_ids)) {
    s <- draw_ids[ii]
    tau_s <- fit_obj$mcmc$samples_tau[s, ]
    lower <- c(-Inf, tau_s)[c_star]
    upper <- c(tau_s, Inf)[c_star]
    u_mc <- matrix(
      rtruncnorm_vec(
        mean = 0,
        sd = 1,
        lower = rep(lower, each = n_latent),
        upper = rep(upper, each = n_latent)
      ),
      nrow = n_latent,
      ncol = N
    )
    pred <- gp_integrated_mean_state_1d(
      x_train = x_train,
      u_train = fit_obj$mcmc$samples_u[s, ],
      y_train = fit_obj$data$y,
      x_star = x_star,
      u_mc = u_mc,
      logtheta = fit_obj$mcmc$samples_logtheta[s, ],
      sigma2_eps = fit_obj$mcmc$samples_sigma2[s],
      kernel = kernel_spec$name,
      matern_nu = kernel_spec$matern_nu,
      return_cov = isTRUE(joint) && isTRUE(include_process_uncertainty)
    )
    conditional_means[ii, ] <-
      fit_obj$data$y_center + fit_obj$data$y_scale * pred$mean
    conditional_vars[ii, ] <- fit_obj$data$y_scale^2 * pred$var
    for (rr in seq_len(n_per_draw)) {
      row_id <- row_id + 1L
      value_std <- if (!isTRUE(include_process_uncertainty)) {
        pred$mean
      } else if (isTRUE(joint)) {
        as.numeric(rmvnorm_psd(1L, pred$mean, pred$cov))
      } else {
        pred$mean + sqrt(pred$var) * rnorm(N)
      }
      out[row_id, ] <- fit_obj$data$y_center + fit_obj$data$y_scale * value_std
    }
  }
  attr(out, "target") <- "m(x,c)"
  attr(out, "n_latent") <- n_latent
  attr(out, "joint") <- isTRUE(joint)
  attr(out, "include_process_uncertainty") <- isTRUE(include_process_uncertainty)
  if (isTRUE(return_components)) {
    attr(out, "conditional_means") <- conditional_means
    attr(out, "conditional_vars") <- conditional_vars
  }
  out
}


sample_eiv_u_given_c_1d <- function(c_new,
                                     fit_obj,
                                     draw_ids = NULL,
                                     n_per_draw = 1L,
                                     scale = c("auto", "raw", "model"),
                                     seed = NULL) {
  scale <- match.arg(scale)
  anchored <- isTRUE(fit_obj$data$latent_scale_anchored)
  if (scale == "raw" && !anchored) {
    stop("scale='raw' requires calibration data that anchor latent U.")
  }
  if (scale == "auto") scale <- if (anchored) "raw" else "model"
  if (!is.null(seed)) set.seed(as.integer(seed))
  c_new <- as.integer(c_new)
  if (length(c_new) < 1L || anyNA(c_new) ||
      any(c_new < 1L | c_new > fit_obj$data$m)) {
    stop("c_new contains invalid ordinal codes.")
  }
  n_saved <- nrow(fit_obj$mcmc$samples_tau)
  if (is.null(draw_ids)) draw_ids <- seq_len(n_saved)
  draw_ids <- as.integer(draw_ids)
  if (length(draw_ids) < 1L || anyNA(draw_ids) ||
      any(draw_ids < 1L | draw_ids > n_saved)) {
    stop("draw_ids must index saved posterior draws.")
  }
  n_per_draw <- as.integer(n_per_draw)
  if (length(n_per_draw) != 1L || is.na(n_per_draw) || n_per_draw < 1L) {
    stop("n_per_draw must be a positive integer.")
  }
  out <- array(
    NA_real_,
    dim = c(length(draw_ids) * n_per_draw, length(c_new), 1L),
    dimnames = list(NULL, NULL, fit_obj$data$U_names)
  )
  row_id <- 0L
  for (s in draw_ids) {
    tau_s <- fit_obj$mcmc$samples_tau[s, ]
    lower <- c(-Inf, tau_s)[c_new]
    upper <- c(tau_s, Inf)[c_new]
    for (rr in seq_len(n_per_draw)) {
      row_id <- row_id + 1L
      out[row_id, , 1L] <- rtruncnorm_vec(0, 1, lower, upper)
    }
  }
  if (scale == "raw") {
    center <- if (is.null(fit_obj$data$U_center)) 0 else fit_obj$data$U_center
    scale_value <- if (is.null(fit_obj$data$U_scale)) 1 else fit_obj$data$U_scale
    out[, , 1L] <- center + scale_value * out[, , 1L]
  }
  attr(out, "source") <- "prospective"
  attr(out, "latent_sampler") <- "exact_truncated_normal"
  attr(out, "scale") <- if (!anchored && scale == "model") "working" else scale
  attr(out, "latent_scale_anchored") <- anchored
  out
}



# ---- core/00_study2_functions.R ----
############################################################
## 00_study2_functions.R
##
## General-purpose fully Bayesian ordinal-probit EIV-GP code,
## plus Study II synthetic data utilities.
############################################################


require_study2_reporting_packages <- function(packages, context) {
  missing <- packages[
    !vapply(packages, requireNamespace, logical(1), quietly = TRUE)
  ]
  if (length(missing) > 0L) {
    stop(
      "The ", context, " requires optional package(s): ",
      paste(missing, collapse = ", ")
    )
  }
  invisible(TRUE)
}


############################################################
## Configurations
############################################################

study2_config_settings <- function(config = c("quick", "balanced", "thorough")) {
  config <- match.arg(config)
  
  if (config == "quick") {
    return(list(
      config = config,
      n_test = 150L,
      n_rep = 2L,
      
      rep_n_iter = 800L,
      rep_burn = 300L,
      rep_thin = 1L,
      rep_n_chains = 2L,
      
      mc_n_iter = 700L,
      mc_burn = 300L,
      mc_thin = 1L,
      mc_n_chains = 1L,
      
      n_pred_draw = 120L,
      n_density_draw = 200L,
      predictive_latent_sampler = "minimax_tilting",
      diagnostic_n_new_latent_gibbs = 100L,
      rejection_max_batches = 1000L,
      n_oracle_pool = 50000L,
      n_starts_learned = 2L,
      preset = "fast"
    ))
  }
  
  if (config == "balanced") {
    return(list(
      config = config,
      n_test = 400L,
      n_rep = 10L,
      
      rep_n_iter = 6000L,
      rep_burn = 2000L,
      rep_thin = 1L,
      rep_n_chains = 12L,
      
      mc_n_iter = 2500L,
      mc_burn = 800L,
      mc_thin = 1L,
      mc_n_chains = 12L,
      
      n_pred_draw = 400L,
      n_density_draw = 800L,
      predictive_latent_sampler = "minimax_tilting",
      diagnostic_n_new_latent_gibbs = 100L,
      rejection_max_batches = 1000L,
      n_oracle_pool = 150000L,
      n_starts_learned = 4L,
      preset = "balanced"
    ))
  }
  
  list(
    config = config,
    n_test = 500L,
    n_rep = 50L,
    
    rep_n_iter = 6000L,
    rep_burn = 2000L,
    rep_thin = 1L,
    rep_n_chains = 12L,
    
    mc_n_iter = 5000L,
    mc_burn = 1500L,
    mc_thin = 1L,
    mc_n_chains = 12L,
    
    n_pred_draw = 600L,
    n_density_draw = 1200L,
    predictive_latent_sampler = "minimax_tilting",
    diagnostic_n_new_latent_gibbs = 100L,
    rejection_max_batches = 1000L,
    n_oracle_pool = 300000L,
    n_starts_learned = 6L,
    preset = "thorough"
  )
}


make_default_control_ordprobit <- function(n,
                                           n_mis,
                                           preset = c("fast", "balanced", "thorough"),
                                           d = 2L) {
  preset <- match.arg(preset)
  
  if (preset == "fast") {
    theta_update_every <- 8L
    n_u_blocks_per_iter <- 1L
    u_block_size <- min(n_mis, max(8L, ceiling(sqrt(max(n_mis, 1L)))))
    global_u_every <- 50L
    joint_local_blocks_per_iter <- 1L
    joint_global_every <- 10L
    joint_collapsed_every <- 10L
    marginal_measurement_every <- 5L
  }
  
  if (preset == "balanced") {
    theta_update_every <- 4L
    n_u_blocks_per_iter <- 2L
    u_block_size <- min(n_mis, max(10L, ceiling(1.5 * sqrt(max(n_mis, 1L)))))
    global_u_every <- 25L
    joint_local_blocks_per_iter <- 2L
    joint_global_every <- 5L
    joint_collapsed_every <- 5L
    marginal_measurement_every <- 2L
  }
  
  if (preset == "thorough") {
    theta_update_every <- 2L
    n_u_blocks_per_iter <- 3L
    u_block_size <- min(n_mis, max(12L, ceiling(2.0 * sqrt(max(n_mis, 1L)))))
    global_u_every <- 15L
    joint_local_blocks_per_iter <- 3L
    joint_global_every <- 3L
    joint_collapsed_every <- 3L
    marginal_measurement_every <- 1L
  }
  
  list(
    preset = preset,
    
    ## Main MCMC knobs.
    theta_update_every = theta_update_every,
    n_u_blocks_per_iter = n_u_blocks_per_iter,
    u_block_size = u_block_size,
    global_u_every = global_u_every,

    ## Exact partially collapsed/interwoven transition frequencies.
    joint_local_blocks_per_iter = joint_local_blocks_per_iter,
    joint_theta_every = 1L,
    joint_global_every = joint_global_every,
    joint_collapsed_every = joint_collapsed_every,
    marginal_measurement_every = marginal_measurement_every,
    ## Optional prior-Gaussianized joint loading/cut-point row ESS. Zero
    ## retains the scalar marginal measurement update for paired ablations.
    joint_measurement_every = 0L,
    ## Optional score-conditioned loading/missing-U transport, disabled
    ## unless requested. This is an exact reparameterized posterior move.
    loading_transport_every = 0L,
    ## NULL chooses Schur at n >= 120; TRUE/FALSE explicitly force the route.
    ## Both paths evaluate the identical dense-GP target.
    gp_use_block_schur = NULL,
    
    ## Ordinal-probit update frequencies.
    score_update_every = 1L,
    A_update_every = 1L,
    tau_update_every = 1L,
    
    ## Slice sampler knobs for GP hyperparameters.
    theta_slice_width_init = rep(0.8, 1L + d + 10L), # trimmed inside fitter
    adapt_theta_width = TRUE,
    adapt_every = 100L,
    adapt_window = 500L,
    theta_width_min = 0.15,
    theta_width_max = 2.50,
    
    ## Measurement-model prior and constraints.
    tau_bound = 6,
    s_A = 2.5,
    
    ## Elliptical-slice safeguard.
    max_ess_try = 5000L
  )
}


as_numeric_matrix_strict <- function(x,
                                     name,
                                     nrow_expected = NULL,
                                     ncol_expected = NULL,
                                     allow_na = FALSE) {
  if (is.null(dim(x))) x <- matrix(x, ncol = 1L)

  if (is.data.frame(x)) {
    is_numeric_column <- vapply(
      x,
      function(z) is.numeric(z) || is.integer(z) || is.logical(z),
      logical(1)
    )
    if (any(!is_numeric_column)) {
      stop(
        name,
        " must contain only numeric columns; encode factors before fitting."
      )
    }
  } else if (!(is.numeric(x) || is.integer(x) || is.logical(x))) {
    stop(name, " must be numeric.")
  }

  out <- as.matrix(x)
  storage.mode(out) <- "double"

  if (ncol(out) < 1L) stop(name, " must contain at least one column.")
  mixedgp_validate_column_names(out, name)
  if (!is.null(nrow_expected) && nrow(out) != nrow_expected) {
    stop(name, " must have ", nrow_expected, " rows.")
  }
  if (!is.null(ncol_expected) && ncol(out) != ncol_expected) {
    stop(name, " must have ", ncol_expected, " columns.")
  }
  if (isTRUE(allow_na)) {
    if (any(is.infinite(out))) stop(name, " contains infinite values.")
  } else if (any(!is.finite(out))) {
    stop(name, " contains missing or non-finite values.")
  }

  out
}


pairwise_sqdist <- function(a, b = NULL) {
  a <- as.matrix(a)
  if (is.null(b)) b <- a
  b <- as.matrix(b)
  
  aa <- rowSums(a^2)
  bb <- rowSums(b^2)
  
  out <- outer(aa, bb, "+") - 2 * tcrossprod(a, b)
  pmax(out, 0)
}


solve_chol <- function(U, b) {
  backsolve(U, forwardsolve(t(U), b))
}


rmvnorm_chol <- function(n, mean, Sigma) {
  mean <- as.numeric(mean)
  d <- length(mean)
  U <- safe_chol(Sigma)
  Z <- matrix(rnorm(n * d), n, d)
  sweep(Z %*% U, 2, mean, "+")
}


rmvnorm_psd <- function(n, mean, Sigma, tolerance = NULL) {
  mean <- as.numeric(mean)
  Sigma <- as.matrix(Sigma)
  if (!identical(dim(Sigma), c(length(mean), length(mean)))) {
    stop("Sigma has incompatible dimensions in rmvnorm_psd().")
  }
  Sigma <- 0.5 * (Sigma + t(Sigma))
  eig <- eigen(Sigma, symmetric = TRUE)
  scale <- max(1, max(abs(eig$values)))
  if (is.null(tolerance)) {
    tolerance <- 100 * length(mean) * .Machine$double.eps * scale
  }
  if (min(eig$values) < -tolerance) {
    stop("Predictive covariance is not positive semidefinite within tolerance.")
  }
  values <- pmax(eig$values, 0)
  factor <- diag(sqrt(values), nrow = length(values)) %*% t(eig$vectors)
  Z <- matrix(rnorm(as.integer(n) * length(mean)), as.integer(n), length(mean))
  out <- sweep(Z %*% factor, 2, mean, "+")
  attr(out, "psd_tolerance") <- tolerance
  attr(out, "clipped_eigenvalues") <- sum(eig$values < 0)
  out
}


rmvnorm_rows_common <- function(mean_mat, Sigma) {
  mean_mat <- as.matrix(mean_mat)
  n <- nrow(mean_mat)
  d <- ncol(mean_mat)
  
  U <- safe_chol(Sigma)
  Z <- matrix(rnorm(n * d), n, d)
  
  mean_mat + Z %*% U
}


rtruncnorm_vec <- function(mean, sd, lower, upper) {
  n <- max(length(mean), length(lower), length(upper))
  mean <- rep(mean, length.out = n)
  sd <- rep(sd, length.out = n)
  lower <- rep(lower, length.out = n)
  upper <- rep(upper, length.out = n)

  if (any(!is.finite(mean)) || any(!is.finite(sd)) || any(sd <= 0) ||
      anyNA(lower) || anyNA(upper) || any(lower >= upper)) {
    stop("Invalid parameters in rtruncnorm_vec().")
  }

  a <- (lower - mean) / sd
  b <- (upper - mean) / sd
  z <- numeric(n)

  ## Tail-specific probability scales avoid catastrophic cancellation when
  ## both ordinary CDF values round to zero or one.
  log_uniform_between <- function(log_lo, log_hi, size) {
    ratio <- exp(log_lo - log_hi)
    log_hi + log(ratio + runif(size) * (1 - ratio))
  }

  right <- a >= 0
  if (any(right)) {
    log_lo <- pnorm(b[right], lower.tail = FALSE, log.p = TRUE)
    log_hi <- pnorm(a[right], lower.tail = FALSE, log.p = TRUE)
    log_prob <- log_uniform_between(log_lo, log_hi, sum(right))
    z[right] <- qnorm(log_prob, lower.tail = FALSE, log.p = TRUE)
  }

  left <- !right & b <= 0
  if (any(left)) {
    log_lo <- pnorm(a[left], log.p = TRUE)
    log_hi <- pnorm(b[left], log.p = TRUE)
    log_prob <- log_uniform_between(log_lo, log_hi, sum(left))
    z[left] <- qnorm(log_prob, log.p = TRUE)
  }

  middle <- !right & !left
  if (any(middle)) {
    p_lo <- pnorm(a[middle])
    p_hi <- pnorm(b[middle])
    prob <- p_lo + runif(sum(middle)) * (p_hi - p_lo)
    z[middle] <- qnorm(prob)
  }

  out <- mean + sd * z
  tolerance <- 64 * .Machine$double.eps * pmax(1, abs(lower), abs(upper))
  violates <- !is.finite(out) |
    out < lower - tolerance | out > upper + tolerance
  if (any(violates)) {
    stop("Numerical failure in tail-stable truncated-normal sampling.")
  }
  pmin(pmax(out, lower), upper)
}


rtruncnorm_one <- function(mean, sd, lower = -Inf, upper = Inf) {
  rtruncnorm_vec(mean, sd, lower, upper)[1]
}


maximin_lhs_nd <- function(n,
                           d,
                           lower = -1,
                           upper = 1,
                           n_starts = 50L) {
  n <- as.integer(n)
  d <- as.integer(d)
  n_starts <- as.integer(n_starts)
  if (n < 1L || d < 1L || n_starts < 1L) {
    stop("n, d, and n_starts must be positive integers.")
  }

  random_lhs <- function() {
    X <- matrix(NA_real_, n, d)
    for (j in seq_len(d)) {
      z <- (seq_len(n) - runif(n)) / n
      X[, j] <- sample(z)
    }
    X
  }

  best <- NULL
  best_min_distance <- -Inf
  for (start in seq_len(n_starts)) {
    candidate <- random_lhs()
    min_distance <- if (n == 1L) {
      Inf
    } else {
      min(stats::dist(candidate))
    }
    if (min_distance > best_min_distance) {
      best <- candidate
      best_min_distance <- min_distance
    }
  }

  lower + (upper - lower) * best
}


pattern_key <- function(C) {
  C <- as.matrix(C)
  apply(C, 1, paste, collapse = "_")
}


combine_chain_arrays <- function(arr_list) {
  dims <- dim(arr_list[[1]])
  total <- sum(vapply(arr_list, function(a) dim(a)[1], integer(1)))
  
  out <- array(NA_real_, dim = c(total, dims[-1]))
  
  pos <- 0L
  for (a in arr_list) {
    ns <- dim(a)[1]
    idx <- pos + seq_len(ns)
    
    if (length(dims) == 3L) {
      out[idx, , ] <- a
    } else if (length(dims) == 4L) {
      out[idx, , , ] <- a
    } else if (length(dims) == 2L) {
      out[idx, ] <- a
    } else {
      out[idx] <- a
    }
    
    pos <- pos + ns
  }
  
  out
}


############################################################
## Tau helpers for proxy-specific level counts
############################################################

tau_names_from_mvec <- function(m_vec) {
  unlist(
    lapply(seq_along(m_vec), function(j) {
      paste0("tau[", j, ",", seq_len(m_vec[j] - 1L), "]")
    }),
    use.names = FALSE
  )
}


flatten_tau <- function(tau_list) {
  unlist(tau_list, use.names = FALSE)
}


unflatten_tau <- function(tau_vec, m_vec) {
  tau_list <- vector("list", length(m_vec))
  pos <- 0L
  
  for (j in seq_along(m_vec)) {
    len <- m_vec[j] - 1L
    
    if (len > 0L) {
      tau_list[[j]] <- tau_vec[pos + seq_len(len)]
      pos <- pos + len
    } else {
      tau_list[[j]] <- numeric(0)
    }
  }
  
  tau_list
}


############################################################
## Slice sampler
############################################################

bounded_slice_update <- function(x0, logf, w = 1,
                                 lower = -Inf, upper = Inf,
                                 max_steps_out = 50,
                                 max_iter = 200,
                                 fail_on_limit = FALSE) {
  if (upper <= lower) {
    if (isTRUE(fail_on_limit)) stop("Slice interval is empty.")
    return(list(x = x0, n_eval = 0L))
  }
  
  eps <- 1e-12
  if (is.finite(lower)) x0 <- max(x0, lower + eps)
  if (is.finite(upper)) x0 <- min(x0, upper - eps)
  
  f0 <- logf(x0)
  n_eval <- 1L
  
  if (!is.finite(f0)) {
    if (isTRUE(fail_on_limit)) stop("Current slice state has non-finite density.")
    return(list(x = x0, n_eval = n_eval))
  }
  
  logy <- f0 + log(runif(1))
  
  L <- x0 - runif(1) * w
  R <- L + w
  
  if (is.finite(lower)) L <- max(L, lower)
  if (is.finite(upper)) R <- min(R, upper)
  
  J <- floor(runif(1) * max_steps_out)
  K <- max_steps_out - 1L - J
  
  while (J > 0 && (!is.finite(lower) || L > lower)) {
    fL <- logf(L)
    n_eval <- n_eval + 1L
    
    if (!is.finite(fL) || fL <= logy) break
    
    L <- L - w
    if (is.finite(lower)) L <- max(L, lower)
    J <- J - 1L
  }
  
  while (K > 0 && (!is.finite(upper) || R < upper)) {
    fR <- logf(R)
    n_eval <- n_eval + 1L
    
    if (!is.finite(fR) || fR <= logy) break
    
    R <- R + w
    if (is.finite(upper)) R <- min(R, upper)
    K <- K - 1L
  }
  
  for (iter in seq_len(max_iter)) {
    x1 <- runif(1, L, R)
    f1 <- logf(x1)
    n_eval <- n_eval + 1L
    
    if (is.finite(f1) && f1 >= logy) {
      return(list(x = x1, n_eval = n_eval))
    }
    
    if (x1 < x0) {
      L <- x1
    } else {
      R <- x1
    }
  }
  
  if (isTRUE(fail_on_limit)) {
    stop("Slice shrinkage exceeded max_iter.")
  }
  list(x = x0, n_eval = n_eval)
}


score_ordinal_matrix <- function(C, level_scores) {
  C <- as.matrix(C)
  n <- nrow(C)
  q <- ncol(C)
  
  Z <- matrix(NA_real_, n, q)
  
  for (j in seq_len(q)) {
    Z[, j] <- level_scores[[j]][C[, j]]
  }
  
  Z
}


initialize_tau_ord <- function(C, m_vec, tau_bound = 6) {
  C <- as.matrix(C)
  q <- ncol(C)
  
  tau <- vector("list", q)
  
  for (j in seq_len(q)) {
    m_j <- m_vec[j]
    
    if (m_j <= 1L) {
      tau[[j]] <- numeric(0)
      next
    }
    
    counts <- tabulate(C[, j], nbins = m_j)
    
    ## Add small smoothing to handle empty levels at initialization.
    counts_s <- counts + 0.5
    probs <- cumsum(counts_s)[seq_len(m_j - 1L)] / sum(counts_s)
    probs <- pmin(pmax(probs, 0.02), 0.98)
    
    tau_j <- qnorm(probs)
    tau_j <- pmin(pmax(tau_j, -tau_bound + 1e-3), tau_bound - 1e-3)
    
    tau[[j]] <- tau_j
  }
  
  tau
}


init_U_from_ordinal <- function(C,
                                d,
                                m_vec,
                                U_obs_full = NULL,
                                calib_idx = integer(0)) {
  C <- as.matrix(C)
  n <- nrow(C)
  q <- ncol(C)
  
  level_scores <- empirical_normal_scores(C, m_vec = m_vec)
  Z <- score_ordinal_matrix(C, level_scores)
  
  Zs <- scale(Z)
  Zs[!is.finite(Zs)] <- 0
  
  pc <- prcomp(Zs, center = FALSE, scale. = FALSE)
  
  U_work <- matrix(0, n, d)
  
  k <- min(d, ncol(pc$x))
  
  if (k > 0L) {
    U_work[, seq_len(k)] <- pc$x[, seq_len(k), drop = FALSE]
  }
  
  if (d > k) {
    U_work[, (k + 1L):d] <- matrix(rnorm(n * (d - k)), n, d - k)
  }
  
  U_work <- scale(U_work)
  U_work[!is.finite(U_work)] <- 0
  
  calib_idx <- sort(as.integer(calib_idx))
  
  if (length(calib_idx) >= d + 1L && !is.null(U_obs_full)) {
    X_cal <- cbind(1, U_work[calib_idx, , drop = FALSE])
    Y_cal <- U_obs_full[calib_idx, , drop = FALSE]
    
    ridge <- diag(c(0, rep(1e-4, d)), d + 1L)
    
    Beta <- solve(crossprod(X_cal) + ridge, crossprod(X_cal, Y_cal))
    
    U_work <- cbind(1, U_work) %*% Beta
  }
  
  if (length(calib_idx) > 0L && !is.null(U_obs_full)) {
    U_work[calib_idx, ] <- U_obs_full[calib_idx, ]
  }
  
  U_work
}


initialize_A_ord <- function(C,
                             U,
                             tau,
                             m_vec,
                             ident = c("lower_triangular", "none")) {
  ident <- match.arg(ident)
  
  C <- as.matrix(C)
  U <- as.matrix(U)
  
  n <- nrow(C)
  q <- ncol(C)
  d <- ncol(U)
  
  A <- matrix(0, q, d)
  Z <- matrix(NA_real_, n, q)
  
  for (j in seq_len(q)) {
    tau_j <- tau[[j]]
    
    lower <- c(-Inf, tau_j)[C[, j]]
    upper <- c(tau_j, Inf)[C[, j]]
    
    denom <- pnorm(upper) - pnorm(lower)
    numer <- dnorm(lower) - dnorm(upper)
    
    Z[, j] <- numer / pmax(denom, .Machine$double.eps)
  }
  
  for (j in seq_len(q)) {
    fit <- try(lm.fit(x = cbind(1, U), y = Z[, j]), silent = TRUE)
    
    if (!inherits(fit, "try-error")) {
      A[j, ] <- fit$coefficients[-1]
    } else {
      A[j, ] <- rnorm(d, 0, 0.2)
    }
  }
  
  if (ident == "lower_triangular") {
    if (q < d) {
      stop("lower_triangular identification requires q >= d.")
    }
    
    for (j in seq_len(d)) {
      if (j < d) {
        A[j, (j + 1L):d] <- 0
      }
      A[j, j] <- max(abs(A[j, j]), 0.2)
    }
  }
  
  A[!is.finite(A)] <- 0
  A
}


assert_ordprobit_state <- function(C,
                                   U,
                                   S,
                                   A,
                                   tau,
                                   m_vec,
                                   ident,
                                   calib_idx = integer(0),
                                   U_obs = NULL,
                                   logtheta = NULL,
                                   sigma2_eps = NULL,
                                   tolerance = 1e-10) {
  C <- as.matrix(C)
  U <- as.matrix(U)
  S <- as.matrix(S)
  A <- as.matrix(A)
  n <- nrow(C)
  q <- ncol(C)
  d <- ncol(U)
  if (!identical(dim(S), c(n, q)) || nrow(U) != n ||
      !identical(dim(A), c(q, d)) || length(tau) != q ||
      any(!is.finite(c(U, S, A)))) {
    stop("Ordinal-probit sampler state has incompatible or nonfinite values.")
  }
  for (j in seq_len(q)) {
    tau_j <- as.numeric(tau[[j]])
    if (length(tau_j) != m_vec[j] - 1L || any(!is.finite(tau_j)) ||
        is.unsorted(tau_j, strictly = TRUE)) {
      stop("Ordinal-probit cut points are not finite and strictly ordered.")
    }
    lower <- c(-Inf, tau_j)[C[, j]]
    upper <- c(tau_j, Inf)[C[, j]]
    if (any(S[, j] < lower - tolerance | S[, j] > upper + tolerance)) {
      stop("Ordinal-probit scores left their observed-category rectangles.")
    }
  }
  if (ident == "lower_triangular") {
    for (j in seq_len(d)) {
      if (A[j, j] <= 0) stop("A constrained diagonal loading is not positive.")
      if (j < d && any(abs(A[j, (j + 1L):d]) > tolerance)) {
        stop("A violates the lower-triangular identification constraint.")
      }
    }
  }
  if (length(calib_idx) > 0L) {
    if (is.null(U_obs) ||
        max(abs(U[calib_idx, , drop = FALSE] -
                  U_obs[calib_idx, , drop = FALSE])) > tolerance) {
      stop("A calibrated latent input moved away from its observed value.")
    }
  }
  if (!is.null(logtheta) && any(!is.finite(logtheta))) {
    stop("A GP hyperparameter is nonfinite.")
  }
  if (!is.null(sigma2_eps) &&
      (length(sigma2_eps) != 1L || !is.finite(sigma2_eps) || sigma2_eps <= 0)) {
    stop("The observation-noise variance is not positive and finite.")
  }
  invisible(TRUE)
}


update_tau_ord <- function(tau, S, C, m_vec, tau_bound = 6) {
  tau_new <- tau
  q <- ncol(C)
  
  for (j in seq_len(q)) {
    m_j <- m_vec[j]
    
    if (m_j <= 1L) next
    
    for (r in seq_len(m_j - 1L)) {
      tau_j <- tau_new[[j]]
      
      L <- max(
        c(
          -tau_bound,
          if (r > 1L) tau_j[r - 1L] else -Inf,
          S[C[, j] <= r, j]
        ),
        na.rm = TRUE
      )
      
      U <- min(
        c(
          tau_bound,
          if (r < m_j - 1L) tau_j[r + 1L] else Inf,
          S[C[, j] > r, j]
        ),
        na.rm = TRUE
      )
      
      if (!is.finite(L) || !is.finite(U) || L >= U) {
        stop("Cut-point Gibbs update encountered an empty conditional interval.")
      }
      tau_j[r] <- runif(1, L, U)
      
      tau_new[[j]] <- tau_j
    }
  }
  
  tau_new
}


rmvnorm_single_truncated_last_positive <- function(mean, Sigma) {
  mean <- as.numeric(mean)
  r <- length(mean)
  
  if (r == 1L) {
    return(rtruncnorm_one(mean[1], sqrt(Sigma[1, 1]), lower = 0, upper = Inf))
  }
  
  idx_minus <- seq_len(r - 1L)
  idx_last <- r
  
  mu_l <- mean[idx_last]
  var_l <- Sigma[idx_last, idx_last]
  
  beta_l <- rtruncnorm_one(mu_l, sqrt(var_l), lower = 0, upper = Inf)
  
  mu_m <- mean[idx_minus]
  Sigma_mm <- Sigma[idx_minus, idx_minus, drop = FALSE]
  Sigma_ml <- Sigma[idx_minus, idx_last, drop = FALSE]
  
  cond_mean <- mu_m + as.numeric(Sigma_ml) / var_l * (beta_l - mu_l)
  cond_cov <- Sigma_mm - Sigma_ml %*% t(Sigma_ml) / var_l
  
  beta_m <- as.numeric(rmvnorm_chol(1, cond_mean, cond_cov))
  
  c(beta_m, beta_l)
}


update_A_ord <- function(S,
                         U,
                         s_A = 2.5,
                         ident = c("lower_triangular", "none")) {
  ident <- match.arg(ident)
  
  S <- as.matrix(S)
  U <- as.matrix(U)
  
  q <- ncol(S)
  d <- ncol(U)
  
  if (ident == "lower_triangular" && q < d) {
    stop("lower_triangular identification requires q >= d.")
  }
  
  A <- matrix(0, q, d)
  prior_prec <- 1 / s_A^2
  
  for (j in seq_len(q)) {
    if (ident == "lower_triangular" && j <= d) {
      active <- seq_len(j)
      Xj <- U[, active, drop = FALSE]
      
      Vj <- solve(crossprod(Xj) + prior_prec * diag(length(active)))
      mj <- as.numeric(Vj %*% crossprod(Xj, S[, j]))
      
      beta <- rmvnorm_single_truncated_last_positive(mj, Vj)
      
      A[j, active] <- beta
      
      if (j < d) {
        A[j, (j + 1L):d] <- 0
      }
    } else {
      Xj <- U
      Vj <- solve(crossprod(Xj) + prior_prec * diag(d))
      mj <- as.numeric(Vj %*% crossprod(Xj, S[, j]))
      
      A[j, ] <- rmvnorm_chol(1, mj, Vj)
    }
  }
  
  A
}


latent_reference_params <- function(S, A) {
  S <- as.matrix(S)
  A <- as.matrix(A)
  
  d <- ncol(A)
  
  V <- solve(diag(d) + crossprod(A))
  M <- S %*% A %*% V
  
  list(mean = M, V = V)
}


############################################################
## General GP likelihood and prediction
############################################################

a_eps0 <- 2

b_eps0 <- 0.05


make_gp_prior <- function(p, d) {
  npar <- 1L + p + d
  
  list(
    mean = c(log(3), rep(log(0.5), p + d)),
    sd = c(1.5, rep(1.5, p + d)),
    lower = c(log(0.05), rep(log(1e-4), p + d)),
    upper = c(log(100), rep(log(100), p + d))
  )
}


log_prior_logtheta_gp <- function(logtheta, gp_prior) {
  if (length(logtheta) != length(gp_prior$mean)) return(-Inf)
  
  if (any(logtheta < gp_prior$lower) || any(logtheta > gp_prior$upper)) {
    return(-Inf)
  }
  
  sum(dnorm(
    logtheta,
    mean = gp_prior$mean,
    sd = gp_prior$sd,
    log = TRUE
  ))
}


weighted_sqdist_general <- function(X_a,
                                    U_a,
                                    theta_x,
                                    theta_u,
                                    X_b = NULL,
                                    U_b = NULL) {
  X_a <- as.matrix(X_a)
  U_a <- as.matrix(U_a)
  if (is.null(X_b)) X_b <- X_a
  if (is.null(U_b)) U_b <- U_a
  X_b <- as.matrix(X_b)
  U_b <- as.matrix(U_b)

  if (nrow(X_a) != nrow(U_a) || nrow(X_b) != nrow(U_b)) {
    stop("X and U must have the same number of rows within each data set.")
  }
  if (ncol(X_a) != ncol(X_b) || ncol(U_a) != ncol(U_b)) {
    stop("Training and prediction inputs have incompatible dimensions.")
  }
  if (length(theta_x) != ncol(X_a) || length(theta_u) != ncol(U_a)) {
    stop("Kernel parameter dimensions do not match X and U.")
  }

  D2 <- matrix(0, nrow(X_a), nrow(X_b))
  for (j in seq_len(ncol(X_a))) {
    D2 <- D2 + theta_x[j] * pairwise_sqdist(
      X_a[, j, drop = FALSE],
      X_b[, j, drop = FALSE]
    )
  }
  for (k in seq_len(ncol(U_a))) {
    D2 <- D2 + theta_u[k] * pairwise_sqdist(
      U_a[, k, drop = FALSE],
      U_b[, k, drop = FALSE]
    )
  }

  pmax(D2, 0)
}


gp_corr_general <- function(X,
                            U,
                            logtheta,
                            kernel = "se",
                            matern_nu = 2.5) {
  X <- as.matrix(X)
  U <- as.matrix(U)
  
  p <- ncol(X)
  d <- ncol(U)
  
  expected_len <- 1L + p + d
  
  if (length(logtheta) != expected_len) {
    stop("logtheta has wrong length.")
  }
  
  rho <- exp(logtheta[1])
  theta_x <- exp(logtheta[1L + seq_len(p)])
  theta_u <- exp(logtheta[1L + p + seq_len(d)])
  
  D2 <- weighted_sqdist_general(X, U, theta_x, theta_u)
  kernel_spec <- normalize_gp_kernel(kernel, matern_nu)
  
  list(
    R = kernel_from_weighted_sqdist(
      D2,
      kernel = kernel_spec$name,
      matern_nu = kernel_spec$matern_nu
    ),
    rho = rho,
    theta_x = theta_x,
    theta_u = theta_u,
    kernel = kernel_spec
  )
}


gp_state_general <- function(y,
                             X,
                             U,
                             logtheta,
                             sigma2_eps,
                             kernel = "se",
                             matern_nu = 2.5) {
  n <- length(y)
  
  cc <- gp_corr_general(X, U, logtheta, kernel, matern_nu)
  A_mat <- diag(n) + cc$rho^2 * cc$R
  
  Uchol <- safe_chol(A_mat)
  Ainv_y <- solve_chol(Uchol, y)
  
  logdetA <- 2 * sum(log(diag(Uchol)))
  quad <- sum(y * Ainv_y)
  
  loglik <- -0.5 * (
    n * log(2 * base::pi * sigma2_eps) +
      logdetA +
      quad / sigma2_eps
  )
  
  list(
    loglik = loglik,
    R = cc$R,
    A = A_mat,
    cholA = Uchol,
    Ainv_y = Ainv_y,
    logdetA = logdetA,
    quad = quad
  )
}


theta_logpost_integrated_general <- function(y,
                                             X,
                                             U,
                                             logtheta,
                                             gp_prior,
                                             kernel = "se",
                                             matern_nu = 2.5) {
  lp <- log_prior_logtheta_gp(logtheta, gp_prior)
  if (!is.finite(lp)) return(-Inf)

  lp + gp_loglik_integrated_general(
    y, X, U, logtheta,
    kernel = kernel,
    matern_nu = matern_nu
  )
}


gp_loglik_integrated_general <- function(y,
                                         X,
                                         U,
                                         logtheta,
                                         kernel = "se",
                                         matern_nu = 2.5) {
  n <- length(y)
  
  cc <- gp_corr_general(X, U, logtheta, kernel, matern_nu)
  A_mat <- diag(n) + cc$rho^2 * cc$R
  
  ## A_mat is positive definite under the stated covariance model. A failure
  ## here is therefore a numerical failure, not an out-of-support proposal;
  ## abort rather than silently changing the transition kernel.
  Uchol <- safe_chol(A_mat)
  
  Ainv_y <- solve_chol(Uchol, y)
  
  logdetA <- 2 * sum(log(diag(Uchol)))
  quad <- sum(y * Ainv_y)
  
  -0.5 * logdetA -
    (a_eps0 + n / 2) * log(b_eps0 + 0.5 * quad)
}


gp_loglik_from_moments_general <- function(state, sigma2_eps = NULL) {
  if (is.null(sigma2_eps)) {
    return(-0.5 * state$logdet -
      (a_eps0 + state$n / 2) * log(b_eps0 + 0.5 * state$quad))
  }
  -0.5 * (state$n * log(2 * base::pi * sigma2_eps) +
            state$logdet + state$quad / sigma2_eps)
}


## Per-chain exact likelihood evaluator. The accepted cache contains sufficient
## GP moments; changing scores/loadings/thresholds does not invalidate it.
## A fixed-theta local proposal changes only the selected covariance rows.
make_gp_evaluator_general <- function(y, X, calib_idx = integer(0),
                                      U_obs = NULL, kernel = "se",
                                      matern_nu = 2.5,
                                      use_block_schur = TRUE) {
  n <- length(y)
  accepted <- NULL
  counts <- c(gp_dense_evaluations = 0L, gp_block_evaluations = 0L,
              gp_cache_hits = 0L, gp_full_cholesky = 0L,
              gp_block_setup_cholesky = 0L, gp_schur_cholesky = 0L)
  increment <- function(name) counts[name] <<- counts[name] + 1L
  assert_calibration <- function(U) {
    if (length(calib_idx) &&
        (is.null(U_obs) ||
         any(U[calib_idx, , drop = FALSE] != U_obs[calib_idx, , drop = FALSE]))) {
      stop("A calibrated latent input changed before a target evaluation.")
    }
    invisible(TRUE)
  }
  validate_block <- function(block_idx) {
    if (anyNA(block_idx) || anyDuplicated(block_idx) ||
        any(block_idx < 1L | block_idx > n) ||
        any(block_idx %in% calib_idx)) {
      stop("A latent update block contains an invalid or calibrated row.")
    }
    invisible(TRUE)
  }
  correlation <- function(Xa, Ua, Xb = Xa, Ub = Ua, logtheta) {
    p <- ncol(X)
    kernel_from_weighted_sqdist(
      weighted_sqdist_general(Xa, Ua,
        exp(logtheta[1L + seq_len(p)]),
        exp(logtheta[1L + p + seq_len(ncol(Ua))]), Xb, Ub),
      kernel, matern_nu)
  }
  dense <- function(U, logtheta) {
    increment("gp_dense_evaluations")
    B <- diag(n) + exp(2 * logtheta[1L]) *
      correlation(X, U, logtheta = logtheta)
    increment("gp_full_cholesky")
    R <- safe_chol(B)
    list(n = n, logdet = 2 * sum(log(diag(R))),
         quad = sum(y * solve_chol(R, y)), U = U, logtheta = logtheta)
  }
  prepare_block <- function(U, logtheta, block_idx) {
    assert_calibration(U)
    validate_block(block_idx)
    b <- length(block_idx)
    if (!isTRUE(use_block_schur) || b == 0L || b >= n / 2) return(NULL)
    rest <- setdiff(seq_len(n), block_idx)
    rho2 <- exp(2 * logtheta[1L])
    B_rest <- diag(length(rest)) + rho2 * correlation(
      X[rest, , drop = FALSE], U[rest, , drop = FALSE], logtheta = logtheta)
    increment("gp_block_setup_cholesky")
    R_rest <- safe_chol(B_rest)
    alpha_rest <- solve_chol(R_rest, y[rest])
    list(idx = block_idx, rest = rest, U_rest = U[rest, , drop = FALSE],
         R_rest = R_rest, alpha_rest = alpha_rest,
         quad_rest = sum(y[rest] * alpha_rest),
         logdet_rest = 2 * sum(log(diag(R_rest))),
         logtheta = logtheta, rho2 = rho2)
  }
  evaluate <- function(U, logtheta, block = NULL) {
    assert_calibration(U)
    if (!is.null(accepted) && identical(U, accepted$U) &&
        identical(logtheta, accepted$logtheta)) {
      increment("gp_cache_hits")
      return(accepted)
    }
    if (is.null(block)) return(dense(U, logtheta))
    if (!identical(logtheta, block$logtheta) ||
        !identical(U[block$rest, , drop = FALSE], block$U_rest)) {
      stop("A Schur likelihood cache was used after its fixed inputs changed.")
    }
    increment("gp_block_evaluations")
    idx <- block$idx
    cross <- block$rho2 * correlation(
      X[idx, , drop = FALSE], U[idx, , drop = FALSE],
      X[block$rest, , drop = FALSE], block$U_rest, logtheta)
    diagonal <- diag(length(idx)) + block$rho2 * correlation(
      X[idx, , drop = FALSE], U[idx, , drop = FALSE], logtheta = logtheta)
    W <- forwardsolve(t(block$R_rest), t(cross))
    H <- diagonal - crossprod(W)
    increment("gp_schur_cholesky")
    R_H <- safe_chol(0.5 * (H + t(H)))
    residual <- y[idx] - as.numeric(cross %*% block$alpha_rest)
    list(n = n,
      logdet = block$logdet_rest + 2 * sum(log(diag(R_H))),
      quad = block$quad_rest + sum(residual * solve_chol(R_H, residual)),
      U = U, logtheta = logtheta)
  }
  accept <- function(state) {
    assert_calibration(state$U)
    if (!is.finite(state$logdet) || !is.finite(state$quad) || state$quad < 0) {
      stop("Cannot cache non-finite or invalid accepted GP moments.")
    }
    accepted <<- state
    invisible(NULL)
  }
  list(evaluate = evaluate, accept = accept, prepare_block = prepare_block,
       assert_calibration = assert_calibration, validate_block = validate_block,
       snapshot = function() accepted,
       restore = function(state) { accepted <<- state; invisible(NULL) },
       counts = function() counts)
}


logspace_sub <- function(log_x, log_y) {
  log_x <- as.numeric(log_x)
  log_y <- as.numeric(log_y)
  if (length(log_x) != length(log_y)) {
    stop("log_x and log_y must have the same length.")
  }

  out <- rep(-Inf, length(log_x))
  only_x <- is.finite(log_x) & !is.finite(log_y)
  out[only_x] <- log_x[only_x]
  both <- is.finite(log_x) & is.finite(log_y) & log_y < log_x
  out[both] <- log_x[both] + log(-expm1(log_y[both] - log_x[both]))
  out
}


log_normal_interval_prob <- function(lower, upper) {
  lower <- as.numeric(lower)
  upper <- as.numeric(upper)
  if (length(lower) != length(upper)) {
    stop("lower and upper must have the same length.")
  }

  out <- rep(-Inf, length(lower))
  full <- is.infinite(lower) & lower < 0 & is.infinite(upper) & upper > 0
  out[full] <- 0
  valid <- is.finite(lower) | is.finite(upper)
  valid <- valid & lower < upper
  if (!any(valid)) return(out)

  right <- valid & lower >= 0
  if (any(right)) {
    out[right] <- logspace_sub(
      pnorm(lower[right], lower.tail = FALSE, log.p = TRUE),
      pnorm(upper[right], lower.tail = FALSE, log.p = TRUE)
    )
  }

  left <- valid & lower < 0
  if (any(left)) {
    out[left] <- logspace_sub(
      pnorm(upper[left], log.p = TRUE),
      pnorm(lower[left], log.p = TRUE)
    )
  }

  out
}


ordinal_loglik_marginal <- function(C, U, A, tau) {
  C <- as.matrix(C)
  U <- as.matrix(U)
  A <- as.matrix(A)

  if (nrow(C) != nrow(U) || ncol(C) != nrow(A) || ncol(U) != ncol(A)) {
    stop("Incompatible dimensions in ordinal_loglik_marginal().")
  }

  ans <- 0
  for (j in seq_len(ncol(C))) {
    mu <- as.numeric(U %*% A[j, ])
    tau_j <- tau[[j]]
    lower <- c(-Inf, tau_j)[C[, j]] - mu
    upper <- c(tau_j, Inf)[C[, j]] - mu
    log_prob <- log_normal_interval_prob(lower, upper)
    if (any(!is.finite(log_prob))) return(-Inf)
    ans <- ans + sum(log_prob)
  }

  ans
}


sample_sigma2_eps_general <- function(y,
                                      X,
                                      U,
                                      logtheta,
                                      kernel = "se",
                                      matern_nu = 2.5,
                                      gp_evaluator = NULL) {
  n <- length(y)
  if (!is.null(gp_evaluator)) {
    state <- gp_evaluator$evaluate(U, logtheta)
    gp_evaluator$accept(state)
    return(1 / rgamma(1, shape = a_eps0 + n / 2,
                     rate = b_eps0 + 0.5 * state$quad))
  }
  
  cc <- gp_corr_general(X, U, logtheta, kernel, matern_nu)
  A_mat <- diag(n) + cc$rho^2 * cc$R
  
  Uchol <- safe_chol(A_mat)
  Ainv_y <- solve_chol(Uchol, y)
  
  quad <- sum(y * Ainv_y)
  
  shape <- a_eps0 + n / 2
  rate <- b_eps0 + 0.5 * quad
  
  1 / rgamma(1, shape = shape, rate = rate)
}


update_logtheta_slice_general <- function(y,
                                          X,
                                          U,
                                          logtheta,
                                          theta_slice_width,
                                          gp_prior,
                                          kernel = "se",
                                          matern_nu = 2.5,
                                          gp_evaluator = NULL) {
  n_eval_total <- 0L
  candidate_gp <- NULL
  
  for (j in seq_along(logtheta)) {
    logf_j <- function(val) {
      lt <- logtheta
      lt[j] <- val
      if (!is.null(gp_evaluator)) {
        lp <- log_prior_logtheta_gp(lt, gp_prior)
        if (!is.finite(lp)) return(-Inf)
        candidate_gp <<- gp_evaluator$evaluate(U, lt)
        return(lp + gp_loglik_from_moments_general(candidate_gp))
      }
      theta_logpost_integrated_general(
        y, X, U, lt, gp_prior,
        kernel = kernel,
        matern_nu = matern_nu
      )
    }
    
    ans <- bounded_slice_update(
      x0 = logtheta[j],
      logf = logf_j,
      w = theta_slice_width[j],
      lower = gp_prior$lower[j],
      upper = gp_prior$upper[j],
      max_steps_out = 30L,
      max_iter = 100L,
      fail_on_limit = TRUE
    )
    
    logtheta[j] <- ans$x
    if (!is.null(gp_evaluator)) gp_evaluator$accept(candidate_gp)
    n_eval_total <- n_eval_total + ans$n_eval
  }
  
  list(logtheta = logtheta, n_eval = n_eval_total)
}


############################################################
## Latent U update
############################################################

update_U_ess_block_general <- function(y,
                                       X,
                                       U_curr,
                                       S,
                                       A,
                                       logtheta,
                                       sigma2_eps,
                                       block_idx,
                                       max_try = 200L,
                                       kernel = "se",
                                       matern_nu = 2.5,
                                       gp_evaluator = NULL,
                                       score_reference = NULL) {
  block_idx <- as.integer(block_idx)
  if (!is.null(gp_evaluator)) gp_evaluator$validate_block(block_idx)
  
  if (length(block_idx) == 0L) {
    return(list(U = U_curr, n_eval = 0L, accepted = TRUE))
  }
  
  ref <- if (is.null(score_reference)) latent_reference_params(S, A) else score_reference
  
  Mref <- ref$mean
  Vref <- ref$V
  
  U_block <- U_curr[block_idx, , drop = FALSE]
  M_block <- Mref[block_idx, , drop = FALSE]
  block_cache <- if (!is.null(gp_evaluator)) {
    gp_evaluator$prepare_block(U_curr, logtheta, block_idx)
  } else NULL
  candidate_gp <- NULL
  
  loglik_fun <- function(U_block_prop) {
    U_prop <- U_curr
    U_prop[block_idx, ] <- U_block_prop
    if (!is.null(gp_evaluator)) {
      candidate_gp <<- gp_evaluator$evaluate(U_prop, logtheta, block_cache)
      return(gp_loglik_from_moments_general(candidate_gp, sigma2_eps))
    }
    gp_state_general(
      y, X, U_prop, logtheta, sigma2_eps,
      kernel = kernel,
      matern_nu = matern_nu
    )$loglik
  }
  
  loglik_cur <- loglik_fun(U_block)
  
  if (!is.finite(loglik_cur)) {
    stop("Current state has a non-finite density in latent ESS update.")
  }
  
  Nu <- rmvnorm_rows_common(
    mean_mat = matrix(0, nrow = length(block_idx), ncol = ncol(U_curr)),
    Sigma = Vref
  )
  
  Zcur <- U_block - M_block
  
  logy <- loglik_cur + log(runif(1))
  
  angle <- runif(1, 0, 2 * base::pi)
  angle_min <- angle - 2 * base::pi
  angle_max <- angle
  
  n_eval <- 1L
  
  for (try_id in seq_len(max_try)) {
    U_prop_block <- M_block + Zcur * cos(angle) + Nu * sin(angle)
    
    loglik_prop <- loglik_fun(U_prop_block)
    
    n_eval <- n_eval + 1L
    
    if (is.finite(loglik_prop) && loglik_prop >= logy) {
      U_new <- U_curr
      U_new[block_idx, ] <- U_prop_block
      if (!is.null(gp_evaluator)) gp_evaluator$accept(candidate_gp)
      
      return(list(U = U_new, n_eval = n_eval, accepted = TRUE))
    }
    
    if (angle < 0) {
      angle_min <- angle
    } else {
      angle_max <- angle
    }
    
    angle <- runif(1, angle_min, angle_max)
  }
  
  stop("Latent elliptical-slice update exceeded max_try.")
}


update_U_theta_ess_integrated_general <- function(y,
                                                  X,
                                                  C,
                                                  U_curr,
                                                  logtheta_curr,
                                                  gp_prior,
                                                  block_idx,
                                                  reference = c("score", "prior"),
                                                  S = NULL,
                                                  A = NULL,
                                                  tau = NULL,
                                                  update_theta = TRUE,
                                                  max_try = 5000L,
                                                  kernel = "se",
                                                  matern_nu = 2.5,
                                                  gp_evaluator = NULL,
                                                  score_reference = NULL) {
  reference <- match.arg(reference)
  block_idx <- as.integer(block_idx)
  d <- ncol(U_curr)
  if (!is.null(gp_evaluator)) gp_evaluator$validate_block(block_idx)

  if (reference == "score" && length(block_idx) > 0L) {
    if (is.null(S) || is.null(A)) {
      stop("S and A are required for a score-reference update.")
    }
    ref <- if (is.null(score_reference)) latent_reference_params(S, A) else score_reference
    M_block <- ref$mean[block_idx, , drop = FALSE]
    V_block <- ref$V
  } else if (reference == "prior") {
    if (is.null(A) || is.null(tau)) {
      stop("A and tau are required for a prior-reference update.")
    }
    M_block <- matrix(0, nrow = length(block_idx), ncol = d)
    V_block <- diag(d)
  } else {
    M_block <- matrix(numeric(0), 0L, d)
    V_block <- diag(d)
  }

  U_block <- U_curr[block_idx, , drop = FALSE]
  theta_mean <- if (isTRUE(update_theta)) gp_prior$mean else numeric(0)
  theta_sd <- if (isTRUE(update_theta)) gp_prior$sd else numeric(0)
  block_cache <- if (!is.null(gp_evaluator) && !isTRUE(update_theta)) {
    gp_evaluator$prepare_block(U_curr, logtheta_curr, block_idx)
  } else NULL
  candidate_gp <- NULL

  residual_log_density <- function(U_block_prop, logtheta_prop) {
    if (any(logtheta_prop < gp_prior$lower) ||
        any(logtheta_prop > gp_prior$upper)) {
      return(-Inf)
    }

    U_prop <- U_curr
    if (length(block_idx) > 0L) {
      U_prop[block_idx, ] <- U_block_prop
    }

    out <- if (is.null(gp_evaluator)) {
      gp_loglik_integrated_general(y, X, U_prop, logtheta_prop,
        kernel = kernel, matern_nu = matern_nu)
    } else {
      candidate_gp <<- gp_evaluator$evaluate(U_prop, logtheta_prop, block_cache)
      gp_loglik_from_moments_general(candidate_gp)
    }
    if (!is.finite(out)) return(-Inf)

    if (reference == "prior") {
      out <- out + ordinal_loglik_marginal(C, U_prop, A, tau)
    }

    out
  }

  loglik_cur <- residual_log_density(U_block, logtheta_curr)
  if (!is.finite(loglik_cur)) {
    stop("Current state has a non-finite density in joint ESS update.")
  }

  Nu_U <- if (length(block_idx) > 0L) {
    rmvnorm_rows_common(
      mean_mat = matrix(0, nrow = length(block_idx), ncol = d),
      Sigma = V_block
    )
  } else {
    matrix(numeric(0), nrow = 0L, ncol = d)
  }
  Nu_theta <- if (isTRUE(update_theta)) {
    rnorm(length(logtheta_curr), 0, theta_sd)
  } else {
    numeric(0)
  }

  Z_U <- U_block - M_block
  Z_theta <- if (isTRUE(update_theta)) logtheta_curr - theta_mean else numeric(0)
  logy <- loglik_cur + log(runif(1))

  angle <- runif(1, 0, 2 * base::pi)
  angle_min <- angle - 2 * base::pi
  angle_max <- angle
  n_eval <- 1L

  for (try_id in seq_len(max_try)) {
    U_prop_block <- M_block + Z_U * cos(angle) + Nu_U * sin(angle)
    logtheta_prop <- if (isTRUE(update_theta)) {
      theta_mean + Z_theta * cos(angle) + Nu_theta * sin(angle)
    } else {
      logtheta_curr
    }

    loglik_prop <- residual_log_density(U_prop_block, logtheta_prop)
    n_eval <- n_eval + 1L

    if (is.finite(loglik_prop) && loglik_prop >= logy) {
      U_new <- U_curr
      if (length(block_idx) > 0L) {
        U_new[block_idx, ] <- U_prop_block
      }
      if (!is.null(gp_evaluator)) gp_evaluator$accept(candidate_gp)

      return(list(
        U = U_new,
        logtheta = logtheta_prop,
        n_eval = n_eval,
        accepted = TRUE,
        angle = angle,
        reference = reference
      ))
    }

    if (angle < 0) {
      angle_min <- angle
    } else {
      angle_max <- angle
    }
    angle <- runif(1, angle_min, angle_max)
  }

  stop("Joint elliptical-slice update exceeded max_try.")
}


## Exact Gaussian reference for one measurement row. The original prior is
## Gaussian on free loadings (half-normal on an identified diagonal) and
## uniform on -B < tau_1 < ... < tau_K < B. Independent Beta(1,K-r+1)
## stick fractions generate exactly the latter ordered-uniform distribution.
## No logistic-normal cut-point prior is introduced by this reparameterization.
measurement_row_active <- function(j, d, ident) {
  if (ident == "lower_triangular" && j <= d) seq_len(j) else seq_len(d)
}


measurement_row_to_gaussian <- function(A_row, tau_row, j, ident,
                                         s_A = 2.5, tau_bound = 6) {
  active <- measurement_row_active(j, length(A_row), ident)
  z_A <- A_row[active] / s_A
  if (ident == "lower_triangular" && j <= length(A_row)) {
    if (A_row[j] <= 0) stop("A diagonal must be strictly positive.")
    ## A^2/s_A^2 is chi-square(1). This avoids subtracting two nearly
    ## equal normal probabilities when a positive diagonal is near zero.
    lower_tail <- A_row[j] / s_A < 1
    z_A[j] <- qnorm(pchisq((A_row[j] / s_A)^2, df = 1,
      lower.tail = lower_tail, log.p = TRUE), lower.tail = lower_tail, log.p = TRUE)
  }
  K <- length(tau_row)
  if (K && (any(!is.finite(tau_row)) ||
            any(diff(c(-tau_bound, tau_row, tau_bound)) <= 0))) {
    stop("Cut points must be strictly ordered inside their prior bounds.")
  }
  z_tau <- numeric(K)
  previous <- -tau_bound
  for (r in seq_len(K)) {
    fraction <- (tau_row[r] - previous) / (tau_bound - previous)
    log_remainder <- if (fraction < 0.5) log1p(-fraction) else {
      log(tau_bound - tau_row[r]) - log(tau_bound - previous)
    }
    z_tau[r] <- qnorm((K - r + 1L) * log_remainder,
                      lower.tail = FALSE, log.p = TRUE)
    previous <- tau_row[r]
  }
  z <- c(z_A, z_tau)
  if (any(!is.finite(z))) stop("Current measurement row has no finite Gaussian coordinates.")
  if (is.null(measurement_row_from_gaussian(z, j, length(A_row), K, ident,
                                           s_A, tau_bound))) {
    stop("Current measurement row cannot be represented without boundary rounding.")
  }
  z
}


measurement_row_from_gaussian <- function(z, j, d, K, ident,
                                           s_A = 2.5, tau_bound = 6) {
  active <- measurement_row_active(j, d, ident)
  k_A <- length(active)
  if (length(z) != k_A + K || any(!is.finite(z))) return(NULL)
  A_row <- numeric(d)
  A_row[active] <- s_A * z[seq_len(k_A)]
  log_jacobian <- k_A * log(s_A)
  if (ident == "lower_triangular" && j <= d) {
    lower_tail <- z[j] < 0
    A_row[j] <- s_A * sqrt(qchisq(pnorm(z[j], lower.tail = lower_tail, log.p = TRUE),
      df = 1, lower.tail = lower_tail, log.p = TRUE))
    if (!is.finite(A_row[j]) || A_row[j] <= 0) return(NULL)
    log_jacobian <- log_jacobian + dnorm(z[j], log = TRUE) -
      dnorm(A_row[j] / s_A, log = TRUE) - log(2)
  }
  tau_row <- numeric(K)
  previous <- -tau_bound
  for (r in seq_len(K)) {
    fraction <- -expm1(pnorm(z[k_A + r], lower.tail = FALSE, log.p = TRUE) /
                       (K - r + 1L))
    tau_row[r] <- previous + (tau_bound - previous) * fraction
    previous <- tau_row[r]
  }
  if (any(!is.finite(A_row)) || (K &&
      any(diff(c(-tau_bound, tau_row, tau_bound)) <= 0))) return(NULL)
  if (K) log_jacobian <- log_jacobian +
    sum(dnorm(z[k_A + seq_len(K)], log = TRUE)) -
    lgamma(K + 1L) + K * log(2 * tau_bound)
  list(A = A_row, tau = tau_row, log_jacobian = log_jacobian)
}


ordinal_row_loglik_marginal <- function(C_row, U, A_row, tau_row) {
  mu <- as.numeric(U %*% A_row)
  lower <- c(-Inf, tau_row)[C_row] - mu
  upper <- c(tau_row, Inf)[C_row] - mu
  lp <- log_normal_interval_prob(lower, upper)
  if (any(!is.finite(lp))) -Inf else sum(lp)
}


update_measurement_joint_ess <- function(C, U, A, tau, m_vec,
                                         s_A = 2.5, tau_bound = 6,
                                         ident = c("lower_triangular", "none"),
                                         max_try = 5000L) {
  ident <- match.arg(ident)
  n_eval <- 0L
  angles <- numeric(nrow(A))
  for (j in seq_len(nrow(A))) {
    K <- m_vec[j] - 1L
    z <- measurement_row_to_gaussian(A[j, ], tau[[j]], j, ident, s_A, tau_bound)
    current <- ordinal_row_loglik_marginal(C[, j], U, A[j, ], tau[[j]])
    n_eval <- n_eval + 1L
    if (!is.finite(current)) stop("Current joint measurement density is non-finite.")
    nu <- rnorm(length(z))
    log_height <- current + log(runif(1))
    angle <- runif(1, 0, 2 * base::pi)
    lower <- angle - 2 * base::pi
    upper <- angle
    accepted <- FALSE
    for (attempt in seq_len(max_try)) {
      candidate <- measurement_row_from_gaussian(z * cos(angle) + nu * sin(angle),
        j, ncol(A), K, ident, s_A, tau_bound)
      proposed <- if (is.null(candidate)) -Inf else {
        ordinal_row_loglik_marginal(C[, j], U, candidate$A, candidate$tau)
      }
      n_eval <- n_eval + 1L
      if (is.finite(proposed) && proposed >= log_height) {
        A[j, ] <- candidate$A
        tau[[j]] <- candidate$tau
        angles[j] <- angle
        accepted <- TRUE
        break
      }
      if (angle < 0) lower <- angle else upper <- angle
      angle <- runif(1, lower, upper)
    }
    if (!accepted) stop("Joint measurement ESS exceeded max_try.")
  }
  list(A = A, tau = tau, n_eval = n_eval, n_rows = nrow(A), angles = angles)
}


## Score-conditioned Gaussian noncentering of the missing latent rows. With
## P=I+A'A and R=chol(P), U=M+W R^{-T}, M=S A P^{-1}; W has iid N(0,1)
## entries in the measurement model. Calibration is never transported.
loading_transport_reference <- function(S, A) {
  precision_chol <- safe_chol(diag(ncol(A)) + crossprod(A))
  mean <- t(solve_chol(precision_chol, t(S %*% A)))
  list(mean = mean, precision_chol = precision_chol)
}


update_loading_transport_ess <- function(y, X, U_curr, S, A_curr,
                                          logtheta, miss_idx, s_A = 2.5,
                                          ident = c("lower_triangular", "none"),
                                          max_try = 5000L, kernel = "se",
                                          matern_nu = 2.5, gp_evaluator = NULL) {
  ident <- match.arg(ident)
  miss_idx <- as.integer(miss_idx)
  if (anyNA(miss_idx) || anyDuplicated(miss_idx) ||
      any(miss_idx < 1L | miss_idx > nrow(U_curr))) {
    stop("miss_idx must index distinct missing rows.")
  }
  if (!is.null(gp_evaluator)) gp_evaluator$validate_block(miss_idx)
  q <- nrow(A_curr)
  d <- ncol(A_curr)
  row_length <- vapply(seq_len(q), function(j)
    length(measurement_row_active(j, d, ident)), integer(1))
  row_end <- cumsum(row_length)
  row_start <- row_end - row_length + 1L
  z <- unlist(lapply(seq_len(q), function(j)
    measurement_row_to_gaussian(A_curr[j, ], numeric(0), j, ident, s_A)),
    use.names = FALSE)
  ref <- loading_transport_reference(S, A_curr)
  W <- (U_curr[miss_idx, , drop = FALSE] - ref$mean[miss_idx, , drop = FALSE]) %*%
    t(ref$precision_chol)
  candidate_gp <- NULL
  evaluate <- function(A, U) {
    gp <- if (is.null(gp_evaluator)) {
      gp_loglik_integrated_general(y, X, U, logtheta, kernel, matern_nu)
    } else {
      candidate_gp <<- gp_evaluator$evaluate(U, logtheta)
      gp_loglik_from_moments_general(candidate_gp)
    }
    gp + loading_transport_score_loglik(S, U, A, miss_idx)
  }
  current <- evaluate(A_curr, U_curr)
  if (!is.finite(current)) stop("Current loading-transport density is non-finite.")
  log_height <- current + log(runif(1))
  nu <- rnorm(length(z))
  angle <- runif(1, 0, 2 * base::pi)
  lower <- angle - 2 * base::pi
  upper <- angle
  n_eval <- 1L
  for (attempt in seq_len(max_try)) {
    z_prop <- z * cos(angle) + nu * sin(angle)
    A_prop <- A_curr
    valid <- TRUE
    for (j in seq_len(q)) {
      row <- measurement_row_from_gaussian(z_prop[row_start[j]:row_end[j]],
        j, d, 0L, ident, s_A)
      if (is.null(row)) {
        valid <- FALSE
        break
      }
      A_prop[j, ] <- row$A
    }
    U_prop <- U_curr
    proposed <- -Inf
    if (valid) {
      ref_prop <- loading_transport_reference(S, A_prop)
      if (length(miss_idx)) {
        U_prop[miss_idx, ] <- ref_prop$mean[miss_idx, , drop = FALSE] +
          t(backsolve(ref_prop$precision_chol, t(W)))
      }
      proposed <- evaluate(A_prop, U_prop)
    }
    n_eval <- n_eval + 1L
    if (is.finite(proposed) && proposed >= log_height) {
      if (!is.null(gp_evaluator)) gp_evaluator$accept(candidate_gp)
      return(list(A = A_prop, U = U_prop, n_eval = n_eval, angle = angle))
    }
    if (angle < 0) lower <- angle else upper <- angle
    angle <- runif(1, lower, upper)
  }
  stop("Loading transport ESS exceeded max_try.")
}


update_measurement_marginal_slice <- function(C,
                                              U,
                                              A,
                                              tau,
                                              m_vec,
                                              s_A = 2.5,
                                              tau_bound = 6,
                                              ident = c("lower_triangular", "none"),
                                              A_width = 0.5,
                                              tau_width = 0.5) {
  ident <- match.arg(ident)
  q <- nrow(A)
  d <- ncol(A)
  n_eval <- 0L

  active_A <- lapply(seq_len(q), function(j) {
    if (ident == "lower_triangular" && j <= d) seq_len(j) else seq_len(d)
  })

  for (j in seq_len(q)) {
    for (k in active_A[[j]]) {
      lower <- if (ident == "lower_triangular" && j == k && j <= d) 0 else -Inf
      logf <- function(value) {
        A_prop <- A
        A_prop[j, k] <- value
        ordinal_loglik_marginal(C, U, A_prop, tau) -
          0.5 * sum(A_prop^2) / s_A^2
      }
      ans <- bounded_slice_update(
        x0 = A[j, k],
        logf = logf,
        w = A_width,
        lower = lower,
        upper = Inf,
        max_steps_out = 50L,
        max_iter = 500L,
        fail_on_limit = TRUE
      )
      A[j, k] <- ans$x
      n_eval <- n_eval + ans$n_eval
    }
  }

  for (j in seq_len(q)) {
    if (m_vec[j] <= 1L) next
    for (r in seq_len(m_vec[j] - 1L)) {
      lower <- max(-tau_bound, if (r > 1L) tau[[j]][r - 1L] else -Inf)
      upper <- min(tau_bound, if (r < m_vec[j] - 1L) tau[[j]][r + 1L] else Inf)
      logf <- function(value) {
        tau_prop <- tau
        tau_prop[[j]][r] <- value
        ordinal_loglik_marginal(C, U, A, tau_prop)
      }
      ans <- bounded_slice_update(
        x0 = tau[[j]][r],
        logf = logf,
        w = tau_width,
        lower = lower,
        upper = upper,
        max_steps_out = 50L,
        max_iter = 500L,
        fail_on_limit = TRUE
      )
      tau[[j]][r] <- ans$x
      n_eval <- n_eval + ans$n_eval
    }
  }

  list(A = A, tau = tau, n_eval = n_eval)
}


############################################################
## MCMC diagnostics
############################################################

split_rhat <- function(chain_list) {
  chain_list <- lapply(chain_list, function(x) as.numeric(x[is.finite(x)]))
  lens <- vapply(chain_list, length, integer(1))
  n0 <- min(lens)
  
  if (length(chain_list) < 2 || n0 < 20) return(NA_real_)
  
  chain_list <- lapply(chain_list, function(x) tail(x, n0))
  
  n_half <- floor(n0 / 2)
  if (n_half < 10) return(NA_real_)
  
  split_mat <- do.call(
    rbind,
    lapply(chain_list, function(x) {
      rbind(
        x[seq_len(n_half)],
        x[(n0 - n_half + 1):n0]
      )
    })
  )
  
  n_split <- ncol(split_mat)
  
  chain_means <- rowMeans(split_mat)
  chain_vars <- apply(split_mat, 1, var)
  
  W <- mean(chain_vars)
  B <- n_split * var(chain_means)
  
  if (!is.finite(W) || W <= 0) return(NA_real_)
  
  var_hat <- ((n_split - 1) / n_split) * W + B / n_split
  
  sqrt(var_hat / W)
}


aligned_chain_matrix <- function(chain_list) {
  chain_list <- lapply(chain_list, as.numeric)
  lens <- vapply(chain_list, length, integer(1))
  n0 <- min(lens)
  if (length(chain_list) < 2L || n0 < 4L) return(NULL)
  do.call(cbind, lapply(chain_list, function(x) tail(x, n0)))
}


rank_rhat <- function(chain_list) {
  mat <- aligned_chain_matrix(chain_list)
  if (is.null(mat)) return(NA_real_)
  if (requireNamespace("posterior", quietly = TRUE)) {
    return(as.numeric(posterior::rhat(mat)))
  }
  split_rhat(chain_list)
}


bulk_ess <- function(chain_list) {
  mat <- aligned_chain_matrix(chain_list)
  if (is.null(mat)) return(NA_real_)
  if (requireNamespace("posterior", quietly = TRUE)) {
    return(as.numeric(posterior::ess_bulk(mat)))
  }
  sum(vapply(chain_list, ess_ips, numeric(1)))
}


tail_ess <- function(chain_list) {
  mat <- aligned_chain_matrix(chain_list)
  if (is.null(mat)) return(NA_real_)
  if (requireNamespace("posterior", quietly = TRUE)) {
    return(as.numeric(posterior::ess_tail(mat)))
  }
  NA_real_
}


ess_ips <- function(x, max_lag = NULL) {
  x <- as.numeric(x)
  x <- x[is.finite(x)]
  
  n <- length(x)
  if (n < 5) return(NA_real_)
  
  sx <- sd(x)
  if (!is.finite(sx)) return(NA_real_)
  if (sx == 0) return(n)
  
  if (is.null(max_lag)) {
    max_lag <- min(n - 1L, 1000L)
  } else {
    max_lag <- min(n - 1L, as.integer(max_lag))
  }
  
  if (max_lag < 1L) return(n)
  
  ac <- tryCatch(
    as.numeric(stats::acf(
      x,
      lag.max = max_lag,
      plot = FALSE,
      demean = TRUE
    )$acf),
    error = function(e) NA_real_
  )
  
  if (length(ac) <= 1 || all(!is.finite(ac))) return(NA_real_)
  
  ac <- ac[-1]
  ac <- ac[is.finite(ac)]
  
  if (length(ac) == 0) return(n)
  
  mm <- floor(length(ac) / 2)
  
  if (mm >= 1) {
    pair_sums <- ac[2 * seq_len(mm) - 1] + ac[2 * seq_len(mm)]
    first_nonpos <- which(pair_sums <= 0)[1]
    
    if (is.na(first_nonpos)) {
      use_lag <- 2 * mm
    } else {
      use_lag <- 2 * (first_nonpos - 1)
    }
  } else {
    use_lag <- ifelse(ac[1] > 0, 1, 0)
  }
  
  if (use_lag > 0) {
    tau_int <- 1 + 2 * sum(ac[seq_len(use_lag)])
  } else {
    tau_int <- 1
  }
  
  tau_int <- max(tau_int, 1)
  
  ess <- n / tau_int
  min(max(ess, 1), n)
}


############################################################
## Fully Bayesian ordinal-probit EIV-GP sampler
############################################################

fit_eivgp_ordprobit_fb <- function(...) {
  if (!exists("fit_eivgp_ordprobit_v030", mode = "function")) {
    stop("Install and load the current eivGP package before fitting.")
  }
  fit_eivgp_ordprobit_v030(...)
}


sample_u_given_c_ordprobit <- function(C_new,
                                       A,
                                       tau,
                                       m_vec,
                                       n_gibbs = 100L,
                                       U_init = NULL) {
  C_new <- as.matrix(C_new)
  A <- as.matrix(A)
  
  N <- nrow(C_new)
  q <- ncol(C_new)
  d <- ncol(A)
  
  if (nrow(A) != q || length(m_vec) != q || length(tau) != q) {
    stop("m_vec has wrong length.")
  }
  if (any(C_new < 1L) ||
      any(C_new > matrix(rep(m_vec, each = N), N, q))) {
    stop("C_new contains invalid ordinal codes.")
  }
  n_gibbs <- as.integer(n_gibbs)
  if (length(n_gibbs) != 1L || is.na(n_gibbs) || n_gibbs < 1L) {
    stop("n_gibbs must be a positive integer.")
  }
  
  if (is.null(U_init)) {
    U <- matrix(rnorm(N * d), N, d)
  } else {
    U <- as_numeric_matrix_strict(
      U_init,
      "U_init",
      nrow_expected = N,
      ncol_expected = d
    )
  }
  
  S <- matrix(0, N, q)
  V <- solve(diag(d) + crossprod(A))
  
  for (gg in seq_len(n_gibbs)) {
    for (j in seq_len(q)) {
      mu <- as.numeric(U %*% A[j, ])
      
      tau_j <- tau[[j]]
      
      lower <- c(-Inf, tau_j)[C_new[, j]]
      upper <- c(tau_j, Inf)[C_new[, j]]
      
      S[, j] <- rtruncnorm_vec(mu, 1, lower, upper)
    }
    
    M <- S %*% A %*% V
    U <- rmvnorm_rows_common(M, V)
  }
  
  U
}


sample_u_given_c_ordprobit_rejection <- function(
    C_new,
    A,
    tau,
    m_vec,
    batch_size = NULL,
    max_batches = 1000L) {
  C_new <- as.matrix(C_new)
  storage.mode(C_new) <- "integer"
  A <- as.matrix(A)

  N <- nrow(C_new)
  q <- ncol(C_new)
  d <- ncol(A)
  if (N < 1L || nrow(A) != q || length(m_vec) != q || length(tau) != q) {
    stop("Incompatible inputs in sample_u_given_c_ordprobit_rejection().")
  }
  if (anyNA(C_new) || any(C_new < 1L) ||
      any(C_new > matrix(rep(m_vec, each = N), N, q))) {
    stop("C_new contains invalid ordinal codes.")
  }
  if (is.null(batch_size)) batch_size <- max(5000L, 20L * N)
  batch_size <- as.integer(batch_size)
  max_batches <- as.integer(max_batches)
  if (is.na(batch_size) || batch_size < 1L ||
      is.na(max_batches) || max_batches < 1L) {
    stop("batch_size and max_batches must be positive integers.")
  }

  target_keys <- pattern_key(C_new)
  unique_keys <- unique(target_keys)
  required <- table(factor(target_keys, levels = unique_keys))
  accepted <- setNames(vector("list", length(unique_keys)), unique_keys)
  accepted_n <- setNames(integer(length(unique_keys)), unique_keys)
  proposal_hits <- setNames(numeric(length(unique_keys)), unique_keys)
  total_candidates <- 0
  batches_used <- 0L

  for (batch_id in seq_len(max_batches)) {
    batches_used <- batch_id
    U_batch <- matrix(rnorm(batch_size * d), batch_size, d)
    S_batch <- U_batch %*% t(A) + matrix(rnorm(batch_size * q), batch_size, q)
    C_batch <- matrix(NA_integer_, batch_size, q)
    for (j in seq_len(q)) {
      C_batch[, j] <- findInterval(
        S_batch[, j], c(-Inf, tau[[j]], Inf), rightmost.closed = TRUE
      )
    }
    candidate_keys <- pattern_key(C_batch)
    total_candidates <- total_candidates + batch_size

    for (key in unique_keys) {
      hit <- which(candidate_keys == key)
      proposal_hits[[key]] <- proposal_hits[[key]] + length(hit)
      remaining <- required[[key]] - accepted_n[[key]]
      if (remaining <= 0L) next
      if (length(hit) == 0L) next
      take <- hit[seq_len(min(remaining, length(hit)))]
      accepted[[key]][[length(accepted[[key]]) + 1L]] <-
        U_batch[take, , drop = FALSE]
      accepted_n[[key]] <- accepted_n[[key]] + length(take)
    }

    if (all(accepted_n >= required)) break
  }

  telemetry <- data.frame(
    pattern = unique_keys,
    requested_draws = as.integer(required),
    retained_draws = as.integer(accepted_n),
    proposal_hits = as.numeric(proposal_hits),
    shortfall = pmax(as.integer(required) - as.integer(accepted_n), 0L),
    complete = as.integer(accepted_n) >= as.integer(required),
    total_candidates = rep.int(total_candidates, length(unique_keys)),
    batches_used = rep.int(batches_used, length(unique_keys)),
    empirical_pattern_probability =
      as.numeric(proposal_hits) / total_candidates,
    retained_fraction_of_candidates =
      as.numeric(accepted_n) / total_candidates,
    stringsAsFactors = FALSE
  )

  if (any(accepted_n < required)) {
    missing <- unique_keys[accepted_n < required]
    message <- paste0(
      "Exact rejection sampling did not fill ordinal pattern(s): ",
      paste(missing, collapse = ", "), ". Use the exact default ",
      "latent_sampler = 'minimax_tilting', or increase ",
      "rejection_max_batches. Finite Gibbs is diagnostic only."
    )
    stop(structure(
      list(
        message = message,
        call = NULL,
        pattern_telemetry = telemetry,
        candidate_draws = total_candidates,
        retained_draws = sum(accepted_n),
        overall_acceptance = sum(accepted_n) / total_candidates
      ),
      class = c("mixedgp_rejection_sampling_error", "error", "condition")
    ))
  }

  out <- matrix(NA_real_, N, d)
  for (key in unique_keys) {
    rows <- which(target_keys == key)
    draws <- do.call(rbind, accepted[[key]])
    out[rows, ] <- draws[seq_along(rows), , drop = FALSE]
  }
  attr(out, "latent_sampler") <- "exact_rejection"
  attr(out, "candidate_draws") <- total_candidates
  attr(out, "overall_acceptance") <- N / total_candidates
  attr(out, "pattern_telemetry") <- telemetry
  attr(out, "rejection_telemetry") <- telemetry
  out
}


checked_minimax_tilting_draw <- function(draw_function) {
  low_acceptance_warnings <- character(0)
  draws <- withCallingHandlers(draw_function(), warning = function(w) {
    message <- conditionMessage(w)
    if (grepl("^Acceptance prob(ability|\\.) smaller than", message)) {
      ## A completed accept--reject run remains exact when acceptance is low.
      ## Record this performance warning; never suppress solver/degeneracy
      ## warnings, which may invalidate the envelope or change the target.
      low_acceptance_warnings <<- c(low_acceptance_warnings, message)
      invokeRestart("muffleWarning")
    }
    stop("Minimax tilting cannot be certified exact after this warning: ",
         message, call. = FALSE)
  })
  list(draws = draws, low_acceptance_warnings = low_acceptance_warnings)
}


sample_u_given_c_ordprobit_minimax <- function(C_new,
                                                A,
                                                tau,
                                                m_vec) {
  if (!requireNamespace("TruncatedNormal", quietly = TRUE)) {
    stop(
      "latent_sampler='minimax_tilting' requires the CRAN package ",
      "TruncatedNormal.",
      call. = FALSE
    )
  }
  C_new <- as.matrix(C_new)
  storage.mode(C_new) <- "integer"
  A <- as.matrix(A)
  N <- nrow(C_new)
  q <- ncol(C_new)
  d <- ncol(A)
  m_vec <- as.integer(m_vec)
  if (N < 1L || nrow(A) != q || length(m_vec) != q || length(tau) != q ||
      anyNA(C_new) || any(!is.finite(A)) ||
      any(C_new < 1L) ||
      any(C_new > matrix(rep(m_vec, each = N), N, q))) {
    stop("Incompatible inputs in sample_u_given_c_ordprobit_minimax().")
  }
  for (j in seq_len(q)) {
    if (length(tau[[j]]) != m_vec[j] - 1L ||
        any(!is.finite(tau[[j]])) ||
        is.unsorted(tau[[j]], strictly = TRUE)) {
      stop("tau contains an invalid cutpoint vector.")
    }
  }

  ## If S = A U + epsilon with U ~ N(0,I_d) and epsilon ~ N(0,I_q),
  ## then S ~ N(0, A A' + I_q). Use the minimax-tilted accept--reject sampler
  ## for S in its ordinal rectangle; on successful completion this draw and the
  ## following Gaussian U | S update are exact.
  Sigma_S <- tcrossprod(A) + diag(q)
  precision_U_given_S <- diag(d) + crossprod(A)
  chol_precision <- safe_chol(precision_U_given_S)
  V <- solve_chol(chol_precision, diag(d))
  B <- V %*% t(A)
  V <- 0.5 * (V + t(V))
  chol_V <- safe_chol(V)

  target_keys <- pattern_key(C_new)
  unique_keys <- unique(target_keys)
  out <- matrix(NA_real_, N, d)
  telemetry <- vector("list", length(unique_keys))
  for (ii in seq_along(unique_keys)) {
    key <- unique_keys[ii]
    rows <- which(target_keys == key)
    c_key <- C_new[rows[1L], ]
    lower <- upper <- numeric(q)
    for (j in seq_len(q)) {
      lower[j] <- c(-Inf, tau[[j]])[c_key[j]]
      upper[j] <- c(tau[[j]], Inf)[c_key[j]]
    }
    checked <- checked_minimax_tilting_draw(function() {
      TruncatedNormal::rtmvnorm(n = length(rows), mu = rep(0, q),
        sigma = Sigma_S, lb = lower, ub = upper)
    })
    S_draw <- checked$draws
    if (length(S_draw) != length(rows) * q || any(!is.finite(S_draw))) {
      stop("Minimax tilting returned malformed or non-finite score draws.")
    }
    S_draw <- matrix(S_draw, nrow = length(rows), ncol = q)
    lower_mat <- matrix(lower, nrow = length(rows), ncol = q, byrow = TRUE)
    upper_mat <- matrix(upper, nrow = length(rows), ncol = q, byrow = TRUE)
    if (any(S_draw < lower_mat | S_draw > upper_mat)) {
      stop("Minimax tilting returned a score outside its ordinal rectangle.")
    }
    U_mean <- S_draw %*% t(B)
    U_draw <- U_mean + matrix(rnorm(length(rows) * d), length(rows), d) %*%
      chol_V
    if (any(!is.finite(U_draw))) {
      stop("The exact Gaussian U | S update returned a non-finite draw.")
    }
    out[rows, ] <- U_draw
    telemetry[[ii]] <- data.frame(
      pattern = key,
      n_draw = length(rows),
      exactness_status = "completed_without_target_validity_warnings",
      low_acceptance_warning_count = length(checked$low_acceptance_warnings),
      low_acceptance_warnings = paste(unique(checked$low_acceptance_warnings),
                                      collapse = "; "),
      stringsAsFactors = FALSE
    )
  }
  attr(out, "latent_sampler") <- "exact_minimax_tilting"
  attr(out, "pattern_telemetry") <- do.call(rbind, telemetry)
  out
}


sample_u_given_c_ordprobit_dispatch <- function(
    C_new,
    A,
    tau,
    m_vec,
    latent_sampler = c("minimax_tilting", "rejection", "gibbs"),
    n_gibbs = 100L,
    rejection_batch_size = NULL,
    rejection_max_batches = 1000L) {
  latent_sampler <- match.arg(latent_sampler)
  if (latent_sampler == "minimax_tilting") {
    return(sample_u_given_c_ordprobit_minimax(C_new, A, tau, m_vec))
  }
  if (latent_sampler == "rejection") {
    return(sample_u_given_c_ordprobit_rejection(
      C_new = C_new,
      A = A,
      tau = tau,
      m_vec = m_vec,
      batch_size = rejection_batch_size,
      max_batches = rejection_max_batches
    ))
  }
  out <- sample_u_given_c_ordprobit(
    C_new = C_new,
    A = A,
    tau = tau,
    m_vec = m_vec,
    n_gibbs = n_gibbs
  )
  attr(out, "latent_sampler") <- "finite_gibbs_diagnostic"
  attr(out, "gibbs_sweeps") <- as.integer(n_gibbs)
  out
}


sample_eiv_test_y_ordprobit_fb <- function(X_test_raw,
                                           C_test,
                                           fit_obj,
                                           draw_ids = NULL,
                                           n_per_draw = 1L,
                                           n_new_latent_gibbs = 100L,
                                           latent_sampler = c(
                                             "minimax_tilting", "rejection", "gibbs"
                                           ),
                                           rejection_batch_size = NULL,
                                           rejection_max_batches = 1000L,
                                           U_test_obs = NULL,
                                           U_input_scale = c("raw", "model"),
                                           predictive_target = c("response", "latent"),
                                           joint = FALSE,
                                           seed = NULL) {
  if (!is.null(mixedgp_budget()) && !isTRUE(getOption("eivGP.prediction_block",FALSE))) {
    args <- lapply(as.list(match.call(expand.dots=TRUE))[-1L],eval,envir=parent.frame())
    return(mixedgp_predict_blocks(args,FUN=sample_eiv_test_y_ordprobit_fb,object_name="fit_obj"))
  }

  predictive_target <- match.arg(predictive_target)
  latent_sampler <- match.arg(latent_sampler)
  U_input_scale <- match.arg(U_input_scale)
  if (!is.null(seed)) set.seed(seed)
  
  samples_U <- fit_obj$mcmc$samples_U
  samples_A <- fit_obj$mcmc$samples_A
  samples_tau <- fit_obj$mcmc$samples_tau
  samples_logtheta <- fit_obj$mcmc$samples_logtheta
  samples_sigma2 <- fit_obj$mcmc$samples_sigma2
  
  X_train <- fit_obj$data$X
  y_train <- fit_obj$data$y
  
  y_center <- fit_obj$data$y_center
  y_scale <- fit_obj$data$y_scale
  
  X_center <- fit_obj$data$X_center
  X_scale <- fit_obj$data$X_scale
  
  m_vec <- fit_obj$data$m_vec
  d <- fit_obj$data$d
  q <- fit_obj$data$q
  p <- fit_obj$data$p
  n_train <- nrow(X_train)

  fit_level_maps <- fit_obj$data$C_level_maps
  if (is.null(fit_level_maps)) {
    fit_level_maps <- lapply(m_vec, seq_len)
  }
  fit_C_names <- fit_obj$data$C_names
  if (is.null(fit_C_names)) fit_C_names <- colnames(fit_obj$data$C_ord)
  fit_U_center <- fit_obj$data$U_center
  if (is.null(fit_U_center)) fit_U_center <- rep(0, d)
  fit_U_scale <- fit_obj$data$U_scale
  if (is.null(fit_U_scale)) fit_U_scale <- rep(1, d)
  kernel_spec <- kernel_spec_from_fit(fit_obj)

  ordinal_test <- prepare_ordinal_matrix(
    C_test,
    m_vec = m_vec,
    level_maps = fit_level_maps,
    expected_names = fit_C_names,
    name = "C_test"
  )
  C_test <- ordinal_test$C
  n_test <- nrow(C_test)

  X_test_names <- if (is.null(dim(X_test_raw))) names(X_test_raw) else colnames(X_test_raw)
  X_test_raw <- as_numeric_matrix_strict(
    X_test_raw,
    "X_test_raw",
    nrow_expected = n_test,
    ncol_expected = p
  )
  train_x_names <- colnames(fit_obj$data$X_raw)
  if (!is.null(X_test_names) && all(nzchar(X_test_names)) &&
      !is.null(train_x_names)) {
    if (!setequal(X_test_names, train_x_names)) {
      stop("X_test_raw column names do not match the training X columns.")
    }
    X_test_raw <- X_test_raw[, train_x_names, drop = FALSE]
  }

  U_test_complete <- integer(0)
  U_test_standardized <- matrix(NA_real_, n_test, d)
  if (!is.null(U_test_obs)) {
    U_test_names <- if (is.null(dim(U_test_obs))) names(U_test_obs) else colnames(U_test_obs)
    U_test_obs <- as_numeric_matrix_strict(
      U_test_obs,
      "U_test_obs",
      nrow_expected = n_test,
      ncol_expected = d,
      allow_na = TRUE
    )
    if (!is.null(U_test_names) && all(nzchar(U_test_names)) &&
        !is.null(fit_obj$data$U_names)) {
      if (!setequal(U_test_names, fit_obj$data$U_names)) {
        stop("U_test_obs column names do not match the training U columns.")
      }
      U_test_obs <- U_test_obs[, fit_obj$data$U_names, drop = FALSE]
    }
    observed_per_row <- rowSums(!is.na(U_test_obs))
    if (any(observed_per_row > 0L & observed_per_row < d)) {
      stop(
        "Partially observed rows of U_test_obs are not supported; provide all ",
        "d coordinates for a row or leave the entire row missing."
      )
    }
    U_test_complete <- which(observed_per_row == d)
    if (length(U_test_complete) > 0L) {
      U_test_standardized[U_test_complete, ] <- if (U_input_scale == "raw") {
        sweep(
          sweep(
            U_test_obs[U_test_complete, , drop = FALSE],
            2,
            fit_U_center,
            "-"
          ),
          2,
          fit_U_scale,
          "/"
        )
      } else {
        U_test_obs[U_test_complete, , drop = FALSE]
      }
    }
  }

  if (is.null(draw_ids)) draw_ids <- seq_len(nrow(samples_logtheta))
  draw_ids <- as.integer(draw_ids)
  if (length(draw_ids) < 1L || anyNA(draw_ids) ||
      any(draw_ids < 1L | draw_ids > nrow(samples_logtheta))) {
    stop("draw_ids contains invalid posterior-draw indices.")
  }
  n_per_draw <- as.integer(n_per_draw)
  if (length(n_per_draw) != 1L || is.na(n_per_draw) || n_per_draw < 1L) {
    stop("n_per_draw must be a positive integer.")
  }
  
  X_test <- sweep(sweep(X_test_raw, 2, X_center, "-"), 2, X_scale, "/")
  n_draw <- length(draw_ids) * n_per_draw
  
  out <- matrix(NA_real_, nrow = n_draw, ncol = n_test)
  conditional_means <- matrix(NA_real_, nrow = n_draw, ncol = n_test)
  conditional_vars <- matrix(NA_real_, nrow = n_draw, ncol = n_test)
  latent_acceptance <- rep(NA_real_, n_draw)
  latent_candidates <- integer(n_draw)
  
  row_id <- 0L
  
  for (s in draw_ids) {
    A_s <- matrix(samples_A[s, , ], nrow = q, ncol = d)
    tau_s <- unflatten_tau(samples_tau[s, ], m_vec)
    
    U_train_s <- matrix(samples_U[s, , ], nrow = n_train, ncol = d)
    
    for (rr in seq_len(n_per_draw)) {
      row_id <- row_id + 1L
      
      missing_test <- setdiff(seq_len(n_test), U_test_complete)
      U_star <- U_test_standardized
      if (length(missing_test) > 0L) {
        U_missing <- sample_u_given_c_ordprobit_dispatch(
          C_new = C_test[missing_test, , drop = FALSE],
          A = A_s,
          tau = tau_s,
          m_vec = m_vec,
          latent_sampler = latent_sampler,
          n_gibbs = n_new_latent_gibbs,
          rejection_batch_size = rejection_batch_size,
          rejection_max_batches = rejection_max_batches
        )
        U_star[missing_test, ] <- U_missing
        if (latent_sampler == "rejection") {
          latent_acceptance[row_id] <- attr(U_missing, "overall_acceptance")
          latent_candidates[row_id] <- attr(U_missing, "candidate_draws")
        }
      }
      if (length(U_test_complete) > 0L) {
        U_star[U_test_complete, ] <-
          U_test_standardized[U_test_complete, , drop = FALSE]
      }
      
      pred <- gp_predict_draw_general(
        X_train = X_train,
        U_train = U_train_s,
        y_train = y_train,
        X_star = X_test,
        U_star = U_star,
        logtheta = samples_logtheta[s, ],
        sigma2_eps = samples_sigma2[s],
        noisy = predictive_target == "response",
        kernel = kernel_spec$name,
        matern_nu = kernel_spec$matern_nu,
        return_cov = isTRUE(joint)
      )
      
      if (isTRUE(joint)) {
        y_std <- as.numeric(rmvnorm_psd(1L, pred$mean, pred$cov))
      } else {
        y_std <- pred$mean + sqrt(pred$var) * rnorm(n_test)
      }
      out[row_id, ] <- y_center + y_scale * y_std
      conditional_means[row_id, ] <- y_center + y_scale * pred$mean
      conditional_vars[row_id, ] <- y_scale^2 * pred$var
    }
  }
  
  attr(out, "predictive_target") <- predictive_target
  attr(out, "joint") <- isTRUE(joint)
  attr(out, "latent_sampler") <- latent_sampler
  attr(out, "latent_input_scale") <- U_input_scale
  attr(out, "rejection_acceptance") <- latent_acceptance
  attr(out, "rejection_candidate_draws") <- latent_candidates
  attr(out, "conditional_means") <- conditional_means
  attr(out, "conditional_vars") <- conditional_vars
  attr(out, "mixture_components") <-
    "posterior and latent-input Monte Carlo Gaussian components"
  out
}


posterior_u_draws <- function(fit_obj,
                              rows = NULL,
                              scale = c("auto", "raw", "model")) {
  scale <- match.arg(scale)
  anchored <- isTRUE(fit_obj$data$latent_scale_anchored)
  if (scale == "raw" && !anchored) {
    stop("scale='raw' requires calibration data that anchor latent U.")
  }
  if (scale == "auto") scale <- if (anchored) "raw" else "model"
  draws <- fit_obj$mcmc$samples_U
  n <- dim(draws)[2]
  d <- dim(draws)[3]

  if (is.null(rows)) rows <- seq_len(n)
  rows <- as.integer(rows)
  if (length(rows) < 1L || anyNA(rows) || any(rows < 1L | rows > n)) {
    stop("rows contains invalid training-row indices.")
  }
  draws <- draws[, rows, , drop = FALSE]

  if (scale == "raw") {
    U_center <- fit_obj$data$U_center
    if (is.null(U_center)) U_center <- rep(0, d)
    U_scale <- fit_obj$data$U_scale
    if (is.null(U_scale)) U_scale <- rep(1, d)
    for (k in seq_len(d)) {
      draws[, , k] <- U_center[k] + U_scale[k] * draws[, , k]
    }
  }

  draw_dimnames <- dimnames(draws)
  if (is.null(draw_dimnames)) draw_dimnames <- vector("list", 3L)
  draw_dimnames[[2]] <- as.character(rows)
  draw_dimnames[[3]] <- fit_obj$data$U_names
  dimnames(draws) <- draw_dimnames
  attr(draws, "scale") <- if (!anchored && scale == "model") "working" else scale
  attr(draws, "latent_scale_anchored") <- anchored
  draws
}


sample_eiv_f_given_xu_fb <- function(X_test_raw,
                                     U_test_raw,
                                     fit_obj,
                                     draw_ids = NULL,
                                     n_per_draw = 1L,
                                     include_process_uncertainty = TRUE,
                                     U_input_scale = c("raw", "model"),
                                     joint = FALSE,
                                     seed = NULL) {
  if (!is.null(mixedgp_budget()) && !isTRUE(getOption("eivGP.prediction_block",FALSE))) {
    args <- lapply(as.list(match.call(expand.dots=TRUE))[-1L],eval,envir=parent.frame())
    return(mixedgp_predict_blocks(args,FUN=sample_eiv_f_given_xu_fb,object_name="fit_obj"))
  }

  U_input_scale <- match.arg(U_input_scale)
  if (!is.null(seed)) set.seed(seed)

  p <- fit_obj$data$p
  d <- fit_obj$data$d
  X_input_names <- if (is.null(dim(X_test_raw))) names(X_test_raw) else colnames(X_test_raw)
  U_input_names <- if (is.null(dim(U_test_raw))) names(U_test_raw) else colnames(U_test_raw)

  X_test_raw <- as_numeric_matrix_strict(
    X_test_raw,
    "X_test_raw",
    ncol_expected = p
  )
  U_test_raw <- as_numeric_matrix_strict(
    U_test_raw,
    "U_test_raw",
    nrow_expected = nrow(X_test_raw),
    ncol_expected = d
  )

  train_x_names <- colnames(fit_obj$data$X_raw)
  if (!is.null(X_input_names) && all(nzchar(X_input_names)) &&
      !is.null(train_x_names)) {
    if (!setequal(X_input_names, train_x_names)) {
      stop("X_test_raw column names do not match the training X columns.")
    }
    X_test_raw <- X_test_raw[, train_x_names, drop = FALSE]
  }
  train_u_names <- fit_obj$data$U_names
  if (!is.null(U_input_names) && all(nzchar(U_input_names)) &&
      !is.null(train_u_names)) {
    if (!setequal(U_input_names, train_u_names)) {
      stop("U_test_raw column names do not match the training U columns.")
    }
    U_test_raw <- U_test_raw[, train_u_names, drop = FALSE]
  }

  U_center <- fit_obj$data$U_center
  if (is.null(U_center)) U_center <- rep(0, d)
  U_scale <- fit_obj$data$U_scale
  if (is.null(U_scale)) U_scale <- rep(1, d)
  X_test <- sweep(
    sweep(X_test_raw, 2, fit_obj$data$X_center, "-"),
    2,
    fit_obj$data$X_scale,
    "/"
  )
  U_test <- if (U_input_scale == "raw") {
    sweep(sweep(U_test_raw, 2, U_center, "-"), 2, U_scale, "/")
  } else {
    U_test_raw
  }

  samples_logtheta <- fit_obj$mcmc$samples_logtheta
  if (is.null(draw_ids)) draw_ids <- seq_len(nrow(samples_logtheta))
  draw_ids <- as.integer(draw_ids)
  if (length(draw_ids) < 1L || anyNA(draw_ids) ||
      any(draw_ids < 1L | draw_ids > nrow(samples_logtheta))) {
    stop("draw_ids contains invalid posterior-draw indices.")
  }
  n_per_draw <- as.integer(n_per_draw)
  if (length(n_per_draw) != 1L || is.na(n_per_draw) || n_per_draw < 1L) {
    stop("n_per_draw must be a positive integer.")
  }

  kernel_spec <- kernel_spec_from_fit(fit_obj)
  n_train <- nrow(fit_obj$data$X)
  out <- matrix(NA_real_, length(draw_ids) * n_per_draw, nrow(X_test))

  row_id <- 0L
  for (ii in seq_along(draw_ids)) {
    s <- draw_ids[ii]
    U_train_s <- matrix(
      fit_obj$mcmc$samples_U[s, , ],
      nrow = n_train,
      ncol = d
    )
    pred <- gp_predict_draw_general(
      X_train = fit_obj$data$X,
      U_train = U_train_s,
      y_train = fit_obj$data$y,
      X_star = X_test,
      U_star = U_test,
      logtheta = samples_logtheta[s, ],
      sigma2_eps = fit_obj$mcmc$samples_sigma2[s],
      noisy = FALSE,
      kernel = kernel_spec$name,
      matern_nu = kernel_spec$matern_nu,
      return_cov = isTRUE(joint) && isTRUE(include_process_uncertainty)
    )

    for (rr in seq_len(n_per_draw)) {
      row_id <- row_id + 1L
      if (!isTRUE(include_process_uncertainty)) {
        f_std <- pred$mean
      } else if (isTRUE(joint)) {
        f_std <- as.numeric(rmvnorm_psd(1L, pred$mean, pred$cov))
      } else {
        f_std <- pred$mean + sqrt(pred$var) * rnorm(nrow(X_test))
      }
      out[row_id, ] <- fit_obj$data$y_center + fit_obj$data$y_scale * f_std
    }
  }

  attr(out, "joint") <- isTRUE(joint)
  attr(out, "include_process_uncertainty") <- isTRUE(include_process_uncertainty)
  attr(out, "latent_input_scale") <- U_input_scale
  out
}


gp_integrated_mean_state_general <- function(X_train,
                                             U_train,
                                             y_train,
                                             X_star,
                                             U_mc,
                                             logtheta,
                                             sigma2_eps,
                                             kernel = "se",
                                             matern_nu = 2.5,
                                             return_cov = FALSE) {
  X_train <- as.matrix(X_train)
  U_train <- as.matrix(U_train)
  X_star <- as.matrix(X_star)
  U_mc <- as.array(U_mc)
  n <- nrow(X_train)
  N <- nrow(X_star)
  p <- ncol(X_train)
  d <- ncol(U_train)
  if (length(dim(U_mc)) != 3L || dim(U_mc)[2L] != N || dim(U_mc)[3L] != d) {
    stop("U_mc must be an M by N by d array.")
  }
  M <- dim(U_mc)[1L]
  if (nrow(U_train) != n || length(y_train) != n || ncol(X_star) != p ||
      M < 1L || any(!is.finite(c(X_train, U_train, y_train, X_star, U_mc)))) {
    stop("Incompatible or nonfinite inputs in gp_integrated_mean_state_general().")
  }

  kernel_spec <- normalize_gp_kernel(kernel, matern_nu)
  rho <- exp(logtheta[1L])
  theta_x <- exp(logtheta[1L + seq_len(p)])
  theta_u <- exp(logtheta[1L + p + seq_len(d)])
  D2_train <- weighted_sqdist_general(X_train, U_train, theta_x, theta_u)
  R_train <- kernel_from_weighted_sqdist(
    D2_train, kernel_spec$name, kernel_spec$matern_nu
  )
  C_train <- sigma2_eps * (rho^2 * R_train + diag(n))
  chol_train <- safe_chol(C_train)
  alpha <- solve_chol(chol_train, y_train)

  K_bar <- matrix(NA_real_, N, n)
  prior_diag <- numeric(N)
  for (i in seq_len(N)) {
    U_i <- matrix(U_mc[, i, ], nrow = M, ncol = d)
    X_i <- X_star[rep(i, M), , drop = FALSE]
    D2_cross <- weighted_sqdist_general(
      X_i, U_i, theta_x, theta_u,
      X_b = X_train, U_b = U_train
    )
    K_cross <- rho^2 * sigma2_eps * kernel_from_weighted_sqdist(
      D2_cross, kernel_spec$name, kernel_spec$matern_nu
    )
    K_bar[i, ] <- colMeans(K_cross)

    D2_ii <- weighted_sqdist_general(X_i, U_i, theta_x, theta_u)
    prior_diag[i] <- mean(
      rho^2 * sigma2_eps * kernel_from_weighted_sqdist(
        D2_ii, kernel_spec$name, kernel_spec$matern_nu
      )
    )
  }

  conditional_mean <- as.numeric(K_bar %*% alpha)
  v <- forwardsolve(t(chol_train), t(K_bar))
  var_tolerance <- 100 * (n + N + M) * .Machine$double.eps *
    max(1, prior_diag)
  if (!isTRUE(return_cov)) {
    conditional_var <- prior_diag - colSums(v^2)
    if (any(conditional_var < -var_tolerance)) {
      stop("Integrated GP conditional variance is materially negative.")
    }
    conditional_var <- pmax(conditional_var, 0)
    return(list(mean = conditional_mean, var = conditional_var))
  }

  prior_cov <- matrix(NA_real_, N, N)
  diag(prior_cov) <- prior_diag
  if (N > 1L) {
    for (i in seq_len(N - 1L)) {
      U_i <- matrix(U_mc[, i, ], nrow = M, ncol = d)
      X_i <- X_star[rep(i, M), , drop = FALSE]
      for (j in (i + 1L):N) {
        U_j <- matrix(U_mc[, j, ], nrow = M, ncol = d)
        X_j <- X_star[rep(j, M), , drop = FALSE]
        D2_ij <- weighted_sqdist_general(
          X_i, U_i, theta_x, theta_u,
          X_b = X_j, U_b = U_j
        )
        value <- mean(
          rho^2 * sigma2_eps * kernel_from_weighted_sqdist(
            D2_ij, kernel_spec$name, kernel_spec$matern_nu
          )
        )
        prior_cov[i, j] <- value
        prior_cov[j, i] <- value
      }
    }
  }
  conditional_cov <- prior_cov - crossprod(v)
  conditional_cov <- 0.5 * (conditional_cov + t(conditional_cov))
  if (any(diag(conditional_cov) < -var_tolerance)) {
    stop("Integrated GP conditional covariance has a materially negative diagonal.")
  }
  diag(conditional_cov) <- pmax(diag(conditional_cov), 0)
  list(
    mean = conditional_mean,
    var = diag(conditional_cov),
    cov = conditional_cov
  )
}


sample_eiv_m_given_xc_fb <- function(X_test_raw,
                                      C_test,
                                      fit_obj,
                                      draw_ids = NULL,
                                      n_per_draw = 1L,
                                      n_latent = 256L,
                                      include_process_uncertainty = TRUE,
                                      joint = FALSE,
                                      return_components = FALSE,
                                      latent_sampler = c(
                                        "minimax_tilting", "rejection", "gibbs"
                                      ),
                                      n_new_latent_gibbs = 100L,
                                      rejection_batch_size = NULL,
                                      rejection_max_batches = 1000L,
                                      seed = NULL) {
  if (!is.null(mixedgp_budget()) && !isTRUE(getOption("eivGP.prediction_block",FALSE))) {
    args <- lapply(as.list(match.call(expand.dots=TRUE))[-1L],eval,envir=parent.frame())
    return(mixedgp_predict_blocks(args,FUN=sample_eiv_m_given_xc_fb,object_name="fit_obj"))
  }

  latent_sampler <- match.arg(latent_sampler)
  if (!is.null(seed)) set.seed(as.integer(seed))
  p <- fit_obj$data$p
  d <- fit_obj$data$d
  q <- fit_obj$data$q
  m_vec <- fit_obj$data$m_vec
  X_input_names <- if (is.null(dim(X_test_raw))) names(X_test_raw) else colnames(X_test_raw)
  X_test_raw <- as_numeric_matrix_strict(
    X_test_raw, "X_test_raw", ncol_expected = p
  )
  train_x_names <- colnames(fit_obj$data$X_raw)
  if (!is.null(X_input_names) && all(nzchar(X_input_names)) &&
      !is.null(train_x_names)) {
    if (!setequal(X_input_names, train_x_names)) {
      stop("X_test_raw column names do not match the training X columns.")
    }
    X_test_raw <- X_test_raw[, train_x_names, drop = FALSE]
  }
  level_maps <- fit_obj$data$C_level_maps
  if (is.null(level_maps)) level_maps <- lapply(m_vec, seq_len)
  C_names <- fit_obj$data$C_names
  if (is.null(C_names)) C_names <- colnames(fit_obj$data$C_ord)
  C_test <- prepare_ordinal_matrix(
    C_test,
    m_vec = m_vec,
    level_maps = level_maps,
    expected_names = C_names,
    name = "C_test"
  )$C
  N <- nrow(C_test)
  if (nrow(X_test_raw) != N || ncol(C_test) != q) {
    stop("X_test_raw and C_test have incompatible dimensions.")
  }
  X_test <- sweep(
    sweep(X_test_raw, 2L, fit_obj$data$X_center, "-"),
    2L, fit_obj$data$X_scale, "/"
  )

  n_saved <- nrow(fit_obj$mcmc$samples_logtheta)
  if (is.null(draw_ids)) draw_ids <- seq_len(n_saved)
  draw_ids <- as.integer(draw_ids)
  if (length(draw_ids) < 1L || anyNA(draw_ids) ||
      any(draw_ids < 1L | draw_ids > n_saved)) {
    stop("draw_ids contains invalid posterior-draw indices.")
  }
  n_per_draw <- as.integer(n_per_draw)
  n_latent <- as.integer(n_latent)
  if (length(n_per_draw) != 1L || is.na(n_per_draw) || n_per_draw < 1L ||
      length(n_latent) != 1L || is.na(n_latent) || n_latent < 1L) {
    stop("n_per_draw and n_latent must be positive integers.")
  }

  kernel_spec <- kernel_spec_from_fit(fit_obj)
  n_train <- nrow(fit_obj$data$X)
  out <- matrix(NA_real_, length(draw_ids) * n_per_draw, N)
  conditional_means <- matrix(NA_real_, length(draw_ids), N)
  conditional_vars <- matrix(NA_real_, length(draw_ids), N)
  acceptance <- numeric(length(draw_ids))
  row_id <- 0L
  C_repeated <- C_test[rep(seq_len(N), each = n_latent), , drop = FALSE]
  for (ii in seq_along(draw_ids)) {
    s <- draw_ids[ii]
    A_s <- matrix(fit_obj$mcmc$samples_A[s, , ], nrow = q, ncol = d)
    tau_s <- unflatten_tau(fit_obj$mcmc$samples_tau[s, ], m_vec)
    U_flat <- sample_u_given_c_ordprobit_dispatch(
      C_new = C_repeated,
      A = A_s,
      tau = tau_s,
      m_vec = m_vec,
      latent_sampler = latent_sampler,
      n_gibbs = n_new_latent_gibbs,
      rejection_batch_size = rejection_batch_size,
      rejection_max_batches = rejection_max_batches
    )
    acceptance[ii] <- if (latent_sampler == "rejection") {
      attr(U_flat, "overall_acceptance")
    } else {
      NA_real_
    }
    U_mc <- array(NA_real_, dim = c(n_latent, N, d))
    for (k in seq_len(d)) {
      U_mc[, , k] <- matrix(U_flat[, k], nrow = n_latent, ncol = N)
    }
    U_train_s <- matrix(
      fit_obj$mcmc$samples_U[s, , ], nrow = n_train, ncol = d
    )
    pred <- gp_integrated_mean_state_general(
      X_train = fit_obj$data$X,
      U_train = U_train_s,
      y_train = fit_obj$data$y,
      X_star = X_test,
      U_mc = U_mc,
      logtheta = fit_obj$mcmc$samples_logtheta[s, ],
      sigma2_eps = fit_obj$mcmc$samples_sigma2[s],
      kernel = kernel_spec$name,
      matern_nu = kernel_spec$matern_nu,
      return_cov = isTRUE(joint) && isTRUE(include_process_uncertainty)
    )
    conditional_means[ii, ] <-
      fit_obj$data$y_center + fit_obj$data$y_scale * pred$mean
    conditional_vars[ii, ] <- fit_obj$data$y_scale^2 * pred$var
    for (rr in seq_len(n_per_draw)) {
      row_id <- row_id + 1L
      value_std <- if (!isTRUE(include_process_uncertainty)) {
        pred$mean
      } else if (isTRUE(joint)) {
        as.numeric(rmvnorm_psd(1L, pred$mean, pred$cov))
      } else {
        pred$mean + sqrt(pred$var) * rnorm(N)
      }
      out[row_id, ] <- fit_obj$data$y_center + fit_obj$data$y_scale * value_std
    }
  }
  attr(out, "target") <- "m(x,c)"
  attr(out, "n_latent") <- n_latent
  attr(out, "joint") <- isTRUE(joint)
  attr(out, "include_process_uncertainty") <- isTRUE(include_process_uncertainty)
  attr(out, "latent_sampler") <- latent_sampler
  attr(out, "rejection_acceptance") <- acceptance
  if (isTRUE(return_components)) {
    attr(out, "conditional_means") <- conditional_means
    attr(out, "conditional_vars") <- conditional_vars
  }
  out
}


sample_eiv_u_given_c_ordprobit <- function(C_new,
                                            fit_obj,
                                            draw_ids = NULL,
                                            n_per_draw = 1L,
                                            scale = c("auto", "raw", "model"),
                                            latent_sampler = c(
                                              "minimax_tilting", "rejection", "gibbs"
                                            ),
                                            n_new_latent_gibbs = 100L,
                                            rejection_batch_size = NULL,
                                            rejection_max_batches = 1000L,
                                            seed = NULL) {
  scale <- match.arg(scale)
  anchored <- isTRUE(fit_obj$data$latent_scale_anchored)
  if (scale == "raw" && !anchored) {
    stop("scale='raw' requires calibration data that anchor latent U.")
  }
  if (scale == "auto") scale <- if (anchored) "raw" else "model"
  latent_sampler <- match.arg(latent_sampler)
  if (!is.null(seed)) set.seed(as.integer(seed))
  m_vec <- fit_obj$data$m_vec
  q <- fit_obj$data$q
  d <- fit_obj$data$d
  level_maps <- fit_obj$data$C_level_maps
  if (is.null(level_maps)) level_maps <- lapply(m_vec, seq_len)
  C_names <- fit_obj$data$C_names
  if (is.null(C_names)) C_names <- colnames(fit_obj$data$C_ord)
  C_new <- prepare_ordinal_matrix(
    C_new,
    m_vec = m_vec,
    level_maps = level_maps,
    expected_names = C_names,
    name = "C_new"
  )$C
  if (ncol(C_new) != q) stop("C_new has the wrong number of ordinal inputs.")
  n_saved <- nrow(fit_obj$mcmc$samples_logtheta)
  if (is.null(draw_ids)) draw_ids <- seq_len(n_saved)
  draw_ids <- as.integer(draw_ids)
  if (length(draw_ids) < 1L || anyNA(draw_ids) ||
      any(draw_ids < 1L | draw_ids > n_saved)) {
    stop("draw_ids contains invalid posterior-draw indices.")
  }
  n_per_draw <- as.integer(n_per_draw)
  if (length(n_per_draw) != 1L || is.na(n_per_draw) || n_per_draw < 1L) {
    stop("n_per_draw must be a positive integer.")
  }

  N <- nrow(C_new)
  out <- array(
    NA_real_,
    dim = c(length(draw_ids) * n_per_draw, N, d),
    dimnames = list(NULL, NULL, fit_obj$data$U_names)
  )
  acceptance <- numeric(length(draw_ids) * n_per_draw)
  row_id <- 0L
  for (s in draw_ids) {
    A_s <- matrix(fit_obj$mcmc$samples_A[s, , ], nrow = q, ncol = d)
    tau_s <- unflatten_tau(fit_obj$mcmc$samples_tau[s, ], m_vec)
    for (rr in seq_len(n_per_draw)) {
      row_id <- row_id + 1L
      U_new <- sample_u_given_c_ordprobit_dispatch(
        C_new = C_new,
        A = A_s,
        tau = tau_s,
        m_vec = m_vec,
        latent_sampler = latent_sampler,
        n_gibbs = n_new_latent_gibbs,
        rejection_batch_size = rejection_batch_size,
        rejection_max_batches = rejection_max_batches
      )
      out[row_id, , ] <- U_new
      acceptance[row_id] <- if (latent_sampler == "rejection") {
        attr(U_new, "overall_acceptance")
      } else {
        NA_real_
      }
    }
  }
  if (scale == "raw") {
    center <- fit_obj$data$U_center
    if (is.null(center)) center <- rep(0, d)
    scale_value <- fit_obj$data$U_scale
    if (is.null(scale_value)) scale_value <- rep(1, d)
    for (k in seq_len(d)) {
      out[, , k] <- center[k] + scale_value[k] * out[, , k]
    }
  }
  attr(out, "source") <- "prospective"
  attr(out, "latent_sampler") <- latent_sampler
  attr(out, "rejection_acceptance") <- acceptance
  attr(out, "scale") <- if (!anchored && scale == "model") "working" else scale
  attr(out, "latent_scale_anchored") <- anchored
  out
}


############################################################
## Study II synthetic data-generating mechanism
############################################################

STUDY2_DESIGN_TAG <- paste(
  "study2-manuscript-v10-tailstable",
  "q4-d2",
  "A-fixed",
  "balanced-4-level",
  "sigma0.12",
  "exact-minimax-interwoven",
  "common-random-numbers",
  sep = "_"
)


study2_scenario_spec <- function(
    scenario = c(
      "primary",
      "latent_additive_control",
      "high_uncertainty",
      "logistic_misspec"
    )) {
  scenario <- match.arg(scenario)

  switch(
    scenario,
    primary = list(
      scenario = scenario,
      response_interactions = TRUE,
      score_error = "gaussian",
      lambda = 1.0
    ),
    latent_additive_control = list(
      scenario = scenario,
      response_interactions = FALSE,
      score_error = "gaussian",
      lambda = 1.0
    ),
    high_uncertainty = list(
      scenario = scenario,
      response_interactions = TRUE,
      score_error = "gaussian",
      lambda = 1.5
    ),
    logistic_misspec = list(
      scenario = scenario,
      response_interactions = TRUE,
      score_error = "logistic",
      lambda = 1.0
    )
  )
}


make_study2_balanced_tau <- function(A, lambda, score_error, m = 4L) {
  m <- as.integer(m)
  if (length(m) != 1L || is.na(m) || m < 2L) {
    stop("m must be one integer at least two.")
  }
  probs <- seq_len(m - 1L) / m

  if (score_error == "gaussian") {
    score_sd <- sqrt(rowSums(A^2) + lambda^2)
    return(outer(score_sd, qnorm(probs)))
  }

  ## For the robustness setting, each score is the convolution of a centered
  ## Gaussian loading term and a variance-one logistic error. Compute fixed
  ## marginal quantiles numerically so that changing the error distribution is
  ## not confounded with category imbalance.
  logistic_scale <- sqrt(3) / base::pi

  tau <- t(vapply(rowSums(A^2), function(var_normal) {
    cdf_score <- function(t) {
      if (var_normal <= .Machine$double.eps) {
        return(plogis(t, scale = lambda * logistic_scale))
      }

      integrate(
        function(z) {
          plogis(
            (t - sqrt(var_normal) * z) / (lambda * logistic_scale)
          ) * dnorm(z)
        },
        lower = -9,
        upper = 9,
        rel.tol = 1e-9
      )$value
    }

    vapply(probs, function(prob) {
      if (prob == 0.5) return(0)
      uniroot(
        function(t) cdf_score(t) - prob,
        interval = c(-10, 10),
        tol = 1e-9
      )$root
    }, numeric(1))
  }, numeric(length(probs))))

  tau
}


make_study2_loading_matrix <- function(q = 4L) {
  q <- as.integer(q)
  if (length(q) != 1L || is.na(q) || q < 2L || q > 6L) {
    stop("The publication proxy-dimension design requires 2 <= q <= 6.")
  }

  ## The first two rows identify the two latent axes. Later rows are mixed or
  ## redundant measurements of the same state. Designs with q = 2, 4, and 6
  ## are therefore nested rather than unrelated DGMs.
  A_master <- rbind(
    c(1.70, 0.00),
    c(0.20, 1.70),
    c(1.20, 0.70),
    c(0.70, 1.20),
    c(1.45, -0.45),
    c(-0.45, 1.45)
  )
  A_master[seq_len(q), , drop = FALSE]
}


make_study2_true_params <- function(scenario = "primary", q = 4L, m = 4L) {
  spec <- study2_scenario_spec(scenario)
  q <- as.integer(q)
  m <- as.integer(m)
  A <- make_study2_loading_matrix(q)
  Omega <- diag(q)
  tau <- make_study2_balanced_tau(
    A = A,
    lambda = spec$lambda,
    score_error = spec$score_error,
    m = m
  )

  c(
    spec,
    list(
      A = A,
      Omega = Omega,
      tau = tau,
      q = q,
      m = m,
      d = 2L,
      sigma_eps = 0.12
    )
  )
}


f0_2d <- function(X, U, scenario = "primary") {
  X <- as.matrix(X)
  U <- as.matrix(U)
  spec <- study2_scenario_spec(scenario)
  
  x1 <- X[, 1]
  x2 <- X[, 2]
  u1 <- U[, 1]
  u2 <- U[, 2]
  
  f_main <- 0.65 * tanh(1.10 * u1) +
    0.50 * tanh(0.90 * u2) +
    0.30 * x1 -
    0.25 * x2 +
    0.10 * x1 * x2

  if (!isTRUE(spec$response_interactions)) {
    return(f_main)
  }

  f_main +
    0.30 * tanh(0.60 * u1 * u2) +
    0.25 * x1 * tanh(u1) -
    0.20 * x2 * tanh(u2)
}


simulate_study2_data <- function(n = 120,
                                 n_test = 400,
                                 seed = NULL,
                                 scenario = "primary",
                                 sigma_eps = NULL,
                                 q = 4L,
                                 m = 4L) {

  n <- as.integer(n)
  n_test <- as.integer(n_test)
  if (n < 1L || n_test < 1L) {
    stop("n and n_test must be positive integers.")
  }
  
  seed <- if (is.null(seed)) sample.int(.Machine$integer.max, 1L) else as.integer(seed)
  pars <- make_study2_true_params(scenario = scenario, q = q, m = m)
  A <- pars$A
  tau <- pars$tau
  lambda <- pars$lambda
  d <- pars$d

  if (is.null(sigma_eps)) sigma_eps <- pars$sigma_eps
  
  ## Generate exactly once.  The experiment must retain naturally sparse
  ## realized patterns rather than condition on a favorable category count.
  set.seed(seed + 1L)
  X <- maximin_lhs_nd(n, d = 2L, lower = -1, upper = 1)
  set.seed(seed + 2L)
  U <- matrix(rnorm(n * d), n, d)
  set.seed(seed + 3L)
  eps_y <- rnorm(n, 0, sigma_eps)
  set.seed(seed + 4L)
  X_test <- matrix(runif(n_test * 2L, -1, 1), n_test, 2L)
  set.seed(seed + 5L)
  U_test <- matrix(rnorm(n_test * d), n_test, d)
  set.seed(seed + 6L)
  eps_y_test <- rnorm(n_test, 0, sigma_eps)

  set.seed(seed + 7L)
  Zeta <- draw_study2_score_errors(n, pars)
  S <- U %*% t(A) + lambda * Zeta
  C <- ordinal_from_scores_matrix_tau(S, tau)
  
  f <- f0_2d(X, U, scenario = scenario)
  y <- f + eps_y

  set.seed(seed + 8L)
  Zeta_test <- draw_study2_score_errors(n_test, pars)
  S_test <- U_test %*% t(A) + lambda * Zeta_test
  C_test <- ordinal_from_scores_matrix_tau(S_test, tau)
  
  f_test <- f0_2d(X_test, U_test, scenario = scenario)
  y_test <- f_test + eps_y_test
  
  list(
    train = list(
      X = X,
      U = U,
      C = C,
      y = y,
      f = f
    ),
    test = list(
      X = X_test,
      U = U_test,
      C = C_test,
      y = y_test,
      f = f_test
    ),
    true_params = pars,
    sigma_eps = sigma_eps,
    scenario = scenario
  )
}


classify_study2_pattern_frequency <- function(C_train,
                                              C_test,
                                              common_min = 5L) {
  train_key <- pattern_key(C_train)
  test_key <- pattern_key(C_test)
  train_counts <- table(train_key)
  n_train_pattern <- as.integer(train_counts[test_key])
  n_train_pattern[is.na(n_train_pattern)] <- 0L

  stratum <- ifelse(
    n_train_pattern == 0L,
    "unobserved",
    ifelse(n_train_pattern < common_min, "rare", "common")
  )

  data.frame(
    pattern = test_key,
    n_train_pattern = n_train_pattern,
    pattern_stratum = factor(
      stratum,
      levels = c("common", "rare", "unobserved")
    ),
    stringsAsFactors = FALSE
  )
}


make_stratified_calibration_sets_2d <- function(C,
                                                calib_grid = c(0, 10, 25, 50, 80),
                                                seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  
  C <- as.matrix(C)
  n <- nrow(C)
  calib_grid <- as.integer(calib_grid)
  if (ncol(C) < 2L) {
    stop("Study II calibration stratification requires at least two proxies.")
  }
  if (any(calib_grid < 0L | calib_grid > n)) {
    stop("Calibration sizes must lie between 0 and n.")
  }

  ## The first two proxies primarily anchor u1 and u2.  Randomize units within
  ## each observed (c1,c2) cell, randomize the order of the observed cells, and
  ## cycle through cells.  Prefixes of the resulting ordering give nested
  ## calibration sets without consulting y or the latent truth.
  cell <- interaction(C[, 1], C[, 2], drop = TRUE, lex.order = TRUE)
  idx_by_stratum <- split(seq_len(n), cell)
  idx_by_stratum <- idx_by_stratum[sample(seq_along(idx_by_stratum))]
  idx_by_stratum <- lapply(idx_by_stratum, function(idx) {
    if (length(idx) <= 1L) idx else sample(idx, size = length(idx), replace = FALSE)
  })
  
  ordering <- integer(0)
  
  repeat {
    added <- FALSE
    
    for (ss in seq_along(idx_by_stratum)) {
      if (length(idx_by_stratum[[ss]]) > 0) {
        ordering <- c(ordering, idx_by_stratum[[ss]][1])
        idx_by_stratum[[ss]] <- idx_by_stratum[[ss]][-1]
        added <- TRUE
      }
    }
    
    if (!added) break
  }
  
  out <- lapply(calib_grid, function(k) {
    if (k == 0) integer(0) else sort(ordering[seq_len(k)])
  })
  stopifnot(
    length(ordering) == n,
    length(unique(ordering)) == n,
    all(ordering %in% seq_len(n)),
    all(vapply(out, function(idx) length(idx) == length(unique(idx)), logical(1)))
  )
  
  names(out) <- as.character(calib_grid)
  out
}


make_nested_calibration_sets_2d <- function(C,
                                            calib_grid,
                                            seed = NULL,
                                            scheme = c("random", "anchor_stratified")) {
  scheme <- match.arg(scheme)
  C <- as.matrix(C)
  n <- nrow(C)
  calib_grid <- as.integer(calib_grid)
  if (anyNA(calib_grid) || any(calib_grid < 0L | calib_grid > n)) {
    stop("Calibration sizes must lie between zero and the training size.")
  }
  if (scheme == "anchor_stratified") {
    return(make_stratified_calibration_sets_2d(C, calib_grid, seed))
  }

  if (!is.null(seed)) set.seed(seed)
  ordering <- sample.int(n)
  out <- lapply(calib_grid, function(k) {
    if (k == 0L) integer(0) else sort(ordering[seq_len(k)])
  })
  names(out) <- as.character(calib_grid)
  out
}


############################################################
## Oracle prediction for Study II simulation
############################################################

study2_oracle_pool_provenance <- function(true_params) {
  list(
    schema_version = "study2-oracle-pool-1.0.0",
    scenario = as.character(true_params$scenario),
    score_error = as.character(true_params$score_error),
    lambda = as.numeric(true_params$lambda),
    A = unname(as.matrix(true_params$A)),
    tau = unname(as.matrix(true_params$tau)),
    Omega = unname(as.matrix(true_params$Omega)),
    d = as.integer(true_params$d),
    q = as.integer(true_params$q),
    m = as.integer(true_params$m)
  )
}


make_oracle_pool_2d <- function(true_params,
                                n_pool = 200000L,
                                seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  
  A <- true_params$A
  tau <- true_params$tau
  lambda <- true_params$lambda
  
  d <- true_params$d
  U <- matrix(rnorm(n_pool * d), n_pool, d)
  Zeta <- draw_study2_score_errors(n_pool, true_params)
  
  S <- U %*% t(A) + lambda * Zeta
  C <- ordinal_from_scores_matrix_tau(S, tau)
  
  key <- pattern_key(C)
  split_idx <- split(seq_len(n_pool), key)
  
  list(
    U = U,
    C = C,
    key = key,
    split_idx = split_idx,
    provenance = study2_oracle_pool_provenance(true_params)
  )
}


sample_oracle_u_rejection <- function(c_star,
                                      true_params,
                                      n_draw,
                                      max_batches = 100L,
                                      batch_size = 20000L) {
  A <- true_params$A
  tau <- true_params$tau
  lambda <- true_params$lambda
  
  d <- true_params$d
  U_keep <- matrix(NA_real_, 0, d)
  key_star <- paste(c_star, collapse = "_")
  
  for (bb in seq_len(max_batches)) {
    U <- matrix(rnorm(batch_size * d), batch_size, d)
    Zeta <- draw_study2_score_errors(batch_size, true_params)
    S <- U %*% t(A) + lambda * Zeta
    C <- ordinal_from_scores_matrix_tau(S, tau)
    
    idx <- which(pattern_key(C) == key_star)
    
    if (length(idx) > 0) {
      U_keep <- rbind(U_keep, U[idx, , drop = FALSE])
    }
    
    if (nrow(U_keep) >= n_draw) break
  }
  
  if (nrow(U_keep) < n_draw) {
    stop(
      "Oracle rejection sampling produced only ", nrow(U_keep), " of ",
      n_draw, " required draws for ordinal pattern ",
      key_star,
      ". Increase max_batches or batch_size; do not bootstrap a short pool."
    )
  }
  
  U_keep[seq_len(n_draw), , drop = FALSE]
}


sample_oracle_u_patterns_shared_2d <- function(
    C,
    true_params,
    n_per_pattern,
    oracle_pool = NULL,
    n_pool = 200000L,
    batch_size = 200000L,
    max_candidates = 20000000L,
    seed = NULL) {
  if (!is.null(seed)) set.seed(as.integer(seed))
  C <- as.matrix(C)
  if (nrow(C) < 1L || ncol(C) != true_params$q || anyNA(C) ||
      any(!is.finite(C)) || any(C != floor(C)) ||
      any(C < 1L | C > true_params$m)) {
    stop("C is incompatible with the Study II oracle model.")
  }

  keys <- unique(pattern_key(C))
  if (length(n_per_pattern) == 1L) {
    required <- setNames(rep(as.integer(n_per_pattern), length(keys)), keys)
  } else {
    if (is.null(names(n_per_pattern)) || !all(keys %in% names(n_per_pattern))) {
      stop("A vector n_per_pattern must be named for every requested pattern.")
    }
    required <- as.integer(n_per_pattern[keys])
    names(required) <- keys
  }
  if (anyNA(required) || any(required < 1L)) {
    stop("Every requested oracle pattern must require at least one draw.")
  }

  n_pool <- as.integer(n_pool)
  batch_size <- as.integer(batch_size)
  max_candidates <- as.numeric(max_candidates)
  if (is.na(n_pool) || n_pool < 1L || is.na(batch_size) || batch_size < 1L ||
      !is.finite(max_candidates) || max_candidates < batch_size) {
    stop("Oracle rejection budgets are invalid.")
  }
  if (is.null(oracle_pool)) {
    oracle_pool <- make_oracle_pool_2d(true_params, n_pool = n_pool)
  }
  validate_oracle_pool_2d(oracle_pool, true_params)

  accepted <- setNames(vector("list", length(keys)), keys)
  accepted_n <- setNames(integer(length(keys)), keys)
  proposal_hits <- setNames(integer(length(keys)), keys)
  initial_candidates <- nrow(oracle_pool$U)
  initial_pool_hits <- setNames(integer(length(keys)), keys)
  for (key in keys) {
    idx <- oracle_pool$split_idx[[key]]
    if (is.null(idx)) idx <- integer(0)
    initial_pool_hits[[key]] <- length(idx)
    proposal_hits[[key]] <- length(idx)
    if (length(idx) > 0L) {
      take <- idx[seq_len(min(length(idx), required[[key]]))]
      accepted[[key]] <- list(oracle_pool$U[take, , drop = FALSE])
      accepted_n[[key]] <- length(take)
    }
  }

  additional_candidates <- 0
  while (any(accepted_n < required)) {
    if (additional_candidates + batch_size > max_candidates) {
      missing <- keys[accepted_n < required]
      detail <- paste0(
        missing, " (", accepted_n[missing], "/", required[missing], ")"
      )
      stop(
        "Shared exact oracle rejection reached max_candidates before filling ",
        "pattern(s): ", paste(detail, collapse = ", "), "."
      )
    }
    U_batch <- matrix(
      rnorm(batch_size * true_params$d), batch_size, true_params$d
    )
    Zeta_batch <- draw_study2_score_errors(batch_size, true_params)
    S_batch <- U_batch %*% t(true_params$A) +
      true_params$lambda * Zeta_batch
    C_batch <- ordinal_from_scores_matrix_tau(S_batch, true_params$tau)
    batch_split <- split(seq_len(batch_size), pattern_key(C_batch))
    additional_candidates <- additional_candidates + batch_size

    for (key in keys) {
      idx <- batch_split[[key]]
      if (is.null(idx) || length(idx) == 0L) next
      proposal_hits[[key]] <- proposal_hits[[key]] + length(idx)
      need <- required[[key]] - accepted_n[[key]]
      if (need <= 0L) next
      take <- idx[seq_len(min(need, length(idx)))]
      accepted[[key]][[length(accepted[[key]]) + 1L]] <-
        U_batch[take, , drop = FALSE]
      accepted_n[[key]] <- accepted_n[[key]] + length(take)
    }
  }

  out <- setNames(lapply(keys, function(key) {
    do.call(rbind, accepted[[key]])[
      seq_len(required[[key]]), , drop = FALSE
    ]
  }), keys)
  attr(out, "rejection_telemetry") <- data.frame(
    pattern = keys,
    required = as.integer(required[keys]),
    retained = as.integer(accepted_n[keys]),
    initial_pool_hits = as.integer(initial_pool_hits[keys]),
    initial_candidates = initial_candidates,
    additional_candidates = additional_candidates,
    total_candidates = initial_candidates + additional_candidates,
    empirical_acceptance = as.integer(proposal_hits[keys]) /
      (initial_candidates + additional_candidates),
    stringsAsFactors = FALSE
  )
  out
}


oracle_m0_2d <- function(X,
                         C,
                         true_params,
                         oracle_pool = NULL,
                         n_latent = 5000L,
                         batch_size = 200000L,
                         max_candidates = 20000000L,
                         seed = NULL) {
  X <- as.matrix(X)
  C <- as.matrix(C)
  n_latent <- as.integer(n_latent)
  if (nrow(X) != nrow(C) || ncol(C) != true_params$q ||
      n_latent < 1L || any(!is.finite(X)) || anyNA(C)) {
    stop("Incompatible inputs in oracle_m0_2d().")
  }
  keys <- pattern_key(C)
  unique_keys <- unique(keys)
  pattern_draws <- sample_oracle_u_patterns_shared_2d(
    C = C,
    true_params = true_params,
    n_per_pattern = n_latent,
    oracle_pool = oracle_pool,
    n_pool = max(200000L, 20L * n_latent),
    batch_size = batch_size,
    max_candidates = max_candidates,
    seed = seed
  )

  out <- numeric(nrow(X))
  truth_mcse <- numeric(nrow(X))
  nested_half_difference <- numeric(nrow(X))
  for (key in unique_keys) {
    rows <- which(keys == key)
    U_draw <- pattern_draws[[key]]
    for (i in rows) {
      X_rep <- X[rep(i, n_latent), , drop = FALSE]
      f_value <- f0_2d(X_rep, U_draw, scenario = true_params$scenario)
      out[i] <- mean(f_value)
      truth_mcse[i] <- if (n_latent > 1L) {
        stats::sd(f_value) / sqrt(n_latent)
      } else {
        NA_real_
      }
      n_half <- floor(n_latent / 2L)
      nested_half_difference[i] <- if (n_half > 0L) {
        abs(mean(f_value[seq_len(n_half)]) - out[i])
      } else {
        NA_real_
      }
    }
  }
  attr(out, "n_latent_truth") <- n_latent
  attr(out, "rejection_telemetry") <-
    attr(pattern_draws, "rejection_telemetry")
  attr(out, "truth_diagnostics") <- data.frame(
    evaluation_row = seq_len(nrow(X)),
    pattern = keys,
    n_latent = n_latent,
    mcse = truth_mcse,
    nested_half_difference = nested_half_difference,
    stringsAsFactors = FALSE
  )
  out
}


summarize_mean_recovery_2d <- function(draw_mat,
                                       m_true,
                                       method,
                                       rep_id,
                                       n_calib,
                                       scenario,
                                       valid_function_draws = FALSE) {
  draw_mat <- as.matrix(draw_mat)
  m_true <- as.numeric(m_true)
  if (ncol(draw_mat) != length(m_true) || any(!is.finite(m_true))) {
    stop("draw_mat and m_true are incompatible in summarize_mean_recovery_2d().")
  }
  conditional_means <- attr(draw_mat, "conditional_means")
  post_mean <- if (!is.null(conditional_means)) {
    colMeans(conditional_means)
  } else {
    colMeans(draw_mat)
  }
  lo <- hi <- rep(NA_real_, length(m_true))
  if (isTRUE(valid_function_draws)) {
    qs <- apply(
      draw_mat, 2L, stats::quantile,
      probs = c(0.025, 0.975), names = FALSE
    )
    lo <- qs[1L, ]
    hi <- qs[2L, ]
  }
  error <- post_mean - m_true
  data.frame(
    rep = rep_id,
    scenario = scenario,
    n_calib = n_calib,
    method = method,
    n_eval = length(m_true),
    RMSE = sqrt(mean(error^2)),
    MAE = mean(abs(error)),
    Bias = mean(error),
    Coverage95 = if (isTRUE(valid_function_draws)) {
      mean(m_true >= lo & m_true <= hi)
    } else {
      NA_real_
    },
    Width95 = if (isTRUE(valid_function_draws)) mean(hi - lo) else NA_real_,
    uncertainty_target = if (isTRUE(valid_function_draws)) "m(x,c)" else "none",
    stringsAsFactors = FALSE
  )
}


sample_oracle_test_y_2d <- function(X_test,
                                    C_test,
                                    true_params,
                                    sigma_eps,
                                    n_draw = 1000L,
                                    oracle_pool = NULL,
                                    n_pool = 200000L,
                                    batch_size = 200000L,
                                    max_candidates = 20000000L,
                                    seed = NULL) {
  X_test <- as.matrix(X_test)
  C_test <- as.matrix(C_test)
  n_draw <- as.integer(n_draw)
  if (nrow(X_test) != nrow(C_test) || ncol(C_test) != true_params$q ||
      n_draw < 1L || length(sigma_eps) != 1L || !is.finite(sigma_eps) ||
      sigma_eps < 0 || any(!is.finite(X_test)) || anyNA(C_test)) {
    stop("Incompatible inputs in sample_oracle_test_y_2d().")
  }

  n_test <- nrow(C_test)
  out <- matrix(NA_real_, nrow = n_draw, ncol = n_test)
  conditional_means <- matrix(NA_real_, nrow = n_draw, ncol = n_test)
  conditional_vars <- matrix(
    sigma_eps^2, nrow = n_draw, ncol = n_test
  )
  keys_test <- pattern_key(C_test)
  pattern_counts <- table(factor(keys_test, levels = unique(keys_test)))
  required <- n_draw * as.integer(pattern_counts)
  names(required) <- names(pattern_counts)
  pattern_draws <- sample_oracle_u_patterns_shared_2d(
    C = C_test,
    true_params = true_params,
    n_per_pattern = required,
    oracle_pool = oracle_pool,
    n_pool = n_pool,
    batch_size = batch_size,
    max_candidates = max_candidates,
    seed = seed
  )
  cursor <- setNames(integer(length(required)), names(required))

  for (j in seq_len(n_test)) {
    key_j <- keys_test[j]
    draw_idx <- cursor[[key_j]] + seq_len(n_draw)
    U_draw <- pattern_draws[[key_j]][draw_idx, , drop = FALSE]
    cursor[[key_j]] <- cursor[[key_j]] + n_draw

    X_rep <- matrix(
      rep(X_test[j, ], each = n_draw),
      nrow = n_draw,
      ncol = ncol(X_test)
    )
    f_draw <- f0_2d(
      X_rep,
      U_draw,
      scenario = true_params$scenario
    )
    
    out[, j] <- f_draw + rnorm(n_draw, 0, sigma_eps)
    conditional_means[, j] <- f_draw
  }

  attr(out, "latent_sampler") <- "exact_shared_rejection"
  attr(out, "rejection_telemetry") <-
    attr(pattern_draws, "rejection_telemetry")
  attr(out, "conditional_means") <- conditional_means
  attr(out, "conditional_vars") <- conditional_vars
  attr(out, "mixture_components") <- "DGM conditional on Monte Carlo U*"
  out
}


## Internal small-data GP engine used by the appendix Full-U, PI, and CC
## benchmarks as well as the archived deterministic-embedding pilots.  Every
## fit uses reproducible multistart optimization and fails closed unless at
## least one L-BFGS-B run reports convergence.
gp_mle_fit <- function(X,
                       y,
                       n_starts = 5L,
                       seed = 1L,
                       maxit = 500L) {
  X <- as.matrix(X)
  storage.mode(X) <- "double"
  y <- as.numeric(y)

  n <- nrow(X)
  d <- ncol(X)
  if (n < 2L || d < 1L || length(y) != n || any(!is.finite(X)) ||
      any(!is.finite(y))) {
    stop("X and y must contain compatible finite GP training data.")
  }
  n_starts <- as.integer(n_starts)
  seed <- as.integer(seed)
  maxit <- as.integer(maxit)
  if (length(n_starts) != 1L || is.na(n_starts) || n_starts < 1L ||
      length(seed) != 1L || is.na(seed) ||
      length(maxit) != 1L || is.na(maxit) || maxit < 0L) {
    stop("n_starts, seed, and maxit must be finite nonnegative integers, with n_starts >= 1.")
  }
  
  Dlist <- lapply(seq_len(d), function(j) {
    pairwise_sqdist(X[, j, drop = FALSE])
  })
  
  nll <- function(par) {
    log_sigma2 <- par[1]
    log_rho <- par[2]
    log_theta <- par[-c(1, 2)]
    
    sigma2 <- exp(log_sigma2)
    rho <- exp(log_rho)
    theta <- exp(log_theta)
    
    Rexp <- matrix(0, n, n)
    
    for (j in seq_len(d)) {
      Rexp <- Rexp + theta[j] * Dlist[[j]]
    }
    
    R <- exp(-Rexp)
    A_mat <- rho^2 * R + diag(n)
    
    Uchol <- try(safe_chol(A_mat), silent = TRUE)
    if (inherits(Uchol, "try-error")) return(1e20)
    
    Ainv_y <- solve_chol(Uchol, y)
    
    logdetA <- 2 * sum(log(diag(Uchol)))
    quad <- sum(y * Ainv_y)
    
    0.5 * (n * log(2 * base::pi * sigma2) + logdetA + quad / sigma2)
  }
  
  lower <- c(log(1e-5), log(0.05), rep(log(1e-4), d))
  upper <- c(log(5), log(100), rep(log(100), d))

  caller_rng <- mixedgp_rng_state()
  on.exit(mixedgp_restore_rng_state(caller_rng), add = TRUE)
  set.seed(seed)
  initial <- matrix(NA_real_, nrow = n_starts, ncol = length(lower))
  initial[1L, ] <- c(log(0.05), log(3), rep(log(0.5), d))
  if (n_starts > 1L) {
    for (ss in 2:n_starts) {
      initial[ss, ] <- stats::runif(length(lower), lower, upper)
    }
  }

  attempts <- lapply(seq_len(n_starts), function(ss) {
    tryCatch(
      stats::optim(
        par = initial[ss, ],
        fn = nll,
        method = "L-BFGS-B",
        lower = lower,
        upper = upper,
        control = list(maxit = maxit)
      ),
      error = function(e) e
    )
  })
  parameter_names <- c(
    "log_sigma2", "log_rho", paste0("log_theta", seq_len(d))
  )
  initial_report <- as.data.frame(initial, stringsAsFactors = FALSE)
  names(initial_report) <- paste0("initial_", parameter_names)
  final_report <- as.data.frame(
    matrix(NA_real_, nrow = n_starts, ncol = length(parameter_names)),
    stringsAsFactors = FALSE
  )
  names(final_report) <- paste0("final_", parameter_names)
  report <- cbind(
    data.frame(
      start = seq_len(n_starts),
      n_starts = rep.int(n_starts, n_starts),
      optimizer_seed = rep.int(seed, n_starts),
      maxit = rep.int(maxit, n_starts),
      stringsAsFactors = FALSE
    ),
    initial_report,
    final_report,
    data.frame(
      convergence = rep.int(NA_integer_, n_starts),
      objective = rep.int(Inf, n_starts),
      function_evaluations = rep.int(NA_integer_, n_starts),
      gradient_evaluations = rep.int(NA_integer_, n_starts),
      message = rep.int("", n_starts),
      stringsAsFactors = FALSE
    )
  )
  for (ss in seq_along(attempts)) {
    attempt <- attempts[[ss]]
    if (inherits(attempt, "error")) {
      report$message[ss] <- conditionMessage(attempt)
    } else {
      report$convergence[ss] <- as.integer(attempt$convergence)
      report$objective[ss] <- as.numeric(attempt$value)
      if (!is.null(attempt$counts)) {
        function_count <- attempt$counts[["function"]]
        gradient_count <- attempt$counts[["gradient"]]
        if (length(function_count) == 1L) {
          report$function_evaluations[ss] <- as.integer(function_count)
        }
        if (length(gradient_count) == 1L) {
          report$gradient_evaluations[ss] <- as.integer(gradient_count)
        }
      }
      if (length(attempt$par) == length(parameter_names)) {
        report[ss, paste0("final_", parameter_names)] <- attempt$par
      }
      attempt_message <- attempt$message
      report$message[ss] <- if (
        is.null(attempt_message) || length(attempt_message) == 0L
      ) {
        ""
      } else {
        as.character(attempt_message[[1L]])
      }
    }
  }
  valid <- which(
    report$convergence == 0L & is.finite(report$objective) &
      report$objective < 1e19
  )
  if (length(valid) == 0L) {
    failure <- structure(
      list(
        message = paste0(
          "Internal Study II GP MLE did not converge from any of the ",
          n_starts, " prespecified starts."
        ),
        call = sys.call(),
        optimizer_attempts = report,
        optimizer_control = list(
          n_starts = n_starts, seed = seed, maxit = maxit
        )
      ),
      class = c("mixedgp_gp_mle_error", "error", "condition")
    )
    stop(failure)
  }
  best <- valid[which.min(report$objective[valid])]
  opt <- attempts[[best]]

  list(
    par = opt$par,
    value = opt$value,
    convergence = opt$convergence,
    optimizer_attempts = report,
    selected_start = best,
    optimizer_control = list(
      n_starts = n_starts, seed = seed, maxit = maxit
    ),
    X = X,
    y = y
  )
}


gp_mle_fit_learned_embedding_ord <- function(X,
                                             C_ord,
                                             y,
                                             m_vec,
                                             n_starts = 4L) {
  X <- as.matrix(X)
  C_ord <- as.matrix(C_ord)
  
  n <- nrow(X)
  p <- ncol(X)
  q <- ncol(C_ord)
  
  n_feat <- p + q
  n_embed <- sum(pmax(m_vec - 2L, 0L))
  
  nll <- function(par) {
    log_sigma2 <- par[1]
    log_rho <- par[2]
    log_theta <- par[3:(2 + n_feat)]
    
    sigma2 <- exp(log_sigma2)
    rho <- exp(log_rho)
    theta <- exp(log_theta)
    
    embed_par <- if (n_embed > 0L) {
      par[(3 + n_feat):length(par)]
    } else {
      numeric(0)
    }
    
    Z <- matrix(NA_real_, n, q)
    
    pos <- 0L
    for (j in seq_len(q)) {
      mj <- m_vec[j]
      
      if (mj > 2L) {
        a_j <- embed_par[pos + seq_len(mj - 2L)]
        pos <- pos + (mj - 2L)
      } else {
        a_j <- numeric(0)
      }
      
      z_scores <- make_monotone_scores(a_j, mj)
      Z[, j] <- z_scores[C_ord[, j]]
    }
    
    X_aug <- cbind(X, Z)
    
    Rexp <- matrix(0, n, n)
    
    for (j in seq_len(n_feat)) {
      Rexp <- Rexp + theta[j] * pairwise_sqdist(X_aug[, j, drop = FALSE])
    }
    
    R <- exp(-Rexp)
    A_mat <- rho^2 * R + diag(n)
    
    Uchol <- try(safe_chol(A_mat), silent = TRUE)
    if (inherits(Uchol, "try-error")) return(1e20)
    
    Ainv_y <- solve_chol(Uchol, y)
    
    logdetA <- 2 * sum(log(diag(Uchol)))
    quad <- sum(y * Ainv_y)
    
    0.5 * (n * log(2 * base::pi * sigma2) + logdetA + quad / sigma2)
  }
  
  lower <- c(
    log(1e-5),
    log(0.05),
    rep(log(1e-4), n_feat),
    rep(-6, n_embed)
  )
  
  upper <- c(
    log(5),
    log(100),
    rep(log(100), n_feat),
    rep(6, n_embed)
  )
  
  make_init <- function() {
    c(
      log(0.05),
      log(3),
      rep(log(0.5), n_feat),
      rnorm(n_embed, 0, 0.5)
    )
  }
  
  starts <- replicate(n_starts, make_init(), simplify = FALSE)
  
  opts <- lapply(starts, function(init) {
    optim(
      par = init,
      fn = nll,
      method = "L-BFGS-B",
      lower = lower,
      upper = upper,
      control = list(maxit = 500)
    )
  })
  
  best_id <- which.min(vapply(opts, function(z) z$value, numeric(1)))
  opt <- opts[[best_id]]
  
  embed_par <- if (n_embed > 0L) {
    opt$par[(3 + n_feat):length(opt$par)]
  } else {
    numeric(0)
  }
  
  z_scores_list <- vector("list", q)
  
  pos <- 0L
  for (j in seq_len(q)) {
    mj <- m_vec[j]
    
    if (mj > 2L) {
      a_j <- embed_par[pos + seq_len(mj - 2L)]
      pos <- pos + (mj - 2L)
    } else {
      a_j <- numeric(0)
    }
    
    z_scores_list[[j]] <- make_monotone_scores(a_j, mj)
  }
  
  Z <- matrix(NA_real_, n, q)
  for (j in seq_len(q)) {
    Z[, j] <- z_scores_list[[j]][C_ord[, j]]
  }
  
  list(
    par = opt$par[1:(2 + n_feat)],
    embed_par = embed_par,
    z_scores_list = z_scores_list,
    value = opt$value,
    convergence = opt$convergence,
    X_aug = cbind(X, Z),
    X = X,
    C_ord = C_ord,
    y = y,
    m_vec = m_vec
  )
}


fit_embedding_baselines_ord <- function(X_raw,
                                        y_raw,
                                        C_ord,
                                        m_vec,
                                        n_starts_learned = 4L) {
  X_raw <- as.matrix(X_raw)
  C_ord <- as.matrix(C_ord)
  
  X_center <- colMeans(X_raw)
  X_scale <- apply(X_raw, 2, sd)
  X_scale[!is.finite(X_scale) | X_scale <= 0] <- 1
  
  X <- sweep(sweep(X_raw, 2, X_center, "-"), 2, X_scale, "/")
  
  y_center <- mean(y_raw)
  y_scale <- sd(y_raw)
  
  y <- as.numeric((y_raw - y_center) / y_scale)
  
  q <- ncol(C_ord)
  
  Z_gauss <- matrix(NA_real_, nrow(C_ord), q)
  
  for (j in seq_len(q)) {
    Z_gauss[, j] <- qnorm(C_ord[, j] / (m_vec[j] + 1))
  }
  
  fit_gauss <- gp_mle_fit(
    X = cbind(X, Z_gauss),
    y = y
  )
  
  z_cm_scores <- lapply(seq_len(q), function(j) {
    conditional_mean_scores_ord(C_ord[, j], m_vec[j])
  })
  
  Z_cm <- matrix(NA_real_, nrow(C_ord), q)
  
  for (j in seq_len(q)) {
    Z_cm[, j] <- z_cm_scores[[j]][C_ord[, j]]
  }
  
  fit_cm <- gp_mle_fit(
    X = cbind(X, Z_cm),
    y = y
  )
  
  fit_learned <- gp_mle_fit_learned_embedding_ord(
    X = X,
    C_ord = C_ord,
    y = y,
    m_vec = m_vec,
    n_starts = n_starts_learned
  )
  
  list(
    X_center = X_center,
    X_scale = X_scale,
    y_center = y_center,
    y_scale = y_scale,
    fit_gauss = fit_gauss,
    fit_cm = fit_cm,
    fit_learned = fit_learned,
    z_cm_scores = z_cm_scores,
    m_vec = m_vec
  )
}


############################################################
## Predictive scoring
############################################################

crps_sample_one <- function(draws, y) {
  draws <- sort(as.numeric(draws))
  S <- length(draws)
  
  term1 <- mean(abs(draws - y))
  
  weights <- 2 * seq_len(S) - S - 1
  mean_abs_pair <- 2 * sum(weights * draws) / S^2
  
  term1 - 0.5 * mean_abs_pair
}


crps_sample_matrix <- function(draw_mat, y) {
  vapply(seq_along(y), function(j) {
    crps_sample_one(draw_mat[, j], y[j])
  }, numeric(1))
}


normal_mixture_nlpd <- function(y_true,
                                component_means,
                                component_vars) {
  y_true <- as.numeric(y_true)
  if (is.null(component_means) || is.null(component_vars)) {
    missing <- c(
      if (is.null(component_means)) "conditional_means" else NULL,
      if (is.null(component_vars)) "conditional_vars" else NULL
    )
    return(list(
      value = NA_real_,
      pointwise = rep(NA_real_, length(y_true)),
      reason = paste0(
        "predictive normal-mixture components unavailable: ",
        paste(missing, collapse = " and ")
      )
    ))
  }

  component_means <- as.matrix(component_means)
  component_vars <- as.matrix(component_vars)
  if (ncol(component_means) != length(y_true) ||
      ncol(component_vars) != length(y_true)) {
    stop("Predictive mixture components must have one column per test case.")
  }
  n_component <- max(nrow(component_means), nrow(component_vars))
  if (nrow(component_means) == 1L && n_component > 1L) {
    component_means <- component_means[
      rep(1L, n_component), , drop = FALSE
    ]
  }
  if (nrow(component_vars) == 1L && n_component > 1L) {
    component_vars <- component_vars[
      rep(1L, n_component), , drop = FALSE
    ]
  }
  if (nrow(component_means) != n_component ||
      nrow(component_vars) != n_component || n_component < 1L) {
    stop("Predictive component means and variances have incompatible rows.")
  }
  if (any(!is.finite(component_means)) ||
      any(!is.finite(component_vars)) || any(component_vars < 0)) {
    stop("Predictive normal-mixture components must have finite, nonnegative variances.")
  }
  if (any(component_vars == 0)) {
    return(list(
      value = NA_real_,
      pointwise = rep(NA_real_, length(y_true)),
      reason = paste0(
        "continuous NLPD undefined: predictive mixture contains ",
        "zero-variance components"
      )
    ))
  }

  pointwise <- vapply(seq_along(y_true), function(j) {
    log_component <- stats::dnorm(
      y_true[j],
      mean = component_means[, j],
      sd = sqrt(component_vars[, j]),
      log = TRUE
    )
    largest <- max(log_component)
    if (is.infinite(largest) && largest < 0) return(Inf)
    -(largest + log(mean(exp(log_component - largest))))
  }, numeric(1L))
  list(value = mean(pointwise), pointwise = pointwise, reason = "")
}


study2_latent_anchor_from_fit <- function(fit) {
  if (!is.list(fit) || !is.list(fit$data)) {
    stop("fit must contain a data component.")
  }
  rank <- fit$data$latent_anchor_rank
  required_rank <- fit$data$latent_anchor_required_rank
  anchored <- fit$data$latent_scale_anchored
  if (!is.null(rank) && !is.null(required_rank) && !is.null(anchored)) {
    return(list(
      anchored = isTRUE(anchored),
      affine_rank = as.integer(rank),
      required_rank = as.integer(required_rank),
      n_calibrated = length(fit$data$calib_idx)
    ))
  }
  if (is.null(fit$data$U_obs) || is.null(fit$data$calib_idx) ||
      is.null(fit$data$d)) {
    stop("fit does not contain enough information to assess latent anchoring.")
  }
  mixedgp_latent_anchor_status(
    fit$data$U_obs,
    fit$data$calib_idx,
    d = fit$data$d
  )
}


study2_latent_imputation_status <- function(fit,
                                            method,
                                            target,
                                            rep_id,
                                            n_calib,
                                            scenario,
                                            n_units) {
  target <- match.arg(
    target,
    c("training_missing_U", "prospective_U_given_C")
  )
  n_units <- as.integer(n_units)
  if (length(n_units) != 1L || is.na(n_units) || n_units < 0L) {
    stop("n_units must be one nonnegative integer.")
  }
  anchor <- study2_latent_anchor_from_fit(fit)
  if (!isTRUE(anchor$anchored)) {
    status <- "not_identified"
    reason <- paste0(
      "Raw-coordinate latent error is not identified: calibration affine ",
      "rank is ", anchor$affine_rank, " but ", anchor$required_rank,
      " is required."
    )
  } else if (n_units == 0L) {
    status <- "not_applicable"
    reason <- if (target == "training_missing_U") {
      "Every training latent input is calibrated."
    } else {
      "No prospective units were supplied."
    }
  } else {
    status <- "eligible"
    reason <- ""
  }
  data.frame(
    rep = as.integer(rep_id),
    scenario = as.character(scenario),
    n_calib = as.integer(n_calib),
    method = as.character(method),
    target = target,
    status = status,
    task_eligible = identical(status, "eligible"),
    reason = reason,
    n_units = n_units,
    latent_scale = if (isTRUE(anchor$anchored)) "physical" else "working",
    calibration_affine_rank = as.integer(anchor$affine_rank),
    calibration_required_rank = as.integer(anchor$required_rank),
    stringsAsFactors = FALSE
  )
}


summarize_study2_latent_draws <- function(draws,
                                          U_true,
                                          method,
                                          target,
                                          rep_id,
                                          n_calib,
                                          scenario) {
  target <- match.arg(
    target,
    c("training_missing_U", "prospective_U_given_C")
  )
  U_true <- as.matrix(U_true)
  draw_dim <- dim(draws)
  if (length(draw_dim) != 3L || draw_dim[1L] < 1L ||
      draw_dim[2L] != nrow(U_true) || draw_dim[3L] != ncol(U_true)) {
    stop("draws must have dimensions draw by unit by latent coordinate.")
  }
  if (nrow(U_true) < 1L || ncol(U_true) < 1L ||
      any(!is.finite(U_true)) || any(!is.finite(draws))) {
    stop("Latent draws and truth must be nonempty and finite.")
  }

  out <- lapply(seq_len(ncol(U_true)), function(j) {
    draw_j <- draws[, , j, drop = FALSE]
    dim(draw_j) <- c(draw_dim[1L], draw_dim[2L])
    post_mean <- colMeans(draw_j)
    lo <- apply(draw_j, 2L, stats::quantile, probs = 0.025, names = FALSE)
    hi <- apply(draw_j, 2L, stats::quantile, probs = 0.975, names = FALSE)
    err <- post_mean - U_true[, j]

    data.frame(
      rep = as.integer(rep_id),
      scenario = as.character(scenario),
      n_calib = as.integer(n_calib),
      method = as.character(method),
      target = target,
      coordinate = paste0("u", j),
      n_units = nrow(U_true),
      n_missing = if (target == "training_missing_U") {
        nrow(U_true)
      } else {
        NA_integer_
      },
      Bias = mean(err),
      RMSE = sqrt(mean(err^2)),
      MAE = mean(abs(err)),
      Coverage95 = mean(U_true[, j] >= lo & U_true[, j] <= hi),
      Width95 = mean(hi - lo),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, out)
}


############################################################
## MCMC review helpers for ordinal-probit EIV-GP
############################################################

classify_parameter_block <- function(parameter) {
  out <- rep("other", length(parameter))
  
  out[grepl("^sigma|^rho|^theta", parameter)] <- "GP/noise"
  out[grepl("^A\\[|^A[0-9]", parameter)] <- "loading A"
  out[grepl("^tau\\[|^tau", parameter)] <- "cutpoints"
  out[grepl("^U\\[", parameter)] <- "latent U"
  
  out
}


summarise_rhat_block <- function(rhat_df) {
  rhat_df |>
    dplyr::filter(is.finite(rhat)) |>
    dplyr::group_by(block) |>
    dplyr::summarise(
      n = dplyr::n(),
      median = median(rhat),
      q90 = quantile(rhat, 0.90),
      q95 = quantile(rhat, 0.95),
      max = max(rhat),
      prop_gt_1.01 = mean(rhat > 1.01),
      prop_gt_1.05 = mean(rhat > 1.05),
      prop_gt_1.10 = mean(rhat > 1.10),
      .groups = "drop"
    )
}


summarise_ess_block <- function(ess_df) {
  ess_df |>
    dplyr::filter(is.finite(ess)) |>
    dplyr::group_by(block) |>
    dplyr::summarise(
      n = dplyr::n(),
      median = median(ess),
      q10 = quantile(ess, 0.10),
      q25 = quantile(ess, 0.25),
      min = min(ess),
      .groups = "drop"
    )
}


make_mcmc_review_ordprobit <- function(fit,
                                       label = "fit",
                                       out_dir = NULL,
                                       save_plots = TRUE,
                                       save_csvs = TRUE) {
  stopifnot(!is.null(fit$diagnostics))
  required <- c("dplyr", "tidyr")
  if (isTRUE(save_plots)) required <- c(required, "ggplot2")
  require_study2_reporting_packages(required, "MCMC review")
  
  rhat_parts <- list()
  
  if (!is.null(fit$diagnostics$rhat_hyper)) {
    rhat_parts$hyper <- fit$diagnostics$rhat_hyper |>
      dplyr::mutate(block = "GP/noise")
  }
  
  if (!is.null(fit$diagnostics$rhat_A)) {
    rhat_parts$A <- fit$diagnostics$rhat_A |>
      dplyr::mutate(block = "loading A")
  }
  
  if (!is.null(fit$diagnostics$rhat_tau)) {
    rhat_parts$tau <- fit$diagnostics$rhat_tau |>
      dplyr::mutate(block = "cutpoints")
  }
  
  if (!is.null(fit$diagnostics$rhat_U)) {
    rhat_parts$U <- fit$diagnostics$rhat_U |>
      dplyr::mutate(block = "latent U")
  }
  
  rhat_df <- dplyr::bind_rows(rhat_parts)
  
  if (!"block" %in% names(rhat_df)) {
    rhat_df$block <- classify_parameter_block(rhat_df$parameter)
  }
  
  rhat_summary <- summarise_rhat_block(rhat_df)
  
  ess_df <- fit$diagnostics$ess_key
  
  if (!is.null(ess_df) && nrow(ess_df) > 0L) {
    ess_df$block <- classify_parameter_block(ess_df$parameter)
    ess_summary <- summarise_ess_block(ess_df)
  } else {
    ess_summary <- data.frame()
  }
  
  chain_stats <- fit$mcmc$chain_stats
  
  if (!is.null(chain_stats) && nrow(chain_stats) > 0L) {
    chain_stats <- chain_stats |>
      dplyr::mutate(
        u_ess_accept_rate = u_ess_accept_total / pmax(u_ess_total, 1),
        global_u_accept_rate = global_u_accept_total / pmax(global_u_total, 1),
        u_eval_per_update = u_ess_eval_total / pmax(u_ess_total, 1),
        global_u_eval_per_update = global_u_eval_total / pmax(global_u_total, 1),
        theta_eval_per_update = theta_eval_total / pmax(theta_update_total, 1)
      )
  }
  
  if (!is.null(out_dir)) {
    dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
    
    if (save_csvs) {
      utils::write.csv(
        rhat_df,
        file.path(out_dir, paste0(label, "_rhat_all.csv")),
        row.names = FALSE
      )
      
      utils::write.csv(
        rhat_summary,
        file.path(out_dir, paste0(label, "_rhat_summary_by_block.csv")),
        row.names = FALSE
      )
      
      if (!is.null(ess_df) && nrow(ess_df) > 0L) {
        utils::write.csv(
          ess_df,
          file.path(out_dir, paste0(label, "_ess_all.csv")),
          row.names = FALSE
        )
        
        utils::write.csv(
          ess_summary,
          file.path(out_dir, paste0(label, "_ess_summary_by_block.csv")),
          row.names = FALSE
        )
      }
      
      if (!is.null(chain_stats) && nrow(chain_stats) > 0L) {
        utils::write.csv(
          chain_stats,
          file.path(out_dir, paste0(label, "_chain_update_stats.csv")),
          row.names = FALSE
        )
      }
    }
  }
  
  plots <- list()
  
  if (nrow(rhat_df) > 0L) {
    if (nrow(rhat_df) > 0L) {
      rhat_plot_df <- rhat_df |>
        dplyr::filter(is.finite(rhat))
      
      rhat_ref <- data.frame(
        xintercept = c(1.01, 1.05, 1.10),
        ref = factor(
          c("1.01", "1.05", "1.10"),
          levels = c("1.01", "1.05", "1.10")
        )
      )
      
      plots$rhat_hist <- ggplot2::ggplot(
        rhat_plot_df,
        ggplot2::aes(x = rhat)
      ) +
        ggplot2::geom_histogram(
          bins = 40,
          fill = "steelblue",
          color = "white",
          alpha = 0.85
        ) +
        ggplot2::geom_vline(
          data = rhat_ref,
          ggplot2::aes(
            xintercept = xintercept,
            color = ref,
            linetype = ref
          ),
          inherit.aes = FALSE,
          linewidth = 0.6
        ) +
        ggplot2::scale_color_manual(
          values = c(
            "1.01" = "grey30",
            "1.05" = "orange3",
            "1.10" = "red3"
          ),
          name = "Reference"
        ) +
        ggplot2::scale_linetype_manual(
          values = c(
            "1.01" = "dotted",
            "1.05" = "dashed",
            "1.10" = "solid"
          ),
          name = "Reference"
        ) +
        ggplot2::facet_wrap(~block, scales = "free_y") +
        ggplot2::labs(
          x = expression(hat(R)),
          y = "Count",
          title = paste0("MCMC review: R-hat distribution by block, ", label),
          subtitle = "Vertical lines: 1.01, 1.05, 1.10"
        ) +
        ggplot2::theme(legend.position = "bottom")
    }
  }
  
  if (!is.null(ess_df) && nrow(ess_df) > 0L) {
    plots$ess_box <- ggplot2::ggplot(
      ess_df |> dplyr::filter(is.finite(ess), ess > 0),
      ggplot2::aes(x = block, y = ess, fill = block)
    ) +
      ggplot2::geom_boxplot(alpha = 0.80, outlier.alpha = 0.45) +
      ggplot2::scale_y_log10() +
      ggplot2::labs(
        x = NULL,
        y = "ESS, log scale",
        title = paste0("MCMC review: ESS distribution by block, ", label)
      ) +
      ggplot2::theme(legend.position = "none")
  }
  
  if (!is.null(chain_stats) && nrow(chain_stats) > 0L) {
    
    chain_stat_cols <- c(
      "u_ess_accept_rate",
      "global_u_accept_rate",
      "u_eval_per_update",
      "global_u_eval_per_update",
      "theta_eval_per_update"
    )
    
    chain_stat_cols <- intersect(chain_stat_cols, names(chain_stats))
    
    chain_long <- chain_stats |>
      dplyr::select(
        dplyr::all_of(c("chain", chain_stat_cols))
      ) |>
      tidyr::pivot_longer(
        cols = -chain,
        names_to = "statistic",
        values_to = "value"
      )
    
    plots$chain_update_stats <- ggplot2::ggplot(
      chain_long, ggplot2::aes(x = factor(chain), y = value)
    ) +
      ggplot2::geom_col(fill = "steelblue", alpha = 0.85) +
      ggplot2::facet_wrap(~statistic, scales = "free_y", ncol = 3) +
      ggplot2::labs(
        x = "Chain",
        y = NULL,
        title = paste0("MCMC review: chain update statistics, ", label)
      )
  }
  
  if (!is.null(out_dir) && save_plots) {
    if (!is.null(plots$rhat_hist)) {
      ggplot2::ggsave(
        file.path(out_dir, paste0(label, "_rhat_hist_by_block.pdf")),
        plots$rhat_hist,
        width = 10,
        height = 6
      )
    }
    
    if (!is.null(plots$ess_box)) {
      ggplot2::ggsave(
        file.path(out_dir, paste0(label, "_ess_box_by_block.pdf")),
        plots$ess_box,
        width = 8,
        height = 5
      )
    }
    
    if (!is.null(plots$chain_update_stats)) {
      ggplot2::ggsave(
        file.path(out_dir, paste0(label, "_chain_update_stats.pdf")),
        plots$chain_update_stats,
        width = 10,
        height = 6
      )
    }
  }
  
  list(
    rhat_df = rhat_df,
    rhat_summary = rhat_summary,
    ess_df = ess_df,
    ess_summary = ess_summary,
    chain_stats = chain_stats,
    plots = plots
  )
}


# Deterministic reference integration for the two-dimensional simulation DGP.
# This is a numerical oracle, not an exact latent sampler or fitted posterior.
oracle_m0_quadrature_2d <- function(X, C, true_params,
    orders = c(40L, 80L, 160L, 240L), tolerance = 1e-5) {
  X <- as.matrix(X); C <- as.matrix(C)
  if (nrow(X) != nrow(C) || nrow(X) < 1L || ncol(C) != true_params$q ||
      true_params$d != 2L || any(!is.finite(X)) || anyNA(C) ||
      any(C != floor(C) | C < 1 | C > true_params$m) ||
      !isTRUE(all.equal(unname(true_params$Omega), diag(true_params$q))) ||
      !true_params$score_error %in% c("gaussian", "logistic") ||
      !is.finite(true_params$lambda) || true_params$lambda <= 0) {
    stop("Quadrature requires valid inputs and the independent-score, two-dimensional Study II DGP.")
  }
  if (length(orders) < 2L || anyNA(orders) || any(orders < 2 | orders != floor(orders)) ||
      any(diff(orders) <= 0) || length(tolerance) != 1L || !is.finite(tolerance) || tolerance <= 0)
    stop("Invalid quadrature convergence controls.")
  keys <- pattern_key(C); previous <- NULL
  # Stable interval probability: subtract survival probabilities in the right tail.
  cdf <- if (true_params$score_error == "gaussian") stats::pnorm else stats::plogis
  scale <- true_params$lambda * if (true_params$score_error == "logistic") sqrt(3)/pi else 1
  for (order in orders) {
    jacobi <- matrix(0, order, order)
    jacobi[cbind(seq_len(order-1L), 2:order)] <- sqrt(seq_len(order-1L))
    eig <- eigen(jacobi + t(jacobi), symmetric = TRUE)
    ix <- expand.grid(seq_len(order), seq_len(order))
    U <- cbind(eig$values[ix[,1]], eig$values[ix[,2]])
    weights <- eig$vectors[1,]^2
    log_prior_weights <- log(weights[ix[,1]]) + log(weights[ix[,2]])
    location <- U %*% t(true_params$A)
    current <- probability <- numeric(nrow(X))
    for (key in unique(keys)) {
      rows <- which(keys == key); c <- C[rows[1],]; logw <- log_prior_weights
      for (j in seq_len(ncol(C))) {
        cuts <- c(-Inf, true_params$tau[j,], Inf)
        lower <- (cuts[c[j]] - location[,j])/scale
        upper <- (cuts[c[j]+1L] - location[,j])/scale
        prob <- ifelse(lower > 0, cdf(lower, lower.tail=FALSE)-cdf(upper, lower.tail=FALSE),
                       cdf(upper)-cdf(lower))
        logw <- logw + log(pmax(prob, 0))
      }
      peak <- max(logw)
      if (!is.finite(peak)) stop("Quadrature has no positive mass for pattern ", key)
      w <- exp(logw-peak); mass <- sum(w); w <- w/mass
      probability[rows] <- exp(peak)*mass
      for (i in rows) current[i] <- sum(w*f0_2d(X[rep(i,nrow(U)),,drop=FALSE],U,
                                               scenario=true_params$scenario))
    }
    if (!is.null(previous)) {
      difference <- abs(current-previous)
      # Check both the normalized expectation and the normalizing probability.
      probability_change <- abs(probability-previous_probability)/pmax(probability, .Machine$double.xmin)
      if (all(is.finite(current)) && max(difference) <= tolerance && max(probability_change) <= tolerance) {
        attr(current,"truth_diagnostics") <- data.frame(evaluation_row=seq_len(nrow(X)),
          pattern=keys, n_latent=NA_integer_, mcse=NA_real_, nested_half_difference=NA_real_,
          method="gauss_hermite", quadrature_order=order, quadrature_difference=difference,
          pattern_probability=probability, probability_relative_difference=probability_change)
        attr(current,"rejection_telemetry") <- data.frame()
        return(current)
      }
    }
    previous <- current; previous_probability <- probability
  }
  stop("Oracle quadrature did not meet the declared convergence tolerance; retain an unavailable truth score.")
}


# Failed reference evaluation must never select which datasets get fitted.
# Preserve the caller's RNG even when the reference task throws after drawing.
mixedgp_reference_task <- function(stage, computation) {
  had_seed <- exists(".Random.seed", envir=.GlobalEnv, inherits=FALSE)
  if (had_seed) saved_seed <- get(".Random.seed", envir=.GlobalEnv)
  on.exit(if (had_seed) assign(".Random.seed", saved_seed, envir=.GlobalEnv)
          else if (exists(".Random.seed", envir=.GlobalEnv, inherits=FALSE))
            rm(".Random.seed", envir=.GlobalEnv), add=TRUE)
  tryCatch(list(value=force(computation), status=data.frame(stage=stage,status="success",message="")),
    error=function(e) list(value=NULL,status=data.frame(stage=stage,status="failed",message=conditionMessage(e))))
}



# ---- core/00_sampler_v030.R ----
## eivGP 0.3.1: exact dense-GP, variance-collapsed, score-marginal sampler.
## This module is data-generic. Study designs and run budgets belong to callers.

mixedgp_v030_logadd <- function(a, b) {
  m <- pmax(a, b)
  out <- m + log1p(exp(-abs(a - b)))
  out[is.infinite(m) & m < 0] <- -Inf
  out
}


mixedgp_v030_priors <- function(priors = list(), p, d, m_vec) {
  defaults <- list(
    variance_shape = 3, variance_rate = 2, signal_shape = c(13, 3),
    log_theta_x_mean = rep(log(0.5), p), log_theta_x_sd = rep(1.5, p),
    u_dictionary = rep(list(exp(log(0.5) + 1.5 *
      stats::qnorm((seq_len(5L) - 0.5) / 5))), d),
    u_weights = NULL, category_alpha = lapply(m_vec, function(m) rep(2, m)),
    loading_sd = 2.5
  )
  if (!is.list(priors) || (length(priors) && (is.null(names(priors)) ||
      any(!nzchar(names(priors))) || anyDuplicated(names(priors))))) {
    stop("priors must be a uniquely named list.")
  }
  bad <- setdiff(names(priors), names(defaults))
  if (length(bad)) stop("Unknown prior fields: ", paste(bad, collapse = ", "))
  out <- defaults
  out[names(priors)] <- priors
  positive <- function(x, n, name) {
    if (!is.numeric(x) || length(x) != n || any(!is.finite(x)) || any(x <= 0))
      stop(name, " must contain ", n, " positive finite value(s).")
    as.numeric(x)
  }
  out$variance_shape <- positive(out$variance_shape, 1L, "variance_shape")
  out$variance_rate <- positive(out$variance_rate, 1L, "variance_rate")
  out$signal_shape <- positive(out$signal_shape, 2L, "signal_shape")
  out$loading_sd <- positive(out$loading_sd, 1L, "loading_sd")
  if (length(out$log_theta_x_mean) == 1L) out$log_theta_x_mean <- rep(out$log_theta_x_mean, p)
  if (length(out$log_theta_x_sd) == 1L) out$log_theta_x_sd <- rep(out$log_theta_x_sd, p)
  if (!is.numeric(out$log_theta_x_mean) || length(out$log_theta_x_mean) != p ||
      any(!is.finite(out$log_theta_x_mean))) stop("log_theta_x_mean must be finite and have length p or 1.")
  out$log_theta_x_sd <- positive(out$log_theta_x_sd, p, "log_theta_x_sd")
  grid <- out$u_dictionary
  if (is.numeric(grid) && is.null(dim(grid))) grid <- rep(list(grid), d)
  if (is.matrix(grid)) {
    if (!is.numeric(grid) || ncol(grid) != d || !nrow(grid) ||
        any(!is.finite(grid)) || any(grid <= 0)) stop("u_dictionary matrix must be positive, with d columns.")
    if (anyDuplicated(as.data.frame(grid))) stop("u_dictionary vectors must be unique.")
    out$dictionary_type <- "vector"
    w <- out$u_weights
    if (is.null(w)) w <- rep(1, nrow(grid))
    w <- positive(w, nrow(grid), "u_weights")
    out$u_weights <- w / sum(w)
  } else {
    if (!is.list(grid) || length(grid) != d) stop("u_dictionary must be a vector, a list of d grids, or an M by d matrix.")
    for (k in seq_len(d)) {
      grid[[k]] <- positive(grid[[k]], length(grid[[k]]), "u_dictionary")
      if (!length(grid[[k]]) || anyDuplicated(grid[[k]])) stop("Coordinate grids must be nonempty and unique.")
    }
    out$dictionary_type <- "product"
    w <- out$u_weights
    if (is.null(w)) w <- lapply(grid, function(g) rep(1, length(g)))
    if (d == 1L && is.numeric(w)) w <- list(w)
    if (!is.list(w) || length(w) != d) stop("u_weights must be a list with d entries.")
    out$u_weights <- lapply(seq_len(d), function(k) {
      a <- positive(w[[k]], length(grid[[k]]), "u_weights")
      a / sum(a)
    })
  }
  out$u_dictionary <- grid
  if (length(m_vec) == 1L && is.numeric(out$category_alpha))
    out$category_alpha <- list(out$category_alpha)
  if (!is.list(out$category_alpha) || length(out$category_alpha) != length(m_vec))
    stop("category_alpha must supply one positive vector per ordinal item.")
  out$category_alpha <- lapply(seq_along(m_vec), function(j)
    positive(out$category_alpha[[j]], m_vec[j], "category_alpha"))
  out
}


mixedgp_v030_control <- function(control = list(), n_mis) {
  out <- list(u_block_size = 8L, cross_every = 10L, ess_max_steps = 10000L,
              dictionary_mode = "conditional", max_dictionary_size = 125L,
              gp_block_schur = TRUE, threshold_update = "gibbs",
              sampler = "ESS", hmc = list(step=.05, leapfrog=c(8L,24L), target_accept=.8))
  if (!is.list(control) || (length(control) && (is.null(names(control)) ||
      any(!nzchar(names(control))) || anyDuplicated(names(control)))))
    stop("sampler_control must be a uniquely named list.")
  bad <- setdiff(names(control), c(names(out), "hybrid_every"))
  if (length(bad)) stop("Unknown sampler_control fields: ", paste(bad, collapse = ", "))
  out <- utils::modifyList(out, control)
  for (nm in c("u_block_size", "ess_max_steps", "max_dictionary_size"))
    out[[nm]] <- mixedgp_as_integer_strict(out[[nm]], nm, 1L, 1L)
  out$cross_every <- mixedgp_as_integer_strict(out$cross_every, "cross_every", 0L, 1L)
  out$dictionary_mode <- match.arg(out$dictionary_mode, c("conditional", "marginal"))
  out$threshold_update <- match.arg(out$threshold_update, c("gibbs", "ess"))
  out$sampler <- match.arg(out$sampler, c("ESS", "HMC", "HYBRID"))
  if (out$sampler == "HYBRID") {
    out$hybrid_every <- mixedgp_as_integer_strict(
      if (is.null(out$hybrid_every)) 20L else out$hybrid_every, "hybrid_every", 1L, 1L)
  } else if (!is.null(out$hybrid_every)) {
    stop("hybrid_every is only used with sampler = 'HYBRID'.")
  }
  if (!is.list(out$hmc) || !setequal(names(out$hmc), c("step","leapfrog","target_accept"))) stop("Invalid HMC controls.")
  if (length(out$hmc$step)!=1L || !is.finite(out$hmc$step) || out$hmc$step<=0 ||
      length(out$hmc$target_accept)!=1L || !is.finite(out$hmc$target_accept) || out$hmc$target_accept<=0 || out$hmc$target_accept>=1) stop("Invalid HMC step or target acceptance.")
  out$hmc$leapfrog <- mixedgp_as_integer_strict(out$hmc$leapfrog,"HMC leapfrog",1L,2L)
  if(out$hmc$leapfrog[1]>out$hmc$leapfrog[2])stop("HMC leapfrog bounds must be ordered.")
  if(out$sampler %in% c("HMC", "HYBRID") && out$dictionary_mode!="conditional")stop("HMC and HYBRID require conditional dictionary updates; use ESS for marginal mode.")
  out$gp_block_schur <- mixedgp_validate_flag(out$gp_block_schur, "gp_block_schur")
  out
}


mixedgp_v030_beta_forward <- function(e, shape) {
  lower <- e <= 0
  prob <- stats::pnorm(abs(e), lower.tail = FALSE, log.p = TRUE)
  ans <- numeric(length(e))
  ans[lower] <- stats::qbeta(prob[lower], shape[1L], shape[2L], log.p = TRUE)
  ans[!lower] <- stats::qbeta(prob[!lower], shape[1L], shape[2L],
                             lower.tail = FALSE, log.p = TRUE)
  if (any(!is.finite(ans)) || any(ans <= 0 | ans >= 1))
    stop("Beta transformation exceeded floating-point resolution; no boundary clipping is applied.")
  ans
}


mixedgp_v030_stick_forward <- function(e, alpha) {
  if (length(e) != length(alpha) - 1L || any(!is.finite(e)) ||
      any(!is.finite(alpha) | alpha <= 0)) stop("Invalid Gaussian stick coordinates.")
  v <- vapply(seq_along(e), function(k)
    mixedgp_v030_beta_forward(e[k], c(alpha[k], sum(alpha[-seq_len(k)]))), numeric(1))
  log_rem <- c(0, cumsum(log1p(-v)))
  pi <- exp(c(log(v) + log_rem[seq_along(v)], tail(log_rem, 1L)))
  if (any(!is.finite(pi)) || any(pi <= 0)) stop("Category probabilities underflowed.")
  pi / sum(pi)
}


mixedgp_v030_stick_inverse <- function(pi, alpha) {
  if (length(pi) != length(alpha) || any(!is.finite(pi) | pi <= 0) ||
      abs(sum(pi) - 1) > 1e-10 || any(!is.finite(alpha) | alpha <= 0))
    stop("Invalid category probabilities.")
  vapply(seq_len(length(pi) - 1L), function(k) {
    v <- pi[k] / sum(pi[k:length(pi)])
    a <- alpha[k]; b <- sum(alpha[-seq_len(k)])
    lp <- stats::pbeta(v, a, b, log.p = TRUE)
    if (lp < log(0.5)) stats::qnorm(lp, log.p = TRUE) else
      stats::qnorm(stats::pbeta(v, a, b, lower.tail = FALSE, log.p = TRUE),
                   lower.tail = FALSE, log.p = TRUE)
  }, numeric(1))
}


mixedgp_v030_cutpoints <- function(pi, scale = 1) {
  b <- vapply(seq_len(length(pi) - 1L), function(k) {
    lo <- sum(pi[seq_len(k)]); hi <- sum(pi[-seq_len(k)])
    if (lo <= hi) stats::qnorm(log(lo), log.p = TRUE) else
      stats::qnorm(log(hi), lower.tail = FALSE, log.p = TRUE)
  }, numeric(1))
  c(-Inf, scale * b, Inf)
}


## Inverse-CDF sampling from a truncated Beta, using both log tails. No
## rejection tuning or clipping: an unrepresentable interval is an error.
mixedgp_v031_truncated_beta <- function(lower, upper, shape, uniform = stats::runif(1L)) {
  if (length(lower) != 1L || length(upper) != 1L ||
      any(!is.finite(c(lower, upper, shape, uniform))) ||
      length(shape) != 2L || any(shape <= 0) || lower < 0 || upper > 1 ||
      lower >= upper || length(uniform) != 1L || uniform <= 0 || uniform >= 1)
    stop("Invalid truncated Beta interval, shapes, or uniform variate.")
  lp <- mixedgp_v030_logadd(log1p(-uniform) + stats::pbeta(lower, shape[1L], shape[2L], log.p = TRUE),
                            log(uniform) + stats::pbeta(upper, shape[1L], shape[2L], log.p = TRUE))
  lq <- mixedgp_v030_logadd(log1p(-uniform) + stats::pbeta(lower, shape[1L], shape[2L], lower.tail = FALSE, log.p = TRUE),
                            log(uniform) + stats::pbeta(upper, shape[1L], shape[2L], lower.tail = FALSE, log.p = TRUE))
  ans <- if (lp <= lq) stats::qbeta(lp, shape[1L], shape[2L], log.p = TRUE) else
    stats::qbeta(lq, shape[1L], shape[2L], lower.tail = FALSE, log.p = TRUE)
  if (!is.finite(ans) || ans <= lower || ans >= upper)
    stop("Truncated Beta draw exceeded floating-point resolution.")
  ans
}


## Conditional on all U, adjacent Dirichlet spacings have a Beta split.
## Observed categories restrict each split to the interval between the
## largest U with C <= k and smallest U with C > k. Use every subject,
## including calibration; empty categories are allowed. No GP likelihood.
mixedgp_v031_threshold_gibbs <- function(state, ctx) {
  if (ctx$measurement != "threshold") stop("Cutoff Gibbs requires threshold measurement.")
  pi <- mixedgp_v030_decode(state, ctx)$pi[[1L]]
  alpha <- ctx$priors$category_alpha[[1L]]
  u <- state$U[, 1L]; categories <- ctx$C[, 1L]
  for (k in seq_len(length(pi) - 1L)) {
    tau <- mixedgp_v030_cutpoints(pi)
    left <- tau[k]; right <- tau[k + 2L]
    low <- max(c(left, u[categories <= k]))
    high <- min(c(right, u[categories > k]))
    if (!is.finite(tau[k + 1L]) || low >= high)
      stop("Empty conditional cutoff support; check U and ordinal categories.")
    logmass <- log_normal_interval_prob(left, right)
    ## Work with the smaller side of the current split to avoid loss of
    ## precision when a cutpoint is deep in the upper normal tail.
    flip <- pi[k] > pi[k + 1L]
    fraction <- function(a, b) if (a == b) 0 else
      exp(log_normal_interval_prob(a, b) - logmass)
    lower <- if (flip) fraction(high, right) else fraction(left, low)
    upper <- if (flip) fraction(low, right) else fraction(left, high)
    split <- mixedgp_v031_truncated_beta(lower, upper,
      if (flip) alpha[c(k + 1L, k)] else alpha[c(k, k + 1L)])
    total <- pi[k] + pi[k + 1L]
    pi[c(k, k + 1L)] <- total * if (flip) c(1 - split, split) else c(split, 1 - split)
  }
  state$e_pi[[1L]] <- mixedgp_v030_stick_inverse(pi, alpha)
  if (!is.finite(mixedgp_v030_measurement_loglik(state, mixedgp_v030_decode(state, ctx), ctx)))
    stop("Cutoff Gibbs roundtrip exceeded floating-point resolution.")
  previous <- ctx$counts$threshold_gibbs_updates
  ctx$counts$threshold_gibbs_updates <- (if (is.null(previous)) 0 else previous) + length(pi) - 1L
  state
}


mixedgp_v030_interval_forward <- function(z, lower, upper) {
  n <- max(length(z), length(lower), length(upper))
  z <- rep_len(z, n); lower <- rep_len(lower, n); upper <- rep_len(upper, n)
  if (any(!is.finite(z)) || anyNA(lower) || anyNA(upper) || any(lower >= upper))
    stop("Invalid interval transform arguments.")
  mass <- log_normal_interval_prob(lower, upper)
  lp <- mixedgp_v030_logadd(stats::pnorm(lower, log.p = TRUE),
                           mass + stats::pnorm(z, log.p = TRUE))
  lq <- mixedgp_v030_logadd(stats::pnorm(upper, lower.tail = FALSE, log.p = TRUE),
                           mass + stats::pnorm(z, lower.tail = FALSE, log.p = TRUE))
  ans <- ifelse(lp <= lq, stats::qnorm(lp, log.p = TRUE),
                stats::qnorm(lq, lower.tail = FALSE, log.p = TRUE))
  if (any(!is.finite(ans)) || any(ans <= lower | ans >= upper))
    stop("Interval transformation exceeded floating-point resolution.")
  ans
}


mixedgp_v030_interval_inverse <- function(u, lower, upper) {
  n <- max(length(u), length(lower), length(upper))
  u <- rep_len(u, n); lower <- rep_len(lower, n); upper <- rep_len(upper, n)
  if (any(!is.finite(u)) || anyNA(lower) || anyNA(upper) ||
      any(u <= lower | u >= upper)) stop("Inputs must lie strictly inside their intervals.")
  mass <- log_normal_interval_prob(lower, upper)
  lp <- log_normal_interval_prob(lower, u) - mass
  lq <- log_normal_interval_prob(u, upper) - mass
  ans <- ifelse(lp <= lq, stats::qnorm(lp, log.p = TRUE),
                stats::qnorm(lq, lower.tail = FALSE, log.p = TRUE))
  if (any(!is.finite(ans))) stop("Inverse interval transformation exceeded floating-point resolution.")
  ans
}


mixedgp_v030_collapsed <- function(y, R, r, a = 3, b = 2) {
  n <- length(y)
  if (!is.matrix(R) || !identical(dim(R), c(n, n)) ||
      any(!is.finite(R)) || any(!is.finite(y)) ||
      length(r) != 1L || !is.finite(r) || r <= 0 || r >= 1 ||
      !is.finite(a) || !is.finite(b) || a <= 0 || b <= 0)
    stop("Invalid collapsed GP arguments.")
  B <- r * R
  diag(B) <- diag(B) + 1 - r
  ch <- tryCatch(chol(B), error = function(e) NULL)
  if (is.null(ch)) stop("Dense GP factorization failed; the covariance was not modified.")
  v <- forwardsolve(t(ch), y)
  Q <- sum(v * v)
  ld <- 2 * sum(log(diag(ch)))
  list(loglik = -0.5 * ld - (a + n / 2) * log(b + Q / 2),
       Q = Q, logdet = ld, chol = ch, B = B)
}


mixedgp_v030_ess <- function(v, loglik, max_steps = 10000L) {
  if (!length(v)) return(list(value = v, loglik = loglik(v), evaluations = 1L, angle = 0))
  if (any(!is.finite(v))) stop("ESS initial coordinates must be finite.")
  current <- loglik(v)
  if (length(current) != 1L || !is.finite(current))
    stop("ESS initial state has nonfinite log likelihood.")
  nu <- stats::rnorm(length(v))
  height <- current + log(stats::runif(1L))
  theta <- stats::runif(1L, 0, 2 * pi)
  lo <- theta - 2 * pi; hi <- theta
  for (step in seq_len(max_steps)) {
    proposal <- v * cos(theta) + nu * sin(theta)
    lp <- loglik(proposal)
    if (length(lp) != 1L || is.na(lp) || lp == Inf)
      stop("ESS likelihood evaluation is undefined.")
    if (lp >= height) return(list(value = proposal, loglik = lp,
                                  evaluations = step + 1L,
                                  angle = atan2(sin(theta), cos(theta))))
    if (theta < 0) lo <- theta else hi <- theta
    theta <- stats::runif(1L, lo, hi)
  }
  stop("ESS search budget exhausted; no successful transition was produced.")
}


mixedgp_v030_decode <- function(state, ctx) {
  prior <- ctx$priors
  pi <- lapply(seq_len(ctx$q), function(j)
    mixedgp_v030_stick_forward(state$e_pi[[j]], prior$category_alpha[[j]]))
  A <- matrix(0, ctx$q, ctx$d)
  if (ctx$measurement == "probit") {
    for (j in seq_len(ctx$q)) {
      active <- if (ctx$ident == "lower_triangular") seq_len(min(j, ctx$d)) else seq_len(ctx$d)
      A[j, active] <- prior$loading_sd * state$e_A[j, active]
      if (ctx$ident == "lower_triangular" && j <= ctx$d) {
        ## Half-normal quantile via log survival: P(|Z|>a)=2 P(Z>a).
        lq <- stats::pnorm(state$e_A[j, j], lower.tail = FALSE, log.p = TRUE) - log(2)
        A[j, j] <- prior$loading_sd * stats::qnorm(lq, lower.tail = FALSE, log.p = TRUE)
      }
    }
  }
  tau <- lapply(seq_len(ctx$q), function(j)
    mixedgp_v030_cutpoints(pi[[j]], if (ctx$measurement == "threshold") 1 else sqrt(1 + sum(A[j, ]^2))))
  theta_u <- if (prior$dictionary_type == "product")
    vapply(seq_len(ctx$d), function(k) prior$u_dictionary[[k]][state$J[k]], numeric(1)) else
    as.numeric(prior$u_dictionary[state$J, ])
  r <- mixedgp_v030_beta_forward(state$e_r, prior$signal_shape)
  list(pi = pi, tau = tau, A = A, r = r, theta_u = theta_u)
}


mixedgp_v030_corr <- function(state, dec, ctx) {
  weights <- c(exp(state$logtheta_x), dec$theta_u)
  if (any(!is.finite(weights))) stop("Kernel coefficients exceeded floating-point resolution.")
  Z <- cbind(ctx$X, state$U)
  D <- matrix(0, ctx$n, ctx$n)
  for (j in seq_along(weights)) D <- D + weights[j] * outer(Z[, j], Z[, j], "-")^2
  kernel_from_weighted_sqdist(D, ctx$kernel, ctx$matern_nu)
}


mixedgp_v030_measurement_loglik <- function(state, dec, ctx, rows = seq_len(ctx$n),
                                           items = seq_len(ctx$q)) {
  if (ctx$measurement == "threshold") {
    c <- ctx$C[rows, 1L]; u <- state$U[rows, 1L]; tau <- dec$tau[[1L]]
    return(if (all(u > tau[c] & u <= tau[c + 1L])) 0 else -Inf)
  }
  ans <- 0
  for (j in items) {
    mu <- as.numeric(state$U[rows, , drop = FALSE] %*% dec$A[j, ])
    c <- ctx$C[rows, j]
    ans <- ans + sum(log_normal_interval_prob(dec$tau[[j]][c] - mu,
                                              dec$tau[[j]][c + 1L] - mu))
  }
  ans
}


mixedgp_v030_dictionary <- function(priors, limit = Inf) {
  if (priors$dictionary_type == "vector") {
    if (nrow(priors$u_dictionary) > limit) stop("Dictionary exceeds max_dictionary_size for marginal updates.")
    return(list(indices = matrix(seq_len(nrow(priors$u_dictionary)), ncol = 1L),
                logweights = log(priors$u_weights)))
  }
  sizes <- lengths(priors$u_dictionary)
  if (sum(log(sizes)) > log(limit) + 1e-12)
    stop("Product dictionary exceeds max_dictionary_size; use conditional coordinate updates.")
  ind <- as.matrix(do.call(expand.grid, lapply(sizes, seq_len)))
  lw <- numeric(nrow(ind))
  for (k in seq_along(sizes)) lw <- lw + log(priors$u_weights[[k]][ind[, k]])
  list(indices = ind, logweights = lw)
}


mixedgp_v030_gp <- function(state, ctx, prepared = NULL) {
  dec <- mixedgp_v030_decode(state, ctx)
  key <- list(U = state$U, logtheta_x = state$logtheta_x, J = state$J, r = dec$r)
  if (!is.null(ctx$cache$key) && identical(ctx$cache$key, key)) return(ctx$cache$value)
  R <- mixedgp_v030_corr(state, dec, ctx)
  value <- NULL
  if (!is.null(prepared) && !identical(prepared$key, list(
      U = state$U[prepared$J, , drop = FALSE], logtheta_x = state$logtheta_x,
      J = state$J, r = dec$r))) prepared <- NULL
  if (!is.null(prepared)) {
    value <- tryCatch({
      B <- dec$r * R; diag(B) <- diag(B) + 1 - dec$r
      I <- prepared$I; J <- prepared$J
      cross <- B[I, J, drop = FALSE]
      z <- forwardsolve(t(prepared$chol), t(cross))
      S <- B[I, I, drop = FALSE] - crossprod(z)
      ch <- chol(S)
      e <- ctx$y[I] - as.numeric(cross %*% prepared$solve_y)
      v <- forwardsolve(t(ch), e)
      Q <- prepared$Q + sum(v * v)
      ld <- prepared$logdet + 2 * sum(log(diag(ch)))
      ll <- -0.5 * ld - (ctx$priors$variance_shape + ctx$n / 2) *
        log(ctx$priors$variance_rate + Q / 2)
      if (any(!is.finite(c(Q, ld, ll))) || Q < 0) stop("Invalid Schur calculation.")
      list(loglik = ll, Q = Q, logdet = ld, B = B, chol = NULL)
    }, error = function(e) NULL)
    if (is.null(value)) {
      ctx$counts$gp_block_fallbacks <- ctx$counts$gp_block_fallbacks + 1L
    } else ctx$counts$gp_block_evaluations <- ctx$counts$gp_block_evaluations + 1L
  }
  if (is.null(value)) {
    ctx$counts$gp_full_factorizations <- ctx$counts$gp_full_factorizations + 1L
    value <- mixedgp_v030_collapsed(ctx$y, R, dec$r,
                                    ctx$priors$variance_shape, ctx$priors$variance_rate)
  }
  ctx$cache$key <- key; ctx$cache$value <- value
  value
}


mixedgp_v030_prepare_block <- function(state, ctx, rows) {
  if (!ctx$control$gp_block_schur || ctx$control$dictionary_mode == "marginal" ||
      !length(rows) || length(rows) == ctx$n) return(NULL)
  base <- mixedgp_v030_gp(state, ctx)
  J <- setdiff(seq_len(ctx$n), rows)
  ctx$counts$gp_block_setups <- ctx$counts$gp_block_setups + 1L
  tryCatch({
    ch <- chol(base$B[J, J, drop = FALSE])
    v <- forwardsolve(t(ch), ctx$y[J])
    list(I = rows, J = J, chol = ch, solve_y = backsolve(ch, v),
         Q = sum(v * v), logdet = 2 * sum(log(diag(ch))),
         key = list(U = state$U[J, , drop = FALSE], logtheta_x = state$logtheta_x,
                    J = state$J, r = mixedgp_v030_decode(state, ctx)$r))
  }, error = function(e) NULL)
}


mixedgp_v030_loggp <- function(state, ctx, prepared = NULL) {
  if (ctx$control$dictionary_mode == "conditional")
    return(mixedgp_v030_gp(state, ctx, prepared)$loglik)
  dict <- ctx$dictionary
  lp <- vapply(seq_len(nrow(dict$indices)), function(h) {
    st <- state; st$J <- as.integer(dict$indices[h, ])
    mixedgp_v030_gp(st, ctx)$loglik + dict$logweights[h]
  }, numeric(1))
  max(lp) + log(sum(exp(lp - max(lp))))
}


mixedgp_v030_update_dictionary <- function(state, ctx) {
  prior <- ctx$priors
  draw_one <- function(indices, weights) {
    lp <- vapply(seq_len(nrow(indices)), function(h) {
      st <- state; st$J <- as.integer(indices[h, ])
      mixedgp_v030_gp(st, ctx)$loglik + log(weights[h])
    }, numeric(1))
    if (any(!is.finite(lp))) stop("Nonfinite dictionary likelihood.")
    ctx$counts$dictionary_evaluations <- ctx$counts$dictionary_evaluations + length(lp)
    probs <- exp(lp - max(lp))
    sample.int(length(lp), 1L, prob = probs)
  }
  before <- state$J
  if (prior$dictionary_type == "vector") {
    state$J <- as.integer(draw_one(matrix(seq_len(nrow(prior$u_dictionary)), ncol = 1L),
                                    prior$u_weights))
  } else {
    for (k in seq_len(ctx$d)) {
      sz <- length(prior$u_dictionary[[k]])
      choices <- matrix(rep(state$J, each = sz), nrow = sz)
      choices[, k] <- seq_len(sz)
      state$J[k] <- draw_one(choices, prior$u_weights[[k]])
    }
  }
  ctx$counts$dictionary_switches <- ctx$counts$dictionary_switches + sum(state$J != before)
  state
}


mixedgp_v030_initial_state <- function(ctx) {
  prior <- ctx$priors
  state <- list(U = matrix(stats::rnorm(ctx$n * ctx$d), ctx$n, ctx$d),
                e_pi = lapply(ctx$m_vec, function(m) stats::rnorm(m - 1L, sd = 0.5)),
                e_A = matrix(stats::rnorm(ctx$q * ctx$d, sd = 0.5), ctx$q, ctx$d),
                logtheta_x = stats::rnorm(ctx$p, prior$log_theta_x_mean,
                                          prior$log_theta_x_sd / 2),
                e_r = stats::rnorm(1L, sd = 0.5))
  state$J <- if (prior$dictionary_type == "product")
    vapply(seq_len(ctx$d), function(k) sample.int(length(prior$u_dictionary[[k]]),
                                                 1L, prob = prior$u_weights[[k]]), integer(1)) else
    sample.int(nrow(prior$u_dictionary), 1L, prob = prior$u_weights)
  if (length(ctx$observed)) state$U[ctx$observed, ] <- ctx$U_obs[ctx$observed, , drop = FALSE]
  if (ctx$measurement == "threshold") {
    ## Calibration-compatible increasing cutpoints, including unobserved levels.
    ## Work between consecutive distinct calibrated classes and interpolate all
    ## intervening boundaries; this avoids arbitrary raw-cutpoint bounds.
    obs <- ctx$observed
    if (length(obs)) {
      cls <- sort(unique(ctx$C[obs, 1L]))
      mins <- vapply(cls, function(c) min(ctx$U_obs[obs[ctx$C[obs, 1L] == c], 1L]), numeric(1))
      maxs <- vapply(cls, function(c) max(ctx$U_obs[obs[ctx$C[obs, 1L] == c], 1L]), numeric(1))
      if (length(cls) > 1L && any(head(maxs, -1L) >= tail(mins, -1L)))
        stop("Calibration values are incompatible with deterministic ordinal thresholds.")
      tau <- numeric(ctx$m_vec[1L] - 1L)
      left <- seq_len(cls[1L] - 1L)
      if (length(left)) tau[left] <- mins[1L] - rev(seq_along(left))
      if (length(cls) > 1L) for (j in seq_len(length(cls) - 1L)) {
        ids <- seq.int(cls[j], cls[j + 1L] - 1L)
        tau[ids] <- maxs[j] + (mins[j + 1L] - maxs[j]) *
          seq_along(ids) / (length(ids) + 1)
      }
      right <- if (tail(cls, 1L) < ctx$m_vec[1L])
        seq.int(tail(cls, 1L), ctx$m_vec[1L] - 1L) else integer(0)
      if (length(right)) tau[right] <- tail(maxs, 1L) + seq_along(right)
      edges <- c(-Inf, tau, Inf)
      pi <- exp(log_normal_interval_prob(head(edges, -1L), tail(edges, -1L)))
      if (any(pi <= 0) || any(!is.finite(pi))) stop("Calibrated cutpoints exceed numerical prior support.")
      pi <- pi / sum(pi)
      state$e_pi[[1L]] <- mixedgp_v030_stick_inverse(pi, prior$category_alpha[[1L]])
    }
    dec <- mixedgp_v030_decode(state, ctx)
    mis <- ctx$missing
    if (length(mis)) {
      c <- ctx$C[mis, 1L]
      state$U[mis, 1L] <- mixedgp_v030_interval_forward(stats::rnorm(length(mis)),
                                                       dec$tau[[1L]][c], dec$tau[[1L]][c + 1L])
    }
  }
  state
}


mixedgp_v030_record_ess <- function(ctx, result, block) {
  nm <- paste0(block, "_ess_total")
  en <- paste0(block, "_ess_eval_total")
  ctx$counts[[nm]] <- ctx$counts[[nm]] + 1L
  ctx$counts[[en]] <- ctx$counts[[en]] + result$evaluations
  ctx$counts$abs_angle_sum <- ctx$counts$abs_angle_sum + abs(result$angle)
  ctx$counts$angle_count <- ctx$counts$angle_count + 1L
}


mixedgp_v030_sweep_base <- function(state, ctx, iteration) {
  mis <- ctx$missing
  blocks <- if (length(mis)) split(mis[sample.int(length(mis))],
                ceiling(seq_along(mis) / ctx$control$u_block_size)) else list()
  if (isTRUE(ctx$control$sampler == "HMC")) state <- mixedgp_hmc_joint(state,ctx,iteration,ctx$burn)
  else for (rows in blocks) {
    dec <- mixedgp_v030_decode(state, ctx)
    prep <- mixedgp_v030_prepare_block(state, ctx, rows)
    if (ctx$measurement == "threshold") {
      c <- ctx$C[rows, 1L]
      lo <- dec$tau[[1L]][c]; hi <- dec$tau[[1L]][c + 1L]
      z <- mixedgp_v030_interval_inverse(state$U[rows, 1L], lo, hi)
      propose <- function(v) {
        st <- state
        st$U[rows, 1L] <- mixedgp_v030_interval_forward(v, lo, hi)
        st
      }
      residual <- function(v) mixedgp_v030_loggp(propose(v), ctx, prep)
    } else {
      z <- as.numeric(state$U[rows, , drop = FALSE])
      propose <- function(v) {
        st <- state; st$U[rows, ] <- matrix(v, length(rows), ctx$d); st
      }
      residual <- function(v) {
        st <- propose(v)
        mixedgp_v030_loggp(st, ctx, prep) +
          mixedgp_v030_measurement_loglik(st, dec, ctx, rows = rows)
      }
    }
    update <- mixedgp_v030_ess(z, residual, ctx$control$ess_max_steps)
    state <- propose(update$value)
    mixedgp_v030_record_ess(ctx, update, "u")
  }

  if (ctx$measurement == "threshold" && ctx$control$threshold_update == "gibbs") {
    state <- mixedgp_v031_threshold_gibbs(state, ctx)
  } else for (j in seq_len(ctx$q)) {
    active <- if (ctx$ident == "lower_triangular") seq_len(min(j, ctx$d)) else seq_len(ctx$d)
    z <- state$e_pi[[j]]
    if (ctx$measurement == "probit") z <- c(state$e_A[j, active], z)
    propose <- function(v) {
      st <- state
      if (ctx$measurement == "probit") {
        st$e_A[j, active] <- v[seq_along(active)]
        v <- v[-seq_along(active)]
      }
      st$e_pi[[j]] <- v
      st
    }
    residual <- function(v) {
      st <- propose(v)
      mixedgp_v030_measurement_loglik(st, mixedgp_v030_decode(st, ctx), ctx, items = j)
    }
    update <- mixedgp_v030_ess(z, residual, ctx$control$ess_max_steps)
    state <- propose(update$value)
    mixedgp_v030_record_ess(ctx, update, "measurement")
  }

  if (ctx$control$cross_every > 0L && iteration %% ctx$control$cross_every == 0L &&
      length(mis)) {
    if (ctx$measurement == "threshold") {
      dec <- mixedgp_v030_decode(state, ctx)
      c <- ctx$C[mis, 1L]
      zmis <- mixedgp_v030_interval_inverse(state$U[mis, 1L],
                                             dec$tau[[1L]][c], dec$tau[[1L]][c + 1L])
      propose <- function(v) {
        st <- state; st$e_pi[[1L]] <- v
        dd <- mixedgp_v030_decode(st, ctx)
        st$U[mis, 1L] <- mixedgp_v030_interval_forward(zmis, dd$tau[[1L]][c],
                                                      dd$tau[[1L]][c + 1L])
        st
      }
      residual <- function(v) {
        st <- propose(v); dd <- mixedgp_v030_decode(st, ctx)
        ll <- mixedgp_v030_measurement_loglik(st, dd, ctx, rows = ctx$observed)
        if (!is.finite(ll)) return(-Inf)
        mixedgp_v030_loggp(st, ctx) + sum(log(dd$pi[[1L]][c]))
      }
      z <- state$e_pi[[1L]]
    } else {
      ## Cycle measurement rows deterministically; the selected subject block is
      ## random independently of its values. Every coordinate of a subject moves.
      j <- 1L + ((iteration %/% ctx$control$cross_every - 1L) %% ctx$q)
      rows <- blocks[[1L]]
      prep <- mixedgp_v030_prepare_block(state, ctx, rows)
      active <- if (ctx$ident == "lower_triangular") seq_len(min(j, ctx$d)) else seq_len(ctx$d)
      na <- length(active); np <- length(state$e_pi[[j]])
      z <- c(state$e_A[j, active], state$e_pi[[j]], as.numeric(state$U[rows, , drop = FALSE]))
      propose <- function(v) {
        st <- state
        st$e_A[j, active] <- v[seq_len(na)]
        st$e_pi[[j]] <- v[na + seq_len(np)]
        st$U[rows, ] <- matrix(v[-seq_len(na + np)], length(rows), ctx$d)
        st
      }
      residual <- function(v) {
        st <- propose(v); dd <- mixedgp_v030_decode(st, ctx)
        ## Only item j changes for all subjects; only the selected subjects
        ## change for other items. Count this affected-factor union once.
        mixedgp_v030_loggp(st, ctx, prep) +
          mixedgp_v030_measurement_loglik(st, dd, ctx, items = j) +
          mixedgp_v030_measurement_loglik(st, dd, ctx, rows = rows,
                                         items = setdiff(seq_len(ctx$q), j))
      }
    }
    update <- mixedgp_v030_ess(z, residual, ctx$control$ess_max_steps)
    state <- propose(update$value)
    mixedgp_v030_record_ess(ctx, update, "cross")
  }

  if (!isTRUE(ctx$control$sampler == "HMC")) {
  prior <- ctx$priors
  z <- c((state$logtheta_x - prior$log_theta_x_mean) / prior$log_theta_x_sd, state$e_r)
  propose <- function(v) {
    st <- state
    st$logtheta_x <- prior$log_theta_x_mean + prior$log_theta_x_sd * v[seq_len(ctx$p)]
    st$e_r <- v[ctx$p + 1L]
    st
  }
  update <- mixedgp_v030_ess(z, function(v) mixedgp_v030_loggp(propose(v), ctx),
                             ctx$control$ess_max_steps)
  state <- propose(update$value)
  mixedgp_v030_record_ess(ctx, update, "theta")
  }
  ## In marginal mode ALL preceding blocks are J-marginal; restore J before
  ## conditional recovery/retention. Product conditional mode uses coordinate Gibbs.
  if (ctx$control$dictionary_mode == "marginal") {
    dict <- ctx$dictionary
    lp <- vapply(seq_len(nrow(dict$indices)), function(h) {
      st <- state; st$J <- as.integer(dict$indices[h, ])
      mixedgp_v030_gp(st, ctx)$loglik + dict$logweights[h]
    }, numeric(1))
    h <- sample.int(length(lp), 1L, prob = exp(lp - max(lp)))
    old <- state$J
    state$J <- as.integer(dict$indices[h, ])
    ctx$counts$dictionary_evaluations <- ctx$counts$dictionary_evaluations + length(lp)
    ctx$counts$dictionary_switches <- ctx$counts$dictionary_switches + sum(state$J != old)
  } else state <- mixedgp_v030_update_dictionary(state, ctx)
  state
}


# Fixed, state-independent composition; keep pure-kernel RNG streams unchanged.
mixedgp_v030_sweep <- function(state, ctx, iteration) {
  if (!identical(ctx$control$sampler, "HYBRID"))
    return(mixedgp_v030_sweep_base(state, ctx, iteration))
  on.exit(ctx$control$sampler <- "HYBRID")
  ctx$control$sampler <- "HMC"
  state <- mixedgp_v030_sweep_base(state, ctx, iteration)
  if (iteration %% ctx$control$hybrid_every == 0L) {
    ctx$control$sampler <- "ESS"
    state <- mixedgp_v030_sweep_base(state, ctx, iteration)
    state$hybrid_extra_sweeps <- if (is.null(state$hybrid_extra_sweeps)) 1L else
      state$hybrid_extra_sweeps + 1L
  }
  state
}


mixedgp_v030_fit <- function(X, y, C, U_obs, m_vec,
                             measurement = c("threshold", "probit"), ident = "none",
                             kernel = "se", matern_nu = 2.5,
                             n_iter = 1750L, burn = 500L, n_chains = 4L, seed = 1L,
                             parallel_chains = TRUE, n_cores = NULL,
                             priors = list(), sampler_control = list(),
                             store_scores = FALSE, .resume = NULL) {
  measurement <- match.arg(measurement)
  ident <- match.arg(ident, c("none", "lower_triangular"))
  spec <- normalize_gp_kernel(kernel, matern_nu)
  X <- as_numeric_matrix_strict(X, "X")
  U_obs <- as_numeric_matrix_strict(U_obs, "U_obs", nrow_expected = nrow(X), allow_na = TRUE)
  n <- nrow(X); d <- ncol(U_obs); p <- ncol(X)
  C <- as.matrix(C)
  if (!is.numeric(y) || length(y) != n || any(!is.finite(y)) ||
      !is.numeric(C) || nrow(C) != n || anyNA(C) ||
      length(m_vec) != ncol(C) || any(C != floor(C))) stop("Invalid sampler data.")
  for (j in seq_along(m_vec)) if (m_vec[j] < 2 || any(C[, j] < 1 | C[, j] > m_vec[j]))
    stop("Ordinal codes must lie within the declared categories.")
  count <- rowSums(is.finite(U_obs))
  if (any(count != 0 & count != d)) stop("Partially observed U rows are not supported.")
  observed <- which(count == d); missing_rows <- setdiff(seq_len(n), observed)
  if (measurement == "threshold" && (d != 1L || ncol(C) != 1L))
    stop("Deterministic threshold measurement requires one latent input and one ordinal item.")
  if (ident == "lower_triangular" && ncol(C) < d) stop("Triangular loadings require q >= d.")
  for (nm in c("n_iter", "n_chains", "seed", "burn")) {
    value <- mixedgp_as_integer_strict(get(nm), nm, if (nm %in% c("burn", "seed")) 0L else 1L, 1L)
    assign(nm, value)
  }
  if (n_iter <= burn) stop("n_iter must exceed burn.")
  priors <- mixedgp_v030_priors(priors, p, d, m_vec)
  control <- mixedgp_v030_control(sampler_control, length(missing_rows))
  if(control$sampler %in% c("HMC", "HYBRID") && !(spec$name=="se" || spec$matern_nu %in% c(1.5,2.5)))
    stop("HMC supports SE and Matern 3/2 or 5/2 kernels; select ESS for other smoothness values.")
  dictionary <- if (control$dictionary_mode == "marginal")
    mixedgp_v030_dictionary(priors, control$max_dictionary_size) else NULL
  if (!is.null(.resume) && !identical(mixedgp_v030_control(.resume$control, length(missing_rows)), control))
    stop("Continuation cannot change sampler controls, including the sampler, hybrid frequency or HMC tuning.")
  start <- if (is.null(.resume)) 0L else .resume$iteration
  if (!is.null(.resume) && (!identical(.resume$version, 2L) ||
      !identical(.resume$sampler_version, "0.3.1") || length(.resume$states) != n_chains ||
      .resume$burn != burn || start >= n_iter)) stop("Incompatible 0.3.1 continuation state.")
  seeds <- as.integer((as.double(seed) + 10000 * seq_len(n_chains)) %% .Machine$integer.max)
  backend <- mixedgp_parallel_backend(if (parallel_chains) n_cores else 1L)
  workers <- min(n_chains, backend$cores)
  use_fork <- parallel_chains && backend$backend == "fork" && workers > 1L
  if (!use_fork) workers <- 1L
  disk <- mixedgp_sampler_disk(list(X=X, y=y, C=C, U_obs=U_obs, m_vec=m_vec,
    measurement=measurement, ident=ident, kernel=spec, n_iter=n_iter, burn=burn,
    n_chains=n_chains, seed=seed, priors=priors, control=control,
    store_scores=store_scores, resume=.resume))
  if (!is.null(disk)) on.exit(unlink(file.path(disk, "LOCK"), recursive=TRUE), add=TRUE)
  checkpoint_enabled <- is.finite(getOption("eivGP.checkpoint_seconds",Inf)) ||
    is.finite(getOption("eivGP.checkpoint_iterations",Inf))
  run_chain <- function(cc) {
    ctx <- new.env(parent = environment())
    ctx$X <- X; ctx$y <- y; ctx$C <- C; ctx$U_obs <- U_obs
    ctx$n <- n; ctx$p <- p; ctx$q <- ncol(C); ctx$d <- d; ctx$m_vec <- m_vec
    ctx$observed <- observed; ctx$missing <- missing_rows
    ctx$measurement <- measurement; ctx$ident <- ident
    ctx$kernel <- spec$name; ctx$matern_nu <- spec$matern_nu
    ctx$priors <- priors; ctx$control <- control; ctx$dictionary <- dictionary; ctx$burn <- burn
    ctx$cache <- new.env(parent = emptyenv())
    ctx$counts <- new.env(parent = emptyenv())
    names <- c("u_ess_total", "u_ess_eval_total",
      "measurement_ess_total", "measurement_ess_eval_total", "cross_ess_total",
      "threshold_gibbs_updates",
      "cross_ess_eval_total", "theta_ess_total", "theta_ess_eval_total",
      "dictionary_evaluations", "dictionary_switches", "gp_full_factorizations",
      "gp_block_setups", "gp_block_evaluations", "gp_block_fallbacks",
      "abs_angle_sum", "angle_count")
    for (name in names) ctx$counts[[name]] <- 0
    if (is.null(.resume)) {
      set.seed(seeds[cc])
      state <- mixedgp_v030_initial_state(ctx)
    } else {
      cp <- .resume$states[[cc]]
      mixedgp_restore_rng_state(cp$rng)
      state <- cp$state
      if (!is.null(cp$cache)) {
        ctx$cache$key <- cp$cache$key
        ctx$cache$value <- cp$cache$value
      }
    }
    initial <- state
    dec <- mixedgp_v030_decode(state, ctx)
    if (!is.finite(mixedgp_v030_measurement_loglik(state, dec, ctx)))
      stop("Initial measurement likelihood is nonfinite.")
    if (is.null(.resume)) mixedgp_v030_gp(state, ctx)
    keep <- n_iter - max(burn, start)
    draws <- list(
      samples_U = array(NA_real_, c(keep, n, d)),
      samples_A = if (measurement == "probit") array(NA_real_, c(keep, ncol(C), d)) else NULL,
      samples_S = if (store_scores && measurement == "probit") array(NA_real_, c(keep, n, ncol(C))) else NULL,
      samples_tau = matrix(NA_real_, keep, sum(m_vec - 1L)),
      samples_logtheta = matrix(NA_real_, keep, 1L + p + d),
      samples_sigma2 = numeric(keep), samples_V = numeric(keep), samples_r = numeric(keep),
      samples_J = matrix(NA_integer_, keep, if (priors$dictionary_type == "product") d else 1L),
      samples_pi = matrix(NA_real_, keep, sum(m_vec))
    )
    colnames(draws$samples_tau) <- if (measurement == "threshold") paste0("tau", seq_len(m_vec - 1L)) else
      tau_names_from_mvec(m_vec)
    colnames(draws$samples_logtheta) <- c("log_rho", paste0("log_theta_x", seq_len(p)),
                                         paste0("log_theta_u", seq_len(d)))
    saved <- 0L
    chain_start <- start
    elapsed_before <- 0
    clock_start <- proc.time()[["elapsed"]]
    path <- if (!is.null(disk)) file.path(disk, sprintf("chain-%02d.rds", cc)) else NULL
    if (checkpoint_enabled && !is.null(path) && file.exists(path)) {
      cp <- readRDS(path)
      if (!identical(cp$schema, 1L) || cp$iteration < start || cp$iteration > n_iter)
        stop("Invalid sampler disk checkpoint: ", path)
      state <- cp$state; initial <- cp$initial
      draws <- mixedgp_checkpoint_restore_draws(cp$draws, draws); saved <- cp$saved
      chain_start <- cp$iteration; elapsed_before <- cp$elapsed_seconds
      list2env(cp$counts, ctx$counts)
      ctx$cache$key <- cp$cache$key; ctx$cache$value <- cp$cache$value
      mixedgp_restore_rng_state(cp$rng)
      message("Resuming chain ", cc, " at iteration ", chain_start, "/", n_iter)
    }
    last_save <- last_progress <- proc.time()[["elapsed"]]
    last_iteration <- chain_start
    checkpoint <- function(iteration, force=FALSE) {
      if (is.null(path)) return(invisible(NULL))
      now <- proc.time()[["elapsed"]]
      elapsed <- elapsed_before + now - clock_start
      saving <- mixedgp_checkpoint_due(now,last_save,iteration-last_iteration,force)
      if (saving) {
        mixedgp_checkpoint_write(list(schema=1L, iteration=iteration, state=state,
          initial=initial, draws=mixedgp_checkpoint_prefix(draws, saved), saved=saved, counts=as.list(ctx$counts),
          cache=list(key=ctx$cache$key,value=ctx$cache$value), rng=mixedgp_rng_state(),
          elapsed_seconds=elapsed), path)
        last_save <<- proc.time()[["elapsed"]]; last_iteration <<- iteration
      }
      if (force || saving || now-last_progress >= 60) {
        elapsed <- elapsed_before + proc.time()[["elapsed"]] - clock_start
        # Rate is local to this chain; queued chains and prediction are not included.
        remaining <- if (iteration > start) elapsed/(iteration-start)*(n_iter-iteration) else NA_real_
        mixedgp_timing_record(sub("[.]rds$", "-timing.csv", path), data.frame(
          updated_at=format(Sys.time(),tz="UTC",usetz=TRUE), chain=cc,
          iteration=iteration, target_iteration=n_iter, checkpoint_iteration=last_iteration,
          phase=if(iteration==n_iter) "complete" else if(iteration<=burn) "warmup" else "sampling",
          elapsed_seconds=elapsed, remaining_seconds=remaining,
          estimated_finish_utc=if(is.finite(remaining)) format(Sys.time()+remaining,tz="UTC",usetz=TRUE) else NA_character_))
        last_progress <<- proc.time()[["elapsed"]]
      }
      invisible(NULL)
    }
    if (chain_start < n_iter) for (iteration in seq.int(chain_start + 1L, n_iter)) {
      state <- mixedgp_v030_sweep(state, ctx, iteration)
      if (iteration <= burn) { checkpoint(iteration); next }
      saved <- saved + 1L
      dec <- mixedgp_v030_decode(state, ctx)
      if (length(observed) && !isTRUE(all.equal(unname(state$U[observed, , drop = FALSE]),
            unname(U_obs[observed, , drop = FALSE]), tolerance = 0)))
        stop("Calibrated inputs changed.")
      if (!is.finite(mixedgp_v030_measurement_loglik(state, dec, ctx)))
        stop("Retained state has invalid ordinal likelihood.")
      gp <- mixedgp_v030_gp(state, ctx)
      V <- 1 / stats::rgamma(1L, shape = priors$variance_shape + n / 2,
                             rate = priors$variance_rate + gp$Q / 2)
      if (!is.finite(V) || V <= 0) stop("Conditional variance recovery failed.")
      draws$samples_U[saved, , ] <- state$U
      if (measurement == "probit") draws$samples_A[saved, , ] <- dec$A
      if (store_scores && measurement == "probit") {
        S <- matrix(NA_real_, n, ncol(C))
        for (j in seq_len(ncol(C))) {
          mu <- as.numeric(state$U %*% dec$A[j, ])
          S[, j] <- mu + mixedgp_v030_interval_forward(stats::rnorm(n),
            dec$tau[[j]][C[, j]] - mu, dec$tau[[j]][C[, j] + 1L] - mu)
        }
        draws$samples_S[saved, , ] <- S
      }
      draws$samples_tau[saved, ] <- unlist(lapply(dec$tau, function(t) t[-c(1L, length(t))]), use.names = FALSE)
      draws$samples_logtheta[saved, ] <- c(0.5 * (log(dec$r) - log1p(-dec$r)),
                                           state$logtheta_x, log(dec$theta_u))
      draws$samples_sigma2[saved] <- V * (1 - dec$r)
      draws$samples_V[saved] <- V; draws$samples_r[saved] <- dec$r
      draws$samples_J[saved, ] <- state$J
      draws$samples_pi[saved, ] <- unlist(dec$pi, use.names = FALSE)
      checkpoint(iteration)
    }
    checkpoint(n_iter, force=TRUE)
    counters <- as.list(ctx$counts)
    counters$theta_eval_total <- counters$theta_ess_eval_total
    counters$mean_abs_angle <- counters$abs_angle_sum / max(1, counters$angle_count)
    stats <- as.data.frame(c(list(chain = cc, seed = seeds[cc], saved = saved), counters))
    c(draws, list(chain_id = cc, seed = seeds[cc], initial_state = initial,
      checkpoint = list(state = state, rng = mixedgp_rng_state(),
                        cache = list(key = ctx$cache$key, value = ctx$cache$value)), stats = stats))
  }
  timing <- system.time({
    chains <- mixedgp_parallel_lapply(as.list(seq_len(n_chains)), run_chain,
      n_cores = workers, seeds = seeds, mc.preschedule = FALSE)
  })
  list(chains = chains, priors = priors, control = control,
       iteration = n_iter, burn = burn, sampler_version = "0.3.1",
       time_seconds = unname(timing["elapsed"]),
       parallel_backend = if (use_fork) "fork" else "serial", parallel_cores = workers,
       restart_storage = if(checkpoint_enabled && !is.null(disk) &&
         length(list.files(disk,"^chain-[0-9]+[.]rds$"))) mixedgp_checkpoint_reference(disk) else NULL)
}



# ---- core/00_sampler_v030_api.R ----
## General-purpose adapters for the 0.3.1 collapsed sampler.
## Existing prediction helpers use an algebraically equivalent covariance:
## sigma2 = V * (1-r), rho^2 = r/(1-r).

mixedgp_v030_validate_arguments <- function(n_iter, burn, thin, n_chains,
                                           seed, parallel_chains, n_cores,
                                           verbose, preset, dots) {
  mixedgp_validate_named_dots(dots)
  if (length(dots)) {
    stop("Unsupported 0.3.1 sampler argument(s): ", paste(names(dots), collapse = ", "),
         ". Configure the new model with priors and sampler_control; historical sampler controls are incompatible.")
  }
  n_iter <- mixedgp_as_integer_strict(n_iter, "n_iter", 1L, 1L)
  burn <- mixedgp_as_integer_strict(burn, "burn", 0L, 1L)
  thin <- mixedgp_as_integer_strict(thin, "thin", 1L, 1L)
  n_chains <- mixedgp_as_integer_strict(n_chains, "n_chains", 1L, 1L)
  seed <- mixedgp_as_integer_strict(seed, "seed", 0L, 1L)
  if (n_iter <= burn) stop("n_iter must exceed burn.")
  if (thin != 1L) stop("Every post-warm-up draw is retained; thin must be 1.")
  parallel_chains <- mixedgp_validate_flag(parallel_chains, "parallel_chains")
  verbose <- mixedgp_validate_flag(verbose, "verbose")
  if (!is.null(n_cores)) n_cores <- mixedgp_as_integer_strict(n_cores, "n_cores", 1L, 1L)
  if (!is.null(preset)) {
    if (!is.character(preset) || length(preset) != 1L || is.na(preset) ||
        !preset %in% c("fast", "balanced", "robust", "thorough")) {
      stop("preset must be NULL or a recognized historical budget label.")
    }
    warning("preset is a deprecated budget label; 0.3.1 uses explicit n_iter, burn, priors, and sampler_control.",
            call. = FALSE)
  }
  list(n_iter = n_iter, burn = burn, n_chains = n_chains, seed = seed,
       parallel_chains = parallel_chains, n_cores = n_cores, verbose = verbose)
}


mixedgp_v030_standardize_xy <- function(X, y) {
  y <- as_numeric_matrix_strict(y, "y_raw")
  if (ncol(y) != 1L || nrow(y) < 2L) {
    stop("y_raw must contain at least two finite univariate observations.")
  }
  y <- as.numeric(y[, 1L])
  X <- as_numeric_matrix_strict(X, "X_raw", nrow_expected = length(y))
  if (is.null(colnames(X))) colnames(X) <- paste0("x", seq_len(ncol(X)))
  X_center <- colMeans(X)
  X_scale <- apply(X, 2L, stats::sd)
  y_center <- mean(y)
  y_scale <- stats::sd(y)
  if (any(!is.finite(X_scale) | X_scale <= 0)) {
    stop("Every numeric predictor must have positive finite variation.")
  }
  if (!is.finite(y_scale) || y_scale <= 0) stop("y_raw must have positive finite variation.")
  list(X_raw = X, X = sweep(sweep(X, 2L, X_center, "-"), 2L, X_scale, "/"),
       X_center = X_center, X_scale = X_scale, y_raw = y,
       y = (y - y_center) / y_scale, y_center = y_center, y_scale = y_scale)
}


mixedgp_v030_calibration <- function(U_obs, n, d, calib_idx = NULL,
                                     standardize_U = FALSE) {
  standardize_U <- mixedgp_validate_flag(standardize_U, "standardize_U")
  explicit_indices <- !is.null(calib_idx)
  if (explicit_indices) {
    calib_idx <- if (length(calib_idx)) mixedgp_as_integer_strict(calib_idx, "calib_idx", 1L) else integer(0)
    if (anyDuplicated(calib_idx) || any(calib_idx > n)) stop("Invalid calib_idx.")
    calib_idx <- sort(calib_idx)
  }
  full <- matrix(NA_real_, n, d, dimnames = list(NULL, paste0("u", seq_len(d))))
  if (!is.null(U_obs)) {
    U_obs <- as_numeric_matrix_strict(U_obs, "U_obs", ncol_expected = d, allow_na = TRUE)
    if (!is.null(colnames(U_obs))) colnames(full) <- colnames(U_obs)
    if (nrow(U_obs) == n) {
      count <- rowSums(is.finite(U_obs))
      if (any(count > 0L & count < d)) stop("Each U_obs row must contain every latent coordinate or only NA.")
      observed <- which(count == d)
      if (explicit_indices && !all(calib_idx %in% observed)) {
        stop("Every row selected by calib_idx must be completely observed in U_obs.")
      }
      if (!explicit_indices) calib_idx <- observed
      ## Low-level callers may supply a complete calibration matrix and select
      ## its observed rows explicitly. Values outside that selection are hidden.
      full[calib_idx, ] <- U_obs[calib_idx, , drop = FALSE]
    } else if (explicit_indices && nrow(U_obs) == length(calib_idx) && all(is.finite(U_obs))) {
      full[calib_idx, ] <- U_obs
    } else stop("U_obs must have n rows, or one complete row per calib_idx.")
  } else if (explicit_indices && length(calib_idx)) {
    stop("calib_idx requires observed latent values.")
  }
  if (is.null(calib_idx)) calib_idx <- integer(0)
  anchor <- mixedgp_latent_anchor_status(full, calib_idx, d)
  center <- rep(0, d)
  scale <- rep(1, d)
  if (standardize_U && length(calib_idx)) {
    center <- colMeans(full[calib_idx, , drop = FALSE])
    scale <- apply(full[calib_idx, , drop = FALSE], 2L, stats::sd)
    scale[!is.finite(scale) | scale <= 0] <- 1
  }
  model <- full
  if (length(calib_idx)) {
    model[calib_idx, ] <- sweep(sweep(full[calib_idx, , drop = FALSE],
                                    2L, center, "-"), 2L, scale, "/")
  }
  list(U_obs = model, U_obs_raw = full, U_center = center, U_scale = scale,
       U_names = colnames(full), standardize_U = standardize_U,
       calib_idx = calib_idx, miss_idx = setdiff(seq_len(n), calib_idx),
       latent_scale_anchored = anchor$anchored,
       latent_anchor_rank = anchor$affine_rank,
       latent_anchor_required_rank = anchor$required_rank)
}


mixedgp_v030_join_draws <- function(parts) {
  if (is.null(parts[[1L]])) return(NULL)
  if (length(dim(parts[[1L]])) == 3L) return(combine_chain_arrays(parts))
  if (is.matrix(parts[[1L]])) return(do.call(rbind, parts))
  unlist(parts, use.names = FALSE)
}


mixedgp_v030_append_chains <- function(chains, resume) {
  if (is.null(resume)) return(chains)
  if (length(resume$old_chains) != length(chains)) stop("Checkpoint chain count changed.")
  for (i in seq_along(chains)) {
    old <- resume$old_chains[[i]]
    for (name in names(old)[startsWith(names(old), "samples_")]) {
      if (is.null(old[[name]])) next
      if (is.null(chains[[i]][[name]])) stop("Checkpoint draw fields changed: ", name)
      chains[[i]][[name]] <- mixedgp_v030_join_draws(list(old[[name]], chains[[i]][[name]]))
    }
    if (!is.null(old$initial_state)) chains[[i]]$initial_state <- old$initial_state
    counters <- setdiff(names(old$stats), c("chain", "seed", "mean_abs_angle"))
    for (name in counters) chains[[i]]$stats[[name]] <-
      old$stats[[name]] + chains[[i]]$stats[[name]]
    if (all(c("abs_angle_sum", "angle_count") %in% names(chains[[i]]$stats))) {
      chains[[i]]$stats$mean_abs_angle <- chains[[i]]$stats$abs_angle_sum /
        max(1, chains[[i]]$stats$angle_count)
    } else if ("mean_abs_angle" %in% names(old$stats)) {
      before <- resume$iteration
      added <- resume$target_iteration - before
      chains[[i]]$stats$mean_abs_angle <- (old$stats$mean_abs_angle * before +
        chains[[i]]$stats$mean_abs_angle * added) / (before + added)
    }
  }
  chains
}


mixedgp_v030_raw_series <- function(fit) {
  engine <- mixedgp_fit_engine(fit)
  raw <- if (engine == "univariate") mixedgp_study1_raw_series(fit) else mixedgp_study2_raw_series(fit)
  ## These compatibility quantities are derived from V, r and dictionary states.
  raw <- raw[!grepl("^(rho$|sigma_epsilon$|theta_u)", names(raw))]
  chains <- fit$mcmc$samples_by_chain
  raw <- c(list(V = chains$V, r = chains$r), raw)
  dictionary <- fit$priors$u_dictionary
  for (k in seq_len(ncol(chains$J[[1L]]))) {
    size <- if (is.matrix(dictionary)) nrow(dictionary) else length(dictionary[[k]])
    if (size > 1L) raw[[paste0("J[", k, "]")]] <- lapply(chains$J, function(z) as.numeric(z[, k]))
    ## A fixed singleton is known from the model; any unvisited nontrivial
    ## dictionary state remains present and cannot be silently screened out.
    if (size > 1L) for (j in seq_len(size)) {
      raw[[paste0("dictionary_occupancy[", k, ",", j, "]")]] <-
        lapply(chains$J, function(z) as.numeric(z[, k] == j))
    }
  }
  if (length(chains$pi) && !is.null(chains$pi[[1L]])) {
    for (j in seq_len(ncol(chains$pi[[1L]]))) raw[[paste0("category_probability[", j, "]")]] <-
      lapply(chains$pi, function(z) z[, j])
  }
  raw
}


mixedgp_v030_checkpoint_hash <- function(object) {
  cp <- object$checkpoint
  cp$signature <- NULL
  path <- tempfile("eivgp-v030-checkpoint-", fileext = ".rds")
  on.exit(unlink(path), add = TRUE)
  saveRDS(list(data = object$data, control = object$control, priors = object$priors,
               kernel = object$kernel, mcmc = object$mcmc, checkpoint = cp,
               interface = object$interface, model_specification = object$model_specification),
          path, compress = FALSE, version = 2L)
  unname(tools::md5sum(path))
}


mixedgp_v030_fit_object <- function(result, data, engine, kernel, matern_nu,
                                    resume = NULL) {
  if (!identical(result$sampler_version, "0.3.1")) stop("Unexpected sampler version.")
  if (!is.null(resume)) resume$target_iteration <- result$iteration
  chains <- result$chains
  m_vec <- if (engine == "univariate") data$m else data$m_vec
  tau_names <- if (engine == "univariate") paste0("tau", seq_len(m_vec - 1L)) else
    unlist(lapply(seq_along(m_vec), function(j) paste0("tau[", j, ",", seq_len(m_vec[j] - 1L), "]")),
           use.names = FALSE)
  for (i in seq_along(chains)) colnames(chains[[i]]$samples_tau) <- tau_names
  if (engine == "univariate") for (i in seq_along(chains)) {
    dims <- dim(chains[[i]]$samples_U)
    chains[[i]]$samples_u <- matrix(chains[[i]]$samples_U, dims[1L], dims[2L])
  }
  chains <- mixedgp_v030_append_chains(chains, resume)
  fields <- c(if (engine == "univariate") "u" else c("U", "S", "A"),
              "tau", "logtheta", "sigma2", "V", "r", "J", "pi")
  by_chain <- setNames(lapply(fields, function(name) {
    lapply(chains, function(chain) chain[[paste0("samples_", name)]])
  }), fields)
  mcmc <- setNames(lapply(by_chain, mixedgp_v030_join_draws), paste0("samples_", fields))
  mcmc$samples_by_chain <- by_chain
  mcmc$chain_initial <- lapply(chains, function(chain) chain$initial_state)
  lengths <- vapply(by_chain$logtheta, nrow, integer(1L))
  mcmc$mcmc_draw_info <- data.frame(
    chain = rep(seq_along(chains), times = lengths),
    draw_within_chain = unlist(lapply(lengths, seq_len), use.names = FALSE))
  mcmc$chain_stats <- do.call(rbind, lapply(chains, function(chain) chain$stats))
  fit <- list(data = data, kernel = list(name = kernel, matern_nu = matern_nu),
              priors = result$priors, gp_prior = result$priors, control = result$control,
              sampler_strategy = switch(result$control$sampler, HMC = "collapsed_dictionary_hmc", HYBRID = "collapsed_dictionary_hybrid", "collapsed_dictionary_ess"),
              sampler_version = "0.3.1", mcmc = mcmc,
              checkpoint = list(version = 2L, sampler_version = "0.3.1",
                iteration = result$iteration, burn = result$burn, thin = 1L,
                control = result$control,
                states = lapply(chains, function(chain) chain$checkpoint)))
  raw <- mixedgp_v030_raw_series(fit)
  table <- suppressWarnings(mixedgp_summarize_diagnostic_series(raw))
  get_rows <- function(pattern) table[grepl(pattern, table$parameter), , drop = FALSE]
  hyper <- get_rows("^(V$|r$|theta_x|J\\[|dictionary_occupancy)")
  tau <- get_rows("^tau")
  latent <- get_rows("^[uU]\\[")
  latent$global_index <- as.integer(sub("^[uU]\\[([0-9]+).*$", "\\1", latent$parameter))
  if (engine == "multivariate") {
    latent$coord <- as.integer(sub("^U\\[[0-9]+,([0-9]+)\\]$", "\\1", latent$parameter))
  }
  finite_summary <- function(x, f) if (length(x) && all(is.finite(x))) f(x) else NA_real_
  diagnostic_summary <- data.frame(
    sampler_version = "0.3.1", sampler_strategy = switch(result$control$sampler, HMC = "collapsed_dictionary_hmc", HYBRID = "collapsed_dictionary_hybrid", "collapsed_dictionary_ess"),
    kernel = kernel, matern_nu = matern_nu, n_chains = length(chains),
    parallel_backend = result$parallel_backend, parallel_cores = result$parallel_cores,
    n_iter = result$iteration, burn = result$burn, thin = 1L,
    saved_per_chain = lengths[1L], total_saved_draws = sum(lengths),
    max_rhat_hyper = finite_summary(hyper$rhat, max),
    max_rhat_tau = finite_summary(tau$rhat, max),
    min_ess_key = finite_summary(table$ess, min),
    covariance_jitter = 0, forms_explicit_covariance_inverse = FALSE,
    time_seconds = as.numeric(result$time_seconds))
  diagnostic_summary[[if (engine == "univariate") "max_rhat_missing_u" else "max_rhat_missing_U"]] <-
    finite_summary(latent$rhat, max)
  diagnostic_summary[[if (engine == "univariate") "median_rhat_missing_u" else "median_rhat_missing_U"]] <-
    finite_summary(latent$rhat, stats::median)
  stats <- mcmc$chain_stats
  for (name in intersect(c("gp_full_factorizations", "gp_block_setups", "gp_block_evaluations",
                           "dictionary_evaluations", "dictionary_switches"), names(stats))) {
    diagnostic_summary[[name]] <- sum(stats[[name]])
  }
  fit$diagnostics <- list(rhat_hyper = hyper, rhat_tau = tau, ess_key = table,
                          summary = diagnostic_summary)
  fit$diagnostics[[if (engine == "univariate") "rhat_u" else "rhat_U"]] <- latent
  if (engine == "multivariate") {
    fit$diagnostics$rhat_A <- get_rows("^A\\[")
    fit$diagnostics$summary$max_rhat_A <- finite_summary(fit$diagnostics$rhat_A$rhat, max)
  }
  if (!is.null(result$restart_storage)) attr(fit,"eivgp_restart_storage") <- result$restart_storage
  mixedgp_hmc_attach_diagnostics(fit)
}


fit_eivgp_1d_v030 <- function(x_raw, y_raw, c_ord, u_true = NULL,
                              calib_idx = integer(0), m = 6L, tau_true = NULL,
                              n_iter = 1750L, burn = 500L, thin = 1L,
                              n_chains = 4L, preset = NULL, seed = 1L,
                              parallel_chains = TRUE, n_cores = NULL,
                              verbose = FALSE, kernel = c("se", "matern"),
                              matern_nu = 2.5, u_obs = NULL, priors = list(),
                              sampler_control = list(), .resume = NULL, ...) {
  args <- mixedgp_v030_validate_arguments(n_iter, burn, thin, n_chains, seed,
    parallel_chains, n_cores, verbose, preset, list(...))
  xy <- mixedgp_v030_standardize_xy(x_raw, y_raw)
  n <- length(xy$y)
  m <- mixedgp_as_integer_strict(m, "m", 2L, 1L)
  C <- prepare_ordinal_matrix(c_ord, m_vec = m, name = "c_ord")$C
  if (nrow(C) != n || ncol(C) != 1L) stop("c_ord must have one ordinal code per observation.")
  kernel <- match.arg(kernel)
  if (is.null(u_obs) && length(calib_idx)) {
    if (is.null(u_true) || length(u_true) != n) stop("Supply calibration values through u_obs.")
    u_obs <- rep(NA_real_, n)
    u_obs[calib_idx] <- u_true[calib_idx]
  }
  calibration <- mixedgp_v030_calibration(u_obs, n, 1L,
    if (length(calib_idx)) calib_idx else NULL)
  args$verbose <- NULL
  if (verbose) cat("Running collapsed dictionary EIV-GP (threshold measurement).\n")
  result <- do.call(mixedgp_v030_fit, c(list(X = xy$X, y = xy$y, C = C,
    U_obs = calibration$U_obs, m_vec = m, measurement = "threshold", ident = "none",
    kernel = kernel, matern_nu = matern_nu, priors = priors,
    sampler_control = sampler_control, store_scores = FALSE, .resume = .resume), args))
  data <- list(x_raw = xy$X_raw, x = xy$X, x_center = xy$X_center, x_scale = xy$X_scale,
    y_raw = xy$y_raw, y = xy$y, y_center = xy$y_center, y_scale = xy$y_scale,
    c_ord = as.integer(C[, 1L]), u_obs = as.numeric(calibration$U_obs),
    u_true = u_true, tau_true = tau_true, m = m, p_x = ncol(xy$X),
    predictor_names = colnames(xy$X), theta_spec = list(
      p_x = ncol(xy$X), x_index = 1L + seq_len(ncol(xy$X)),
      u_index = 2L + ncol(xy$X)))
  data <- c(data, calibration[setdiff(names(calibration), "U_obs")])
  mixedgp_v030_fit_object(result, data, "univariate", kernel, matern_nu, .resume)
}


fit_eivgp_ordprobit_v030 <- function(X_raw, y_raw, C_ord, U_obs = NULL,
                                     calib_idx = NULL, U_true_eval = NULL,
                                     d = 2L, m_vec = NULL,
                                     ident = c("lower_triangular", "none"),
                                     n_iter = 1750L, burn = 500L, thin = 1L,
                                     n_chains = 4L, preset = NULL, seed = 1L,
                                     parallel_chains = TRUE, n_cores = NULL,
                                     verbose = FALSE, kernel = c("se", "matern"),
                                     matern_nu = 2.5, standardize_U = FALSE,
                                     ordinal_levels = NULL, store_scores = FALSE,
                                     priors = list(), sampler_control = list(),
                                     .resume = NULL, ...) {
  args <- mixedgp_v030_validate_arguments(n_iter, burn, thin, n_chains, seed,
    parallel_chains, n_cores, verbose, preset, list(...))
  xy <- mixedgp_v030_standardize_xy(X_raw, y_raw)
  n <- length(xy$y)
  if (missing(d) && !is.null(U_obs)) d <- ncol(as.matrix(U_obs))
  d <- mixedgp_as_integer_strict(d, "d", 1L, 1L)
  ordinal <- prepare_ordinal_matrix(C_ord, m_vec = m_vec, level_maps = ordinal_levels, name = "C_ord")
  if (nrow(ordinal$C) != n) stop("C_ord must have one row per observation.")
  ident <- match.arg(ident)
  if (ident == "lower_triangular" && ncol(ordinal$C) < d) stop("lower_triangular requires q >= d.")
  kernel <- match.arg(kernel)
  store_scores <- mixedgp_validate_flag(store_scores, "store_scores")
  calibration <- mixedgp_v030_calibration(U_obs, n, d, calib_idx, standardize_U)
  args$verbose <- NULL
  if (verbose) cat("Running collapsed dictionary EIV-GP (ordinal-probit measurement).\n")
  result <- do.call(mixedgp_v030_fit, c(list(X = xy$X, y = xy$y, C = ordinal$C,
    U_obs = calibration$U_obs, m_vec = ordinal$m_vec, measurement = "probit", ident = ident,
    kernel = kernel, matern_nu = matern_nu, priors = priors,
    sampler_control = sampler_control, store_scores = store_scores, .resume = .resume), args))
  truth <- if (is.null(U_true_eval)) NULL else as_numeric_matrix_strict(
    U_true_eval, "U_true_eval", nrow_expected = n, ncol_expected = d)
  truth_model <- if (is.null(truth)) NULL else sweep(sweep(truth, 2L,
    calibration$U_center, "-"), 2L, calibration$U_scale, "/")
  data <- c(xy, calibration, list(C_ord = ordinal$C, C_level_maps = ordinal$level_maps,
    C_names = colnames(ordinal$C), U_true_eval = truth_model, U_true_eval_raw = truth,
    m_vec = ordinal$m_vec, q = ncol(ordinal$C), d = d, p = ncol(xy$X), ident = ident))
  mixedgp_v030_fit_object(result, data, "multivariate", kernel, matern_nu, .resume)
}



# ---- core/03_study2_published_competitors.R ----
############################################################
## 03_study2_published_competitors.R
##
## Common, auditable competitor layer for Studies I and II.
## The historical filename is retained so older entry points still work.
##
## Publication rules:
##   * use only the public implementation attached to the cited method;
##   * never relabel a home-built approximation as a published method;
##   * record unavailable and failed fits explicitly;
##   * return predictive draws for the noisy-response target Y* | X*, C*;
##   * retain the package's fitted latent/predictive mean separately, so that
##     recovery of m(x,c) is not scored by averaging Monte Carlo-noisy Y draws.
############################################################

`%||%` <- function(x, y) if (is.null(x)) y else x


mixedgp_competitor_registry <- function() {
  data.frame(
    method = c("UC-GP", "LVGP", "EzGP"),
    reference = c(
      "Qian et al. (2008); Zhou et al. (2011)",
      "Zhang et al. (2020)",
      "Xiao et al. (2021)"
    ),
    package = c("kergp", "LVGP", "EzGP"),
    implementation = c(
      "kergp::q1Symm",
      "LVGP::LVGP_fit and LVGP::LVGP_predict",
      "EzGP::EzGP_fit and EzGP::EzGP_predict"
    ),
    stringsAsFactors = FALSE
  )
}


## Backward-compatible name used by existing Study II scripts.
study2_competitor_registry <- mixedgp_competitor_registry


mixedgp_competitor_preflight <- function(
    methods = mixedgp_competitor_registry()$method,
    strict = FALSE) {
  registry <- mixedgp_competitor_registry()
  unknown <- setdiff(methods, registry$method)
  if (length(unknown) > 0L) {
    stop("Unknown published competitor(s): ", paste(unknown, collapse = ", "))
  }

  out <- registry[match(methods, registry$method), , drop = FALSE]
  out$available <- vapply(
    out$package,
    requireNamespace,
    logical(1),
    quietly = TRUE
  )
  out$version <- vapply(seq_len(nrow(out)), function(i) {
    if (!out$available[i]) return(NA_character_)
    as.character(utils::packageVersion(out$package[i]))
  }, character(1))

  if (isTRUE(strict) && any(!out$available)) {
    stop(
      "Publication preflight failed. Install the missing public package(s): ",
      paste(out$package[!out$available], collapse = ", "),
      ". No substitute implementation will be used."
    )
  }
  out
}


mixedgp_competitor_status <- function(method,
                                      status,
                                      optimization_status = "",
                                      optimization_attempts = "",
                                      message = "",
                                      warnings = "",
                                      elapsed = NA_real_,
                                      implementation = "",
                                      fitted_noise_variance = NA_real_) {
  data.frame(
    method = method,
    status = status,
    optimization_status = optimization_status,
    optimization_attempts = optimization_attempts,
    message = message,
    warnings = warnings,
    elapsed_seconds = elapsed,
    implementation = implementation,
    fitted_noise_variance = fitted_noise_variance,
    stringsAsFactors = FALSE
  )
}


study2_competitor_status <- mixedgp_competitor_status


format_mixedgp_optimizer_attempts <- function(attempt) {
  if (is.null(attempt) || nrow(attempt) == 0L) return("")
  if ("budget_name" %in% names(attempt)) return(paste(
    paste0("retry=",attempt$retry,",seed=",attempt$seed,",",attempt$budget_name,"=",attempt$budget,
      ",elapsed=",signif(attempt$elapsed_seconds,5),",status=",attempt$status,
      ",message=",attempt$message),collapse=" | "))

  fmt_num <- function(x, digits = 6L) {
    ifelse(
      is.finite(x),
      trimws(formatC(x, digits = digits, format = "fg")),
      "NA"
    )
  }
  paste(
    paste0(
      "retry=", attempt$retry,
      ",seed=", attempt$seed,
      ",limits=", attempt$max_iter_ini, "/", attempt$max_iter_lat,
      ",elapsed=", fmt_num(attempt$elapsed_seconds, 4L),
      ",status=", attempt$status,
      ",flag=", ifelse(is.na(attempt$selected_flag), "NA", attempt$selected_flag),
      ",objective=", fmt_num(attempt$reported_objective)
    ),
    collapse = " | "
  )
}


make_mixedgp_competitor_data <- function(X, C, m_vec, y = NULL) {
  checked <- validate_mixedgp_inputs(X, C, m_vec)
  X <- checked$X
  C <- checked$C
  m_vec <- checked$m_vec

  out <- as.data.frame(X)
  names(out) <- paste0("x", seq_len(ncol(X)))
  for (j in seq_len(ncol(C))) {
    out[[paste0("c", j)]] <- factor(C[, j], levels = seq_len(m_vec[j]))
  }
  if (!is.null(y)) {
    if (length(y) != nrow(out) || any(!is.finite(y))) {
      stop("y must be finite and have one value per training row.")
    }
    out$y <- as.numeric(y)
  }
  out
}


make_mixedgp_numeric_matrix <- function(X, C, m_vec) {
  checked <- validate_mixedgp_inputs(X, C, m_vec)
  out <- cbind(checked$X, checked$C)
  storage.mode(out) <- "double"
  colnames(out) <- c(
    paste0("x", seq_len(ncol(checked$X))),
    paste0("c", seq_len(ncol(checked$C)))
  )
  out
}


mixedgp_adapter_ucgp_once <- function(X_train,
                                 y_train,
                                 C_train,
                                 X_test,
                                 C_test,
                                 m_vec,
                                 n_draw,
                                 seed,
                                 n_starts = 8L) {
  if (!requireNamespace("kergp", quietly = TRUE)) {
    stop("Package 'kergp' is not installed.")
  }
  checked <- validate_mixedgp_inputs(X_train, C_train, m_vec)
  p <- ncol(checked$X)
  q <- ncol(checked$C)
  train <- make_mixedgp_competitor_data(
    checked$X, checked$C, checked$m_vec, y_train
  )
  test <- make_mixedgp_competitor_data(X_test, C_test, checked$m_vec)

  kernel_env <- new.env(parent = parent.frame())
  kernel_env$kx <- kergp::kGauss(d = p)
  if (!identical(kergp::inputNames(kernel_env$kx), paste0("x", seq_len(p)))) {
    stop("Unexpected quantitative-input names returned by kergp::kGauss().")
  }
  for (j in seq_len(q)) {
    kernel_env[[paste0("kc", j)]] <- kergp::q1Symm(
      factor = factor(seq_len(checked$m_vec[j]),
                      levels = seq_len(checked$m_vec[j])),
      input = paste0("c", j),
      cov = "corr"
    )
  }
  kernel_formula <- stats::as.formula(
    paste0(
      "~ ",
      paste(c("kx()", paste0("kc", seq_len(q), "()")), collapse = " * ")
    ),
    env = kernel_env
  )
  covariance <- kergp::covComp(formula = kernel_formula, where = kernel_env)

  set.seed(seed)
  multistart_fit <- kergp::gp(
    y ~ 1,
    data = train,
    cov = covariance,
    noise = TRUE,
    multistart = as.integer(n_starts),
    trace = 0
  )
  report <- multistart_fit$MLE$report
  finite_rows <- apply(report$par, 1L, function(z) all(is.finite(z)))
  eligible <- which(report$convergence & is.finite(report$logLik) & finite_rows)
  if (length(eligible) == 0L) {
    stop("No kergp likelihood start converged to a finite fit.")
  }
  selected <- eligible[which.max(report$logLik[eligible])]

  ## kergp's multistart routine selects the largest likelihood even when that
  ## start has a nonconvergence flag. Reconstruct the public-package GP object
  ## at the best *converged* reported solution instead of accepting that flag
  ## or discarding converged alternatives.
  selected_par <- as.numeric(report$par[selected, ])
  covariance_par <- methods::getMethod("coef", "covComp")(covariance)
  n_covariance_par <- length(covariance_par)
  covariance_selected <- kergp::`coef<-`(
    covariance,
    value = selected_par[seq_len(n_covariance_par)]
  )
  noise_var <- selected_par[n_covariance_par + 1L]
  fit <- kergp::gp(
    y ~ 1,
    data = train,
    cov = covariance_selected,
    estim = FALSE,
    varNoise = noise_var
  )
  fit$MLE <- multistart_fit$MLE
  fit$logLik <- as.numeric(report$logLik[selected])
  fit$selected_multistart_index <- selected
  fit$selected_multistart_report <- data.frame(
    start = seq_len(nrow(report$par)),
    logLik = as.numeric(report$logLik),
    convergence = as.logical(report$convergence)
  )
  pred <- stats::predict(
    fit,
    newdata = test,
    type = "UK",
    seCompute = TRUE,
    covCompute = FALSE,
    lightReturn = TRUE
  )
  noise_var <- as.numeric(fit$varNoise)
  if (length(noise_var) != 1L || !is.finite(noise_var) || noise_var < 0) {
    stop("kergp returned an invalid fitted noise variance.")
  }
  pred_var_y <- predictive_variance_diagonal(
    as.numeric(pred$sd)^2 + noise_var,
    nrow(test),
    context = "kergp"
  )
  latent_mean <- as.numeric(pred$mean)

  predictive_draws <- sample_independent_predictive_marginals(
    latent_mean, pred_var_y, n_draw, seed + 10000L
  )
  predictive_draws <- attach_predictive_normal_components(
    predictive_draws, latent_mean, pred_var_y
  )

  list(
    draws = predictive_draws,
    ## The fitted response nugget is zero mean, so the latent-function and
    ## noisy-response predictive means coincide.  Keep both names explicit:
    ## the former is used for m(x,c), whereas draws above target Y*.
    latent_mean = latent_mean,
    predictive_mean = latent_mean,
    predictive_variance = pred_var_y,
    fit = fit,
    fitted_noise_variance = noise_var,
    implementation = paste0(
      "kergp ", utils::packageVersion("kergp"),
      "::q1Symm; best of ", sum(report$convergence),
      " converged fits among ", n_starts, " likelihood starts"
    )
  )
}


# Run the author package's unchanged optimizer through our bounded executor.
# Only parallel dispatch bindings are scoped locally; no namespace is modified.
mixedgp_lvgp_fit <- function(..., parallel=FALSE) {
  if(!isTRUE(parallel) || is.null(mixedgp_budget()))
    return(LVGP::LVGP_fit(...,parallel=FALSE))
  fit <- LVGP::LVGP_fit
  dispatch <- new.env(parent=environment(fit))
  dispatch$detectCores <- function(...)mixedgp_budget_cores()
  dispatch$makeCluster <- function(spec,...)as.integer(spec)
  dispatch$stopCluster <- function(cl)invisible(NULL)
  dispatch$parLapply <- function(cl,X,fun,...)mixedgp_stage_lapply(X,fun,...,workers=cl)
  environment(fit) <- dispatch
  fit(...,parallel=TRUE)
}

mixedgp_adapter_lvgp <- function(X_train,
                                 y_train,
                                 C_train,
                                 X_test,
                                 C_test,
                                 m_vec,
                                 n_draw,
                                 seed,
                                 n_starts = 8L,
                                 max_retries = 3L,
                                 dim_z = 2L,
                                 max_iter_ini = 100L,
                                 max_iter_lat = 20L,
                                 rescue_iter_ini = 300L,
                                 rescue_iter_lat = 100L,
                                 max_elapsed_seconds = Inf,
                                 parallel = FALSE) {
  if (!requireNamespace("LVGP", quietly = TRUE)) {
    stop("Package 'LVGP' is not installed.")
  }
  checked <- validate_mixedgp_inputs(X_train, C_train, m_vec)
  p <- ncol(checked$X)
  q <- ncol(checked$C)
  train_x <- make_mixedgp_numeric_matrix(
    checked$X, checked$C, checked$m_vec
  )
  test_x <- make_mixedgp_numeric_matrix(X_test, C_test, checked$m_vec)

  ## LVGP_fit() selects the smallest objective without requiring its optim
  ## convergence code to be zero. Because the returned fit object is tied to
  ## that selected solution, do not reconstruct an undocumented alternative by
  ## hand. Use a prespecified sequence of complete author-package multistart
  ## fits, prefer a converged return, and—only if every rescue still reports a
  ## nonzero code—retain the best finite package-selected fit with an explicit
  ## optimizer warning. This follows the public package more faithfully than
  ## discarding every fit that the authors' prediction routine accepts.
  n_starts <- as.integer(n_starts)
  max_retries <- as.integer(max_retries)
  if (n_starts < 1L || max_retries < 1L) {
    stop("n_starts and max_retries must be positive integers.")
  }
  if (length(max_elapsed_seconds) != 1L ||
      is.na(max_elapsed_seconds) || max_elapsed_seconds <= 0) {
    stop("max_elapsed_seconds must be a positive number or Inf.")
  }
  candidate_fits <- vector("list", max_retries)
  attempt <- data.frame(
    retry = seq_len(max_retries),
    seed = as.integer(seed) + 1000L * (seq_len(max_retries) - 1L),
    max_iter_ini = c(
      as.integer(max_iter_ini),
      rep(as.integer(rescue_iter_ini), max_retries - 1L)
    ),
    max_iter_lat = c(
      as.integer(max_iter_lat),
      rep(as.integer(rescue_iter_lat), max_retries - 1L)
    ),
    elapsed_seconds = NA_real_, selected_flag = NA_integer_, converged = FALSE,
    status = "not_run",
    reported_objective = NA_real_, selection_objective = Inf, message = "",
    stringsAsFactors = FALSE
  )
  for (ss in seq_len(max_retries)) {
    attempt_start <- proc.time()[["elapsed"]]
    candidate <- tryCatch(
      local({
        if (is.finite(max_elapsed_seconds)) {
          setTimeLimit(
            elapsed = as.numeric(max_elapsed_seconds), transient = TRUE
          )
          on.exit(
            setTimeLimit(cpu = Inf, elapsed = Inf, transient = FALSE),
            add = TRUE
          )
        }
        mixedgp_lvgp_fit(
          X = train_x,
          Y = as.numeric(y_train),
          ind_qual = p + seq_len(q),
          dim_z = as.integer(dim_z),
          n_opt = n_starts,
          max_iter_ini = attempt$max_iter_ini[ss],
          max_iter_lat = attempt$max_iter_lat[ss],
          seed = as.integer(attempt$seed[ss]),
          progress = FALSE,
          parallel = isTRUE(parallel),
          noise = TRUE
        )
      }),
      error = function(e) e
    )
    attempt$elapsed_seconds[ss] <-
      proc.time()[["elapsed"]] - attempt_start
    if (inherits(candidate, "error")) {
      attempt$message[ss] <- conditionMessage(candidate)
      attempt$status[ss] <- if (grepl(
        "elapsed time limit|reached elapsed", attempt$message[ss],
        ignore.case = TRUE
      )) "timeout" else "error"
      next
    }
    best_try <- which.min(unlist(candidate$optim_hist$obj_sol_best))
    best_flag <- as.numeric(
      candidate$optim_hist$flag_sol_unique[[best_try]][1L]
    )
    finite_package_fit <- is.finite(candidate$fit_detail$min_n_log_l) &&
      is.finite(candidate$fit_detail$raw_min_eig)
    valid <- identical(best_flag, 0) && finite_package_fit
    attempt$selected_flag[ss] <- best_flag
    attempt$converged[ss] <- valid
    attempt$status[ss] <- if (valid) "converged" else "nonconverged"
    attempt$reported_objective[ss] <- candidate$fit_detail$min_n_log_l
    attempt$selection_objective[ss] <- if (finite_package_fit) {
      candidate$fit_detail$min_n_log_l
    } else {
      Inf
    }
    if (!valid) {
      attempt$message[ss] <- paste0(
        "Package-selected optim convergence code ", best_flag,
        "; extended-limit rescue attempted when available."
      )
    }
    if (finite_package_fit) candidate_fits[[ss]] <- candidate
    if (valid) break
  }
  finite_attempt <- is.finite(attempt$selection_objective)
  if (!any(finite_attempt)) {
    attempted <- attempt[!is.na(attempt$selected_flag), , drop = FALSE]
    detail <- if (nrow(attempted) == 0L) {
      paste(unique(attempt$message[nzchar(attempt$message)]), collapse = " | ")
    } else {
      paste0(
        "seed=", attempted$seed,
        ", limits=", attempted$max_iter_ini, "/", attempted$max_iter_lat,
        ", flag=", attempted$selected_flag,
        collapse = "; "
      )
    }
    failure <- simpleError(paste0(
      "LVGP did not return a converged finite-likelihood fit after the ",
      "prespecified default and extended-limit attempts. ", detail
    ))
    failure$optimizer_attempts <- attempt
    stop(failure)
  }
  converged_attempt <- which(attempt$converged)
  best_retry <- if (length(converged_attempt) > 0L) {
    converged_attempt[which.min(attempt$selection_objective[converged_attempt])]
  } else {
    which.min(attempt$selection_objective)
  }
  fit <- candidate_fits[[best_retry]]
  accepted_package_nonconvergence <- !isTRUE(attempt$converged[best_retry])
  if (accepted_package_nonconvergence) {
    warning(
      "LVGP returned a finite author-package-selected solution after the ",
      "prespecified attempt sequence, but optim reported code ",
      attempt$selected_flag[best_retry],
      ". The package-selected fit is retained and flagged; no internal ",
      "replacement or hand-edited latent embedding is used.",
      call. = FALSE
    )
  }
  pred <- LVGP::LVGP_predict(test_x, fit, MSE_on = 1)
  latent_var <- predictive_variance_diagonal(
    pred$MSE, nrow(test_x), context = "LVGP"
  )

  ## LVGP represents the fitted nugget on the internally standardized response
  ## correlation scale. Convert it back to the original response scale before
  ## adding it to the latent-mean MSE returned by LVGP_predict().
  response_range <- as.numeric(fit$data$Y_max - fit$data$Y_min)
  noise_var <- as.numeric(
    fit$fit_detail$sigma2 * fit$fit_detail$nug_opt * response_range^2
  )
  if (length(noise_var) != 1L || !is.finite(noise_var) || noise_var < 0) {
    stop("LVGP returned an invalid fitted noise variance.")
  }
  pred_var_y <- predictive_variance_diagonal(
    latent_var + noise_var,
    nrow(test_x),
    context = "LVGP noisy prediction"
  )
  latent_mean <- as.numeric(pred$Y_hat)

  predictive_draws <- sample_independent_predictive_marginals(
    latent_mean, pred_var_y, n_draw, seed + 10000L
  )
  predictive_draws <- attach_predictive_normal_components(
    predictive_draws, latent_mean, pred_var_y
  )

  list(
    draws = predictive_draws,
    latent_mean = latent_mean,
    predictive_mean = latent_mean,
    predictive_variance = pred_var_y,
    fit = fit,
    optimizer_attempts = attempt,
    optimization_status = if (accepted_package_nonconvergence) {
      if (max_retries > 1L) {
        "package_selected_nonconverged_after_rescue"
      } else {
        "package_selected_nonconverged"
      }
    } else if (best_retry == 1L) {
      "default_converged"
    } else {
      "rescued_converged"
    },
    fitted_noise_variance = noise_var,
    implementation = paste0(
      "LVGP ", utils::packageVersion("LVGP"),
      "::LVGP_fit(noise=TRUE, dim_z=", dim_z,
      ", n_opt=", n_starts, "); author-package-selected finite fit after ",
      "a prespecified sequence of default and extended-limit attempts; ",
      "nonzero optim codes are retained only after rescue exhaustion and ",
      "are explicitly flagged"
    )
  )
}


fit_ezgp_at_tau <- function(train_x,
                            y_train,
                            p,
                            q,
                            m_vec,
                            tau,
                            maxeval) {
  EzGP::EzGP_fit(
    X = train_x,
    Y = as.numeric(y_train),
    p = p,
    q = q,
    m = m_vec,
    tau = tau,
    maxeval = as.integer(maxeval)
  )
}


select_ezgp_tau_cv <- function(train_x,
                              y_train,
                              p,
                              q,
                              m_vec,
                              seed,
                              tau_fractions,
                              cv_folds,
                              maxeval,
                              cv_score = c("nlpd", "mse")) {
  cv_score <- match.arg(cv_score)
  n <- nrow(train_x)
  cv_folds <- max(2L, min(as.integer(cv_folds), n))
  set.seed(seed)
  fold_id <- sample(rep(seq_len(cv_folds), length.out = n))
  response_var <- max(stats::var(y_train), .Machine$double.eps)
  tau_grid <- unique(as.numeric(tau_fractions) * response_var)
  cv_loss <- rep(Inf, length(tau_grid))

  for (tt in seq_along(tau_grid)) {
    fold_loss <- rep(NA_real_, cv_folds)
    for (fold in seq_len(cv_folds)) {
      train_id <- fold_id != fold
      validation_id <- !train_id
      ans <- tryCatch({
        fit <- fit_ezgp_at_tau(
          train_x[train_id, , drop = FALSE],
          y_train[train_id], p, q, m_vec, tau_grid[tt], maxeval
        )
        pred <- EzGP::EzGP_predict(
          train_x[validation_id, , drop = FALSE], fit,
          MSE_on = as.integer(cv_score == "nlpd")
        )
        error <- as.numeric(pred$Y_hat) - y_train[validation_id]
        if (cv_score == "mse") {
          mean(error^2)
        } else {
          latent_var <- predictive_variance_diagonal(
            pred$MSE, sum(validation_id), context = "EzGP cross-validation"
          )
          predictive_var <- predictive_variance_diagonal(
            latent_var + tau_grid[tt],
            sum(validation_id),
            context = "EzGP cross-validation noisy prediction"
          )
          if (any(predictive_var <= 0)) {
            stop("EzGP cross-validation returned zero predictive variance.")
          }
          mean(
            0.5 * log(2 * base::pi * predictive_var) +
              0.5 * error^2 / predictive_var
          )
        }
      }, error = function(e) NA_real_)
      fold_loss[fold] <- ans
    }
    if (all(is.finite(fold_loss))) cv_loss[tt] <- mean(fold_loss)
  }
  if (!any(is.finite(cv_loss))) {
    stop("EzGP training-only cross-validation failed for every nugget value.")
  }
  best <- which.min(cv_loss)
  list(
    tau = tau_grid[best], tau_grid = tau_grid, cv_loss = cv_loss,
    cv_score = cv_score
  )
}


mixedgp_adapter_ezgp_once <- function(X_train,
                                 y_train,
                                 C_train,
                                 X_test,
                                 C_test,
                                 m_vec,
                                 n_draw,
                                 seed,
                                 tau_fractions = c(1e-6, 0.0025, 0.01, 0.04, 0.16),
                                 cv_folds = 3L,
                                 maxeval = 100L,
                                 cv_score = c("nlpd", "mse")) {
  cv_score <- match.arg(cv_score)
  if (!requireNamespace("EzGP", quietly = TRUE)) {
    stop("Package 'EzGP' is not installed.")
  }
  checked <- validate_mixedgp_inputs(X_train, C_train, m_vec)
  p <- ncol(checked$X)
  q <- ncol(checked$C)
  train_x <- make_mixedgp_numeric_matrix(
    checked$X, checked$C, checked$m_vec
  )
  test_x <- make_mixedgp_numeric_matrix(X_test, C_test, checked$m_vec)

  tau_selection <- select_ezgp_tau_cv(
    train_x = train_x,
    y_train = as.numeric(y_train),
    p = p,
    q = q,
    m_vec = checked$m_vec,
    seed = seed,
    tau_fractions = tau_fractions,
    cv_folds = cv_folds,
    maxeval = maxeval,
    cv_score = cv_score
  )
  fit <- fit_ezgp_at_tau(
    train_x, y_train, p, q, checked$m_vec,
    tau_selection$tau, maxeval
  )
  if (length(fit$param) == 0L || any(!is.finite(fit$param))) {
    stop("EzGP returned non-finite fitted parameters.")
  }
  pred <- EzGP::EzGP_predict(test_x, fit, MSE_on = 1)
  latent_var <- predictive_variance_diagonal(
    pred$MSE, nrow(test_x), context = "EzGP"
  )
  pred_var_y <- predictive_variance_diagonal(
    latent_var + tau_selection$tau,
    nrow(test_x),
    context = "EzGP noisy prediction"
  )
  latent_mean <- as.numeric(pred$Y_hat)

  predictive_draws <- sample_independent_predictive_marginals(
    latent_mean, pred_var_y, n_draw, seed + 10000L
  )
  predictive_draws <- attach_predictive_normal_components(
    predictive_draws, latent_mean, pred_var_y
  )

  list(
    draws = predictive_draws,
    latent_mean = latent_mean,
    predictive_mean = latent_mean,
    predictive_variance = pred_var_y,
    fit = fit,
    fitted_noise_variance = tau_selection$tau,
    tuning = tau_selection,
    implementation = paste0(
      "EzGP ", utils::packageVersion("EzGP"),
      "::EzGP_fit; nugget selected by ", cv_folds,
      "-fold training-only ", toupper(cv_score),
      " CV over prespecified variance fractions"
    )
  )
}


run_published_mixedgp_competitors <- function(
    X_train,
    y_train,
    C_train,
    X_test,
    C_test,
    n_draw,
    seed,
    m_vec = NULL,
    methods = mixedgp_competitor_registry()$method,
    strict = FALSE,
    controls = list()) {
  adapters <- list(
    `UC-GP` = mixedgp_adapter_ucgp,
    LVGP = mixedgp_adapter_lvgp,
    EzGP = mixedgp_adapter_ezgp
  )
  unknown <- setdiff(methods, names(adapters))
  if (length(unknown) > 0L) {
    stop("Unknown published competitor(s): ", paste(unknown, collapse = ", "))
  }
  checked <- validate_mixedgp_inputs(X_train, C_train, m_vec)
  validate_mixedgp_inputs(X_test, C_test, checked$m_vec)

  if (length(methods) == 0L) {
    return(list(
      draws = list(), means = list(), latent_means = list(),
      predictive_means = list(), predictive_variances = list(), fits = list(),
      status = mixedgp_competitor_status(
        method = character(0), status = character(0),
        optimization_status = character(0),
        optimization_attempts = character(0), message = character(0),
        warnings = character(0), elapsed = numeric(0),
        implementation = character(0), fitted_noise_variance = numeric(0)
      )
    ))
  }

  draws <- list()
  latent_means <- list()
  predictive_means <- list()
  predictive_variances <- list()
  fits <- list()
  status <- list()
  for (j in seq_along(methods)) {
    method <- methods[j]
    start <- proc.time()[3]
    common_args <- list(
      X_train = checked$X,
      y_train = as.numeric(y_train),
      C_train = checked$C,
      X_test = as.matrix(X_test),
      C_test = as.matrix(C_test),
      m_vec = checked$m_vec,
      n_draw = as.integer(n_draw),
      seed = as.integer(seed + j)
    )
    method_args <- controls[[method]] %||% list()
    duplicated_args <- intersect(names(common_args), names(method_args))
    if (length(duplicated_args) > 0L) {
      stop(
        "Do not override common competitor arguments in controls: ",
        paste(duplicated_args, collapse = ", ")
      )
    }
    captured_warnings <- character(0)
    ans <- tryCatch(
      withCallingHandlers(
        do.call(adapters[[method]], c(common_args, method_args)),
        warning = function(w) {
          captured_warnings <<- c(captured_warnings, conditionMessage(w))
          invokeRestart("muffleWarning")
        }
      ),
      error = function(e) e
    )
    elapsed <- proc.time()[3] - start
    warning_text <- paste(
      unique(trimws(gsub("[[:space:]]+", " ", captured_warnings))),
      collapse = " | "
    )

    if (inherits(ans, "error")) {
      failure_attempts <- ans$optimizer_attempts %||% NULL
      failure_status <- if (!is.null(failure_attempts) &&
          any(failure_attempts$status == "timeout")) {
        "timeout"
      } else if (!is.null(failure_attempts) &&
                 any(failure_attempts$status == "nonconverged")) {
        "nonconverged"
      } else {
        "failed"
      }
      status[[method]] <- mixedgp_competitor_status(
        method = method,
        status = "unavailable_or_failed",
        optimization_status = failure_status,
        optimization_attempts = format_mixedgp_optimizer_attempts(
          failure_attempts
        ),
        message = conditionMessage(ans),
        warnings = warning_text,
        elapsed = elapsed
      )
      if (isTRUE(strict)) stop(ans)
      next
    }

    ans <- validate_mixedgp_competitor_result(
      ans, n_draw = n_draw, n_test = nrow(X_test)
    )
    draws[[method]] <- ans$draws
    latent_means[[method]] <- ans$latent_mean
    predictive_means[[method]] <- ans$predictive_mean
    predictive_variances[[method]] <- ans$predictive_variance %||% NA_real_
    fits[[method]] <- ans$fit %||% NULL
    status[[method]] <- mixedgp_competitor_status(
      method = method,
      status = "success",
      optimization_status = ans$optimization_status %||% "converged",
      optimization_attempts = format_mixedgp_optimizer_attempts(
        ans$optimizer_attempts %||% NULL
      ),
      warnings = warning_text,
      elapsed = elapsed,
      implementation = ans$implementation %||% "public package",
      fitted_noise_variance = ans$fitted_noise_variance %||% NA_real_
    )
  }

  list(
    draws = draws,
    ## `means` is a concise alias for callers interested in m(x,c).  The two
    ## explicit names make the estimand transparent and remain identical for
    ## the zero-mean response-noise models fitted by these packages.
    means = latent_means,
    latent_means = latent_means,
    predictive_means = predictive_means,
    predictive_variances = predictive_variances,
    fits = fits,
    status = do.call(rbind, status)
  )
}


run_study2_published_competitors <- run_published_mixedgp_competitors

run_study1_published_competitors <- run_published_mixedgp_competitors


## Backward-compatible adapter name used in older validation code.
# Prespecified rescue sequence: stop at the first valid fit. Test responses are
# never supplied to an adapter and cannot influence retry or model selection.
mixedgp_retry_adapter <- function(adapter, args, budget_name, budgets, change_seed=TRUE) {
  attempts <- list()
  for (i in seq_along(budgets)) {
    trial <- args; trial[[budget_name]] <- budgets[i]
    trial$seed <- as.integer(args$seed + if (change_seed) (i-1L)*1000L else 0L)
    start <- proc.time()[3]
    result <- tryCatch(do.call(adapter,trial),error=function(e)e)
    ok <- !inherits(result,"error")
    attempts[[i]] <- data.frame(retry=i, seed=trial$seed, budget_name=budget_name,
      budget=budgets[i], elapsed_seconds=unname(proc.time()[3]-start),
      status=if(ok)"converged" else "failed", message=if(ok)"" else conditionMessage(result))
    if (ok) {
      result$optimizer_attempts <- do.call(rbind,attempts)
      result$optimization_status <- if(i==1L)"default_converged" else "rescued_converged"
      result$fit$recovery_attempts <- result$optimizer_attempts
      return(result)
    }
  }
  failure <- simpleError(paste0("Prespecified optimizer rescue exhausted: ",conditionMessage(result)))
  failure$optimizer_attempts <- do.call(rbind,attempts)
  stop(failure)
}


mixedgp_adapter_ucgp <- function(X_train,y_train,C_train,X_test,C_test,m_vec,n_draw,seed,
    n_starts=8L, rescue_starts=c(32L,64L)) {
  budgets <- c(n_starts,rescue_starts)
  if(anyNA(budgets)||any(!is.finite(budgets))||any(budgets<1|budgets!=floor(budgets)))stop("Invalid UC-GP start budgets.")
  mixedgp_retry_adapter(mixedgp_adapter_ucgp_once,
    list(X_train=X_train,y_train=y_train,C_train=C_train,X_test=X_test,C_test=C_test,
      m_vec=m_vec,n_draw=n_draw,seed=seed),"n_starts",as.integer(budgets))
}


mixedgp_adapter_ezgp <- function(X_train,y_train,C_train,X_test,C_test,m_vec,n_draw,seed,
    tau_fractions=c(1e-6,.0025,.01,.04,.16), cv_folds=3L, maxeval=100L,
    cv_score=c("nlpd","mse"), rescue_maxeval=c(300L,1000L)) {
  budgets <- c(maxeval,rescue_maxeval)
  if(anyNA(budgets)||any(!is.finite(budgets))||any(budgets<1|budgets!=floor(budgets)))stop("Invalid EzGP evaluation budgets.")
  mixedgp_retry_adapter(mixedgp_adapter_ezgp_once,
    list(X_train=X_train,y_train=y_train,C_train=C_train,X_test=X_test,C_test=C_test,
      m_vec=m_vec,n_draw=n_draw,seed=seed,tau_fractions=tau_fractions,cv_folds=cv_folds,
      cv_score=match.arg(cv_score)),"maxeval",as.integer(budgets),change_seed=FALSE)
}


study2_adapter_ucgp_qian <- mixedgp_adapter_ucgp




# ---- core/00_public_api.R ----
############################################################
## Package-ready public fitting interface
############################################################

mixedgp_encode_ordinal <- function(C, m_vec = NULL, ordinal_levels = NULL) {
  C_df <- if (is.data.frame(C)) C else as.data.frame(C, check.names = FALSE)
  prepared <- prepare_ordinal_matrix(
    C_df,
    m_vec = m_vec,
    level_maps = ordinal_levels,
    name = "C"
  )
  names(C_df) <- prepared$column_names
  list(
    C = prepared$C,
    raw = C_df,
    level_maps = prepared$level_maps,
    m_vec = prepared$m_vec
  )
}


mixedgp_validate_named_dots <- function(dots, argument = "...") {
  if (!is.list(dots)) stop(argument, " must be supplied as a list internally.")
  if (length(dots) > 0L &&
      (is.null(names(dots)) || anyNA(names(dots)) ||
       any(!nzchar(names(dots))))) {
    stop("Every argument passed through ", argument, " must be named.")
  }
  if (length(dots) > 0L && anyDuplicated(names(dots))) {
    stop("Arguments passed through ", argument, " must have unique names.")
  }
  invisible(dots)
}


mixedgp_do_call <- function(FUN, common, dots) {
  mixedgp_validate_named_dots(dots)
  duplicate <- intersect(names(common), names(dots))
  if (length(duplicate) > 0L) {
    stop("Do not duplicate public arguments through ...: ",
         paste(duplicate, collapse = ", "))
  }
  if (".resume" %in% names(dots)) stop(".resume is an internal argument.")
  arguments <- c(common, dots)
  result <- do.call(FUN, arguments)
  if (is.list(result) && !is.null(result$checkpoint)) {
    result$checkpoint$arguments <- arguments
  }
  result
}


mixedgp_resolve_loading_rule <- function(ident, n_calibrated, anchored) {
  # Affine anchoring is a sufficient geometric check, not a necessity theorem
  # for this already restricted Gaussian/diagonal-ARD model. In particular,
  # sparse calibration must not silently impose scientific zeros in A.
  if (is.null(ident)) {
    if (n_calibrated > 0L && !isTRUE(anchored)) {
      stop(
        "Sparse/rank-deficient calibration requires an explicit ident choice: ",
        "use ident='none' for unrestricted loadings in the supplied coordinates, ",
        "or ident='lower_triangular' only when its loading restrictions are ",
        "scientifically intended. Full affine rank is a sufficient anchoring ",
        "check, not a necessary condition for fitting unrestricted loadings."
      )
    }
    ident <- if (isTRUE(anchored)) "none" else "lower_triangular"
  }
  ident <- match.arg(ident, c("lower_triangular", "none"))
  if (ident == "none" && n_calibrated == 0L) {
    stop(
      "The public interface currently requires some calibration for ",
      "ident='none'. This is a supported-workflow restriction, not an ",
      "identification theorem; use the explicit lower_triangular working ",
      "model without calibration."
    )
  }
  ident
}


mixedgp_resolve_latent_scale <- function(object,
                                         scale = c("auto", "raw", "model"),
                                         argument = "scale") {
  scale <- match.arg(scale)
  anchored <- isTRUE(object$data$latent_scale_anchored)
  if (scale == "raw" && !anchored) {
    stop(
      argument,
      "='raw' requires calibration data that anchor the latent coordinates; ",
      "use 'model' for the identified working coordinates."
    )
  }
  internal <- if (scale == "auto") {
    if (anchored) "raw" else "model"
  } else {
    scale
  }
  list(
    internal = internal,
    label = if (!anchored && internal == "model") "working" else internal,
    anchored = anchored
  )
}


#' Fit the EIV-GP for ordinal mixed inputs
#'
#' This is the stable package-facing interface for both the univariate and
#' multivariate ordinal-input samplers. Independent chains use forked R
#' processes on macOS/Linux and automatically fall back to serial execution.
#'
#' @param X Numeric training design matrix.
#' @param y Numeric response vector.
#' @param C Ordinal inputs as factors or positive integer codes.
#' @param U_obs Optional calibrated latent inputs; unavailable rows are `NA`.
#' @param latent_dim Latent dimension. Inferred from informative `U_obs` when
#'   possible. It must be supplied explicitly for the multivariate engine when
#'   no calibrated latent input is available.
#' @param engine The scientific measurement model: `"univariate"` uses
#'   deterministic thresholding of one latent input, whereas `"multivariate"`
#'   uses the noisy ordinal-probit model. It must be selected explicitly.
#' @param m_vec Number of levels for each ordinal input.
#' @param ordinal_levels Optional ordered level labels, one vector per ordinal
#'   input. Unordered factors require this argument.
#' @param kernel `"se"` or `"matern"`.
#' @param matern_nu Smoothness of the Matérn kernel.
#' @param ident Multivariate loading structure. The default is `"none"` for
#'   full-affine-rank calibration and `"lower_triangular"` without calibration.
#'   Sparse/rank-deficient calibration requires an explicit choice. Triangular
#'   loadings together with diagonal ARD restrict the joint statistical model;
#'   they are not merely a change of coordinates. Explicit `"none"` requires
#'   at least one calibrated row in the supported public workflow.
#' @param standardize_U Whether calibrated latent inputs are standardized before
#'   fitting. The default is `TRUE` when calibration is available.
#' @param store_scores Whether to retain ordinal-probit score draws.
#' @param parallel Whether to run independent MCMC chains in parallel.
#' @param n_cores Number of forked workers; `NULL` uses package settings.
#' @param priors Named list of 0.3.1 prior settings. These include
#'   `variance_shape`, `variance_rate`, `signal_shape`, `log_theta_x_mean`,
#'   `log_theta_x_sd`, `u_dictionary`, `u_weights`, `category_alpha`, and
#'   `loading_sd`. `u_dictionary` is a list of positive coordinate grids or a
#'   positive matrix with one complete latent-kernel parameter vector per row.
#' @param sampler Main continuous update: `"ESS"` (default) for blocked
#'   elliptical slice sampling, `"HMC"` for joint Hamiltonian Monte Carlo, or
#'   `"HYBRID"` for joint HMC plus a complete ESS sweep every 20 iterations.
#'   Set `sampler_control$hybrid_every` to a positive integer to change that
#'   fixed schedule; it cannot change during continuation.
#'   All retain measurement and dictionary updates. HMC and HYBRID support
#'   conditional dictionary updates with SE or Matern 3/2 and 5/2 kernels.
#' @param sampler_control Named list controlling the collapsed sampler:
#'   `u_block_size`, `cross_every`, `ess_max_steps`, `dictionary_mode`,
#'   `max_dictionary_size`, `gp_block_schur`, `threshold_update`, and `hmc`.
#'   The `hmc` list accepts `step` (0.05), `leapfrog` (c(8L,24L)), and
#'   `target_accept` (0.8). Step size and a diagonal metric adapt only during
#'   warmup; continuation preserves the frozen adaptation and RNG state. Deterministic-threshold models
#'   default to exact coordinate cutoff Gibbs updates. Use
#'   `threshold_update = "gibbs"` (default), or `"ess"` to retain the previous joint cutoff slice update for
#'   comparisons. This setting does not change the ordinal-probit updates.
#'   `gp_block_schur` also accelerates the ordinal-probit joint move exactly.
#'   `u_block_size` counts subjects (default 8), with every latent coordinate
#'   of a selected subject updated together. Smaller blocks require more
#'   setups per sweep; assess effective samples per second when tuning it.
#'   The default dictionary update is conditional;
#'   `dictionary_mode = "marginal"` integrates the finite dictionary in the
#'   continuous latent-input updates.
#' @param ess_min Optional strict minimum scientific-target bulk and tail ESS.
#'   NULL keeps the fixed iteration budget. A positive number enables adaptive
#'   saved-state extensions after the initial budget, without new warmup.
#' @param ess_panel Optional list with X, C, and optional U in original units,
#'   plus n_latent and seed. Supply a fixed scientific evaluation panel. If
#'   omitted, five evenly spaced training rows monitor observable targets.
#'   At least four chains are required. ESS is a precision rule, not proof of
#'   convergence; inspect the full-history scientific traces.
#' @param ... Additional sampler arguments such as `n_iter`, `burn`, and
#'   `n_chains`. Every post-warm-up draw is retained: `thin`, if supplied, must
#'   equal one. Defaults are 1750 transitions, 500 warm-up transitions, and
#'   four chains. Version 0.3.1 uses a collapsed total-variance update, a
#'   continuous signal fraction and X-kernel parameters, a finite U-kernel
#'   dictionary, and ESS or joint HMC updates. Historical sampler
#'   settings such as `noise_strategy`, `sampler_strategy`, and
#'   `control_overrides` are rejected.
#' @return An EIV-GP fitted object.
#' @details X and y are always centered and scaled using training means and
#'   standard deviations. Variance priors and stored `samples_V` refer to the
#'   standardized response: the signal covariance is \eqn{V r R}, and noise
#'   variance is \eqn{V(1-r)}. Multiply either variance by `fit$data$y_scale^2`
#'   to express it in response units. The defaults are independent
#'   \eqn{V\sim IG(3,2)} (density proportional to \eqn{v^{-4}\exp(-2/v)}) and
#'   \eqn{r\sim Beta(13,3)}; r is a pointwise variance fraction, not realized
#'   regression R-squared. Kernel coefficients refer to model-scaled inputs.
#'   In particular, `u_dictionary` contains positive inverse squared-distance
#'   coefficients, not lengths: the SE correlation is
#'   \eqn{\exp\{-\sum_k\theta_k(z_k-z'_k)^2\}}.
#'
#'   Both engines use a fixed standard-normal latent population law,
#'   independent of X, in model coordinates. Standardization uses calibrated
#'   training rows and is plug-in preprocessing: its estimated moments are not
#'   additional posterior parameters, and scaling does not remove correlations.
#'   The multivariate score residual covariance is fixed to the identity.
#'   These restrictions are recorded in `model_specification`. The affine-rank
#'   flag used to permit raw-scale output is a conservative software reporting
#'   policy; it neither proves statistical identification nor checks coverage.
#' @export
fit_eivgp <- function(X,
                      y,
                      C,
                      U_obs = NULL,
                      latent_dim = NULL,
                      engine = NULL,
                      m_vec = NULL,
                      ordinal_levels = NULL,
                      kernel = c("se", "matern"),
                      matern_nu = 2.5,
                      ident = NULL,
                      standardize_U = NULL,
                      store_scores = FALSE,
                      parallel = TRUE,
                      n_cores = NULL,
                      priors = list(),
                      sampler_control = list(),
                      sampler = c("ESS", "HMC", "HYBRID"),
                      ess_min = NULL, ess_panel = NULL,
                      ...) {
  mixedgp_validate_ess_request(ess_min, ess_panel)
  fit_call <- match.call()
  sampler <- match.arg(sampler)
  if (!is.null(sampler_control$sampler) && sampler_control$sampler != sampler) stop("Conflicting sampler options.")
  sampler_control$sampler <- sampler
  caller_rng <- mixedgp_rng_state()
  on.exit(mixedgp_restore_rng_state(caller_rng), add = TRUE)
  parallel <- mixedgp_validate_flag(parallel, "parallel")
  store_scores <- mixedgp_validate_flag(store_scores, "store_scores")
  standardize_U <- mixedgp_validate_flag(
    standardize_U, "standardize_U", allow_null = TRUE
  )
  if (is.null(engine)) {
    stop(
      "Select engine explicitly: 'univariate' for deterministic thresholding ",
      "or 'multivariate' for the noisy ordinal-probit measurement model."
    )
  }
  engine <- match.arg(engine, c("univariate", "multivariate"))
  kernel <- match.arg(kernel)
  X <- as_numeric_matrix_strict(X, "X")
  y_mat <- as_numeric_matrix_strict(y, "y")
  if (ncol(y_mat) != 1L) stop("y must be a univariate response.")
  y <- as.numeric(y_mat[, 1L])
  encoded <- mixedgp_encode_ordinal(C, m_vec, ordinal_levels)
  C_mat <- encoded$C
  if (nrow(X) != length(y) || nrow(C_mat) != length(y)) {
    stop("X, y, and C must have the same number of observations.")
  }
  m_vec <- encoded$m_vec
  U_input_names <- if (is.null(U_obs) || is.null(dim(U_obs))) {
    NULL
  } else {
    colnames(U_obs)
  }
  U_checked <- if (is.null(U_obs)) {
    NULL
  } else {
    as_numeric_matrix_strict(
      U_obs, "U_obs", nrow_expected = length(y), allow_na = TRUE
    )
  }

  if (engine == "multivariate" && is.null(latent_dim)) {
    if (is.null(U_checked)) {
      stop(
        "latent_dim must be supplied explicitly for the multivariate engine ",
        "when no calibrated U_obs is available."
      )
    }
    finite_probe <- is.finite(U_checked)
    if (length(finite_probe) == 0L || !any(finite_probe)) {
      stop(
        "latent_dim must be supplied explicitly for the multivariate engine ",
        "when U_obs contains no calibration values."
      )
    }
    latent_dim <- ncol(U_checked)
  }
  dots <- list(...)
  if(!is.null(ess_min) && !is.null(dots$n_chains) && dots$n_chains<4L)
    stop("Scientific ESS requires at least four chains.")
  mixedgp_validate_named_dots(dots)
  if ("thin" %in% names(dots) &&
      !identical(mixedgp_as_integer_strict(dots$thin, "thin", 1L, 1L), 1L)) {
    stop("The public workflow keeps every post-warm-up draw; thin must be 1.")
  }
  dots$thin <- 1L
  if ("seed" %in% names(dots)) {
    dots$seed <- mixedgp_as_integer_strict(
      dots$seed, "seed", min_value = 0L, length_expected = 1L
    )
  }

  if (engine == "univariate") {
    if (ncol(C_mat) != 1L || length(m_vec) != 1L) {
      stop("The univariate engine requires exactly one ordinal input.")
    }
    u_obs <- if (is.null(U_checked)) {
      rep(NA_real_, length(y))
    } else {
      if (ncol(U_checked) != 1L) {
        stop("For the univariate engine, U_obs must have one column.")
      }
      as.numeric(U_checked[, 1L])
    }
    if (length(u_obs) != length(y)) {
      stop("For the univariate engine, U_obs must have one value per row.")
    }
    if (!is.null(ident) || isTRUE(store_scores)) {
      stop("ident and store_scores apply only to the multivariate engine.")
    }
    calibrated <- which(is.finite(u_obs))
    if (is.null(standardize_U)) standardize_U <- length(calibrated) > 0L
    u_center <- 0
    u_scale <- 1
    if (isTRUE(standardize_U) && length(calibrated) > 0L) {
      u_center <- mean(u_obs[calibrated])
      u_scale <- stats::sd(u_obs[calibrated])
      if (!is.finite(u_scale) || u_scale <= 0) u_scale <- 1
    }
    u_obs_raw <- u_obs
    if (length(calibrated) > 0L) {
      u_obs[calibrated] <- (u_obs[calibrated] - u_center) / u_scale
    }
    fit <- mixedgp_do_call(
      fit_eivgp_1d_v030,
      list(
        x_raw = X, y_raw = y, c_ord = C_mat[, 1L],
        u_obs = u_obs, calib_idx = which(is.finite(u_obs)), m = m_vec,
        kernel = kernel, matern_nu = matern_nu,
        priors = priors, sampler_control = sampler_control,
        parallel_chains = isTRUE(parallel), n_cores = n_cores
      ),
      dots
    )
    fit$data$u_obs_raw <- u_obs_raw
    fit$data$U_center <- u_center
    fit$data$U_scale <- u_scale
    fit$data$standardize_U <- isTRUE(standardize_U)
    fit$data$U_names <- U_input_names
    if (is.null(fit$data$U_names) || length(fit$data$U_names) != 1L ||
        !nzchar(fit$data$U_names)) {
      fit$data$U_names <- "u1"
    }
  } else {
    latent_dim <- mixedgp_as_integer_strict(
      latent_dim, "latent_dim", min_value = 1L, length_expected = 1L
    )
    U_full <- if (is.null(U_checked)) {
      matrix(NA_real_, length(y), latent_dim)
    } else {
      U_checked
    }
    if (!identical(dim(U_full), c(length(y), latent_dim))) {
      stop("U_obs must be an n by latent_dim matrix.")
    }
    partly_observed <- apply(U_full, 1L, function(z) any(is.finite(z)) &&
      !all(is.finite(z)))
    if (any(partly_observed)) {
      stop("Each calibrated row of U_obs must observe every latent coordinate.")
    }
    calib_idx <- which(apply(U_full, 1L, function(z) all(is.finite(z))))
    anchor_status <- mixedgp_latent_anchor_status(
      U_full, calib_idx, d = latent_dim
    )
    ident <- mixedgp_resolve_loading_rule(
      ident, length(calib_idx), anchor_status$anchored
    )
    if (ident == "lower_triangular") {
      warning(
        "Triangular loadings and diagonal ARD jointly restrict the model. ",
        "The loading zeros/sign constraint are not a cost-free identification ",
        "convention, particularly for physically calibrated U.",
        call. = FALSE
      )
    }
    if (length(calib_idx) == 0L) {
      warning(
        "No calibration data anchor the latent coordinates; using model-defined ",
        "working coordinates, and raw/physical-scale latent outputs remain ",
        "unavailable.",
        call. = FALSE
      )
    }
    if (is.null(standardize_U)) standardize_U <- length(calib_idx) > 0L
    fit <- mixedgp_do_call(
      fit_eivgp_ordprobit_v030,
      list(
        X_raw = X, y_raw = y, C_ord = encoded$raw, U_obs = U_full,
        calib_idx = calib_idx,
        d = latent_dim, m_vec = m_vec,
        ordinal_levels = encoded$level_maps,
        kernel = kernel, matern_nu = matern_nu,
        ident = ident, standardize_U = isTRUE(standardize_U),
        store_scores = isTRUE(store_scores),
        priors = priors, sampler_control = sampler_control,
        parallel_chains = isTRUE(parallel), n_cores = n_cores
      ),
      dots
    )
  }
  fit$data$C_names <- colnames(C_mat)
  fit$data$C_level_maps <- encoded$level_maps
  fit$call <- fit_call
  fit$interface <- list(
    engine = engine,
    parallel = isTRUE(parallel),
    n_cores = fit$diagnostics$summary$parallel_cores,
    ordinal_level_maps = encoded$level_maps
  )
  fit$model_specification <- list(
    latent_population = "N_d(0, I_d) in model coordinates, independent of X",
    latent_population_moments = "fixed; calibrated standardization is plug-in",
    standardize_U = isTRUE(fit$data$standardize_U),
    latent_kernel_metric = "diagonal ARD",
    latent_kernel_prior = "finite dictionary",
    variance_parameterization = "V { r R + (1-r) I }; V total variance, r signal fraction",
    sampler_version = "0.3.1",
    sampler = if (fit$control$sampler == "HYBRID")
      paste0("collapsed total variance, finite kernel dictionary, joint HMC plus a complete ESS sweep every ", fit$control$hybrid_every, " iterations") else if (fit$control$sampler == "HMC")
      "collapsed total variance, finite kernel dictionary, joint HMC within Gibbs" else
      "collapsed total variance, finite kernel dictionary, blocked elliptical slice",
    loading_structure = if (engine == "multivariate") ident else NA_character_,
    ordinal_residual_covariance = if (engine == "multivariate") "I_q" else NA_character_,
    prospective_convention = "fixed_fit",
    raw_scale_policy = "conservative full-affine-rank calibration gate"
  )
  class(fit) <- unique(c("eivgp_fit", class(fit)))
  fit$checkpoint$signature <- mixedgp_v030_checkpoint_hash(fit)
  mixedgp_finish_ess(fit, ess_min, ess_panel, parallel, n_cores)
}


mixedgp_first_value <- function(x, default = NA) {
  if (is.null(x) || length(x) == 0L) return(default)
  x[[1L]]
}


mixedgp_diagnostic_value <- function(diagnostics, name, default = NA) {
  if (is.null(diagnostics) || !name %in% names(diagnostics) ||
      length(diagnostics[[name]]) == 0L) {
    return(default)
  }
  diagnostics[[name]][[1L]]
}


mixedgp_fit_engine <- function(object) {
  engine <- object$interface$engine
  if (is.null(engine)) {
    engine <- if (!is.null(object$data$d)) "multivariate" else "univariate"
  }
  engine <- as.character(engine[[1L]])
  if (!engine %in% c("univariate", "multivariate")) {
    stop("Malformed eivgp_fit object: unknown engine.")
  }
  engine
}


mixedgp_fit_dimensions <- function(object, engine) {
  data <- object$data
  X <- if (engine == "univariate") {
    if (!is.null(data$x_raw)) data$x_raw else data$x
  } else {
    if (!is.null(data$X_raw)) data$X_raw else data$X
  }
  n <- if (!is.null(X)) nrow(as.matrix(X)) else length(data$y)
  p <- if (!is.null(X)) ncol(as.matrix(X)) else if (engine == "univariate") {
    mixedgp_first_value(data$p_x, NA_integer_)
  } else {
    mixedgp_first_value(data$p, NA_integer_)
  }
  q <- if (engine == "univariate") {
    1L
  } else if (!is.null(data$q)) {
    data$q
  } else {
    ncol(as.matrix(data$C_ord))
  }
  d <- if (engine == "univariate") {
    1L
  } else if (!is.null(data$d)) {
    data$d
  } else {
    dim(object$mcmc$samples_U)[3L]
  }
  c(
    n = as.integer(mixedgp_first_value(n, NA_integer_)),
    p = as.integer(mixedgp_first_value(p, NA_integer_)),
    q = as.integer(mixedgp_first_value(q, NA_integer_)),
    d = as.integer(mixedgp_first_value(d, NA_integer_))
  )
}


#' Summarize a fitted EIV-GP model
#'
#' The summary reports the fitted measurement-model engine, problem dimensions,
#' calibration rank and latent-scale interpretation, kernel, MCMC configuration,
#' and the convergence summaries already computed by the sampler.
#'
#' @param object An object returned by [fit_eivgp()].
#' @param ... Reserved for future summary controls.
#' @return An object of class `summary_eivgp_fit`.
#' @export
summary.eivgp_fit <- function(object, ...) {
  if (!inherits(object, "eivgp_fit")) {
    stop("object must be returned by fit_eivgp().")
  }
  engine <- mixedgp_fit_engine(object)
  dimensions <- mixedgp_fit_dimensions(object, engine)
  calib_idx <- object$data$calib_idx
  if (is.null(calib_idx)) calib_idx <- integer(0)
  n_calibrated <- length(calib_idx)
  anchored <- isTRUE(object$data$latent_scale_anchored)
  required_rank <- mixedgp_first_value(
    object$data$latent_anchor_required_rank,
    dimensions[["d"]] + 1L
  )
  affine_rank <- mixedgp_first_value(
    object$data$latent_anchor_rank,
    if (n_calibrated == 0L) 0L else NA_integer_
  )
  calibration_status <- if (n_calibrated == 0L) {
    "none"
  } else if (anchored) {
    "affine_anchored"
  } else {
    "rank_deficient"
  }
  interpretation <- if (n_calibrated == 0L) {
    paste(
      "No calibration data were observed. The observed-input mean m(x,c) and",
      "prediction given (x,c) remain model-defined; latent coordinates,",
      "f(x,u), and latent-input imputations are defined on the model's",
      "working scale, so raw/physical-scale latent claims are unavailable."
    )
  } else if (!anchored) {
    paste0(
      "Calibration data inform U but have affine rank ", affine_rank, " of ",
      required_rank, ". The conservative software reporting gate retains ",
      "model-scale outputs. This rank check is sufficient for affine anchoring, ",
      "not a necessary-condition theorem for statistical identification."
    )
  } else {
    paste0(
      "The calibration design has full affine rank ", affine_rank, " of ",
      required_rank, "; the software permits raw-coordinate output. ",
      "This geometric check does not establish likelihood identification, ",
      "representative calibration, or response-surface coverage."
    )
  }

  diagnostics <- object$diagnostics$summary
  if (is.null(diagnostics)) diagnostics <- data.frame()
  if (is.data.frame(diagnostics) && nrow(diagnostics) > 1L) {
    diagnostics <- diagnostics[1L, , drop = FALSE]
  }
  convergence_names <- names(diagnostics)[
    grepl("rhat|ess", names(diagnostics), ignore.case = TRUE) &
      !grepl("accept|eval|total", names(diagnostics), ignore.case = TRUE)
  ]
  convergence <- lapply(
    convergence_names,
    function(name) mixedgp_diagnostic_value(diagnostics, name)
  )
  names(convergence) <- convergence_names
  mcmc_names <- c(
    "sampler_strategy", "n_chains", "parallel_backend", "parallel_cores",
    "n_iter", "burn", "thin", "saved_per_chain", "total_saved_draws",
    "time_seconds"
  )
  mcmc <- lapply(
    mcmc_names,
    function(name) mixedgp_diagnostic_value(diagnostics, name)
  )
  names(mcmc) <- mcmc_names
  mcmc$sampler_version <- mixedgp_first_value(object$sampler_version, NA_character_)
  mcmc$noise_strategy <- if (isTRUE(object$sampler_version %in% c("0.3.0", "0.3.1")))
    "collapsed total variance" else mixedgp_first_value(object$control$noise_strategy, NA_character_)
  mcmc$dictionary_mode <- mixedgp_first_value(object$control$dictionary_mode, NA_character_)
  mcmc$u_block_size <- mixedgp_first_value(object$control$u_block_size, NA_integer_)

  structure(
    list(
      call = object$call,
      engine = engine,
      dimensions = dimensions,
      calibration = list(
        n = as.integer(n_calibrated),
        status = calibration_status,
        affine_rank = as.integer(affine_rank),
        required_rank = as.integer(required_rank),
        anchored = anchored,
        output_scale = if (anchored) "raw" else "working",
        interpretation = interpretation
      ),
      identification = mixedgp_first_value(object$data$ident, NA_character_),
      model_specification = object$model_specification,
      kernel = list(
        name = mixedgp_first_value(object$kernel$name, "unknown"),
        matern_nu = mixedgp_first_value(object$kernel$matern_nu, NA_real_),
        label = mixedgp_kernel_label(object$kernel)
      ),
      mcmc = mcmc,
      convergence = convergence,
      diagnostic_tables = object$diagnostics[
        setdiff(names(object$diagnostics), "summary")
      ]
    ),
    class = c("summary_eivgp_fit", "list")
  )
}


#' Print a fitted EIV-GP model
#'
#' @param x An object returned by [fit_eivgp()].
#' @param ... Reserved for future print controls.
#' @return `x`, invisibly.
#' @export
print.eivgp_fit <- function(x, ...) {
  fit_summary <- summary.eivgp_fit(x)
  dims <- fit_summary$dimensions
  calibration <- fit_summary$calibration
  cat("<eivgp_fit>\n")
  cat(" Engine: ", fit_summary$engine, "\n", sep = "")
  cat(
    " Data: ", dims[["n"]], " rows; ", dims[["p"]], " numeric, ",
    dims[["q"]], " ordinal, and ", dims[["d"]], " latent dimension(s)\n",
    sep = ""
  )
  cat(
    " Calibration: ", calibration$n, " rows; affine rank ",
    calibration$affine_rank, "/", calibration$required_rank, "; ",
    calibration$output_scale, " latent scale\n",
    sep = ""
  )
  cat(" Kernel: ", fit_summary$kernel$label, "\n", sep = "")
  saved <- fit_summary$mcmc$total_saved_draws
  if (length(saved) == 1L && is.finite(saved)) {
    cat(" Posterior draws: ", format(saved, scientific = FALSE), "\n", sep = "")
  }
  invisible(x)
}


#' @rdname summary.eivgp_fit
#' @param x A fitted-model summary returned by `summary()`.
#' @return `x`, invisibly, for the print method.
#' @export
print.summary_eivgp_fit <- function(x, ...) {
  dims <- x$dimensions
  calibration <- x$calibration
  cat("Ordinal mixed-input EIV-GP fit\n")
  cat(" Engine: ", x$engine, "\n", sep = "")
  cat(
    " Training design: n = ", dims[["n"]], ", p = ", dims[["p"]],
    ", q = ", dims[["q"]], ", d = ", dims[["d"]], "\n",
    sep = ""
  )
  cat(
    " Calibration: ", calibration$n, " rows; affine rank ",
    calibration$affine_rank, "/", calibration$required_rank,
    "; status = ", calibration$status, "\n",
    sep = ""
  )
  if (!is.na(x$identification)) {
    cat(" Loading structure: ", x$identification, "\n", sep = "")
  }
  cat(" Kernel: ", x$kernel$label, "\n", sep = "")
  if (!is.null(x$model_specification)) {
    cat(" Latent population: ", x$model_specification$latent_population,
        "\n", sep = "")
    cat(" Calibration standardization: ",
        if (isTRUE(x$model_specification$standardize_U)) "plug-in" else "off",
        "; prospective inference: fixed fit\n", sep = "")
  }

  mcmc <- x$mcmc
  if (is.finite(mcmc$n_chains) && is.finite(mcmc$n_iter)) {
    cat(
      " MCMC: ", mcmc$n_chains, " chain(s) x ", mcmc$n_iter,
      " iterations; burn = ", mcmc$burn, ", thin = ", mcmc$thin,
      "; saved = ", mcmc$total_saved_draws, "\n",
      sep = ""
    )
  }
  convergence <- unlist(x$convergence, use.names = TRUE)
  rhat <- convergence[
    grepl("rhat", names(convergence), ignore.case = TRUE) &
      is.finite(convergence)
  ]
  ess <- convergence[
    grepl("ess", names(convergence), ignore.case = TRUE) &
      is.finite(convergence)
  ]
  if (length(rhat) > 0L || length(ess) > 0L) {
    pieces <- character(0)
    if (length(rhat) > 0L) {
      pieces <- c(pieces, paste0("largest reported R-hat = ", format(max(rhat))))
    }
    if (length(ess) > 0L) {
      pieces <- c(pieces, paste0("smallest reported ESS = ", format(min(ess))))
    }
    cat(" Convergence: ", paste(pieces, collapse = "; "), "\n", sep = "")
  }
  cat(" Latent-scale interpretation: ", calibration$interpretation, "\n", sep = "")
  invisible(x)
}


#' Fit the unrestricted-correlation GP competitor
#' @inheritParams fit_lvgp
#' @export
fit_ucgp <- function(X,
                     y,
                     C,
                     new_X,
                     new_C,
                     m_vec = NULL,
                     n_draw = 1000L,
                     seed = 1L,
                     ...) {
  n_draw <- mixedgp_as_integer_strict(
    n_draw, "n_draw", min_value = 1L, length_expected = 1L
  )
  seed <- mixedgp_as_integer_strict(
    seed, "seed", min_value = 0L, length_expected = 1L
  )
  caller_rng <- mixedgp_rng_state()
  on.exit(mixedgp_restore_rng_state(caller_rng), add = TRUE)
  mixedgp_do_call(
    mixedgp_adapter_ucgp,
    list(
      X_train = X, y_train = y, C_train = C,
      X_test = new_X, C_test = new_C, m_vec = m_vec,
      n_draw = n_draw, seed = seed
    ),
    list(...)
  )
}


#' Fit the published latent-variable GP competitor
#'
#' @param X Numeric training design matrix.
#' @param y Numeric response vector.
#' @param C Integer-coded ordinal training inputs.
#' @param new_X Numeric test design matrix.
#' @param new_C Integer-coded ordinal test inputs.
#' @param m_vec Number of levels for each ordinal input.
#' @param n_draw Number of noisy-response predictive draws.
#' @param seed Reproducible fitting seed.
#' @param ... Method-specific controls passed to the audited wrapper.
#' @return A list containing the public-package fit and predictive draws.
#' @export
fit_lvgp <- function(X,
                     y,
                     C,
                     new_X,
                     new_C,
                     m_vec = NULL,
                     n_draw = 1000L,
                     seed = 1L,
                     ...) {
  n_draw <- mixedgp_as_integer_strict(
    n_draw, "n_draw", min_value = 1L, length_expected = 1L
  )
  seed <- mixedgp_as_integer_strict(
    seed, "seed", min_value = 0L, length_expected = 1L
  )
  caller_rng <- mixedgp_rng_state()
  on.exit(mixedgp_restore_rng_state(caller_rng), add = TRUE)
  mixedgp_do_call(
    mixedgp_adapter_lvgp,
    list(
      X_train = X, y_train = y, C_train = C,
      X_test = new_X, C_test = new_C, m_vec = m_vec,
      n_draw = n_draw, seed = seed
    ),
    list(...)
  )
}


#' Fit the published easy-to-interpret GP competitor
#' @inheritParams fit_lvgp
#' @export
fit_ezgp <- function(X,
                     y,
                     C,
                     new_X,
                     new_C,
                     m_vec = NULL,
                     n_draw = 1000L,
                     seed = 1L,
                     ...) {
  n_draw <- mixedgp_as_integer_strict(
    n_draw, "n_draw", min_value = 1L, length_expected = 1L
  )
  seed <- mixedgp_as_integer_strict(
    seed, "seed", min_value = 0L, length_expected = 1L
  )
  caller_rng <- mixedgp_rng_state()
  on.exit(mixedgp_restore_rng_state(caller_rng), add = TRUE)
  mixedgp_do_call(
    mixedgp_adapter_ezgp,
    list(
      X_train = X, y_train = y, C_train = C,
      X_test = new_X, C_test = new_C, m_vec = m_vec,
      n_draw = n_draw, seed = seed
    ),
    list(...)
  )
}


#' Fit one audited published mixed-input GP competitor
#' @param method One of `"UC-GP"`, `"LVGP"`, or `"EzGP"`.
#' @param ... Arguments passed to the method-specific single-call wrapper.
#' @details UC-GP tries its original start budget before the fixed
#'   `rescue_starts = c(32L, 64L)` sequence after failure. EzGP similarly uses
#'   `rescue_maxeval = c(300L, 1000L)` with unchanged training CV folds.
#'   Both stop at the first valid fit and record every attempt; no test responses
#'   enter fitting or selection. Set the rescue vector to `integer(0)` to disable it.
#' @export
fit_mixedgp_competitor <- function(method, ...) {
  method <- match.arg(method, c("UC-GP", "LVGP", "EzGP"))
  dots <- list(...)
  mixedgp_validate_named_dots(dots)
  FUN <- switch(method, `UC-GP` = fit_ucgp, LVGP = fit_lvgp,
                EzGP = fit_ezgp)
  do.call(FUN, dots)
}



# ---- core/00_diagnostics.R ----
## Shared export protocol; never recomputes diagnostics or changes their gates.
mixedgp_write_diagnostic_tables <- function(summary, detail, output_dir,
                                            prefix = "eivgp", overwrite = FALSE) {
  if (!is.data.frame(summary) || !is.data.frame(detail))
    stop("summary and detail must be data frames.")
  if (!is.character(prefix) || length(prefix) != 1L || is.na(prefix) ||
      !grepl("^[A-Za-z0-9_-]+$", prefix)) stop("Invalid export prefix.")
  if (!is.logical(overwrite) || length(overwrite) != 1L || is.na(overwrite))
    stop("overwrite must be TRUE or FALSE.")
  if (!is.character(output_dir) || length(output_dir) != 1L ||
      is.na(output_dir) || !nzchar(output_dir)) stop("Invalid output_dir.")
  if (!"parameter" %in% names(detail)) {
    if (nrow(detail)) stop("detail must contain parameter.")
    detail$parameter <- character()
  }
  role <- if ("diagnostic_role" %in% names(detail)) as.character(detail$diagnostic_role) else
    if ("group" %in% names(detail)) as.character(detail$group) else
      if ("target" %in% names(detail)) ifelse(detail$target == "raw_coordinate", "raw",
        ifelse(detail$target == "ordinal_measurement_invariant", "invariant", "target")) else
          rep("unclassified", nrow(detail))
  role[!is.na(role) & role == "raw_and_training_imputation"] <- "raw"
  role[!is.na(role) & role == "measurement_invariant"] <- "invariant"
  role[!is.na(role) & role == "scientific_panel"] <- "target"
  role[!is.na(role) & role == "user_functional"] <- "additional"
  detail$diagnostic_role <- role
  raw <- !is.na(role) & role == "raw"
  target <- !is.na(role) & role %in% c("target", "invariant", "additional")
  tables <- list(mcmc_diagnostics = summary, mcmc_parameter_diagnostics = detail[raw, , drop = FALSE],
    mcmc_target_diagnostics = detail[target, , drop = FALSE], mcmc_diagnostic_details = detail)
  paths <- file.path(output_dir, paste0(prefix, "_", names(tables), ".csv"))
  if (!overwrite && any(file.exists(paths))) stop("Diagnostic export exists; use overwrite = TRUE explicitly.")
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  for (i in seq_along(tables)) {
    tmp <- tempfile(".diagnostics-", tmpdir = output_dir)
    tryCatch({
      utils::write.csv(tables[[i]], tmp, row.names = FALSE, na = "NA")
      if (!file.rename(tmp, paths[i])) stop("Cannot publish diagnostic export: ", paths[i])
    }, finally = if (file.exists(tmp)) unlink(tmp))
  }
  invisible(setNames(paths, names(tables)))
}


## Release diagnostics must retain nonfinite entries: dropping one changes the
## proposition being checked. Only known fixed model coordinates are excluded
## when constructing the series, never because their diagnostics are inconvenient.
mixedgp_diagnostic_table_pass <- function(x, rhat_limit = 1.01,
                                           bulk_ess_limit = 400,
                                           tail_ess_limit = 400,
                                           expected_parameters = NULL,
                                           min_chains = 4L,
                                           mcse_limit = NULL,
                                           mcse_sd_ratio_limit = NULL) {
  required <- c("parameter", "rhat", "ess_bulk", "ess_tail", "n_chains")
  if (!is.data.frame(x) || nrow(x) == 0L ||
      !all(required %in% names(x)) || anyNA(x$parameter) ||
      anyDuplicated(x$parameter)) return(FALSE)
  if (!is.null(expected_parameters) &&
      !setequal(x$parameter, expected_parameters)) return(FALSE)
  if ("draw_window_complete" %in% names(x) &&
      !all(!is.na(x$draw_window_complete) & x$draw_window_complete)) return(FALSE)
  for (item in list(c("mcse", "mcse_limit"),
                    c("mcse_sd_ratio", "mcse_sd_ratio_limit"))) {
    limit <- get(item[2L])
    if (!is.null(limit)) {
      if (!is.numeric(limit) || length(limit) != 1L ||
          !is.finite(limit) || limit <= 0) stop(item[2L], " must be positive and finite.")
      if (!item[1L] %in% names(x) || any(!is.finite(x[[item[1L]]])) ||
          any(x[[item[1L]]] > limit)) return(FALSE)
    }
  }
  values <- as.matrix(x[, c("rhat", "ess_bulk", "ess_tail", "n_chains")])
  all(is.finite(values)) && all(x$n_chains >= min_chains) &&
    all(x$rhat <= rhat_limit) && all(x$ess_bulk >= bulk_ess_limit) &&
    all(x$ess_tail >= tail_ess_limit)
}


mixedgp_summarize_diagnostic_series <- function(series, rhat_limit = 1.01,
                                                bulk_ess_limit = 400,
                                                tail_ess_limit = 400,
                                                mcse_limit = NULL,
                                                mcse_sd_ratio_limit = NULL) {
  if (!requireNamespace("posterior", quietly = TRUE)) {
    stop("Release diagnostics require the posterior package.")
  }
  if (length(series) == 0L) return(data.frame())
  do.call(rbind, mixedgp_stage_lapply(names(series), function(parameter) {
    chains <- lapply(series[[parameter]], as.numeric)
    lens <- lengths(chains)
    values <- unlist(chains, use.names = FALSE)
    valid <- length(chains) >= 2L && length(unique(lens)) == 1L &&
      all(lens >= 4L) && all(is.finite(values))
    mat <- if (valid) do.call(cbind, chains) else NULL
    rhat <- if (valid) as.numeric(posterior::rhat(mat)) else NA_real_
    eb <- if (valid) as.numeric(posterior::ess_bulk(mat)) else NA_real_
    et <- if (valid) as.numeric(posterior::ess_tail(mat)) else NA_real_
    em <- if (valid) as.numeric(posterior::ess_mean(mat)) else NA_real_
    mcse <- if (valid) as.numeric(posterior::mcse_mean(mat)) else NA_real_
    row <- data.frame(
      parameter = parameter, n_chains = length(chains),
      min_draws_per_chain = if (length(lens)) min(lens) else 0L,
      max_draws_per_chain = if (length(lens)) max(lens) else 0L,
      mean = if (length(values)) mean(values) else NA_real_,
      posterior_sd = if (length(values) > 1L) stats::sd(values) else NA_real_,
      rhat = rhat, ess_bulk = eb, ess_tail = et, ess_mean = em,
      draw_window_complete = !identical(attr(series, "draw_window_complete"), FALSE) &&
        !identical(attr(series[[parameter]], "draw_window_complete"), FALSE),
      stringsAsFactors = FALSE
    )
    row$ess <- min(eb, et)
    # Rank-normalized bulk ESS is a mixing screen, not the ESS of a raw mean.
    # Each series is the actual functional being averaged (including second
    # moments and conditional variances), so use its mean-specific MCSE.
    row$mcse <- mcse
    row$mcse_sd_ratio <- if (is.finite(row$posterior_sd) && row$posterior_sd > 0)
      mcse / row$posterior_sd else NA_real_
    row$target_pass <- mixedgp_diagnostic_table_pass(
      row, rhat_limit, bulk_ess_limit, tail_ess_limit,
      mcse_limit = mcse_limit, mcse_sd_ratio_limit = mcse_sd_ratio_limit
    )
    row
  }))
}


## Inf means the entire reporting window. An explicitly truncated window is
## useful for pilot costs, but its diagnostics cannot certify all retained draws.
mixedgp_diagnostic_draw_indices <- function(n_saved, max_draws_per_chain = Inf) {
  if (!is.numeric(max_draws_per_chain) || length(max_draws_per_chain) != 1L ||
      is.na(max_draws_per_chain) || max_draws_per_chain < 1 ||
      (is.finite(max_draws_per_chain) &&
       max_draws_per_chain != floor(max_draws_per_chain))) {
    stop("max_draws_per_chain must be a positive integer or Inf (all draws).")
  }
  tail(seq_len(n_saved), min(n_saved, max_draws_per_chain))
}


mixedgp_mark_diagnostic_window <- function(series, complete) {
  attr(series, "draw_window_complete") <- isTRUE(complete)
  # Per-parameter provenance survives ordinary outer-list c() and [ subsetting.
  for (name in names(series)) attr(series[[name]], "draw_window_complete") <- isTRUE(complete)
  series
}


mixedgp_study2_invariant_series <- function(fit) {
  chains <- fit$mcmc$samples_by_chain
  q <- fit$data$q
  d <- fit$data$d
  series <- list()
  score_sd <- lapply(chains$A, function(a) {
    out <- matrix(NA_real_, dim(a)[1L], q)
    for (j in seq_len(q)) {
      aj <- matrix(a[, j, , drop = FALSE], ncol = d)
      out[, j] <- sqrt(1 + rowSums(aj^2))
    }
    out
  })
  if (q > 1L) for (pair in combn(seq_len(q), 2L, simplify = FALSE)) {
    j <- pair[1L]
    k <- pair[2L]
    name <- paste0("ordinal_score_cor[", j, ",", k, "]")
    series[[name]] <- lapply(seq_along(chains$A), function(id) {
      a <- chains$A[[id]]
      aj <- matrix(a[, j, , drop = FALSE], ncol = d)
      ak <- matrix(a[, k, , drop = FALSE], ncol = d)
      rowSums(aj * ak) / (score_sd[[id]][, j] * score_sd[[id]][, k])
    })
  }
  for (j in seq_len(q)) for (r in seq_len(fit$data$m_vec[j] - 1L)) {
    name <- paste0("tau[", j, ",", r, "]")
    column <- match(name, colnames(chains$tau[[1L]]))
    if (is.na(column)) stop("Missing diagnostic cutpoint: ", name)
    series[[paste0("standardized_", name)]] <- lapply(
      seq_along(chains$tau), function(id) {
        chains$tau[[id]][, column] / score_sd[[id]][, j]
      }
    )
  }
  series
}


mixedgp_study2_raw_series <- function(fit) {
  chains <- fit$mcmc$samples_by_chain
  series <- list(sigma_epsilon = lapply(chains$sigma2, sqrt))
  theta_names <- c("rho", paste0("theta_x", seq_len(fit$data$p)),
                   paste0("theta_u", seq_len(fit$data$d)))
  for (k in seq_along(theta_names)) {
    series[[theta_names[k]]] <- lapply(chains$logtheta, function(z) exp(z[, k]))
  }
  for (j in seq_len(fit$data$q)) for (k in seq_len(fit$data$d)) {
    fixed <- identical(fit$data$ident, "lower_triangular") &&
      j <= fit$data$d && k > j
    if (!fixed) series[[paste0("A[", j, ",", k, "]")]] <-
      lapply(chains$A, function(a) a[, j, k])
  }
  for (k in seq_len(ncol(chains$tau[[1L]]))) {
    series[[colnames(chains$tau[[1L]])[k]]] <-
      lapply(chains$tau, function(z) z[, k])
  }
  for (i in fit$data$miss_idx) for (k in seq_len(fit$data$d)) {
    series[[paste0("U[", i, ",", k, "]")]] <-
      lapply(chains$U, function(z) z[, i, k])
  }
  series
}


mixedgp_study1_raw_series <- function(fit) {
  chains <- fit$mcmc$samples_by_chain
  p <- ncol(fit$data$x)
  names_theta <- c("rho", if (p == 1L) "theta_x" else paste0("theta_x", seq_len(p)),
                   "theta_u")
  series <- list(sigma_epsilon = lapply(chains$sigma2, sqrt))
  for (k in seq_along(names_theta)) {
    series[[names_theta[k]]] <- lapply(chains$logtheta, function(z) exp(z[, k]))
  }
  for (k in seq_len(ncol(chains$tau[[1L]]))) {
    series[[paste0("tau", k)]] <- lapply(chains$tau, function(z) z[, k])
  }
  for (i in fit$data$miss_idx) {
    series[[paste0("u[", i, "]")]] <- lapply(chains$u, function(z) z[, i])
  }
  series
}


mixedgp_study1_raw_diagnostics <- function(fit, rhat_limit = 1.01,
                                            bulk_ess_limit = 400,
                                            tail_ess_limit = 400) {
  mixedgp_summarize_diagnostic_series(
    mixedgp_study1_raw_series(fit), rhat_limit, bulk_ess_limit, tail_ess_limit
  )
}


mixedgp_study1_target_series <- function(fit, X, C, U = NULL,
                                         max_draws_per_chain = Inf,
                                         n_latent = 64L, seed = 481517L) {
  caller_rng <- mixedgp_rng_state()
  on.exit(mixedgp_restore_rng_state(caller_rng), add = TRUE)
  X <- as.matrix(X)
  C <- as.integer(C)
  G <- nrow(X)
  M <- as.integer(n_latent)
  if (G < 1L || length(C) != G || anyNA(C) ||
      any(C < 1L | C > fit$data$m) || M < 2L) stop("Invalid Study I diagnostic panel.")
  X_std <- sweep(sweep(X, 2L, fit$data$x_center, "-"),
                 2L, fit$data$x_scale, "/")
  colnames(X_std) <- colnames(fit$data$x)
  calibrated <- length(fit$data$calib_idx) > 0L
  U_std <- if (!is.null(U)) {
    uc <- if (is.null(fit$data$U_center)) 0 else fit$data$U_center
    us <- if (is.null(fit$data$U_scale)) 1 else fit$data$U_scale
    (as.numeric(U) - uc) / us
  } else NULL
  if (!is.null(U_std) && (length(U_std) != G || any(!is.finite(U_std)))) {
    stop("U does not match the Study I diagnostic panel.")
  }
  nms <- c(paste0("m_conditional_mean[", seq_len(G), "]"),
           paste0("m_conditional_variance[", seq_len(G), "]"),
           paste0("predictive_second_moment[", seq_len(G), "]"))
  if (calibrated) nms <- c(nms,
    paste0("prospective_U_mean[", seq_len(G), "]"),
    paste0("prospective_U_second_moment[", seq_len(G), "]"))
  if (calibrated && !is.null(U_std)) nms <- c(nms,
    paste0("f_conditional_mean[", seq_len(G), "]"),
    paste0("f_conditional_variance[", seq_len(G), "]"))
  chains <- fit$mcmc$samples_by_chain
  series <- setNames(lapply(nms, function(z) vector("list", length(chains$u))), nms)
  attr(series, "draw_window_complete") <- TRUE
  kernel <- kernel_spec_from_fit_1d(fit)
  repeated <- rep(seq_len(G), each = M)
  for (id in seq_along(chains$u)) {
    n_saved <- nrow(chains$u[[id]])
    keep <- mixedgp_diagnostic_draw_indices(n_saved, max_draws_per_chain)
    if (length(keep) < n_saved) attr(series, "draw_window_complete") <- FALSE
    values <- matrix(NA_real_, length(keep), length(nms), dimnames = list(NULL, nms))
    for (ii in seq_along(keep)) {
      s <- keep[ii]
      tau <- chains$tau[[id]][s, ]
      # Common integration randomness avoids injecting fresh outcome/GP noise.
      # This finite-M approximation still requires an integration sensitivity check.
      set.seed(seed)
      u_mc <- matrix(rtruncnorm_vec(0, 1,
        c(-Inf, tau)[C[repeated]], c(tau, Inf)[C[repeated]]), M, G)
      args_gp <- list(x_train = fit$data$x, u_train = chains$u[[id]][s, ],
        y_train = fit$data$y, logtheta = chains$logtheta[[id]][s, ],
        sigma2_eps = chains$sigma2[[id]][s], kernel = kernel$name,
        matern_nu = kernel$matern_nu)
      mp <- do.call(gp_integrated_mean_state_1d,
                    c(args_gp, list(x_star = X_std, u_mc = u_mc)))
      values[ii, paste0("m_conditional_mean[", seq_len(G), "]")] <-
        fit$data$y_center + fit$data$y_scale * mp$mean
      values[ii, paste0("m_conditional_variance[", seq_len(G), "]")] <-
        fit$data$y_scale^2 * mp$var
      yp <- do.call(gp_predict_draw, c(args_gp, list(
        x_star = X_std[repeated, , drop = FALSE], u_star = as.numeric(u_mc),
        noisy = TRUE)))
      mu <- fit$data$y_center + fit$data$y_scale * yp$mean
      values[ii, paste0("predictive_second_moment[", seq_len(G), "]")] <-
        colMeans(matrix(fit$data$y_scale^2 * yp$var + mu^2, M, G))
      if (calibrated) {
        uc <- if (is.null(fit$data$U_center)) 0 else fit$data$U_center
        us <- if (is.null(fit$data$U_scale)) 1 else fit$data$U_scale
        u_physical <- uc + us * u_mc
        values[ii, paste0("prospective_U_mean[", seq_len(G), "]")] <- colMeans(u_physical)
        values[ii, paste0("prospective_U_second_moment[", seq_len(G), "]")] <- colMeans(u_physical^2)
      }
      if (calibrated && !is.null(U_std)) {
        fp <- do.call(gp_predict_draw, c(args_gp, list(x_star = X_std, u_star = U_std)))
        values[ii, paste0("f_conditional_mean[", seq_len(G), "]")] <-
          fit$data$y_center + fit$data$y_scale * fp$mean
        values[ii, paste0("f_conditional_variance[", seq_len(G), "]")] <-
          fit$data$y_scale^2 * fp$var
      }
    }
    for (name in nms) series[[name]][[id]] <- values[, name]
  }
  mixedgp_mark_diagnostic_window(series, attr(series, "draw_window_complete"))
}


mixedgp_study2_diagnostic_rows <- function(X, C, max_points = 5L) {
  X <- as.matrix(X)
  C <- as.matrix(C)
  keys <- apply(C, 1L, paste, collapse = ":")
  counts <- sort(table(keys), decreasing = TRUE)
  patterns <- names(counts)[unique(round(seq(
    1L, length(counts), length.out = min(3L, length(counts))
  )))]
  rows <- vapply(patterns, function(key) which(keys == key)[1L], integer(1))
  for (j in seq_len(ncol(X))) rows <- c(rows, which.min(X[, j]), which.max(X[, j]))
  rows <- unique(c(rows, seq_len(nrow(X))))
  head(rows, min(as.integer(max_points), nrow(X)))
}


## Fixed common integration randomness makes each signature a deterministic
## function of the retained state. Fresh observation or GP noise must not make
## a poorly mixing parameter chain appear to have a large effective sample size.
mixedgp_study2_target_series <- function(fit, X, C, U = NULL,
                                          max_draws_per_chain = Inf,
                                          n_latent = 64L, seed = 481516L) {
  caller_rng <- mixedgp_rng_state()
  on.exit(mixedgp_restore_rng_state(caller_rng), add = TRUE)
  X <- as.matrix(X)
  C <- as.matrix(C)
  G <- nrow(X)
  d <- fit$data$d
  M <- as.integer(n_latent)
  if (G < 1L || M < 2L || nrow(C) != G) stop("Invalid target diagnostic panel.")
  X_std <- sweep(sweep(X, 2L, fit$data$X_center, "-"),
                 2L, fit$data$X_scale, "/")
  repeated <- rep(seq_len(G), each = M)
  calibrated <- length(fit$data$calib_idx) > 0L
  U_std <- if (!is.null(U)) {
    uc <- if (is.null(fit$data$U_center)) rep(0, d) else fit$data$U_center
    us <- if (is.null(fit$data$U_scale)) rep(1, d) else fit$data$U_scale
    sweep(sweep(as.matrix(U), 2L, uc, "-"), 2L, us, "/")
  } else NULL
  if (!is.null(U_std) && !identical(dim(U_std), c(G, d))) {
    stop("U does not match the target diagnostic panel.")
  }
  by_chain <- fit$mcmc$samples_by_chain
  names_target <- c(paste0("m_conditional_mean[", seq_len(G), "]"),
                    paste0("m_conditional_variance[", seq_len(G), "]"),
                    paste0("predictive_second_moment[", seq_len(G), "]"))
  if (calibrated && !is.null(U_std)) names_target <- c(
    names_target, paste0("f_conditional_mean[", seq_len(G), "]"),
    paste0("f_conditional_variance[", seq_len(G), "]")
  )
  if (calibrated) for (k in seq_len(d)) names_target <- c(
    names_target, paste0("prospective_U_mean[", seq_len(G), ",", k, "]"),
    paste0("prospective_U_second_moment[", seq_len(G), ",", k, "]")
  )
  series <- setNames(lapply(names_target, function(z) vector("list", length(by_chain$U))),
                     names_target)
  attr(series, "draw_window_complete") <- TRUE
  kernel <- kernel_spec_from_fit(fit)
  for (id in seq_along(by_chain$U)) {
    n_saved <- dim(by_chain$U[[id]])[1L]
    keep <- mixedgp_diagnostic_draw_indices(n_saved, max_draws_per_chain)
    if (length(keep) < n_saved) attr(series, "draw_window_complete") <- FALSE
    values <- matrix(NA_real_, length(keep), length(names_target),
                     dimnames = list(NULL, names_target))
    for (ii in seq_along(keep)) {
      s <- keep[ii]
      tau <- unflatten_tau(by_chain$tau[[id]][s, ], fit$data$m_vec)
      A <- matrix(by_chain$A[[id]][s, , ], fit$data$q, d)
      set.seed(seed)
      u_mc <- sample_u_given_c_ordprobit_minimax(C[repeated, , drop = FALSE],
                                                A, tau, fit$data$m_vec)
      U_train <- matrix(by_chain$U[[id]][s, , ], nrow(fit$data$X), d)
      pred <- gp_predict_draw_general(
        X_train = fit$data$X, U_train = U_train, y_train = fit$data$y,
        X_star = X_std[repeated, , drop = FALSE], U_star = u_mc,
        logtheta = by_chain$logtheta[[id]][s, ],
        sigma2_eps = by_chain$sigma2[[id]][s], noisy = TRUE,
        kernel = kernel$name, matern_nu = kernel$matern_nu, return_cov = FALSE
      )
      mu <- fit$data$y_center + fit$data$y_scale * pred$mean
      variance <- fit$data$y_scale^2 * pred$var
      mp <- gp_integrated_mean_state_general(
        X_train = fit$data$X, U_train = U_train, y_train = fit$data$y,
        X_star = X_std, U_mc = array(u_mc, c(M, G, d)),
        logtheta = by_chain$logtheta[[id]][s, ],
        sigma2_eps = by_chain$sigma2[[id]][s], kernel = kernel$name,
        matern_nu = kernel$matern_nu
      )
      values[ii, paste0("m_conditional_mean[", seq_len(G), "]")] <-
        fit$data$y_center + fit$data$y_scale * mp$mean
      values[ii, paste0("m_conditional_variance[", seq_len(G), "]")] <-
        fit$data$y_scale^2 * mp$var
      values[ii, paste0("predictive_second_moment[", seq_len(G), "]")] <-
        colMeans(matrix(variance + mu^2, M, G))
      if (calibrated && !is.null(U_std)) {
        fp <- gp_predict_draw_general(
          X_train = fit$data$X, U_train = U_train, y_train = fit$data$y,
          X_star = X_std, U_star = U_std,
          logtheta = by_chain$logtheta[[id]][s, ],
          sigma2_eps = by_chain$sigma2[[id]][s], noisy = FALSE,
          kernel = kernel$name, matern_nu = kernel$matern_nu, return_cov = FALSE
        )
        values[ii, paste0("f_conditional_mean[", seq_len(G), "]")] <-
          fit$data$y_center + fit$data$y_scale * fp$mean
        values[ii, paste0("f_conditional_variance[", seq_len(G), "]")] <-
          fit$data$y_scale^2 * fp$var
      }
      if (calibrated) for (k in seq_len(d)) {
        uc <- if (is.null(fit$data$U_center)) 0 else fit$data$U_center[k]
        us <- if (is.null(fit$data$U_scale)) 1 else fit$data$U_scale[k]
        u_physical <- uc + us * u_mc[, k]
        values[ii, paste0("prospective_U_mean[", seq_len(G), ",", k, "]")] <-
          colMeans(matrix(u_physical, M, G))
        values[ii, paste0("prospective_U_second_moment[", seq_len(G), ",", k, "]")] <-
          colMeans(matrix(u_physical^2, M, G))
      }
    }
    for (name in names_target) series[[name]][[id]] <- values[, name]
  }
  mixedgp_mark_diagnostic_window(series, attr(series, "draw_window_complete"))
}



# ---- core/00_mcmc_workflow.R ----
## Public, explicit diagnostic-and-extension workflow. No automatic stopping
## or convergence claim: finite-chain diagnostics can miss common-mode failures.

#' Extend the existing posterior chains
#'
#' @param object An unmodified fit returned by [fit_eivgp()] in version 0.3.1,
#'   or by a previous call to this function. Earlier sampler checkpoints are
#'   incompatible with this transition revision. Version 0.3.0 targets the
#'   same posterior, but reproduce/continue those chains with that version,
#'   or start a fresh fit to use the new transitions.
#' @param n_iter Additional transitions per chain, not additional saved draws.
#'   Every post-warm-up transition is retained; there is no thinning.
#' @param parallel Whether to run chains in parallel. NULL keeps the fit setting.
#' @param n_cores Worker count. NULL keeps the fit setting.
#' @param verbose Whether to print sampler progress.
#' @return An updated `eivgp_fit`, containing both old and new retained draws.
#' @details The terminal state, random-number stream, update schedule,
#'   prior dictionary, and sampler settings are preserved. No new
#'   warm-up is discarded and no model or sampler options can be changed here.
#'   Checkpoints survive `saveRDS()`/`readRDS()`; use the same package version
#'   and numerical environment for reproducible continuation. Older fits lacking
#'   checkpoints must be refitted. Diagnostics in the returned fit are recomputed
#'   over all retained draws; call [diagnose_eivgp()] again for target diagnostics.
#'   Running longer is not a remedy for unidentified coordinates or chains
#'   trapped in different modes. This is continuation of a completed call, not
#'   periodic crash recovery during an interrupted fit.
#' @export
continue_eivgp <- function(object, n_iter, parallel = NULL, n_cores = NULL,
                           verbose = FALSE) {
  if (!inherits(object, "eivgp_fit")) stop("object must be an eivgp_fit.")
  caller_rng <- mixedgp_rng_state()
  on.exit(mixedgp_restore_rng_state(caller_rng), add = TRUE)
  n_iter <- mixedgp_as_integer_strict(n_iter, "n_iter", 1L, 1L)
  verbose <- mixedgp_validate_flag(verbose, "verbose")
  cp <- object$checkpoint
  if (is.null(cp) || is.null(cp$arguments)) {
    stop("This fit has no supported continuation checkpoint; refit with version 0.3.1.")
  }
  if (!identical(cp$version, 2L) || !identical(cp$sampler_version, "0.3.1") ||
      !identical(object$sampler_version, "0.3.1")) {
    stop("Incompatible sampler checkpoint: continuation requires a version 0.3.1 fit; earlier fits must be refitted.")
  }
  if (!identical(cp$control, object$control)) stop("The fit controls were modified; cannot continue.")
  if (!identical(cp$signature, mixedgp_v030_checkpoint_hash(object))) {
    stop("The fit data, draws, or checkpoint were modified; cannot safely continue.")
  }
  if (cp$thin != 1L) stop("Thinned historical fits cannot recover discarded draws; refit without thinning.")
  total <- mixedgp_as_integer_strict(as.double(cp$iteration) + n_iter,
                                   "total iterations", 2L, 1L)
  expected <- floor((cp$iteration - cp$burn) / cp$thin)
  if (floor((total - cp$burn) / cp$thin) <= expected) {
    stop("The extension must reach at least one new retained draw per chain.")
  }
  by_chain <- object$mcmc$samples_by_chain
  n_chains <- length(cp$states)
  if (!n_chains || length(by_chain$logtheta) != n_chains ||
      any(vapply(by_chain$logtheta, nrow, integer(1)) != expected) ||
      nrow(object$mcmc$chain_stats) != n_chains) {
    stop("The retained chain history does not match its checkpoint.")
  }
  cp$old_chains <- lapply(seq_len(n_chains), function(i) {
    out <- lapply(by_chain, function(chains) chains[[i]])
    names(out) <- paste0("samples_", names(out))
    out$stats <- object$mcmc$chain_stats[i, , drop = FALSE]
    if (!is.null(object$mcmc$chain_initial)) out$initial_state <- object$mcmc$chain_initial[[i]]
    out
  })
  if (is.null(parallel)) parallel <- object$interface$parallel
  parallel <- mixedgp_validate_flag(parallel, "parallel")
  if (is.null(n_cores)) n_cores <- object$interface$n_cores
  if (!is.null(n_cores)) n_cores <- mixedgp_as_integer_strict(n_cores, "n_cores", 1L, 1L)
  arguments <- cp$arguments
  arguments$n_iter <- total
  arguments$parallel_chains <- parallel
  arguments$n_cores <- n_cores
  arguments$verbose <- verbose
  engine <- mixedgp_fit_engine(object)
  sampler <- if (engine == "univariate") fit_eivgp_1d_v030 else fit_eivgp_ordprobit_v030
  updated <- do.call(sampler, c(arguments, list(.resume = cp)))
  updated$checkpoint$arguments <- arguments
  updated$data <- object$data
  updated$call <- object$call
  updated$model_specification <- object$model_specification
  updated$interface <- object$interface
  updated$interface$parallel <- parallel
  updated$interface$n_cores <- updated$diagnostics$summary$parallel_cores
  prior_time <- object$diagnostics$summary$time_seconds
  updated$diagnostics$summary$time_seconds <- prior_time + updated$diagnostics$summary$time_seconds
  updated$continuation_history <- rbind(object$continuation_history, data.frame(
    from_iteration = cp$iteration, to_iteration = total,
    added_draws_per_chain = floor((total - cp$burn) / cp$thin) - expected
  ))
  class(updated) <- class(object)
  updated$checkpoint$signature <- mixedgp_v030_checkpoint_hash(updated)
  updated
}


mixedgp_diagnostic_status <- function(tab, rhat_limit, bulk_ess_limit,
                                      tail_ess_limit, mcse_sd_ratio_limit) {
  valid <- is.finite(tab$rhat) & is.finite(tab$ess_bulk) & is.finite(tab$ess_tail) &
    is.finite(tab$mcse_sd_ratio) & tab$n_chains >= 4L & tab$draw_window_complete
  status <- rep("insufficient_diagnostics", nrow(tab))
  status[valid] <- "screen_passed"
  precision <- tab$ess_bulk < bulk_ess_limit | tab$ess_tail < tail_ess_limit |
    tab$mcse_sd_ratio > mcse_sd_ratio_limit
  status[valid & precision] <- "insufficient_precision"
  status[valid & tab$rhat > rhat_limit] <- "poor_exploration"
  status
}


#' Diagnose posterior exploration and Monte Carlo precision
#'
#' @param object An `eivgp_fit`.
#' @param X,C Diagnostic panel in the same units/ordinal labels as fitting data.
#'   Supply both or neither. By default at most five training-input locations
#'   are chosen deterministically; these are locations for evaluating m(x,c),
#'   not a claim to cover the response surface.
#' @param U Optional exact latent inputs at the panel locations, on the original
#'   calibration scale, to additionally monitor f(x,u). Requires the package's
#'   calibrated-scale reporting gate to be satisfied.
#' @param n_latent Integration points per location and posterior state.
#' @param seed Seed for common integration randomness, not for the fitted chains.
#' @param rhat_limit Maximum rank-normalized split R-hat (default 1.01).
#' @param bulk_ess_limit,tail_ess_limit Minimum bulk/tail ESS (default 400 each).
#' @param mcse_sd_ratio_limit Maximum MCSE of each monitored mean relative to
#'   its posterior standard deviation (default 0.05). This is a screening
#'   tolerance, not a decision-specific absolute accuracy guarantee.
#' @param additional_series Optional uniquely named list of scientific
#'   functionals. Each entry must be a list of one numeric vector per chain,
#'   aligned with every retained draw in that chain. For example supply a
#'   scientifically important contrast. No thinning or truncation is allowed.
#' @return An `eivgp_diagnostics` with per-functional `table`, overall `status`,
#'   `recommendation`, diagnostic `panel`, and `settings`.
#' @details All retained draws are checked, using at least four chains. Reports
#'   include free parameters, missing training U coordinates and their squares,
#'   multivariate measurement invariants, and panel-specific conditional means,
#'   conditional GP variances, and predictive second moments. These quantities
#'   avoid adding independent predictive noise that can disguise stuck chains.
#'   Mean MCSE is computed with `posterior::mcse_mean`, not bulk ESS.
#'   Common finite integration randomness makes target signatures deterministic
#'   functions of chain states, but integration error is not included in MCSE:
#'   repeat with larger `n_latent` and another integration seed for sensitivity.
#'   Quantile MCSE, simultaneous bands, and decision losses require separate
#'   checks. Passing this screen proves neither convergence, identification,
#'   nor model adequacy. Uncalibrated latent-coordinate failures may reflect
#'   symmetries; the raw and invariant panels are reported separately, without
#'   silently dismissing an unsuccessful raw-coordinate check.
#' @export
diagnose_eivgp <- function(object, X = NULL, C = NULL, U = NULL,
                           n_latent = 64L, seed = 481517L,
                           rhat_limit = 1.01, bulk_ess_limit = 400,
                           tail_ess_limit = 400, mcse_sd_ratio_limit = 0.05,
                           additional_series = list()) {
  if (!inherits(object, "eivgp_fit")) stop("object must be an eivgp_fit.")
  for (name in c("rhat_limit", "bulk_ess_limit", "tail_ess_limit", "mcse_sd_ratio_limit")) {
    value <- get(name)
    if (!is.numeric(value) || length(value) != 1L || !is.finite(value) || value <= 0) {
      stop(name, " must be positive and finite.")
    }
  }
  if (rhat_limit < 1) stop("rhat_limit must be at least 1.")
  n_latent <- mixedgp_as_integer_strict(n_latent, "n_latent", 2L, 1L)
  seed <- mixedgp_as_integer_strict(seed, "seed", 0L, 1L)
  engine <- mixedgp_fit_engine(object)
  if (xor(is.null(X), is.null(C))) stop("Supply both X and C, or neither.")
  if (is.null(X)) {
    X <- if (engine == "univariate") object$data$x_raw else object$data$X_raw
    C <- if (engine == "univariate") matrix(object$data$c_ord, ncol = 1L) else object$data$C_ord
    rows <- mixedgp_study2_diagnostic_rows(X, C)
    X <- X[rows, , drop = FALSE]
    C <- C[rows, , drop = FALSE]
    if (!is.null(U)) stop("Supply explicit X and C when supplying U.")
  } else {
    X <- as_numeric_matrix_strict(X, "X")
    encoded <- mixedgp_encode_ordinal(C,
      if (engine == "univariate") object$data$m else object$data$m_vec,
      object$data$C_level_maps)
    C <- encoded$C
  }
  expected_p <- if (engine == "univariate") ncol(object$data$x) else object$data$p
  expected_q <- if (engine == "univariate") 1L else object$data$q
  if (ncol(X) != expected_p || nrow(X) != nrow(C) || ncol(C) != expected_q ||
      nrow(X) < 1L || any(!is.finite(X))) stop("Invalid diagnostic panel dimensions or values.")
  if (!is.null(U)) {
    if (!isTRUE(object$data$latent_scale_anchored)) {
      stop("Physical-scale surface diagnostics require adequate calibration anchoring.")
    }
    U <- as_numeric_matrix_strict(U, "U", nrow_expected = nrow(X))
    if (ncol(U) != if (engine == "univariate") 1L else object$data$d) {
      stop("U has the wrong latent dimension.")
    }
  }
  raw <- if (isTRUE(object$sampler_version %in% c("0.3.0", "0.3.1"))) mixedgp_v030_raw_series(object) else
    if (engine == "univariate") mixedgp_study1_raw_series(object) else mixedgp_study2_raw_series(object)
  latent_names <- names(raw)[grepl("^[uU]\\[", names(raw))]
  for (name in latent_names) raw[[paste0(name, "^2")]] <- lapply(raw[[name]], function(z) z^2)
  invariant <- if (engine == "multivariate") mixedgp_study2_invariant_series(object) else list()
  targets <- if (engine == "univariate") {
    mixedgp_study1_target_series(object, X, C, U, n_latent = n_latent, seed = seed)
  } else {
    mixedgp_study2_target_series(object, X, C, U, n_latent = n_latent, seed = seed)
  }
  series <- c(raw, invariant, targets)
  mixedgp_validate_named_dots(additional_series, "additional_series")
  if (length(intersect(names(series), names(additional_series)))) {
    stop("additional_series names must not duplicate built-in functionals.")
  }
  lens <- vapply(object$mcmc$samples_by_chain$logtheta, nrow, integer(1))
  for (name in names(additional_series)) {
    z <- additional_series[[name]]
    if (!is.list(z) || length(z) != length(lens) ||
        !identical(as.integer(lengths(z)), as.integer(lens)) ||
        !all(vapply(z, function(v) is.numeric(v) && is.null(dim(v)), logical(1)))) {
      stop("Each additional functional must match all retained draws in every chain.")
    }
  }
  series <- c(series, additional_series)
  tab <- mixedgp_summarize_diagnostic_series(series, rhat_limit, bulk_ess_limit,
    tail_ess_limit, mcse_sd_ratio_limit = mcse_sd_ratio_limit)
  tab$group <- rep(c("raw_and_training_imputation", "measurement_invariant",
                     "scientific_panel", "user_functional"),
                   c(length(raw), length(invariant), length(targets), length(additional_series)))
  tab$status <- mixedgp_diagnostic_status(tab, rhat_limit, bulk_ess_limit,
                                        tail_ess_limit, mcse_sd_ratio_limit)
  priority <- c("poor_exploration", "insufficient_diagnostics", "insufficient_precision", "screen_passed")
  status <- priority[priority %in% tab$status][1L]
  recommendation <- switch(status,
    poor_exploration = paste("Inspect trace plots, chain-specific targets, calibration and latent symmetries.",
      "Do not assume extending alone will resolve poor exploration."),
    insufficient_diagnostics = paste("Use at least four chains and more retained draws; inspect nonfinite",
      "or constant functionals. A diagnostic pass cannot be issued."),
    insufficient_precision = paste("Chain agreement meets the R-hat screen but precision is insufficient.",
      "Consider continue_eivgp(), then recompute this report on all draws."),
    screen_passed = paste("The monitored functionals pass the specified screen. Check integration sensitivity,",
      "decision-specific MC error and model adequacy before scientific reporting."))
  structure(list(status = status, recommendation = recommendation, table = tab,
    panel = list(X = X, C = C, U = U),
    settings = list(rhat_limit = rhat_limit, bulk_ess_limit = bulk_ess_limit,
      tail_ess_limit = tail_ess_limit, mcse_sd_ratio_limit = mcse_sd_ratio_limit,
      n_latent = n_latent, seed = seed, full_retained_window = TRUE,
      physical_scale_anchored = isTRUE(object$data$latent_scale_anchored))),
    class = "eivgp_diagnostics")
}


#' Export an existing diagnostic report without refitting
#'
#' @param object An `eivgp_diagnostics` returned by [diagnose_eivgp()].
#' @param output_dir Destination directory.
#' @param prefix File prefix, containing letters, digits, underscores or hyphens.
#' @param overwrite Explicitly allow replacement of existing exports.
#' @return Invisibly, named paths to summary, raw-parameter, target and complete
#'   diagnostic CSV files and an RDS containing the full diagnostic report.
#' @details Undefined diagnostics are retained as NA, not interpreted as passes.
#'   The RDS preserves settings, panel locations and recommendations. No fitting,
#'   draw selection, thinning or automatic chain extension is performed.
#' @export
write_diagnostics_eivgp <- function(object, output_dir, prefix = "eivgp", overwrite = FALSE) {
  if (!inherits(object, "eivgp_diagnostics")) stop("object must be an eivgp_diagnostics.")
  if (!is.character(prefix) || length(prefix) != 1L || is.na(prefix) ||
      !grepl("^[A-Za-z0-9_-]+$", prefix)) stop("Invalid export prefix.")
  rds <- file.path(output_dir, paste0(prefix, "_diagnostics.rds"))
  if (!isTRUE(overwrite) && file.exists(rds)) stop("Diagnostic export exists; use overwrite = TRUE explicitly.")
  paths <- mixedgp_write_diagnostic_tables(
    data.frame(status = object$status, recommendation = object$recommendation),
    object$table, output_dir, prefix, overwrite)
  tmp <- tempfile(".diagnostics-", tmpdir = output_dir)
  tryCatch({
    saveRDS(object, tmp)
    if (!file.rename(tmp, rds)) stop("Cannot publish diagnostic report.")
  }, finally = if (file.exists(tmp)) unlink(tmp))
  invisible(c(paths, report = rds))
}


#' @export
print.eivgp_diagnostics <- function(x, ...) {
  cat("EIV-GP diagnostic screen:", x$status, "\n")
  print(table(x$table$group, x$table$status))
  cat(x$recommendation, "\n")
  invisible(x)
}

# Extend the internal Study I format without manufacturing a public-fit signature.
mixedgp_continue_study1_fit <- function(object, n_iter, parallel=FALSE, n_cores=1L) {
  if (inherits(object,"eivgp_fit")) return(continue_eivgp(object,n_iter,parallel,n_cores))
  cp <- object$checkpoint
  if (!identical(object$sampler_version,"0.3.1") || !identical(cp$version,2L) ||
      !identical(cp$sampler_version,"0.3.1") || !identical(cp$thin,1L) ||
      !identical(cp$control,object$control) || mixedgp_fit_engine(object)!="univariate")
    stop("Unsupported or modified Study I continuation state.")
  n_iter <- mixedgp_as_integer_strict(n_iter,"n_iter",1L,1L)
  by <- object$mcmc$samples_by_chain
  expected <- cp$iteration-cp$burn
  nc <- length(cp$states)
  if (nc<1L || any(vapply(by$u,nrow,integer(1))!=expected) || length(by$u)!=nc ||
      nrow(object$mcmc$chain_stats)!=nc || any(vapply(cp$states,function(s)is.null(s$rng)||is.null(s$state),logical(1))))
    stop("Incomplete retained history or terminal chain states.")
  cp$old_chains <- lapply(seq_len(nc),function(i) {
    out<-lapply(by,`[[`,i);names(out)<-paste0("samples_",names(out))
    out$stats<-object$mcmc$chain_stats[i,,drop=FALSE]
    out$initial_state<-object$mcmc$chain_initial[[i]];out
  })
  rng<-mixedgp_rng_state();on.exit(mixedgp_restore_rng_state(rng),add=TRUE)
  pr<-object$priors;pr$dictionary_type<-NULL
  d<-object$data
  result<-mixedgp_v030_fit(X=d$x,y=d$y,C=matrix(d$c_ord,ncol=1L),
    U_obs=matrix(d$u_obs,ncol=1L),m_vec=d$m,measurement="threshold",ident="none",
    kernel=object$kernel$name,matern_nu=object$kernel$matern_nu,
    n_iter=cp$iteration+n_iter,burn=cp$burn,n_chains=nc,
    seed=as.integer(object$mcmc$chain_stats$seed[1]-10000L),
    parallel_chains=parallel,n_cores=n_cores,priors=pr,sampler_control=object$control,.resume=cp)
  updated<-mixedgp_v030_fit_object(result,d,"univariate",object$kernel$name,object$kernel$matern_nu,cp)
  updated$diagnostics$summary$time_seconds<-object$diagnostics$summary$time_seconds+updated$diagnostics$summary$time_seconds
  updated$continuation_history<-rbind(object$continuation_history,data.frame(
    from_iteration=cp$iteration,to_iteration=cp$iteration+n_iter,added_draws_per_chain=n_iter))
  updated
}

mixedgp_scientific_ess_pass <- function(tab,ess_min=100) {
  is.data.frame(tab) && nrow(tab)>0L && all(c("ess_bulk","ess_tail","n_chains","draw_window_complete") %in% names(tab)) &&
    isTRUE(all(is.finite(tab$ess_bulk) & is.finite(tab$ess_tail) & tab$ess_bulk>ess_min &
          tab$ess_tail>ess_min & tab$n_chains>=4L & tab$draw_window_complete))
}

#' Continue saved chains until scientific-target ESS exceeds a threshold
#'
#' @param object A public eivgp_fit or supported internal Study I or Study II fit.
#' @param X,C,U Fixed scientific diagnostic panel, in original units.
#' @param ess_min Strict lower bound for both bulk and tail ESS (default 100).
#' @param chunk_iter Additional transitions per chain between checks.
#' @param adaptive Estimate the next extension from the weakest ESS, bounded
#'   between chunk_iter and four times chunk_iter; always reassess afterwards.
#' @param max_extensions Optional operational limit; Inf continues until success.
#' @param n_latent,seed Fixed integration settings for target functionals.
#' @param parallel,n_cores Chain execution controls.
#' @param checkpoint_file Optional RDS path, atomically saved after each extension
#'   and diagnostic check. Contains the fit and full target series; pass its fit
#'   back to resume after interruption. Terminal states remain in successful fits.
#' @return A list with fit, target_series, diagnostics, history, and success.
#' @details All retained draws count, including the original warm-up exclusion.
#'   No fresh warm-up, thinning or chain restart is performed. R-hat and MCSE
#'   remain reported but do not control this precision stopping rule. Success is
#'   not proof of convergence. Undefined ESS never passes. Inspect the traces.
#' @export
continue_eivgp_until_ess <- function(object,X,C,U=NULL,ess_min=100,chunk_iter=5000L,
    max_extensions=Inf,n_latent=64L,seed=481517L,parallel=FALSE,n_cores=1L,
    checkpoint_file=NULL, adaptive=FALSE) {
  if(length(ess_min)!=1L || !is.finite(ess_min) || ess_min<=0)stop("ess_min must be positive and finite.")
  chunk_iter<-mixedgp_as_integer_strict(chunk_iter,"chunk_iter",1L,1L)
  if(length(max_extensions)!=1L || is.na(max_extensions) || max_extensions<0 ||
     (is.finite(max_extensions)&&max_extensions!=floor(max_extensions)))stop("Invalid max_extensions.")
  if(length(object$mcmc$samples_by_chain$logtheta)<4L)stop("Scientific ESS requires at least four chains.")
  rng<-mixedgp_rng_state();on.exit(mixedgp_restore_rng_state(rng),add=TRUE)
  state<-list(fit=object,target_series=NULL,diagnostics=NULL,history=data.frame(),success=FALSE)
  save_state<-function()if(!is.null(checkpoint_file))mixedgp_checkpoint_write(mixedgp_pack_artifact(state),checkpoint_file)
  state$settings<-list(ess_min=ess_min,chunk_iter=chunk_iter,n_latent=n_latent,seed=seed,X=X,C=C,U=U)
  save_state()
  for(k in seq_len(.Machine$integer.max)) {
    state$target_series<-mixedgp_continuation_target_series(state$fit,X,C,U,n_latent,seed,parallel,n_cores)
    state$diagnostics<-mixedgp_summarize_diagnostic_series(state$target_series)
    state$success<-mixedgp_scientific_ess_pass(state$diagnostics,ess_min)
    state$settings<-list(ess_min=ess_min,chunk_iter=chunk_iter,n_latent=n_latent,seed=seed,X=X,C=C,U=U)
    state$history<-rbind(state$history,data.frame(iteration=state$fit$checkpoint$iteration,
      min_bulk_ess=min(state$diagnostics$ess_bulk),min_tail_ess=min(state$diagnostics$ess_tail),
      max_rhat=max(state$diagnostics$rhat),success=state$success))
    save_state()
    message("Scientific targets at iteration ",state$fit$checkpoint$iteration,": min bulk ESS=",
      signif(tail(state$history$min_bulk_ess,1),5),", min tail ESS=",signif(tail(state$history$min_tail_ess,1),5))
    if(state$success || k>max_extensions)break
    extra <- if(isTRUE(adaptive)) mixedgp_ess_extension_size(state$diagnostics,
      state$fit$checkpoint$iteration-state$fit$checkpoint$burn,ess_min,chunk_iter) else chunk_iter
    state$fit<-if(inherits(state$fit,"eivgp_fit"))continue_eivgp(state$fit,extra,parallel,n_cores) else
      if(mixedgp_fit_engine(state$fit)=="multivariate")mixedgp_continue_study2_fit(state$fit,extra,parallel,n_cores) else
      mixedgp_continue_study1_fit(state$fit,extra,parallel,n_cores)
    state$target_series<-NULL;state$diagnostics<-NULL;state$success<-FALSE
    save_state()
  }
  state
}

mixedgp_continuation_target_series <- function(fit,X,C,U,n_latent,seed,parallel,n_cores) {
  fun<-if(mixedgp_fit_engine(fit)=="univariate")mixedgp_study1_target_series else mixedgp_study2_target_series
  if(!isTRUE(parallel)||is.null(n_cores)||n_cores<=1L)
    return(fun(fit,X=X,C=C,U=U,n_latent=n_latent,seed=seed,max_draws_per_chain=Inf))
  nc<-length(fit$mcmc$samples_by_chain$logtheta)
  parts<-mixedgp_parallel_lapply(seq_len(nc),function(i) {
    one<-fit
    one$mcmc$samples_by_chain<-lapply(fit$mcmc$samples_by_chain,function(chains)chains[i])
    fun(one,X=X,C=C,U=U,n_latent=n_latent,seed=seed,max_draws_per_chain=Inf)
  },n_cores=min(nc,n_cores),mc.preschedule=FALSE)
  series<-setNames(lapply(names(parts[[1]]),function(name)lapply(parts,function(p)p[[name]][[1]])),names(parts[[1]]))
  mixedgp_mark_diagnostic_window(series,TRUE)
}

mixedgp_continue_study2_fit <- function(object, n_iter, parallel=FALSE, n_cores=1L) {
  if (inherits(object,"eivgp_fit")) return(continue_eivgp(object,n_iter,parallel,n_cores))
  cp <- object$checkpoint
  if (!identical(object$sampler_version,"0.3.1") || !identical(cp$version,2L) ||
      !identical(cp$sampler_version,"0.3.1") || !identical(cp$thin,1L) ||
      !identical(cp$control,object$control) || mixedgp_fit_engine(object)!="multivariate")
    stop("Unsupported or modified Study II continuation state.")
  n_iter <- mixedgp_as_integer_strict(n_iter,"n_iter",1L,1L)
  by <- object$mcmc$samples_by_chain
  expected <- cp$iteration-cp$burn
  nc <- length(cp$states)
  if (nc<1L || any(vapply(by$U,function(x)dim(x)[1L],integer(1))!=expected) || length(by$U)!=nc ||
      nrow(object$mcmc$chain_stats)!=nc || any(vapply(cp$states,function(s)is.null(s$rng)||is.null(s$state),logical(1))))
    stop("Incomplete retained history or terminal chain states.")
  cp$old_chains <- lapply(seq_len(nc),function(i) {
    out<-lapply(by,`[[`,i);names(out)<-paste0("samples_",names(out))
    out$stats<-object$mcmc$chain_stats[i,,drop=FALSE]
    out$initial_state<-object$mcmc$chain_initial[[i]];out
  })
  rng<-mixedgp_rng_state();on.exit(mixedgp_restore_rng_state(rng),add=TRUE)
  pr<-object$priors;pr$dictionary_type<-NULL
  d<-object$data
  result<-mixedgp_v030_fit(X=d$X,y=d$y,C=d$C_ord,
    U_obs=d$U_obs,m_vec=d$m_vec,measurement="probit",ident=d$ident,
    kernel=object$kernel$name,matern_nu=object$kernel$matern_nu,
    n_iter=cp$iteration+n_iter,burn=cp$burn,n_chains=nc,
    seed=as.integer(object$mcmc$chain_stats$seed[1]-10000L),
    parallel_chains=parallel,n_cores=n_cores,priors=pr,sampler_control=object$control,.resume=cp)
  updated<-mixedgp_v030_fit_object(result,d,"multivariate",object$kernel$name,object$kernel$matern_nu,cp)
  updated$diagnostics$summary$time_seconds<-object$diagnostics$summary$time_seconds+updated$diagnostics$summary$time_seconds
  updated$continuation_history<-rbind(object$continuation_history,data.frame(
    from_iteration=cp$iteration,to_iteration=cp$iteration+n_iter,added_draws_per_chain=n_iter))
  updated
}


mixedgp_validate_ess_request <- function(ess_min, panel=NULL) {
  if(!is.null(ess_min) && (!is.numeric(ess_min)||length(ess_min)!=1L||
      !is.finite(ess_min)||ess_min<=0))stop("ess_min must be NULL or one finite positive number.")
  if(!is.null(panel) && (!is.list(panel)||is.null(panel$X)||is.null(panel$C)))
    stop("ess_panel must contain X and C, with optional U, n_latent and seed.")
  invisible(NULL)
}

mixedgp_ess_extension_size <- function(tab, retained, ess_min, chunk_iter=5000L) {
  weakest<-suppressWarnings(min(c(tab$ess_bulk,tab$ess_tail)))
  if(!is.finite(weakest)||weakest<=0)return(as.integer(chunk_iter))
  estimate<-ceiling(retained*max(0,1.2*ess_min/weakest-1))
  as.integer(min(4*chunk_iter,max(chunk_iter,estimate)))
}

mixedgp_finish_ess <- function(fit, ess_min=NULL, ess_panel=NULL, parallel=FALSE,
                               n_cores=1L, checkpoint_file=NULL) {
  mixedgp_validate_ess_request(ess_min,ess_panel)
  if(is.null(ess_min))return(fit)
  if(is.null(ess_panel)) {
    d<-fit$data;one<-mixedgp_fit_engine(fit)=="univariate"
    X<-if(one)d$x_raw else d$X_raw;C<-if(one)matrix(d$c_ord,ncol=1L)else d$C_ord
    rows<-unique(as.integer(round(seq(1,nrow(X),length.out=min(5,nrow(X))))))
    ess_panel<-list(X=X[rows,,drop=FALSE],C=C[rows,,drop=FALSE])
  }
  z<-continue_eivgp_until_ess(fit,ess_panel$X,ess_panel$C,ess_panel$U,
    ess_min=ess_min,n_latent=if(is.null(ess_panel$n_latent))64L else ess_panel$n_latent,
    seed=if(is.null(ess_panel$seed))481517L else ess_panel$seed,
    parallel=parallel,n_cores=n_cores,checkpoint_file=checkpoint_file,adaptive=TRUE)
  if(!isTRUE(z$success))stop("Scientific ESS target was not reached; retained state is available.")
  fit<-z$fit
  fit$scientific_ess<-z[setdiff(names(z),'fit')]
  fit
}

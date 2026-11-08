# Disk-backed posterior histories. Only this small index is replaced at commit;
# immutable draw blocks and the original legacy fit survive interrupted writes.
mixedgp_disk_directory <- function(path) paste0(path, ".draws")

summary.eivgp_disk_fit <- function(object, ...) {
  cp <- object$checkpoint; nc <- length(cp$states)
  object$diagnostics <- list(summary = data.frame(sampler_strategy = object$sampler_strategy,
    n_chains = nc, n_iter = cp$iteration, burn = cp$burn, thin = cp$thin,
    saved_per_chain = cp$iteration - cp$burn,
    total_saved_draws = nc * (cp$iteration - cp$burn),
    parallel_backend = if (is.null(object$disk_history$parallel_backend)) 'serial' else object$disk_history$parallel_backend,
    parallel_cores = if (is.null(object$disk_history$parallel_cores)) 1L else object$disk_history$parallel_cores,
    time_seconds = object$disk_history$elapsed_seconds))
  result <- summary.eivgp_fit(object, ...)
  result$storage <- 'disk'
  result
}

print.eivgp_disk_fit <- function(x, ...) {
  cp <- x$checkpoint
  cat('<eivgp_disk_fit>\n', length(cp$states), ' chains; iteration ', cp$iteration,
      '; burn ', cp$burn, '; ', length(cp$states) * (cp$iteration - cp$burn),
      ' retained draws stored on disk.\nUse read_fit(path) for full-history analysis.\n', sep = '')
  invisible(x)
}

# A short call to the unchanged sampler performs the original warmup and at
# most one retained block. All subsequent sampling uses terminal states only.
mixedgp_fit_disk <- function(arguments, dots, path, chunk_iter, call, ess_min, ess_panel) {
  mixedgp_validate_named_dots(dots)
  if (!is.character(path) || length(path) != 1L || is.na(path) || !nzchar(path))
    stop("Disk storage requires one nonempty checkpoint path.")
  chunk_iter <- mixedgp_as_integer_strict(chunk_iter, "chunk_iter", 1L, 1L)
  target <- mixedgp_as_integer_strict(if (is.null(dots$n_iter)) 1750L else dots$n_iter, "n_iter", 1L, 1L)
  burn <- mixedgp_as_integer_strict(if (is.null(dots$burn)) 500L else dots$burn, "burn", 0L, 1L)
  if (target <= burn) stop("n_iter must exceed burn.")
  if (!is.null(ess_min) && !is.null(dots$n_chains) && dots$n_chains < 4L)
    stop("Scientific ESS requires at least four chains.")
  path <- path.expand(path)
  if (file.exists(path) || dir.exists(mixedgp_disk_directory(path)))
    stop("Checkpoint path already exists; use continue_eivgp_disk() to resume it: ", path)
  rng <- mixedgp_rng_state(); on.exit(mixedgp_restore_rng_state(rng), add = TRUE)
  first <- as.integer(min(as.double(target), as.double(burn) + chunk_iter))
  dots$n_iter <- first
  prefix <- do.call(fit_eivgp, c(arguments, dots))
  prefix$call <- call
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  save_fit(prefix, path)
  mixedgp_retire_saved_checkpoints(prefix, path)
  rm(prefix); gc(FALSE)
  path <- normalizePath(path)
  index <- continue_eivgp_disk(path, target - first, chunk_iter = chunk_iter,
    parallel = arguments$parallel, n_cores = arguments$n_cores,
    verbose = isTRUE(dots$verbose))
  if (!is.null(ess_min)) index <- mixedgp_finish_ess_disk(path, ess_min, ess_panel,
    arguments$parallel, arguments$n_cores, chunk_iter)
  index
}

# Use the established scientific-ESS check and extension policy. Diagnostic
# materialization is separate from sampling, and is released before forking.
mixedgp_finish_ess_disk <- function(path, ess_min, panel = NULL, parallel = FALSE,
                                    n_cores = 1L, chunk_iter = 1000L) {
  index <- mixedgp_disk_index(path)
  if (length(index$checkpoint$states) < 4L) stop("Scientific ESS requires at least four chains.")
  if (is.null(panel)) {
    d <- index$data; one <- mixedgp_fit_engine(index) == 'univariate'
    X <- if (one) d$x_raw else d$X_raw
    C <- if (one) matrix(d$c_ord, ncol = 1L) else d$C_ord
    rows <- unique(as.integer(round(seq(1, nrow(X), length.out = min(5, nrow(X))))))
    panel <- list(X = X[rows, , drop = FALSE], C = C[rows, , drop = FALSE])
  }
  history <- data.frame()
  repeat {
    fit <- read_fit(path, fields = setdiff(index$disk_history$fields, 'S'))
    state <- continue_eivgp_until_ess(fit, panel$X, panel$C, panel$U, ess_min = ess_min,
      n_latent = if (is.null(panel$n_latent)) 64L else panel$n_latent,
      seed = if (is.null(panel$seed)) 481517L else panel$seed,
      parallel = parallel, n_cores = n_cores, max_extensions = 0L, adaptive = TRUE)
    history <- rbind(history, state$history)
    if (isTRUE(state$success)) {
      result <- state[setdiff(names(state), 'fit')]; result$history <- history
      # The in-memory ESS loop removes/re-adds these fields after extension.
      if (nrow(history) > 1L)
        result <- result[c('history', 'success', 'settings', 'target_series', 'diagnostics')]
      return(mixedgp_disk_save_scientific(index, result, path))
    }
    extra <- mixedgp_ess_extension_size(state$diagnostics,
      index$checkpoint$iteration - index$checkpoint$burn, ess_min, 5000L)
    rm(fit, state); gc(FALSE)
    index <- continue_eivgp_disk(path, extra, chunk_iter = chunk_iter,
                                parallel = parallel, n_cores = n_cores)
  }
}

mixedgp_disk_save_scientific <- function(index, result, path) {
  name <- sprintf('scientific-ess-%09d.rds', index$checkpoint$iteration)
  destination <- file.path(mixedgp_disk_directory(path), name)
  mixedgp_atomic_rds(result, destination)
  index$scientific_ess <- result[setdiff(names(result), 'target_series')]
  index$disk_history$scientific_ess <- list(file = name, iteration = index$checkpoint$iteration,
    md5 = unname(tools::md5sum(destination)))
  mixedgp_disk_commit(index, path)
}

# Legacy conversion is deliberately isolated and sequential: historical arrays
# must not be inherited by a later pool of application fitting workers.
mixedgp_prepare_disk_fits <- function(paths, chunk_iter = 1000L) {
  for (path in unique(paths)) if (file.exists(path) &&
      !file.exists(file.path(mixedgp_disk_directory(path), 'index.rds'))) {
    cat('Preparing disk history: ', path, '\n', sep = ''); flush.console()
    code <- paste0('a<-commandArgs(TRUE);.libPaths(strsplit(a[3],.Platform$path.sep,fixed=TRUE)[[1]]);',
      'eivGP::continue_eivgp_disk(a[1],0L,chunk_iter=as.integer(a[2]),parallel=FALSE,n_cores=1L)')
    status <- system2(file.path(R.home('bin'), 'Rscript'),
      shQuote(c('--vanilla', '-e', code, normalizePath(path), as.character(chunk_iter),
                paste(unique(c(dirname(getNamespaceInfo(asNamespace('eivGP'), 'path')),
                               .libPaths())), collapse = .Platform$path.sep))),
      env = c('OPENBLAS_NUM_THREADS=1', 'OMP_NUM_THREADS=1', 'MKL_NUM_THREADS=1'))
    if (status != 0L) stop('Disk-history preparation failed; original checkpoint retained: ', path)
  }
  invisible(paths)
}

mixedgp_disk_checkpoint <- function(path, chunk_iter = 1000L) {
  if (file.exists(file.path(mixedgp_disk_directory(path), 'index.rds')))
    mixedgp_disk_index(path) else continue_eivgp_disk(path, 0L, chunk_iter = chunk_iter)
}

mixedgp_application_disk_ess <- function(path, ess_min, panel, cores = 1L,
                                         chunk_iter = 1000L, trace_path = NULL, label = NULL) {
  ess_min <- mixedgp_parse_ess_min(ess_min)
  index <- mixedgp_disk_index(path)
  settings <- index$scientific_ess$settings
  cached <- !is.null(index$disk_history$scientific_ess) &&
    identical(index$disk_history$scientific_ess$iteration, index$checkpoint$iteration) &&
    identical(settings$X, panel$X) && identical(settings$C, panel$C) && identical(settings$U, panel$U) &&
    identical(settings$n_latent, panel$n_latent) && identical(settings$seed, panel$seed) &&
    identical(settings$ess_min, ess_min) && (is.null(ess_min) || isTRUE(index$scientific_ess$success))
  if (!cached && !is.null(ess_min)) {
    mixedgp_finish_ess_disk(path, ess_min, panel, cores > 1L, cores, chunk_iter)
    index <- mixedgp_disk_index(path)
  }
  fit <- read_fit(path, fields = setdiff(index$disk_history$fields, 'S'))
  if (!cached && is.null(ess_min)) {
    fit <- mixedgp_application_ess(fit, NULL, panel, NULL, cores)
    mixedgp_disk_save_scientific(index, fit$scientific_ess, path)
  }
  if (!is.null(trace_path)) mixedgp_save_scientific_trace(fit, trace_path, label)
  mixedgp_atomic_rds(list(iteration = fit$checkpoint$iteration, scientific_ess = fit$scientific_ess),
    file.path(dirname(path), 'scientific-diagnostics.rds'))
  invisible(fit$checkpoint$iteration)
}

mixedgp_disk_analysis_fit <- function(path) {
  if (!file.exists(file.path(mixedgp_disk_directory(path), 'index.rds')))
    return(mixedgp_read_artifact(path))
  index <- mixedgp_disk_index(path)
  read_fit(path, fields = setdiff(index$disk_history$fields, 'S'))
}

# Copy an archive without materializing its draw history, e.g. Ocean's separate
# continuation output. Immutable blocks may be hard-linked on the same device.
mixedgp_disk_copy <- function(source, destination) {
  if (file.exists(destination)) stop('Continuation destination already exists: ', destination)
  directory <- mixedgp_disk_directory(destination)
  if (dir.exists(directory)) {
    # Recover a completed directory publication interrupted before its pointer.
    index <- mixedgp_disk_index(destination, verify = TRUE)
    mixedgp_atomic_rds(index, destination)
    return(index)
  }
  dir.create(dirname(destination), recursive = TRUE, showWarnings = FALSE)
  if (!file.exists(file.path(mixedgp_disk_directory(source), 'index.rds'))) {
    if (!file.link(source, destination) && !file.copy(source, destination))
      stop('Cannot copy legacy checkpoint: ', source)
    index <- continue_eivgp_disk(destination, 0L)
    index$disk_history$copied_from_iteration <- index$checkpoint$iteration
    return(mixedgp_disk_commit(index, destination))
  }
  index <- mixedgp_disk_index(source, verify = TRUE)
  index$disk_history$copied_from_iteration <- index$checkpoint$iteration
  index$disk_history$signature <- mixedgp_disk_signature(index)
  staging <- tempfile('.disk-clone-', tmpdir = dirname(destination))
  dir.create(staging)
  on.exit(unlink(staging, recursive = TRUE), add = TRUE)
  files <- list.files(mixedgp_disk_directory(source), full.names = TRUE)
  for (file in files[basename(files) != 'index.rds']) {
    target <- file.path(staging, basename(file))
    if (!file.link(file, target) && !file.copy(file, target)) stop('Cannot copy disk archive: ', file)
  }
  mixedgp_atomic_rds(index, file.path(staging, 'index.rds'))
  if (!file.rename(staging, directory)) stop('Cannot publish continuation archive: ', directory)
  mixedgp_atomic_rds(index, destination)
  index
}

mixedgp_saved_iteration <- function(path) {
  if (file.exists(file.path(mixedgp_disk_directory(path), "index.rds")))
    return(mixedgp_disk_index(path)$checkpoint$iteration)
  mixedgp_read_artifact(path)$checkpoint$iteration
}

mixedgp_disk_signature <- function(index) {
  index$disk_history$signature <- NULL
  mixedgp_object_fingerprint(index)
}

mixedgp_disk_index <- function(path, verify = FALSE) {
  directory <- mixedgp_disk_directory(path)
  index <- readRDS(file.path(directory, "index.rds"))
  if (!inherits(index, "eivgp_disk_fit") ||
      !identical(index$disk_history$schema, "eivgp-disk-v1") ||
      !identical(index$disk_history$signature, mixedgp_disk_signature(index)))
    stop("Invalid or modified disk continuation index: ", path)
  cp <- index$checkpoint
  if (!identical(cp$version, 2L) || !identical(cp$sampler_version, "0.3.1") ||
      !identical(cp$control, index$control) || !identical(cp$thin, 1L))
    stop("Incompatible disk continuation state.")
  blocks <- index$disk_history$blocks
  for (chain in seq_along(cp$states)) {
    rows <- blocks[blocks$chain == chain, , drop = FALSE]
    rows <- rows[order(rows$from), , drop = FALSE]
    if (!nrow(rows) || rows$from[1L] != cp$burn + 1L || tail(rows$to, 1L) != cp$iteration ||
        any(rows$to < rows$from) || any(tail(rows$from, -1L) != head(rows$to, -1L) + 1L))
      stop("Disk history has missing or overlapping draw blocks.")
  }
  all_blocks <- rbind(blocks, index$disk_history$hmc_blocks)
  files <- file.path(directory, all_blocks$file)
  if (any(!file.exists(files)) || any(file.info(files)$size != all_blocks$bytes))
    stop("Disk history block is missing or changed.")
  if (verify && !identical(unname(tools::md5sum(files)), all_blocks$md5))
    stop("Disk history block checksum differs.")
  index
}

mixedgp_disk_commit <- function(index, path) {
  index$disk_history$signature <- mixedgp_disk_signature(index)
  # The archive index is authoritative if interrupted between these two writes.
  mixedgp_atomic_rds(index, file.path(mixedgp_disk_directory(path), "index.rds"))
  mixedgp_atomic_rds(index, path)
  index
}

mixedgp_disk_slice <- function(x, rows) {
  if (is.null(x)) return(NULL)
  if (length(dim(x)) == 3L) return(x[rows, , , drop = FALSE])
  if (is.matrix(x)) return(x[rows, , drop = FALSE])
  x[rows]
}

mixedgp_disk_write_block <- function(draws, chain, from, to, directory, prefix = "chain") {
  name <- sprintf("%s%02d_%09d_%09d.rds", prefix, chain, from, to)
  path <- file.path(directory, name)
  mixedgp_atomic_rds(draws, path)
  data.frame(chain = as.integer(chain), from = as.integer(from), to = as.integer(to),
             file = name, bytes = unname(file.info(path)$size), md5 = unname(tools::md5sum(path)),
             stringsAsFactors = FALSE)
}

mixedgp_disk_import <- function(path, chunk_iter) {
  if (file.exists(file.path(mixedgp_disk_directory(path), "index.rds")))
    return(mixedgp_disk_index(path, verify = TRUE))
  fit <- mixedgp_read_artifact(path)
  if (!inherits(fit, "eivgp_fit") || !identical(fit$checkpoint$version, 2L) ||
      !identical(fit$checkpoint$sampler_version, "0.3.1") || !identical(fit$checkpoint$thin, 1L) ||
      !identical(fit$checkpoint$control, fit$control) ||
      !identical(fit$checkpoint$signature, mixedgp_v030_checkpoint_hash(fit)))
    stop("Unsupported or modified legacy fit; disk import refused.")
  cp <- fit$checkpoint
  retained <- cp$iteration - cp$burn
  fields <- names(fit$mcmc$samples_by_chain)
  if (retained < 1L || !length(cp$states) ||
      any(vapply(fit$mcmc$samples_by_chain$logtheta, nrow, integer(1)) != retained))
    stop("Incomplete legacy retained history.")
  directory <- mixedgp_disk_directory(path)
  dir.create(directory, recursive = TRUE, showWarnings = FALSE)
  original <- file.path(directory, "legacy-fit.rds")
  if (file.exists(original)) {
    if (!identical(unname(tools::md5sum(original)), unname(tools::md5sum(path))))
      stop("Legacy fit changed during an unfinished disk import.")
  } else if (!file.link(path, original) && !file.copy(path, original))
    stop("Cannot preserve the original legacy fit.")
  blocks <- list()
  for (chain in seq_along(cp$states)) for (from in seq.int(1L, retained, by = chunk_iter)) {
    to <- min(retained, from + chunk_iter - 1L)
    draws <- lapply(fit$mcmc$samples_by_chain, function(chains)
      mixedgp_disk_slice(chains[[chain]], seq.int(from, to)))
    blocks[[length(blocks) + 1L]] <- mixedgp_disk_write_block(draws, chain,
      cp$burn + from, cp$burn + to, directory)
  }
  traces <- list()
  for (chain in seq_along(cp$states)) {
    history <- fit$checkpoint$states[[chain]]$state$hmc$history
    if (!is.null(history)) {
      traces[[length(traces) + 1L]] <- mixedgp_disk_write_block(history, chain, 1L,
        cp$iteration, directory, prefix = "hmc")
      fit$checkpoint$states[[chain]]$state$hmc$history <- history[FALSE, , drop = FALSE]
    }
  }
  elapsed <- fit$diagnostics$summary$time_seconds
  parallel_backend <- fit$diagnostics$summary$parallel_backend
  parallel_cores <- fit$diagnostics$summary$parallel_cores
  fit$mcmc <- fit$mcmc[c("chain_stats", "chain_initial")]
  fit$diagnostics <- NULL
  fit$scientific_ess <- NULL # Recomputed for the completed history during analysis.
  fit$checkpoint$signature <- NULL
  fit$disk_history <- list(schema = "eivgp-disk-v1", fields = fields,
    blocks = do.call(rbind, blocks), hmc_blocks = if (length(traces)) do.call(rbind, traces) else
      do.call(rbind, blocks)[FALSE, ], base_iteration = cp$iteration,
    legacy_md5 = unname(tools::md5sum(original)), elapsed_seconds = elapsed,
    parallel_backend = parallel_backend, parallel_cores = parallel_cores)
  class(fit) <- unique(c("eivgp_disk_fit", class(fit)))
  mixedgp_disk_commit(fit, path)
}

mixedgp_disk_add_stats <- function(old, new) {
  counters <- setdiff(names(old), c("chain", "seed", "mean_abs_angle"))
  for (name in counters) new[[name]] <- old[[name]] + new[[name]]
  new$mean_abs_angle <- new$abs_angle_sum / pmax(1, new$angle_count)
  new
}

#' Continue MCMC with immutable draw blocks on disk
#'
#' @param path Path to an existing public fit saved with save_fit().
#' @param n_iter Additional iterations per chain. Zero imports a legacy fit only.
#' @param chunk_iter Maximum new draws per chain held in RAM (default 1000).
#' @param parallel,n_cores Chain execution controls; sampler tuning is unchanged.
#' @param verbose Print checkpoint progress.
#' @param progress Optional function receiving the small committed index.
#' @return Invisibly, the small disk-backed fit index. The checkpoint at path is
#'   updated atomically; read_fit(path) reconstructs its full history on demand.
#' @details The original checkpoint is preserved in path.draws/legacy-fit.rds.
#'   New sampling retains only complete terminal states and bounded draw buffers;
#'   it uses the same transition function, RNG states, absolute iteration numbers,
#'   warmup exclusion, priors and controls. No draws are thinned or discarded.
#'   An interrupted, uncommitted block is recomputed from the last committed RNG
#'   state. Full history reads and diagnostics are explicit analysis operations.
#'   Only one writer may own a fit; application executors supply that lock.
#' @export
continue_eivgp_disk <- function(path, n_iter, chunk_iter = 1000L, parallel = TRUE,
                               n_cores = NULL, verbose = FALSE, progress = NULL) {
  path <- normalizePath(path, mustWork = TRUE)
  n_iter <- mixedgp_as_integer_strict(n_iter, "n_iter", 0L, 1L)
  chunk_iter <- mixedgp_as_integer_strict(chunk_iter, "chunk_iter", 1L, 1L)
  parallel <- mixedgp_validate_flag(parallel, "parallel")
  verbose <- mixedgp_validate_flag(verbose, "verbose")
  if (!is.null(progress) && !is.function(progress)) stop("progress must be a function.")
  rng <- mixedgp_rng_state(); on.exit(mixedgp_restore_rng_state(rng), add = TRUE)
  old <- options(eivGP.checkpoint_root = NULL, eivGP.checkpoint_iterations = Inf,
                 eivGP.checkpoint_seconds = Inf)
  on.exit(options(old), add = TRUE)
  index <- mixedgp_disk_import(path, chunk_iter)
  gc(FALSE) # Release the one-time legacy import before any chain workers fork.
  target <- mixedgp_as_integer_strict(as.double(index$checkpoint$iteration) + n_iter,
                                      "target iterations", 1L, 1L)
  directory <- mixedgp_disk_directory(path)
  while (index$checkpoint$iteration < target) {
    cp <- index$checkpoint
    before <- cp$iteration
    to <- as.integer(min(target, (floor(before / chunk_iter) + 1L) * chunk_iter))
    cp$signature <- NULL; cp$old_chains <- NULL; cp$state_only <- TRUE
    arguments <- cp$arguments
    arguments$n_iter <- to; arguments$parallel_chains <- parallel
    arguments$n_cores <- n_cores; arguments$verbose <- verbose
    sampler <- if (mixedgp_fit_engine(index) == "univariate") fit_eivgp_1d_v030 else fit_eivgp_ordprobit_v030
    result <- do.call(sampler, c(arguments, list(.resume = cp)))
    blocks <- vector("list", length(result$chains)); traces <- list()
    for (chain in seq_along(result$chains)) {
      draws <- setNames(lapply(index$disk_history$fields, function(field)
        result$chains[[chain]][[paste0("samples_", field)]]), index$disk_history$fields)
      blocks[[chain]] <- mixedgp_disk_write_block(draws, chain, before + 1L, to, directory)
      history <- result$chains[[chain]]$checkpoint$state$hmc$history
      if (!is.null(history)) {
        traces[[length(traces) + 1L]] <- mixedgp_disk_write_block(history, chain,
          before + 1L, to, directory, prefix = "hmc")
        result$chains[[chain]]$checkpoint$state$hmc$history <- history[FALSE, , drop = FALSE]
      }
    }
    index$disk_history$blocks <- rbind(index$disk_history$blocks, do.call(rbind, blocks))
    if (length(traces)) index$disk_history$hmc_blocks <-
      rbind(index$disk_history$hmc_blocks, do.call(rbind, traces))
    index$disk_history$elapsed_seconds <- index$disk_history$elapsed_seconds + result$time_seconds
    index$checkpoint$states <- lapply(result$chains, function(chain) chain$checkpoint)
    index$checkpoint$iteration <- to
    index$checkpoint$arguments <- arguments
    index$scientific_ess <- NULL
    index$disk_history$scientific_ess <- NULL
    index$mcmc$chain_stats <- mixedgp_disk_add_stats(index$mcmc$chain_stats,
      do.call(rbind, lapply(result$chains, function(chain) chain$stats)))
    index$interface$parallel <- parallel
    index$interface$n_cores <- result$parallel_cores
    index$disk_history$parallel_backend <- result$parallel_backend
    index$disk_history$parallel_cores <- result$parallel_cores
    index$continuation_history <- rbind(index$continuation_history, data.frame(
      from_iteration = before, to_iteration = to, added_draws_per_chain = to - before))
    rm(result, draws, blocks, cp, traces, history); gc(FALSE)
    index <- mixedgp_disk_commit(index, path)
    if (verbose) message("Committed disk history through iteration ", to, "/", target)
    if (!is.null(progress)) progress(index)
  }
  invisible(index)
}

# Full-history reconstruction is intentionally separate from the sampling path.
mixedgp_disk_materialize <- function(path, draw_ids = NULL, fields = NULL) {
  index <- mixedgp_disk_index(path, verify = TRUE)
  cp <- index$checkpoint
  full_fields <- index$disk_history$fields
  if (is.null(fields)) fields <- index$disk_history$fields
  if (!all(fields %in% index$disk_history$fields) || anyDuplicated(fields) ||
      !all(setdiff(index$disk_history$fields, 'S') %in% fields))
    stop("Only auxiliary S draws may be omitted from a diagnostic fit.")
  nc <- length(cp$states); retained <- cp$iteration - cp$burn
  if (is.null(draw_ids)) ids <- seq_len(retained * nc) else {
    ids <- mixedgp_as_integer_strict(draw_ids, "draw_ids", 1L)
    if (anyDuplicated(ids) || any(ids > retained * nc)) stop("Invalid or repeated draw IDs.")
    ids <- sort(ids)
    if (any(tabulate((ids - 1L) %/% retained + 1L, nbins = nc) < 2L))
      stop("Select at least two draws per chain for fit diagnostics.")
  }
  wanted <- lapply(seq_len(nc), function(chain) ids[ids > (chain - 1L) * retained &
    ids <= chain * retained] - (chain - 1L) * retained)
  blocks <- index$disk_history$blocks
  chains <- vector("list", nc)
  for (chain in seq_len(nc)) {
    traces <- index$disk_history$hmc_blocks
    traces <- traces[traces$chain == chain, , drop = FALSE]
    if (nrow(traces)) {
      traces <- traces[order(traces$from), , drop = FALSE]
      if (traces$from[1L] != 1L || tail(traces$to, 1L) != cp$iteration ||
          any(tail(traces$from, -1L) != head(traces$to, -1L) + 1L)) stop("Incomplete HMC trace history.")
      cp$states[[chain]]$state$hmc$history <- do.call(rbind, lapply(traces$file, function(file)
        readRDS(file.path(mixedgp_disk_directory(path), file))))
    }
    output <- NULL
    rows <- blocks[blocks$chain == chain, , drop = FALSE]
    for (j in seq_len(nrow(rows))) {
      row <- rows[j, ]
      local <- wanted[[chain]] + cp$burn
      keep <- which(local >= row$from & local <= row$to)
      if (!length(keep)) next
      block <- readRDS(file.path(mixedgp_disk_directory(path), row$file))[fields]
      if (is.null(output)) output <- lapply(block, function(x) {
        if (is.null(x)) return(NULL)
        dims <- dim(x); n <- length(wanted[[chain]])
        if (is.null(dims)) return(vector(typeof(x), n))
        dims[1L] <- n
        value <- array(vector(typeof(x), prod(dims)), dims)
        dn <- dimnames(x)
        if (!is.null(dn)) { dn[1L] <- list(NULL); dimnames(value) <- dn }
        value
      })
      for (field in names(output)) {
        x <- block[[field]]
        if (is.null(x)) next
        source <- local[keep] - row$from + 1L
        if (length(dim(x)) == 3L) output[[field]][keep, , ] <- x[source, , , drop = FALSE]
        else if (is.matrix(x)) output[[field]][keep, ] <- x[source, , drop = FALSE]
        else output[[field]][keep] <- x[source]
      }
    }
    names(output) <- paste0("samples_", names(output))
    if (mixedgp_fit_engine(index) == 'univariate') {
      u <- output$samples_u
      output$samples_U <- array(u, c(nrow(u), ncol(u), 1L))
    }
    output$stats <- index$mcmc$chain_stats[chain, , drop = FALSE]
    output$initial_state <- index$mcmc$chain_initial[[chain]]
    output$checkpoint <- cp$states[[chain]]
    chains[[chain]] <- output
  }
  result <- list(chains = chains, priors = index$priors, control = index$control,
    iteration = cp$iteration, burn = cp$burn, sampler_version = "0.3.1",
    time_seconds = index$disk_history$elapsed_seconds,
    parallel_backend = if (is.null(index$disk_history$parallel_backend)) "serial" else index$disk_history$parallel_backend,
    parallel_cores = if (is.null(index$disk_history$parallel_cores)) 1L else index$disk_history$parallel_cores)
  built <- mixedgp_v030_fit_object(result, index$data, mixedgp_fit_engine(index),
                                   index$kernel$name, index$kernel$matern_nu)
  index$mcmc <- built$mcmc; index$diagnostics <- built$diagnostics
  scientific <- index$disk_history$scientific_ess
  if (!is.null(scientific) && identical(scientific$iteration, cp$iteration)) {
    file <- file.path(mixedgp_disk_directory(path), scientific$file)
    if (!identical(unname(tools::md5sum(file)), scientific$md5)) stop("Scientific ESS archive changed.")
    index$scientific_ess <- readRDS(file)
  }
  index$checkpoint <- cp
  index$disk_history <- NULL
  class(index) <- setdiff(class(index), "eivgp_disk_fit")
  if (is.null(draw_ids) && identical(fields, full_fields))
    index$checkpoint$signature <- mixedgp_v030_checkpoint_hash(index)
  else {
    index$mcmc$mcmc_draw_info$draw_within_chain <- unlist(wanted, use.names = FALSE)
    index$draw_selection <- if (is.null(draw_ids)) NULL else ids
    index$draw_selection_fields <- fields
    index$checkpoint$signature <- NULL
    attr(index, "full_retained_history") <- is.null(draw_ids)
  }
  index
}

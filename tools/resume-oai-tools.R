# Shared command-line recovery for the frozen OAI method groups.
resume_oai_process <- function(pid) {
  path <- file.path('/proc', as.character(pid), 'stat')
  if (!file.exists(path)) return(NULL)
  fields <- tryCatch(strsplit(sub('^.*\\) ', '', readLines(path, warn = FALSE)), ' +')[[1L]],
                     error = function(e) NULL)
  if (length(fields) < 20L) return(NULL)
  list(pid = as.integer(pid), state = fields[1L], group = as.integer(fields[3L]),
       session = as.integer(fields[4L]), start = fields[20L])
}

resume_oai_record_session <- function(path) {
  owner <- resume_oai_process(Sys.getpid())
  if (is.null(owner) || owner$group != owner$pid || owner$session != owner$pid)
    stop('Resume worker must start in its own process session.', call. = FALSE)
  launcher_path <- paste0(path, '.launcher')
  if (file.exists(launcher_path)) owner$launcher <- readRDS(launcher_path)
  tmp <- tempfile(tmpdir = dirname(path))
  on.exit(unlink(tmp), add = TRUE)
  saveRDS(owner, tmp)
  if (!file.rename(tmp, path)) stop('Cannot record isolated resume session.', call. = FALSE)
}

resume_oai_cleanup_session <- function(path) {
  if (!file.exists(path)) return(invisible(integer()))
  owner <- readRDS(path)
  leader <- resume_oai_process(owner$pid)
  # Never signal a recycled PID or the launcher's own process group.
  if (owner$group != owner$pid || owner$session != owner$pid ||
      owner$pid == resume_oai_process(Sys.getpid())$group ||
      (!is.null(leader) && !identical(leader$start, owner$start)))
    stop('Resume session identity changed; worker cleanup refused.', call. = FALSE)
  members <- function() {
    ids <- list.files('/proc', '^[0-9]+$')
    info <- lapply(ids, resume_oai_process)
    vapply(Filter(function(x) !is.null(x) && x$state != 'Z' &&
      x$group == owner$group && x$session == owner$session, info), `[[`, integer(1), 'pid')
  }
  active <- members()
  if (!length(active)) return(invisible(active))
  tools::pskill(active, tools::SIGTERM)
  for (i in seq_len(20L)) {
    Sys.sleep(0.05)
    if (!length(members())) break
  }
  remaining <- members()
  if (length(remaining)) {
    tools::pskill(remaining, tools::SIGKILL)
    for (i in seq_len(20L)) {
      Sys.sleep(0.05)
      if (!length(members())) break
    }
  }
  cat('Stopped ', length(active), ' remaining workers from this resume session.\n', sep = '')
  invisible(active)
}

resume_oai_launch_worker <- function(script, args) {
  setsid <- Sys.which('setsid')
  if (!nzchar(setsid)) stop('Linux setsid is required for reliable worker cleanup.', call. = FALSE)
  path <- tempfile('oai-resume-session-', fileext = '.rds')
  status_path <- paste0(path, '.status'); pending_status <- paste0(status_path, '.tmp')
  runner <- paste0(path, '.sh'); launcher_path <- paste0(path, '.launcher')
  saveRDS(resume_oai_process(Sys.getpid()), launcher_path)
  previous <- Sys.getenv('EIVGP_OAI_SESSION_FILE', unset = NA_character_)
  Sys.setenv(EIVGP_OAI_SESSION_FILE = path)
  on.exit({
    resume_oai_cleanup_session(path)
    unlink(c(path, status_path, pending_status, runner, launcher_path))
    if (is.na(previous)) Sys.unsetenv('EIVGP_OAI_SESSION_FILE') else
      Sys.setenv(EIVGP_OAI_SESSION_FILE = previous)
  }, add = TRUE)
  # A synchronous system2() ignores SIGINT while waiting. Poll an asynchronous
  # shell instead, so Ctrl+C unwinds this function and runs the session cleanup.
  command <- paste(shQuote(c(setsid, '--wait', file.path(R.home('bin'), 'Rscript'),
                             '--vanilla', script, args)), collapse = ' ')
  writeLines(c('#!/bin/sh', command, 'resume_status=$?',
    paste("printf '%s\\n' \"$resume_status\" >", shQuote(pending_status)),
    paste('mv -f --', shQuote(pending_status), shQuote(status_path))), runner)
  launched <- system2('/bin/sh', shQuote(runner), wait = FALSE)
  if (launched != 0L) stop('Cannot launch the isolated resume worker.', call. = FALSE)
  started <- proc.time()[['elapsed']]; dead_polls <- 0L
  while (!file.exists(status_path)) {
    if (file.exists(path)) {
      owner <- readRDS(path)
      leader <- resume_oai_process(owner$pid)
      if (is.null(leader) || leader$state %in% c('Z', 'X')) dead_polls <- dead_polls + 1L
      if (dead_polls > 40L) stop('Resume worker exited without recording its status.', call. = FALSE)
    } else if (proc.time()[['elapsed']] - started > 60)
      stop('Resume worker did not initialize within 60 seconds.', call. = FALSE)
    Sys.sleep(0.05)
  }
  status <- suppressWarnings(as.integer(readLines(status_path, warn = FALSE)))
  if (length(status) != 1L || is.na(status)) stop('Invalid resume exit status.', call. = FALSE)
  status
}

resume_oai_process_output <- function(pid) {
  # Resolve only our known resume entrypoints; unknown executors remain guarded.
  bytes <- tryCatch(readBin(file.path('/proc', as.character(pid), 'cmdline'),
                            'raw', n = 1048576L), error = function(e) raw())
  ends <- which(bytes == as.raw(0))
  if (!length(ends)) return(NULL)
  starts <- c(1L, head(ends, -1L) + 1L)
  args <- vapply(seq_along(ends), function(i) if (starts[i] < ends[i])
    rawToChar(bytes[seq.int(starts[i], ends[i] - 1L)]) else '', character(1))
  script <- sub('^--file=', '', grep('^--file=', args, value = TRUE))
  if (length(script) != 1L || !basename(script) %in%
      c('resume-oai-eivGP.r', 'resume-oai-competitors.r')) return(NULL)
  cwd <- Sys.readlink(file.path('/proc', as.character(pid), 'cwd'))
  absolute <- function(path) {
    path <- path.expand(path)
    if (!startsWith(path, '/')) {
      if (is.na(cwd) || !nzchar(cwd)) return(NULL)
      path <- file.path(cwd, path)
    }
    if (!dir.exists(path)) return(NULL)
    normalizePath(path)
  }
  output <- sub('^--output=', '', grep('^--output=', args, value = TRUE))
  if (length(output) > 1L) return(NULL)
  if (length(output)) return(absolute(output))
  if (!startsWith(script, '/')) {
    if (is.na(cwd) || !nzchar(cwd)) return(NULL)
    script <- file.path(cwd, script)
  }
  absolute(file.path(dirname(script), '..', 'results',
    if (basename(script) == 'resume-oai-eivGP.r') 'oai' else 'oai_competitors'))
}

resume_oai_recover_locks <- function(out, group) {
  fail <- function(...) stop(..., call. = FALSE)
  if (!dir.exists('/proc')) fail('Stale-lock verification requires Linux /proc.')
  legacy <- file.path(out, c('.application.lock', '.adni-application.lock', '.oai-application.lock'))
  modern <- file.path(out, '_meta/locks', paste0(group, '.lock'))
  locks <- c(legacy[dir.exists(legacy)], modern[dir.exists(modern)])
  live <- function(pid) {
    if (!is.finite(pid) || pid <= 0) return(FALSE)
    process <- resume_oai_process(pid)
    # A killed worker can remain in /proc until its parent reaps it. It can
    # no longer execute or write checkpoints. Unreadable entries stay guarded.
    if (is.null(process)) return(dir.exists(file.path('/proc', as.character(pid))))
    !process$state %in% c('Z', 'X', 'x')
  }
  recorded <- integer()
  owners <- list()
  for (lock in locks) {
    path <- file.path(lock, if (identical(lock, modern)) 'owner.json' else 'pid')
    snapshot <- readLines(path, warn = FALSE)
    if (identical(lock, modern)) {
      owner <- tryCatch(jsonlite::fromJSON(paste(snapshot, collapse = '\n')), error = function(e) NULL)
      if (is.null(owner) || length(owner$host) != 1L || is.na(owner$host) ||
          !identical(owner$host, unname(Sys.info()[['nodename']])))
        fail('Cannot verify executor lock host: ', lock)
      if (length(owner$token) != 1L || is.na(owner$token) || !nzchar(owner$token))
        fail('Cannot verify executor lock token: ', lock)
      pid <- suppressWarnings(as.numeric(owner$pid))
    } else pid <- suppressWarnings(as.numeric(snapshot))
    if (length(pid) != 1L || !is.finite(pid) || pid < 1 || pid != floor(pid) || pid > .Machine$integer.max)
      fail('Cannot verify lock owner PID: ', lock)
    recorded <- c(recorded, pid)
    owners[[lock]] <- list(path = path, snapshot = snapshot, pid = pid)
  }
  # A dead coordinator can leave orphaned workers. Check its phase records and
  # matching R processes as well; the other method group may run independently.
  files <- list.files(file.path(out, 'logs'), '^fit_status[.]csv$', recursive = TRUE, full.names = TRUE)
  if (!any(dir.exists(legacy))) files <- files[grepl(paste0('/logs/', group, '_'), files, fixed = TRUE)]
  for (f in files) {
    tab <- read.csv(f, stringsAsFactors = FALSE)
    recorded <- c(recorded, tab$worker_pid[tab$state == 'RUNNING'])
  }
  active <- unique(recorded[vapply(recorded, live, logical(1))])
  if (length(active)) fail('Recorded OAI process IDs still exist: ', paste(active, collapse = ', '))
  ancestors <- Sys.getpid()
  repeat {
    status <- readLines(file.path('/proc', tail(ancestors, 1L), 'status'), warn = FALSE)
    parent <- as.integer(sub('^PPid:\\s*', '', grep('^PPid:', status, value = TRUE)))
    if (length(parent) != 1L || is.na(parent) || parent <= 0L || parent %in% ancestors) break
    ancestors <- c(ancestors, parent)
  }
  # Asynchronous system2() reparents its shell, so the verified launcher need
  # not be an ancestor of the isolated worker in Linux's current process tree.
  session_file <- Sys.getenv('EIVGP_OAI_SESSION_FILE')
  if (nzchar(session_file) && file.exists(session_file)) {
    launcher <- readRDS(session_file)$launcher
    if (!is.null(launcher)) {
      current <- resume_oai_process(launcher$pid)
      if (!is.null(current) && identical(current$start, launcher$start))
        ancestors <- c(ancestors, launcher$pid)
    }
  }
  processes <- system2('ps', c('-eo', 'pid=,comm=,args='), stdout = TRUE)
  pattern <- if (any(dir.exists(legacy))) 'oai' else
    paste0('(resume-oai-', group, '[.]r|oai-', group, '[.]R|/oai[.]R)')
  other <- processes[grepl('^ *[0-9]+ +R(script)? +', processes) & grepl(pattern, processes)]
  pids <- as.integer(sub('^ *([0-9]+).*$', '\\1', other))
  # A legacy lock excludes both method groups only within its own output root.
  # A known resume process writing elsewhere is independent; unknown command
  # lines still block recovery, as do all recorded live owners/workers above.
  independent <- vapply(pids, function(pid) {
    root <- resume_oai_process_output(pid)
    !is.null(root) && !identical(root, normalizePath(out))
  }, logical(1))
  other <- other[!pids %in% ancestors & !independent]
  if (length(other)) fail('Other OAI R processes are present:\n', paste(other, collapse = '\n'))
  # Re-read every owner before moving anything, so a newly acquired lock cannot
  # be mistaken for the stale lock inspected above. Preserve the owner record.
  for (lock in locks) {
    owner <- owners[[lock]]
    if (!dir.exists(lock) || !identical(readLines(owner$path, warn = FALSE), owner$snapshot) || live(owner$pid))
      fail('Lock ownership changed during recovery: ', lock)
  }
  archived <- character()
  for (lock in locks) {
    dest <- file.path(out, '_meta/stale-locks', paste0(basename(lock), '-',
                     format(Sys.time(), '%Y%m%d_%H%M%OS6'), '-', Sys.getpid()))
    dir.create(dirname(dest), recursive = TRUE, showWarnings = FALSE)
    if (!file.rename(lock, dest)) fail('Cannot archive stale lock: ', lock)
    cat('Archived stale lock: ', dest, '\n', sep = '')
    archived <- c(archived, dest)
  }
  invisible(archived)
}

resume_oai_available_gib <- function() {
  path <- '/proc/meminfo'
  if (!file.exists(path)) return(NA_real_)
  line <- grep('^MemAvailable:', readLines(path, warn = FALSE), value = TRUE)
  if (length(line) != 1L) return(NA_real_)
  suppressWarnings(as.numeric(strsplit(trimws(line), '\\s+')[[1L]][2L]) / 1024^2)
}

resume_oai_migration_capacity <- function(cores, sizes_gib, available_gib) {
  # Compressed checkpoint size is only a heuristic for import RAM. Keep ample
  # headroom and never import in the parent process.
  per_fit <- max(8, 4 * sizes_gib)
  memory_slots <- if (is.finite(available_gib))
    max(0, floor((available_gib - 16) / per_fit)) else 1L
  list(workers = as.integer(min(cores, length(sizes_gib), memory_slots)),
       per_fit_gib = per_fit, reserve_gib = 16)
}

resume_oai_migrate <- function(paths, cores, chunk_iter, out) {
  paths <- unique(normalizePath(paths, mustWork = TRUE))
  if (!length(paths)) return(invisible(NULL))
  capacity <- resume_oai_migration_capacity(cores, file.info(paths)$size / 1024^3,
                                          resume_oai_available_gib())
  if (capacity$workers < 1L)
    stop('Insufficient available RAM for legacy import with 16 GiB headroom; checkpoints unchanged.', call. = FALSE)
  cat('Converting ', length(paths), ' legacy histories with up to ', capacity$workers,
      ' isolated import workers (within --cores=', cores, '); estimated RAM allowance ',
      round(capacity$per_fit_gib, 1), ' GiB per worker plus 16 GiB headroom.\n', sep = '')
  flush.console()
  logdir <- file.path(out, 'logs', paste0('eivGP_migration_',
    format(Sys.time(), '%Y%m%d_%H%M%OS6'), '_', Sys.getpid()))
  dir.create(logdir, recursive = TRUE, showWarnings = FALSE)
  status_path <- file.path(logdir, 'fit_status.csv')
  status <- data.frame(path = paths, state = 'PENDING', worker_pid = NA_integer_,
                       error = '', stringsAsFactors = FALSE)
  publish <- function() {
    tmp <- tempfile(tmpdir = logdir)
    write.csv(status, tmp, row.names = FALSE)
    if (!file.rename(tmp, status_path)) stop('Cannot publish import worker status.')
  }
  active <- list()
  on.exit({
    if (length(active)) {
      for (job in active) tools::pskill(job$pid)
      parallel::mccollect(active, wait = TRUE)
      status$state[status$state == 'RUNNING'] <- 'STOPPED'
      publish()
    }
  }, add = TRUE)
  publish()
  while (any(status$state == 'PENDING') || length(active)) {
    while (any(status$state == 'PENDING') && length(active) < capacity$workers) {
      available <- resume_oai_available_gib()
      # Account for active imports which may not have allocated their RAM yet.
      if (is.finite(available) && available < capacity$reserve_gib +
          (length(active) + 1L) * capacity$per_fit_gib) {
        if (!length(active)) stop('Available RAM dropped below the legacy-import allowance; retry later.', call. = FALSE)
        break
      }
      i <- which(status$state == 'PENDING')[1L]
      path <- paths[i]
      job <- parallel::mcparallel({
        tryCatch({
          eivGP::continue_eivgp_disk(path, 0L, chunk_iter = chunk_iter,
                                    parallel = FALSE, n_cores = 1L)
          list(ok = TRUE, error = '')
        }, error = function(e) list(ok = FALSE, error = conditionMessage(e)))
      }, mc.set.seed = FALSE, silent = TRUE)
      active[[as.character(job$pid)]] <- job
      status$state[i] <- 'RUNNING'; status$worker_pid[i] <- job$pid
      publish()
      cat('Migrating saved history [PID ', job$pid, ']: ', path, '\n', sep = '')
      flush.console()
    }
    done <- parallel::mccollect(active, wait = FALSE)
    if (length(done)) for (pid in names(done)) {
      i <- match(as.integer(pid), status$worker_pid)
      result <- done[[pid]]
      active[[pid]] <- NULL
      ok <- is.list(result) && isTRUE(result$ok)
      status$state[i] <- if (ok) 'DONE' else 'FAILED'
      status$error[i] <- if (ok) '' else if (is.list(result)) result$error else
        'Import worker exited without a successful result.'
      publish()
      if (ok) {
        cat('Migrated saved history: ', paths[i], '\n', sep = ''); flush.console()
      }
    }
    if (any(status$state == 'FAILED')) {
      i <- which(status$state == 'FAILED')[1L]
      stop('Disk-history migration failed; original fit preserved: ', paths[i],
           '\n', status$error[i], call. = FALSE)
    }
    if (length(active)) Sys.sleep(0.1)
  }
  invisible(status)
}

resume_oai_main <- function(group, script, args = commandArgs(TRUE)) {
  group <- match.arg(group, c("eivGP", "competitors"))
  fail <- function(...) stop(..., call. = FALSE)
  if (any(args %in% c('--help', '-h'))) {
    cat('Usage: Rscript tools/resume-oai-', group, '.r --action=plan|run|report --cores=N\n', sep = '')
    cat('       [--profile=paper] [--data-root=PATH] [--output=PATH]\n',
        'Defaults: action=run, cores=4, frozen settings and inputs from the result root.\n',
        if (group == 'eivGP') 'EIV-GP: cores >= 4; 4/8/12 cores permit 1/2/3 four-chain fits.\n' else
          'Competitors: cores >= 1; each resident optimizer job uses one core.\n',
        'plan verifies code/settings/inputs; run resumes incomplete work;\n',
        'report renders saved results without fitting. --plan=true is also accepted.\n',
        'Run mode records the implementation update, archives confirmed stale\n',
        'workflow locks, reuses saved fits and completes missing work.\n', sep = '')
  } else {
    script <- normalizePath(script)
    repo <- normalizePath(file.path(dirname(script), '..'))
    opts <- list(cores = '4', action = 'run', plan = 'false', profile = NULL, data_root = NULL,
                 output = file.path(repo, 'results', if (group == 'eivGP') 'oai' else 'oai_competitors'))
    seen <- character()
    for (arg in args) {
      if (!grepl('^--(action|cores|plan|profile|data-root|output)=.+$', arg)) fail('Unknown option: ', arg, '; see --help')
      key <- gsub('-', '_', sub('^--([^=]+)=.*$', '\\1', arg))
      if (key %in% seen) fail('Repeated option: ', key)
      seen <- c(seen, key)
      opts[[key]] <- sub('^[^=]+=', '', arg)
    }
    cores <- suppressWarnings(as.numeric(opts$cores))
    minimum <- if (group == 'eivGP') 4L else 1L
    if (length(cores) != 1L || !is.finite(cores) || cores < minimum ||
        cores != floor(cores) || cores > .Machine$integer.max)
      fail('--cores must be an integer >= ', minimum, '.')
    if (!opts$plan %in% c('true', 'false')) fail('--plan must be true or false.')
    if (!opts$action %in% c('plan', 'run', 'report')) fail('--action must be plan, run or report.')
    if (all(c('action', 'plan') %in% seen)) fail('Use --action or --plan, not both.')
    if ('plan' %in% seen) opts$action <- if (opts$plan == 'true') 'plan' else 'run'
    plan <- opts$action == 'plan'
    cores <- as.integer(cores)
    out <- normalizePath(path.expand(opts$output), mustWork = TRUE)
    # Set thread limits BEFORE the worker starts R and loads numerical libraries.
    if (!identical(Sys.getenv('EIVGP_OAI_RESUME_WORKER'), script)) {
      Sys.setenv(OMP_NUM_THREADS = '1', OPENBLAS_NUM_THREADS = '1', MKL_NUM_THREADS = '1',
                 VECLIB_MAXIMUM_THREADS = '1', BLIS_NUM_THREADS = '1', EIVGP_OAI_RESUME_WORKER = script)
      status <- resume_oai_launch_worker(script, args)
      if (status != 0L) fail('OAI resume exited with status ', status, '. Check the error above; checkpoints retained.')
    } else {
      session_file <- Sys.getenv('EIVGP_OAI_SESSION_FILE')
      if (nzchar(session_file)) resume_oai_record_session(session_file)
      version <- read.dcf(file.path(repo, 'eivGP/DESCRIPTION'), fields = 'Version')[1L]
      prepared_library <- file.path(repo, 'results', paste0('eivgp-library-', version))
      if (!dir.exists(file.path(prepared_library, 'eivGP')))
        prepared_library <- file.path(repo, 'results/oai-resume-library')
      if (!nzchar(Sys.getenv('EIVGP_R_LIBRARY')) && dir.exists(file.path(prepared_library, 'eivGP')))
        Sys.setenv(EIVGP_R_LIBRARY = prepared_library)
      .libPaths(c(Sys.getenv('EIVGP_R_LIBRARY'), file.path(repo, '.r-library'), .libPaths()))
      source(file.path(repo, 'tools/oai-tools.R'), local = TRUE)
      api <- oai_load_api(repo)
      cores <- min(cores, getFromNamespace('mixedgp_available_cores', 'eivGP')())
      if (cores < minimum) fail('At least ', minimum, ' available cores are required.')
      manifest <- file.path(out, '_meta/application_manifest.json')
      identity <- file.path(out, '_meta/workflow_identity.rds')
      if (!all(file.exists(c(manifest, identity)))) fail('Missing frozen manifest/identity in ', out)
      bundle <- api$oai_screen_manifest(manifest, api)
      m <- bundle$manifest
      if (!is.null(opts$profile) && !identical(opts$profile, m$profile))
        fail('--profile must match the saved profile: ', m$profile)
      if (!is.null(opts$data_root)) {
        supplied <- api$oai_plan_inputs(normalizePath(path.expand(opts$data_root), mustWork = TRUE))
        if (!identical(unname(tools::md5sum(supplied$files)),
                       unname(tools::md5sum(bundle$inputs[[m$cases[[1]]$id]]$files))))
          fail('--data-root differs from the frozen input snapshot.')
      }
      if (!identical(normalizePath(m$output_root), out) || length(m$cases) != 1L ||
          !identical(m$trial_id, 'oai') || !identical(as.integer(m$chains), 4L))
        fail('This script requires a single-cohort, four-chain OAI manifest.')
      modules <- api$eivgp_module_identity()
      paths <- file.path(repo, 'eivGP', ifelse(grepl('[.]R$', names(modules)), 'R', 'src'), names(modules))
      if (!identical(unname(modules), unname(tools::md5sum(paths))))
        fail('Installed eivGP differs from this checkout. Build/install the current package before resuming.')
      code <- api$oai_source_files(repo)
      current <- list(schema = 'oai-upstream-screening-v2', trial = m$trial_id, modules = modules,
        application = setNames(unname(tools::md5sum(code)), substring(code, nchar(repo) + 2L)),
        burn = m$burn, stages = unlist(m$stages), smoke = isTRUE(m$smoke), ess_min = m$ess_min,
        chains = as.integer(m$chains), samplers = m$samplers, baseline_profile = m$baseline_profile,
        draws_per_chain = m$draws_per_chain, plan = m$cases[[1]]$id,
        inputs = unname(tools::md5sum(bundle$inputs[[m$cases[[1]]$id]]$files)))
      # The user authorized the optimized implementation of the same algorithm.
      # Only code/schema provenance may differ; statistical settings and inputs
      # must still match exactly. Keep the original identity/receipts intact.
      check_identity <- function(saved, now) {
        if (length(saved$schema) != 1L || !saved$schema %in%
            c('adni-upstream-screening-v2', 'oai-upstream-screening-v2')) fail('Unsupported saved workflow identity.')
        omit <- c('schema', 'modules', 'application')
        a <- saved[setdiff(names(saved), omit)]
        b <- now[setdiff(names(now), omit)]
        if (!identical(a, b)) fail('Frozen model settings or inputs differ: ', paste(all.equal(a, b), collapse = '; '))
        invisible(saved)
      }
      saved <- readRDS(identity)
      check_identity(saved, current)
      task_directory <- function(task) {
        if (task$kind == 'cc') api$oai_path(out, 'CC', task$cell, 'baselines',
                                          sprintf('fold%02d', task$fold), task$method)
        else api$oai_path(out, 'baselines', sprintf('fold%02d', task$fold), task$method)
      }
      cached_files <- function(task) file.path(task_directory(task),
        c('fit.rds', 'prediction_draws.rds', 'metrics.csv', 'Y_predictions.csv'))
      if (group == 'competitors') {
        # Older runs had one OI receipt per fold, whereas the current executor
        # schedules each method independently. Reuse its complete saved outputs
        # even when a method-level receipt has not yet been written.
        baseline <- api$baseline_cell
        api$baseline_cell <- function(task, inputs, out, s, cc = FALSE) {
          if (all(file.exists(cached_files(task)))) return(invisible('REUSED'))
          baseline(task, inputs, out, s, cc)
        }
      }
      if (group == 'eivGP' && utils::packageVersion('eivGP') < '0.4.0.9021')
        fail('Install eivGP >= 0.4.0.9021 for the default disk-backed application workflow.')
      changed <- union(names(saved$modules), names(modules))
      changed <- changed[!vapply(changed, function(n)
        identical(unname(saved$modules[n]), unname(modules[n])), logical(1))]
      cat('Verified frozen settings and inputs; current eivGP ', as.character(utils::packageVersion('eivGP')), '\n',
          'Updated modules: ', if (length(changed)) paste(changed, collapse = ', ') else 'none', '\n',
          'Output: ', out, '\n',
          if (group == 'eivGP') paste0('Target: ', paste(unlist(m$stages), collapse = ', '),
            ' iterations; burn: ', m$burn, '; chains: 4; simultaneous fits: ', cores %/% 4L) else
            paste0('Competitor profile: ', m$baseline_profile, '; simultaneous optimizer jobs: ', cores), '\n', sep = '')
      flush.console()
      if (group == 'eivGP') {
        cat('Disk history: at most ', m$checkpoint_step,
            ' new draws per chain in RAM; one fit at a time during full-history analysis.\n', sep = '')
        flush.console()
      }
      launcher <- c(script, file.path(repo, 'tools/resume-oai-tools.R'))
      launcher_sources <- setNames(unname(tools::md5sum(launcher)), substring(launcher, nchar(repo) + 2L))
      api$oai_workflow_identity <- function(path, id) {
        original <- readRDS(path)
        check_identity(original, id)
        record <- list(schema = 'oai-resume-implementation-v1', group = group, original_identity = original,
          current_identity = id, launcher_sources = launcher_sources,
          package_version = as.character(utils::packageVersion('eivGP')),
          statistical_settings_and_inputs_unchanged = TRUE,
          reason = 'User-authorized resumption with the more efficient implementation of the same algorithm.')
        dest <- file.path(dirname(path), 'implementation_updates',
                          paste0(group, '-', digest::digest(list(id, launcher_sources)), '.rds'))
        if (!file.exists(dest)) {
          record$started_at <- format(Sys.time(), tz = 'UTC', usetz = TRUE)
          api$atomic_rds(record, dest)
        }
        cat('Implementation provenance: ', dest, '\n', sep = ''); flush.console()
        original
      }
      if (plan) {
        if (group == 'competitors') {
          tasks <- api$oai_screen_tasks(m, group, 'BASELINE')
          tasks$state <- vapply(seq_len(nrow(tasks)), function(i) {
            task <- tasks[i, , drop = FALSE]
            receipt <- api$oai_screen_receipt(out, task, 'fit')
            if (api$oai_screen_verified(receipt, saved)) 'VERIFIED'
            else if (all(file.exists(cached_files(task)))) 'SAVED_FIT'
            else 'INCOMPLETE'
          }, '')
          print(tasks[, c('fold', 'cell', 'method', 'state')], row.names = FALSE)
          cat('Reusable tasks: ', sum(tasks$state != 'INCOMPLETE'), '/', nrow(tasks),
              '; incomplete: ', sum(tasks$state == 'INCOMPLETE'),
              '; resident optimizer limit: ', min(cores, nrow(tasks)), '\n', sep = '')
        } else api$oai_screen_run(repo, list(manifest = manifest, cores = cores,
                                           plan = TRUE, reuse_only = FALSE), group)
        cat('Plan only: no fitting, provenance writes or lock changes.\n')
        return(invisible(bundle))
      } else resume_oai_recover_locks(out, group)
      # Use the frozen manifest directly, so newer CLI defaults cannot change
      # checkpoint spacing, priors, samplers, masks or the target budget.
      if (opts$action == 'run')
        api$oai_execution_lock(out, group, function() {
          if (group == 'eivGP') {
            # Old do.call() provenance serialized a versioned package namespace.
            # Validate it with the original library; leave completed outputs intact.
            for (fold in unlist(m$folds)) for (cell in unlist(m$cells)) {
              path <- api$oai_plan_fit_path(out, list(fold = fold, cell = cell))
              if (file.exists(file.path(paste0(path, '.draws'), 'index.rds')))
                getFromNamespace('mixedgp_disk_prepare_compatibility', 'eivGP')(
                  path, file.path(repo, 'results/oai-resume-library'))
            }
            # Each isolated import worker exits before sampling begins, so
            # historical arrays never reside in the coordinator or its chains.
            paths <- character()
            for (fold in unlist(m$folds)) for (cell in unlist(m$cells)) {
              path <- api$oai_plan_fit_path(out, list(fold = fold, cell = cell))
              task <- api$oai_screen_tasks(m, 'eivGP', paste0('S', tail(unlist(m$stages), 1L)))
              task <- task[task$fold == fold & task$cell == cell, , drop = FALSE]
              receipt <- api$oai_screen_receipt(out, task, 'evaluate')
              if (api$oai_screen_verified(receipt, saved)) next
              if (file.exists(path) && !file.exists(file.path(paste0(path, '.draws'), 'index.rds'))) {
                paths <- c(paths, path)
              }
            }
            resume_oai_migrate(paths, cores, as.integer(m$checkpoint_step), out)
            cat('All legacy history conversions complete; resuming workflow (completed fits are reused).\n')
            flush.console()
          }
          api$oai_screen_run(repo, list(manifest = manifest, cores = cores,
            plan = FALSE, reuse_only = FALSE, lock_held = TRUE), group)
        })
      settings <- api$oai_settings(profile = m$profile, cores = max(4L, cores),
        stages = as.integer(unlist(m$stages)), burn = as.integer(m$burn),
        folds = as.integer(unlist(m$folds)), cells = unlist(m$cells),
        first_checkpoint = as.integer(m$first_checkpoint), checkpoint_step = as.integer(m$checkpoint_step),
        draws_per_chain = as.integer(m$draws_per_chain), baseline_profile = m$baseline_profile, ess_min = m$ess_min,
        sampler_o0 = m$samplers$O0, sampler_o20 = m$samplers$O20, sampler_o50 = m$samplers$O50)
      settings$cores <- cores
      api$oai_execution_lock(out, group, function() api$oai_render_report(repo, m$cases[[1]]$input,
        out, settings, report_dir = api$oai_path(out, 'reports', group), group = group, filename = 'index.html'))
    }
  }
}

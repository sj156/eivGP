# Author: Xinye Chen. OAI application; all application functions live in this flat tools file.
.oai_repository <- local({
  files<-as.character(unlist(lapply(sys.frames(),function(frame) {
    c(if(is.character(frame$ofile))frame$ofile,if(is.character(frame$file))frame$file)
  }),use.names=FALSE))
  files<-files[basename(files)=='oai-tools.R']
  path<-if(length(files))tail(files,1L)else sub('^--file=','',grep('^--file=',commandArgs(FALSE),value=TRUE)[1])
  normalizePath(file.path(dirname(normalizePath(path)),'..'))
})
oai_load_api <- function(repo=.oai_repository) {
  repo<-normalizePath(repo)
  lib<-Sys.getenv('EIVGP_R_LIBRARY',file.path(repo,'.r-library'))
  if(nzchar(lib)&&dir.exists(lib)).libPaths(c(normalizePath(lib),.libPaths()))
  if(!requireNamespace('eivGP',quietly=TRUE)||utils::packageVersion('eivGP')<'0.4.0.9014')
    stop('Install eivGP >= 0.4.0.9014 before running this application.')
  ns<-asNamespace('eivGP');e<-new.env(parent=ns)
  evalq({
oai_labels <- function(fields,labels=NULL) {
  values<-unlist(labels,use.names=TRUE)
  vapply(seq_along(fields),function(i) {
    if(fields[i]%in%names(values))as.character(values[[fields[i]]])
    else if(length(values)>=i)as.character(values[i])else fields[i]
  },character(1),USE.NAMES=FALSE)
}

# Task folders and shared directories mirror the Ocean flat result layout.
oai_storage_flat <- function(root) {
  marker<-base::file.path(root,'_meta','storage-layout.rds')
  if(!file.exists(marker))return(FALSE)
  if(!readRDS(marker)$schema%in%c('oai-flat-v1','adni-flat-v1'))stop('Unknown OAI storage schema: ',marker)
  TRUE
}
oai_storage_init <- function(root,stages) {
  if(oai_storage_flat(root))return(invisible(root))
  existing<-if(dir.exists(root))list.files(root,all.files=TRUE,no..=TRUE)else character()
  # Existing unmarked results retain their original physical layout.
  if(length(existing))return(invisible(root))
  dir.create(base::file.path(root,'_meta'),recursive=TRUE,showWarnings=FALSE)
  mixedgp_atomic_rds(list(schema='oai-flat-v1',stages=as.integer(stages)),
    base::file.path(root,'_meta','storage-layout.rds'))
  invisible(root)
}
oai_flat_relative <- function(path,stages) {
  path<-gsub('\\','/',path,fixed=TRUE)
  stage<-regexec('^(S[0-9]+K?)/analysis/(fold[0-9]+)/(O[0-9]+)(/.*)?$',path)
  hit<-regmatches(path,stage)[[1]]
  if(length(hit)) {
    task<-paste0('eiv_',hit[3],'_',hit[4]);suffix<-hit[5]
    if(hit[2]!=paste0('S',tail(stages,1L)))task<-paste0(task,'/stages/',hit[2])
    return(paste0(task,suffix))
  }
  path<-sub('^fits/(fold[0-9]+)/(O[0-9]+)(?=/|$)','eiv_\\1_\\2',path,perl=TRUE)
  cc<-regexec('^CC/(O[0-9]+)/baselines/(fold[0-9]+)/([^/]+)(/.*)?$',path)
  hit<-regmatches(path,cc)[[1]]
  if(length(hit))return(paste0('cc_',hit[3],'_',hit[2],'_',gsub('-','',hit[4]),hit[5]))
  path<-sub('^baselines/(fold[0-9]+)(?=/|$)','baseline_\\1',path,perl=TRUE)
  path<-sub('^receipts(?=/|$)','_meta/tasks',path,perl=TRUE)
  path<-sub('^workflow/logs(?=/|$)','logs',path,perl=TRUE)
  path<-sub('^workflow/reports(?=/|$)','reports/prediction_report',path,perl=TRUE)
  path<-sub('^report-assets/paper(?=/|$)','reports/paper',path,perl=TRUE)
  path<-sub('^report-assets/functionals(?=/|$)','reports/diagnostics/functionals',path,perl=TRUE)
  path<-sub('^report-assets(?=/|$)','reports/diagnostics',path,perl=TRUE)
  path<-sub('^(reports/(eivGP|competitors))/report-assets/paper(?=/|$)','\\1/paper',path,perl=TRUE)
  path<-sub('^(reports/(eivGP|competitors))/report-assets(?=/|$)','\\1/diagnostics',path,perl=TRUE)
  if(grepl('^(application_manifest[.]json|workflow_identity[.]rds|[^/]+_budget[.]txt)$',path))path<-paste0('_meta/',path)
  path
}
oai_path <- function(...,fsep=.Platform$file.sep) {
  paths<-base::file.path(...,fsep=fsep)
  vapply(paths,function(path) {
    candidate<-path
    repeat {
      marker<-base::file.path(candidate,'_meta','storage-layout.rds')
      if(dir.exists(candidate)&&file.exists(marker)) {
        oai_storage_flat(candidate)
        if(identical(path,candidate))return(path)
        rel<-substring(path,nchar(candidate)+2L)
        return(base::file.path(candidate,oai_flat_relative(rel,readRDS(marker)$stages)))
      }
      parent<-dirname(candidate)
      if(identical(parent,candidate)||candidate=='.')return(path)
      candidate<-parent
    }
  },character(1),USE.NAMES=FALSE)
}
oai_workflow_root <- function(out)if(oai_storage_flat(out))out else base::file.path(out,'workflow')
oai_default_output <- function(dataset,profile,repo) {
  base::file.path(repo,'results',if(profile=='paper')dataset else paste(dataset,profile,sep='_'))
}
oai_snapshot_inputs <- function(input,out,files) {
  dest<-oai_path(out,'data');meta<-base::file.path(input,'DATASET.json')
  if(file.exists(meta))files<-c(files,meta)
  for(source in files) {
    target<-base::file.path(dest,oai_input_basename(source))
    if(file.exists(target)) {
      if(!identical(unname(tools::md5sum(source)),unname(tools::md5sum(target))))
        stop('Frozen input snapshot differs: ',target)
    }else {
      dir.create(dirname(target),recursive=TRUE,showWarnings=FALSE)
      if(!file.copy(source,target)||!identical(unname(tools::md5sum(source)),unname(tools::md5sum(target))))
        stop('Could not verify frozen input snapshot: ',target)
    }
  }
  dest
}

oai_source_files <- function(repo) {
  oai_path(repo,'tools',c('oai-tools.R','oai.R','oai-eivGP.R',
    'oai-competitors.R','oai-report.R','oai-prepare-folds.R'))
}
# OAI application: data. Reuses package model and sampler functions.

oai_ezgp_mixed_matrix <- function(X, C, m_vec) {
    X <- as.matrix(X)
    C <- as.matrix(C)
    if (!is.numeric(X) || any(!is.finite(X)))
        stop("X must be a finite numeric matrix.")
    if (nrow(X) != nrow(C) || ncol(C) != length(m_vec))
        stop("X/C dimensions do not match.")
    if (any(!is.finite(C)) || any(C != round(C)))
        stop("C must contain finite integer codes.")
    for (j in seq_along(m_vec)) {
        if (any(C[, j] < 1L | C[, j] > m_vec[j]))
            stop("C contains an out-of-range category.")
    }
    out <- cbind(X, C)
    storage.mode(out) <- "double"
    colnames(out) <- c(paste0("x", seq_len(ncol(X))), paste0("c", seq_len(ncol(C))))
    out
}

oai_ezgp_cross_collisions <- function(train_x, test_x, digits = 5L) {
    train_key <- apply(round(as.matrix(train_x), digits), 1L, paste, collapse = "\r")
    test_key <- apply(round(as.matrix(test_x), digits), 1L, paste, collapse = "\r")
    sum(test_key %in% train_key)
}

oai_ezgp_predict_covariance_consistent <- function(X_new, model) {
    mixedgp_ezgp_predict(X_new, model)
}

oai_ezgp_select_tau_cv <- function(train_x, y_train, p, q, m_vec, seed, tau_fractions, cv_folds = 3L, maxeval = 100L,
    progress = NULL, initial_audit = NULL) {
    n <- nrow(train_x)
    cv_folds <- max(2L, min(as.integer(cv_folds), n))
    set.seed(as.integer(seed))
    fold_id <- sample(rep(seq_len(cv_folds), length.out = n))
    response_var <- max(stats::var(y_train), .Machine$double.eps)
    tau_fractions <- unique(as.numeric(tau_fractions))
    tau_grid <- tau_fractions * response_var
    audit <- list()
    index <- 0L
    for (tt in seq_along(tau_grid)) {
        for (fold in seq_len(cv_folds)) {
            index <- index + 1L
            train_id <- fold_id != fold
            validation_id <- !train_id
            prior <- if (is.data.frame(initial_audit)) {
                initial_audit[initial_audit$tau_index == tt & initial_audit$inner_fold == fold, , drop = FALSE]
            }
            else data.frame()
            if (nrow(prior) == 1L) {
                audit[[index]] <- prior
                if (is.function(progress))
                  progress(do.call(rbind, audit), index, length(tau_grid) * cv_folds)
                next
            }
            started <- proc.time()[["elapsed"]]
            warnings <- character()
            ans <- tryCatch(withCallingHandlers({
                fit <- EzGP::EzGP_fit(X = train_x[train_id, , drop = FALSE], Y = as.numeric(y_train[train_id]),
                  p = p, q = q, m = m_vec, tau = tau_grid[tt], maxeval = as.integer(maxeval))
                if (length(fit$param) == 0L || any(!is.finite(fit$param))) {
                  stop("EzGP returned non-finite fitted parameters.")
                }
                pred <- oai_ezgp_predict_covariance_consistent(train_x[validation_id, , drop = FALSE], fit)
                error <- pred$Y_hat - y_train[validation_id]
                if (any(pred$predictive_MSE <= 0))
                  stop("EzGP CV returned zero predictive variance.")
                list(nlpd = mean(0.5 * log(2 * pi * pred$predictive_MSE) + 0.5 * error^2/pred$predictive_MSE),
                  mse = mean(error^2), message = "")
            }, warning = function(w) {
                warnings <<- c(warnings, conditionMessage(w))
                invokeRestart("muffleWarning")
            }), error = function(e) list(nlpd = NA_real_, mse = NA_real_, message = conditionMessage(e)))
            audit[[index]] <- data.frame(tau_index = tt, tau_fraction = tau_fractions[tt], tau = tau_grid[tt],
                inner_fold = fold, n_train = sum(train_id), n_validation = sum(validation_id), rounded_cross_collisions = oai_ezgp_cross_collisions(train_x[train_id,
                  , drop = FALSE], train_x[validation_id, , drop = FALSE]), nlpd = ans$nlpd, mse = ans$mse, elapsed_seconds = unname(proc.time()[["elapsed"]] -
                  started), message = ans$message, warnings = paste(unique(warnings), collapse = " | "), stringsAsFactors = FALSE)
            if (is.function(progress))
                progress(do.call(rbind, audit), index, length(tau_grid) * cv_folds)
        }
    }
    audit <- do.call(rbind, audit)
    candidate <- do.call(rbind, lapply(seq_along(tau_grid), function(tt) {
        z <- audit[audit$tau_index == tt, , drop = FALSE]
        data.frame(tau_index = tt, tau_fraction = tau_fractions[tt], tau = tau_grid[tt], cv_nlpd = if (all(is.finite(z$nlpd)))
            mean(z$nlpd)
        else Inf, cv_mse = if (all(is.finite(z$mse)))
            mean(z$mse)
        else Inf, successful_folds = sum(is.finite(z$nlpd)), stringsAsFactors = FALSE)
    }))
    if (!any(is.finite(candidate$cv_nlpd))) {
        stop("Covariance-consistent EzGP training-only CV failed for every nugget value.")
    }
    best <- which.min(candidate$cv_nlpd)
    list(tau = candidate$tau[best], selected_index = candidate$tau_index[best], cv_score = "nlpd", fold_id = fold_id,
        audit = audit, candidates = candidate)
}

oai_ezgp_draws <- function(mean, variance, n_draw, seed) {
    if (any(!is.finite(mean)) || any(!is.finite(variance)) || any(variance < 0)) {
        stop("Invalid predictive mean or variance.")
    }
    set.seed(as.integer(seed))
    matrix(stats::rnorm(n_draw * length(mean), rep(mean, each = n_draw), rep(sqrt(variance), each = n_draw)), nrow = n_draw,
        ncol = length(mean))
}

oai_plan_require <- function() {
    for (p in c("jsonlite", "digest", "kergp", "LVGP", "EzGP")) if (!requireNamespace(p, quietly = TRUE))
        stop("Required OAI dependency: ", p)
}

# Flat filenames encode fold and calibration; older prepared trees remain readable.
oai_input_basename <- function(path) {
    name <- basename(path)
    if (name %in% c('train.csv','test.csv','PROTOCOL.json')) {
        normalized <- gsub('\\', '/', path, fixed=TRUE)
        hit <- regmatches(normalized, regexec('(fold[0-9]+)/(O[0-9]+)/(private/)?(train[.]csv|test[.]csv|PROTOCOL[.]json)$', normalized))[[1]]
        if (length(hit)) name <- paste(hit[2], hit[3], hit[5], sep='_')
    }
    name
}
oai_input_files <- function(root,fold,cell) {
    prefix <- sprintf('fold%02d_%s_',fold,cell)
    flat <- base::file.path(root,paste0(prefix,c('PROTOCOL.json','train.csv','test.csv')))
    if (any(file.exists(flat))) return(flat)
    base::file.path(root,sprintf('fold%02d',fold),cell,c('PROTOCOL.json','private/train.csv','private/test.csv'))
}
oai_plan_input_root <- function(path) {
    path <- path.expand(path)
    if (!file.exists(base::file.path(path,'fold01_O0_PROTOCOL.json'))) {
        if (dir.exists(base::file.path(path,'MC'))) path <- base::file.path(path,'MC')
        else if (dir.exists(base::file.path(path,'inputs','MC'))) path <- base::file.path(path,'inputs','MC')
    }
    normalizePath(path, mustWork = TRUE)
}

oai_plan_inputs <- function(data_dir) {
    if (!requireNamespace("jsonlite", quietly = TRUE))
        stop("Install jsonlite to read OAI protocols.")
    root <- oai_plan_input_root(data_dir)
    all <- list()
    audit <- list()
    files <- character()
    held <- list()
    reference <- NULL
    for (f in 1:5) {
        by_cell <- list()
        for (cell in c("O0", "O20", "O50")) {
            dir <- oai_path(root, sprintf("fold%02d", f), cell)
            paths <- oai_input_files(root,f,cell)
            if (!all(file.exists(paths)))
                stop("Missing private fold inputs: ", dir)
            p <- jsonlite::fromJSON(paths[1], simplifyVector = TRUE)
            tr <- read.csv(paths[2], stringsAsFactors = FALSE)
            te <- read.csv(paths[3], stringsAsFactors = FALSE)
            q <- as.integer(p$q)
            d <- as.integer(p$d)
            xc <- grep("^X[0-9]+$", names(tr), value = TRUE)
            xc <- xc[order(as.integer(sub("X", "", xc)))]
            cc <- paste0("C", seq_len(q))
            uc <- paste0("U", seq_len(d))
            needed <- c("RID", "Y", xc, cc, uc, "calibration")
            if (!length(xc) || !all(needed %in% names(tr)) || !all(needed %in% names(te)))
                stop("Invalid fold columns: ", dir)
            if(length(p$X_order)!=length(xc)||p$p!=length(xc))stop('Continuous-X labels/dimensions differ: ',dir)
            if(any(vapply(tr[,xc,drop=FALSE],function(v)length(unique(v))<=2L,TRUE)))
                stop('X must be continuous; fix binary/categorical population attributes before fitting: ',dir)
            if(any(te$calibration!=0L))stop('Test calibration must be zero: ',dir)
            if (anyNA(tr[, needed]) || anyNA(te[, needed]) || any(!is.finite(as.matrix(tr[, setdiff(needed, "RID")]))) ||
                any(!is.finite(as.matrix(te[, setdiff(needed, "RID")]))))
                stop("Non-finite private inputs: ", dir)
            if (anyDuplicated(tr$RID) || anyDuplicated(te$RID) || length(intersect(tr$RID, te$RID)))
                stop("Training/test ID overlap: ", dir)
            if (nrow(tr) != p$training_N || nrow(te) != p$test_N || p$fold != f || p$O != cell || nrow(tr) + nrow(te) !=
                p$N)
                stop("Protocol row counts do not match: ", dir)
            if (!all(tr$calibration %in% 0:1) || sum(tr$calibration) != p$observed_training_U)
                stop("Calibration mask differs: ", dir)
            if (!is.null(p$train_md5) && !identical(tolower(p$train_md5), unname(tools::md5sum(paths[2]))))
                stop("Training checksum differs: ", dir)
            if (!is.null(p$test_md5) && !identical(tolower(p$test_md5), unname(tools::md5sum(paths[3]))))
                stop("Test checksum differs: ", dir)
            if (length(p$U_order) != d || length(p$C_order) != q || length(p$m_vec) != q)
                stop("Protocol coordinate labels differ: ", dir)
            for (j in seq_len(q)) if (!all(tr[[cc[j]]] %in% seq_len(p$m_vec[j])) || !all(te[[cc[j]]] %in% seq_len(p$m_vec[j])) ||
                !all(seq_len(p$m_vec[j]) %in% tr[[cc[j]]]))
                stop("Invalid native ordinal levels: ", dir)
            observed <- tr$calibration == 1L
            available <- as.character(p$test_U_available_RID)
            if (cell == "O0") {
                if (any(observed) || length(available) || p$ident != "lower_triangular" || isTRUE(p$standardize_U))
                  stop("O0 must use unanchored model coordinates.")
            }
            else {
                if (!any(observed) || !isTRUE(p$standardize_U) || p$ident != "none")
                  stop("Calibrated U must use the stated raw-coordinate anchoring.")
                if (!all(available %in% as.character(te$RID)) || anyDuplicated(available))
                  stop("Invalid fixed test-U mask: ", dir)
                if (qr(cbind(1, as.matrix(tr[observed, uc, drop = FALSE])))$rank < d + 1L)
                  stop("Calibration vectors have insufficient affine rank.")
            }
            seed_names <- c("production", "posterior", "prediction", "mixed_prediction", "UC_GP", "LVGP", "EzGP",
                "EzGP_draws")
            if (!all(seed_names %in% names(p$seeds)))
                stop("Missing fixed protocol seeds: ", dir)
            for (n in seed_names) mixedgp_as_integer_strict(p$seeds[[n]], paste("seed", n), 0L, 1L)
            by_cell[[cell]] <- list(train = tr, test = te, protocol = p, x = xc, c = cc, u = uc, available = te$RID %in%
                available)
            files <- c(files, paths)
            audit[[length(audit) + 1L]] <- data.frame(fold = f, cell = cell, n_train = nrow(tr), n_test = nrow(te),
                n_calibration = sum(observed), n_test_U = sum(te$RID %in% available), q = q, d = d)
        }
        a <- by_cell$O0
        for (cell in c("O20", "O50")) {
            b <- by_cell[[cell]]
            if (!identical(a$train[, setdiff(names(a$train), "calibration")], b$train[, setdiff(names(b$train),
                "calibration")]) || !identical(a$test, b$test))
                stop("Cells must share the same frozen fold records.")
        }
        if (any(by_cell$O20$train$calibration > by_cell$O50$train$calibration) || any(by_cell$O20$available & !by_cell$O50$available))
            stop("Calibration/test-U masks must be nested.")
        ids <- sort(as.character(c(a$train$RID, a$test$RID)))
        if (is.null(reference))
            reference <- ids
        else if (!identical(reference, ids))
            stop("Fold populations differ.")
        held[[f]] <- as.character(a$test$RID)
        for (cell in names(by_cell)) all[[paste(f, cell, sep = "/")]] <- by_cell[[cell]]
    }
    if (anyDuplicated(unlist(held)) || !identical(sort(unlist(held)), reference))
        stop("Each person must be held out exactly once.")
    list(root = root, cells = all, audit = do.call(rbind, audit), files = files, N = length(reference))
}

oai_plan_fit_args <- function(z, s) {
    p <- z$protocol
    u <- as.matrix(z$train[, z$u, drop = FALSE])
    colnames(u) <- p$U_order
    u[z$train$calibration == 0L, ] <- NA_real_
    cmat <- as.matrix(z$train[, z$c, drop = FALSE])
    colnames(cmat) <- p$C_order
    # The input contract encodes each native item as1:m; JSON may simplify nested arrays.
    levels <- lapply(as.integer(p$m_vec), function(m) as.numeric(seq_len(m)))
    list(X = as.matrix(z$train[, z$x, drop = FALSE]), y = z$train$Y, C = cmat, U_obs = u, latent_dim = p$d, engine = "multivariate",
        m_vec = as.integer(p$m_vec), ordinal_levels = levels, ident = p$ident, standardize_U = isTRUE(p$standardize_U),
        kernel = "se", store_scores = TRUE, thin = 1L,
        sampler = if (!is.null(s$samplers)) s$samplers[[p$O]] else if (p$O == "O50") "ESS" else "HMC",
        parallel = s$cores > 1L, n_cores = s$cores,
        n_chains = s$chains, seed = p$seeds$production, verbose = TRUE)
}

oai_plan_fit_path <- function(out, task) oai_path(out, "fits", sprintf("fold%02d", task$fold), task$cell, "fit.rds")

oai_plan_read_fit <- function(path) mixedgp_read_artifact(path)

oai_plan_save_fit <- function(fit, path) {
    dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
    tmp <- tempfile(".oai-fit-", tmpdir = dirname(path))
    on.exit(unlink(tmp), add = TRUE)
    packed <- mixedgp_pack_artifact(fit)
    saveRDS(packed, tmp)
    restored <- readRDS(tmp)
    if (!identical(restored, packed))
        stop("OAI checkpoint read-back differs; previous checkpoint retained.")
    rm(restored)
    gc(FALSE)
    if (!file.rename(tmp, path))
        stop("Could not publish OAI checkpoint; previous retained: ", path)
    mixedgp_retire_saved_checkpoints(fit, path)
    invisible(path)
}

oai_plan_progress <- function(fit,path,task,s,event,elapsed) {
    row<-data.frame(updated_at=format(Sys.time(),tz='UTC',usetz=TRUE),fold=task$fold,cell=task$cell,
      event=event,iteration=fit$checkpoint$iteration,target_iteration=s$iter,burn=fit$checkpoint$burn,
      retained_per_chain=fit$checkpoint$iteration-fit$checkpoint$burn,chains=length(fit$checkpoint$states),
      elapsed_seconds=as.numeric(elapsed),fit_bytes=file.info(path)$size,
      checkpoint_file=basename(path),stringsAsFactors=FALSE)
    history<-oai_path(dirname(path),'mcmc_history.csv')
    old<-if(file.exists(history))read.csv(history,stringsAsFactors=FALSE)else NULL
    mixedgp_atomic_csv(rbind(old,row),history)
    mixedgp_atomic_csv(row,oai_path(dirname(path),'mcmc_progress.csv'))
    invisible(row)
}

oai_plan_fit_task <- function(task, inputs, out, s, extend = FALSE) {
    z <- inputs$cells[[paste(task$fold, task$cell, sep = "/")]]
    path <- oai_plan_fit_path(out, task)
    old <- options(eivGP.checkpoint_root = oai_path(dirname(path),'checkpoints'),
      eivGP.checkpoint_iterations=s$checkpoint_step,eivGP.checkpoint_seconds=Inf,
      eivGP.checkpoint_after_seconds=0)
    on.exit(options(old), add = TRUE)
    started<-proc.time()[['elapsed']]
    fit <- if (file.exists(path))
        oai_plan_read_fit(path)
    else NULL
    if(!is.null(fit)) {
      if(!extend && fit$checkpoint$iteration>s$iter && is.null(s$ess_min))
        stop('Cannot refit an earlier stage after continuation; restore its completed evaluation receipt')
      if(fit$checkpoint$burn!=s$burn||fit$control$sampler!=s$samplers[[task$cell]])
        stop('Checkpoint sampler/burn differs')
    }
    target <- s$iter
    request <- oai_path(dirname(path), "continuation_request.rds")
    if (extend) {
        if (is.null(fit))
            stop("Missing initial fit: ", path)
        tab <- oai_plan_diagnostics(fit)
        if (!file.exists(request)) {
            if (isTRUE(tab$summary$core_pass))
                return(invisible("SKIPPED_CONVERGED"))
            mixedgp_atomic_rds(list(from = fit$checkpoint$iteration, target = fit$checkpoint$iteration + s$continuation_iter,
                added = s$continuation_iter), request)
        }
        req <- readRDS(request)
        if (req$added != s$continuation_iter)
            stop("Saved continuation length differs; use its original --continuation-iter.")
        target <- req$target
    }
    else if (file.exists(request))
        target <- readRDS(request)$target
    if (is.null(fit)) {
        first <- min(target, s$first_checkpoint)
        fit <- do.call(fit_eivgp, c(oai_plan_fit_args(z, s), list(n_iter = first, burn = s$burn)))
        oai_plan_save_fit(fit, path)
        oai_plan_progress(fit,path,task,s,'SAVED',proc.time()[['elapsed']]-started)
    }
    else oai_plan_progress(fit,path,task,s,'RESUMED',proc.time()[['elapsed']]-started)
    # Later segments recover from the published full fit. Avoid duplicating its
    # entire retained history in the package's per-call restart identity.
    options(eivGP.checkpoint_root=NULL)
    while (fit$checkpoint$iteration < target) {
        current <- fit$checkpoint$iteration
        next_target <- min(target, (floor(current/s$checkpoint_step) + 1L) * s$checkpoint_step)
        fit <- continue_eivgp(fit, n_iter = next_target - current, parallel = s$cores > 1L, n_cores = s$cores,
            verbose = TRUE)
        oai_plan_save_fit(fit, path)
        oai_plan_progress(fit,path,task,s,'SAVED',proc.time()[['elapsed']]-started)
    }
    if (fit$checkpoint$burn != s$burn || length(fit$checkpoint$states) != s$chains)
        stop("Saved fit settings differ.")
    if(isTRUE(s$final_stage) && length(fit$checkpoint$states)>=4L) {
      X<-as.matrix(z$test[,z$x,drop=FALSE]);C<-as.matrix(z$test[,z$c,drop=FALSE])
      rows<-mixedgp_study2_diagnostic_rows(X,C,5L)
      panel<-list(X=X[rows,,drop=FALSE],C=C[rows,,drop=FALSE],
        U=if(task$cell!="O0")as.matrix(z$test[rows,z$u,drop=FALSE])else NULL,
        n_latent=64L,seed=z$protocol$seeds$prediction)
      fit<-mixedgp_application_ess(fit,s$ess_min,panel,oai_path(dirname(path),"scientific-ess-state.rds"),s$cores)
      oai_plan_save_fit(fit,path)
      mixedgp_save_scientific_trace(fit,oai_path(dirname(path),"scientific-trace-targets"),
        sprintf("OAI | fold %d | %s",task$fold,task$cell))
      unlink(oai_path(dirname(path),"scientific-ess-state.rds"))
    }
    dg <- oai_plan_diagnostics(fit)
    mixedgp_atomic_csv(dg$parameters, oai_path(dirname(path), "diagnostic_parameters.csv"))
    mixedgp_atomic_csv(dg$summary, oai_path(dirname(path), "diagnostic_summary.csv"))
    oai_plan_progress(fit,path,task,s,'FIT_COMPLETE',proc.time()[['elapsed']]-started)
    invisible("DONE")
}

oai_plan_diagnostics <- function(fit) {
    parts <- fit$diagnostics[grep("^rhat_", names(fit$diagnostics))]
    parts <- Filter(function(z) is.data.frame(z) && nrow(z)>0L && all(c("parameter", "rhat", "ess_bulk", "ess_tail") %in% names(z)),
        parts)
    d <- do.call(rbind, lapply(names(parts), function(n) {
        z <- parts[[n]][, c("parameter", "rhat", "ess_bulk", "ess_tail")]
        z$group <- n
        z
    }))
    core <- d[!grepl("^dictionary_occupancy", d$parameter), , drop = FALSE]
    pass <- function(z) nrow(z) > 0L && all(is.finite(z$rhat) & is.finite(z$ess_bulk) & is.finite(z$ess_tail) &
        z$rhat <= 1.01 & z$ess_bulk >= 400 & z$ess_tail >= 400)
    mx <- function(x) if (any(is.finite(x)))
        max(x[is.finite(x)])
    else NA_real_
    mn <- function(x) if (any(is.finite(x)))
        min(x[is.finite(x)])
    else NA_real_
    list(parameters = d, summary = data.frame(iterations = fit$checkpoint$iteration, burn = fit$checkpoint$burn,
        chains = length(fit$checkpoint$states), max_rhat = mx(core$rhat), min_bulk_ESS = mn(core$ess_bulk), min_tail_ESS = mn(core$ess_tail),
        core_nonfinite = sum(!is.finite(core$rhat) | !is.finite(core$ess_bulk) | !is.finite(core$ess_tail)), core_pass = pass(core),
        strict_all_pass = pass(d)))
}

oai_plan_score <- function(draws, point, truth, training_sd) {
    draws <- as.matrix(draws)
    point <- as.numeric(point)
    if (ncol(draws) != length(truth) || any(!is.finite(draws)) || any(!is.finite(point)) || !is.finite(training_sd) ||
        training_sd <= 0)
        stop("Invalid held-out predictive scores.")
    qq <- apply(draws, 2, quantile, c(0.025, 0.975))
    r <- sqrt(mean((point - truth)^2))
    a <- mean(abs(point - truth))
    crps <- mean(vapply(seq_along(truth), function(j) {
        v <- sort(draws[, j])
        n <- length(v)
        mean(abs(v - truth[j])) - sum((2 * seq_len(n) - n - 1) * v)/n^2
    }, 0))
    width <- mean(qq[2, ] - qq[1, ])
    interval <- mean(qq[2, ] - qq[1, ] + 40 * pmax(qq[1, ] - truth, 0) + 40 * pmax(truth - qq[2, ], 0))
    data.frame(n = length(truth), RMSE = r, RMSE_SD = r/training_sd, MAE = a, MAE_SD = a/training_sd, CRPS = crps,
        CRPS_SD = crps/training_sd, Coverage95 = mean(truth >= qq[1, ] & truth <= qq[2, ]), Width95 = width, Width95_SD = width/training_sd,
        IntervalScore95 = interval, IntervalScore95_SD = interval/training_sd)
}

oai_plan_baseline <- function(task, inputs, out, s) {
    dir <- oai_path(out, "baselines", sprintf("fold%02d", task$fold), task$method)
    path <- oai_path(dir, "fit.rds")
    if (all(file.exists(oai_path(dir, c("fit.rds", "prediction_draws.rds")))))
        return(invisible("DONE"))
    z <- inputs$cells[[paste(task$fold, "O50", sep = "/")]]
    tr <- z$train
    te <- z$test
    p <- z$protocol
    mu <- colMeans(tr[, z$x, drop = FALSE])
    ss <- vapply(tr[, z$x, drop = FALSE], sd, 0)
    if (any(ss <= 0))
        stop("Constant OAI training X.")
    sx <- function(d) sweep(sweep(as.matrix(d[, z$x, drop = FALSE]), 2, mu, "-"), 2, ss, "/")
    xc <- sx(tr)
    xs <- sx(te)
    yc <- mean(tr$Y)
    ys <- sd(tr$Y)
    yt <- (tr$Y - yc)/ys
    ct <- as.matrix(tr[, z$c, drop = FALSE])
    cs <- as.matrix(te[, z$c, drop = FALSE])
    mvec <- as.integer(p$m_vec)
    if (task$method != "EzGP") {
        controls <- if (s$profile == "paper") mixedgp_competitor_controls()[[task$method]] else
            if (task$method == "LVGP") list(n_starts=1L,max_retries=1L,max_iter_ini=10L,max_iter_lat=3L,
                rescue_iter_ini=10L,rescue_iter_lat=3L,max_elapsed_seconds=120,parallel=TRUE) else
                list(n_starts=1L,rescue_starts=integer())
        if(task$method == "LVGP") controls$parallel <- FALSE
        ans <- do.call(fit_mixedgp_competitor, c(list(method = task$method, X = xc, y = yt, C = ct, new_X = xs,
            new_C = cs, m_vec = mvec, n_draw = 4L * s$draws_per_chain, seed = p$seeds[[if (task$method == "UC-GP") "UC_GP" else "LVGP"]]),
            controls))
        point <- yc + ys * as.numeric(ans$predictive_mean)
        draws <- yc + ys * as.matrix(ans$draws)
        model <- ans
    }
    else {
        tx <- oai_ezgp_mixed_matrix(xc, ct, mvec)
        vx <- oai_ezgp_mixed_matrix(xs, cs, mvec)
        sel <- oai_ezgp_select_tau_cv(tx, yt, ncol(xc), length(z$c), mvec, p$seeds$EzGP, if (s$profile == "paper")
            c(1e-06, 0.0025, 0.01, 0.04, 0.16)
        else c(0.0025, 0.04), cv_folds = if (s$profile == "paper")
            3L
        else 2L, maxeval = if (s$profile == "paper")
            100L
        else 8L)
        model <- EzGP::EzGP_fit(tx, yt, p = ncol(xc), q = length(z$c), m = mvec, tau = sel$tau, maxeval = if (s$profile ==
            "paper")
            100L
        else 8L)
        pp <- oai_ezgp_predict_covariance_consistent(vx, model)
        point <- yc + ys * pp$Y_hat
        draws <- yc + ys * oai_ezgp_draws(pp$Y_hat, pp$predictive_MSE, 4L * s$draws_per_chain, p$seeds$EzGP_draws)
        mixedgp_atomic_csv(sel$audit, oai_path(dir, "tau_cv_audit.csv"))
        mixedgp_atomic_csv(sel$candidates, oai_path(dir, "tau_cv_candidates.csv"))
    }
    mixedgp_atomic_rds(list(model = model, center = mu, scale = ss, y_center = yc, y_scale = ys), path)
    mm <- cbind(data.frame(fold = task$fold, cell = "OI", method = task$method, protocol = "OI"), oai_plan_score(draws,
        point, te$Y, ys))
    mixedgp_atomic_csv(mm, oai_path(dir, "metrics.csv"))
    mixedgp_atomic_rds(list(draws = draws, point = point), oai_path(dir, "prediction_draws.rds"))
    invisible("DONE")
}

oai_plan_analyze_cell <- function(task, inputs, out, s) {
    path <- oai_plan_fit_path(out, task)
    if (!file.exists(path))
        stop("Missing saved OAI fit: ", path)
    fit <- oai_plan_read_fit(path)
    z <- inputs$cells[[paste(task$fold, task$cell, sep = "/")]]
    p <- z$protocol
    tr <- z$train
    te <- z$test
    dir <- oai_path(out, "analysis", sprintf("fold%02d", task$fold), task$cell)
    dir.create(dir, recursive = TRUE, showWarnings = FALSE)
    n <- nrow(fit$mcmc$samples_tau)
    cn <- n%/%s$chains
    set.seed(p$seeds$posterior)
    ids <- sort(unlist(lapply(0:(s$chains - 1L), function(i) i * cn + sample.int(cn, min(s$draws_per_chain, cn)))))
    ct <- as.matrix(te[, z$c, drop = FALSE])
    colnames(ct) <- p$C_order
    yd <- predict_eivgp(fit, new_X = as.matrix(te[, z$x, drop = FALSE]), new_C = ct, target = "response", draw_ids = ids,
        n_per_draw = 1L, joint = FALSE, seed = p$seeds$prediction)
    conditional <- attr(yd, "conditional_means", exact = TRUE)
    point <- if (is.null(conditional))
        colMeans(as.matrix(yd))
    else colMeans(conditional)
    score <- function(draws, point, protocol) cbind(data.frame(fold = task$fold, cell = task$cell, method = "EIV-GP",
        protocol = protocol), oai_plan_score(draws, point, te$Y, sd(tr$Y)))
    metrics <- score(yd, point, "OI")
    predictions <- data.frame(row = seq_len(nrow(te)), observed_Y = te$Y, predicted_Y = point)
    mixedgp_atomic_csv(predictions, oai_path(dir, "Y_predictions.csv"))
    if (any(z$available)) {
        if (!isTRUE(fit$data$latent_scale_anchored))
            stop("Unanchored fit cannot receive raw test U.")
        sel <- which(z$available)
        u <- as.matrix(te[sel, z$u, drop = FALSE])
        colnames(u) <- p$U_order
        exact <- predict_eivgp(fit, new_X = as.matrix(te[sel, z$x, drop = FALSE]), new_C = ct[sel, , drop = FALSE],
            new_U = u, new_U_scale = "raw", target = "response", draw_ids = ids, n_per_draw = 1L, joint = FALSE,
            seed = p$seeds$mixed_prediction)
        mixed <- as.matrix(yd)
        mixed[, sel] <- as.matrix(exact)
        cm <- attr(exact, "conditional_means", exact = TRUE)
        mp <- point
        mp[sel] <- if (is.null(cm))
            colMeans(as.matrix(exact))
        else colMeans(cm)
        metrics <- rbind(metrics, score(mixed, mp, "MI"))
        mixedgp_atomic_csv(data.frame(row = seq_len(nrow(te)), observed_Y = te$Y, predicted_Y = mp, U_available = z$available),
            oai_path(dir, "mixed_Y_predictions.csv"))
    }
    mixedgp_atomic_csv(metrics, oai_path(dir, "prediction_metrics.csv"))
    dg <- oai_plan_diagnostics(fit)
    mixedgp_atomic_csv(cbind(fold = task$fold, cell = task$cell, dg$summary), oai_path(dir, "convergence.csv"))
    mixedgp_atomic_csv(dg$parameters, oai_path(dir, "diagnostic_parameters.csv"))
    scale <- if (any(tr$calibration == 1L))
        "raw"
    else "model"
    imp <- impute_eivgp(fit, rows = seq_len(nrow(tr)), draw_ids = ids, scale = scale)
    mixedgp_atomic_rds(list(draws = imp, scale = scale, draw_ids = ids, coordinates = p$U_order), oai_path(dir,
        "U_draws.rds"))
    if (scale == "raw") {
        plots <- tables <- list()
        hidden <- which(tr$calibration == 0)
        cal <- which(tr$calibration == 1)
        common <- inputs$cells[[paste(task$fold, "O50", sep = "/")]]$train$calibration == 0
        for (j in seq_along(z$u)) {
            ud <- matrix(imp[, , j], nrow = length(ids))
            up <- colMeans(ud)
            qq <- apply(ud, 2, quantile, c(0.025, 0.975))
            truth <- tr[[z$u[j]]]
            sd0 <- sd(truth)
            tables[[j]] <- cbind(data.frame(fold = task$fold, cell = task$cell, coordinate = p$U_order[j]), oai_plan_score(ud[,
                hidden, drop = FALSE], up[hidden], truth[hidden], sd0), calibration_mean_RMSE = sqrt(mean((mean(truth[cal]) -
                truth[hidden])^2)))
            plots[[j]] <- data.frame(fold = task$fold, cell = task$cell, coordinate = p$U_order[j], row = seq_len(nrow(tr)),
                calibration = tr$calibration, common_hidden = common, true_U = truth, posterior_mean = up, q025 = qq[1,
                  ], q975 = qq[2, ])
        }
        mixedgp_atomic_csv(do.call(rbind, tables), oai_path(dir, "U_metrics.csv"))
        mixedgp_atomic_csv(do.call(rbind, plots), oai_path(dir, "U_predictions.csv"))
    }
    if (task$fold == s$folds[1]) {
        series <- mixedgp_v030_raw_series(fit)
        parameters <- head(grep("^(V$|r$|theta_x|tau|A\\[)", names(series), value = TRUE), 6L)
        if (length(parameters) && length(fit$checkpoint$states) >= 2L)
            ggplot2::ggsave(oai_path(dir, "mcmc_trace.png"), plot_mcmc(fit, parameters), width = 10, height = 7,
                dpi = 130)
        else if (length(fit$checkpoint$states) < 2L)
            writeLines("Single-chain workflow check: multi-chain trace diagnostics are unavailable; convergence is not established.",
                oai_path(dir, "MCMC_TRACE_UNAVAILABLE.txt"))
    }
    invisible("DONE")
}

# OAI application: evaluation. Reuses package model and sampler functions.

oai_save_plot <- function(filename, plot, width, height, dpi = 130, ...) {
    input <- paste0(filename, ".plot.rds")
    saveRDS(plot, input)
    code <- paste0("source(", deparse(oai_path(repository_dir,"tools/oai.R")),
      ");oai_load_api()$oai_render_plot(commandArgs(TRUE))")
    args <- c('-e', shQuote(code), shQuote(input), shQuote(filename), width, height, dpi)
    rc <- system2(oai_path(R.home("bin"), "Rscript"), args, stdout = paste0(filename, ".render.log"), stderr = paste0(filename,
        ".render.log"))
    stopifnot(rc == 0L, file.exists(filename), file.info(filename)$size > 1000)
    unlink(input)
    invisible(filename)
}

oai_replace_plot_device <- function(x) {
    if (is.call(x) && identical(x, quote(ggplot2::ggsave)))
        return(as.name("oai_save_plot"))
    if (is.call(x))
        for (i in seq_along(x)) x[[i]] <- oai_replace_plot_device(x[[i]])
    x
}

oai_reference_points <- function(calibration, targets, c_names, u_name) {
    global <- rep(mean(calibration[[u_name]]), nrow(targets))
    c1 <- tapply(calibration[[u_name]], calibration$C1, mean)
    one <- as.numeric(c1[as.character(targets$C1)])
    one[!is.finite(one)] <- global[!is.finite(one)]
    key <- function(d) apply(d[, c_names, drop = FALSE], 1, paste, collapse = "/")
    joint <- tapply(calibration[[u_name]], key(calibration), mean)
    jp <- as.numeric(joint[key(targets)])
    fallback <- !is.finite(jp)
    jp[fallback] <- one[fallback]
    list(global = global, C1 = one, joint = jp, fallback = sum(fallback))
}

oai_analyze_cell <- function(task, inputs, out, s, stage) {
    z <- inputs$cells[[paste(task$fold, task$cell, sep = "/")]]
    tr <- z$train
    te <- z$test
    p <- z$protocol
    stageout <- oai_path(out, stage)
    dir <- oai_path(stageout, "analysis", sprintf("fold%02d", task$fold), task$cell)
    dir.create(dir, recursive = TRUE, showWarnings = FALSE)
    cache <- new.env(parent = emptyenv())
    env <- new.env(parent = environment(api$oai_plan_analyze_cell))
    env$oai_plan_fit_path <- function(ignored, task) api$oai_plan_fit_path(out, task)
    env$oai_plan_read_fit <- function(path) {
        cache$fit <- api$oai_plan_read_fit(path)
        stopifnot(cache$fit$checkpoint$iteration == s$iter)
        cache$fit
    }
    env$predict_eivgp <- function(...) {
        a <- list(...)
        r <- eivGP::predict_eivgp(...)
        if (is.null(a$new_U)) {
            stopifnot(is.null(a$new_y))
            cache$ids <- a$draw_ids
            cache$pure <- r
            atomic_rds(list(draws = r, draw_ids = a$draw_ids), oai_path(dir, "Y_draws.rds"))
        }
        else {
            cache$exact <- r
            atomic_rds(list(draws = r, draw_ids = a$draw_ids), oai_path(dir, "known_U_Y_draws.rds"))
        }
        r
    }
    env$oai_save_plot <- oai_save_plot
    fun <- api$oai_plan_analyze_cell
    body(fun) <- oai_replace_plot_device(body(fun))
    environment(fun) <- env
    fun(task, inputs, stageout, s)
    fit <- cache$fit
    stopifnot(fit$data$p == length(z$x), fit$data$q == p$q, fit$data$d == length(z$u), all(is.na(fit$data$U_obs_raw[tr$calibration ==
        0, , drop = FALSE])))
    met <- read.csv(oai_path(dir, "prediction_metrics.csv"))
    met$stage <- stage
    met$iterations <- s$iter
    met$sampler <- fit$control$sampler
    if (exists("exact", cache, inherits = FALSE)) {
        ids <- which(z$available)
        cm <- attr(cache$exact, "conditional_means", exact = TRUE)
        stopifnot(!is.null(cm))
        extra <- cbind(data.frame(fold = task$fold, cell = task$cell, method = "EIV-GP", protocol = "CC-known"),
            api$oai_plan_score(cache$exact, colMeans(cm), te$Y[ids], sd(tr$Y)))
        extra$stage <- stage
        extra$iterations <- s$iter
        extra$sampler <- fit$control$sampler
        met <- rbind(met, extra)
    }
    atomic_csv(met, oai_path(dir, "prediction_metrics.csv"))
    dg <- read.csv(oai_path(dir, "convergence.csv"))
    dg$stage <- stage
    dg$sampler <- fit$control$sampler
    if (!is.null(fit$diagnostics$hmc))
        atomic_csv(fit$diagnostics$hmc, oai_path(dir, "HMC_DIAGNOSTICS.csv"))
    atomic_csv(dg, oai_path(dir, "convergence.csv"))
    cm <- attr(cache$pure, "conditional_means", exact = TRUE)
    per <- nrow(cm)%/%s$chains
    target_diagnostics <- do.call(rbind, lapply(seq_len(ncol(cm)), function(j) {
        a <- matrix(cm[, j], per, s$chains)
        data.frame(RID = te$RID[j], predicted_Y = mean(a), selected_Rhat = posterior::rhat(a), selected_ESS = posterior::ess_bulk(a),
            selected_MCSE = posterior::mcse_mean(a))
    }))
    atomic_csv(target_diagnostics, oai_path(dir, "PREDICTION_FUNCTION_DIAGNOSTICS.csv"))
    if (task$cell != "O0") {
        cal <- tr[tr$calibration == 1, , drop = FALSE]
        hidden <- which(tr$calibration == 0)
        ui <- read.csv(oai_path(dir, "U_predictions.csv"))
        um <- read.csv(oai_path(dir, "U_metrics.csv"))
        for (j in seq_along(z$u)) {
            truth <- tr[[z$u[j]]]
            ref <- oai_reference_points(cal, tr[hidden, , drop = FALSE], z$c, z$u[j])
            row <- which(um$coordinate == p$U_order[j])
            um$C1_mean_RMSE[row] <- sqrt(mean((ref$C1 - truth[hidden])^2))
            um$joint_C_mean_RMSE[row] <- sqrt(mean((ref$joint - truth[hidden])^2))
            um$joint_C_fallback_N[row] <- ref$fallback
            ui$C1[ui$coordinate == p$U_order[j]] <- tr$C1
            hit <- ui$coordinate == p$U_order[j]
            ui$true_std[hit] <- (ui$true_U[hit] - mean(truth))/sd(truth)
            ui$mean_std[hit] <- (ui$posterior_mean[hit] - mean(truth))/sd(truth)
            ui$lo_std[hit] <- (ui$q025[hit] - mean(truth))/sd(truth)
            ui$hi_std[hit] <- (ui$q975[hit] - mean(truth))/sd(truth)
        }
        um$stage <- stage
        atomic_csv(um, oai_path(dir, "U_metrics.csv"))
        atomic_csv(ui, oai_path(dir, "U_predictions.csv"))
        ct <- as.matrix(te[, z$c, drop = FALSE])
        colnames(ct) <- p$C_order
        prospective <- impute_eivgp(fit, new_C = ct, draw_ids = cache$ids, scale = "raw", seed = p$seeds$prospective_U)
        stopifnot(identical(dim(prospective), c(length(cache$ids), nrow(te), length(z$u))))
        pm <- list()
        pp <- list()
        for (j in seq_along(z$u)) {
            ud <- matrix(prospective[, , j], nrow = length(cache$ids))
            point <- colMeans(ud)
            truth <- te[[z$u[j]]]
            ref <- oai_reference_points(cal, te, z$c, z$u[j])
            qq <- apply(ud, 2, quantile, c(0.025, 0.975))
            pm[[j]] <- cbind(data.frame(stage = stage, fold = task$fold, cell = task$cell, coordinate = p$U_order[j]),
                api$oai_plan_score(ud, point, truth, sd(tr[[z$u[j]]])), calibration_mean_RMSE = sqrt(mean((ref$global -
                  truth)^2)), C1_mean_RMSE = sqrt(mean((ref$C1 - truth)^2)), joint_C_mean_RMSE = sqrt(mean((ref$joint -
                  truth)^2)), joint_C_fallback_N = ref$fallback)
            pp[[j]] <- data.frame(RID = te$RID, coordinate = p$U_order[j], true_U = truth, posterior_mean = point,
                q025 = qq[1, ], q975 = qq[2, ])
        }
        atomic_csv(do.call(rbind, pm), oai_path(dir, "PROSPECTIVE_U_metrics.csv"))
        atomic_csv(do.call(rbind, pp), oai_path(dir, "PROSPECTIVE_U_predictions.csv"))
        ui$role <- factor(ifelse(ui$calibration == 1, "Observed calibration", "Masked training U"), levels = c("Masked training U",
            "Observed calibration"))
        plot <- ggplot(ui, aes(true_std, mean_std, colour = factor(C1))) + geom_abline(slope = 1, intercept = 0,
            linetype = 2, colour = "grey50") + geom_linerange(data = ui[ui$calibration == 0, ], aes(ymin = lo_std,
            ymax = hi_std), alpha = 0.15, linewidth = 0.3) + geom_point(aes(shape = role), size = 1.2, alpha = 0.65) +
            scale_shape_manual(values = c(16, 1)) + facet_wrap(~coordinate) + theme_bw(base_size = 11) + theme(legend.position = "bottom") +
            labs(x = "Reference U (training SD units; evaluation only)", y = "Posterior mean (training SD units)",
                colour = p$C_order[1], shape = NULL, title = paste(p$plan_id, "fold", task$fold, task$cell, stage),
                subtitle = paste(fit$control$sampler, "; CORE pass:", dg$core_pass))
        oai_save_plot(oai_path(dir, "training_U_imputation_std.png"), plot, width = 10, height = 5.3, dpi = 140,
            bg = "white")
    }
    else writeLines("O0 is unanchored: no physical-U RMSE or truth-based alignment.", oai_path(dir, "O0_SCALE_LIMITATION.txt"))
    atomic_json(list(state = "PASS", stage = stage, iterations = s$iter, primary_no_test_U_Y = TRUE, hidden_training_U_masked = TRUE,
        calibration_excluded_from_U_scoring = TRUE, prospective_U_uses_C_not_test_Y = TRUE, sampler = fit$control$sampler),
        oai_path(dir, "EVALUATION_VALIDATION.json"))
    rm(fit, cache)
    gc(FALSE)
    invisible("DONE")
}

oai_baseline_cell <- function(task, inputs, out, s, cc = FALSE) {
    source <- inputs$cells[[paste(task$fold, if (cc)
        task$cell
    else "O50", sep = "/")]]
    tr <- source$train
    te <- source$test
    root <- if (cc)
        oai_path(out, "CC", task$cell)
    else out
    inp <- inputs
    if (cc) {
        z <- source
        z$train <- tr[tr$calibration == 1, , drop = FALSE]
        z$test <- te
        for (j in seq_along(z$u)) {
            name <- paste0("X", length(z$x) + j)
            z$train[[name]] <- z$train[[z$u[j]]]
            z$test[[name]] <- z$test[[z$u[j]]]
        }
        z$x <- paste0("X", seq_len(length(z$x) + length(z$u)))
        for (n in c("UC_GP", "LVGP", "EzGP", "EzGP_draws")) z$protocol$seeds[[n]] <- as.integer((as.double(z$protocol$seeds[[n]]) +
            if (task$cell == "O20") 17041 else 35023)%%(2^30 - 1) + 1)
        inp$cells[[paste(task$fold, "O50", sep = "/")]] <- z
    }
    api$oai_plan_baseline(task, inp, root, s)
    dir <- oai_path(root, "baselines", sprintf("fold%02d", task$fold), task$method)
    saved <- readRDS(oai_path(dir, "prediction_draws.rds"))
    selected <- if (cc)
        which(source$available)
    else seq_len(nrow(te))
    mm <- cbind(data.frame(stage = "BASELINE", fold = task$fold, cell = if (cc)
        task$cell
    else "OI", method = paste0(task$method, if (cc)
        "-CC"
    else ""), protocol = if (cc)
        "CC-known"
    else "OI"), api$oai_plan_score(saved$draws[, selected, drop = FALSE], saved$point[selected], te$Y[selected],
        sd(tr$Y)))
    mm$training_N <- if (cc)
        sum(tr$calibration == 1)
    else nrow(tr)
    atomic_csv(mm, oai_path(dir, "metrics.csv"))
    atomic_csv(data.frame(RID = te$RID, observed_Y = te$Y, predicted_Y = saved$point, scored = seq_len(nrow(te)) %in%
        selected), oai_path(dir, "Y_predictions.csv"))
    invisible("DONE")
}

# OAI application: execution. Reuses package model and sampler functions.

oai_with_lock <- function(path,FUN,wait=FALSE,timeout=60) {
  dir.create(dirname(path),recursive=TRUE,showWarnings=FALSE)
  started<-proc.time()[['elapsed']];announced<-FALSE
  repeat {
    if(dir.create(path,showWarnings=FALSE))break
    if(!dir.exists(path))stop('OAI cannot create lock; check directory permissions and free space: ',path)
    owner<-base::file.path(path,'owner.json')
    info<-tryCatch(if(file.exists(owner))jsonlite::fromJSON(owner)else NULL,error=function(e)NULL)
    detail<-if(is.null(info))'owner record not yet available'else paste('PID',info$pid,'on',info$host)
    if(!wait||proc.time()[['elapsed']]-started>=timeout)
      stop('OAI lock is held (',detail,'): ',path,
        '. Check the owning process; remove a stale lock only after its workers have stopped.')
    if(!announced){cat('[OAI] Waiting briefly for shared initialization: ',detail,'\n',sep='');flush.console();announced<-TRUE}
    Sys.sleep(.1)
  }
  token<-paste(Sys.getpid(),format(Sys.time(),'%Y%m%d%H%M%OS6'),sep='-')
  on.exit({
    owner<-base::file.path(path,'owner.json')
    saved<-tryCatch(jsonlite::fromJSON(owner),error=function(e)NULL)
    if(!is.null(saved)&&identical(saved$token,token))unlink(path,recursive=TRUE)
  },add=TRUE)
  atomic_json(list(pid=Sys.getpid(),host=Sys.info()[['nodename']],token=token,
    started_at=format(Sys.time(),tz='UTC',usetz=TRUE)),base::file.path(path,'owner.json'))
  FUN()
}

oai_check_legacy_lock <- function(root) {
  old<-unlist(lapply(unique(c(root,oai_workflow_root(root))),function(path)
    base::file.path(path,c('.application.lock','.adni-application.lock','.oai-application.lock'))),use.names=FALSE)
  if(any(dir.exists(old)))stop('OAI found a lock from the earlier exclusive executor: ',old[dir.exists(old)][1],
    '. Finish or stop that run before using the new executor; confirm all workers have stopped before removing a stale lock.')
  invisible(NULL)
}

oai_execution_lock <- function(root,group,FUN) {
  oai_check_legacy_lock(root)
  if(group=='report')return(oai_execution_lock(root,'eivGP',function()
    oai_execution_lock(root,'competitors',FUN)))
  oai_with_lock(base::file.path(root,'_meta','locks',paste0(group,'.lock')),FUN)
}

oai_publish_progress <- function(root,group,phase,status) {
  active<-status[status$state=='RUNNING',,drop=FALSE]
  record<-list(application='OAI',group=group,phase=phase,pid=Sys.getpid(),
    updated_at=format(Sys.time(),tz='UTC',usetz=TRUE),completed=sum(status$state%in%c('DONE','REUSED')),
    failed=sum(status$state=='FAILED'),total=nrow(status),active_jobs=active)
  atomic_json(record,base::file.path(root,'_meta',paste0(group,'_status.json')))
  lines<-c(paste('OAI |',group,'|',phase,'|',record$updated_at),
    sprintf('Completed %d/%d; failed %d; running %d',record$completed,record$total,record$failed,nrow(active)),
    if(nrow(active))apply(active[,c('fold','cell','method','state'),drop=FALSE],1,paste,collapse=' | ')else 'No job currently running.',
    'Detailed task logs and phase status CSVs are under logs/.')
  path<-base::file.path(root,paste0('PROGRESS_',group,'.txt'))
  tmp<-tempfile(tmpdir=dirname(path));on.exit(unlink(tmp),add=TRUE);writeLines(lines,tmp)
  if(!file.rename(tmp,path))stop('OAI cannot publish progress: ',path)
}

oai_workflow_identity <- function(path,id) {
  if(!file.exists(path)){atomic_rds(id,path);return(id)}
  saved<-readRDS(path)
  if(identical(saved,id))return(saved)
  a<-saved;b<-id;a$application<-b$application<-NULL;a$schema<-b$schema<-NULL
  if(length(saved$schema)!=1L||!saved$schema%in%c('adni-upstream-screening-v2','oai-upstream-screening-v2')||!identical(a,b))
    stop('OAI saved model settings, package modules or inputs differ; use a new output root: ',dirname(path))
  # Preserve the original fitting identity and receipts. Record only the new
  # execution source; sampler/package/input identities must match exactly.
  dest<-base::file.path(dirname(path),'executor_updates',paste0(digest::digest(id$application),'.rds'))
  if(!file.exists(dest))atomic_rds(list(schema='oai-executor-update-v1',
    original_identity=saved,current_execution_source=id$application,
    package_and_statistical_identity_unchanged=TRUE),dest)
  cat('[OAI] Reusing the original saved fitting identity; execution-source update recorded.\n');flush.console()
  saved
}

oai_output_check <- function(path,repo) {
  path<-normalizePath(path.expand(path),mustWork=FALSE);repo<-normalizePath(repo)
  if(startsWith(paste0(path,'/'),paste0(repo,'/'))&&
     !startsWith(paste0(path,'/'),paste0(repo,'/results/')))
    stop('Choose ignored results/ or an external output directory; never write to versioned source/data.')
  invisible(path)
}

oai_screen_args <- function(args) {
  p<-list(cores=5L,plan=FALSE,reuse_only=FALSE)
  for(a in args){
    if(!grepl('^--[^=]+=.+$',a))stop('Use --name=value: ',a)
    k<-gsub('-','_',sub('^--([^=]+)=.*$','\\1',a));v<-sub('^[^=]+=','',a)
    if(!k%in%c('manifest','cores','plan','reuse_only'))stop('Unknown screening option: ',k)
    if(k%in%c('plan','reuse_only')){if(!v%in%c('true','false'))stop(k,' must be true or false');v<-v=='true'}
    if(k%in%c('cores')){n<-suppressWarnings(as.numeric(v));if(!is.finite(n)||n<1||n!=floor(n)||n>.Machine$integer.max)stop('cores must be a positive integer');v<-as.integer(n)}
    if(k%in%names(p)&&!k%in%c('cores','plan','reuse_only'))stop('Repeated option: ',k)
    p[[k]]<-v
  }
  if(is.null(p$manifest))stop('Supply --manifest=/external/SCREENING.json')
  p
}

oai_screen_settings <- function(m,iter,cores) {
  smoke<-isTRUE(m$smoke)
  get<-function(n,default)if(is.null(m[[n]]))default else m[[n]]
  list(profile=get('baseline_profile',if(smoke)'pilot'else'paper'),run_profile=get('profile',if(smoke)'quick'else'pilot'),
    action='fit',folds=as.integer(unlist(m$folds)),cells=unlist(m$cells),
    cores=as.integer(cores),chains=as.integer(get('chains',4L)),iter=as.integer(iter),burn=as.integer(m$burn),
    first_checkpoint=as.integer(min(iter,get('first_checkpoint',if(smoke)iter else m$burn+250L))),
    checkpoint_step=as.integer(get('checkpoint_step',1000L)),
    ess_min=if(iter==max(unlist(m$stages)))mixedgp_parse_ess_min(m$ess_min)else NULL,
    final_stage=iter==max(unlist(m$stages)),
    draws_per_chain=as.integer(get('draws_per_chain',if(smoke)5L else 250L)),protocol=m$trial_id,
    samplers=get('samplers',list(O0='HMC',O20='HMC',O50='ESS')))
}

oai_screen_manifest <- function(path,api) {
  m<-jsonlite::fromJSON(normalizePath(path),simplifyVector=FALSE)
  if(!m$schema%in%c('oai-screening-v2','oai-application-v1','adni-screening-v2','adni-application-v1')||!length(m$cases))stop('Unsupported OAI manifest')
  if(m$schema=='oai-screening-v2'&&!isTRUE(m$smoke)&&
    (m$burn!=500L||!identical(as.integer(unlist(m$folds)),1:5)))stop('Production screening requires all five folds and burn500')
  m$ess_min<-mixedgp_parse_ess_min(m$ess_min)
  if(!is.null(m$ess_min) && !is.null(m$chains) && m$chains<4L)stop('Scientific ESS requires four chains.')
  stages<-unlist(m$stages)
  if(!length(stages)||any(!is.finite(stages))||any(stages!=as.integer(stages))||any(diff(stages)<=0)||
    length(m$burn)!=1L||!is.finite(m$burn)||m$burn<0||m$burn!=as.integer(m$burn)||min(stages)<=m$burn)stop('Invalid stage budgets')
  if(!length(m$folds)||any(!unlist(m$folds)%in%1:5)||anyDuplicated(unlist(m$folds)))stop('Invalid folds')
  if(!is.null(m$chains)&&(length(m$chains)!=1L||m$chains<1||m$chains!=as.integer(m$chains)))stop('Invalid chain count')
  if(!all(unlist(m$cells)%in%c('O0','O20','O50')))stop('Invalid calibration cells')
  if(anyDuplicated(vapply(m$cases,`[[`,'','id')))stop('Repeated Plan identity')
  for(k in names(m$frozen_inputs))if(!identical(digest::digest(file=k,algo='sha256'),m$frozen_inputs[[k]]))stop('Frozen input changed: ',k)
  inputs<-list()
  for(z in m$cases){
    if(!grepl('^[A-Za-z0-9][A-Za-z0-9_-]*$',z$id))stop('Invalid dataset identity')
    oai_output_check(z$output,api$repository_dir)
    inp<-api$oai_plan_inputs(z$input)
    for(cell in inp$cells)if(cell$protocol$plan_id!=z$id||any(cell$test$calibration!=0L))stop('Dataset identity or test mask differs')
    inputs[[z$id]]<-inp
  }
  if(length(inputs)>1L)for(i in 2:length(inputs))for(key in names(inputs[[1]]$cells)){
    a<-inputs[[1]]$cells[[key]];b<-inputs[[i]]$cells[[key]]
    if(!identical(a$train$RID,b$train$RID)||!identical(a$test$RID,b$test$RID)||
       !identical(a$train[,c('Y','X1','X2','U1','U2','calibration')],b$train[,c('Y','X1','X2','U1','U2','calibration')])||
       !identical(a$available,b$available))stop('Paired people, anchors or masks differ')
  }
  list(manifest=m,inputs=inputs)
}

oai_screen_tasks <- function(m,group,stage='') {
  rows<-list();add<-function(plan,fold,cell,method,kind)rows[[length(rows)+1L]]<<-data.frame(plan=plan,fold=as.integer(fold),cell=cell,method=method,kind=kind,stage=stage)
  for(z in m$cases)for(f in unlist(m$folds)){
    if(group=='eivGP')for(cell in unlist(m$cells))add(z$id,f,cell,'EIV-GP','eiv')
    if(group=='competitors'){
      for(method in c('UC-GP','LVGP','EzGP'))add(z$id,f,'OI',method,'baseline')
      for(cell in intersect(unlist(m$cells),c('O20','O50')))for(method in c('UC-GP','LVGP','EzGP'))add(z$id,f,cell,method,'cc')
    }
  }
  if(!length(rows))return(data.frame())
  x<-do.call(rbind,rows);x$id<-paste(x$plan,x$stage,x$kind,x$fold,x$cell,gsub('[+]','_',x$method),sep='__');x
}

oai_screen_receipt <- function(out,task,phase)oai_path(out,'receipts',paste0(phase,'_',task$id,'.rds'))

oai_screen_verified <- function(path,identity) {
  if(!file.exists(path))return(FALSE)
  r<-readRDS(path);if(!identical(r$identity,identity))stop('Receipt identity differs: ',path)
  if(length(r$files)&&(!all(file.exists(r$files))||!identical(unname(tools::md5sum(r$files)),r$md5)))stop('Completed output changed: ',path)
  TRUE
}

oai_resident_workers <- function(cores,phase,n_tasks,chains=4L) {
  phase<-match.arg(phase,c('fit','evaluate'))
  cores<-oai_integer(cores,'cores');chains<-oai_integer(chains,'chains')
  n_tasks<-oai_integer(n_tasks,'n_tasks',0L)
  # One resident job per chain-sized group in either phase; the shared pool
  # lends unused CPU slots to other jobs, including a complete group.
  as.integer(min(n_tasks,mixedgp_stage_workers(cores,phase,n_tasks,chains)))
}

oai_screen_phase <- function(tasks,phase,cores,api,FUN,status_dir,chains=4L,root=NULL,group='eivGP') {
  dir.create(status_dir,recursive=TRUE,showWarnings=FALSE)
  status<-tasks;status$state<-'PENDING';status$detail<-'';status$worker_pid<-NA_integer_;status$seconds<-NA_real_
  outer<-oai_resident_workers(cores,phase,nrow(tasks),chains)
  publish<-function() {
    api$atomic_csv(status,oai_path(status_dir,paste0(phase,'_status.csv')))
    if(!is.null(root))oai_publish_progress(root,group,phase,status)
  }
  publish()
  cat(sprintf('[OAI %s] Phase %s: %d jobs; %d cores; up to %d resident jobs.\n',group,phase,nrow(tasks),cores,outer));flush.console()
  # Match the public Ocean executor: resident jobs share the full active-worker pool.
  api$mixedgp_budget_yield(function()api$queue$mixedgp_schedule(seq_len(nrow(tasks)),rep(1L,nrow(tasks)),outer,
    function(i,allocated)api$mixedgp_budget_worker(function(){
      logpath<-oai_path(status_dir,paste0(phase,'_',make.names(tasks$id[i]),'.log'));con<-file(logpath,'at')
      old_message<-getConnection(sink.number(type='message'));sink(con);sink(con,type='message')
      on.exit({sink(old_message,type='message');sink();close(con)},add=TRUE);t<-proc.time()[3]
      result<-NULL
      err<-tryCatch({result<-FUN(tasks[i,,drop=FALSE],min(cores,chains));NULL},error=function(e)conditionMessage(e))
      if(!is.null(err))cat('FAILED:',err,'\n')
      list(state=if(!is.null(err))'FAILED'else if(identical(result,'REUSED'))'REUSED'else'DONE',
        detail=if(is.null(err))''else err,seconds=as.numeric(proc.time()[3]-t))
    }),
    function(i,n,pid,active){
      status$state[i]<<-'RUNNING';status$worker_pid[i]<<-pid;publish()
      cat(sprintf('[OAI %s] START %s | fold %d | %s | %s | PID %d\n',group,phase,tasks$fold[i],tasks$cell[i],tasks$method[i],pid));flush.console()
    },
    function(i,r,active){if(!is.list(r)||is.null(r$state))r<-list(state='FAILED',detail='Worker exited without a result',seconds=NA_real_)
      status$state[i]<<-r$state;status$detail[i]<<-r$detail;status$seconds[i]<<-r$seconds;publish()
      cat(sprintf('[OAI %s] %s %s | fold %d | %s | %s | %.1f sec %s\n',group,r$state,phase,tasks$fold[i],tasks$cell[i],tasks$method[i],r$seconds,r$detail));flush.console()
    }))
  status
}

oai_screen_collect <- function(m,api,group='report') {
  pred<-conv<-u<-pros<-list()
  for(z in m$cases){
    for(mark in if(group=='competitors')character()else list.files(z$output,'^EVALUATION_VALIDATION[.]json$',recursive=TRUE,full.names=TRUE)){
      dir<-dirname(mark)
      for(pair in list(c('prediction_metrics.csv','pred'),c('convergence.csv','conv'),c('U_metrics.csv','u'),c('PROSPECTIVE_U_metrics.csv','pros'))){
        p<-oai_path(dir,pair[1]);if(!file.exists(p))next
        val<-cbind(plan=z$id,read.csv(p));dest<-get(pair[2]);dest[[length(dest)+1L]]<-val;assign(pair[2],dest)
      }
    }
    files<-if(group=='eivGP')character()else list.files(z$output,'^metrics[.]csv$',recursive=TRUE,full.names=TRUE)
    for(p in files){x<-read.csv(p);if('stage'%in%names(x)&&all(x$stage=='BASELINE'))pred[[length(pred)+1L]]<-cbind(plan=z$id,x)}
  }
  bind<-function(x)if(length(x))api$queue$mixedgp_bind(x)else data.frame()
  tables<-list(prediction=bind(pred),convergence=bind(conv),training_U=bind(u),prospective_U=bind(pros))
  dest<-oai_path(m$output_root,'reports',if(group=='report')'prediction_report'else group,'progress')
  dir.create(dest,recursive=TRUE,showWarnings=FALSE)
  for(k in names(tables))if(nrow(tables[[k]]))api$atomic_csv(tables[[k]],oai_path(dest,paste0(k,'_fold_metrics.csv')))
  p<-tables$prediction;summary<-list()
  if(nrow(p))for(ix in split(seq_len(nrow(p)),interaction(p$plan,p$stage,p$cell,p$method,p$protocol,drop=TRUE))){
    x<-p[ix,,drop=FALSE];v<-x$RMSE_SD
    summary[[length(summary)+1L]]<-data.frame(x[1,c('plan','stage','cell','method','protocol')],completed_folds=length(unique(x$fold)),
       RMSE_SD_mean=mean(v),RMSE_SD_SD=if(length(v)>1)sd(v)else NA_real_)
  }
  five<-bind(summary);if(nrow(five))api$atomic_csv(five,oai_path(dest,'fivefold_prediction.csv'))
  stability<-list()
  for(kind in names(tables)[names(tables)%in%c('prediction','training_U','prospective_U')]){
    tab<-tables[[kind]];if(!nrow(tab)||!'stage'%in%names(tab))next
    stage_names<-paste0('S',unlist(m$stages));if(length(stage_names)<2L)next
    x<-tab[tab$stage==stage_names[1],,drop=FALSE];y<-tab[tab$stage==tail(stage_names,1),,drop=FALSE]
    keys<-if(kind=='prediction')c('plan','fold','cell','method','protocol')else c('plan','fold','cell','coordinate')
    if(nrow(x)&&nrow(y)){
      z<-merge(x[,c(keys,'RMSE_SD'),drop=FALSE],y[,c(keys,'RMSE_SD'),drop=FALSE],by=keys,suffixes=c('_first','_last'))
      z$first_stage<-stage_names[1];z$last_stage<-tail(stage_names,1)
      z$target<-kind;z$delta_last_minus_first<-z$RMSE_SD_last-z$RMSE_SD_first;stability[[kind]]<-z
    }
  }
  stable<-bind(stability);if(nrow(stable))api$atomic_csv(stable,oai_path(dest,'stage_stability.csv'))
  lines<-c('# OAI screening: saved results',paste('Author: Xinye Chen. Trial',m$trial_id),
    'Partial fold counts and CORE failures remain explicit. Budget completion does not establish convergence.',
    if(nrow(five))capture.output(print(five,row.names=FALSE))else'No completed scores yet.')
  writeLines(lines,oai_path(dest,'REPORT.md'))
  html<-c('<!doctype html><meta charset="utf-8"><title>OAI screening</title>',
    '<h1>OAI screening: saved results</h1><p>Author: Xinye Chen. Partial folds and CORE failures remain explicit; completing a budget does not establish convergence.</p>',
    api$queue$mixedgp_table_html(five,'Response prediction: fold means and sample SD'),
    api$queue$mixedgp_table_html(tables$convergence,'Full-chain CORE diagnostics'),
    api$queue$mixedgp_table_html(stable,'Paired stage stability; descriptive screening only'))
  writeLines(html,oai_path(dest,'index.html'));invisible(tables)
}

oai_screen_run <- function(repo,p,group) {
  api<-oai_screen_api(repo);bundle<-oai_screen_manifest(p$manifest,api);m<-bundle$manifest;inputs<-bundle$inputs
  if(length(p$run_folds)) {
    selected<-oai_integer_list(p$run_folds,'run-folds')
    if(anyDuplicated(selected)||any(!selected%in%unlist(m$folds)))stop('run-folds must be distinct members of the frozen folds')
    m$folds<-selected
  }
  cases<-setNames(m$cases,vapply(m$cases,`[[`,'','id'));api$B<-dirname(normalizePath(p$manifest))
  if(p$plan){
    for(g in if(group=='report')c('eivGP','competitors')else group){x<-oai_screen_tasks(m,g);print(x,row.names=FALSE)}
    cat('Public shared CPU pool:',p$cores,'; fit resident limit:',oai_resident_workers(p$cores,'fit',100L,if(is.null(m$chains))4L else m$chains),
      '; evaluation resident limit:',oai_resident_workers(p$cores,'evaluate',100L,if(is.null(m$chains))4L else m$chains),'; full-chain reservations/backfill disabled.\n');return(invisible(bundle))
  }
  if(group=='eivGP' && !isTRUE(p$reuse_only) &&
     (if(is.null(m$chains))4L else m$chains)!=mixedgp_application_allocation(p$cores)$chains)
    stop('Saved/manifest chains must equal min(--cores,4); use a matching core budget or a new output.')
  root<-m$output_root
  oai_output_check(root,repo)
  if(!isTRUE(p$lock_held))return(oai_execution_lock(root,group,function(){
    p$lock_held<-TRUE;oai_screen_run(repo,p,group)
  }))
  cat('[OAI ',group,'] BEGIN | output: ',root,'\n',sep='');flush.console()
  code<-oai_source_files(repo)
  base_id<-list(schema='oai-upstream-screening-v2',trial=m$trial_id,modules=api$eivgp_module_identity(),
    application=setNames(unname(tools::md5sum(code)),substring(code,nchar(repo)+2L)),burn=m$burn,stages=unlist(m$stages),smoke=isTRUE(m$smoke),ess_min=m$ess_min,
    chains=if(is.null(m$chains))4L else as.integer(m$chains),
    samplers=if(is.null(m$samplers))list(O0='HMC',O20='HMC',O50='ESS')else m$samplers,
    baseline_profile=if(is.null(m$baseline_profile))if(isTRUE(m$smoke))'pilot'else'paper'else m$baseline_profile,
    draws_per_chain=if(is.null(m$draws_per_chain))if(isTRUE(m$smoke))5L else 250L else m$draws_per_chain)
  identities<-list()
  oai_with_lock(base::file.path(dirname(root),paste0('.',basename(root),'.initialization.lock')),function(){
    for(z in m$cases){
      id<-c(base_id,list(plan=z$id,inputs=unname(tools::md5sum(inputs[[z$id]]$files))))
      path<-oai_path(z$output,'workflow_identity.rds')
      identities[[z$id]]<<-oai_workflow_identity(path,id)
    }
  },wait=TRUE)
  required<-list()
  for(g in if(group=='report')c('eivGP','competitors')else group){
    for(st in if(g=='eivGP')paste0('S',unlist(m$stages))else'BASELINE'){
      tasks<-oai_screen_tasks(m,g,st)
      for(i in seq_len(nrow(tasks))){t<-tasks[i,];out<-cases[[t$plan]]$output
        required[[length(required)+1L]]<-list(path=oai_screen_receipt(out,t,if(g=='eivGP')'evaluate'else'fit'),identity=identities[[t$plan]])}
    }
  }
  if(group=='report'||p$reuse_only){
    if(!all(vapply(required,function(r)oai_screen_verified(r$path,r$identity),logical(1))))stop('Saved-only prerequisites incomplete; no fitting or continuation started')
    return(oai_screen_collect(m,api,group))
  }
  api$oai_plan_require();job<-oai_path(root,'logs',paste0(group,'_',format(Sys.time(),'%Y%m%d_%H%M%S')))
  old<-options(eivGP.core_events=oai_path(job,'core-events'));on.exit(options(old),add=TRUE)
  close<-api$mixedgp_budget_enter(p$cores);on.exit(close(),add=TRUE,after=FALSE);cores<-api$mixedgp_budget_cores()
  writeLines(paste('Effective shared worker budget',cores,'; public resident-fit scheduling.'),oai_path(root,paste0(group,'_budget.txt')))
  stages<-if(group=='eivGP')unlist(m$stages)else 0L
  for(iter in stages){
    stage<-if(group=='eivGP')paste0('S',iter)else'BASELINE';tasks<-oai_screen_tasks(m,group,stage)
    work<-function(t,n,evaluate=FALSE){
      out<-cases[[t$plan]]$output;inp<-inputs[[t$plan]];id<-identities[[t$plan]];s<-oai_screen_settings(m,if(iter)iter else max(unlist(m$stages)),n)
      phase<-if(evaluate)'evaluate'else'fit';receipt<-oai_screen_receipt(out,t,phase)
      if(oai_screen_verified(receipt,id))return(invisible('REUSED'))
      if(t$kind=='eiv'){
        path<-api$oai_plan_fit_path(out,t)
        if(evaluate){
          if(!file.exists(path))stop('Stage checkpoint unavailable for evaluation')
          actual<-api$oai_plan_read_fit(path)$checkpoint$iteration
          if(is.null(s$ess_min) && actual!=iter || !is.null(s$ess_min) && actual<iter)stop('Stage checkpoint unavailable for evaluation')
          s$iter<-actual
          api$analyze_cell(t,inp,out,s,stage);folder<-oai_path(out,stage,'analysis',sprintf('fold%02d',t$fold),t$cell)
        }else{
          api$oai_plan_fit_task(t,inp,out,s,extend=FALSE);folder<-NULL
        }
      }else{
        methods<-if(t$kind=='baseline_group')c('UC-GP','LVGP','EzGP')else t$method
        failures<-character()
        for(method in methods){one<-t;one$method<-method;tryCatch(api$baseline_cell(one,inp,out,s,cc=t$kind=='cc'),error=function(e)failures<<-c(failures,paste(method,conditionMessage(e))))}
        if(length(failures))stop(paste(failures,collapse='; '))
        folder<-if(t$kind=='cc')oai_path(out,'CC',t$cell,'baselines',sprintf('fold%02d',t$fold),t$method)else
          oai_path(out,'baselines',sprintf('fold%02d',t$fold),t$method)
      }
      files<-if(is.null(folder))character()else sort(list.files(folder,recursive=TRUE,full.names=TRUE))
      api$atomic_rds(list(identity=id,task=t,files=files,md5=unname(tools::md5sum(files)),
        execution_source=base_id$application,finished=as.character(Sys.time())),receipt)
    }
    nchains<-if(group=='competitors')1L else if(is.null(m$chains))4L else as.integer(m$chains)
    status<-oai_screen_phase(tasks,'fit',cores,api,function(t,n)work(t,n),oai_path(job,stage),nchains,root,group)
    if(group=='eivGP'&&any(status$state=='FAILED'))stop('EIV stage incomplete; checkpoints retained. Inspect ',job)
    if(group=='eivGP'){
      status<-oai_screen_phase(tasks,'evaluate',cores,api,function(t,n)work(t,n,TRUE),oai_path(job,stage),nchains,root,group)
      if(any(status$state=='FAILED'))stop('EIV evaluation incomplete; continuation not started. Inspect ',job)
    }
    oai_screen_collect(m,api,group)
    if(any(status$state=='FAILED'))stop('OAI comparator tasks incomplete; completed results retained. Inspect ',job)
  }
  oai_publish_progress(root,group,'FITTING_COMPLETE',status)
  cat('[OAI ',group,'] FITTING AND EVALUATION COMPLETE; rendering its saved-results report.\n',sep='');flush.console()
  tables<-oai_screen_collect(m,api,group);attr(tables,'phase_status')<-status
  invisible(tables)
}

# OAI application: functionals. Reuses package model and sampler functions.

oai_trace_diagnostic <- function(values,chain) {
  groups<-split(values,chain);n<-min(lengths(groups));x<-do.call(cbind,lapply(groups,head,n))
  data.frame(states_per_chain=n,chains=ncol(x),mean=mean(values),posterior_SD=sd(values),
    Rhat=posterior::rhat(x),bulk_ESS=posterior::ess_bulk(x),tail_ESS=posterior::ess_tail(x),MCSE=posterior::mcse_mean(x))
}

oai_baseline_mean <- function(api,dir,X,C,mvec,method) {
  a<-readRDS(oai_path(dir,'fit.rds'));x<-sweep(sweep(X,2,a$center,'-'),2,a$scale,'/')
  if(method=='UC-GP') {
    nd<-getFromNamespace('make_mixedgp_competitor_data','eivGP')(x,C,mvec)
    mean<-as.numeric(stats::predict(a$model$fit,newdata=nd,type='UK',seCompute=FALSE,lightReturn=TRUE)$mean)
  }else if(method=='LVGP')mean<-as.numeric(LVGP::LVGP_predict(cbind(x,C),a$model$fit,MSE_on=0)$Y_hat)
  else mean<-api$oai_ezgp_predict_covariance_consistent(api$oai_ezgp_mixed_matrix(x,C,mvec),a$model)$Y_hat
  a$y_center+a$y_scale*mean
}

oai_saved_functionals <- function(a) {
  library(ggplot2)
  dest<-oai_path(a$report_dir,'report-assets','functionals');dir.create(dest,recursive=TRUE,showWarnings=FALSE)
  traces<-raw_traces<-curves<-surfaces<-precision<-integration<-list()
  code<-oai_source_files(a$repo)
  summarize_draws<-function(x)data.frame(mean=colMeans(x),lo=apply(x,2,quantile,.025),hi=apply(x,2,quantile,.975))
  for(key in names(a$paths)) {
    path<-a$paths[[key]];f<-path$fold;cell<-path$cell;z<-a$input$cells[[key]];p<-z$protocol
    cache_root<-oai_path(dest,'cache',digest::digest(unname(tools::md5sum(code))))
    dir.create(cache_root,recursive=TRUE,showWarnings=FALSE)
    cache<-oai_path(cache_root,paste0('fold',f,'_',cell,'.rds'))
    identity<-list(fit=unname(tools::md5sum(path$fit)),inputs=unname(tools::md5sum(a$input$files)),
      code=unname(tools::md5sum(code)),states=a$settings$functional_draws_per_chain,n_latent=a$settings$n_latent,
      analysis_modules=a$api$eivgp_module_identity())
    if(file.exists(cache)) {
      saved<-readRDS(cache)
      if(!identical(saved$identity,identity))stop('Derived functional cache differs; choose a fresh report directory: ',cache)
      result<-saved$result
    }else {
      cat('Read-only functionals fold',f,cell,'\n');flush.console()
      fit<-a$api$oai_plan_read_fit(path$fit)
      stopifnot(if(is.null(fit$scientific_ess$settings$ess_min))fit$checkpoint$iteration==a$numerics$iterations else fit$checkpoint$iteration>=a$numerics$iterations,fit$checkpoint$burn==a$numerics$burn)
      info<-fit$mcmc$mcmc_draw_info
      ids<-unlist(lapply(split(seq_len(nrow(info)),info$chain),function(v)
        v[unique(round(seq(1,length(v),length.out=min(length(v),a$settings$functional_draws_per_chain))))]),use.names=FALSE)
      chains<-info$chain[ids]
      xref<-vapply(z$x,function(n)median(z$train[[n]]),0)
      patterns<-apply(z$train[,z$c,drop=FALSE],1,paste,collapse='/');mode<-sort(names(table(patterns))[table(patterns)==max(table(patterns))])[1]
      cref<-as.numeric(strsplit(mode,'/',fixed=TRUE)[[1]])
      age<-seq(quantile(z$train[[z$x[1]]],.1),quantile(z$train[[z$x[1]]],.9),length.out=7)
      grid<-expand.grid(x=age,C1=seq_len(p$m_vec[1]))
      xx<-matrix(rep(xref,each=nrow(grid)),ncol=length(xref));xx[,1]<-grid$x
      cc<-matrix(rep(cref,each=nrow(grid)),ncol=length(cref));cc[,1]<-grid$C1;colnames(cc)<-p$C_order
      seed<-as.integer((as.double(p$seeds$posterior)+63001)%%(2^30-1)+1)
      md<-eivGP::predict_eivgp(fit,new_X=xx,new_C=cc,target='mean',draw_ids=ids,n_latent=a$settings$n_latent,
        include_process_uncertainty=TRUE,return_components=TRUE,seed=seed)
      cm<-attr(md,'conditional_means',exact=TRUE);stopifnot(is.matrix(cm))
      cv<-a$convergence[a$convergence$fold==f&a$convergence$cell==cell,]
      method_label<-if(nrow(cv)&&all(cv$core_pass))'EIV-GP'else'EIV-GP*'
      curve<-cbind(data.frame(fold=f,cell=cell,method=method_label,grid),summarize_draws(md))
      target_grade<-min(2L,p$m_vec[1]);target<-which(grid$C1==target_grade)[4L]
      mt<-data.frame(fold=f,cell=cell,target='m(X reference, C profile)',chain=chains,
        iteration=fit$checkpoint$burn+info$draw_within_chain[ids],value=cm[,target])
      pr<-cbind(data.frame(fold=f,cell=cell,target='m(X reference, C profile)',units='OAI_MEM'),oai_trace_diagnostic(cm[,target],chains))
      # Repeat only the fixed target with twice as many integration nodes.
      md2<-eivGP::predict_eivgp(fit,new_X=xx[target,,drop=FALSE],new_C=cc[target,,drop=FALSE],target='mean',
        draw_ids=ids,n_latent=2L*a$settings$n_latent,include_process_uncertainty=FALSE,seed=seed)
      check<-data.frame(fold=f,cell=cell,n_latent=a$settings$n_latent,n_latent_check=2L*a$settings$n_latent,
        mean=mean(cm[,target]),mean_check=mean(md2),difference=mean(md2)-mean(cm[,target]))
      refcell<-a$input$cells[[paste(f,'O50',sep='/')]]
      hidden<-which(refcell$train$calibration==0L&refcell$train$C1==target_grade)
      if(!length(hidden))hidden<-which(refcell$train$calibration==0L)
      for(j in seq_along(z$u)) {
        v<-rowMeans(matrix(fit$mcmc$samples_U[,hidden,j],nrow=nrow(info)))
        anchored<-isTRUE(fit$data$latent_scale_anchored)
        if(anchored)v<-fit$data$U_center[j]+fit$data$U_scale[j]*v
        name<-paste0('Mean hidden ',p$U_order[j])
        mt<-rbind(mt,data.frame(fold=f,cell=cell,target=name,chain=info$chain,
          iteration=fit$checkpoint$burn+info$draw_within_chain,value=v))
        pr<-rbind(pr,cbind(data.frame(fold=f,cell=cell,target=name,
          units=if(anchored)'physical U'else'unanchored model U'),oai_trace_diagnostic(v,info$chain)))
      }
      surface<-data.frame()
      if(cell!='O0') {
        cal<-z$train[z$train$calibration==1,,drop=FALSE]
        ug<-lapply(z$u,function(n)seq(quantile(cal[[n]],.1),quantile(cal[[n]],.9),length.out=11))
        if(length(z$u)>=2L)uv<-expand.grid(U1=ug[[1]],U2=ug[[2]])else uv<-data.frame(U1=ug[[1]],U2=0)
        uu<-matrix(rep(vapply(z$u,function(n)median(cal[[n]]),0),each=nrow(uv)),ncol=length(z$u))
        uu[,1]<-uv$U1;if(length(z$u)>=2L)uu[,2]<-uv$U2;colnames(uu)<-p$U_order
        x<-matrix(rep(xref,each=nrow(uv)),ncol=length(xref))
        fd<-eivGP::predict_eivgp(fit,new_X=x,new_U=uu,new_U_scale='raw',target='surface',draw_ids=ids,
          include_process_uncertainty=TRUE,seed=seed+3000L)
        surface<-cbind(data.frame(fold=f,cell=cell,uv),summarize_draws(fd))
      }
      raw<-oai_bind(lapply(c('V','sigma2','r'),function(n){
        v<-fit$mcmc[[paste0('samples_',n)]]
        if(is.null(v))return(data.frame())
        data.frame(fold=f,cell=cell,target=n,chain=info$chain,iteration=fit$checkpoint$burn+info$draw_within_chain,value=as.numeric(v))
      }))
      result<-list(curves=curve,surfaces=surface,traces=mt,raw_traces=raw,precision=pr,integration=check,
        reference=list(X=xref,C=cref,target_C1=target_grade,common_hidden_N=length(hidden)))
      saveRDS(list(identity=identity,result=result),cache)
      rm(fit,md,md2);gc(FALSE)
    }
    curves[[key]]<-result$curves;surfaces[[key]]<-result$surfaces;traces[[key]]<-result$traces
    raw_traces[[key]]<-result$raw_traces;precision[[key]]<-result$precision;integration[[key]]<-result$integration
  }
  # Native observable-mean curves on the same fixed training-defined slices.
  native<-list()
  for(f in a$settings$folds) {
    z<-a$input$cells[[paste(f,'O50',sep='/')]];p<-z$protocol
    grid<-curves[[paste(f,a$settings$cells[1],sep='/')]][,c('x','C1')]
    xref<-vapply(z$x,function(n)median(z$train[[n]]),0)
    patterns<-apply(z$train[,z$c,drop=FALSE],1,paste,collapse='/');mode<-sort(names(table(patterns))[table(patterns)==max(table(patterns))])[1]
    cref<-as.numeric(strsplit(mode,'/',fixed=TRUE)[[1]])
    xx<-matrix(rep(xref,each=nrow(grid)),ncol=length(xref));xx[,1]<-grid$x
    cc<-matrix(rep(cref,each=nrow(grid)),ncol=length(cref));cc[,1]<-grid$C1
    for(method in c('UC-GP','LVGP','EzGP')) {
      status<-a$method_status[a$method_status$fold==f&a$method_status$cell=='OI'&a$method_status$method==method,]
      if(!nrow(status)||status$state=='MISSING')next
      value<-oai_baseline_mean(a$api,status$source_directory,xx,cc,p$m_vec,method)
      native[[length(native)+1L]]<-data.frame(fold=f,cell='native OI',
        method=paste0(method,if(status$state=='OPTIMIZER_FLAGGED')'*'else''),grid,mean=value,lo=NA_real_,hi=NA_real_)
    }
  }
  out<-list(curves=oai_bind(c(curves,native)),surfaces=oai_bind(surfaces),traces=oai_bind(traces),
    raw_traces=oai_bind(raw_traces),precision=oai_bind(precision),integration=oai_bind(integration))
  for(n in names(out))if(nrow(out[[n]]))write.csv(out[[n]],oai_path(dest,paste0(n,'.csv')),row.names=FALSE)
  reference_fold<-a$settings$representative_fold
  plot_curve<-function(x)ggplot(x,aes(x,mean,colour=interaction(method,cell),group=interaction(method,cell)))+
    geom_ribbon(aes(ymin=lo,ymax=hi,fill=interaction(method,cell)),alpha=.06,colour=NA,na.rm=TRUE)+
    geom_line(linewidth=.65)+facet_wrap(~C1,ncol=2)+theme_bw()+
    labs(x=a$protocol$X_order[1],y='m(x,c): OAI_MEM',colour='Method / calibration',fill='Method / calibration')
  out$mean_plot<-plot_curve(out$curves[out$curves$fold==reference_fold,])
  out$mean_all_plot<-plot_curve(out$curves)+facet_grid(fold~C1)
  if(nrow(out$surfaces)) {
    out$surface_plot<-ggplot(out$surfaces[out$surfaces$fold==reference_fold,],aes(U1,U2,fill=mean))+
      geom_raster()+geom_contour(aes(z=mean),colour='white',alpha=.6)+facet_wrap(~cell)+
      scale_fill_viridis_c()+theme_bw()+labs(x=a$protocol$U_order[1],y=a$protocol$U_order[min(2,a$protocol$d)],fill='f(x,u)')
    out$surface_all_plot<-ggplot(out$surfaces,aes(U1,U2,fill=mean))+geom_raster()+
      facet_grid(fold~cell)+scale_fill_viridis_c()+theme_bw()+labs(x=a$protocol$U_order[1],y=a$protocol$U_order[min(2,a$protocol$d)],fill='f(x,u)')
  }
  traceplot<-function(d)ggplot(d[d$fold==reference_fold,],aes(iteration,value,colour=factor(chain)))+
    geom_line(alpha=.6,linewidth=.25)+facet_grid(target~cell,scales='free_y')+theme_bw()+
    labs(x='Saved iteration',y=NULL,colour='Chain')
  out$trace_plot<-traceplot(out$traces);out$raw_trace_plot<-traceplot(out$raw_traces)
  out
}

# OAI application: render. Reuses package model and sampler functions.

oai_render_plot <- function(args=commandArgs(TRUE)) {
  stopifnot(length(args)==5L)
  plot<-readRDS(args[1])
  ggplot2::ggsave(args[2],plot,width=as.numeric(args[3]),height=as.numeric(args[4]),
    dpi=as.numeric(args[5]),bg='white')
  stopifnot(file.exists(args[2]),file.info(args[2])$size>1000)
  invisible(args[2])
}

# OAI application: reporting. Reuses package model and sampler functions.

oai_table_html <- function(d,caption='',note='') {
  if(!is.data.frame(d)||!nrow(d))return('<p>No saved result available for this quantity.</p>')
  for(j in seq_along(d))if(is.numeric(d[[j]]))d[[j]]<-round(d[[j]],4)
  # Reuse Ocean's padding, rules, header shading, escaping and bold-mean markup.
  mixedgp_table_html(d,caption,note)
}

oai_bind <- function(parts) {
  parts<-Filter(function(x)is.data.frame(x)&&nrow(x)>0L,parts)
  if(!length(parts))return(data.frame())
  cols<-unique(unlist(lapply(parts,names)))
  do.call(rbind,lapply(parts,function(x){for(n in setdiff(cols,names(x)))x[[n]]<-NA;x[,cols,drop=FALSE]}))
}

oai_csv <- function(path,required=TRUE) {
  if(!file.exists(path)){if(required)stop('Missing saved report prerequisite: ',path);return(data.frame())}
  read.csv(path,stringsAsFactors=FALSE,check.names=FALSE)
}

oai_summary <- function(d,groups,metrics) {
  if(!nrow(d))return(data.frame())
  split_key<-interaction(d[,groups,drop=FALSE],drop=TRUE,lex.order=TRUE)
  oai_bind(lapply(split(d,split_key),function(x){
    out<-x[1,groups,drop=FALSE];out$folds<-length(unique(x$fold))
    for(k in metrics){v<-as.numeric(x[[k]]);out[[paste0(k,'_mean')]]<-mean(v);out[[paste0(k,'_SD')]]<-if(length(v)>1)sd(v)else NA_real_}
    out
  }))
}

oai_saved_source <- function(output_dir,source=NULL) {
  if(!is.null(source)&&nzchar(source)) {
    m<-jsonlite::fromJSON(source,simplifyVector=FALSE)
    if(!m$schema%in%c('oai-saved-report-v1','adni-saved-report-v1'))stop('Unsupported saved-report source')
    return(m)
  }
  p<-oai_path(output_dir,'application_manifest.json')
  if(!file.exists(p))stop('Saved report requires application_manifest.json or --report-source=FILE')
  m<-jsonlite::fromJSON(p,simplifyVector=FALSE)
  list(schema='oai-saved-report-v1',output=output_dir,input=m$cases[[1]]$input,
    stage=paste0('S',tail(unlist(m$stages),1)),fit_package_version=m$package_version,fit_profile=m$profile,
    fitting_manifest=p,input_sha256=m$frozen_inputs)
}

oai_report_prepare <- function(repo,data_dir,output_dir,settings,report_source=NULL,report_dir=output_dir,group='report') {
  for(p in c('jsonlite','digest','ggplot2','patchwork','posterior'))if(!requireNamespace(p,quietly=TRUE))stop('Report dependency: ',p)
  api<-oai_load_api(repo);api$application$oai_output_check(report_dir,repo)
  m<-oai_saved_source(output_dir,report_source)
  fitting<-if(!is.null(m$fitting_manifest))jsonlite::fromJSON(m$fitting_manifest,simplifyVector=FALSE)else list()
  for(n in setdiff(c('folds','cells'),attr(settings,'explicit_selection')))if(length(fitting[[n]]))
    settings[[n]]<-if(n=='folds')as.integer(unlist(fitting[[n]]))else unlist(fitting[[n]],use.names=FALSE)
  for(n in c('draws_per_chain','baseline_profile'))if(!is.null(fitting[[n]]))settings[[n]]<-fitting[[n]]
  if(length(fitting$samplers))for(cell in names(fitting$samplers))
    settings[[paste0('sampler_',tolower(cell))]]<-fitting$samplers[[cell]]
  if(!is.null(m$fit_profile))settings$profile<-m$fit_profile
  if(!is.null(m$input))data_dir<-m$input
  if(!is.null(m$output))output_dir<-m$output
  input<-api$oai_plan_inputs(data_dir)
  if(length(m$input_sha256))for(p in names(m$input_sha256))
    if(digest::digest(file=p,algo='sha256')!=m$input_sha256[[p]])stop('Frozen report input changed: ',p)
  # Historical fitting identities are read, not rewritten to match reporting code.
  identity<-oai_path(output_dir,'workflow_identity.rds')
  if(file.exists(identity)) {
    id<-readRDS(identity)
    if(!is.null(id$inputs)&&!identical(unname(tools::md5sum(input$files)),id$inputs))stop('Saved fit input identity differs')
  }
  if(is.null(m$stage)) {
    stages<-list.dirs(output_dir,recursive=FALSE,full.names=FALSE)
    stages<-stages[grepl('^S[0-9]+K?$',stages)]
    amount<-as.numeric(gsub('K','',sub('^S','',stages)))*ifelse(grepl('K$',stages),1000,1)
    if(!length(stages))stop('No saved EIV stages')
    m$stage<-stages[which.max(amount)]
  }
  pred<-conv<-up<-um<-pp<-pm<-parameters<-statuses<-paths<-list()
  final<-m$stage
  selected<-expand.grid(fold=settings$folds,cell=settings$cells,stringsAsFactors=FALSE)
  want_eiv<-group!='competitors';want_baseline<-group!='eivGP'
  for(i in seq_len(nrow(selected))) {
    f<-selected$fold[i];cell<-selected$cell[i];z<-input$cells[[paste(f,cell,sep='/')]]
    dir<-oai_path(output_dir,final,'analysis',sprintf('fold%02d',f),cell)
    if(want_eiv) {
      d<-oai_csv(oai_path(dir,'prediction_metrics.csv'));d<-d[d$protocol%in%c('OI','CC-known'),,drop=FALSE]
      pred[[length(pred)+1L]]<-d
      cv<-oai_csv(oai_path(dir,'convergence.csv'));conv[[length(conv)+1L]]<-cv
      par<-oai_csv(oai_path(dir,'diagnostic_parameters.csv'));par$fold<-f;par$cell<-cell
      parameters[[length(parameters)+1L]]<-par
      paths[[paste(f,cell,sep='/')]]<-list(fold=f,cell=cell,dir=dir,
        fit=oai_path(output_dir,'fits',sprintf('fold%02d',f),cell,'fit.rds'))
      if(!file.exists(paths[[paste(f,cell,sep='/')]]$fit))stop('Missing saved EIV fit; report never fits')
      if(cell!='O0') {
        u<-oai_csv(oai_path(dir,'U_predictions.csv'));u$fold<-f;u$cell<-cell
        up[[length(up)+1L]]<-u;um[[length(um)+1L]]<-oai_csv(oai_path(dir,'U_metrics.csv'))
        p<-oai_csv(oai_path(dir,'PROSPECTIVE_U_predictions.csv'));p$fold<-f;p$cell<-cell
        pp[[length(pp)+1L]]<-p;pm[[length(pm)+1L]]<-oai_csv(oai_path(dir,'PROSPECTIVE_U_metrics.csv'))
      }
    }
  }
  if(want_baseline)for(f in settings$folds)for(cell in c('OI',intersect(settings$cells,c('O20','O50'))))for(method in c('UC-GP','LVGP','EzGP')) {
    dir<-if(cell=='OI')oai_path(output_dir,'baselines',sprintf('fold%02d',f),method)else
      oai_path(output_dir,'CC',cell,'baselines',sprintf('fold%02d',f),method)
    version<-m$fit_package_version
    if(length(m$baseline_overrides))for(b in m$baseline_overrides)
      if(b$fold==f&&b$cell==cell&&b$method==method){dir<-b$directory;version<-b$package_version}
    valid<-all(file.exists(oai_path(dir,c('metrics.csv','fit.rds','prediction_draws.rds'))))
    state<-'MISSING';optimization<-'not_available'
    if(valid) {
      d<-oai_csv(oai_path(dir,'metrics.csv'));d$method<-paste0(method,if(cell=='OI')''else'-CC');d$cell<-cell
      d$stage<-'BASELINE';pred[[length(pred)+1L]]<-d
      a<-readRDS(oai_path(dir,'fit.rds'));optimization<-a$model$optimization_status
      if(is.null(optimization))optimization<-'NOT_EXPOSED'
      state<-if(grepl('nonconverged',optimization))'OPTIMIZER_FLAGGED'else'AVAILABLE'
      rm(a)
    }
    statuses[[length(statuses)+1L]]<-data.frame(fold=f,cell=cell,method=paste0(method,if(cell=='OI')''else'-CC'),
      state=state,optimization_status=optimization,fit_package_version=if(is.null(version))'recorded in saved provenance'else version,
      source_directory=dir)
  }
  p<-input$cells[[1]]$protocol
  z<-input$cells[['1/O50']];cohort<-rbind(z$train,z$test);cohort<-cohort[order(cohort$RID),]
  metadata_path<-oai_path(data_dir,'DATASET.json')
  meta<-if(file.exists(metadata_path))jsonlite::fromJSON(metadata_path,simplifyVector=FALSE)else list()
  meta<-modifyList(meta,if(is.null(m$metadata))list()else m$metadata)
  if(is.null(meta))meta<-list()
  if(is.null(meta$title))meta$title<-'OAI cartilage morphology and fixed-location joint-space width'
  if(is.null(meta$population))meta$population<-if(is.null(p$population))'Population specified in private preparation records'else p$population
  if(is.null(meta$calibration_design))meta$calibration_design<-if(is.null(p$calibration_design))
    'C1-proportional stratified sampling, fixed seeded within-stratum order; O20 nested in O50'else
    paste(p$calibration_design,'; strata:',paste(unlist(p$stratification_columns),collapse=', '))
  analysis<-list(repo=repo,api=api,source=m,input=input,output_dir=output_dir,report_dir=report_dir,
    settings=settings,group=group,stage=final,metadata=meta,protocol=p,cohort=cohort,
    prediction=oai_bind(pred),convergence=oai_bind(conv),parameters=oai_bind(parameters),
    U_predictions=oai_bind(up),U_metrics=oai_bind(um),prospective_predictions=oai_bind(pp),
    prospective_metrics=oai_bind(pm),method_status=oai_bind(statuses),paths=paths)
  if(nrow(analysis$convergence)) {
    nums<-unique(analysis$convergence[,c('iterations','burn','chains')])
    nums<-mixedgp_application_numerics(nums)
    analysis$numerics<-nums
    analysis$settings$chains<-as.integer(nums$chains)
    analysis$settings$iter<-as.integer(nums$iterations)
    analysis$settings$burn<-as.integer(nums$burn)
  }else analysis$numerics<-data.frame(iterations=NA,burn=NA,chains=NA)
  if(!is.null(m$fitting_manifest)) {
    stages<-paste0('S',unlist(fitting$stages))
    analysis$settings$stages<-as.integer(unlist(fitting$stages))
  }else {
    stages<-paste0('S',unlist(m$evaluation_iterations))
    if(length(m$evaluation_iterations))analysis$settings$stages<-as.integer(unlist(m$evaluation_iterations))
  }
  # Historical folders may use explicit stage labels rather than S<iteration>.
  if(length(m$stage_directories))stages<-unlist(m$stage_directories)
  stability<-list()
  for(st in stages)for(i in seq_len(nrow(selected))) {
    path<-oai_path(output_dir,st,'analysis',sprintf('fold%02d',selected$fold[i]),selected$cell[i],'prediction_metrics.csv')
    d<-oai_csv(path,FALSE)
    if(nrow(d)){d<-d[d$protocol=='OI',];d$stage<-st;stability[[length(stability)+1L]]<-d}
  }
  analysis$stage_scores<-oai_bind(stability)
  oai_report_exports(analysis)
  analysis$plots<-oai_report_plots(analysis)
  if(want_eiv) {
    analysis$functionals<-oai_saved_functionals(analysis)
  }
  analysis$artifact_manifest<-oai_paper_assets(analysis)
  analysis
}

oai_report_exports <- function(a) {
  dest<-oai_path(a$report_dir,'report-assets');dir.create(dest,recursive=TRUE,showWarnings=FALSE)
  for(n in c('prediction','convergence','parameters','U_metrics','prospective_metrics','method_status','stage_scores'))
    if(nrow(a[[n]]))write.csv(a[[n]],oai_path(dest,paste0(n,'.csv')),row.names=FALSE)
  a$api$atomic_json(list(author='Xinye Chen',saved_only=TRUE,new_MCMC=FALSE,new_optimization=FALSE,
    fitting_source=a$source,stage=a$stage,saved_fitting_numerics=a$numerics,report_settings=a$settings,
    input_md5=as.list(tools::md5sum(a$input$files)),analysis_modules=a$api$eivgp_module_identity(),
    reporting_source=as.list(tools::md5sum(c(oai_path(a$repo,'paper/applications/oai.Rmd'),oai_source_files(a$repo))))),
    oai_path(dest,'REPORT_PROVENANCE.json'))
}

oai_metric_table <- function(a,protocol='OI',full=FALSE) {
  d<-a$prediction[a$prediction$protocol==protocol,,drop=FALSE]
  metrics<-c('RMSE','RMSE_SD','MAE_SD','CRPS_SD','Coverage95','Width95_SD','IntervalScore95_SD')
  tab<-oai_summary(d,c('cell','method'),metrics)
  if(!nrow(tab))return('<p>No available saved scores for this protocol.</p>')
  tab$CORE<-vapply(seq_len(nrow(tab)),function(i){
    if(tab$method[i]!='EIV-GP')return('N/A')
    q<-a$convergence[a$convergence$cell==tab$cell[i],];paste0(sum(q$core_pass),'/',nrow(q))},'')
  method<-tab$method
  for(i in seq_len(nrow(tab))) {
    flagged<-if(method[i]=='EIV-GP')tab$CORE[i]!=paste0(tab$folds[i],'/',tab$folds[i])else
      any(a$method_status$cell==tab$cell[i]&a$method_status$method==method[i]&a$method_status$state=='OPTIMIZER_FLAGGED')
    if(flagged)method[i]<-paste0(method[i],'*')
  }
  result<-data.frame(Method=method,Calibration=tab$cell,
    Test=if(protocol=='OI')'(X,C)'else ifelse(tab$method=='EIV-GP','(X,U)','(X,U,C)'),
    Folds=tab$folds,check.names=FALSE)
  shown<-if(full)metrics else c('RMSE','RMSE_SD','CRPS_SD','Coverage95','Width95_SD')
  labels<-c(RMSE='RMSE',RMSE_SD='RMSE/sY',MAE_SD='MAE/sY',CRPS_SD='CRPS/sY',
    Coverage95='Coverage (%)',Width95_SD='Width/sY',IntervalScore95_SD='Interval score/sY')
  best_mask<-matrix(FALSE,nrow(result),length(shown),dimnames=list(NULL,unname(labels[shown])))
  for(k in shown) {
    scale<-if(k=='Coverage95')100 else 1
    values<-sprintf(if(k=='Coverage95')'%.2f (%.2f)'else'%.4f (%.4f)',
      scale*tab[[paste0(k,'_mean')]],scale*tab[[paste0(k,'_SD')]])
    if(k%in%c('RMSE','RMSE_SD','MAE_SD','CRPS_SD','IntervalScore95_SD')) {
      strata<-if(protocol=='OI')rep('shared',nrow(tab))else tab$cell
      for(g in unique(strata)) {
        ix<-which(strata==g);best<-min(tab[[paste0(k,'_mean')]][ix])
        hit<-ix[tab[[paste0(k,'_mean')]][ix]==best]
        # Unconverged methods remain eligible; the asterisk records their status.
        best_mask[hit,labels[[k]]]<-TRUE
      }
    }
    result[[labels[[k]]]]<-values
  }
  result$CORE<-tab$CORE
  attr(result,'best')<-best_mask
  oai_table_html(result,if(protocol=='OI')'Ordinal-input prediction: fold mean (sample SD)'else
    'Complete-case prediction: fold mean (sample SD)')
}

oai_u_table <- function(d) {
  if(!nrow(d))return(data.frame())
  tab<-oai_summary(d,c('cell','coordinate'),c('RMSE','RMSE_SD','MAE_SD','CRPS_SD','Coverage95','Width95_SD','IntervalScore95_SD'))
  tab$gain_vs_C1_percent<-vapply(seq_len(nrow(tab)),function(i){
    v<-d[d$cell==tab$cell[i]&d$coordinate==tab$coordinate[i],]
    100*(1-sqrt(mean(v$RMSE^2))/sqrt(mean(v$C1_mean_RMSE^2)))},0)
  out<-tab[,c('cell','coordinate','folds')]
  for(k in c('RMSE','RMSE_SD','MAE_SD','CRPS_SD','Coverage95','Width95_SD','IntervalScore95_SD'))
    out[[k]]<-sprintf('%.4f (%.4f)',tab[[paste0(k,'_mean')]],tab[[paste0(k,'_SD')]])
  out$gain_vs_C1_percent<-sprintf('%+.2f',tab$gain_vs_C1_percent);out
}

oai_target_table <- function(a,all_folds=FALSE) {
  if(is.null(a$functionals))return(data.frame())
  d<-a$functionals$precision
  if(!all_folds)d<-d[d$fold==a$settings$representative_fold,,drop=FALSE]
  z<-a$input$cells[[1]];target<-rep('m(x,c)',nrow(d))
  for(j in seq_along(z$u)) {
    name<-paste0('Mean hidden ',a$protocol$U_order[j]);at<-d$target==name
    label<-if(grepl('FDG',a$protocol$U_order[j]))'FDG'else if(grepl('HIP',a$protocol$U_order[j]))'HIP1000/ICV'else a$protocol$U_order[j]
    target[at]<-paste0(label,' mean (',ifelse(d$cell[at]=='O0','model','raw'),')')
  }
  d$Target<-target
  d<-d[order(match(d$target,unique(d$target)),match(d$cell,c('O0','O20','O50')),d$fold),,drop=FALSE]
  out<-data.frame(Target=d$Target,Calibration=d$cell,`Mean (SD)`=sprintf('%.4f (%.4f)',d$mean,d$posterior_SD),
    `R-hat`=sprintf('%.3f',d$Rhat),`Bulk ESS`=sprintf('%.1f',d$bulk_ESS),`Tail ESS`=sprintf('%.1f',d$tail_ESS),
    MCSE=ifelse(abs(d$MCSE)<1e-4,formatC(d$MCSE,format='e',digits=2),sprintf('%.4f',d$MCSE)),check.names=FALSE)
  if(all_folds)out<-data.frame(Fold=d$fold,out,check.names=FALSE)
  out
}

oai_integration_table <- function(a,all_folds=FALSE) {
  if(is.null(a$functionals))return(data.frame())
  d<-a$functionals$integration
  if(!all_folds)d<-d[d$fold==a$settings$representative_fold,,drop=FALSE]
  data.frame(Fold=d$fold,Calibration=d$cell,Nodes=d$n_latent,`Check nodes`=d$n_latent_check,
    `Mean at nodes`=sprintf('%.5f',d$mean),`Mean at check nodes`=sprintf('%.5f',d$mean_check),
    Difference=sprintf('%+.5f',d$difference),check.names=FALSE)
}

oai_report_plots <- function(a) {
  library(ggplot2);library(patchwork)
  z<-a$input$cells[['1/O50']];d<-a$cohort;p<-a$protocol;plots<-list()
  long<-oai_bind(lapply(seq_along(z$u),function(j)data.frame(U=d[[z$u[j]]],Y=d$Y,C1=factor(d$C1),coordinate=p$U_order[j])))
  plots$cohort<-ggplot(long,aes(U,Y,colour=C1))+geom_point(alpha=.55,size=1.3)+
    facet_wrap(~coordinate,scales='free_x')+labs(x='Observed U reference',y='Response Y',colour=p$C_order[1])+
    theme_minimal(base_size=11)+theme(legend.position='bottom')
  xlong<-oai_bind(lapply(seq_along(z$x),function(j)data.frame(X=d[[z$x[j]]],Y=d$Y,C1=factor(d$C1),coordinate=p$X_order[j])))
  plots$covariates<-ggplot(xlong,aes(X,Y,colour=C1))+geom_point(alpha=.55,size=1.3)+
    facet_wrap(~coordinate,scales='free_x')+labs(x='Observed covariate',y='Response Y',colour=p$C_order[1])+
    theme_minimal(base_size=11)+theme(legend.position='bottom')
  shown<-seq_len(min(2L,length(z$c)))
  cu<-oai_bind(lapply(shown,function(k)oai_bind(lapply(seq_along(z$u),function(j)
    data.frame(grade=factor(d[[z$c[k]]],levels=seq_len(p$m_vec[k])),U=d[[z$u[j]]],
      item=paste0('C',k,': ',p$C_order[k]),coordinate=p$U_order[j])))))
  plots$C_U<-ggplot(cu,aes(grade,U))+geom_boxplot(width=.55,outlier.shape=NA,fill='grey90',colour='grey35')+
    geom_point(aes(colour=grade),position=position_jitter(width=.12,height=0,seed=6203L),alpha=.45,size=.9)+
    facet_wrap(~item+coordinate,ncol=length(z$u),scales='free_y',labeller=label_wrap_gen(36))+
    labs(x='Ordered model code',y='Measured U reference (raw units)')+
    theme_minimal(base_size=11)+theme(legend.position='none')
  counts<-oai_bind(lapply(seq_along(z$c),function(j){t<-as.data.frame(table(factor(d[[z$c[j]]],levels=seq_len(p$m_vec[j]))));
    data.frame(item=p$C_order[j],grade=t$Var1,N=t$Freq)}))
  plots$C_counts<-ggplot(counts,aes(grade,N,fill=grade))+geom_col()+facet_wrap(~item,ncol=3,scales='free_y')+
    labs(x='Ordered model code',y='Participants')+theme_bw()+theme(legend.position='none')
  u<-a$U_predictions
  if(nrow(u)) {
    u$role<-ifelse(u$calibration==1,'Observed calibration','Masked training U')
    u$C1<-factor(u$C1)
    draw_u<-function(x,all=FALSE){
      ggplot(x,aes(true_std,mean_std,colour=C1))+geom_abline(slope=1,intercept=0,linetype=2,colour='grey50')+
        geom_linerange(data=x[x$calibration==0,],aes(ymin=lo_std,ymax=hi_std),alpha=.15)+
        geom_point(aes(shape=role),size=1,alpha=.6)+scale_shape_manual(values=c('Masked training U'=16,'Observed calibration'=1))+
        facet_grid(if(all)fold~coordinate else cell~coordinate,scales='free')+
        labs(x='Reference U (training SD units; evaluation only)',y='Posterior mean (training SD units)',colour=p$C_order[1],shape=NULL)+theme_bw()
    }
    plots$imputation<-draw_u(u[u$fold==a$settings$representative_fold,])
    for(cell in unique(u$cell))plots[[paste0('imputation_',cell)]]<-draw_u(u[u$cell==cell,],TRUE)
    hidden<-u[u$calibration==0,]
    density<-oai_bind(list(data.frame(hidden[,c('cell','coordinate')],value=hidden$true_std,distribution='Reference'),
      data.frame(hidden[,c('cell','coordinate')],value=hidden$mean_std,distribution='Posterior means')))
    plots$input_distribution<-ggplot(density,aes(value,colour=distribution,fill=distribution))+
      geom_density(alpha=.12)+facet_grid(cell~coordinate,scales='free')+
      labs(x='Training SD units',y='Density',colour=NULL,fill=NULL)+theme_bw()
    class_rows<-oai_bind(lapply(split(hidden,interaction(hidden$fold,hidden$cell,hidden$coordinate,hidden$C1,drop=TRUE)),function(x)
      data.frame(fold=x$fold[1],cell=x$cell[1],coordinate=x$coordinate[1],C1=as.character(x$C1[1]),N=nrow(x),
        RMSE_SD=sqrt(mean((x$mean_std-x$true_std)^2)),Coverage95=mean(x$true_U>=x$q025&x$true_U<=x$q975))))
    plots$class_table<-oai_summary(class_rows,c('cell','coordinate','C1'),c('RMSE_SD','Coverage95'))
    common<-hidden[as.character(hidden$common_hidden)%in%c('TRUE','True','1'),]
    common_rows<-oai_bind(lapply(split(common,interaction(common$fold,common$cell,common$coordinate,drop=TRUE)),function(x)
      data.frame(fold=x$fold[1],cell=x$cell[1],coordinate=x$coordinate[1],N=nrow(x),RMSE_SD=sqrt(mean((x$mean_std-x$true_std)^2)))))
    plots$common_hidden<-oai_summary(common_rows,c('cell','coordinate'),'RMSE_SD')
    plots$support<-oai_bind(lapply(split(hidden,interaction(hidden$cell,hidden$coordinate,drop=TRUE)),function(x)
      data.frame(cell=x$cell[1],coordinate=x$coordinate[1],fold_rows=nrow(x),lower_nonpositive=sum(x$q025<=0),
        mean_nonpositive=sum(x$posterior_mean<=0),interval_outside_C_grade='Not applicable: C is not an imaging bin')))
  }
  if(nrow(u))plots<-modifyList(plots,oai_input_distributions(a))
  plots
}

oai_input_distributions <- function(a) {
  library(ggplot2);parts<-list()
  for(key in names(a$paths)) {
    path<-a$paths[[key]];if(path$cell=='O0'||path$fold!=a$settings$representative_fold)next
    z<-a$input$cells[[key]]
    hidden<-which(a$input$cells[[paste(path$fold,'O50',sep='/')]]$train$calibration==0L)
    saved<-readRDS(oai_path(path$dir,'U_draws.rds'))
    if(!identical(saved$scale,'raw'))stop('Raw saved U draws required for physical distribution checks')
    ids<-unique(round(seq(1,dim(saved$draws)[1],length.out=min(128,dim(saved$draws)[1]))))
    for(j in seq_along(z$u)) {
      ref<-z$train[[z$u[j]]];mu<-mean(ref);ss<-sd(ref)
      draw<-matrix(saved$draws[ids,hidden,j],nrow=length(ids))
      tag<-data.frame(cell=path$cell,coordinate=z$protocol$U_order[j])
      parts[[length(parts)+1L]]<-data.frame(tag,C1=rep(as.character(z$train$C1[hidden]),each=length(ids)),
        value=as.vector((draw-mu)/ss),distribution='Posterior mixture')
      parts[[length(parts)+1L]]<-data.frame(tag,C1=as.character(z$train$C1[hidden]),value=(ref[hidden]-mu)/ss,distribution='Reference')
      means<-colMeans(matrix(saved$draws[,hidden,j],ncol=length(hidden)))
      parts[[length(parts)+1L]]<-data.frame(tag,C1=as.character(z$train$C1[hidden]),value=(means-mu)/ss,distribution='Posterior means')
    }
  }
  d<-oai_bind(parts);if(!nrow(d))return(list())
  d$C1<-factor(d$C1);reference<-d[d$distribution=='Reference',];post<-d[d$distribution!='Reference',]
  draw_distribution<-function(class=FALSE) {
    g<-ggplot()+geom_histogram(data=reference,aes(value,after_stat(density)),bins=15,
      fill='grey80',colour='white')+geom_rug(data=reference,aes(value),alpha=.3)+
      geom_density(data=post,aes(value,colour=distribution,linetype=distribution),linewidth=.7)+
      scale_linetype_manual(values=c('Posterior mixture'='solid','Posterior means'='dashed'))+
      theme_minimal(base_size=11)+theme(legend.position='bottom')+
      labs(x='Common-hidden U (training SD units)',y='Density',colour=NULL,linetype=NULL)
    if(class)g+facet_grid(cell+coordinate~C1,scales='free')else g+facet_grid(cell~coordinate,scales='free')
  }
  list(input_distribution=draw_distribution(),class_distribution=draw_distribution(TRUE))
}

oai_paper_assets <- function(a) {
  dest<-oai_path(a$report_dir,'report-assets','paper');dir.create(dest,recursive=TRUE,showWarnings=FALSE)
  plots<-list(fig_oai_cohort=a$plots$cohort,fig_oai_C_support=a$plots$C_counts,
    fig_oai_Y_X=a$plots$covariates,fig_oai_C_U=a$plots$C_U,
    fig_oai_target_traces=a$functionals$trace_plot,fig_oai_mean_functions=a$functionals$mean_plot,
    fig_oai_response_surfaces=a$functionals$surface_plot,fig_oai_U_imputation=a$plots$imputation,
    fig_oai_U_distributions=a$plots$input_distribution,fig_oai_C1_U_distributions=a$plots$class_distribution)
  manifest<-list()
  for(name in names(plots)) {
    file<-paste0(name,'.pdf');available<-!is.null(plots[[name]])
    if(available)ggplot2::ggsave(oai_path(dest,file),plots[[name]],width=12,height=if(grepl('support|traces|C1',name))9 else 6,
      device=grDevices::pdf,limitsize=FALSE)
    manifest[[length(manifest)+1L]]<-data.frame(label=name,file=file,status=if(available)'available'else'not_applicable',
      N=a$input$N,profile=a$settings$profile,MD5=if(available)unname(tools::md5sum(oai_path(dest,file)))else'')
  }
  oi<-a$prediction[a$prediction$protocol=='OI',,drop=FALSE]
  cc<-a$prediction[a$prediction$protocol=='CC-known',,drop=FALSE]
  metrics<-c('RMSE','RMSE_SD','MAE_SD','CRPS_SD','Coverage95','Width95_SD','IntervalScore95_SD')
  tables<-list(tab_oai_design=a$input$audit,
    tab_oai_C_definitions=oai_C_dictionary(a$protocol$C_order,a$protocol$m_vec,a$metadata$C_labels),
    tab_oai_prediction_OI=oai_summary(oi,c('cell','method'),metrics),
    tab_oai_prediction_CC=oai_summary(cc,c('cell','method'),metrics),tab_oai_training_U=oai_u_table(a$U_metrics),
    tab_oai_prospective_U=oai_u_table(a$prospective_metrics),tab_oai_convergence=a$convergence)
  for(name in names(tables)) {
    d<-tables[[name]];available<-nrow(d)>0L;file<-paste0(name,'.csv')
    if(available) {
      write.csv(d,oai_path(dest,file),row.names=FALSE)
      writeLines(oai_table_html(d,name),oai_path(dest,paste0(name,'.html')))
    }
    manifest[[length(manifest)+1L]]<-data.frame(label=name,file=file,status=if(available)'available'else'unavailable',
      N=a$input$N,profile=a$settings$profile,MD5=if(available)unname(tools::md5sum(oai_path(dest,file)))else'')
  }
  out<-oai_bind(manifest);write.csv(out,oai_path(dest,'artifact_manifest.csv'),row.names=FALSE);out
}

# OAI application: settings. Reuses package model and sampler functions.

oai_resolve_path <- function(path,repo) {
  path<-path.expand(path)
  if(!grepl('^(/|[A-Za-z]:[/\\\\])',path))path<-oai_path(repo,path)
  base<-path;suffix<-character()
  while(!file.exists(base)&&dirname(base)!=base) {
    suffix<-c(basename(base),suffix);base<-dirname(base)
  }
  resolved<-normalizePath(base,mustWork=TRUE)
  if(length(suffix))oai_path(resolved,paste(suffix,collapse='/'))else resolved
}

oai_cli_params <- function(args) {
  allowed<-c('cores','profile','dataset','data_root','data_dir','output','output_dir','config','manifest',
    'chains','iter','burn','stages','folds','cells','draws_per_chain','first_checkpoint','checkpoint_step',
    'baseline_profile','sampler_o0','sampler_o20','sampler_o50','representative_fold','appendix',
    'plan','reuse_only','report_source','report_dir','functional_draws_per_chain','n_latent','run_id','ess_min',
    'run_folds')
  p<-list()
  for(a in args) {
    if(!grepl('^--[^=]+=.+$',a))stop('Use --name=value: ',a)
    k<-gsub('-','_',sub('^--([^=]+)=.*$','\\1',a));v<-sub('^[^=]+=','',a)
    if(!k%in%allowed)stop('Unknown OAI option: ',k)
    if(k%in%names(p))stop('Repeated option: ',k)
    if(k%in%c('plan','reuse_only','appendix')) {
      if(!v%in%c('true','false'))stop(k,' must be true or false')
      v<-v=='true'
    }
    p[[k]]<-v
  }
  p
}

oai_integer <- function(value,name,minimum=1L) {
  n<-suppressWarnings(as.numeric(value))
  if(length(n)!=1L||!is.finite(n)||n<minimum||n!=floor(n)||n>.Machine$integer.max)
    stop(name,' must be an integer >= ',minimum)
  as.integer(n)
}

oai_integer_list <- function(value,name) {
  if(is.character(value)&&length(value)==1L) {
    if(grepl('^[0-9]+:[0-9]+$',value)) {
      z<-as.integer(strsplit(value,':',fixed=TRUE)[[1]]);value<-seq.int(z[1],z[2])
    }else value<-strsplit(value,',',fixed=TRUE)[[1]]
  }
  unname(vapply(unlist(value),oai_integer,1L,name=name))
}

oai_settings <- function(profile='pilot',cores=4L,...) {
  overrides<-list(...);profile<-match.arg(profile,c('pilot','paper','quick'))
  allocation<-mixedgp_application_allocation(oai_integer(cores,"cores"));cores<-allocation$requested_cores
  if(!is.null(overrides$chains) && oai_integer(overrides$chains,'chains')!=allocation$chains)
    stop('For a new OAI run, chains must equal min(--cores, 4); use a matching core budget.')
  defaults<-switch(profile,paper=list(chains=allocation$chains,iter=50000L,burn=5000L,stages=50000L),
    pilot=list(chains=allocation$chains,iter=1500L,burn=500L,stages=c(1000L,1500L)),
    quick=list(chains=allocation$chains,iter=80L,burn=30L,stages=80L))
  s<-modifyList(c(defaults,list(profile=profile,cores=cores,folds=if(profile=='quick')1L else 1:5,
    cells=c('O50','O20','O0'),draws_per_chain=if(profile=='quick')10L else 250L,
    checkpoint_step=1000L,baseline_profile=if(profile=='quick')'quick'else'paper',
    sampler_o0='HMC',sampler_o20='HMC',sampler_o50='ESS',representative_fold=1L,
    appendix=TRUE,functional_draws_per_chain=64L,n_latent=256L,run_id='local')),overrides)
  s$ess_min<-mixedgp_parse_ess_min(s$ess_min)
  if(!is.null(s$ess_min) && as.numeric(s$chains)<4L)stop('Scientific ESS requires at least four chains.')
  for(n in c('cores','chains','iter','draws_per_chain','checkpoint_step','representative_fold',
             'functional_draws_per_chain','n_latent'))s[[n]]<-oai_integer(s[[n]],n)
  s$burn<-oai_integer(s$burn,'burn',0L);s$folds<-oai_integer_list(s$folds,'folds')
  if(any(!s$folds%in%1:5)||anyDuplicated(s$folds))stop('folds must be distinct members of1:5')
  s$cells<-if(length(s$cells)==1L)strsplit(s$cells,',',fixed=TRUE)[[1]]else unlist(s$cells)
  if(!length(s$cells)||any(!s$cells%in%c('O0','O20','O50'))||anyDuplicated(s$cells))stop('Invalid cells')
  if(!is.null(overrides$stages)) {
    s$stages<-oai_integer_list(overrides$stages,'stages')
    if(!is.null(overrides$iter)&&tail(s$stages,1)!=s$iter)stop('iter must equal the last stage')
    s$iter<-tail(s$stages,1)
  }else if(!is.null(overrides$iter))s$stages<-s$iter
  if(any(diff(s$stages)<=0)||min(s$stages)<=s$burn)stop('Stages must increase and exceed burn')
  s$first_checkpoint<-oai_integer(if(is.null(s$first_checkpoint))min(s$stages[1],s$burn+250L)else s$first_checkpoint,'first_checkpoint')
  if(s$first_checkpoint<=s$burn||s$first_checkpoint>min(s$stages))stop('first-checkpoint must exceed burn and not exceed first stage')
  if(!s$baseline_profile%in%c('paper','quick'))stop('baseline-profile must be paper or quick')
  for(n in c('sampler_o0','sampler_o20','sampler_o50'))if(!s[[n]]%in%c('HMC','ESS'))stop('Unknown sampler')
  attr(s,'explicit_selection')<-intersect(names(overrides),c('folds','cells'))
  s
}

oai_method_run <- function(repo,params,group) {
  if(!is.null(params$config)) {
    config<-jsonlite::fromJSON(oai_resolve_path(params$config,repo),simplifyVector=TRUE)
    params<-modifyList(config,params)
  }
  get<-function(n,default)if(is.null(params[[n]]))default else params[[n]]
  if(any(c('fit_jobs','evaluation_jobs')%in%names(params)))
    stop('Resident-job options have been retired; specify only --cores for resource allocation.')
  if(!is.null(params$manifest)) {
    if(!is.null(params$ess_min)) {
      declared<-jsonlite::fromJSON(oai_resolve_path(params$manifest,repo))$ess_min
      if(!identical(mixedgp_parse_ess_min(params$ess_min),mixedgp_parse_ess_min(declared)))
        stop('--ess-min must match the supplied manifest; set ess_min in that manifest.')
    }
    env<-oai_load_api(repo)
    return(env$oai_screen_run(repo,list(manifest=oai_resolve_path(params$manifest,repo),
      cores=oai_integer(get('cores',4L),'cores'),plan=isTRUE(params$plan),reuse_only=isTRUE(params$reuse_only),run_folds=params$run_folds),group))
  }
  fields<-c('chains','iter','burn','stages','folds','cells','draws_per_chain','first_checkpoint','checkpoint_step',
    'baseline_profile','sampler_o0','sampler_o20','sampler_o50','representative_fold','appendix',
    'functional_draws_per_chain','n_latent','run_id','ess_min')
  s<-do.call(oai_settings,c(list(profile=get('profile','pilot'),cores=get('cores',4L)),params[intersect(names(params),fields)]))
  if(length(params$run_folds)) {
    selected<-oai_integer_list(params$run_folds,'run-folds')
    if(anyDuplicated(selected)||any(!selected%in%s$folds))stop('run-folds must be distinct members of the frozen folds')
    params$run_folds<-selected
  }
  dataset<-get('dataset','oai')
  if(!grepl('^[A-Za-z0-9][A-Za-z0-9_-]*$',dataset))stop('Invalid dataset name')
  default_input<-oai_default_input(dataset,repo)
  input<-oai_resolve_path(get('data_root',get('data_dir',default_input)),repo)
  out<-oai_resolve_path(get('output',get('output_dir',oai_default_output(dataset,s$profile,repo))),repo)
  if(group=='competitors' && file.exists(oai_path(out,'application_manifest.json'))) {
    saved<-jsonlite::fromJSON(oai_path(out,'application_manifest.json'))
    s$chains<-as.integer(saved$chains)
  }
  if(isTRUE(params$plan)&&(!dir.exists(input)||!all(file.exists(oai_input_files(oai_plan_input_root(input),1L,'O0'))))) {
    cat('Dataset:',dataset,'; data:',input,'\n');print(s[c('profile','chains','iter','burn','cores')])
    cat('Supply private inputs before fitting. No computation started.\n');return(invisible(s))
  }
  report_source<-if(!is.null(params$report_source))oai_resolve_path(params$report_source,repo)else NULL
  saved_only<-group=='report'||isTRUE(params$reuse_only)
  if(saved_only&&!isTRUE(params$plan)) {
    if(length(params$run_folds)) {
      selected<-oai_integer_list(params$run_folds,'run-folds')
      if(anyDuplicated(selected)||any(!selected%in%s$folds))stop('run-folds must be distinct members of the frozen folds')
      s$folds<-selected;attr(s,'explicit_selection')<-union(attr(s,'explicit_selection'),'folds')
      if(!s$representative_fold%in%s$folds)s$representative_fold<-s$folds[1L]
    }
    default_report<-if(group=='report')out else oai_path(out,'reports',group)
    return(oai_execution_lock(oai_workflow_root(out),group,function()oai_render_report(repo,input,out,s,report_source,
      report_dir=oai_resolve_path(get('report_dir',default_report),repo),group=group,
      filename=if(group!='report'&&oai_storage_flat(out))'index.html'else'oai.html')))
  }
  api<-oai_load_api(repo);inp<-api$oai_plan_inputs(input);p<-inp$cells[[1]]$protocol
  api$application$oai_output_check(out,repo)
  if(!isTRUE(params$plan))oai_check_legacy_lock(out)
  initialize<-function() {
  input<-input;inp<-inp;s<-s
  if(!isTRUE(params$plan)) {
    oai_storage_init(out,s$stages)
    if(oai_storage_flat(out)) {
      input<-oai_snapshot_inputs(input,out,inp$files)
      inp<-api$oai_plan_inputs(input)
    }
  }
  manifest<-oai_path(out,'application_manifest.json')
  if(file.exists(manifest)) {
    prior<-jsonlite::fromJSON(manifest,simplifyVector=FALSE)
    # Existing runs retain their checkpoint cadence unless explicitly changed.
    for(n in c('first_checkpoint','checkpoint_step'))if(is.null(params[[n]])&&!is.null(prior[[n]]))s[[n]]<-as.integer(prior[[n]])
    if(group=='competitors')s$chains<-as.integer(prior$chains)
  }else if(group=='competitors')s$chains<-4L
  files<-setNames(lapply(inp$files,function(x)digest::digest(file=x,algo='sha256')),inp$files)
  m<-list(schema='oai-application-v1',author='Xinye Chen',profile=s$profile,
    package_version=as.character(packageVersion('eivGP')),
    baseline_profile=s$baseline_profile,trial_id=dataset,smoke=s$profile=='quick',
    output_root=oai_workflow_root(out),cores=s$cores,chains=s$chains,stages=s$stages,burn=s$burn,
    first_checkpoint=s$first_checkpoint,checkpoint_step=s$checkpoint_step,draws_per_chain=s$draws_per_chain,ess_min=s$ess_min,
    samplers=list(O0=s$sampler_o0,O20=s$sampler_o20,O50=s$sampler_o50),folds=s$folds,cells=s$cells,
    cases=list(list(id=p$plan_id,input=input,output=out,C=p$C_order,N=inp$N,q=p$q)),frozen_inputs=files)
  if(!isTRUE(params$plan)) {
    if(file.exists(manifest)) {
      a<-m;b<-prior;a$cores<-b$cores<-NULL
      if(b$schema=='adni-application-v1')b$schema<-'oai-application-v1'
      canonical<-function(x)jsonlite::fromJSON(jsonlite::toJSON(x,auto_unbox=TRUE,null='null'),simplifyVector=FALSE)
      if(!identical(canonical(a),canonical(b)))stop('Saved OAI settings/inputs differ; select a new output directory. ',
        paste(all.equal(canonical(a),canonical(b)),collapse='; '))
    }else api$atomic_json(m,manifest)
  }
  list(m=m,input=input,s=s,manifest=manifest)
  }
  prepared<-if(isTRUE(params$plan))initialize()else oai_with_lock(
    base::file.path(dirname(out),paste0('.',basename(out),'.initialization.lock')),initialize,wait=TRUE)
  m<-prepared$m;input<-prepared$input;s<-prepared$s;manifest<-prepared$manifest
  if(isTRUE(params$plan)) {
    e<-oai_load_api(repo)
    plan_manifest<-m
    if(length(params$run_folds)) {
      selected<-oai_integer_list(params$run_folds,'run-folds')
      if(anyDuplicated(selected)||any(!selected%in%s$folds))stop('run-folds must be distinct members of the frozen folds')
      plan_manifest$folds<-selected
    }
    tasks<-e$oai_screen_tasks(plan_manifest,group)
    print(tasks,row.names=FALSE);print(s[c('profile','chains','iter','burn','cores')])
    cat('Resident fit jobs:',oai_resident_workers(s$cores,'fit',nrow(tasks),if(group=='competitors')1L else s$chains),
      '; resident evaluation jobs:',oai_resident_workers(s$cores,'evaluate',nrow(tasks),if(group=='competitors')1L else s$chains),'\n')
    return(invisible(m))
  }
  oai_execution_lock(oai_workflow_root(out),group,function(){
  e<-oai_load_api(repo)
  result<-e$oai_screen_run(repo,list(manifest=manifest,cores=s$cores,plan=FALSE,reuse_only=FALSE,
    run_folds=params$run_folds,lock_held=TRUE),group)
  if(length(params$run_folds)) {
    s$folds<-oai_integer_list(params$run_folds,'run-folds')
    attr(s,'explicit_selection')<-union(attr(s,'explicit_selection'),'folds')
    if(!s$representative_fold%in%s$folds)s$representative_fold<-s$folds[1L]
  }
  phase_status<-attr(result,'phase_status')
  if(!is.null(phase_status))oai_publish_progress(oai_workflow_root(out),group,'REPORTING',phase_status)
  tryCatch(oai_render_report(repo,input,out,s,report_dir=oai_path(out,'reports',group),group=group,
    filename=if(oai_storage_flat(out))'index.html'else'oai.html'),error=function(e){
      if(!is.null(phase_status))oai_publish_progress(oai_workflow_root(out),group,'REPORT_FAILED',phase_status)
      stop(e)
    })
  if(!is.null(phase_status))oai_publish_progress(oai_workflow_root(out),group,'COMPLETE',phase_status)
  cat('[OAI ',group,'] COMPLETE | report: ',oai_path(out,'reports',group,'index.html'),'\n',sep='');flush.console()
  invisible(result)
  })
}

# Author: Xinye Chen. Base-R C/X-only folds and nested stratified observation masks.
oai_seed <- function(seed,label) {
  h<-digest::digest(paste(seed,label,sep='/'),algo='sha256',serialize=FALSE)
  as.integer(strtoi(substr(h,1,7),16L)%%(2^30-1)+1)
}
oai_nested_masks <- function(d,strata='C1',seed,fractions=c(.2,.5)) {
  key<-do.call(paste,c(d[,strata,drop=FALSE],sep='/'));groups<-split(seq_len(nrow(d)),key)
  ranked<-lapply(names(groups),function(g){set.seed(oai_seed(seed,g));v<-groups[[g]];v[sample.int(length(v))]})
  names(ranked)<-names(groups);counts<-lengths(groups);previous<-integer(length(groups));out<-list()
  for(fraction in fractions) {
    n<-as.integer(floor(nrow(d)*fraction+.5));ideal<-n*counts/sum(counts)
    quota<-pmax(previous,pmin(counts,floor(ideal)))
    while(sum(quota)<n){candidates<-which(quota<counts);j<-candidates[which.max(ideal[candidates]-quota[candidates])];quota[j]<-quota[j]+1L}
    while(sum(quota)>n){candidates<-which(quota>previous);j<-candidates[which.max(quota[candidates]-ideal[candidates])];quota[j]<-quota[j]-1L}
    selected<-unlist(Map(function(v,k)head(v,k),ranked,quota),use.names=FALSE)
    out[[paste0('O',round(100*fraction))]]<-sort(d$RID[selected]);previous<-quota
  }
  stopifnot(all(out$O20%in%out$O50));out
}
oai_balanced_folds <- function(d,C,X,seed) {
  n<-nrow(d);features<-list()
  for(c in C)for(level in sort(unique(d[[c]]))) {
    v<-as.numeric(d[[c]]==level);features[[length(features)+1L]]<-v/sqrt(max(1,sum(v)))
  }
  for(x in X) {
    v<-d[[x]];features[[length(features)+1L]]<-.15*(v-mean(v))/sd(v)
    edges<-unique(quantile(v,c(.2,.4,.6,.8)));bins<-findInterval(v,edges)
    for(level in unique(bins)){z<-as.numeric(bins==level);features[[length(features)+1L]]<-.5*z/sqrt(sum(z))}
  }
  features<-do.call(cbind,features);cap<-rep(n%/%5,5);if(n%%5)cap[seq_len(n%%5)]<-cap[seq_len(n%%5)]+1L
  target<-outer(cap/n,colSums(features));totals<-matrix(0,5,ncol(features));sizes<-integer(5);fold<-integer(n)
  rarity<-apply(sapply(C,function(c)as.integer(table(d[[c]])[as.character(d[[c]])])),1,min)
  set.seed(seed);tie<-runif(n);order<-order(rarity,tie)
  for(i in order) {
    available<-which(sizes<cap)
    cost<-vapply(available,function(k)sum(2*(totals[k,]-target[k,])*features[i,]+features[i,]^2),0)
    best<-available[which.min(cost)];fold[i]<-best;sizes[best]<-sizes[best]+1L;totals[best,]<-totals[best,]+features[i,]
  }
  fold
}
oai_prepare_folds <- function(cohort_file,metadata_file,output_dir) {
  if(!requireNamespace('jsonlite',quietly=TRUE)||!requireNamespace('digest',quietly=TRUE))stop('Install jsonlite and digest')
  if(dir.exists(oai_path(output_dir,'MC'))||
     (dir.exists(output_dir)&&any(grepl('^fold[0-9]+_',list.files(output_dir)))))
    stop('Prepared inputs already exist; never overwrite or reroll frozen folds')
  meta<-jsonlite::fromJSON(metadata_file);d<-read.csv(cohort_file,stringsAsFactors=FALSE)
  for(name in c('fold_seed','calibration_seed','test_seed','model_seed')) {
    v<-meta[[name]]
    if(length(v)!=1L||!is.finite(v)||v<0||v!=floor(v)||v>2^30-1)stop('Supply a fixed integer seed: ',name)
  }
  C<-paste0('C',seq_along(meta$C_order));X<-paste0('X',seq_along(meta$X_order));U<-paste0('U',seq_along(meta$U_order))
  fields<-c('RID','Y',X,U,C)
  if(!all(fields%in%names(d))||anyNA(d[,fields])||any(!is.finite(as.matrix(d[,setdiff(fields,'RID')]))))stop('Invalid prepared cohort columns')
  if(anyDuplicated(d$RID)||nrow(d)<10)stop('Expected unique participant IDs and sufficient rows')
  if(length(meta$m_vec)!=length(C)||length(U)<1||length(C)<length(U))stop('Invalid coordinate dimensions')
  for(j in seq_along(C))if(any(!d[[C[j]]]%in%seq_len(meta$m_vec[j])))stop('C must contain native integer grades')
  for(x in X)if(length(unique(d[[x]]))<=2||sd(d[[x]])<=0)stop('X must be continuous; fix categorical population attributes before preparing data')
  strata<-meta$stratification_columns;if(is.null(strata))strata<-'C1'
  if(any(!strata%in%C))stop('Stratification uses native C only')
  d<-d[order(d$RID),fields]
  original<-read.csv(cohort_file)
  if('fold'%in%names(original)) {
    d$fold<-original$fold[match(d$RID,original$RID)]
    if(any(!d$fold%in%1:5)||!all(1:5%in%d$fold))stop('Invalid supplied folds')
  }else d$fold<-oai_balanced_folds(d,C,X,meta$fold_seed)
  for(f in 1:5)for(j in seq_along(C))if(!all(seq_len(meta$m_vec[j])%in%d[d$fold!=f,C[j]]))
    stop('Every declared native grade must be represented in each training fold; do not silently merge levels or reroll')
  for(f in 1:5) {
    tr<-d[d$fold!=f,fields];te<-d[d$fold==f,fields];tr$calibration<-te$calibration<-0L
    cal<-oai_nested_masks(tr,strata,oai_seed(meta$calibration_seed,paste0('fold',f)))
    mask<-oai_nested_masks(te,strata,oai_seed(meta$test_seed,paste0('fold',f)))
    for(cell in c('O20','O50'))if(qr(cbind(1,as.matrix(tr[tr$RID%in%cal[[cell]],U,drop=FALSE])))$rank<length(U)+1L)
      stop('The frozen calibration vectors have insufficient affine rank; no outcome-guided mask reroll is performed')
    for(cell in c('O0','O20','O50')) {
      dir.create(output_dir,recursive=TRUE,showWarnings=FALSE)
      paths<-base::file.path(output_dir,paste0(sprintf('fold%02d_%s_',f,cell),c('PROTOCOL.json','train.csv','test.csv')))
      train<-tr;train$calibration<-as.integer(train$RID%in%cal[[cell]])
      write.csv(train,paths[2],row.names=FALSE);write.csv(te,paths[3],row.names=FALSE)
      seed_names<-c('production','posterior','prediction','mixed_prediction','prospective_U','UC_GP','LVGP','EzGP','EzGP_draws')
      seeds<-setNames(lapply(seed_names,function(n)oai_seed(meta$model_seed,paste(f,if(n%in%c('UC_GP','LVGP','EzGP','EzGP_draws'))'baseline'else cell,n,sep='/'))),seed_names)
      p<-list(author='Xinye Chen',plan_id=meta$dataset_id,fold=f,O=cell,N=nrow(d),training_N=nrow(tr),test_N=nrow(te),
        observed_training_U=sum(train$calibration),p=length(X),q=length(C),d=length(U),X_order=meta$X_order,U_order=meta$U_order,
        C_order=meta$C_order,m_vec=meta$m_vec,C_levels=lapply(meta$m_vec,seq_len),
        ident=if(cell=='O0')'lower_triangular'else'none',standardize_U=cell!='O0',test_U_available_RID=if(cell=='O0')integer()else mask[[cell]],
        seeds=seeds,stratification_columns=strata,calibration_design='Proportional stratified sampling; fixed seeded order within stratum; nested20/50',
        fold_design=if('fold'%in%names(original))'Supplied frozen folds'else'C/X-only balanced folds; no Y or U used',
        train_md5=unname(tools::md5sum(paths[2])),test_md5=unname(tools::md5sum(paths[3])))
      jsonlite::write_json(p,paths[1],auto_unbox=TRUE,pretty=TRUE)
    }
  }
  jsonlite::write_json(meta,oai_path(output_dir,'DATASET.json'),auto_unbox=TRUE,pretty=TRUE)
  invisible(d[,c('RID','fold')])
}


  },e)
  e$repository_dir<-repo;e$api<-e;e$application<-e
  e$queue<-setNames(lapply(c('mixedgp_stage_workers','mixedgp_schedule','mixedgp_bind','mixedgp_table_html','mixedgp_competitor_controls'),function(n)get(n,ns)),c('mixedgp_stage_workers','mixedgp_schedule','mixedgp_bind','mixedgp_table_html','mixedgp_competitor_controls'))
  for(n in c('mixedgp_budget_yield','mixedgp_budget_worker','mixedgp_budget_enter',
    'mixedgp_budget_cores','eivgp_module_identity'))e[[n]]<-get(n,ns)
  e$atomic_json<-function(value,path) {
    dir.create(dirname(path),recursive=TRUE,showWarnings=FALSE)
    p<-tempfile(tmpdir=dirname(path));on.exit(unlink(p))
    jsonlite::write_json(value,p,auto_unbox=TRUE,pretty=TRUE,null='null',na='null')
    stopifnot(file.rename(p,path))
  }
  e$atomic_csv<-get('mixedgp_atomic_csv',ns);e$atomic_rds<-get('mixedgp_atomic_rds',ns)
  e$analyze_cell<-e$oai_analyze_cell;e$baseline_cell<-e$oai_baseline_cell
  e$oai_load_api<-e$oai_screen_api<-e$oai_workflow_api<-function(repo=NULL)e
  e$run_oai_application<-function(repository_dir=repo,params=list(),group=c('report','eivGP','competitors'))
    e$oai_method_run(repository_dir,params,match.arg(group))
  oai_configure_api(e,repo)
}

oai_configure_api <- function(e,repo) {
  example<-file.path(repo,'data/oai/DATASET.example.json')
  e$oai_metadata<-if(file.exists(example))jsonlite::fromJSON(example,simplifyVector=FALSE)else list()
  e$oai_default_input<-function(dataset,repo)base::file.path(repo,'data',if(dataset=='oai')'oai'else paste0('oai_',dataset))
  source_files<-e$oai_source_files
  e$oai_source_files<-function(repo)unique(c(source_files(repo),base::file.path(repo,'tools',
    c('oai-tools.R','oai.R','oai-eivGP.R','oai-competitors.R','oai-report.R','oai-prepare-folds.R')),
    base::file.path(repo,'paper/applications/oai.Rmd')))
  settings<-e$oai_settings
  e$oai_settings<-function(profile='pilot',cores=4L,...) {
    overrides<-list(...)
    if(profile=='pilot') {
      if(!any(c('iter','stages')%in%names(overrides)))overrides$stages<-1500L
      if(is.null(overrides$checkpoint_step))overrides$checkpoint_step<-500L
    }
    if(profile=='paper'&&is.null(overrides$checkpoint_step))overrides$checkpoint_step<-5000L
    s<-do.call(settings,c(list(profile=profile,cores=cores),overrides))
    if(profile=='paper'&&is.null(overrides$first_checkpoint))
      s$first_checkpoint<-as.integer(min(s$stages[1],(s$burn%/%s$checkpoint_step+1L)*s$checkpoint_step))
    s
  }
  read_inputs<-e$oai_plan_inputs
  e$oai_plan_inputs<-function(data_dir) {
    x<-read_inputs(data_dir)
    candidates<-c(base::file.path(data_dir,'DATASET.json'),base::file.path(dirname(data_dir),'DATASET.json'))
    meta<-candidates[file.exists(candidates)]
    if(length(meta))e$oai_metadata<-jsonlite::fromJSON(meta[1],simplifyVector=FALSE)
    x
  }
  run<-e$oai_method_run
  e$oai_method_run<-function(repo,params,group) {
    if(is.null(params$dataset))params$dataset<-'oai'
    run(repo,params,group)
  }
  table_html<-e$oai_table_html
  e$oai_table_html<-function(x,caption='',note='') {
    table_labels<-c(CORE='Parameter checks passed',core_pass='Parameter checks passed',
      core_nonfinite='Undefined monitored diagnostics',strict_all_pass='All parameter checks passed')
    if(is.data.frame(x)) {
      at<-match(names(x),names(table_labels));names(x)[!is.na(at)]<-unname(table_labels[at[!is.na(at)]])
    }
    table_html(x,caption,note)
  }
  e$oai_mask_notice<-function(a=NULL) {
    meta<-if(is.null(a))e$oai_metadata else a$metadata
    selected<-isTRUE(meta$test_U_selection$outcome_selected)||
      grepl('post.hoc|outcome.selected',paste(unlist(meta$test_U_design),collapse=' '),ignore.case=TRUE)
    if(!selected)return('')
    paste('The supplied test-U masks were selected after inspecting held-out results.',
      'CC comparisons are exploratory and conditional on that selection; they are not unbiased confirmatory results.',
      'The workflow preserves the supplied masks and performs no performance-based mask search.')
  }
  e$oai_C_dictionary<-function(items,levels,labels=NULL) {
    encoding<-e$oai_metadata$grade_encoding
    native<-vapply(seq_along(items),function(i) {
      map<-encoding[[paste0('C',i)]]
      if(length(map)==levels[i])paste(unlist(map,use.names=FALSE),collapse=' / ')else 'Supplied native-grade mapping'
    },character(1))
    data.frame(Model_column=paste0('C',seq_along(items)),OAI_field=items,
      Item_meaning=oai_labels(items,labels),Clinical_grades=native,
      Model_codes=vapply(levels,function(n)paste(seq_len(n),collapse='/'),character(1)),check.names=FALSE)
  }
  save_plot<-e$oai_save_plot
  e$oai_save_plot<-function(filename,plot,...) {
    if(grepl('training_U',filename,fixed=TRUE))plot<-plot+ggplot2::labs(
      colour=oai_labels(unlist(e$oai_metadata$C_order),e$oai_metadata$C_labels)[1])
    save_plot(filename,plot,...)
  }
  plots<-e$oai_report_plots
  e$oai_report_plots<-function(a) {
    display<-a;fold<-as.integer(unlist(a$metadata$imputation_display_fold))[1]
    if(length(fold)&&!is.na(fold)&&fold%in%a$settings$folds)display$settings$representative_fold<-fold
    p<-plots(display)
    y<-if(length(a$metadata$Y_label))as.character(a$metadata$Y_label)else as.character(a$metadata$Y)
    if(!length(y))y<-'Response Y'
    u<-if(a$protocol$d==1L)oai_labels(a$protocol$U_order,a$metadata$U_labels)else 'Measured U reference'
    c1<-oai_labels(a$protocol$C_order,a$metadata$C_labels)[1]
    p$cohort<-p$cohort+ggplot2::labs(x=u,y=y,colour=c1)
    p$covariates<-p$covariates+ggplot2::labs(y=y,colour=c1)
    p$C_U<-p$C_U+ggplot2::labs(x='Ordered model code (see native-grade dictionary)',y=u)
    p$C_counts<-p$C_counts+ggplot2::labs(x='Ordered model code (see native-grade dictionary)')
    p
  }
  functionals<-e$oai_saved_functionals
  e$oai_saved_functionals<-function(a) {
    f<-functionals(a)
    y<-if(length(a$metadata$Y_label))as.character(a$metadata$Y_label)else 'Response Y'
    f$precision$units[f$precision$units=='OAI_MEM']<-y
    x<-oai_labels(a$protocol$X_order,a$metadata$X_labels)[1]
    f$mean_plot<-f$mean_plot+ggplot2::labs(x=x,y=paste('m(x,c):',y))
    f$mean_all_plot<-f$mean_all_plot+ggplot2::labs(x=x,y=paste('m(x,c):',y))
    if(a$protocol$d==1L&&nrow(f$surfaces)) {
      surface<-function(d)ggplot2::ggplot(d,ggplot2::aes(U1,mean,colour=cell,fill=cell))+
        ggplot2::geom_ribbon(ggplot2::aes(ymin=lo,ymax=hi),alpha=.15,colour=NA)+
        ggplot2::geom_line(linewidth=.7)+ggplot2::theme_bw()+
        ggplot2::labs(x=oai_labels(a$protocol$U_order,a$metadata$U_labels),y=paste('f(x,u):',y),colour='Calibration',fill='Calibration')
      f$surface_plot<-surface(f$surfaces[f$surfaces$fold==a$settings$representative_fold,])+ggplot2::facet_wrap(~cell)
      f$surface_all_plot<-surface(f$surfaces)+ggplot2::facet_grid(fold~cell)
    }
    utils::write.csv(f$precision,e$oai_path(a$report_dir,'report-assets/functionals/precision.csv'),row.names=FALSE)
    f
  }
  exports<-e$oai_report_exports
  e$oai_report_exports<-function(a) {
    exports(a)
    file<-e$oai_path(a$report_dir,'report-assets/REPORT_PROVENANCE.json')
    record<-jsonlite::fromJSON(file,simplifyVector=FALSE)
    record$application<-'OAI';record$test_U_design<-a$metadata$test_U_design
    record$test_U_selection<-a$metadata$test_U_selection
    record$test_U_selection_notice<-e$oai_mask_notice(a)
    record$imputation_display_fold<-if(length(a$metadata$imputation_display_fold)&&
      as.integer(a$metadata$imputation_display_fold)%in%a$settings$folds)as.integer(a$metadata$imputation_display_fold)else a$settings$representative_fold
    record$reporting_source<-as.list(tools::md5sum(e$oai_source_files(a$repo)))
    e$atomic_json(record,file)
  }
  assets<-e$oai_paper_assets
  e$oai_paper_assets<-function(a) {
    m<-assets(a);folder<-e$oai_path(a$report_dir,'report-assets/paper')
    for(i in seq_len(nrow(m))) {
      old<-m$file[i];new<-gsub('_oai_','_oai_',old,fixed=TRUE)
      if(old!=new&&file.exists(base::file.path(folder,old)))stopifnot(file.rename(base::file.path(folder,old),base::file.path(folder,new)))
      oh<-sub('[.]csv$','.html',old);nh<-sub('[.]csv$','.html',new)
      if(oh!=nh&&file.exists(base::file.path(folder,oh)))stopifnot(file.rename(base::file.path(folder,oh),base::file.path(folder,nh)))
      m$file[i]<-new;m$label[i]<-gsub('_oai_','_oai_',m$label[i],fixed=TRUE)
      if(file.exists(base::file.path(folder,new)))m$MD5[i]<-unname(tools::md5sum(base::file.path(folder,new)))
    }
    utils::write.csv(m,base::file.path(folder,'artifact_manifest.csv'),row.names=FALSE);m
  }
  e$oai_render_report<-function(repo,data_dir,output_dir,settings,report_source=NULL,
      report_dir=output_dir,group='report',rmd=NULL,filename='oai.html') {
    if(!requireNamespace('rmarkdown',quietly=TRUE)||!rmarkdown::pandoc_available())stop('HTML requires rmarkdown and Pandoc')
    a<-e$oai_report_prepare(repo,data_dir,output_dir,settings,report_source,report_dir,group)
    env<-new.env(parent=globalenv());env$.oai_precomputed_analysis<-a;env$.oai_repository_dir<-repo
    if(is.null(rmd))rmd<-base::file.path(repo,'paper/applications/oai.Rmd')
    if(filename=='oai.html')filename<-'oai.html'
    rmarkdown::render(rmd,output_file=filename,output_dir=report_dir,intermediates_dir=report_dir,
      params=list(profile=settings$profile,action='report',cores=settings$cores,data_dir=data_dir,output_dir=output_dir,
        report_source=if(is.null(report_source))''else report_source,representative_fold=settings$representative_fold,appendix=settings$appendix),envir=env,quiet=TRUE)
    cat('Saved OAI HTML: ',base::file.path(report_dir,filename),'\n',sep='')
    invisible(list(html=base::file.path(report_dir,filename),analysis=a))
  }
  e$oai_settings<-e$oai_settings;e$oai_prepare_folds<-e$oai_prepare_folds
  e$oai_method_run<-e$oai_method_run;e$oai_default_output<-e$oai_default_output
  e$oai_load_api<-function(repo=NULL)e
  e
}

oai_workflow_api <- oai_screen_api <- oai_load_api
# Export application functions when sourced, retaining the fresh-environment loader.
local({
  api<-oai_load_api()
  bindings<-as.list(api)
  bindings[c('oai_load_api','oai_screen_api','oai_workflow_api')]<-NULL
  list2env(bindings,envir=parent.env(environment()))
})

oai_settings <- function(profile='pilot',cores=4L,...)oai_load_api()$oai_settings(profile,cores,...)
oai_prepare_folds <- function(cohort_file,metadata_file,output_dir)oai_load_api()$oai_prepare_folds(cohort_file,metadata_file,output_dir)
oai_method_main <- function(group,args=commandArgs(TRUE)) {
  if(identical(args,'--help')) {
    cat('Rscript tools/oai-',group,'.R --profile=paper --cores=4\n',
      '  --profile=pilot|paper|quick; --dataset=NAME; --data-root=PATH; --output=PATH\n',
      '  --cores=N is the sole resource control; new fits use min(N,4) chains.\n',
      '  --plan=true; --reuse-only=true; --iter=N --burn=N --stages=N,N\n',
      '  --folds=1:5 --cells=O0,O20,O50; --ess-min=100; --report-source=FILE\n',
      '  --run-folds=1 executes a subset of the frozen folds in the same output.\n',
      '  --checkpoint-step=5000; chain states and draws are saved for recovery.\n',
      '  EIV and competitors may run concurrently; each --cores budget is independent.\n',
      '  Run the combined report after both fitting groups finish.\n',sep='')
    return(invisible(NULL))
  }
  Sys.setenv(OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1',VECLIB_MAXIMUM_THREADS='1')
  cat('[OAI ',group,'] Initializing inputs and saved settings.\n',sep='');flush.console()
  api<-oai_load_api(.oai_repository)
  tryCatch(api$oai_method_run(.oai_repository,api$oai_cli_params(args),group),error=function(e){
    cat('[OAI ',group,'] FAILED: ',conditionMessage(e),'\n',sep='');flush.console();stop(e)
  })
}

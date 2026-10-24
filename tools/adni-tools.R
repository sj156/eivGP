# ADNI application functions. All implementation stays in this flat tools file.
.adni_repository <- local({
  files<-as.character(unlist(lapply(sys.frames(),function(frame) {
    c(if(is.character(frame$ofile))frame$ofile,if(is.character(frame$file))frame$file)
  }),use.names=FALSE))
  files<-files[basename(files)=='adni-tools.R']
  path<-if(length(files))tail(files,1L)else sub('^--file=','',grep('^--file=',commandArgs(FALSE),value=TRUE)[1])
  normalizePath(file.path(dirname(normalizePath(path)),'..'))
})
adni_load_api <- function(repo=.adni_repository) {
  repo<-normalizePath(repo)
  lib<-Sys.getenv('EIVGP_R_LIBRARY',file.path(repo,'.r-library'))
  if(nzchar(lib)&&dir.exists(lib)).libPaths(c(normalizePath(lib),.libPaths()))
  if(!requireNamespace('eivGP',quietly=TRUE)||utils::packageVersion('eivGP')<'0.4.0.9014')
    stop('Install eivGP >= 0.4.0.9014 before running this application.')
  ns<-asNamespace('eivGP');e<-new.env(parent=ns)
  evalq({
# Task folders and shared directories mirror the Ocean flat result layout.
adni_storage_flat <- function(root) {
  marker<-base::file.path(root,'_meta','storage-layout.rds')
  if(!file.exists(marker))return(FALSE)
  if(!identical(readRDS(marker)$schema,'adni-flat-v1'))stop('Unknown ADNI storage schema: ',marker)
  TRUE
}
adni_storage_init <- function(root,stages) {
  if(adni_storage_flat(root))return(invisible(root))
  existing<-if(dir.exists(root))list.files(root,all.files=TRUE,no..=TRUE)else character()
  # Existing unmarked results retain their original physical layout.
  if(length(existing))return(invisible(root))
  dir.create(base::file.path(root,'_meta'),recursive=TRUE,showWarnings=FALSE)
  mixedgp_atomic_rds(list(schema='adni-flat-v1',stages=as.integer(stages)),
    base::file.path(root,'_meta','storage-layout.rds'))
  invisible(root)
}
adni_flat_relative <- function(path,stages) {
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
adni_path <- function(...,fsep=.Platform$file.sep) {
  paths<-base::file.path(...,fsep=fsep)
  vapply(paths,function(path) {
    candidate<-path
    repeat {
      marker<-base::file.path(candidate,'_meta','storage-layout.rds')
      if(dir.exists(candidate)&&file.exists(marker)) {
        adni_storage_flat(candidate)
        if(identical(path,candidate))return(path)
        rel<-substring(path,nchar(candidate)+2L)
        return(base::file.path(candidate,adni_flat_relative(rel,readRDS(marker)$stages)))
      }
      parent<-dirname(candidate)
      if(identical(parent,candidate)||candidate=='.')return(path)
      candidate<-parent
    }
  },character(1),USE.NAMES=FALSE)
}
adni_workflow_root <- function(out)if(adni_storage_flat(out))out else base::file.path(out,'workflow')
adni_default_output <- function(dataset,profile,repo) {
  base::file.path(repo,'results',if(profile=='paper')dataset else paste(dataset,profile,sep='_'))
}
adni_snapshot_inputs <- function(input,out,files) {
  dest<-adni_path(out,'data');meta<-base::file.path(input,'DATASET.json')
  if(file.exists(meta))files<-c(files,meta)
  for(source in files) {
    target<-base::file.path(dest,adni_input_basename(source))
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

adni_source_files <- function(repo) {
  adni_path(repo,'tools',c('adni-tools.R','adni.R','adni-eivGP.R',
    'adni-competitors.R','adni-report.R','adni-prepare-folds.R'))
}
# ADNI application: data. Reuses package model and sampler functions.

adni_ezgp_mixed_matrix <- function(X, C, m_vec) {
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

adni_ezgp_cross_collisions <- function(train_x, test_x, digits = 5L) {
    train_key <- apply(round(as.matrix(train_x), digits), 1L, paste, collapse = "\r")
    test_key <- apply(round(as.matrix(test_x), digits), 1L, paste, collapse = "\r")
    sum(test_key %in% train_key)
}

adni_ezgp_predict_covariance_consistent <- function(X_new, model) {
    mixedgp_ezgp_predict(X_new, model)
}

adni_ezgp_select_tau_cv <- function(train_x, y_train, p, q, m_vec, seed, tau_fractions, cv_folds = 3L, maxeval = 100L,
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
                pred <- adni_ezgp_predict_covariance_consistent(train_x[validation_id, , drop = FALSE], fit)
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
                inner_fold = fold, n_train = sum(train_id), n_validation = sum(validation_id), rounded_cross_collisions = adni_ezgp_cross_collisions(train_x[train_id,
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

adni_ezgp_draws <- function(mean, variance, n_draw, seed) {
    if (any(!is.finite(mean)) || any(!is.finite(variance)) || any(variance < 0)) {
        stop("Invalid predictive mean or variance.")
    }
    set.seed(as.integer(seed))
    matrix(stats::rnorm(n_draw * length(mean), rep(mean, each = n_draw), rep(sqrt(variance), each = n_draw)), nrow = n_draw,
        ncol = length(mean))
}

adni_plan_require <- function() {
    for (p in c("jsonlite", "digest", "kergp", "LVGP", "EzGP")) if (!requireNamespace(p, quietly = TRUE))
        stop("Required ADNI dependency: ", p)
}

# Flat filenames encode fold and calibration; older prepared trees remain readable.
adni_input_basename <- function(path) {
    name <- basename(path)
    if (name %in% c('train.csv','test.csv','PROTOCOL.json')) {
        normalized <- gsub('\\', '/', path, fixed=TRUE)
        hit <- regmatches(normalized, regexec('(fold[0-9]+)/(O[0-9]+)/(private/)?(train[.]csv|test[.]csv|PROTOCOL[.]json)$', normalized))[[1]]
        if (length(hit)) name <- paste(hit[2], hit[3], hit[5], sep='_')
    }
    name
}
adni_default_input <- function(dataset,repo) {
    base::file.path(repo,'data',if(dataset=='adni')'adni'else paste0('adni_',dataset))
}
adni_input_files <- function(root,fold,cell) {
    prefix <- sprintf('fold%02d_%s_',fold,cell)
    flat <- base::file.path(root,paste0(prefix,c('PROTOCOL.json','train.csv','test.csv')))
    if (any(file.exists(flat))) return(flat)
    base::file.path(root,sprintf('fold%02d',fold),cell,c('PROTOCOL.json','private/train.csv','private/test.csv'))
}
adni_plan_input_root <- function(path) {
    path <- path.expand(path)
    if (!file.exists(base::file.path(path,'fold01_O0_PROTOCOL.json'))) {
        if (dir.exists(base::file.path(path,'MC'))) path <- base::file.path(path,'MC')
        else if (dir.exists(base::file.path(path,'inputs','MC'))) path <- base::file.path(path,'inputs','MC')
    }
    normalizePath(path, mustWork = TRUE)
}

adni_plan_inputs <- function(data_dir) {
    if (!requireNamespace("jsonlite", quietly = TRUE))
        stop("Install jsonlite to read ADNI protocols.")
    root <- adni_plan_input_root(data_dir)
    all <- list()
    audit <- list()
    files <- character()
    held <- list()
    reference <- NULL
    for (f in 1:5) {
        by_cell <- list()
        for (cell in c("O0", "O20", "O50")) {
            dir <- adni_path(root, sprintf("fold%02d", f), cell)
            paths <- adni_input_files(root,f,cell)
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

adni_plan_fit_args <- function(z, s) {
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

adni_plan_fit_path <- function(out, task) adni_path(out, "fits", sprintf("fold%02d", task$fold), task$cell, "fit.rds")

adni_plan_read_fit <- function(path) mixedgp_read_artifact(path)

adni_plan_save_fit <- function(fit, path) {
    dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
    tmp <- tempfile(".adni-fit-", tmpdir = dirname(path))
    on.exit(unlink(tmp), add = TRUE)
    saveRDS(mixedgp_pack_artifact(fit), tmp)
    restored <- mixedgp_read_artifact(tmp)
    if (!identical(restored, fit))
        stop("ADNI checkpoint read-back differs; previous checkpoint retained.")
    rm(restored)
    gc(FALSE)
    if (!file.rename(tmp, path))
        stop("Could not publish ADNI checkpoint; previous retained: ", path)
    invisible(path)
}

adni_plan_fit_task <- function(task, inputs, out, s, extend = FALSE) {
    z <- inputs$cells[[paste(task$fold, task$cell, sep = "/")]]
    path <- adni_plan_fit_path(out, task)
    old <- options(eivGP.checkpoint_root = NULL)
    on.exit(options(old), add = TRUE)
    fit <- if (file.exists(path))
        adni_plan_read_fit(path)
    else NULL
    if(!is.null(fit)) {
      if(!extend && fit$checkpoint$iteration>s$iter && is.null(s$ess_min))
        stop('Cannot refit an earlier stage after continuation; restore its completed evaluation receipt')
      if(fit$checkpoint$burn!=s$burn||fit$control$sampler!=s$samplers[[task$cell]])
        stop('Checkpoint sampler/burn differs')
    }
    target <- s$iter
    request <- adni_path(dirname(path), "continuation_request.rds")
    if (extend) {
        if (is.null(fit))
            stop("Missing initial fit: ", path)
        tab <- adni_plan_diagnostics(fit)
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
        fit <- do.call(fit_eivgp, c(adni_plan_fit_args(z, s), list(n_iter = first, burn = s$burn)))
        adni_plan_save_fit(fit, path)
    }
    while (fit$checkpoint$iteration < target) {
        current <- fit$checkpoint$iteration
        next_target <- min(target, (floor(current/s$checkpoint_step) + 1L) * s$checkpoint_step)
        fit <- continue_eivgp(fit, n_iter = next_target - current, parallel = s$cores > 1L, n_cores = s$cores,
            verbose = TRUE)
        adni_plan_save_fit(fit, path)
    }
    if (fit$checkpoint$burn != s$burn || length(fit$checkpoint$states) != s$chains)
        stop("Saved fit settings differ.")
    if(isTRUE(s$final_stage) && length(fit$checkpoint$states)>=4L) {
      X<-as.matrix(z$test[,z$x,drop=FALSE]);C<-as.matrix(z$test[,z$c,drop=FALSE])
      rows<-mixedgp_study2_diagnostic_rows(X,C,5L)
      panel<-list(X=X[rows,,drop=FALSE],C=C[rows,,drop=FALSE],
        U=if(task$cell!="O0")as.matrix(z$test[rows,z$u,drop=FALSE])else NULL,
        n_latent=64L,seed=z$protocol$seeds$prediction)
      fit<-mixedgp_application_ess(fit,s$ess_min,panel,adni_path(dirname(path),"scientific-ess-state.rds"),s$cores)
      adni_plan_save_fit(fit,path)
      mixedgp_save_scientific_trace(fit,adni_path(dirname(path),"scientific-trace-targets"),
        sprintf("ADNI | fold %d | %s",task$fold,task$cell))
      unlink(adni_path(dirname(path),"scientific-ess-state.rds"))
    }
    dg <- adni_plan_diagnostics(fit)
    mixedgp_atomic_csv(dg$parameters, adni_path(dirname(path), "diagnostic_parameters.csv"))
    mixedgp_atomic_csv(dg$summary, adni_path(dirname(path), "diagnostic_summary.csv"))
    invisible("DONE")
}

adni_plan_diagnostics <- function(fit) {
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

adni_plan_score <- function(draws, point, truth, training_sd) {
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

adni_plan_baseline <- function(task, inputs, out, s) {
    dir <- adni_path(out, "baselines", sprintf("fold%02d", task$fold), task$method)
    path <- adni_path(dir, "fit.rds")
    if (all(file.exists(adni_path(dir, c("fit.rds", "prediction_draws.rds")))))
        return(invisible("DONE"))
    z <- inputs$cells[[paste(task$fold, "O50", sep = "/")]]
    tr <- z$train
    te <- z$test
    p <- z$protocol
    mu <- colMeans(tr[, z$x, drop = FALSE])
    ss <- vapply(tr[, z$x, drop = FALSE], sd, 0)
    if (any(ss <= 0))
        stop("Constant ADNI training X.")
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
        ans <- do.call(fit_mixedgp_competitor, c(list(method = task$method, X = xc, y = yt, C = ct, new_X = xs,
            new_C = cs, m_vec = mvec, n_draw = 4L * s$draws_per_chain, seed = p$seeds[[if (task$method == "UC-GP") "UC_GP" else "LVGP"]]),
            controls))
        point <- yc + ys * as.numeric(ans$predictive_mean)
        draws <- yc + ys * as.matrix(ans$draws)
        model <- ans
    }
    else {
        tx <- adni_ezgp_mixed_matrix(xc, ct, mvec)
        vx <- adni_ezgp_mixed_matrix(xs, cs, mvec)
        sel <- adni_ezgp_select_tau_cv(tx, yt, ncol(xc), length(z$c), mvec, p$seeds$EzGP, if (s$profile == "paper")
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
        pp <- adni_ezgp_predict_covariance_consistent(vx, model)
        point <- yc + ys * pp$Y_hat
        draws <- yc + ys * adni_ezgp_draws(pp$Y_hat, pp$predictive_MSE, 4L * s$draws_per_chain, p$seeds$EzGP_draws)
        mixedgp_atomic_csv(sel$audit, adni_path(dir, "tau_cv_audit.csv"))
        mixedgp_atomic_csv(sel$candidates, adni_path(dir, "tau_cv_candidates.csv"))
    }
    mixedgp_atomic_rds(list(model = model, center = mu, scale = ss, y_center = yc, y_scale = ys), path)
    mm <- cbind(data.frame(fold = task$fold, cell = "OI", method = task$method, protocol = "OI"), adni_plan_score(draws,
        point, te$Y, ys))
    mixedgp_atomic_csv(mm, adni_path(dir, "metrics.csv"))
    mixedgp_atomic_rds(list(draws = draws, point = point), adni_path(dir, "prediction_draws.rds"))
    invisible("DONE")
}

adni_plan_analyze_cell <- function(task, inputs, out, s) {
    path <- adni_plan_fit_path(out, task)
    if (!file.exists(path))
        stop("Missing saved ADNI fit: ", path)
    fit <- adni_plan_read_fit(path)
    z <- inputs$cells[[paste(task$fold, task$cell, sep = "/")]]
    p <- z$protocol
    tr <- z$train
    te <- z$test
    dir <- adni_path(out, "analysis", sprintf("fold%02d", task$fold), task$cell)
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
        protocol = protocol), adni_plan_score(draws, point, te$Y, sd(tr$Y)))
    metrics <- score(yd, point, "OI")
    predictions <- data.frame(row = seq_len(nrow(te)), observed_Y = te$Y, predicted_Y = point)
    mixedgp_atomic_csv(predictions, adni_path(dir, "Y_predictions.csv"))
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
            adni_path(dir, "mixed_Y_predictions.csv"))
    }
    mixedgp_atomic_csv(metrics, adni_path(dir, "prediction_metrics.csv"))
    dg <- adni_plan_diagnostics(fit)
    mixedgp_atomic_csv(cbind(fold = task$fold, cell = task$cell, dg$summary), adni_path(dir, "convergence.csv"))
    mixedgp_atomic_csv(dg$parameters, adni_path(dir, "diagnostic_parameters.csv"))
    scale <- if (any(tr$calibration == 1L))
        "raw"
    else "model"
    imp <- impute_eivgp(fit, rows = seq_len(nrow(tr)), draw_ids = ids, scale = scale)
    mixedgp_atomic_rds(list(draws = imp, scale = scale, draw_ids = ids, coordinates = p$U_order), adni_path(dir,
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
            tables[[j]] <- cbind(data.frame(fold = task$fold, cell = task$cell, coordinate = p$U_order[j]), adni_plan_score(ud[,
                hidden, drop = FALSE], up[hidden], truth[hidden], sd0), calibration_mean_RMSE = sqrt(mean((mean(truth[cal]) -
                truth[hidden])^2)))
            plots[[j]] <- data.frame(fold = task$fold, cell = task$cell, coordinate = p$U_order[j], row = seq_len(nrow(tr)),
                calibration = tr$calibration, common_hidden = common, true_U = truth, posterior_mean = up, q025 = qq[1,
                  ], q975 = qq[2, ])
        }
        mixedgp_atomic_csv(do.call(rbind, tables), adni_path(dir, "U_metrics.csv"))
        mixedgp_atomic_csv(do.call(rbind, plots), adni_path(dir, "U_predictions.csv"))
    }
    if (task$fold == s$folds[1]) {
        series <- mixedgp_v030_raw_series(fit)
        parameters <- head(grep("^(V$|r$|theta_x|tau|A\\[)", names(series), value = TRUE), 6L)
        if (length(parameters) && length(fit$checkpoint$states) >= 2L)
            ggplot2::ggsave(adni_path(dir, "mcmc_trace.png"), plot_mcmc(fit, parameters), width = 10, height = 7,
                dpi = 130)
        else if (length(fit$checkpoint$states) < 2L)
            writeLines("Single-chain workflow check: multi-chain trace diagnostics are unavailable; convergence is not established.",
                adni_path(dir, "MCMC_TRACE_UNAVAILABLE.txt"))
    }
    invisible("DONE")
}

# ADNI application: evaluation. Reuses package model and sampler functions.

adni_save_plot <- function(filename, plot, width, height, dpi = 130, ...) {
    input <- paste0(filename, ".plot.rds")
    saveRDS(plot, input)
    code <- paste0("source(", deparse(adni_path(repository_dir,"tools/adni.R")),
      ");adni_load_api()$adni_render_plot(commandArgs(TRUE))")
    args <- c('-e', shQuote(code), shQuote(input), shQuote(filename), width, height, dpi)
    rc <- system2(adni_path(R.home("bin"), "Rscript"), args, stdout = paste0(filename, ".render.log"), stderr = paste0(filename,
        ".render.log"))
    stopifnot(rc == 0L, file.exists(filename), file.info(filename)$size > 1000)
    unlink(input)
    invisible(filename)
}

adni_replace_plot_device <- function(x) {
    if (is.call(x) && identical(x, quote(ggplot2::ggsave)))
        return(as.name("adni_save_plot"))
    if (is.call(x))
        for (i in seq_along(x)) x[[i]] <- adni_replace_plot_device(x[[i]])
    x
}

adni_reference_points <- function(calibration, targets, c_names, u_name) {
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

adni_analyze_cell <- function(task, inputs, out, s, stage) {
    z <- inputs$cells[[paste(task$fold, task$cell, sep = "/")]]
    tr <- z$train
    te <- z$test
    p <- z$protocol
    stageout <- adni_path(out, stage)
    dir <- adni_path(stageout, "analysis", sprintf("fold%02d", task$fold), task$cell)
    dir.create(dir, recursive = TRUE, showWarnings = FALSE)
    cache <- new.env(parent = emptyenv())
    env <- new.env(parent = environment(api$adni_plan_analyze_cell))
    env$adni_plan_fit_path <- function(ignored, task) api$adni_plan_fit_path(out, task)
    env$adni_plan_read_fit <- function(path) {
        cache$fit <- api$adni_plan_read_fit(path)
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
            atomic_rds(list(draws = r, draw_ids = a$draw_ids), adni_path(dir, "Y_draws.rds"))
        }
        else {
            cache$exact <- r
            atomic_rds(list(draws = r, draw_ids = a$draw_ids), adni_path(dir, "known_U_Y_draws.rds"))
        }
        r
    }
    env$adni_save_plot <- adni_save_plot
    fun <- api$adni_plan_analyze_cell
    body(fun) <- adni_replace_plot_device(body(fun))
    environment(fun) <- env
    fun(task, inputs, stageout, s)
    fit <- cache$fit
    stopifnot(fit$data$p == length(z$x), fit$data$q == p$q, fit$data$d == length(z$u), all(is.na(fit$data$U_obs_raw[tr$calibration ==
        0, , drop = FALSE])))
    met <- read.csv(adni_path(dir, "prediction_metrics.csv"))
    met$stage <- stage
    met$iterations <- s$iter
    met$sampler <- fit$control$sampler
    if (exists("exact", cache, inherits = FALSE)) {
        ids <- which(z$available)
        cm <- attr(cache$exact, "conditional_means", exact = TRUE)
        stopifnot(!is.null(cm))
        extra <- cbind(data.frame(fold = task$fold, cell = task$cell, method = "EIV-GP", protocol = "CC-known"),
            api$adni_plan_score(cache$exact, colMeans(cm), te$Y[ids], sd(tr$Y)))
        extra$stage <- stage
        extra$iterations <- s$iter
        extra$sampler <- fit$control$sampler
        met <- rbind(met, extra)
    }
    atomic_csv(met, adni_path(dir, "prediction_metrics.csv"))
    dg <- read.csv(adni_path(dir, "convergence.csv"))
    dg$stage <- stage
    dg$sampler <- fit$control$sampler
    if (!is.null(fit$diagnostics$hmc))
        atomic_csv(fit$diagnostics$hmc, adni_path(dir, "HMC_DIAGNOSTICS.csv"))
    atomic_csv(dg, adni_path(dir, "convergence.csv"))
    cm <- attr(cache$pure, "conditional_means", exact = TRUE)
    per <- nrow(cm)%/%s$chains
    target_diagnostics <- do.call(rbind, lapply(seq_len(ncol(cm)), function(j) {
        a <- matrix(cm[, j], per, s$chains)
        data.frame(RID = te$RID[j], predicted_Y = mean(a), selected_Rhat = posterior::rhat(a), selected_ESS = posterior::ess_bulk(a),
            selected_MCSE = posterior::mcse_mean(a))
    }))
    atomic_csv(target_diagnostics, adni_path(dir, "PREDICTION_FUNCTION_DIAGNOSTICS.csv"))
    if (task$cell != "O0") {
        cal <- tr[tr$calibration == 1, , drop = FALSE]
        hidden <- which(tr$calibration == 0)
        ui <- read.csv(adni_path(dir, "U_predictions.csv"))
        um <- read.csv(adni_path(dir, "U_metrics.csv"))
        for (j in seq_along(z$u)) {
            truth <- tr[[z$u[j]]]
            ref <- adni_reference_points(cal, tr[hidden, , drop = FALSE], z$c, z$u[j])
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
        atomic_csv(um, adni_path(dir, "U_metrics.csv"))
        atomic_csv(ui, adni_path(dir, "U_predictions.csv"))
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
            ref <- adni_reference_points(cal, te, z$c, z$u[j])
            qq <- apply(ud, 2, quantile, c(0.025, 0.975))
            pm[[j]] <- cbind(data.frame(stage = stage, fold = task$fold, cell = task$cell, coordinate = p$U_order[j]),
                api$adni_plan_score(ud, point, truth, sd(tr[[z$u[j]]])), calibration_mean_RMSE = sqrt(mean((ref$global -
                  truth)^2)), C1_mean_RMSE = sqrt(mean((ref$C1 - truth)^2)), joint_C_mean_RMSE = sqrt(mean((ref$joint -
                  truth)^2)), joint_C_fallback_N = ref$fallback)
            pp[[j]] <- data.frame(RID = te$RID, coordinate = p$U_order[j], true_U = truth, posterior_mean = point,
                q025 = qq[1, ], q975 = qq[2, ])
        }
        atomic_csv(do.call(rbind, pm), adni_path(dir, "PROSPECTIVE_U_metrics.csv"))
        atomic_csv(do.call(rbind, pp), adni_path(dir, "PROSPECTIVE_U_predictions.csv"))
        ui$role <- factor(ifelse(ui$calibration == 1, "Observed calibration", "Masked training U"), levels = c("Masked training U",
            "Observed calibration"))
        plot <- ggplot(ui, aes(true_std, mean_std, colour = factor(C1))) + geom_abline(slope = 1, intercept = 0,
            linetype = 2, colour = "grey50") + geom_linerange(data = ui[ui$calibration == 0, ], aes(ymin = lo_std,
            ymax = hi_std), alpha = 0.15, linewidth = 0.3) + geom_point(aes(shape = role), size = 1.2, alpha = 0.65) +
            scale_shape_manual(values = c(16, 1)) + facet_wrap(~coordinate) + theme_bw(base_size = 11) + theme(legend.position = "bottom") +
            labs(x = "Reference U (training SD units; evaluation only)", y = "Posterior mean (training SD units)",
                colour = "MEMORY3", shape = NULL, title = paste(p$plan_id, "fold", task$fold, task$cell, stage),
                subtitle = paste(fit$control$sampler, "; CORE pass:", dg$core_pass))
        adni_save_plot(adni_path(dir, "training_U_imputation_std.png"), plot, width = 10, height = 5.3, dpi = 140,
            bg = "white")
    }
    else writeLines("O0 is unanchored: no physical-U RMSE or truth-based alignment.", adni_path(dir, "O0_SCALE_LIMITATION.txt"))
    atomic_json(list(state = "PASS", stage = stage, iterations = s$iter, primary_no_test_U_Y = TRUE, hidden_training_U_masked = TRUE,
        calibration_excluded_from_U_scoring = TRUE, prospective_U_uses_C_not_test_Y = TRUE, sampler = fit$control$sampler),
        adni_path(dir, "EVALUATION_VALIDATION.json"))
    rm(fit, cache)
    gc(FALSE)
    invisible("DONE")
}

adni_baseline_cell <- function(task, inputs, out, s, cc = FALSE) {
    source <- inputs$cells[[paste(task$fold, if (cc)
        task$cell
    else "O50", sep = "/")]]
    tr <- source$train
    te <- source$test
    root <- if (cc)
        adni_path(out, "CC", task$cell)
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
    api$adni_plan_baseline(task, inp, root, s)
    dir <- adni_path(root, "baselines", sprintf("fold%02d", task$fold), task$method)
    saved <- readRDS(adni_path(dir, "prediction_draws.rds"))
    selected <- if (cc)
        which(source$available)
    else seq_len(nrow(te))
    mm <- cbind(data.frame(stage = "BASELINE", fold = task$fold, cell = if (cc)
        task$cell
    else "OI", method = paste0(task$method, if (cc)
        "-CC"
    else ""), protocol = if (cc)
        "CC-known"
    else "OI"), api$adni_plan_score(saved$draws[, selected, drop = FALSE], saved$point[selected], te$Y[selected],
        sd(tr$Y)))
    mm$training_N <- if (cc)
        sum(tr$calibration == 1)
    else nrow(tr)
    atomic_csv(mm, adni_path(dir, "metrics.csv"))
    atomic_csv(data.frame(RID = te$RID, observed_Y = te$Y, predicted_Y = saved$point, scored = seq_len(nrow(te)) %in%
        selected), adni_path(dir, "Y_predictions.csv"))
    invisible("DONE")
}

# ADNI application: execution. Reuses package model and sampler functions.

adni_output_check <- function(path,repo) {
  path<-normalizePath(path.expand(path),mustWork=FALSE);repo<-normalizePath(repo)
  if(startsWith(paste0(path,'/'),paste0(repo,'/'))&&
     !startsWith(paste0(path,'/'),paste0(repo,'/results/')))
    stop('Choose ignored results/ or an external output directory; never write to versioned source/data.')
  invisible(path)
}

adni_screen_args <- function(args) {
  p<-list(cores=5L,plan=FALSE,reuse_only=FALSE)
  for(a in args){
    if(!grepl('^--[^=]+=.+$',a))stop('Use --name=value: ',a)
    k<-gsub('-','_',sub('^--([^=]+)=.*$','\\1',a));v<-sub('^[^=]+=','',a)
    if(!k%in%c('manifest','cores','fit_jobs','evaluation_jobs','plan','reuse_only'))stop('Unknown screening option: ',k)
    if(k%in%c('plan','reuse_only')){if(!v%in%c('true','false'))stop(k,' must be true or false');v<-v=='true'}
    if(k%in%c('cores','fit_jobs','evaluation_jobs')){n<-suppressWarnings(as.numeric(v));if(!is.finite(n)||n<1||n!=floor(n)||n>.Machine$integer.max)stop('cores must be a positive integer');v<-as.integer(n)}
    if(k%in%names(p)&&!k%in%c('cores','plan','reuse_only'))stop('Repeated option: ',k)
    p[[k]]<-v
  }
  if(is.null(p$manifest))stop('Supply --manifest=/external/SCREENING.json')
  p
}

adni_screen_settings <- function(m,iter,cores) {
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

adni_screen_manifest <- function(path,api) {
  m<-jsonlite::fromJSON(normalizePath(path),simplifyVector=FALSE)
  if(!m$schema%in%c('adni-screening-v2','adni-application-v1')||!length(m$cases))stop('Unsupported ADNI manifest')
  if(m$schema=='adni-screening-v2'&&!isTRUE(m$smoke)&&
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
    adni_output_check(z$output,api$repository_dir)
    inp<-api$adni_plan_inputs(z$input)
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

adni_screen_tasks <- function(m,group,stage='') {
  rows<-list();add<-function(plan,fold,cell,method,kind)rows[[length(rows)+1L]]<<-data.frame(plan=plan,fold=as.integer(fold),cell=cell,method=method,kind=kind,stage=stage)
  for(z in m$cases)for(f in unlist(m$folds)){
    if(group=='eivGP')for(cell in unlist(m$cells))add(z$id,f,cell,'EIV-GP','eiv')
    if(group=='competitors'){
      add(z$id,f,'OI','UC-GP+LVGP+EzGP','baseline_group')
      for(cell in intersect(unlist(m$cells),c('O20','O50')))for(method in c('UC-GP','LVGP','EzGP'))add(z$id,f,cell,method,'cc')
    }
  }
  if(!length(rows))return(data.frame())
  x<-do.call(rbind,rows);x$id<-paste(x$plan,x$stage,x$kind,x$fold,x$cell,gsub('[+]','_',x$method),sep='__');x
}

adni_screen_receipt <- function(out,task,phase)adni_path(out,'receipts',paste0(phase,'_',task$id,'.rds'))

adni_screen_verified <- function(path,identity) {
  if(!file.exists(path))return(FALSE)
  r<-readRDS(path);if(!identical(r$identity,identity))stop('Receipt identity differs: ',path)
  if(length(r$files)&&(!all(file.exists(r$files))||!identical(unname(tools::md5sum(r$files)),r$md5)))stop('Completed output changed: ',path)
  TRUE
}

adni_resident_workers <- function(cores,phase,n_tasks,chains=4L,limit=NULL) {
  workers<-mixedgp_stage_workers(cores,if(phase=='fit')'fit'else'analyze',n_tasks,chains)
  if(!is.null(limit))workers<-min(workers,adni_integer(limit,paste0(phase,' jobs')))
  as.integer(workers)
}

adni_screen_phase <- function(tasks,phase,cores,api,FUN,status_dir,chains=4L,resident_limit=NULL) {
  dir.create(status_dir,recursive=TRUE,showWarnings=FALSE)
  status<-tasks;status$state<-'PENDING';status$detail<-'';status$worker_pid<-NA_integer_;status$seconds<-NA_real_
  outer<-adni_resident_workers(cores,phase,nrow(tasks),chains,resident_limit)
  publish<-function()api$atomic_csv(status,adni_path(status_dir,paste0(phase,'_status.csv')))
  publish()
  # Match the public Ocean executor: resident jobs share the full active-worker pool.
  api$mixedgp_budget_yield(function()api$queue$mixedgp_schedule(seq_len(nrow(tasks)),rep(1L,nrow(tasks)),outer,
    function(i,allocated)api$mixedgp_budget_worker(function(){
      logpath<-adni_path(status_dir,paste0(phase,'_',make.names(tasks$id[i]),'.log'));con<-file(logpath,'at')
      old_message<-getConnection(sink.number(type='message'));sink(con);sink(con,type='message')
      on.exit({sink(old_message,type='message');sink();close(con)},add=TRUE);t<-proc.time()[3]
      err<-tryCatch({FUN(tasks[i,,drop=FALSE],cores);NULL},error=function(e)conditionMessage(e))
      if(!is.null(err))cat('FAILED:',err,'\n')
      list(state=if(is.null(err))'DONE'else'FAILED',detail=if(is.null(err))''else err,seconds=as.numeric(proc.time()[3]-t))
    }),
    function(i,n,pid,active){status$state[i]<<-'RUNNING';status$worker_pid[i]<<-pid;publish()},
    function(i,r,active){if(!is.list(r)||is.null(r$state))r<-list(state='FAILED',detail='Worker exited without a result',seconds=NA_real_)
      status$state[i]<<-r$state;status$detail[i]<<-r$detail;status$seconds[i]<<-r$seconds;publish()}))
  status
}

adni_screen_collect <- function(m,api) {
  pred<-conv<-u<-pros<-list()
  for(z in m$cases){
    for(mark in list.files(z$output,'^EVALUATION_VALIDATION[.]json$',recursive=TRUE,full.names=TRUE)){
      dir<-dirname(mark)
      for(pair in list(c('prediction_metrics.csv','pred'),c('convergence.csv','conv'),c('U_metrics.csv','u'),c('PROSPECTIVE_U_metrics.csv','pros'))){
        p<-adni_path(dir,pair[1]);if(!file.exists(p))next
        val<-cbind(plan=z$id,read.csv(p));dest<-get(pair[2]);dest[[length(dest)+1L]]<-val;assign(pair[2],dest)
      }
    }
    files<-list.files(z$output,'^metrics[.]csv$',recursive=TRUE,full.names=TRUE)
    for(p in files){x<-read.csv(p);if('stage'%in%names(x)&&all(x$stage=='BASELINE'))pred[[length(pred)+1L]]<-cbind(plan=z$id,x)}
  }
  bind<-function(x)if(length(x))api$queue$mixedgp_bind(x)else data.frame()
  tables<-list(prediction=bind(pred),convergence=bind(conv),training_U=bind(u),prospective_U=bind(pros))
  dest<-if(adni_storage_flat(m$output_root))adni_path(m$output_root,'reports','prediction_report')else adni_path(m$output_root,'reports');dir.create(dest,recursive=TRUE,showWarnings=FALSE)
  for(k in names(tables))if(nrow(tables[[k]]))api$atomic_csv(tables[[k]],adni_path(dest,paste0(k,'_fold_metrics.csv')))
  p<-tables$prediction;summary<-list()
  if(nrow(p))for(ix in split(seq_len(nrow(p)),interaction(p$plan,p$stage,p$cell,p$method,p$protocol,drop=TRUE))){
    x<-p[ix,,drop=FALSE];v<-x$RMSE_SD
    summary[[length(summary)+1L]]<-data.frame(x[1,c('plan','stage','cell','method','protocol')],completed_folds=length(unique(x$fold)),
       RMSE_SD_mean=mean(v),RMSE_SD_SD=if(length(v)>1)sd(v)else NA_real_)
  }
  five<-bind(summary);if(nrow(five))api$atomic_csv(five,adni_path(dest,'fivefold_prediction.csv'))
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
  stable<-bind(stability);if(nrow(stable))api$atomic_csv(stable,adni_path(dest,'stage_stability.csv'))
  lines<-c('# ADNI screening: saved results',paste('Author: Xinye Chen. Trial',m$trial_id),
    'Partial fold counts and CORE failures remain explicit. Budget completion does not establish convergence.',
    if(nrow(five))capture.output(print(five,row.names=FALSE))else'No completed scores yet.')
  writeLines(lines,adni_path(dest,'REPORT.md'))
  html<-c('<!doctype html><meta charset="utf-8"><title>ADNI screening</title>',
    '<h1>ADNI screening: saved results</h1><p>Author: Xinye Chen. Partial folds and CORE failures remain explicit; completing a budget does not establish convergence.</p>',
    api$queue$mixedgp_table_html(five,'Response prediction: fold means and sample SD'),
    api$queue$mixedgp_table_html(tables$convergence,'Full-chain CORE diagnostics'),
    api$queue$mixedgp_table_html(stable,'Paired stage stability; descriptive screening only'))
  writeLines(html,adni_path(dest,'index.html'));invisible(tables)
}

adni_screen_run <- function(repo,p,group) {
  api<-adni_screen_api(repo);bundle<-adni_screen_manifest(p$manifest,api);m<-bundle$manifest;inputs<-bundle$inputs
  cases<-setNames(m$cases,vapply(m$cases,`[[`,'','id'));api$B<-dirname(normalizePath(p$manifest))
  if(p$plan){
    for(g in if(group=='report')c('eivGP','competitors')else group){x<-adni_screen_tasks(m,g);print(x,row.names=FALSE)}
    cat('Public shared CPU pool:',p$cores,'; fit resident limit:',adni_resident_workers(p$cores,'fit',100L,if(is.null(m$chains))4L else m$chains,p$fit_jobs),
      '; evaluation resident limit:',adni_resident_workers(p$cores,'evaluate',100L,limit=p$evaluation_jobs),'; full-chain reservations/backfill disabled.\n');return(invisible(bundle))
  }
  root<-m$output_root
  adni_output_check(root,repo)
  dir.create(root,recursive=TRUE,showWarnings=FALSE);lock<-adni_path(root,if(adni_storage_flat(root))'.application.lock'else'.adni-application.lock')
  if(!dir.create(lock,showWarnings=FALSE))stop('ADNI workflow already active: ',lock)
  writeLines(as.character(Sys.getpid()),adni_path(lock,'pid'));on.exit(unlink(lock,recursive=TRUE),add=TRUE)
  code<-adni_source_files(repo)
  base_id<-list(schema='adni-upstream-screening-v2',trial=m$trial_id,modules=api$eivgp_module_identity(),
    application=setNames(unname(tools::md5sum(code)),substring(code,nchar(repo)+2L)),burn=m$burn,stages=unlist(m$stages),smoke=isTRUE(m$smoke),ess_min=m$ess_min,
    chains=if(is.null(m$chains))4L else as.integer(m$chains),
    samplers=if(is.null(m$samplers))list(O0='HMC',O20='HMC',O50='ESS')else m$samplers,
    baseline_profile=if(is.null(m$baseline_profile))if(isTRUE(m$smoke))'pilot'else'paper'else m$baseline_profile,
    draws_per_chain=if(is.null(m$draws_per_chain))if(isTRUE(m$smoke))5L else 250L else m$draws_per_chain)
  identities<-list()
  for(z in m$cases){
    id<-c(base_id,list(plan=z$id,inputs=unname(tools::md5sum(inputs[[z$id]]$files))))
    identities[[z$id]]<-id;path<-adni_path(z$output,'workflow_identity.rds')
    if(file.exists(path)){if(!identical(readRDS(path),id))stop('Saved workflow identity differs; use a new output root')}else if(group!='report'&&!p$reuse_only)api$atomic_rds(id,path)
  }
  required<-list()
  for(g in if(group=='report')c('eivGP','competitors')else group){
    for(st in if(g=='eivGP')paste0('S',unlist(m$stages))else'BASELINE'){
      tasks<-adni_screen_tasks(m,g,st)
      for(i in seq_len(nrow(tasks))){t<-tasks[i,];out<-cases[[t$plan]]$output
        required[[length(required)+1L]]<-list(path=adni_screen_receipt(out,t,if(g=='eivGP')'evaluate'else'fit'),identity=identities[[t$plan]])}
    }
  }
  if(group=='report'||p$reuse_only){
    if(!all(vapply(required,function(r)adni_screen_verified(r$path,r$identity),logical(1))))stop('Saved-only prerequisites incomplete; no fitting or continuation started')
    return(adni_screen_collect(m,api))
  }
  api$adni_plan_require();job<-adni_path(root,'logs',paste0(group,'_',format(Sys.time(),'%Y%m%d_%H%M%S')))
  old<-options(eivGP.core_events=adni_path(job,'core-events'));on.exit(options(old),add=TRUE)
  close<-api$mixedgp_budget_enter(p$cores);on.exit(close(),add=TRUE,after=FALSE);cores<-api$mixedgp_budget_cores()
  writeLines(paste('Effective shared worker budget',cores,'; public resident-fit scheduling.'),adni_path(root,paste0(group,'_budget.txt')))
  stages<-if(group=='eivGP')unlist(m$stages)else 0L
  for(iter in stages){
    stage<-if(group=='eivGP')paste0('S',iter)else'BASELINE';tasks<-adni_screen_tasks(m,group,stage)
    work<-function(t,n,evaluate=FALSE){
      out<-cases[[t$plan]]$output;inp<-inputs[[t$plan]];id<-identities[[t$plan]];s<-adni_screen_settings(m,if(iter)iter else max(unlist(m$stages)),n)
      phase<-if(evaluate)'evaluate'else'fit';receipt<-adni_screen_receipt(out,t,phase)
      if(adni_screen_verified(receipt,id))return(invisible(NULL))
      if(t$kind=='eiv'){
        path<-api$adni_plan_fit_path(out,t)
        if(evaluate){
          if(!file.exists(path))stop('Stage checkpoint unavailable for evaluation')
          actual<-api$adni_plan_read_fit(path)$checkpoint$iteration
          if(is.null(s$ess_min) && actual!=iter || !is.null(s$ess_min) && actual<iter)stop('Stage checkpoint unavailable for evaluation')
          s$iter<-actual
          api$analyze_cell(t,inp,out,s,stage);folder<-adni_path(out,stage,'analysis',sprintf('fold%02d',t$fold),t$cell)
        }else{
          api$adni_plan_fit_task(t,inp,out,s,extend=FALSE);folder<-NULL
        }
      }else{
        methods<-if(t$kind=='baseline_group')c('UC-GP','LVGP','EzGP')else t$method
        failures<-character()
        for(method in methods){one<-t;one$method<-method;tryCatch(api$baseline_cell(one,inp,out,s,cc=t$kind=='cc'),error=function(e)failures<<-c(failures,paste(method,conditionMessage(e))))}
        if(length(failures))stop(paste(failures,collapse='; '))
        folder<-if(t$kind=='cc')adni_path(out,'CC',t$cell,'baselines',sprintf('fold%02d',t$fold),t$method)else adni_path(out,'baselines',sprintf('fold%02d',t$fold))
      }
      files<-if(is.null(folder))character()else sort(list.files(folder,recursive=TRUE,full.names=TRUE))
      api$atomic_rds(list(identity=id,task=t,files=files,md5=unname(tools::md5sum(files)),finished=as.character(Sys.time())),receipt)
    }
    nchains<-if(is.null(m$chains))4L else as.integer(m$chains)
    status<-adni_screen_phase(tasks,'fit',cores,api,function(t,n)work(t,n),adni_path(job,stage),nchains,p$fit_jobs)
    if(group=='eivGP'&&any(status$state=='FAILED'))stop('EIV stage incomplete; checkpoints retained. Inspect ',job)
    if(group=='eivGP'){
      status<-adni_screen_phase(tasks,'evaluate',cores,api,function(t,n)work(t,n,TRUE),adni_path(job,stage),nchains,p$evaluation_jobs)
      if(any(status$state=='FAILED'))stop('EIV evaluation incomplete; continuation not started. Inspect ',job)
    }else if(any(status$state=='FAILED'))warning('Comparator failures retained; saved results remain available.')
    adni_screen_collect(m,api)
  }
  for(z in m$cases)mixedgp_application_trace_report(z$output,adni_path(z$output,'reports','scientific-traces'))
  invisible(adni_screen_collect(m,api))
}

# ADNI application: functionals. Reuses package model and sampler functions.

adni_trace_diagnostic <- function(values,chain) {
  groups<-split(values,chain);n<-min(lengths(groups));x<-do.call(cbind,lapply(groups,head,n))
  data.frame(states_per_chain=n,chains=ncol(x),mean=mean(values),posterior_SD=sd(values),
    Rhat=posterior::rhat(x),bulk_ESS=posterior::ess_bulk(x),tail_ESS=posterior::ess_tail(x),MCSE=posterior::mcse_mean(x))
}

adni_baseline_mean <- function(api,dir,X,C,mvec,method) {
  a<-readRDS(adni_path(dir,'fit.rds'));x<-sweep(sweep(X,2,a$center,'-'),2,a$scale,'/')
  if(method=='UC-GP') {
    nd<-getFromNamespace('make_mixedgp_competitor_data','eivGP')(x,C,mvec)
    mean<-as.numeric(stats::predict(a$model$fit,newdata=nd,type='UK',seCompute=FALSE,lightReturn=TRUE)$mean)
  }else if(method=='LVGP')mean<-as.numeric(LVGP::LVGP_predict(cbind(x,C),a$model$fit,MSE_on=0)$Y_hat)
  else mean<-api$adni_ezgp_predict_covariance_consistent(api$adni_ezgp_mixed_matrix(x,C,mvec),a$model)$Y_hat
  a$y_center+a$y_scale*mean
}

adni_saved_functionals <- function(a) {
  library(ggplot2)
  dest<-adni_path(a$report_dir,'report-assets','functionals');dir.create(dest,recursive=TRUE,showWarnings=FALSE)
  traces<-raw_traces<-curves<-surfaces<-precision<-integration<-list()
  code<-adni_source_files(a$repo)
  summarize_draws<-function(x)data.frame(mean=colMeans(x),lo=apply(x,2,quantile,.025),hi=apply(x,2,quantile,.975))
  for(key in names(a$paths)) {
    path<-a$paths[[key]];f<-path$fold;cell<-path$cell;z<-a$input$cells[[key]];p<-z$protocol
    cache<-adni_path(dest,paste0('fold',f,'_',cell,'.rds'))
    identity<-list(fit=unname(tools::md5sum(path$fit)),inputs=unname(tools::md5sum(a$input$files)),
      code=unname(tools::md5sum(code)),states=a$settings$functional_draws_per_chain,n_latent=a$settings$n_latent,
      analysis_modules=a$api$eivgp_module_identity())
    if(file.exists(cache)) {
      saved<-readRDS(cache)
      if(!identical(saved$identity,identity))stop('Derived functional cache differs; choose a fresh report directory: ',cache)
      result<-saved$result
    }else {
      cat('Read-only functionals fold',f,cell,'\n');flush.console()
      fit<-a$api$adni_plan_read_fit(path$fit)
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
      pr<-cbind(data.frame(fold=f,cell=cell,target='m(X reference, C profile)',units='ADNI_MEM'),adni_trace_diagnostic(cm[,target],chains))
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
          units=if(anchored)'physical U'else'unanchored model U'),adni_trace_diagnostic(v,info$chain)))
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
      raw<-adni_bind(lapply(c('V','sigma2','r'),function(n){
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
      value<-adni_baseline_mean(a$api,status$source_directory,xx,cc,p$m_vec,method)
      native[[length(native)+1L]]<-data.frame(fold=f,cell='native OI',
        method=paste0(method,if(status$state=='OPTIMIZER_FLAGGED')'*'else''),grid,mean=value,lo=NA_real_,hi=NA_real_)
    }
  }
  out<-list(curves=adni_bind(c(curves,native)),surfaces=adni_bind(surfaces),traces=adni_bind(traces),
    raw_traces=adni_bind(raw_traces),precision=adni_bind(precision),integration=adni_bind(integration))
  for(n in names(out))if(nrow(out[[n]]))write.csv(out[[n]],adni_path(dest,paste0(n,'.csv')),row.names=FALSE)
  reference_fold<-a$settings$representative_fold
  plot_curve<-function(x)ggplot(x,aes(x,mean,colour=interaction(method,cell),group=interaction(method,cell)))+
    geom_ribbon(aes(ymin=lo,ymax=hi,fill=interaction(method,cell)),alpha=.06,colour=NA,na.rm=TRUE)+
    geom_line(linewidth=.65)+facet_wrap(~C1,ncol=2)+theme_bw()+
    labs(x=a$protocol$X_order[1],y='m(x,c): ADNI_MEM',colour='Method / calibration',fill='Method / calibration')
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

# ADNI application: legacy. Reuses package model and sampler functions.

ADNI_CV_FOLDS <- 5L

ADNI_CV_SEED <- 20260914L

adni_make_folds <- function(dat) {
  stopifnot(nrow(dat) == 495L, !anyNA(dat$RID), !anyDuplicated(dat$RID),
            !anyNA(dat$R), all(dat$R %in% 0:1), sum(dat$R == 1L) == 133L)
  old_kind <- RNGkind()
  old_seed <- get0(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  on.exit({ do.call(RNGkind, as.list(old_kind))
    if (is.null(old_seed)) {
      if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
        rm(".Random.seed", envir = .GlobalEnv)
    } else assign(".Random.seed", old_seed, envir = .GlobalEnv)
  }, add = TRUE)
  RNGkind("Mersenne-Twister", "Inversion", "Rejection")
  set.seed(ADNI_CV_SEED)
  z <- dat[order(dat$RID), c("RID", "R")]
  z$fold <- NA_integer_; offset <- 0L
  # Continuing the cycle across strata gives 99 people in every fold and
  # either 26 or 27 observed-CSF participants in each held-out group.
  for (r in 0:1) {
    idx <- which(z$R == r)
    idx <- idx[sample.int(length(idx))]
    z$fold[idx] <- (offset + seq_along(idx) - 1L) %% ADNI_CV_FOLDS + 1L
    offset <- offset + length(idx)
  }
  stopifnot(all(table(z$fold) == 99L), all(table(z$fold[z$R == 1L]) %in% 26:27))
  rownames(z) <- NULL
  z
}

adni_read_cohort <- function(data_dir) {
  path <- adni_path(data_dir, "toledo_adni_cohort_n495.csv")
  if (unname(tools::md5sum(path)) != "d10b980253304620efea9f2a013166a9")
    stop("The cleaned cohort differs from the frozen ADNI cohort: ", path)
  dat <- read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
  required <- c("RID", "R", "y_centiloid", "u_csf_abeta42", "x_age_years",
                "x_female", "x_apoe4_dose", "c_ab42_ab40_3", "c_gfap_3", "diagnosis")
  stopifnot(all(required %in% names(dat)))
  dat
}

adni_core_plan <- function(total_cores, n_folds = 5L, n_chains = 4L,
                           max_fold_workers = 5L, os_type = .Platform$OS.type) {
  valid <- function(x) is.numeric(x) && length(x) == 1L && is.finite(x) && x >= 1 && x == floor(x)
  if (!all(vapply(list(total_cores, n_folds, n_chains, max_fold_workers), valid, logical(1))))
    stop("Core, fold, and chain counts must be positive integers.")
  candidates <- seq_len(min(n_folds, max_fold_workers))
  chains <- pmin(n_chains, floor(total_cores / candidates))
  active <- candidates * chains
  # Maximize concurrent chain workers; ties prefer fewer in-memory folds.
  best <- which.max(active)
  f <- candidates[best]; c <- chains[best]
  if (os_type == "windows") { f <- 1L; c <- 1L }
  list(total_cores = as.integer(total_cores), fold_workers = as.integer(f),
       chain_workers = as.integer(c), active_chain_limit = as.integer(f * c),
       n_folds = as.integer(n_folds), n_chains = as.integer(n_chains))
}

adni_print_core_plan <- function(plan) {
  cat("Total core budget:", plan$total_cores, "| concurrent folds:", plan$fold_workers,
      "| chain workers per fold:", plan$chain_workers,
      "| maximum concurrent chains:", plan$active_chain_limit, "\n")
}

adni_run_jobs <- function(jobs, FUN, plan, status_path) {
  release_budget <- mixedgp_budget_enter(plan$total_cores)
  on.exit(release_budget(),add=TRUE)
  worker <- function(i) {
    fold <- jobs$fold[i]
    dir <- adni_path(jobs$output_root[i], paste0("fold_", fold))
    dir.create(dir, recursive = TRUE, showWarnings = FALSE)
    lock <- adni_path(dir, ".fold-lock")
    if (!dir.create(lock, showWarnings = FALSE))
      return(list(ok = FALSE, fold = fold, message = paste("Fold is locked:", lock)))
    on.exit(unlink(lock, recursive = TRUE), add = TRUE)
    writeLines(paste(Sys.info()[["nodename"]], Sys.getpid()), adni_path(lock, "owner"))
    log <- file(adni_path(dir, "worker.log"), open = "at")
    sink(log); sink(log, type = "message")
    on.exit({ sink(type = "message"); sink(); close(log) }, add = TRUE)
    cat("\nWorker start:", as.character(Sys.time()), "PID", Sys.getpid(), "fold", fold, "\n")
    unlink(adni_path(dir, "WORKER_FAILURE.txt"))
    tryCatch({
      value <- FUN(fold)
      cat("Worker finished:", as.character(Sys.time()), "\n")
      list(ok = TRUE, fold = fold, value = value)
    }, error = function(e) {
      text <- paste("Fold", fold, "failed:", conditionMessage(e))
      writeLines(text, adni_path(dir, "WORKER_FAILURE.txt"))
      cat(text, "\n")
      list(ok = FALSE, fold = fold, message = text)
    })
  }
  results <- mixedgp_parallel_lapply(seq_len(nrow(jobs)),worker,
    n_cores=plan$fold_workers,.preserve_rng=TRUE)
  good <- vapply(results, function(x) is.list(x) && isTRUE(x$ok), logical(1))
  status <- data.frame(fold = jobs$fold, success = good,
    message = vapply(seq_along(results), function(i) {
      x <- results[[i]]
      if (good[i]) "finished" else if (is.list(x) && !is.null(x$message)) x$message else
        "Worker exited without a result; inspect worker.log and any stale fold lock."
    }, character(1)))
  write.csv(status, status_path, row.names = FALSE)
  if (any(!good)) stop("One or more folds failed. Completed fold checkpoints are retained. See ",
                       status_path, ". No combined report generated.")
  lapply(results, `[[`, "value")
}

adni_run_folds <- function(folds, FUN, plan, output_root) {
  jobs <- data.frame(fold = folds, output_root = output_root)
  values <- adni_run_jobs(jobs, function(fold) FUN(fold), plan,
                         adni_path(output_root, "fold_worker_status.csv"))
  setNames(values, as.character(folds))
}

ADNI_CASE_SCHEMA <- "case-study-v5-five-fold"

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
  if (any(!tab$installed)) stop("Install case-study dependencies via setup_real_data.R: ",
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
  ggplot2::ggsave(adni_path(root, "main", "fig1_exploration.pdf"), fig, width = 12, height = 3.8)
  ggplot2::ggsave(adni_path(root, "main", "fig1_exploration.png"), fig, width = 12, height = 3.8, dpi = 150)
  counts <- as.data.frame(table(z$pattern, z$CSF))
  names(counts) <- c("pattern", "CSF_availability", "n")
  adni_write_csv(counts, adni_path(root, "appendix", "pattern_counts.csv"))
  desc <- do.call(rbind, lapply(split(z, z$CSF), function(a) data.frame(
    group = as.character(a$CSF[1]), n = nrow(a), mean_age = mean(a$x_age_years),
    female_fraction = mean(a$x_female), mean_apoe4_dose = mean(a$x_apoe4_dose),
    mean_centiloid = mean(a$y_centiloid), sd_centiloid = sd(a$y_centiloid))))
  adni_write_csv(desc, adni_path(root, "appendix", "cohort_by_CSF_availability.csv"))
}

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

adni_prediction_rows <- function(draws, test, method, scenario, cv_folds, fold,
                                 convergence = NA, point_mean = NULL) {
  z <- adni_draw_summary(draws, test$y_centiloid, point_mean)
  data.frame(cv_folds = cv_folds, fold = fold,
    RID = test$RID, R = test$R,
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
  save_fit(list(key = key, value = out), path)
  out
}

adni_extended_fold <- function(fit, train, test, draw_ids, fold_dir, cv_folds,
                                fold, smoke, diagnostic_passed, available_draws,
                                seed_base = NULL) {
  out_dir <- adni_path(fold_dir, "case-study")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  seed <- if (is.null(seed_base)) as.integer(202680000 + cv_folds) else
    as.integer(seed_base)
  cc <- c("c_ab42_ab40_3", "c_gfap_3")
  xc <- c("x_age_years", "x_female", "x_apoe4_dose")
  Cnew <- adni_make_C(test)
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
    zz <- adni_prediction_rows(dd, test, "EIV-GP", "no_test_CSF", cv_folds, fold,
                               diagnostic_passed)
    mm <- adni_metrics(zz); mm$chain <- k; mm
  }))
  if (nrow(chain_scores)) adni_write_csv(chain_scores, adni_path(out_dir, "prediction_scores_by_chain.csv"))
  rows <- list(adni_prediction_rows(no_csf, test, "EIV-GP", "no_test_CSF", cv_folds,
    fold, diagnostic_passed),
    adni_prediction_rows(available_draws, test, "EIV-GP", "available_CSF",
      cv_folds, fold, diagnostic_passed))
  comparison <- adni_competitors(train, test, nrow(no_csf), seed + 1000L, smoke,
                                adni_path(out_dir, "competitors.rds"))
  for (method in names(comparison$draws)) for (scenario in c("no_test_CSF", "available_CSF")) {
    rows[[length(rows) + 1L]] <- adni_prediction_rows(comparison$draws[[method]], test, method,
      scenario, cv_folds, fold, point_mean = comparison$means[[method]])
  }
  adni_write_csv(do.call(rbind, rows), adni_path(out_dir, "predictions.csv"))
  status <- rbind(data.frame(method = "EIV-GP", status = "success",
    optimization_status = if (diagnostic_passed) "mcmc_gate_passed" else "mcmc_gate_failed",
    elapsed_seconds = NA_real_, message = "Sampling time in PROGRESS.md", warnings = ""), comparison$status)
  status$cv_folds <- cv_folds; status$fold <- fold
  adni_write_csv(status, adni_path(out_dir, "method_status.csv"))
  # Save only the bounded predictive draw subset, not another full MCMC fit.
  save_fit(list(no_test_CSF = no_csf, available_CSF = available_draws,
    draw_ids = draw_ids, RID = test$RID), adni_path(out_dir, "eiv_predictive_draws.rds"))

  observed <- which(test$R == 1L)
  if (length(observed)) {
    # Prospective U|C; neither held-out U nor held-out Y is passed to imputation.
    ud <- eivGP::impute_eivgp(fit, new_C = Cnew[observed, , drop = FALSE],
      draw_ids = draw_ids, scale = "raw", seed = seed + 2000L)
    ud <- matrix(ud[, , 1L], nrow = length(draw_ids), ncol = length(observed))
    scores <- adni_draw_summary(ud, test$u_csf_abeta42[observed])
    scores$probability_nonpositive_CSF <- colMeans(ud <= 0)
    urows <- data.frame(cv_folds = cv_folds, fold = fold,
      RID = test$RID[observed],
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
    if (!inherits(ub, "error")) urows <- rbind(urows, data.frame(cv_folds = cv_folds,
      fold = fold, RID = test$RID[observed],
      method = "LM-CSF: response-free", convergence_passed = NA, ub))
    writeLines(if (inherits(ub, "error")) conditionMessage(ub) else "success",
               adni_path(out_dir, "CSF_benchmark_status.txt"))
    adni_write_csv(urows, adni_path(out_dir, "hidden_CSF_validation.csv"))
  }
  # Draws for missing training U are conditional on observed training Y; no truth score exists.
  missing <- which(train$R == 0L)
  ud <- eivGP::impute_eivgp(fit, rows = missing, draw_ids = draw_ids, scale = "raw")
  umat <- matrix(ud[, , 1L], nrow = length(draw_ids), ncol = length(missing))
  adni_write_csv(data.frame(RID = train$RID[missing],
    mean = colMeans(umat), probability_nonpositive_CSF = colMeans(umat <= 0), lo95 = apply(umat, 2, quantile, .025),
    hi95 = apply(umat, 2, quantile, .975)), adni_path(out_dir, "training_missing_CSF_posterior.csv"))

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
      cv_folds = cv_folds, fold = fold,
      convergence_passed = diagnostic_passed),
      adni_path(out_dir, "CSF_PET_surface.csv"))
  }
  # Include the dictionary states deliberately omitted from the lightweight stopping gate.
  adni_write_csv(fit$diagnostics$rhat_hyper, adni_path(out_dir, "all_hyper_diagnostics.csv"))
  j <- fit$mcmc$samples_by_chain$J
  if (length(j)) {
    occupancy <- do.call(rbind, lapply(seq_along(j), function(k) do.call(rbind,
      lapply(seq_len(ncol(j[[k]])), function(d) {
        tb <- table(j[[k]][, d]); data.frame(chain = k, coordinate = d,
          index = as.integer(names(tb)), fraction = as.numeric(tb) / sum(tb))
      }))))
    adni_write_csv(occupancy, adni_path(out_dir, "dictionary_occupancy.csv"))
  }
  invisible(NULL)
}

eivgp_adni_driver <- function(data_dir,output_dir) {
# ADNI implementation extracted from the former analysis notebook.
# Called by adni_run_analysis(); reader-facing ADNI work is deferred.


## Environment variables override the notebook defaults without editing it:
## EIVGP_CORES=12 or EIVGP_SMOKE_TEST=1. Five folds and four chains are fixed.
env_int <- function(name, default, minimum = 0L) {
  raw <- Sys.getenv(name, unset = as.character(default))
  if (!grepl("^[0-9]+$", raw)) stop(name, " must be an integer.")
  x <- suppressWarnings(as.integer(raw))
  if (length(x) != 1L || is.na(x) || x < minimum) stop(name, " is invalid.")
  x
}
env_flag <- function(name, default = FALSE) {
  x <- tolower(Sys.getenv(name, unset = if (default) "true" else "false"))
  if (!x %in% c("true", "false", "1", "0", "yes", "no")) stop(name, " is invalid.")
  x %in% c("true", "1", "yes")
}

PREPARE_ONLY <- isTRUE(get0("ADNI_PREPARE_ONLY", inherits = FALSE))
CV_FOLDS <- 5L
TOTAL_CORES <- env_int("EIVGP_CORES", env_int("EIVGP_N_CORES", 4L, 1L), 1L)
release_budget <- mixedgp_budget_enter(TOTAL_CORES)
on.exit(release_budget(),add=TRUE)
MAX_FOLD_WORKERS <- env_int("ADNI_MAX_FOLD_WORKERS", 5L, 1L)
SMOKE_TEST <- env_flag("EIVGP_SMOKE_TEST", FALSE)


if (SMOKE_TEST) {
  N_CHAINS <- 4L; BURN <- 2L; INITIAL_ITER <- 8L
  EXTENSION_ITER <- 4L; CHECKPOINT_EVERY <- 3L
  FOLDS_TO_RUN <- if (env_flag("ADNI_SMOKE_ALL_FOLDS", FALSE)) 1:5 else 1L
  PREDICTION_DRAWS <- 8L
} else {
  N_CHAINS <- 4L; BURN <- 4000L; INITIAL_ITER <- 20000L
  EXTENSION_ITER <- 10000L; CHECKPOINT_EVERY <- 2000L
  FOLDS_TO_RUN <- 1:5
  PREDICTION_DRAWS <- 1000L
}

## Frozen sampler choice. A matched ADNI pilot found that block sizes 16 and
## 32 reduced wall time per transition but worsened missing-U mixing enough
## that they did not improve convergence per equal waiting time. Keep the
## package reference schedule for the confirmatory run.
SAMPLER_CONTROL <- list(
  u_block_size = 8L,
  cross_every = 10L,
  ess_max_steps = 10000L,
  dictionary_mode = "conditional",
  max_dictionary_size = 125L,
  gp_block_schur = TRUE,
  threshold_update = "gibbs"
)

R_HAT_LIMIT <- 1.01
BULK_ESS_LIMIT <- 400
TAIL_ESS_LIMIT <- 400
MCSE_SD_RATIO_LIMIT <- 0.10

PROJECT_DIR <- system.file(package="eivGP")
DATA_DIR <- data_dir
CORE_PLAN <- if (PREPARE_ONLY) get("ADNI_RUN_CORE_PLAN", inherits = FALSE) else
  adni_core_plan(TOTAL_CORES, length(FOLDS_TO_RUN), N_CHAINS, MAX_FOLD_WORKERS)
N_CORES <- CORE_PLAN$chain_workers
adni_print_core_plan(CORE_PLAN)
DATA_FILE <- adni_path(DATA_DIR, "toledo_adni_cohort_n495.csv")
OUTPUT_ROOT <- adni_path(output_dir, if (SMOKE_TEST) "smoke" else "five_fold")
FOLD_FILE <- adni_path(OUTPUT_ROOT,"appendix","fold_assignments.csv")

dir.create(OUTPUT_ROOT, recursive = TRUE, showWarnings = FALSE)
RUN_LOCK <- adni_path(OUTPUT_ROOT, ".run-lock")
LOCK_OWNED <- FALSE
if (!dir.create(RUN_LOCK, showWarnings = FALSE)) stop("Run already locked: ", RUN_LOCK,
  ". If a job was killed, verify it has stopped before removing this directory.")
LOCK_OWNED <- TRUE
on.exit(if (LOCK_OWNED && dir.exists(RUN_LOCK)) unlink(RUN_LOCK,recursive=TRUE),add=TRUE)
writeLines(paste(Sys.info()[["nodename"]], Sys.getpid()), adni_path(RUN_LOCK, "owner"))


## Dependency installation is a separate, explicit setup step.
EIVGP_TESTED_COMMIT <- paste(eivgp_module_identity(),collapse=";")
if (!requireNamespace("eivGP", quietly = TRUE)) stop("Run setup_real_data.R before analysis.")
pkg <- utils::packageDescription("eivGP")




for (d in c("main", "appendix")) dir.create(adni_path(OUTPUT_ROOT, d), showWarnings = FALSE)
case_packages <- adni_case_preflight(SMOKE_TEST)
adni_write_csv(as.data.frame(CORE_PLAN), adni_path(OUTPUT_ROOT, "appendix", "resource_plan.csv"))
package_files <- list.files(system.file("R", package = "eivGP"), full.names = TRUE)
adni_write_csv(data.frame(file = basename(package_files), md5 = unname(tools::md5sum(package_files))),
  adni_path(OUTPUT_ROOT, "appendix", "eivGP_code_fingerprints.csv"))
adni_write_csv(case_packages, adni_path(OUTPUT_ROOT, "appendix", "package_versions.csv"))
adni_write_csv(data.frame(schema = ADNI_CASE_SCHEMA, smoke = SMOKE_TEST,
  prediction_draws = PREDICTION_DRAWS, run_mode = "five_fold", cv_folds = CV_FOLDS, split_seed = ADNI_CV_SEED,
  primary_target = "Y | test X,C (no test CSF)",
  competitors = paste(ADNI_METHODS, collapse = ";")),
  adni_path(OUTPUT_ROOT, "appendix", "case_config.csv"))

cat("Project:", PROJECT_DIR, "\n")
cat("Data:", DATA_DIR, "\n")
cat("Five-fold CV | folds:", paste(FOLDS_TO_RUN, collapse = ","), "\n")
cat("Mode:", if (SMOKE_TEST) "SMOKE" else "PRODUCTION", "\n")
cat("eivGP:", as.character(utils::packageVersion("eivGP")), "\n")

dat <- adni_read_cohort(DATA_DIR)
folds <- adni_make_folds(dat)
if (file.exists(FOLD_FILE)) {
  previous <- read.csv(FOLD_FILE)
  if (!isTRUE(all.equal(previous, folds, check.attributes = FALSE)))
    stop("Saved fold assignment differs. Use a separate output directory.")
} else adni_write_csv(folds, FOLD_FILE)
adni_write_csv(as.data.frame(table(folds$fold, folds$R)),
               adni_path(OUTPUT_ROOT, "appendix", "fold_CSF_counts.csv"))

X_COLUMNS <- c("x_age_years", "x_female", "x_apoe4_dose")
C_COLUMNS <- c("c_ab42_ab40_3", "c_gfap_3")
Y_COLUMN <- "y_centiloid"
U_COLUMN <- "u_csf_abeta42"


adni_eda(dat, OUTPUT_ROOT, SMOKE_TEST)



atomic_write_lines <- function(lines, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  tmp <- tempfile(paste0(basename(path), "-"), tmpdir = dirname(path))
  on.exit(if (file.exists(tmp)) unlink(tmp), add = TRUE)
  writeLines(lines, tmp)
  if (!file.rename(tmp, path)) stop("Could not publish ", path)
}

stamp <- function(run_dir, message) {
  line <- sprintf("- [%s] %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"), message)
  cat(line, "\n", file = adni_path(run_dir, "PROGRESS.md"), append = TRUE)
  cat(line, "\n")
  flush.console()
}

set_status <- function(run_dir, fold, stage, detail) {
  atomic_write_lines(c(
    paste0("updated_at=", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
    paste0("cv_folds=", CV_FOLDS), paste0("fit_unit=", fold),
    paste0("stage=", stage), paste0("detail=", detail)
  ), adni_path(run_dir, "CURRENT_STATUS.txt"))
  stamp(run_dir, paste("Five-fold CV", "|", stage, "|", detail))
}

light_diagnostics <- function(fit, fold_dir, stage) {
  common <- c("parameter", "rhat", "ess_bulk", "ess_tail", "mcse_sd_ratio")
  parts <- list(
    hyper = subset(fit$diagnostics$rhat_hyper,
                   !grepl("^(J\\[|dictionary_occupancy)", parameter))[, common],
    loading = fit$diagnostics$rhat_A[, common],
    threshold = fit$diagnostics$rhat_tau[, common],
    missing_U = fit$diagnostics$rhat_U[, common]
  )
  tab <- do.call(rbind, lapply(names(parts), function(g) {
    z <- parts[[g]]; z$group <- g; z
  }))
  rownames(tab) <- NULL
  finite <- with(tab, is.finite(rhat) & is.finite(ess_bulk) &
                   is.finite(ess_tail) & is.finite(mcse_sd_ratio))
  passed <- all(finite) && all(tab$rhat <= R_HAT_LIMIT) &&
    all(tab$ess_bulk >= BULK_ESS_LIMIT) && all(tab$ess_tail >= TAIL_ESS_LIMIT) &&
    all(tab$mcse_sd_ratio <= MCSE_SD_RATIO_LIMIT)
  summary <- data.frame(
    stage = stage,
    n_iter = fit$diagnostics$summary$n_iter,
    saved_per_chain = fit$diagnostics$summary$saved_per_chain,
    max_rhat = if (all(finite)) max(tab$rhat) else NA_real_,
    min_bulk_ess = if (all(finite)) min(tab$ess_bulk) else NA_real_,
    min_tail_ess = if (all(finite)) min(tab$ess_tail) else NA_real_,
    max_mcse_sd_ratio = if (all(finite)) max(tab$mcse_sd_ratio) else NA_real_,
    light_gate_passed = passed,
    stringsAsFactors = FALSE
  )
  write.csv(tab, adni_path(fold_dir, paste0("diagnostics_", stage, "_parameters.csv")),
            row.names = FALSE)
  write.csv(summary, adni_path(fold_dir, paste0("diagnostics_", stage, "_summary.csv")),
            row.names = FALSE)
  list(table = tab, summary = summary, passed = passed)
}

trace_series <- function(fit) {
  out <- list(V = fit$mcmc$samples_by_chain$V, r = fit$mcmc$samples_by_chain$r)
  H <- fit$mcmc$samples_by_chain$logtheta
  for (j in seq_len(ncol(H[[1L]]))) {
    nm <- colnames(H[[1L]])[j]
    out[[nm]] <- lapply(H, function(x) x[, j])
  }
  A <- fit$mcmc$samples_by_chain$A
  for (j in seq_len(dim(A[[1L]])[2L]))
    out[[paste0("A[", j, ",1]")]] <- lapply(A, function(x) x[, j, 1L])
  T <- fit$mcmc$samples_by_chain$tau
  for (j in seq_len(ncol(T[[1L]])))
    out[[colnames(T[[1L]])[j]]] <- lapply(T, function(x) x[, j])
  out
}

draw_trace_panel <- function(series, path, title) {
  png(path, width = 2400, height = 1800, res = 180)
  on.exit(dev.off(), add = TRUE)
  old <- par(mfrow = c(4, 4), mar = c(2.4, 2.7, 2.1, 0.7), oma = c(0, 0, 2, 0))
  on.exit(par(old), add = TRUE)
  cols <- c("#0072B2", "#D55E00", "#009E73", "#CC79A7")
  for (nm in names(series)) {
    chains <- series[[nm]]
    n <- length(chains[[1L]])
    ids <- unique(as.integer(round(seq(1, n, length.out = min(n, 2000L)))))
    ylim <- range(unlist(lapply(chains, function(x) x[ids])), finite = TRUE)
    plot(ids, chains[[1L]][ids], type = "l", col = cols[1L], xlab = "retained draw",
         ylab = "", main = nm, ylim = ylim)
    if (length(chains) > 1L) for (cc in 2:length(chains))
      lines(ids, chains[[cc]][ids], col = cols[cc])
  }
  mtext(title, outer = TRUE, cex = 1.2)
}

draw_convergence_plots <- function(fit, diagnostic, fold_dir, stage) {
  draw_trace_panel(trace_series(fit),
                   adni_path(fold_dir, paste0("trace_key_", stage, ".png")),
                   paste("Five-fold CV", stage, "key traces"))

  udiag <- fit$diagnostics$rhat_U
  worst <- head(udiag$global_index[order(udiag$rhat, decreasing = TRUE)], 6L)
  U <- fit$mcmc$samples_by_chain$U
  useries <- setNames(lapply(worst, function(j) lapply(U, function(x) x[, j, 1L])),
                      paste0("U[", worst, "]"))
  draw_trace_panel(useries,
                   adni_path(fold_dir, paste0("trace_worst_U_", stage, ".png")),
                   paste("Worst missing-U traces by raw R-hat:", stage))

  tab <- diagnostic$table
  png(adni_path(fold_dir, paste0("rhat_overview_", stage, ".png")),
      width = 1800, height = 1100, res = 180)
  cols <- c(hyper = "#0072B2", loading = "#D55E00",
            threshold = "#009E73", missing_U = "#777777")
  plot(seq_len(nrow(tab)), tab$rhat, pch = 16, col = cols[tab$group],
       xlab = "raw parameter (grouped)", ylab = "rank-normalized split R-hat",
       main = paste("Five-fold CV", stage),
       ylim = range(c(1, tab$rhat, R_HAT_LIMIT), finite = TRUE))
  abline(h = R_HAT_LIMIT, col = "red", lty = 2, lwd = 2)
  legend("topright", legend = names(cols), col = cols, pch = 16, bty = "n")
  dev.off()
}

balanced_draw_ids <- function(fit, total = 500L) {
  info <- fit$mcmc$mcmc_draw_info
  chains <- sort(unique(info$chain))
  target <- max(1L, floor(total / length(chains)))
  sort(unique(unlist(lapply(chains, function(cc) {
    a <- which(info$chain == cc)
    a[unique(as.integer(round(seq(1, length(a), length.out = min(target, length(a))))))]
  }), use.names = FALSE)))
}

score_predictions <- function(z) {
  split_z <- c(list(overall = z), split(z, paste0("R", z$R)))
  do.call(rbind, lapply(names(split_z), function(g) {
    x <- split_z[[g]]
    data.frame(
      group = g, n = nrow(x),
      RMSE = sqrt(mean((x$y_true - x$pred_mean)^2)),
      MAE = mean(abs(x$y_true - x$pred_mean)),
      Coverage95 = mean(x$y_true >= x$pred_lo95 & x$y_true <= x$pred_hi95),
      Width95 = mean(x$pred_hi95 - x$pred_lo95),
      stringsAsFactors = FALSE
    )
  }))
}

run_one_fold <- function(fold) {
  fold_dir <- adni_path(OUTPUT_ROOT, paste0("fold_", fold))
  dir.create(fold_dir, recursive = TRUE, showWarnings = FALSE)
  unlink(adni_path(fold_dir, "CASE_STUDY_COMPLETED.txt"))
  set_status(fold_dir, fold, "prepare", "Preparing frozen train/test data.")

  test_ids <- folds$RID[folds$fold == fold]
  train <- dat[!(dat$RID %in% test_ids), , drop = FALSE]
  test <- dat[dat$RID %in% test_ids, , drop = FALSE]
  train <- train[order(train$RID), , drop = FALSE]
  test <- test[order(test$RID), , drop = FALSE]
  stopifnot(nrow(train) == 396L, nrow(test) == 99L)

  fit_seed <- 202660000L + fold
  prediction_seed <- 202670000L + fold
  case_seed <- 202680000L + fold

  X_train <- as.matrix(train[, X_COLUMNS, drop = FALSE])
  C_train <- adni_make_C(train)
  U_train <- matrix(train[[U_COLUMN]], ncol = 1L,
                    dimnames = list(NULL, U_COLUMN))
  U_train[train$R == 0L, 1L] <- NA_real_

  config <- data.frame(
    run_mode = "five_fold", cv_folds = CV_FOLDS, fold = fold,
    assignment_source = "fixed_R_stratified_5fold", assignment_seed = ADNI_CV_SEED,
    data_md5 = unname(tools::md5sum(DATA_FILE)),
    folds_md5 = unname(tools::md5sum(FOLD_FILE)),
    n_train = nrow(train), n_test = nrow(test),
    n_train_observed_U = sum(train$R == 1L),
    n_train_missing_U = sum(train$R == 0L),
    engine = "multivariate", latent_dim = 1L, ident = "none",
    kernel = "se", standardize_U = TRUE,
    n_chains = N_CHAINS, burn = BURN, initial_iter = INITIAL_ITER,
    extension_iter = EXTENSION_ITER, checkpoint_every = CHECKPOINT_EVERY,
    u_block_size = SAMPLER_CONTROL$u_block_size,
    cross_every = SAMPLER_CONTROL$cross_every,
    gp_block_schur = SAMPLER_CONTROL$gp_block_schur,
    dictionary_mode = SAMPLER_CONTROL$dictionary_mode,
    max_dictionary_size = SAMPLER_CONTROL$max_dictionary_size,
    ess_max_steps = SAMPLER_CONTROL$ess_max_steps,
    threshold_update = SAMPLER_CONTROL$threshold_update,
    fit_seed = fit_seed, prediction_seed = prediction_seed, case_seed = case_seed,
    eivgp_version = as.character(utils::packageVersion("eivGP")),
    tested_commit = EIVGP_TESTED_COMMIT,
    installed_commit = if (is.null(pkg$RemoteSha)) "unrecorded" else pkg$RemoteSha,
    stringsAsFactors = FALSE
  )
  config_path <- adni_path(fold_dir, "run_config.csv")
  if (file.exists(adni_path(fold_dir, "fit_latest.rds"))) {
    if (!file.exists(config_path)) stop("Checkpoint has no run configuration.")
    old_config <- read.csv(config_path, stringsAsFactors = FALSE)
    if (!isTRUE(all.equal(old_config, config, check.attributes = FALSE)))
      stop("Checkpoint configuration mismatch; use a separate output directory.")
  }
  write.csv(config, config_path, row.names = FALSE)

  fit_path <- adni_path(fold_dir, "fit_latest.rds")
  if (file.exists(fit_path)) {
    fit <- mixedgp_read_artifact(fit_path)
    stopifnot(
      identical(fit$sampler_version, "0.3.1"),
      identical(fit$control$u_block_size, SAMPLER_CONTROL$u_block_size),
      identical(fit$control$cross_every, SAMPLER_CONTROL$cross_every),
      fit$diagnostics$summary$burn == BURN
    )
    stamp(fold_dir, paste("Loaded checkpoint at iteration", fit$diagnostics$summary$n_iter))
  } else {
    first_target <- min(INITIAL_ITER, BURN + CHECKPOINT_EVERY)
    set_status(fold_dir, fold, "fit", paste("Initial fit to iteration", first_target))
    started <- Sys.time()
    fit <- fit_eivgp(
      X = X_train, y = train[[Y_COLUMN]], C = C_train, U_obs = U_train,
      engine = "multivariate", latent_dim = 1L, ident = "none",
      standardize_U = TRUE, store_scores = FALSE, kernel = "se",
      parallel = TRUE, n_cores = N_CORES, n_chains = N_CHAINS,
      n_iter = first_target, burn = BURN, thin = 1L,
      seed = fit_seed,
      verbose = FALSE, sampler_control = SAMPLER_CONTROL
    )
    save_fit(fit, fit_path)
    elapsed <- as.numeric(difftime(Sys.time(), started, units = "secs"))
    stamp(fold_dir, sprintf("Saved checkpoint at %d; segment %.1f min.",
                            first_target, elapsed / 60))
  }

  extend_to <- function(fit, target) {
    while (fit$diagnostics$summary$n_iter < target) {
      current <- fit$diagnostics$summary$n_iter
      add <- min(CHECKPOINT_EVERY, target - current)
      set_status(fold_dir, fold, "continue",
                 sprintf("Iteration %d -> %d; checkpoint will overwrite fit_latest.rds.",
                         current, current + add))
      started <- Sys.time()
      fit <- continue_eivgp(fit, n_iter = add, parallel = TRUE,
                            n_cores = N_CORES, verbose = FALSE)
      save_fit(fit, fit_path)
      elapsed <- as.numeric(difftime(Sys.time(), started, units = "secs"))
      remaining <- target - fit$diagnostics$summary$n_iter
      eta <- if (remaining > 0L) elapsed / add * remaining / 60 else 0
      stamp(fold_dir, sprintf(
        "Saved iteration %d; last %d took %.1f min; ETA to %d about %.1f min.",
        fit$diagnostics$summary$n_iter, add, elapsed / 60, target, eta
      ))
    }
    fit
  }

  fit <- extend_to(fit, INITIAL_ITER)
  current_iter <- fit$diagnostics$summary$n_iter
  final_stage <- paste0("iter", current_iter)
  diagnostic <- light_diagnostics(fit, fold_dir, final_stage)
  draw_convergence_plots(fit, diagnostic, fold_dir, final_stage)
  if (!diagnostic$passed && current_iter < INITIAL_ITER + EXTENSION_ITER) {
    fit <- extend_to(fit, INITIAL_ITER + EXTENSION_ITER)
    final_stage <- paste0("iter", fit$diagnostics$summary$n_iter)
    diagnostic <- light_diagnostics(fit, fold_dir, final_stage)
    draw_convergence_plots(fit, diagnostic, fold_dir, final_stage)
  }
  stamp(fold_dir, paste(final_stage, "convergence gate:", diagnostic$passed))

  set_status(fold_dir, fold, "prediction", paste("Predicting from", final_stage, "fit."))
  X_test <- as.matrix(test[, X_COLUMNS, drop = FALSE])
  C_test <- adni_make_C(test)
  U_test <- matrix(test[[U_COLUMN]], ncol = 1L,
                   dimnames = list(NULL, U_COLUMN))
  U_test[test$R == 0L, 1L] <- NA_real_
  draw_ids <- balanced_draw_ids(fit, PREDICTION_DRAWS)
  response_draws <- predict_eivgp(
    fit, new_X = X_test, new_C = C_test, new_U = U_test,
    new_U_scale = "raw", target = "response", draw_ids = draw_ids,
    seed = prediction_seed
  )
  pred <- data.frame(
    run_mode = "five_fold",
    cv_folds = CV_FOLDS,
    fold = fold, final_stage = final_stage,
    convergence_passed = diagnostic$passed,
    RID = test$RID, diagnosis = test$diagnosis, R = test$R,
    y_true = test[[Y_COLUMN]], pred_mean = colMeans(response_draws),
    pred_sd = apply(response_draws, 2L, stats::sd),
    pred_lo95 = apply(response_draws, 2L, stats::quantile, probs = 0.025),
    pred_hi95 = apply(response_draws, 2L, stats::quantile, probs = 0.975),
    n_prediction_draws = nrow(response_draws), stringsAsFactors = FALSE
  )
  write.csv(pred, adni_path(fold_dir, "eivgp_predictions.csv"), row.names = FALSE)
  metrics <- score_predictions(pred)
  metrics$convergence_passed <- diagnostic$passed
  metrics$run_mode <- "five_fold"
  metrics$cv_folds <- CV_FOLDS
  metrics$fold <- fold; metrics$final_stage <- final_stage
  write.csv(metrics, adni_path(fold_dir, "eivgp_metrics.csv"), row.names = FALSE)
  atomic_write_lines(c(
    paste0("completed_at=", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
    paste0("cv_folds=", CV_FOLDS),
    paste0("fit_unit=", fold),
    paste0("convergence_passed=", diagnostic$passed),
    paste0("final_stage=", final_stage)
  ), adni_path(fold_dir, "EIV_COMPLETED.txt"))
  set_status(fold_dir, fold, "EIV_complete", paste("EIV prediction completed at", final_stage))
  cat("\n### Fold", fold, "\n\n")
  print(knitr::kable(metrics, digits = 3))
  set_status(fold_dir, fold, "comparison", "Running shared prediction and CSF inference tasks.")
  adni_extended_fold(fit, train, test, draw_ids, fold_dir, CV_FOLDS, fold,
                     SMOKE_TEST, diagnostic$passed, response_draws,
                     seed_base = case_seed)
  atomic_write_lines(c(paste0("completed_at=", Sys.time()),
    paste0("schema=", ADNI_CASE_SCHEMA), paste0("convergence_passed=", diagnostic$passed)),
    adni_path(fold_dir, "CASE_STUDY_COMPLETED.txt"))
  set_status(fold_dir, fold, "complete", "EIV fit, comparison, and additional inference finished; inspect method status.")
  rm(fit, response_draws)
  invisible(gc())
  metrics
}

if (!PREPARE_ONLY) {
cat("Running fit units; progress and errors are written to fold_N/worker.log.\n")
all_metrics <- adni_run_folds(FOLDS_TO_RUN, run_one_fold, CORE_PLAN, OUTPUT_ROOT)
combined_metrics <- do.call(rbind, all_metrics)
write.csv(combined_metrics, adni_path(OUTPUT_ROOT, "all_fold_metrics.csv"), row.names = FALSE)
}

if (!PREPARE_ONLY) {
  main_table <- adni_case_report(OUTPUT_ROOT, folds_expected = FOLDS_TO_RUN,
                                 smoke = SMOKE_TEST)

}

writeLines(capture.output(sessionInfo()), adni_path(OUTPUT_ROOT, "sessionInfo.txt"))
sessionInfo()
if (!PREPARE_ONLY) unlink(RUN_LOCK, recursive = TRUE)

}

run_adni_analysis <- function(data_dir,output_dir,action=c("plan","run","report"),
                              profile=c("paper","quick"),cores=4L) {
  action<-match.arg(action);profile<-match.arg(profile)
  folds<-if(profile=="paper")1:5 else 1L
  if(action=="plan") return(data.frame(profile=profile,fold=folds,action="planned",
    data_available=file.exists(adni_path(data_dir,"toledo_adni_cohort_n495.csv"))))
  target<-adni_path(output_dir,if(profile=="paper")"five_fold" else "smoke")
  if(action=="report") return(adni_case_report(target,folds_expected=folds,smoke=profile=="quick"))
  notebook_with_env(c(EIVGP_SMOKE_TEST=if(profile=="quick")"1" else "0",
    EIVGP_CORES=as.character(cores)),eivgp_adni_driver(data_dir,output_dir))
}

adni_make_C <- function(z) {
  out <- data.frame(
    c_ab42_ab40_3 = ordered(z$c_ab42_ab40_3, levels = 1:3),
    c_gfap_3 = ordered(z$c_gfap_3, levels = 1:3),
    check.names = FALSE
  )
  names(out) <- c("c_ab42_ab40_3", "c_gfap_3")
  out
}

adni_case_report <- function(root, folds_expected = 1:5, smoke = FALSE) {
  stopifnot(!anyDuplicated(folds_expected), all(folds_expected %in% 1:5))
  case_dirs <- adni_path(root, paste0("fold_", folds_expected), "case-study")
  files <- adni_path(case_dirs, "predictions.csv")
  if (!all(file.exists(files))) stop("Missing fold prediction files; reporting requires every requested fold.")
  z <- do.call(rbind, lapply(files, read.csv, stringsAsFactors = FALSE))
  status <- do.call(rbind, lapply(adni_path(case_dirs, "method_status.csv"), read.csv,
                                stringsAsFactors = FALSE))
  adni_write_csv(status, adni_path(root, "appendix", "method_status.csv"))
  stopifnot(!anyDuplicated(z[, c("RID", "method", "scenario")]))
  reference <- z[z$method == "EIV-GP" & z$scenario == "no_test_CSF", c("RID", "fold")]
  for (m in ADNI_METHODS[-1]) {
    a <- z[z$method == m & z$scenario == "no_test_CSF", c("RID", "fold")]
    for (f in unique(a$fold)) if (!setequal(a$RID[a$fold == f], reference$RID[reference$fold == f]))
      stop("Participant mismatch for ", m, " fold ", f)
  }
  if (!smoke) {
    assignment_file <- adni_path(root, "appendix", "fold_assignments.csv")
    if (!file.exists(assignment_file)) stop("Missing saved fold assignments.")
    assignment <- read.csv(assignment_file)
    stopifnot(nrow(assignment) == 495L, !anyDuplicated(assignment$RID),
              all(table(factor(assignment$fold, levels = 1:5)) == 99L))
    expected <- assignment[assignment$fold %in% folds_expected, c("RID", "fold")]
    if (nrow(reference) != nrow(expected) || !setequal(reference$RID, expected$RID) ||
        any(reference$fold != expected$fold[match(reference$RID, expected$RID)]))
      stop("Predictions do not match the saved held-out assignments.")
  }
  coverage <- do.call(rbind, lapply(ADNI_METHODS, function(m) {
    a <- z[z$method == m & z$scenario == "no_test_CSF", ]
    complete_unit <- !smoke && nrow(a) == 495L && setequal(unique(a$fold), 1:5)
    data.frame(method = m, n_folds = length(unique(a$fold)), n_predictions = nrow(a),
      all_requested_folds = setequal(unique(a$fold), folds_expected),
      complete_cv = complete_unit,
      flagged_folds = sum(status$method == m & (status$optimization_status != "converged" &
        status$optimization_status != "mcmc_gate_passed"), na.rm = TRUE))
  }))
  adni_write_csv(coverage, adni_path(root, "appendix", "completeness.csv"))
  # Main table requires all methods on the same fold set. Never compare each method's favorable subset.
  common_folds <- Reduce(intersect, lapply(ADNI_METHODS, function(m)
    unique(z$fold[z$method == m & z$scenario == "no_test_CSF"])))
  main <- list(); detailed <- list(); fold_metrics <- list()
  for (m in ADNI_METHODS) for (scenario in unique(z$scenario)) {
    a <- z[z$method == m & z$scenario == scenario, ]
    if (!nrow(a)) next
    metric <- adni_metrics(a); metric$method <- m; metric$scenario <- scenario
    metric$complete_cv <- coverage$complete_cv[coverage$method == m]
    detailed[[length(detailed) + 1L]] <- metric
    for (f in unique(a$fold)) {
      fm <- adni_metrics(a[a$fold == f, ]); fm$method <- m; fm$scenario <- scenario; fm$fold <- f
      fold_metrics[[length(fold_metrics) + 1L]] <- fm
    }
    if (scenario == "no_test_CSF") {
      matched <- a[a$fold %in% common_folds, ]
      if (nrow(matched)) {
        mm <- adni_metrics(matched)
        mm <- mm[mm$group %in% c("overall", "R0", "R1"), ]
        mm$method <- m; mm$matched_folds <- length(common_folds)
        mm$complete_comparison <- !smoke && all(coverage$complete_cv) &&
          setequal(common_folds, 1:5)
        mm$diagnostic_flag <- coverage$flagged_folds[coverage$method == m] > 0
        main[[length(main) + 1L]] <- mm
      }
    }
  }
  adni_write_csv(do.call(rbind, detailed), adni_path(root, "appendix", "all_scenario_metrics.csv"))
  adni_write_csv(do.call(rbind, fold_metrics), adni_path(root, "appendix", "fold_metrics.csv"))
  main_table <- if (length(main)) do.call(rbind, main) else data.frame(
    status = "No fold has valid predictions from all four methods; no comparative ranking available.")
  if ("method" %in% names(main_table)) main_table <- main_table[, c("method", "group", "n",
    "RMSE", "CRPS", "Coverage95", "Width95", "IntervalScore95", "complete_comparison", "diagnostic_flag")]
  adni_write_csv(main_table, adni_path(root, "main", "table1_prediction.csv"))
  writeLines(c("# ADNI prediction comparison", "",
    if (smoke) "SMOKE OUTPUT - no inferential interpretation." else
      if (!all(coverage$complete_cv) || !setequal(common_folds, 1:5))
        "INCOMPLETE: descriptive common-fit subset only; no full-run superiority claim." else
        "All methods evaluated on the identical held-out participants.", "",
    "Primary task: PET response prediction with no test CSF. Overall, R=0, and R=1 scores are reported.",
    "EIV-GP additionally uses calibrated training CSF; competitors use training X,C,Y. This compares available workflows, not an isolated algorithmic effect.",
    "Diagnostic warnings are retained. Coverage and width are not ranked. See appendix for all failures and settings.", "",
    as.character(knitr::kable(main_table, format = "pipe", digits = 3))),
    adni_path(root, "main", "table1_prediction.md"))

  # Per-participant paired contrasts. Bootstrap conditions on fitted folds, not independent refitted studies.
  paired <- list(); set.seed(202650000L)
  for (m in ADNI_METHODS[-1]) for (group in c("overall", "R0")) {
    e <- z[z$method == "EIV-GP" & z$scenario == "no_test_CSF", ]
    a <- z[z$method == m & z$scenario == "no_test_CSF", ]
    if (group == "R0") { e <- e[e$R == 0, ]; a <- a[a$R == 0, ] }
    merge_cols <- c("cv_folds", "fold", "RID")
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
  if (length(paired)) adni_write_csv(do.call(rbind, paired), adni_path(root, "appendix", "paired_losses.csv"))
  # Observable task shared by every method: held-out response predictions.
  zp <- z[z$scenario == "no_test_CSF", ]; zp$method <- factor(zp$method, levels = ADNI_METHODS)
  p <- ggplot2::ggplot(zp, ggplot2::aes(x = y_true, y = pred_mean)) +
    ggplot2::geom_point(alpha = .45, size = .8) + ggplot2::geom_abline(slope = 1, intercept = 0, linetype = 2) +
    ggplot2::facet_wrap(~method, nrow = 1) + ggplot2::theme_bw() +
    ggplot2::labs(x = "Observed Centiloid", y = "Out-of-fold prediction", title = "Common response prediction task")
  ggplot2::ggsave(adni_path(root, "appendix", "prediction_scatter.pdf"), p, width = 11, height = 3.3)
  uf <- adni_path(case_dirs, "hidden_CSF_validation.csv"); uf <- uf[file.exists(uf)]
  if (length(uf)) {
    uz <- do.call(rbind, lapply(uf, read.csv, stringsAsFactors = FALSE))
    adni_write_csv(uz, adni_path(root, "appendix", "hidden_CSF_predictions.csv"))
    um <- do.call(rbind, lapply(split(uz, uz$method), function(a) data.frame(method = a$method[1],
      n = nrow(a), RMSE = sqrt(mean(a$squared_error)), CRPS = mean(a$CRPS),
      Coverage95 = mean(a$covered95), Width95 = mean(a$Width95),
      mean_probability_nonpositive_CSF = mean(a$probability_nonpositive_CSF))))
    adni_write_csv(um, adni_path(root, "appendix", "hidden_CSF_metrics.csv"))
    e <- uz[uz$method == "EIV-GP: prospective U|C", ]
    pi <- ggplot2::ggplot(e, ggplot2::aes(x = y_true, y = pred_mean)) +
      ggplot2::geom_linerange(ggplot2::aes(ymin = pred_lo95, ymax = pred_hi95), alpha = .17) +
      ggplot2::geom_point(size = 1.2, color = "#0072B2") + ggplot2::geom_abline(slope = 1, intercept = 0, linetype = 2) +
      ggplot2::labs(x = "Held-out observed CSF A-beta42", y = "Prospective U|C mean and 95% interval",
                    title = "A  Masked-CSF validation (observed subset)") + ggplot2::theme_bw()
    sf <- adni_path(root, "fold_1", "case-study", "CSF_PET_surface.csv")
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
      ggplot2::ggsave(adni_path(root, "main", "fig2_CSF_inference.pdf"), combined, width = 10, height = 4.5)
      ggplot2::ggsave(adni_path(root, "main", "fig2_CSF_inference.png"), combined, width = 10, height = 4.5, dpi = 150)
    }
  }
  writeLines(c("# Reporting boundaries", "",
    "Do not claim superiority from smoke, incomplete, or inadequately mixed fits.",
    "Positive paired differences favor EIV-GP. Bootstrap intervals are descriptive conditional-on-fits intervals; they omit refitting and model-selection uncertainty.",
    "Each participant is tested once. Training sets overlap, so fold scores are not independent replications.",
    "The literature competitors do not estimate physical CSF posteriors or a CSF-coordinate surface. Other calibrated models can; these capabilities are not universally unique to EIV-GP.",
    "Background and task specifications are in the notebook. EDA is descriptive and must not drive post-hoc test-set tuning."),
    adni_path(root, "appendix", "INTERPRETATION.md"))
  invisible(main_table)
}

# ADNI application: render. Reuses package model and sampler functions.

adni_render_plot <- function(args=commandArgs(TRUE)) {
  stopifnot(length(args)==5L)
  plot<-readRDS(args[1])
  ggplot2::ggsave(args[2],plot,width=as.numeric(args[3]),height=as.numeric(args[4]),
    dpi=as.numeric(args[5]),bg='white')
  stopifnot(file.exists(args[2]),file.info(args[2])$size>1000)
  invisible(args[2])
}

# ADNI application: reporting. Reuses package model and sampler functions.

adni_table_html <- function(d,caption='',note='') {
  if(!is.data.frame(d)||!nrow(d))return('<p>No saved result available for this quantity.</p>')
  for(j in seq_along(d))if(is.numeric(d[[j]]))d[[j]]<-round(d[[j]],4)
  # Reuse Ocean's padding, rules, header shading, escaping and bold-mean markup.
  mixedgp_table_html(d,caption,note)
}

adni_bind <- function(parts) {
  parts<-Filter(function(x)is.data.frame(x)&&nrow(x)>0L,parts)
  if(!length(parts))return(data.frame())
  cols<-unique(unlist(lapply(parts,names)))
  do.call(rbind,lapply(parts,function(x){for(n in setdiff(cols,names(x)))x[[n]]<-NA;x[,cols,drop=FALSE]}))
}

adni_csv <- function(path,required=TRUE) {
  if(!file.exists(path)){if(required)stop('Missing saved report prerequisite: ',path);return(data.frame())}
  read.csv(path,stringsAsFactors=FALSE,check.names=FALSE)
}

adni_summary <- function(d,groups,metrics) {
  if(!nrow(d))return(data.frame())
  split_key<-interaction(d[,groups,drop=FALSE],drop=TRUE,lex.order=TRUE)
  adni_bind(lapply(split(d,split_key),function(x){
    out<-x[1,groups,drop=FALSE];out$folds<-length(unique(x$fold))
    for(k in metrics){v<-as.numeric(x[[k]]);out[[paste0(k,'_mean')]]<-mean(v);out[[paste0(k,'_SD')]]<-if(length(v)>1)sd(v)else NA_real_}
    out
  }))
}

adni_C_dictionary <- function(items,levels,labels=NULL) {
  meanings<-c(MEMORY1='Remembering shopping items without a list',
    MEMORY2='Remembering recent events',MEMORY3='Recalling conversations a few days later',
    MEMORY4='Remembering where objects were placed',MEMORY5='Repeating stories or questions',
    MEMORY6='Remembering the current date or day',MEMORY7='Remembering already having told someone something',
    MEMORY8='Remembering appointments, meetings or engagements',
    DIVATT3='Concentrating without distraction from the environment')
  if(length(labels))for(n in names(labels))meanings[n]<-as.character(labels[[n]])
  text<-unname(meanings[items]);text[is.na(text)]<-'Description not supplied in dataset metadata'
  data.frame(Model_column=paste0('C',seq_along(items)),ECog_field=items,
    Item_meaning=text,Native_grades=vapply(levels,function(n)paste(seq_len(n),collapse='/'),''),check.names=FALSE)
}

adni_saved_source <- function(output_dir,source=NULL) {
  if(!is.null(source)&&nzchar(source)) {
    m<-jsonlite::fromJSON(source,simplifyVector=FALSE)
    if(!identical(m$schema,'adni-saved-report-v1'))stop('Unsupported saved-report source')
    return(m)
  }
  p<-adni_path(output_dir,'application_manifest.json')
  if(!file.exists(p))stop('Saved report requires application_manifest.json or --report-source=FILE')
  m<-jsonlite::fromJSON(p,simplifyVector=FALSE)
  list(schema='adni-saved-report-v1',output=output_dir,input=m$cases[[1]]$input,
    stage=paste0('S',tail(unlist(m$stages),1)),fit_package_version=m$package_version,fit_profile=m$profile,
    fitting_manifest=p,input_sha256=m$frozen_inputs)
}

adni_report_prepare <- function(repo,data_dir,output_dir,settings,report_source=NULL,report_dir=output_dir,group='report') {
  for(p in c('jsonlite','digest','ggplot2','patchwork','posterior'))if(!requireNamespace(p,quietly=TRUE))stop('Report dependency: ',p)
  api<-adni_load_api(repo);api$application$adni_output_check(report_dir,repo)
  m<-adni_saved_source(output_dir,report_source)
  fitting<-if(!is.null(m$fitting_manifest))jsonlite::fromJSON(m$fitting_manifest,simplifyVector=FALSE)else list()
  for(n in setdiff(c('folds','cells'),attr(settings,'explicit_selection')))if(length(fitting[[n]]))
    settings[[n]]<-if(n=='folds')as.integer(unlist(fitting[[n]]))else unlist(fitting[[n]],use.names=FALSE)
  for(n in c('draws_per_chain','baseline_profile'))if(!is.null(fitting[[n]]))settings[[n]]<-fitting[[n]]
  if(length(fitting$samplers))for(cell in names(fitting$samplers))
    settings[[paste0('sampler_',tolower(cell))]]<-fitting$samplers[[cell]]
  if(!is.null(m$fit_profile))settings$profile<-m$fit_profile
  if(!is.null(m$input))data_dir<-m$input
  if(!is.null(m$output))output_dir<-m$output
  input<-api$adni_plan_inputs(data_dir)
  if(length(m$input_sha256))for(p in names(m$input_sha256))
    if(digest::digest(file=p,algo='sha256')!=m$input_sha256[[p]])stop('Frozen report input changed: ',p)
  # Historical fitting identities are read, not rewritten to match reporting code.
  identity<-adni_path(output_dir,'workflow_identity.rds')
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
    dir<-adni_path(output_dir,final,'analysis',sprintf('fold%02d',f),cell)
    if(want_eiv) {
      d<-adni_csv(adni_path(dir,'prediction_metrics.csv'));d<-d[d$protocol%in%c('OI','CC-known'),,drop=FALSE]
      pred[[length(pred)+1L]]<-d
      cv<-adni_csv(adni_path(dir,'convergence.csv'));conv[[length(conv)+1L]]<-cv
      par<-adni_csv(adni_path(dir,'diagnostic_parameters.csv'));par$fold<-f;par$cell<-cell
      parameters[[length(parameters)+1L]]<-par
      paths[[paste(f,cell,sep='/')]]<-list(fold=f,cell=cell,dir=dir,
        fit=adni_path(output_dir,'fits',sprintf('fold%02d',f),cell,'fit.rds'))
      if(!file.exists(paths[[paste(f,cell,sep='/')]]$fit))stop('Missing saved EIV fit; report never fits')
      if(cell!='O0') {
        u<-adni_csv(adni_path(dir,'U_predictions.csv'));u$fold<-f;u$cell<-cell
        up[[length(up)+1L]]<-u;um[[length(um)+1L]]<-adni_csv(adni_path(dir,'U_metrics.csv'))
        p<-adni_csv(adni_path(dir,'PROSPECTIVE_U_predictions.csv'));p$fold<-f;p$cell<-cell
        pp[[length(pp)+1L]]<-p;pm[[length(pm)+1L]]<-adni_csv(adni_path(dir,'PROSPECTIVE_U_metrics.csv'))
      }
    }
  }
  if(want_baseline)for(f in settings$folds)for(cell in c('OI',intersect(settings$cells,c('O20','O50'))))for(method in c('UC-GP','LVGP','EzGP')) {
    dir<-if(cell=='OI')adni_path(output_dir,'baselines',sprintf('fold%02d',f),method)else
      adni_path(output_dir,'CC',cell,'baselines',sprintf('fold%02d',f),method)
    version<-m$fit_package_version
    if(length(m$baseline_overrides))for(b in m$baseline_overrides)
      if(b$fold==f&&b$cell==cell&&b$method==method){dir<-b$directory;version<-b$package_version}
    valid<-all(file.exists(adni_path(dir,c('metrics.csv','fit.rds','prediction_draws.rds'))))
    state<-'MISSING';optimization<-'not_available'
    if(valid) {
      d<-adni_csv(adni_path(dir,'metrics.csv'));d$method<-paste0(method,if(cell=='OI')''else'-CC');d$cell<-cell
      d$stage<-'BASELINE';pred[[length(pred)+1L]]<-d
      a<-readRDS(adni_path(dir,'fit.rds'));optimization<-a$model$optimization_status
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
  metadata_path<-adni_path(data_dir,'DATASET.json')
  meta<-if(file.exists(metadata_path))jsonlite::fromJSON(metadata_path,simplifyVector=FALSE)else list()
  meta<-modifyList(meta,if(is.null(m$metadata))list()else m$metadata)
  if(is.null(meta))meta<-list()
  if(is.null(meta$title))meta$title<-'ADNI memory: partially observed brain imaging'
  if(is.null(meta$population))meta$population<-if(is.null(p$population))'Population specified in private preparation records'else p$population
  if(is.null(meta$calibration_design))meta$calibration_design<-if(is.null(p$calibration_design))
    'C1-proportional stratified sampling, fixed seeded within-stratum order; O20 nested in O50'else
    paste(p$calibration_design,'; strata:',paste(unlist(p$stratification_columns),collapse=', '))
  analysis<-list(repo=repo,api=api,source=m,input=input,output_dir=output_dir,report_dir=report_dir,
    settings=settings,group=group,stage=final,metadata=meta,protocol=p,cohort=cohort,
    prediction=adni_bind(pred),convergence=adni_bind(conv),parameters=adni_bind(parameters),
    U_predictions=adni_bind(up),U_metrics=adni_bind(um),prospective_predictions=adni_bind(pp),
    prospective_metrics=adni_bind(pm),method_status=adni_bind(statuses),paths=paths)
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
    path<-adni_path(output_dir,st,'analysis',sprintf('fold%02d',selected$fold[i]),selected$cell[i],'prediction_metrics.csv')
    d<-adni_csv(path,FALSE)
    if(nrow(d)){d<-d[d$protocol=='OI',];d$stage<-st;stability[[length(stability)+1L]]<-d}
  }
  analysis$stage_scores<-adni_bind(stability)
  adni_report_exports(analysis)
  analysis$plots<-adni_report_plots(analysis)
  if(want_eiv) {
    analysis$functionals<-adni_saved_functionals(analysis)
  }
  analysis$artifact_manifest<-adni_paper_assets(analysis)
  analysis
}

adni_report_exports <- function(a) {
  dest<-adni_path(a$report_dir,'report-assets');dir.create(dest,recursive=TRUE,showWarnings=FALSE)
  for(n in c('prediction','convergence','parameters','U_metrics','prospective_metrics','method_status','stage_scores'))
    if(nrow(a[[n]]))write.csv(a[[n]],adni_path(dest,paste0(n,'.csv')),row.names=FALSE)
  a$api$atomic_json(list(author='Xinye Chen',saved_only=TRUE,new_MCMC=FALSE,new_optimization=FALSE,
    fitting_source=a$source,stage=a$stage,saved_fitting_numerics=a$numerics,report_settings=a$settings,
    input_md5=as.list(tools::md5sum(a$input$files)),analysis_modules=a$api$eivgp_module_identity(),
    reporting_source=as.list(tools::md5sum(c(adni_path(a$repo,'paper/applications/adni.Rmd'),adni_source_files(a$repo))))),
    adni_path(dest,'REPORT_PROVENANCE.json'))
}

adni_metric_table <- function(a,protocol='OI',full=FALSE) {
  d<-a$prediction[a$prediction$protocol==protocol,,drop=FALSE]
  metrics<-c('RMSE','RMSE_SD','MAE_SD','CRPS_SD','Coverage95','Width95_SD','IntervalScore95_SD')
  tab<-adni_summary(d,c('cell','method'),metrics)
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
  adni_table_html(result,if(protocol=='OI')'Ordinal-input prediction: fold mean (sample SD)'else
    'Complete-case prediction: fold mean (sample SD)')
}

adni_u_table <- function(d) {
  if(!nrow(d))return(data.frame())
  tab<-adni_summary(d,c('cell','coordinate'),c('RMSE','RMSE_SD','MAE_SD','CRPS_SD','Coverage95','Width95_SD','IntervalScore95_SD'))
  tab$gain_vs_C1_percent<-vapply(seq_len(nrow(tab)),function(i){
    v<-d[d$cell==tab$cell[i]&d$coordinate==tab$coordinate[i],]
    100*(1-sqrt(mean(v$RMSE^2))/sqrt(mean(v$C1_mean_RMSE^2)))},0)
  out<-tab[,c('cell','coordinate','folds')]
  for(k in c('RMSE','RMSE_SD','MAE_SD','CRPS_SD','Coverage95','Width95_SD','IntervalScore95_SD'))
    out[[k]]<-sprintf('%.4f (%.4f)',tab[[paste0(k,'_mean')]],tab[[paste0(k,'_SD')]])
  out$gain_vs_C1_percent<-sprintf('%+.2f',tab$gain_vs_C1_percent);out
}

adni_target_table <- function(a,all_folds=FALSE) {
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

adni_integration_table <- function(a,all_folds=FALSE) {
  if(is.null(a$functionals))return(data.frame())
  d<-a$functionals$integration
  if(!all_folds)d<-d[d$fold==a$settings$representative_fold,,drop=FALSE]
  data.frame(Fold=d$fold,Calibration=d$cell,Nodes=d$n_latent,`Check nodes`=d$n_latent_check,
    `Mean at nodes`=sprintf('%.5f',d$mean),`Mean at check nodes`=sprintf('%.5f',d$mean_check),
    Difference=sprintf('%+.5f',d$difference),check.names=FALSE)
}

adni_report_plots <- function(a) {
  library(ggplot2);library(patchwork)
  z<-a$input$cells[['1/O50']];d<-a$cohort;p<-a$protocol;plots<-list()
  long<-adni_bind(lapply(seq_along(z$u),function(j)data.frame(U=d[[z$u[j]]],Y=d$Y,C1=factor(d$C1),coordinate=p$U_order[j])))
  plots$cohort<-ggplot(long,aes(U,Y,colour=C1))+geom_point(alpha=.55,size=1.3)+
    facet_wrap(~coordinate,scales='free_x')+labs(x='Observed imaging reference',y='ADNI memory score',colour=p$C_order[1])+
    theme_minimal(base_size=11)+theme(legend.position='bottom')
  xlong<-adni_bind(lapply(seq_along(z$x),function(j)data.frame(X=d[[z$x[j]]],Y=d$Y,C1=factor(d$C1),coordinate=p$X_order[j])))
  plots$covariates<-ggplot(xlong,aes(X,Y,colour=C1))+geom_point(alpha=.55,size=1.3)+
    facet_wrap(~coordinate,scales='free_x')+labs(x='Observed covariate',y='ADNI memory score',colour=p$C_order[1])+
    theme_minimal(base_size=11)+theme(legend.position='bottom')
  shown<-seq_len(min(2L,length(z$c)))
  cu<-adni_bind(lapply(shown,function(k)adni_bind(lapply(seq_along(z$u),function(j)
    data.frame(grade=factor(d[[z$c[k]]],levels=seq_len(p$m_vec[k])),U=d[[z$u[j]]],
      item=paste0('C',k,': ',p$C_order[k]),coordinate=p$U_order[j])))))
  plots$C_U<-ggplot(cu,aes(grade,U))+geom_boxplot(width=.55,outlier.shape=NA,fill='grey90',colour='grey35')+
    geom_point(aes(colour=grade),position=position_jitter(width=.12,height=0,seed=6203L),alpha=.45,size=.9)+
    facet_wrap(~item+coordinate,ncol=length(z$u),scales='free_y',labeller=label_wrap_gen(36))+
    labs(x='Native ECog grade (higher = greater functional difficulty)',y='Measured imaging reference (raw units)')+
    theme_minimal(base_size=11)+theme(legend.position='none')
  counts<-adni_bind(lapply(seq_along(z$c),function(j){t<-as.data.frame(table(factor(d[[z$c[j]]],levels=seq_len(p$m_vec[j]))));
    data.frame(item=p$C_order[j],grade=t$Var1,N=t$Freq)}))
  plots$C_counts<-ggplot(counts,aes(grade,N,fill=grade))+geom_col()+facet_wrap(~item,ncol=3,scales='free_y')+
    labs(x='Native ECog grade',y='Participants')+theme_bw()+theme(legend.position='none')
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
    density<-adni_bind(list(data.frame(hidden[,c('cell','coordinate')],value=hidden$true_std,distribution='Reference'),
      data.frame(hidden[,c('cell','coordinate')],value=hidden$mean_std,distribution='Posterior means')))
    plots$input_distribution<-ggplot(density,aes(value,colour=distribution,fill=distribution))+
      geom_density(alpha=.12)+facet_grid(cell~coordinate,scales='free')+
      labs(x='Training SD units',y='Density',colour=NULL,fill=NULL)+theme_bw()
    class_rows<-adni_bind(lapply(split(hidden,interaction(hidden$fold,hidden$cell,hidden$coordinate,hidden$C1,drop=TRUE)),function(x)
      data.frame(fold=x$fold[1],cell=x$cell[1],coordinate=x$coordinate[1],C1=as.character(x$C1[1]),N=nrow(x),
        RMSE_SD=sqrt(mean((x$mean_std-x$true_std)^2)),Coverage95=mean(x$true_U>=x$q025&x$true_U<=x$q975))))
    plots$class_table<-adni_summary(class_rows,c('cell','coordinate','C1'),c('RMSE_SD','Coverage95'))
    common<-hidden[as.character(hidden$common_hidden)%in%c('TRUE','True','1'),]
    common_rows<-adni_bind(lapply(split(common,interaction(common$fold,common$cell,common$coordinate,drop=TRUE)),function(x)
      data.frame(fold=x$fold[1],cell=x$cell[1],coordinate=x$coordinate[1],N=nrow(x),RMSE_SD=sqrt(mean((x$mean_std-x$true_std)^2)))))
    plots$common_hidden<-adni_summary(common_rows,c('cell','coordinate'),'RMSE_SD')
    plots$support<-adni_bind(lapply(split(hidden,interaction(hidden$cell,hidden$coordinate,drop=TRUE)),function(x)
      data.frame(cell=x$cell[1],coordinate=x$coordinate[1],fold_rows=nrow(x),lower_nonpositive=sum(x$q025<=0),
        mean_nonpositive=sum(x$posterior_mean<=0),interval_outside_C_grade='Not applicable: C is not an imaging bin')))
  }
  if(nrow(u))plots<-modifyList(plots,adni_input_distributions(a))
  plots
}

adni_input_distributions <- function(a) {
  library(ggplot2);parts<-list()
  for(key in names(a$paths)) {
    path<-a$paths[[key]];if(path$cell=='O0'||path$fold!=a$settings$representative_fold)next
    z<-a$input$cells[[key]]
    hidden<-which(a$input$cells[[paste(path$fold,'O50',sep='/')]]$train$calibration==0L)
    saved<-readRDS(adni_path(path$dir,'U_draws.rds'))
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
  d<-adni_bind(parts);if(!nrow(d))return(list())
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

adni_paper_assets <- function(a) {
  dest<-adni_path(a$report_dir,'report-assets','paper');dir.create(dest,recursive=TRUE,showWarnings=FALSE)
  plots<-list(fig_adni_cohort=a$plots$cohort,fig_adni_C_support=a$plots$C_counts,
    fig_adni_Y_X=a$plots$covariates,fig_adni_C_U=a$plots$C_U,
    fig_adni_target_traces=a$functionals$trace_plot,fig_adni_mean_functions=a$functionals$mean_plot,
    fig_adni_response_surfaces=a$functionals$surface_plot,fig_adni_U_imputation=a$plots$imputation,
    fig_adni_U_distributions=a$plots$input_distribution,fig_adni_C1_U_distributions=a$plots$class_distribution)
  manifest<-list()
  for(name in names(plots)) {
    file<-paste0(name,'.pdf');available<-!is.null(plots[[name]])
    if(available)ggplot2::ggsave(adni_path(dest,file),plots[[name]],width=12,height=if(grepl('support|traces|C1',name))9 else 6,
      device=grDevices::pdf,limitsize=FALSE)
    manifest[[length(manifest)+1L]]<-data.frame(label=name,file=file,status=if(available)'available'else'not_applicable',
      N=a$input$N,profile=a$settings$profile,MD5=if(available)unname(tools::md5sum(adni_path(dest,file)))else'')
  }
  oi<-a$prediction[a$prediction$protocol=='OI',,drop=FALSE]
  cc<-a$prediction[a$prediction$protocol=='CC-known',,drop=FALSE]
  metrics<-c('RMSE','RMSE_SD','MAE_SD','CRPS_SD','Coverage95','Width95_SD','IntervalScore95_SD')
  tables<-list(tab_adni_design=a$input$audit,
    tab_adni_C_definitions=adni_C_dictionary(a$protocol$C_order,a$protocol$m_vec,a$metadata$C_labels),
    tab_adni_prediction_OI=adni_summary(oi,c('cell','method'),metrics),
    tab_adni_prediction_CC=adni_summary(cc,c('cell','method'),metrics),tab_adni_training_U=adni_u_table(a$U_metrics),
    tab_adni_prospective_U=adni_u_table(a$prospective_metrics),tab_adni_convergence=a$convergence)
  for(name in names(tables)) {
    d<-tables[[name]];available<-nrow(d)>0L;file<-paste0(name,'.csv')
    if(available) {
      write.csv(d,adni_path(dest,file),row.names=FALSE)
      writeLines(adni_table_html(d,name),adni_path(dest,paste0(name,'.html')))
    }
    manifest[[length(manifest)+1L]]<-data.frame(label=name,file=file,status=if(available)'available'else'unavailable',
      N=a$input$N,profile=a$settings$profile,MD5=if(available)unname(tools::md5sum(adni_path(dest,file)))else'')
  }
  out<-adni_bind(manifest);write.csv(out,adni_path(dest,'artifact_manifest.csv'),row.names=FALSE);out
}

# ADNI application: settings. Reuses package model and sampler functions.

adni_resolve_path <- function(path,repo) {
  path<-path.expand(path)
  if(!grepl('^(/|[A-Za-z]:[/\\\\])',path))path<-adni_path(repo,path)
  base<-path;suffix<-character()
  while(!file.exists(base)&&dirname(base)!=base) {
    suffix<-c(basename(base),suffix);base<-dirname(base)
  }
  resolved<-normalizePath(base,mustWork=TRUE)
  if(length(suffix))adni_path(resolved,paste(suffix,collapse='/'))else resolved
}

adni_cli_params <- function(args) {
  allowed<-c('fit_jobs','evaluation_jobs','cores','profile','dataset','data_root','data_dir','output','output_dir','config','manifest',
    'chains','iter','burn','stages','folds','cells','draws_per_chain','first_checkpoint','checkpoint_step',
    'baseline_profile','sampler_o0','sampler_o20','sampler_o50','representative_fold','appendix',
    'plan','reuse_only','report_source','report_dir','functional_draws_per_chain','n_latent','run_id','ess_min')
  p<-list()
  for(a in args) {
    if(!grepl('^--[^=]+=.+$',a))stop('Use --name=value: ',a)
    k<-gsub('-','_',sub('^--([^=]+)=.*$','\\1',a));v<-sub('^[^=]+=','',a)
    if(!k%in%allowed)stop('Unknown ADNI option: ',k)
    if(k%in%names(p))stop('Repeated option: ',k)
    if(k%in%c('plan','reuse_only','appendix')) {
      if(!v%in%c('true','false'))stop(k,' must be true or false')
      v<-v=='true'
    }
    p[[k]]<-v
  }
  p
}

adni_integer <- function(value,name,minimum=1L) {
  n<-suppressWarnings(as.numeric(value))
  if(length(n)!=1L||!is.finite(n)||n<minimum||n!=floor(n)||n>.Machine$integer.max)
    stop(name,' must be an integer >= ',minimum)
  as.integer(n)
}

adni_integer_list <- function(value,name) {
  if(is.character(value)&&length(value)==1L) {
    if(grepl('^[0-9]+:[0-9]+$',value)) {
      z<-as.integer(strsplit(value,':',fixed=TRUE)[[1]]);value<-seq.int(z[1],z[2])
    }else value<-strsplit(value,',',fixed=TRUE)[[1]]
  }
  unname(vapply(unlist(value),adni_integer,1L,name=name))
}

adni_settings <- function(profile='pilot',cores=4L,...) {
  overrides<-list(...);profile<-match.arg(profile,c('pilot','paper','quick'))
  defaults<-switch(profile,paper=list(chains=4L,iter=50000L,burn=5000L,stages=50000L),
    pilot=list(chains=4L,iter=1500L,burn=500L,stages=c(1000L,1500L)),
    quick=list(chains=1L,iter=80L,burn=30L,stages=80L))
  s<-modifyList(c(defaults,list(profile=profile,cores=cores,folds=if(profile=='quick')1L else 1:5,
    cells=c('O50','O20','O0'),draws_per_chain=if(profile=='quick')10L else 250L,
    checkpoint_step=1000L,baseline_profile=if(profile=='quick')'quick'else'paper',
    sampler_o0='HMC',sampler_o20='HMC',sampler_o50='ESS',representative_fold=1L,
    appendix=TRUE,functional_draws_per_chain=64L,n_latent=256L,run_id='local')),overrides)
  s$ess_min<-mixedgp_parse_ess_min(s$ess_min)
  if(!is.null(s$ess_min) && as.numeric(s$chains)<4L)stop('Scientific ESS requires at least four chains.')
  for(n in c('cores','chains','iter','draws_per_chain','checkpoint_step','representative_fold',
             'functional_draws_per_chain','n_latent'))s[[n]]<-adni_integer(s[[n]],n)
  s$burn<-adni_integer(s$burn,'burn',0L);s$folds<-adni_integer_list(s$folds,'folds')
  if(any(!s$folds%in%1:5)||anyDuplicated(s$folds))stop('folds must be distinct members of1:5')
  s$cells<-if(length(s$cells)==1L)strsplit(s$cells,',',fixed=TRUE)[[1]]else unlist(s$cells)
  if(!length(s$cells)||any(!s$cells%in%c('O0','O20','O50'))||anyDuplicated(s$cells))stop('Invalid cells')
  if(!is.null(overrides$stages)) {
    s$stages<-adni_integer_list(overrides$stages,'stages')
    if(!is.null(overrides$iter)&&tail(s$stages,1)!=s$iter)stop('iter must equal the last stage')
    s$iter<-tail(s$stages,1)
  }else if(!is.null(overrides$iter))s$stages<-s$iter
  if(any(diff(s$stages)<=0)||min(s$stages)<=s$burn)stop('Stages must increase and exceed burn')
  s$first_checkpoint<-adni_integer(if(is.null(s$first_checkpoint))min(s$stages[1],s$burn+250L)else s$first_checkpoint,'first_checkpoint')
  if(s$first_checkpoint<=s$burn||s$first_checkpoint>min(s$stages))stop('first-checkpoint must exceed burn and not exceed first stage')
  if(!s$baseline_profile%in%c('paper','quick'))stop('baseline-profile must be paper or quick')
  for(n in c('sampler_o0','sampler_o20','sampler_o50'))if(!s[[n]]%in%c('HMC','ESS'))stop('Unknown sampler')
  attr(s,'explicit_selection')<-intersect(names(overrides),c('folds','cells'))
  s
}

adni_method_run <- function(repo,params,group) {
  if(!is.null(params$config)) {
    config<-jsonlite::fromJSON(adni_resolve_path(params$config,repo),simplifyVector=TRUE)
    params<-modifyList(config,params)
  }
  get<-function(n,default)if(is.null(params[[n]]))default else params[[n]]
  fit_jobs<-if(is.null(params$fit_jobs))NULL else adni_integer(params$fit_jobs,'fit-jobs')
  evaluation_jobs<-if(is.null(params$evaluation_jobs))NULL else adni_integer(params$evaluation_jobs,'evaluation-jobs')
  if(!is.null(params$manifest)) {
    if(!is.null(params$ess_min)) {
      declared<-jsonlite::fromJSON(adni_resolve_path(params$manifest,repo))$ess_min
      if(!identical(mixedgp_parse_ess_min(params$ess_min),mixedgp_parse_ess_min(declared)))
        stop('--ess-min must match the supplied manifest; set ess_min in that manifest.')
    }
    env<-adni_load_api(repo)
    return(env$adni_screen_run(repo,list(manifest=adni_resolve_path(params$manifest,repo),
      cores=adni_integer(get('cores',4L),'cores'),fit_jobs=fit_jobs,evaluation_jobs=evaluation_jobs,plan=isTRUE(params$plan),reuse_only=isTRUE(params$reuse_only)),group))
  }
  fields<-c('chains','iter','burn','stages','folds','cells','draws_per_chain','first_checkpoint','checkpoint_step',
    'baseline_profile','sampler_o0','sampler_o20','sampler_o50','representative_fold','appendix',
    'functional_draws_per_chain','n_latent','run_id','ess_min')
  s<-do.call(adni_settings,c(list(profile=get('profile','pilot'),cores=get('cores',4L)),params[intersect(names(params),fields)]))
  dataset<-get('dataset','adni')
  if(!grepl('^[A-Za-z0-9][A-Za-z0-9_-]*$',dataset))stop('Invalid dataset name')
  default_input<-adni_default_input(dataset,repo)
  input<-adni_resolve_path(get('data_root',get('data_dir',default_input)),repo)
  out<-adni_resolve_path(get('output',get('output_dir',adni_default_output(dataset,s$profile,repo))),repo)
  if(isTRUE(params$plan)&&(!dir.exists(input)||!all(file.exists(adni_input_files(adni_plan_input_root(input),1L,'O0'))))) {
    cat('Dataset:',dataset,'; data:',input,'\n');print(s[c('profile','chains','iter','burn','cores')])
    cat('Supply private inputs before fitting. No computation started.\n');return(invisible(s))
  }
  report_source<-if(!is.null(params$report_source))adni_resolve_path(params$report_source,repo)else NULL
  saved_only<-group=='report'||isTRUE(params$reuse_only)
  if(saved_only&&!isTRUE(params$plan)) {
    default_report<-if(group=='report')out else adni_path(out,'reports',group)
    return(adni_render_report(repo,input,out,s,report_source,
      report_dir=adni_resolve_path(get('report_dir',default_report),repo),group=group,
      filename=if(group!='report'&&adni_storage_flat(out))'index.html'else'adni.html'))
  }
  api<-adni_load_api(repo);inp<-api$adni_plan_inputs(input);p<-inp$cells[[1]]$protocol
  api$application$adni_output_check(out,repo)
  if(!isTRUE(params$plan)) {
    adni_storage_init(out,s$stages)
    if(adni_storage_flat(out)) {
      input<-adni_snapshot_inputs(input,out,inp$files)
      inp<-api$adni_plan_inputs(input)
    }
  }
  files<-setNames(lapply(inp$files,function(x)digest::digest(file=x,algo='sha256')),inp$files)
  m<-list(schema='adni-application-v1',author='Xinye Chen',profile=s$profile,
    package_version=as.character(packageVersion('eivGP')),
    baseline_profile=s$baseline_profile,trial_id=dataset,smoke=s$profile=='quick',
    output_root=adni_workflow_root(out),cores=s$cores,chains=s$chains,stages=s$stages,burn=s$burn,
    first_checkpoint=s$first_checkpoint,checkpoint_step=s$checkpoint_step,draws_per_chain=s$draws_per_chain,ess_min=s$ess_min,
    samplers=list(O0=s$sampler_o0,O20=s$sampler_o20,O50=s$sampler_o50),folds=s$folds,cells=s$cells,
    cases=list(list(id=p$plan_id,input=input,output=out,C=p$C_order,N=inp$N,q=p$q)),frozen_inputs=files)
  if(isTRUE(params$plan)) {
    e<-adni_load_api(repo)
    tasks<-e$adni_screen_tasks(m,group)
    print(tasks,row.names=FALSE);print(s[c('profile','chains','iter','burn','cores')])
    cat('Resident fit jobs:',adni_resident_workers(s$cores,'fit',nrow(tasks),s$chains,fit_jobs),
      '; resident evaluation jobs:',adni_resident_workers(s$cores,'evaluate',nrow(tasks),s$chains,evaluation_jobs),'\n')
    return(invisible(m))
  }
  dir.create(out,recursive=TRUE,showWarnings=FALSE);manifest<-adni_path(out,'application_manifest.json')
  if(file.exists(manifest)) {
    prior<-jsonlite::fromJSON(manifest,simplifyVector=FALSE);a<-m;b<-prior;a$cores<-b$cores<-NULL
    canonical<-function(x)jsonlite::fromJSON(jsonlite::toJSON(x,auto_unbox=TRUE,null='null'),simplifyVector=FALSE)
    if(!identical(canonical(a),canonical(b)))
      stop('Saved ADNI settings/inputs differ; select a new output directory. ',
        paste(all.equal(canonical(a),canonical(b)),collapse='; '))
  }else api$atomic_json(m,manifest)
  e<-adni_load_api(repo)
  result<-e$adni_screen_run(repo,list(manifest=manifest,cores=s$cores,fit_jobs=fit_jobs,evaluation_jobs=evaluation_jobs,plan=FALSE,reuse_only=FALSE),group)
  adni_render_report(repo,input,out,s,report_dir=adni_path(out,'reports',group),group=group,
    filename=if(adni_storage_flat(out))'index.html'else'adni.html')
  invisible(result)
}

adni_render_report <- function(repo,data_dir,output_dir,settings,report_source=NULL,
                              report_dir=output_dir,group='report',rmd=NULL,filename='adni.html') {
  if(!requireNamespace('rmarkdown',quietly=TRUE)||!rmarkdown::pandoc_available())
    stop('HTML requires rmarkdown and Pandoc on PATH or RSTUDIO_PANDOC')
  env<-adni_load_api(repo)
  analysis<-env$adni_report_prepare(repo,data_dir,output_dir,settings,report_source,report_dir,group)
  render_env<-new.env(parent=globalenv());render_env$.adni_precomputed_analysis<-analysis
  render_env$.adni_repository_dir<-repo
  if(is.null(rmd))rmd<-adni_path(repo,'paper/applications/adni.Rmd')
  rmarkdown::render(rmd,output_file=filename,output_dir=report_dir,intermediates_dir=report_dir,
    params=list(profile=settings$profile,action='report',cores=settings$cores,data_dir=data_dir,
      output_dir=output_dir,report_source=if(is.null(report_source))''else report_source,
      representative_fold=settings$representative_fold,appendix=settings$appendix),envir=render_env,quiet=TRUE)
  cat('Saved ADNI HTML: ',adni_path(report_dir,filename),'\n',sep='')
  invisible(list(html=adni_path(report_dir,filename),analysis=analysis))
}

# Author: Xinye Chen. Base-R C/X-only folds and nested stratified observation masks.
adni_seed <- function(seed,label) {
  h<-digest::digest(paste(seed,label,sep='/'),algo='sha256',serialize=FALSE)
  as.integer(strtoi(substr(h,1,7),16L)%%(2^30-1)+1)
}
adni_nested_masks <- function(d,strata='C1',seed,fractions=c(.2,.5)) {
  key<-do.call(paste,c(d[,strata,drop=FALSE],sep='/'));groups<-split(seq_len(nrow(d)),key)
  ranked<-lapply(names(groups),function(g){set.seed(adni_seed(seed,g));v<-groups[[g]];v[sample.int(length(v))]})
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
adni_balanced_folds <- function(d,C,X,seed) {
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
adni_prepare_folds <- function(cohort_file,metadata_file,output_dir) {
  if(!requireNamespace('jsonlite',quietly=TRUE)||!requireNamespace('digest',quietly=TRUE))stop('Install jsonlite and digest')
  if(dir.exists(adni_path(output_dir,'MC'))||
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
  }else d$fold<-adni_balanced_folds(d,C,X,meta$fold_seed)
  for(f in 1:5)for(j in seq_along(C))if(!all(seq_len(meta$m_vec[j])%in%d[d$fold!=f,C[j]]))
    stop('Every declared native grade must be represented in each training fold; do not silently merge levels or reroll')
  for(f in 1:5) {
    tr<-d[d$fold!=f,fields];te<-d[d$fold==f,fields];tr$calibration<-te$calibration<-0L
    cal<-adni_nested_masks(tr,strata,adni_seed(meta$calibration_seed,paste0('fold',f)))
    mask<-adni_nested_masks(te,strata,adni_seed(meta$test_seed,paste0('fold',f)))
    for(cell in c('O20','O50'))if(qr(cbind(1,as.matrix(tr[tr$RID%in%cal[[cell]],U,drop=FALSE])))$rank<length(U)+1L)
      stop('The frozen calibration vectors have insufficient affine rank; no outcome-guided mask reroll is performed')
    for(cell in c('O0','O20','O50')) {
      dir.create(output_dir,recursive=TRUE,showWarnings=FALSE)
      paths<-base::file.path(output_dir,paste0(sprintf('fold%02d_%s_',f,cell),c('PROTOCOL.json','train.csv','test.csv')))
      train<-tr;train$calibration<-as.integer(train$RID%in%cal[[cell]])
      write.csv(train,paths[2],row.names=FALSE);write.csv(te,paths[3],row.names=FALSE)
      seed_names<-c('production','posterior','prediction','mixed_prediction','prospective_U','UC_GP','LVGP','EzGP','EzGP_draws')
      seeds<-setNames(lapply(seed_names,function(n)adni_seed(meta$model_seed,paste(f,if(n%in%c('UC_GP','LVGP','EzGP','EzGP_draws'))'baseline'else cell,n,sep='/'))),seed_names)
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
  jsonlite::write_json(meta,adni_path(output_dir,'DATASET.json'),auto_unbox=TRUE,pretty=TRUE)
  invisible(d[,c('RID','fold')])
}

# Author: Xinye Chen. ADNI CLI orchestration; application implementation is in this file.
adni_method_main <- function(group,args=commandArgs(TRUE)) {
  script<-sub('^--file=','',grep('^--file=',commandArgs(FALSE),value=TRUE)[1])
  repo<-normalizePath(adni_path(dirname(normalizePath(script)),'..'))
  api<-adni_load_api(repo)
  if(identical(args,'--help')) {
    cat('Rscript tools/adni-',group,'.R --profile=paper --cores=4\n',
      '  --profile=pilot|paper|quick; --dataset=NAME; --data-root=PATH; --output=PATH\n',
      '  --fit-jobs=N --evaluation-jobs=N (resident-job caps; CPU budget unchanged)\n',
      '  --config=FILE; --iter=N --burn=N --chains=N --stages=N,N\n',
      '  --folds=1:5 --cells=O0,O20,O50; --plan=true; --reuse-only=true; --report-source=FILE\n',
      '  Run EIV-GP and competitors sequentially; report reads saved results.\n',sep='')
    return(invisible(NULL))
  }
  Sys.setenv(OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1',VECLIB_MAXIMUM_THREADS='1')
  api$adni_method_run(repo,api$adni_cli_params(args),group)
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
  e$analyze_cell<-e$adni_analyze_cell;e$baseline_cell<-e$adni_baseline_cell
  e$adni_load_api<-e$adni_screen_api<-e$adni_workflow_api<-function(repo=NULL)e
  e$run_adni_application<-function(repository_dir=repo,params=list(),group=c('report','eivGP','competitors'))
    e$adni_method_run(repository_dir,params,match.arg(group))
  e
}
adni_workflow_api <- adni_screen_api <- adni_load_api
# Export application functions when sourced, retaining the fresh-environment loader.
local({
  api<-adni_load_api()
  bindings<-as.list(api)
  bindings[c('adni_load_api','adni_screen_api','adni_workflow_api')]<-NULL
  list2env(bindings,envir=parent.env(environment()))
})

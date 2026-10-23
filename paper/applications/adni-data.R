# Author: Xinye Chen. ADNI measurement, validation and scoring adapters.
# Models, prediction, optimizers and CPU scheduling remain in the unmodified installed upstream package.
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
    if (!inherits(model, "EzGP model"))
        stop("model must be an EzGP model.")
    if (!requireNamespace("EzGP", quietly = TRUE))
        stop("Package 'EzGP' is not installed.")
    X_new <- as.matrix(X_new)
    X_train <- as.matrix(model$data$X)
    Y <- as.numeric(model$data$Y)
    p <- as.integer(model$data$p)
    q <- as.integer(model$data$q)
    m <- as.integer(model$data$m)
    n <- nrow(X_train)
    tau <- as.numeric(model$data$tau)
    parv <- as.numeric(model$param)
    if (ncol(X_new) != ncol(X_train))
        stop("X_new has the wrong number of columns.")
    if (length(parv) == 0L || any(!is.finite(parv)) || !is.finite(tau) || tau < 0) {
        stop("EzGP fit has invalid parameters or nugget.")
    }
    K <- EzGP:::cov_m(X_train, p, q, m, n, parv, tau)
    Tm <- try(chol(K), silent = TRUE)
    if (inherits(Tm, "try-error"))
        stop("EzGP training covariance is not Cholesky-factorable.")
    invT <- backsolve(Tm, diag(n))
    invc <- invT %*% t(invT)
    one <- matrix(1, nrow = n, ncol = 1L)
    y <- matrix(Y, ncol = 1L)
    mu <- as.numeric((t(one) %*% invc %*% y)/sum(invc))
    npar <- 1L + q + p + p * sum(m)
    psum <- function(x1, x2, par2) sum(-par2 * (x1 - x2)^2)
    covx <- function(w1, w2) {
        par1 <- parv[seq_len(q + 1L)]
        par2 <- parv[(q + 2L):(q + 1L + p)]
        par3 <- parv[(q + 2L + p):npar]
        x1 <- w1[seq_len(p)]
        z1 <- w1[p + seq_len(q)]
        x2 <- w2[seq_len(p)]
        z2 <- w2[p + seq_len(q)]
        value <- par1[1L] * exp(psum(x1, x2, par2))
        for (j in seq_len(q)) {
            if (z1[j] == z2[j]) {
                level <- as.integer(z1[j])
                offset <- (sum(m[seq_len(j)]) - m[j] + level - 1L) * p
                idx <- offset + seq_len(p)
                value <- value + par1[j + 1L] * exp(psum(x1, x2, par3[idx]))
            }
        }
        value
    }
    means <- latent_var <- numeric(nrow(X_new))
    prior_var <- sum(parv[seq_len(q + 1L)])
    centered_y <- y - mu * one
    for (j in seq_len(nrow(X_new))) {
        gamma <- matrix(vapply(seq_len(n), function(i) covx(X_new[j, ], X_train[i, ]), numeric(1)), ncol = 1L)
        means[j] <- as.numeric(mu + t(gamma) %*% invc %*% centered_y)
        latent_var[j] <- as.numeric(prior_var - t(gamma) %*% invc %*% gamma)
    }
    predictive_var <- latent_var + tau
    tolerance <- 100 * (nrow(X_new) + 1L) * .Machine$double.eps * max(1, max(abs(predictive_var)))
    if (any(!is.finite(means)) || any(!is.finite(predictive_var)) || any(predictive_var < -tolerance)) {
        stop("Covariance-consistent EzGP prediction returned invalid predictive variances.")
    }
    predictive_var <- pmax(predictive_var, 0)
    list(Y_hat = means, latent_MSE = latent_var, predictive_MSE = predictive_var, tau = tau)
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

adni_plan_input_root <- function(path) {
    path <- path.expand(path)
    if (dir.exists(file.path(path, "MC")))
        path <- file.path(path, "MC")
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
            dir <- file.path(root, sprintf("fold%02d", f), cell)
            paths <- file.path(dir, c("PROTOCOL.json", "private/train.csv", "private/test.csv"))
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

adni_plan_fit_path <- function(out, task) file.path(out, "fits", sprintf("fold%02d", task$fold), task$cell, "fit.rds")

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
    target <- s$iter
    request <- file.path(dirname(path), "continuation_request.rds")
    if (extend) {
        if (is.null(fit))
            stop("Missing initial fit: ", path)
        tab <- adni_plan_diagnostics(fit)
        if (!file.exists(request)) {
            if (isTRUE(tab$summary$core_pass))
                return(invisible("SKIPPED_CONVERGED"))
            ocean_atomic_rds(list(from = fit$checkpoint$iteration, target = fit$checkpoint$iteration + s$continuation_iter,
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
    dg <- adni_plan_diagnostics(fit)
    ocean_atomic_csv(dg$parameters, file.path(dirname(path), "diagnostic_parameters.csv"))
    ocean_atomic_csv(dg$summary, file.path(dirname(path), "diagnostic_summary.csv"))
    invisible("DONE")
}

adni_plan_diagnostics <- function(fit) {
    parts <- fit$diagnostics[grep("^rhat_", names(fit$diagnostics))]
    parts <- Filter(function(z) is.data.frame(z) && all(c("parameter", "rhat", "ess_bulk", "ess_tail") %in% names(z)),
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
    dir <- file.path(out, "baselines", sprintf("fold%02d", task$fold), task$method)
    path <- file.path(dir, "fit.rds")
    if (all(file.exists(file.path(dir, c("fit.rds", "prediction_draws.rds")))))
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
        controls <- if (s$profile == "paper") ocean_controls()[[task$method]] else
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
        ocean_atomic_csv(sel$audit, file.path(dir, "tau_cv_audit.csv"))
        ocean_atomic_csv(sel$candidates, file.path(dir, "tau_cv_candidates.csv"))
    }
    ocean_atomic_rds(list(model = model, center = mu, scale = ss, y_center = yc, y_scale = ys), path)
    mm <- cbind(data.frame(fold = task$fold, cell = "OI", method = task$method, protocol = "OI"), adni_plan_score(draws,
        point, te$Y, ys))
    ocean_atomic_csv(mm, file.path(dir, "metrics.csv"))
    ocean_atomic_rds(list(draws = draws, point = point), file.path(dir, "prediction_draws.rds"))
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
    dir <- file.path(out, "analysis", sprintf("fold%02d", task$fold), task$cell)
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
    ocean_atomic_csv(predictions, file.path(dir, "Y_predictions.csv"))
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
        ocean_atomic_csv(data.frame(row = seq_len(nrow(te)), observed_Y = te$Y, predicted_Y = mp, U_available = z$available),
            file.path(dir, "mixed_Y_predictions.csv"))
    }
    ocean_atomic_csv(metrics, file.path(dir, "prediction_metrics.csv"))
    dg <- adni_plan_diagnostics(fit)
    ocean_atomic_csv(cbind(fold = task$fold, cell = task$cell, dg$summary), file.path(dir, "convergence.csv"))
    ocean_atomic_csv(dg$parameters, file.path(dir, "diagnostic_parameters.csv"))
    scale <- if (any(tr$calibration == 1L))
        "raw"
    else "model"
    imp <- impute_eivgp(fit, rows = seq_len(nrow(tr)), draw_ids = ids, scale = scale)
    ocean_atomic_rds(list(draws = imp, scale = scale, draw_ids = ids, coordinates = p$U_order), file.path(dir,
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
        ocean_atomic_csv(do.call(rbind, tables), file.path(dir, "U_metrics.csv"))
        ocean_atomic_csv(do.call(rbind, plots), file.path(dir, "U_predictions.csv"))
    }
    if (task$fold == s$folds[1]) {
        series <- mixedgp_v030_raw_series(fit)
        parameters <- head(grep("^(V$|r$|theta_x|tau|A\\[)", names(series), value = TRUE), 6L)
        if (length(parameters) && length(fit$checkpoint$states) >= 2L)
            ggplot2::ggsave(file.path(dir, "mcmc_trace.png"), plot_mcmc(fit, parameters), width = 10, height = 7,
                dpi = 130)
        else if (length(fit$checkpoint$states) < 2L)
            writeLines("Single-chain workflow check: multi-chain trace diagnostics are unavailable; convergence is not established.",
                file.path(dir, "MCMC_TRACE_UNAVAILABLE.txt"))
    }
    invisible("DONE")
}

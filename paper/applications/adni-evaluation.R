# Author: Xinye Chen. Preserve approved held-out and latent-recovery protocols as application code.
pilot_ggsave <- function(filename, plot, width, height, dpi = 130, ...) {
    input <- paste0(filename, ".plot.rds")
    saveRDS(plot, input)
    args <- c(shQuote(file.path(repository_dir, "paper/applications/adni-render.R")), shQuote(input), shQuote(filename), width, height, dpi)
    rc <- system2(file.path(R.home("bin"), "Rscript"), args, stdout = paste0(filename, ".render.log"), stderr = paste0(filename,
        ".render.log"))
    stopifnot(rc == 0L, file.exists(filename), file.info(filename)$size > 1000)
    unlink(input)
    invisible(filename)
}

replace_plot_device <- function(x) {
    if (is.call(x) && identical(x, quote(ggplot2::ggsave)))
        return(as.name("pilot_ggsave"))
    if (is.call(x))
        for (i in seq_along(x)) x[[i]] <- replace_plot_device(x[[i]])
    x
}

reference_points <- function(calibration, targets, c_names, u_name) {
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

analyze_cell <- function(task, inputs, out, s, stage) {
    z <- inputs$cells[[paste(task$fold, task$cell, sep = "/")]]
    tr <- z$train
    te <- z$test
    p <- z$protocol
    stageout <- file.path(out, stage)
    dir <- file.path(stageout, "analysis", sprintf("fold%02d", task$fold), task$cell)
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
            atomic_rds(list(draws = r, draw_ids = a$draw_ids), file.path(dir, "Y_draws.rds"))
        }
        else {
            cache$exact <- r
            atomic_rds(list(draws = r, draw_ids = a$draw_ids), file.path(dir, "known_U_Y_draws.rds"))
        }
        r
    }
    env$pilot_ggsave <- pilot_ggsave
    fun <- api$adni_plan_analyze_cell
    body(fun) <- replace_plot_device(body(fun))
    environment(fun) <- env
    fun(task, inputs, stageout, s)
    fit <- cache$fit
    stopifnot(fit$data$p == length(z$x), fit$data$q == p$q, fit$data$d == length(z$u), all(is.na(fit$data$U_obs_raw[tr$calibration ==
        0, , drop = FALSE])))
    met <- read.csv(file.path(dir, "prediction_metrics.csv"))
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
    atomic_csv(met, file.path(dir, "prediction_metrics.csv"))
    dg <- read.csv(file.path(dir, "convergence.csv"))
    dg$stage <- stage
    dg$sampler <- fit$control$sampler
    if (!is.null(fit$diagnostics$hmc))
        atomic_csv(fit$diagnostics$hmc, file.path(dir, "HMC_DIAGNOSTICS.csv"))
    atomic_csv(dg, file.path(dir, "convergence.csv"))
    cm <- attr(cache$pure, "conditional_means", exact = TRUE)
    per <- nrow(cm)%/%s$chains
    target_diagnostics <- do.call(rbind, lapply(seq_len(ncol(cm)), function(j) {
        a <- matrix(cm[, j], per, s$chains)
        data.frame(RID = te$RID[j], predicted_Y = mean(a), selected_Rhat = posterior::rhat(a), selected_ESS = posterior::ess_bulk(a),
            selected_MCSE = posterior::mcse_mean(a))
    }))
    atomic_csv(target_diagnostics, file.path(dir, "PREDICTION_FUNCTION_DIAGNOSTICS.csv"))
    if (task$cell != "O0") {
        cal <- tr[tr$calibration == 1, , drop = FALSE]
        hidden <- which(tr$calibration == 0)
        ui <- read.csv(file.path(dir, "U_predictions.csv"))
        um <- read.csv(file.path(dir, "U_metrics.csv"))
        for (j in seq_along(z$u)) {
            truth <- tr[[z$u[j]]]
            ref <- reference_points(cal, tr[hidden, , drop = FALSE], z$c, z$u[j])
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
        atomic_csv(um, file.path(dir, "U_metrics.csv"))
        atomic_csv(ui, file.path(dir, "U_predictions.csv"))
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
            ref <- reference_points(cal, te, z$c, z$u[j])
            qq <- apply(ud, 2, quantile, c(0.025, 0.975))
            pm[[j]] <- cbind(data.frame(stage = stage, fold = task$fold, cell = task$cell, coordinate = p$U_order[j]),
                api$adni_plan_score(ud, point, truth, sd(tr[[z$u[j]]])), calibration_mean_RMSE = sqrt(mean((ref$global -
                  truth)^2)), C1_mean_RMSE = sqrt(mean((ref$C1 - truth)^2)), joint_C_mean_RMSE = sqrt(mean((ref$joint -
                  truth)^2)), joint_C_fallback_N = ref$fallback)
            pp[[j]] <- data.frame(RID = te$RID, coordinate = p$U_order[j], true_U = truth, posterior_mean = point,
                q025 = qq[1, ], q975 = qq[2, ])
        }
        atomic_csv(do.call(rbind, pm), file.path(dir, "PROSPECTIVE_U_metrics.csv"))
        atomic_csv(do.call(rbind, pp), file.path(dir, "PROSPECTIVE_U_predictions.csv"))
        ui$role <- factor(ifelse(ui$calibration == 1, "Observed calibration", "Masked training U"), levels = c("Masked training U",
            "Observed calibration"))
        plot <- ggplot(ui, aes(true_std, mean_std, colour = factor(C1))) + geom_abline(slope = 1, intercept = 0,
            linetype = 2, colour = "grey50") + geom_linerange(data = ui[ui$calibration == 0, ], aes(ymin = lo_std,
            ymax = hi_std), alpha = 0.15, linewidth = 0.3) + geom_point(aes(shape = role), size = 1.2, alpha = 0.65) +
            scale_shape_manual(values = c(16, 1)) + facet_wrap(~coordinate) + theme_bw(base_size = 11) + theme(legend.position = "bottom") +
            labs(x = "Reference U (training SD units; evaluation only)", y = "Posterior mean (training SD units)",
                colour = "MEMORY3", shape = NULL, title = paste(p$plan_id, "fold", task$fold, task$cell, stage),
                subtitle = paste(fit$control$sampler, "; CORE pass:", dg$core_pass))
        pilot_ggsave(file.path(dir, "training_U_imputation_std.png"), plot, width = 10, height = 5.3, dpi = 140,
            bg = "white")
    }
    else writeLines("O0 is unanchored: no physical-U RMSE or truth-based alignment.", file.path(dir, "O0_SCALE_LIMITATION.txt"))
    atomic_json(list(state = "PASS", stage = stage, iterations = s$iter, primary_no_test_U_Y = TRUE, hidden_training_U_masked = TRUE,
        calibration_excluded_from_U_scoring = TRUE, prospective_U_uses_C_not_test_Y = TRUE, sampler = fit$control$sampler),
        file.path(dir, "EVALUATION_VALIDATION.json"))
    rm(fit, cache)
    gc(FALSE)
    invisible("DONE")
}

baseline_cell <- function(task, inputs, out, s, cc = FALSE) {
    source <- inputs$cells[[paste(task$fold, if (cc)
        task$cell
    else "O50", sep = "/")]]
    tr <- source$train
    te <- source$test
    root <- if (cc)
        file.path(out, "CC", task$cell)
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
    dir <- file.path(root, "baselines", sprintf("fold%02d", task$fold), task$method)
    saved <- readRDS(file.path(dir, "prediction_draws.rds"))
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
    atomic_csv(mm, file.path(dir, "metrics.csv"))
    atomic_csv(data.frame(RID = te$RID, observed_Y = te$Y, predicted_Y = saved$point, scored = seq_len(nrow(te)) %in%
        selected), file.path(dir, "Y_predictions.csv"))
    invisible("DONE")
}

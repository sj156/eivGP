if (!requireNamespace("ggplot2", quietly = TRUE)) {
  stop("Install ggplot2 to draw the frozen-data figures.")
}
library(ggplot2)
VIS_REP <- 1L
VIS_CALIBRATION <- 50L
manifest_paths <- sort(list.files(DATA_ROOT, pattern = "^manifest\\.rds$",
                                 recursive = TRUE, full.names = TRUE))
frozen_examples <- list()
inventory_rows <- list()
for (manifest_path in manifest_paths) {
  manifest <- readRDS(manifest_path)
  if (!is.data.frame(manifest) ||
      !all(c("study", "rep", "file", "md5") %in% names(manifest))) {
    stop("Malformed manifest: ", manifest_path)
  }
  selected <- manifest[manifest$study == "study2" & manifest$rep == VIS_REP, , drop = FALSE]
  if (!nrow(selected)) next
  if (nrow(selected) != 1L) stop("Ambiguous replication in: ", manifest_path)
  data_path <- file.path(dirname(manifest_path), basename(selected$file))
  if (!file.exists(data_path) ||
      !identical(unname(tools::md5sum(data_path)), as.character(selected$md5))) {
    stop("Missing or checksum-mismatched frozen data: ", data_path)
  }
  frozen <- readRDS(data_path)
  stopifnot(identical(frozen$study, "study2"), frozen$rep_id == VIS_REP,
            is.list(frozen$data), is.list(frozen$calibration_sets))
  setting <- basename(dirname(manifest_path))
  if (setting %in% names(frozen_examples)) stop("Duplicate setting folder: ", setting)
  frozen_examples[[setting]] <- frozen
  inventory_rows[[setting]] <- data.frame(
    setting = setting, replication = frozen$rep_id,
    n_train = nrow(frozen$data$train$X),
    n_test = nrow(frozen$data$test$X),
    saved_calibration_sizes = paste(names(frozen$calibration_sets), collapse = ", "),
    data_seed = frozen$seeds$data, calibration_seed = frozen$seeds$calibration,
    file = basename(data_path), md5 = as.character(selected$md5)
  )
}
if (!length(frozen_examples)) {
  message("No frozen Study II replication-1 datasets found under ", DATA_ROOT,
          ". No data were generated and no plots were fabricated.")
} else {
  FIGURE_ROOT <- if (exists("params") && is.character(params$figure_root) &&
                    length(params$figure_root) == 1L && nzchar(params$figure_root)) {
    params$figure_root
  } else Sys.getenv("EIVGP_FIGURE_ROOT", unset = file.path(ARTIFACT_ROOT, "figures"))
  figure_dir <- file.path(FIGURE_ROOT, "frozen-data", "study2")
  dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
  inventory <- do.call(rbind, inventory_rows)
  rownames(inventory) <- NULL
  print(knitr::kable(inventory[, 1:7], caption = "Frozen datasets actually visualized"))
  write.csv(inventory, file.path(figure_dir, "visualized_dataset_inventory.csv"), row.names = FALSE)
}
save_frozen_plot <- function(plot, name, width = 10, height = 7) {
  print(plot)
  ggplot2::ggsave(file.path(figure_dir, paste0(name, ".png")), plot,
                 width = width, height = height, dpi = 160, bg = "white")
  ggplot2::ggsave(file.path(figure_dir, paste0(name, ".pdf")), plot,
                 width = width, height = height, bg = "white")
}
frozen_theme <- theme_bw(base_size = 11) +
  theme(panel.grid.minor = element_blank(), legend.position = "bottom")

if (length(frozen_examples)) {
  category_rows <- latent_rows <- response_rows <- list()
  for (setting in names(frozen_examples)) {
    frozen <- frozen_examples[[setting]]
    train <- frozen$data$train
    stopifnot(is.matrix(train$X), is.matrix(train$U), is.matrix(train$C),
              ncol(train$U) >= 2L, nrow(train$X) == length(train$y),
              all(is.finite(train$X)), all(is.finite(train$U)),
              all(is.finite(train$y)), all(is.finite(train$f)))
    for (j in seq_len(ncol(train$C))) {
      category_rows[[paste(setting, j)]] <- data.frame(
        setting = setting, proxy = paste0("C", j),
        category = factor(seq_len(frozen$design$m)),
        proportion = tabulate(train$C[, j], nbins = frozen$design$m) / nrow(train$C))
    }
    key <- as.character(VIS_CALIBRATION)
    available <- key %in% names(frozen$calibration_sets)
    ids <- if (available) frozen$calibration_sets[[key]] else integer()
    if (available) stopifnot(length(ids) == VIS_CALIBRATION, !anyDuplicated(ids),
                             all(ids %in% seq_len(nrow(train$X))))
    label <- paste0(setting, if (available) " | calibration=50" else " | calibration=50 unavailable")
    latent_rows[[setting]] <- data.frame(
      U1 = train$U[, 1], U2 = train$U[, 2], f = train$f,
      setting = label, calibrated = seq_len(nrow(train$X)) %in% ids)
    for (j in seq_len(ncol(train$X))) {
      response_rows[[paste(setting, j)]] <- data.frame(
        x = train$X[, j], y = train$y, U1 = train$U[, 1],
        input = paste0("X", j), setting = setting)
    }
  }
  category_data <- do.call(rbind, category_rows)
  latent_data <- do.call(rbind, latent_rows)
  response_data <- do.call(rbind, response_rows)
  p_categories <- ggplot(category_data, aes(category, proxy, fill = proportion)) +
    geom_tile(colour = "white") +
    geom_text(aes(label = sprintf("%.2f", proportion)), size = 3, colour = "black") +
    scale_fill_gradient(low = "#F0F9E8", high = "#7BCCC4", limits = c(0, 1)) +
    facet_wrap(~setting, ncol = 2, scales = "free_y") +
    labs(title = "Study II: frequencies of ordinal proxies",
         subtitle = "Frozen replication 1; training observations",
         x = "Ordinal category", y = "Proxy", fill = "Proportion") + frozen_theme
  save_frozen_plot(p_categories, "proxy_frequencies", 11, 8)
  p_latent <- ggplot(latent_data, aes(U1, U2)) +
    geom_point(aes(colour = f), size = 2) +
    geom_point(data = subset(latent_data, calibrated), shape = 1,
               size = 3.3, colour = "black", stroke = 0.7) +
    scale_colour_viridis_c(option = "C") + facet_wrap(~setting, ncol = 2) +
    labs(title = "Study II: latent coordinates and saved response",
         subtitle = "Black rings mark the saved calibration subset when available",
         x = "Latent U1", y = "Latent U2", colour = "True f(X,U)",
         caption = "X also varies: this scatterplot is not a fixed-X response-surface slice.") +
    frozen_theme
  save_frozen_plot(p_latent, "latent_response_and_calibration", 11, 8)
  p_response <- ggplot(response_data, aes(x, y, colour = U1)) +
    geom_point(size = 1.7, alpha = 0.8) + scale_colour_viridis_c() +
    facet_grid(setting ~ input) +
    labs(title = "Study II: responses versus observed inputs",
         subtitle = "Latent U1 is shown for explanation, not supplied as uncalibrated training information",
         x = "Observed input", y = "Response y", colour = "Latent U1") + frozen_theme
  save_frozen_plot(p_response, "response_by_input", 11, max(5, 2 * length(frozen_examples)))
}

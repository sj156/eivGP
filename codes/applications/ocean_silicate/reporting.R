# Combine the original and continued paper fits, retaining every baseline.
ocean_read_paper_results <- function(output_dir) {
  bind <- function(parts) {
    columns <- unique(unlist(lapply(parts, names)))
    do.call(rbind, lapply(parts, function(x) {
      for (name in setdiff(columns, names(x))) x[[name]] <- NA
      x[, columns, drop = FALSE]
    }))
  }
  metrics <- predictions <- diagnostics <- list()
  for (fold in 1:5) {
    initial <- file.path(output_dir, "full5_class4", sprintf("fold%02d", fold))
    metrics[[paste(fold, "baselines")]] <- read.csv(file.path(initial, "baseline_metrics.csv"))
    for (cell in c("O0", "O20", "O50")) {
      continued <- file.path(output_dir, "continuation_plus20k", sprintf("fold%02d", fold), cell)
      path <- if (file.exists(file.path(continued, "predictive_metrics.csv"))) continued else file.path(initial, cell)
      key <- paste(fold, cell)
      metrics[[key]] <- read.csv(file.path(path, "predictive_metrics.csv"))
      z <- read.csv(file.path(path, "predictions.csv")); z$fold <- fold; z$cell <- cell
      predictions[[key]] <- z
      z <- read.csv(file.path(path, "convergence_summary.csv")); z$fold <- fold; z$cell <- cell
      diagnostics[[key]] <- z
    }
  }
  list(metrics = bind(metrics), predictions = bind(predictions), diagnostics = bind(diagnostics))
}

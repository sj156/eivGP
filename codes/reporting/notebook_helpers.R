# Keep HTML focused on analysis while preserving complete execution logs.
notebook_logged <- function(path, expr) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  connection <- file(path, open = "wt")
  output_depth <- sink.number()
  sink(connection)
  on.exit({
    while (sink.number() > output_depth) sink()
    close(connection)
  }, add = TRUE)
  warnings <- character()
  value <- withCallingHandlers(eval(substitute(expr), envir = parent.frame()),
    message = function(m) { cat(conditionMessage(m)); invokeRestart("muffleMessage") },
    warning = function(w) {
      warnings <<- unique(c(warnings, conditionMessage(w)))
      cat("WARNING:", conditionMessage(w), "\n")
      invokeRestart("muffleWarning")
    })
  attr(value, "notebook_warnings") <- warnings
  value
}

notebook_with_env <- function(values, expr) {
  old <- Sys.getenv(names(values), unset = NA_character_)
  on.exit({
    if (any(!is.na(old))) do.call(Sys.setenv, as.list(old[!is.na(old)]))
    if (any(is.na(old))) Sys.unsetenv(names(old)[is.na(old)])
  }, add = TRUE)
  do.call(Sys.setenv, as.list(values))
  eval(substitute(expr), envir = parent.frame())
}

notebook_cache_identity <- function(directory, identity) {
  path <- file.path(directory, "notebook_identity.rds")
  if (file.exists(path)) {
    if (!identical(readRDS(path), identity)) stop("Inputs or settings changed. Use a new output directory: ", directory)
  } else {
    if (dir.exists(directory) && length(list.files(directory, all.files = TRUE, no.. = TRUE)))
      stop("Existing outputs lack notebook provenance. Choose a fresh output directory: ", directory)
    dir.create(directory, recursive = TRUE, showWarnings = FALSE)
    saveRDS(identity, path)
  }
  invisible(identity)
}

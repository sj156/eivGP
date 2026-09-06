############################################################
## 08_published_competitor_validation.R
##
## Interface and small-data validation for the public-package competitor
## wrappers used in Studies I and II. Installed packages are exercised through
## their real fitting and prediction functions; missing packages are reported
## explicitly and are never replaced by an internal approximation.
############################################################

args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
script_dir <- if (length(file_arg) > 0L) {
  dirname(normalizePath(sub("^--file=", "", file_arg[1L])))
} else {
  getwd()
}

source(file.path(script_dir, "00_project_setup.R"), chdir = TRUE)
source(file.path(script_dir, "03_study2_published_competitors.R"))
source(file.path(script_dir, "00_public_api.R"))

read_bool_env <- function(name, default = FALSE) {
  value <- tolower(Sys.getenv(name, unset = if (default) "true" else "false"))
  value %in% c("1", "true", "yes", "y")
}

set.seed(20260830L)
n_train <- 32L
n_test <- 7L
m_vec <- rep(4L, 4L)

X_train <- cbind(
  x1 = stats::rnorm(n_train),
  x2 = stats::runif(n_train, -1, 1)
)
C_train <- cbind(
  c1 = rep(1:4, length.out = n_train),
  c2 = rep(1:4, each = 2L, length.out = n_train),
  c3 = rep(c(1L, 3L, 2L, 4L), each = 4L, length.out = n_train),
  c4 = rep(4:1, each = 8L, length.out = n_train)
)
y_factor_effect <- matrix(
  c(
    -0.50, 0.10, 0.45, 0.80,
    -0.20, 0.25, 0.40, -0.10,
    -0.15, 0.05, 0.20, 0.35,
    0.30, 0.10, -0.10, -0.25
  ),
  nrow = 4L,
  byrow = TRUE
)
y_train <-
  sin(X_train[, 1L]) +
  0.4 * X_train[, 2L] +
  rowSums(vapply(seq_len(ncol(C_train)), function(j) {
    y_factor_effect[j, C_train[, j]]
  }, numeric(n_train))) +
  stats::rnorm(n_train, sd = 0.08)
X_test <- cbind(
  x1 = seq(-0.8, 0.8, length.out = n_test),
  x2 = seq(0.7, -0.7, length.out = n_test)
)
C_test <- cbind(
  c1 = rep(1:4, length.out = n_test),
  c2 = rep(4:1, length.out = n_test),
  c3 = rep(c(1L, 3L, 2L, 4L), length.out = n_test),
  c4 = rep(c(2L, 4L, 1L, 3L), length.out = n_test)
)

## Data construction must preserve quantitative dimensions and the declared
## factor levels, including levels absent from a particular test set.
checked <- validate_mixedgp_inputs(X_train, C_train, m_vec)
stopifnot(
  identical(dim(checked$X), c(n_train, 2L)),
  identical(dim(checked$C), c(n_train, 4L)),
  identical(checked$m_vec, m_vec)
)
factor_data <- make_mixedgp_competitor_data(X_test, C_test, m_vec)
stopifnot(
  identical(
    names(factor_data),
    c("x1", "x2", "c1", "c2", "c3", "c4")
  ),
  identical(levels(factor_data$c1), as.character(seq_len(m_vec[1L]))),
  identical(levels(factor_data$c2), as.character(seq_len(m_vec[2L])))
)
numeric_data <- make_mixedgp_numeric_matrix(X_test, C_test, m_vec)
stopifnot(
  identical(dim(numeric_data), c(n_test, 6L)),
  all(is.finite(numeric_data))
)

## All wrappers target Y* | X*, C*: the public interface accepts no test U.
wrapper_formals <- names(formals(run_published_mixedgp_competitors))
stopifnot(!any(grepl("(^|_)U($|_)", wrapper_formals)))
stopifnot(all(vapply(
  c("fit_ucgp", "fit_lvgp", "fit_ezgp", "fit_mixedgp_competitor"),
  function(name) is.function(get(name, mode = "function")),
  logical(1)
)))

draw_check <- sample_independent_predictive_marginals(
  mean = seq_len(n_test),
  variance = rep(0.25, n_test),
  n_draw = 11L,
  seed = 1L
)
stopifnot(identical(dim(draw_check), c(11L, n_test)))

bad_levels <- C_test
bad_levels[1L, 1L] <- m_vec[1L] + 1L
bad_level_error <- try(
  validate_mixedgp_inputs(X_test, bad_levels, m_vec),
  silent = TRUE
)
stopifnot(inherits(bad_level_error, "try-error"))

preflight <- mixedgp_competitor_preflight(strict = FALSE)
fit_check <- run_published_mixedgp_competitors(
  X_train = X_train,
  y_train = y_train,
  C_train = C_train,
  X_test = X_test,
  C_test = C_test,
  n_draw = 11L,
  seed = 20260831L,
  m_vec = m_vec,
  strict = FALSE,
  controls = list(
    `UC-GP` = list(n_starts = 8L),
    LVGP = list(
      ## Eight starts are part of the frozen publication specification and
      ## materially improve convergence for the four-factor fit.
      n_starts = 8L,
      max_iter_ini = 100L,
      max_iter_lat = 20L
    ),
    EzGP = list(
      tau_fractions = 0.01,
      cv_folds = 2L,
      maxeval = 30L
    )
  )
)

stopifnot(
  identical(sort(fit_check$status$method), sort(preflight$method)),
  all(fit_check$status$status %in% c("success", "unavailable_or_failed"))
)
for (method in names(fit_check$draws)) {
  component_means <- attr(
    fit_check$draws[[method]], "conditional_means", exact = TRUE
  )
  component_vars <- attr(
    fit_check$draws[[method]], "conditional_vars", exact = TRUE
  )
  stopifnot(
    identical(dim(fit_check$draws[[method]]), c(11L, n_test)),
    all(is.finite(fit_check$draws[[method]])),
    identical(dim(component_means), c(1L, n_test)),
    identical(dim(component_vars), c(1L, n_test)),
    all(component_vars > 0),
    isTRUE(all.equal(
      as.numeric(component_vars),
      as.numeric(fit_check$predictive_variances[[method]])
    ))
  )
}

preflight_for_merge <- preflight
names(preflight_for_merge)[
  names(preflight_for_merge) == "implementation"
] <- "registered_implementation"
status_for_merge <- fit_check$status
names(status_for_merge)[
  names(status_for_merge) == "implementation"
] <- "executed_implementation"
validation <- merge(
  preflight_for_merge,
  status_for_merge,
  by = "method",
  all.x = TRUE,
  sort = FALSE
)
validation <- validation[match(preflight$method, validation$method), ]
out_dir <- normalizePath(file.path(script_dir, "..", "tables"), mustWork = TRUE)
utils::write.csv(
  validation,
  file.path(out_dir, "published_competitor_validation.csv"),
  row.names = FALSE
)

if (read_bool_env("MIXEDGP_REQUIRE_ALL_COMPETITORS", FALSE) &&
    any(validation$status != "success")) {
  failed <- validation$method[validation$status != "success"]
  stop(
    "Strict competitor validation failed for: ",
    paste(failed, collapse = ", "),
    ". See published_competitor_validation.csv."
  )
}

print(validation[, c(
  "method", "package", "available", "version", "status", "message",
  "warnings"
)])
cat("Published-competitor interface checks completed.\n")

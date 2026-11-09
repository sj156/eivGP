# Plot subprocesses must not reload a newer application during an older run.
repo <- Sys.getenv('EIVGP_TEST_REPO'); stopifnot(nzchar(repo))
source(file.path(repo, 'tools/oai.R'))
source(file.path(repo, 'tools/adni.R'))

for (application in c('oai', 'adni')) local({
  app <- application
  api <- get(paste0(app, '_load_api'))(repo)
  save_plot <- get(paste0(app, '_save_plot'), api)
  testthat::test_that(paste(app, 'plot rendering is independent of application bootstrap'), {
    root <- tempfile("plot 'with spaces "); dir.create(root)
    on.exit(unlink(root, recursive = TRUE))
    # This inherited library selection caused the production version mismatch.
    old <- Sys.getenv('EIVGP_R_LIBRARY', unset = NA_character_)
    Sys.setenv(EIVGP_R_LIBRARY = file.path(root, 'older-library'))
    on.exit(if (is.na(old)) Sys.unsetenv('EIVGP_R_LIBRARY') else
      Sys.setenv(EIVGP_R_LIBRARY = old), add = TRUE)
    api$repository_dir <- file.path(root, 'nonexistent-checkout')
    api[[paste0(app, '_load_api')]] <- function(...) stop('Application must not be bootstrapped')
    plot <- ggplot2::ggplot(data.frame(x = 1:4, y = 1:4), ggplot2::aes(x, y)) + ggplot2::geom_point()
    path <- file.path(root, 'rendered plot.png')
    testthat::expect_identical(save_plot(path, plot, width = 4, height = 3), path)
    testthat::expect_gt(file.info(path)$size, 1000)
    testthat::expect_false(file.exists(paste0(path, '.plot.rds')))
    testthat::expect_true(file.exists(paste0(path, '.render.log')))
  })
  testthat::test_that(paste(app, 'rendering failures show the subprocess log and retain the plot'), {
    root <- tempfile(); dir.create(root)
    on.exit(unlink(root, recursive = TRUE))
    path <- file.path(root, 'failed.png')
    testthat::expect_error(save_plot(path, structure(list(), class = 'invalid_plot'), 4, 3),
      'Plot rendering failed.*Log:.*render[.]log')
    testthat::expect_true(file.exists(paste0(path, '.plot.rds')))
    testthat::expect_match(paste(readLines(paste0(path, '.render.log')), collapse = '\n'), 'Error')
  })
})

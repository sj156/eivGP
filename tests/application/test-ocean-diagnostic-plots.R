repo <- Sys.getenv('EIVGP_TEST_REPO'); stopifnot(nzchar(repo))
source(file.path(repo, 'tools/ocean-tools.R'), local = TRUE)

ocean_plot_fixture <- function() {
  u <- expand.grid(C = 1:4, cell = c('O0', 'O20', 'O50'), fold = 1:5)
  u$observed <- u$cell != 'O0' & u$C == 1L
  u$plot_reference <- u$C / 2
  u$plot_estimate <- u$plot_reference + .1
  u$plot_q025 <- u$plot_estimate - .3
  u$plot_q975 <- u$plot_estimate + .3
  m <- expand.grid(C = 1:4, cell = c('O0', 'O20', 'O50'), fold = 1:5,
    pressure = c(600, 650, 700), design = 'candidate_A', stringsAsFactors = FALSE)
  same <- m[m$pressure == 650 & m$cell != 'O0', ]
  same$design <- 'same_x'
  m <- rbind(m, same)
  m$posterior_mean <- m$C + m$pressure / 1000
  m$q025 <- m$posterior_mean - .2
  m$q975 <- m$posterior_mean + .2
  f <- expand.grid(U = seq(25, 175, by = 25), cell = c('O20', 'O50'), fold = 1:5)
  f$posterior_mean <- f$U / 50
  f$q025 <- f$posterior_mean - .2
  f$q975 <- f$posterior_mean + .2
  baseline <- expand.grid(C = 1:4, fold = 1:5,
    model = c('UC-GP (x,c)', 'LVGP (x,c)', 'EzGP (x,c)'), stringsAsFactors = FALSE,
    KEEP.OUT.ATTRS = FALSE)
  baseline$posterior_mean <- baseline$C + .5
  list(u = list(display_points = u), m = m, f = f,
    grids = list(near = data.frame(U = seq(25, 175, by = 25), pressure = 650)),
    baseline = baseline)
}

test_that('EIV diagnostic plots render before any competitor fits exist', {
  d <- ocean_plot_fixture()
  for (baseline in list(NULL, d$baseline[FALSE, ])) {
    plots <- ocean_diag_plots(d$u, d$m, d$f, d$grids, 1L, baseline)
    expect_length(plots, 6L)
    for (name in c('fm_main', 'fm_all_folds')) {
      means <- plots[[name]][[2L]]$data
      expect_setequal(unique(means$model), c('EIV-GP O20', 'EIV-GP O50'))
      path <- tempfile(fileext = '.pdf')
      tryCatch({
        ggplot2::ggsave(path, plots[[name]], width = 11, height = 8)
        expect_gt(file.info(path)$size, 1000)
      }, finally = unlink(path))
    }
  }
})

test_that('available competitor curves remain correct for complete and partial results', {
  d <- ocean_plot_fixture()
  for (baseline in list(d$baseline, d$baseline[d$baseline$fold == 2L, ])) {
    plots <- ocean_diag_plots(d$u, d$m, d$f, d$grids, 1L, baseline)
    for (name in c('fm_main', 'fm_all_folds')) {
      means <- plots[[name]][[2L]]$data
      actual <- means[!grepl('^EIV', means$model), c('C', 'fold', 'model', 'posterior_mean')]
      expected <- if (name == 'fm_main') baseline[baseline$fold == 1L, ] else baseline
      rownames(actual) <- rownames(expected) <- NULL
      expect_equal(actual, expected)
      path <- tempfile(fileext = '.pdf')
      tryCatch({
        ggplot2::ggsave(path, plots[[name]], width = 11, height = 8)
        expect_gt(file.info(path)$size, 1000)
      }, finally = unlink(path))
    }
  }
})

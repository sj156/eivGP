repo<-Sys.getenv('EIVGP_TEST_REPO');stopifnot(nzchar(repo))
source(file.path(repo,'tools/ocean-tools.R'),local=TRUE)

test_that('plot exports create real PNG and PDF artifacts with the raster backend', {
  skip_if_not_installed('ragg')
  root<-tempfile();on.exit(unlink(root,recursive=TRUE),add=TRUE)
  p<-ggplot2::ggplot(data.frame(x=1:3,y=c(1,3,2)),ggplot2::aes(x,y))+ggplot2::geom_line()
  release<-mixedgp_budget_enter(2L);on.exit(release(),add=TRUE)
  mixedgp_save_plot_jobs(list(trace=p),root,function(n)c(4,3))
  paths<-file.path(root,c('trace.png','trace.pdf'))
  expect_true(all(file.exists(paths)))
  expect_true(all(file.info(paths)$size>1000))
  expect_identical(readBin(paths[1L],'raw',n=8L),as.raw(c(137,80,78,71,13,10,26,10)))
})

test_that('a PNG backend that silently fails to open is reported immediately', {
  helper<-mixedgp_png_device
  env<-new.env(parent=environment(helper));environment(helper)<-env
  env$requireNamespace<-function(...)FALSE
  testthat::local_mocked_bindings(png=function(...)invisible(NULL),.package='grDevices')
  expect_error(helper(tempfile(fileext='.png'),4,3),'Could not start the PNG device')
})

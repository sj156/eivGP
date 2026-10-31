test_that("calibration lines require distinct valid positions within each facet", {
  z<-data.frame(cell=c('a','a','b','b','c','c','d','d','e','e','e'),
    method='EIV-GP',n_calib=c(20,50,20,20,20,50,20,50,0,20,50),
    mean=c(1,2,1,2,1,NA,1,Inf,1,NA,2))
  keep<-mixedgp_line_rows(z,c('cell','method'))
  expect_identical(keep,z[z$cell%in%c('a','e'),,drop=FALSE])
  # Two singleton facets must not be connected, even with the same method.
  single<-data.frame(cell=c('a','b'),method='EIV-GP',n_calib=c(20,50),mean=c(1,2))
  expect_equal(nrow(mixedgp_line_rows(single,c('cell','method'))),0L)
  expect_identical(mixedgp_line_rows(z[FALSE,],c('cell','method')),z[FALSE,])
})
test_that("comparison plots keep singleton points and omit only invalid lines", {
  z<-data.frame(cell_id=c('one','two','two','two'),method=c('EIV-GP','EIV-GP','EIV-GP','Oracle'),
    n_calib=c(50,20,50,NA),mean=c(.3,.4,.2,.1),metric='RMSE')
  p<-plot_comparison(z)
  messages<-character()
  grob<-withCallingHandlers(ggplot2::ggplotGrob(p),message=function(m){messages<<-c(messages,conditionMessage(m))})
  expect_false(any(grepl('Each group consists',messages)))
  b<-ggplot2::ggplot_build(p)
  expect_equal(nrow(b$data[[1]]),nrow(z))
  expect_equal(sort(b$data[[1]]$y),sort(z$mean))
  expect_equal(nrow(b$data[[2]]),2L)
  expect_equal(length(unique(b$data[[2]]$PANEL)),1L)
  expect_equal(sort(b$data[[2]]$y),c(.2,.4))
  expect_s3_class(grob,'gtable')
  single<-plot_comparison(z[z$cell_id=='one',])
  expect_message(ggplot2::ggplotGrob(single),NA)
  expect_equal(nrow(ggplot2::ggplot_build(single)$data[[2]]),0L)
})

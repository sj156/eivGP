test_that('target traces preserve chains, retained iteration and full oxygen diagnostics', {
  n<-40L;nc<-4L
  u<-lapply(seq_len(nc),function(c)cbind(seq_len(n)/100+c,seq_len(n)/50+c,rep(c,n)))
  fit<-list(data=list(x_center=c(0,0),x_scale=c(1,1),x=matrix(1,3,2),y=1:3,
    y_center=10,y_scale=2,U_center=50,U_scale=10,miss_idx=1:2),
    mcmc=list(samples_by_chain=list(u=u),samples_u=do.call(rbind,u),
      samples_tau=matrix(rep(c(-1,0,1),each=n*nc),ncol=3),
      samples_logtheta=matrix(0,n*nc,3),samples_sigma2=rep(1,n*nc),
      mcmc_draw_info=data.frame(chain=rep(1:nc,each=n))))
  env<-new.env(parent=environment(ocean_target_trace_fit));fun<-ocean_target_trace_fit;environment(fun)<-env
  env$gp_integrated_mean_state_1d<-function(u_train,...)list(mean=u_train[1])
  env$predict_eivgp<-function(fit,draw_ids,...)matrix(draw_ids,nrow=1)
  spec<-ocean_target_trace_settings();spec$states<-10L;spec$nodes<-8L;spec$check_states<-3L
  train<-data.frame(calib_o50=c(FALSE,FALSE,TRUE),C=c(2,2,2))
  z<-fun(fit,train,1,'O20',burn=5L,spec=spec,
    coordinate_diagnostics=data.frame(parameter=c('u[1]','u[2]','u[3]'),rhat=c(1.1,1.2,5)))
  expect_equal(z$individual,2L)
  expect_equal(z$summary$diagnostic_draws,c(10,10,40,40))
  expect_equal(z$traces$value[z$traces$target=='m' & z$traces$chain==1],10+2*u[[1]][unique(round(seq(1,n,length.out=10))),1])
  expect_equal(range(z$traces$iteration),c(6,45))
  expect_equal(z$integration_check$max_abs_delta,0)
  expect_equal(nrow(z$traces),4*10*4)
  expect_equal(ocean_target_trace_selection()$fold,c(5,1,2,1,2,1,1,1,2))
  expect_s3_class(ocean_target_trace_plot(z$traces,data.frame(fold=1,cell='O20',target='m')),'ggplot')
})

test_that('complete-case manuscript table excludes different test protocols', {
  dest<-tempfile();dir.create(dest);on.exit(unlink(dest,recursive=TRUE))
  d<-data.frame(Method=c('EIV-GP','LVGP(CC)','EIV-GP'),
    Training=c('20% calibration; 400 rows','20% CC: 80 rows','20% calibration; 400 rows'),
    Test=c('20 (x,u)','20 (x,u,c)','20% (x,u) + 80% (x,c)'),Folds=5L,
    RMSE=c('6.598 (2.499)','8.019 (0.953)','8.792 (0.883)'),check.names=FALSE)
  for(k in c('CRPS/sY','Coverage (%)','Width/sY'))d[[k]]<-c('0.1 (0.01)','0.2 (0.01)','0.3 (0.01)')
  write.csv(d,file.path(dest,'tab_ocean_complete_case.csv'),row.names=FALSE)
  ocean_export_manuscript_tables(dest,list(response_sd=40.112538,design=list(N=500)))
  tex<-paste(readLines(file.path(dest,'tab_ocean_complete_case.tex')),collapse='\n')
  expect_match(tex,'tab:ocean-si-complete-case',fixed=TRUE)
  expect_match(tex,'\\textbf{6.598}',fixed=TRUE)
  expect_match(tex,'400 &',fixed=TRUE);expect_match(tex,'80 &',fixed=TRUE)
  expect_false(grepl('8.792',tex,fixed=TRUE))
  expect_equal(read.csv(file.path(dest,'tab_ocean_complete_case.csv'),check.names=FALSE),d)
})

test_that('paper trace axes match the appendix without mislabeling pilot traces', {
  z<-data.frame(fold=1L,cell='O20',target='m',chain=1L,
    iteration=c(5001L,50000L),value=c(1,2))
  selection<-data.frame(fold=1L,cell='O20',target='m')
  p<-ocean_target_trace_plot(z,selection)
  expect_identical(p$labels$x,'Iteration (thousands; burn-in ends at 5)')
  expect_equal(p$scales$get_scales('x')$breaks,c(10,30,50))
  z$iteration<-c(31L,80L)
  expect_identical(ocean_target_trace_plot(z,selection)$labels$x,'Iteration (thousands; retained states)')
})

test_that("summary plots preserve failures, nonfinite values and report thresholds", {
  d<-data.frame(parameter=c("a","b","c"),rhat=c(1.001,1.2,NA_real_),
    ess_bulk=c(1000,30,NA_real_),ess_tail=c(800,NA_real_,NA_real_),
    mcse_sd_ratio=c(.02,.2,NA_real_))
  p<-plot_mixing(d);expect_s3_class(p,"ggplot")
  expect_equal(p$data$max_rhat,1.2);expect_equal(p$data$min_bulk,30)
  q<-plot_mixing(d,"coverage");expect_equal(q$data$count,c(2,1))
  report<-structure(list(table=d,settings=list(rhat_limit=1.05,bulk_ess_limit=200)),class="eivgp_diagnostics")
  expect_equal(ggplot2::ggplot_build(plot_mixing(report))$data[[1]]$xintercept,1.05)
  expect_error(plot_mixing(data.frame(gate_pass=TRUE)),"Diagnostics need")
  expect_error(plot_mixing(d,rhat_limit=.9),"exceed 1")
  d$panel<-NA_character_;expect_equal(plot_mixing(d)$data$n,3)
})

test_that("chain plots use the full window and never silently discard bad draws", {
  a<-1:101; b<-rep(4,80);s<-list(a=list(as.numeric(a),b))
  p<-plot_mcmc(s,"a","running_mean",max_points=10)
  expect_equal(tail(p$data$value[p$data$chain==1],1),mean(a))
  expect_equal(tail(p$data$x[p$data$chain==1],1),101)
  q<-plot_mcmc(s,"a","rank",bins=10);expect_equal(as.numeric(tapply(q$data$value,q$data$chain,sum)),c(1,1))
  q<-plot_mcmc(s,"a","acf");expect_true(all(is.na(q$data$value[q$data$chain==2])))
  expect_match(q$labels$subtitle,"Constant")
  expect_error(plot_mcmc(s,"missing"),"available")
  s$a[[1]][3]<-NA_real_;expect_error(plot_mcmc(s),"silently dropped")
  expect_error(plot_mcmc(1:10),"Supply a fit")
})

test_that("directory reporting de-duplicates detail files and does not refit", {
  root<-tempfile();dir.create(root);on.exit(unlink(root,recursive=TRUE))
  d<-data.frame(parameter="r",rhat=1.02,ess_bulk=100,ess_tail=150)
  write.csv(d,file.path(root,"study1_mcmc_parameter_diagnostics.csv"),row.names=FALSE)
  write.csv(d,file.path(root,"study1_mcmc_diagnostic_details.csv"),row.names=FALSE)
  files<-save_mixing_plots(root,formats="pdf");expect_length(files,3);expect_true(all(file.exists(files)))
  manifest<-read.csv(file.path(root,"mixing-plots/diagnostic_sources.csv"))
  expect_equal(nrow(manifest),1);expect_match(manifest$path,"parameter")
  expect_warning(save_mixing_plots(file.path(root,"mixing-plots")),"No parameter-level")
})

test_that("study default plotting is reporting-only and can be disabled", {
  calls<-character();e<-new.env(parent=environment(run_paper_study));root<-tempfile();dir.create(root);dir.create(file.path(root,"core"));on.exit(unlink(root,recursive=TRUE))
  e$eivgp_numerical_api<-function(...)list(numerical_default_output=function(...)root,numerical_default_data=function(...)root)
  e$study_config<-function(...)list();e$run_numerical_study<-function(...)list()
  e$save_mixing_plots<-function(path,...){calls<<-c(calls,path);"figure.pdf"}
  run<-run_paper_study;environment(run)<-e
  args<-list(component="core",output=root,`data-root`=root,action="report")
  x<-run("study1",args);expect_length(calls,1);expect_identical(attr(x,"mixing_plots"),"figure.pdf")
  calls<-character();args$`mixing-plots`<-"false";run("study2",args);expect_length(calls,0)
  args$`mixing-plots`<-TRUE;args$action<-"plan";run("study2",args);expect_length(calls,0)
  args$`mixing-plots`<-"bogus";expect_error(run("study1",args),"true or false")
})

test_that("core competitor routing never requests EIV fitting", {
  f<-run_paper_study;env<-new.env(parent=environment(f));actions<-character();calls<-list()
  env$eivgp_numerical_api<-function(...)list()
  env$study_config<-function(...)list()
  env$run_numerical_study<-function(config,action,competitors) {actions<<-c(actions,action);list()}
  env$mixedgp_run_competitor_workflow<-function(...) {calls[[length(calls)+1L]]<<-list(...);list()}
  environment(f)<-env
  base<-list(component="core-competitors",output=tempfile(),`data-root`=tempfile(),cores="1")
  f("study1",modifyList(base,list(action="run")),tempdir())
  expect_identical(actions,"data");expect_equal(calls[[1]]$action,"run")
  actions<-character();calls<-list()
  f("study1",modifyList(base,list(action="report")),tempdir())
  expect_length(actions,0);expect_equal(calls[[1]]$action,"export")
})


test_that("Study I complete-case paper defaults preserve the published budget", {
  f<-run_paper_study;env<-new.env(parent=environment(f));seen<-NULL
  env$eivgp_numerical_api<-function(...)list()
  env$mixedgp_study1_cc_workflow<-function(args,...) {seen<<-args;list()}
  environment(f)<-env
  f("study1",list(component="complete-case",action="plan",profile="paper",output=tempfile(),`data-root`=tempfile()),tempdir())
  expect_identical(seen$iter,"12000");expect_identical(seen$warmup,"4000")
  expect_identical(seen$outer,"16");expect_identical(seen$latent,"128")
  expect_identical(seen[["cc-response-latent"]],"64")
})

test_that("publication method partitions preserve the scientific plan", {
  f<-mixedgp_publication_workflow
  base<-list(stage="full",output=tempfile(),`data-root`=tempfile())
  capture.output(all<-f("study2",base,tempdir()))
  capture.output(eiv<-f("study2",modifyList(base,list(methods="FULL,CUT")),tempdir()))
  capture.output(cc<-f("study2",modifyList(base,list(methods="UC-GP-CC,LVGP-CC,EzGP-CC")),tempdir()))
  key<-function(d)paste(d$rep,d$calib,d$method)
  expect_equal(nrow(eiv$jobs),200L);expect_equal(nrow(cc$jobs),300L)
  expect_length(intersect(key(eiv$jobs),key(cc$jobs)),0L)
  expect_setequal(c(key(eiv$jobs),key(cc$jobs)),key(all$jobs))
  for(x in list(eiv,cc)) {
    a<-all$commands$reference;b<-x$commands$reference
    a$methods<-b$methods<-NULL
    expect_identical(a,b)
  }
  capture.output(missing20<-f("study2",modifyList(base,list(methods="UC-GP-CC,LVGP-CC,EzGP-CC",reps="29:50",calib="20")),tempdir()))
  capture.output(missing50<-f("study2",modifyList(base,list(methods="UC-GP-CC,LVGP-CC,EzGP-CC",reps="15:50",calib="50")),tempdir()))
  expect_equal(nrow(missing20$jobs),66L);expect_equal(nrow(missing50$jobs),108L)
  expect_false(dir.exists(base$output))
  expect_error(f("study2",modifyList(base,list(methods="CUT,CUT")),tempdir()),"distinct subset")
  expect_error(f("study2",modifyList(base,list(methods="UC-GP")),tempdir()),"distinct subset")
  expect_error(f("study2",modifyList(base,list(calib="10")),tempdir()),"calib")
})

test_that("competitor-only publication dispatch does not require FULL archives", {
  f<-mixedgp_publication_workflow
  # Clone the dispatcher into an environment whose terminal route cannot fit.
  env<-new.env(parent=environment(f));seen<-list()
  env$run_paper_study<-function(study,args,repo_root) {seen[[length(seen)+1L]]<<-args;invisible(NULL)}
  environment(f)<-env
  root<-tempfile();on.exit(unlink(root,recursive=TRUE))
  a<-list(stage="full",action="run",reps="29",calib="20",methods="UC-GP-CC",output=root,`data-root`=tempfile())
  f("study2",a,tempdir())
  expect_length(seen,1L);expect_equal(seen[[1]]$methods,"UC-GP-CC")
  expect_false("full-fit-map"%in%names(seen[[1]]))
  expect_error(f("study2",modifyList(a,list(methods="FULL")),tempdir()),"full-fit-map")
})

test_that("core competitor routing never requests EIV fitting", {
  f<-run_paper_study;env<-new.env(parent=environment(f));actions<-character();calls<-list()
  env$eivgp_numerical_api<-function(...)list()
  env$study_config<-function(...)list()
  env$run_numerical_study<-function(config,action,competitors) {actions<<-c(actions,action);list()}
  env$mixedgp_run_competitor_workflow<-function(...) {calls[[length(calls)+1L]]<<-list(...);list()}
  environment(f)<-env
  base<-list(component="core-competitors",output=tempfile(),`data-root`=tempfile(),cores="1")
  f("study2",modifyList(base,list(action="run")),tempdir())
  expect_identical(actions,"data");expect_equal(calls[[1]]$action,"run")
  actions<-character();calls<-list()
  f("study2",modifyList(base,list(action="report")),tempdir())
  expect_length(actions,0);expect_equal(calls[[1]]$action,"export")
})

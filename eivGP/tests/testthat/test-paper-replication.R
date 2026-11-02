test_that("paper profiles retain reported cells without superseded component ablations", {
  for(study in c("study1","study2")) {
    cfg<-study_config(study,"paper",tempfile(),tempfile(),cores=12L)
    expect_true(all(vapply(cfg$cells,function(c)!c$run_ablations,logical(1))))
    expect_identical(cfg$mcmc$n_chains,4L)
    expect_identical(cfg$mcmc$n_iter,20000L)
    expect_identical(cfg$mcmc$burn,5000L)
    expect_lte(cfg$parallel$workers*cfg$parallel$chain_workers,12L)
    if(study=="study2") {
      expect_null(cfg$mcmc$ess_min)
      expect_identical(vapply(cfg$cells,`[[`,"","id"),c("primary_q2","primary_q4_calibration","logistic_q4"))
      expect_identical(lapply(cfg$cells,`[[`,"calibration_grid"),list(50L,c(0L,20L,50L,80L),50L))
    } else expect_identical(vapply(cfg$cells,`[[`,0,"heterogeneity_eta"),c(0,1))
  }
})

test_that("paired synthetic designs preserve data and common calibration masks", {
  for(rr in c(1L,17L)) {
    a<-make_study1_synthetic_dataset(rr,n=100L,n_test=100L,scenario="heterogeneity_continuum",threshold_design="balanced",heterogeneity_eta=0,calib_grid=c(0L,20L,50L))
    b<-make_study1_synthetic_dataset(rr,n=100L,n_test=100L,scenario="heterogeneity_continuum",threshold_design="balanced",heterogeneity_eta=1,calib_grid=c(0L,20L,50L))
    expect_identical(a$calibration_sets,b$calibration_sets)
    expect_identical(a$data$train$u,b$data$train$u)
    a<-make_study2_synthetic_dataset(rr,n=100L,n_test=100L,q=2L,calib_grid=50L)
    b<-make_study2_synthetic_dataset(rr,n=100L,n_test=100L,q=4L,calib_grid=c(0L,20L,50L,80L))
    expect_identical(a$data$train$U,b$data$train$U)
    expect_identical(a$data$test$y,b$data$test$y)
    expect_identical(a$calibration_sets[["50"]],b$calibration_sets[["50"]])
  }
})

# Extract the nested cache reader without executing its surrounding simulation.
test_that("Study I cache hits do not depend on Study II globals", {
  assignments<-Filter(function(x)is.call(x)&&identical(x[[1L]],as.name("<-"))&&
    identical(x[[2L]],as.name("run_or_load_study1_replication")),as.list(body(eivgp_study1_driver)))
  expect_length(assignments,1L)
  definition<-assignments[[1L]][[3L]]
  expect_false(is.null(definition))
  env<-new.env(parent=environment(eivgp_study1_driver))
  root<-tempfile();dir.create(root);on.exit(unlink(root,recursive=TRUE))
  env$REP_DIR<-env$STUDY1_DATA_DIR<-root
  env$STUDY1_REUSE_LOCKED_EIV<-FALSE;env$STUDY1_USE_CACHE<-TRUE
  env$STUDY1_DESIGN_TAG<-"test";env$STUDY1_SCENARIO<-"heterogeneity_continuum"
  env$n_train<-env$n_test<-100L;env$m<-6L
  path<-file.path(root,mixedgp_dataset_filename("study1",1L,env$STUDY1_SCENARIO,100L,100L,6L))
  saveRDS(list(frozen=TRUE),path)
  cached<-list(metadata=list(rep=1L,design_tag="test",frozen_data_file=basename(path),frozen_data_md5=unname(tools::md5sum(path))))
  saveRDS(cached,file.path(root,"study1_rep001_test_fresh.rds"))
  env$run_one_study1_replication<-function(...)stop("Unexpected refit")
  f<-eval(definition,env)
  expect_true(attr(f(1L),"eivgp_cache_hit"))
})

test_that("paper lanes select the exact appendix pairs without truncating FULL or CC", {
  repo<-Sys.getenv("EIVGP_REPO_ROOT")
  skip_if(!file.exists(file.path(repo,"paper/numerical/study2-common.R")))
  env<-new.env();sys.source(file.path(repo,"paper/numerical/study2-common.R"),env)
  seen<-list()
  runner<-function(study,args,repo_root) {seen[[length(seen)+1L]]<<-args;list()}
  env$study2_paper_workflow("eivGP",list(output=tempfile()),repo,runner)
  expect_length(seen,4L)
  expect_identical(seen[[1]]$component,"core")
  expect_identical(seen[[2]]$methods,"FULL")
  expect_null(seen[[2]]$reps)
  expect_identical(seen[[3]]$calib,"20")
  expect_identical(as.integer(strsplit(seen[[3]]$reps,",",fixed=TRUE)[[1]]),1:28)
  expect_identical(seen[[4]]$calib,"50")
  expect_identical(as.integer(strsplit(seen[[4]]$reps,",",fixed=TRUE)[[1]]),1:14)
  expect_true(all(vapply(seen,function(a)identical(a$cores,"12"),logical(1))))
  seen<-list()
  env$study2_paper_workflow("competitors",list(output=tempfile()),repo,runner)
  expect_length(seen,2L);expect_null(seen[[2]]$reps)
  expect_identical(seen[[2]]$methods,"UC-GP-CC,LVGP-CC,EzGP-CC")
  seen<-list()
  env$study2_paper_workflow("eivGP",list(action="report",methods="CUT",output=tempfile()),repo,runner)
  expect_length(seen,2L)
  expect_true(all(vapply(seen,function(a)identical(a$action,"report"),logical(1))))
  expect_error(env$study2_paper_workflow("eivGP",list(latent="256"),repo,runner),"Paper fixes")
  expect_error(env$study2_paper_workflow("eivGP",list(latnett="128"),repo,runner),"Unknown")
  expect_error(env$study2_paper_workflow("eivGP",list(cores="2"),repo,runner),"four chains")
  sys.source(file.path(repo,"paper/numerical/study1-common.R"),env)
  seen<-list()
  env$study1_paper_workflow("competitors",list(output=tempfile()),repo,runner)
  expect_length(seen,2L);expect_null(seen[[2]]$reps)
  expect_identical(seen[[2]]$methods,"UC-GP-CC,LVGP-CC,EzGP-CC")
  expect_error(env$study1_paper_workflow("competitors",list(methods="LVGP-CC",seed="2"),repo,runner),"Paper seeds")
})

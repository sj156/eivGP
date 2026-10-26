repo<-Sys.getenv('EIVGP_TEST_REPO');stopifnot(nzchar(repo))
source(file.path(repo,'tools/ocean-tools.R'),local=TRUE)
test_that("flat paths are consistent across incremental builders and legacy readers", {
  root<-tempfile();on.exit(unlink(root,recursive=TRUE),add=TRUE)
  ocean_storage_init(root)
  expect_true(ocean_storage_flat(root))
  expect_equal(ocean_path(root,"full5_class4","fold01","O20","fit.rds"),file.path(root,"eiv_fold01_O20","fit.rds"))
  base<-ocean_path(root,"full5_class4");fold<-ocean_path(base,"fold01")
  expect_equal(ocean_path(fold,"O20","fit.rds"),file.path(root,"eiv_fold01_O20","fit.rds"))
  expect_equal(ocean_path(fold,"baseline_metrics.csv"),file.path(root,"baseline_fold01","baseline_metrics.csv"))
  expect_equal(ocean_path(root,"continuation_plus20k","fold01","O20","fit_plus20k.rds"),file.path(root,"eiv_fold01_O20_continued","fit_plus20k.rds"))
  base<-ocean_path(root,"complete_case_baselines","fits")
  expect_equal(ocean_path(base,"cc_fold01_O20_EzGP","package_output.rds"),file.path(root,"cc_fold01_O20_EzGP","package_output.rds"))
  expect_equal(ocean_path(root,"workflow_identity.rds"),file.path(root,"_meta","workflow_identity.rds"))
  expect_equal(ocean_path(root,"diagnostics","tables","a.csv"),file.path(root,"reports","diagnostics","tables","a.csv"))
  expect_equal(ocean_path(ocean_job_dir(root,"paper_fit"),"tasks","a.log"),file.path(root,"logs","paper_fit","a.log"))
  expect_equal(ocean_path(root,"full5_class4",sprintf("fold%02d",1:2)),file.path(root,paste0("baseline_fold0",1:2)))
  legacy<-tempfile();dir.create(legacy);on.exit(unlink(legacy,recursive=TRUE),add=TRUE)
  saveRDS(list(profile="paper"),file.path(legacy,"workflow_identity.rds"))
  ocean_storage_init(legacy);expect_false(ocean_storage_flat(legacy))
  expect_equal(ocean_path(legacy,"full5_class4","fold01","O20"),file.path(legacy,"full5_class4","fold01","O20"))
})

test_that("migration planning is read only; fixture migration preserves fits and receipts", {
  parent<-tempfile();dir.create(parent);on.exit(unlink(parent,recursive=TRUE),add=TRUE)
  old<-file.path(parent,"old");new<-file.path(parent,"new");dir.create(old)
  identity<-list(profile="paper",numerics=list(chains=4L),code="historical-code")
  saveRDS(identity,file.path(old,"workflow_identity.rds"))
  rel<-"full5_class4/fold01/O20/fit.rds";dir.create(dirname(file.path(old,rel)),recursive=TRUE)
  saveRDS(list(chain_state=1:4),file.path(old,rel))
  hash<-unname(tools::md5sum(file.path(old,rel)))
  dir.create(file.path(old,"workflow_tasks"))
  saveRDS(list(identity=identity,files=rel,hashes=hash),file.path(old,"workflow_tasks","fit.done.rds"))
  before<-tools::md5sum(list.files(old,recursive=TRUE,full.names=TRUE))
  plan<-ocean_migrate_storage(old,new)
  expect_false(dir.exists(new));expect_true("eiv_fold01_O20/fit.rds"%in%plan$destination)
  expect_identical(before,tools::md5sum(list.files(old,recursive=TRUE,full.names=TRUE)))
  dir.create(file.path(old,".workflow.lock"))
  expect_error(ocean_migrate_storage(old,new),"locked")
  unlink(file.path(old,".workflow.lock"),recursive=TRUE)
  ocean_migrate_storage(old,new,dry_run=FALSE)
  expect_true(ocean_storage_flat(new));expect_true(ocean_verify_receipt(new,"fit",identity))
  expect_identical(readRDS(ocean_path(new,"workflow_identity.rds")),identity)
  expect_equal(unname(tools::md5sum(ocean_path(new,rel))),hash)
  expect_identical(before,tools::md5sum(list.files(old,recursive=TRUE,full.names=TRUE)))
  expect_error(ocean_verify_receipt(new,"fit",list(code="different")),"identity mismatch")
  expect_error(ocean_migrate_storage(old,new),"exists")
  expect_error(ocean_migrate_storage(old,file.path(old,"nested")),"outside")
})

test_that("flat continuation and report source selection retain relative task identity", {
  root<-tempfile();on.exit(unlink(root,recursive=TRUE),add=TRUE);ocean_storage_init(root)
  mixedgp_atomic_rds(list(profile="quick"),ocean_path(root,"workflow_identity.rds"))
  for(cell in c("O0","O20","O50")) {
    p<-ocean_path(root,"smoke_class4","fold01",cell)
    mixedgp_atomic_csv(data.frame(rhat=1,ess_bulk=500,ess_tail=500),ocean_path(p,"convergence_parameters.csv"))
    mixedgp_atomic_csv(data.frame(score=1),ocean_path(p,"predictive_metrics.csv"))
  }
  expect_equal(ocean_saved_profile(root),"quick")
  selected<-ocean_final_eiv_sources(root,"quick")
  expect_true(all(selected$available));expect_equal(selected$source,paste0("eiv_fold01_",c("O0","O20","O50")))
  ocean_continuation_iterations(root,"pilot",1,"O20",requested=7L,register=TRUE)
  expect_equal(ocean_continuation_lengths(root,"pilot"),7L)
  expect_error(ocean_continuation_iterations(root,"pilot",1,"O20",requested=8L,register=TRUE),"differs")
})

test_that("flat workflow writes receipts, reuses completed tasks and reports without fitting", {
  root<-tempfile();on.exit(unlink(root,recursive=TRUE),add=TRUE)
  repo<-normalizePath(testthat::test_path("../.."))
  data<-file.path(repo,"data/ocean_silicate/n200")
  skip_if_not(dir.exists(data))
  env<-new.env(parent=environment(ocean_execute));execute<-ocean_execute;environment(execute)<-env
  env$ocean_identity<-function(...)list(profile="pilot",code="fixture")
  calls<-0L
  env$ocean_run_task<-function(task,repository_dir,output_dir,settings,inputs,identity) {
    calls<<-calls+1L
    reader<-new.env(parent=environment(ocean_execute))
    reader$root<-output_dir
    ocean_initial_functions(reader)
    copied<-reader$read_fold(task$fold)
    expect_identical(copied$train$cohort_id,
      inputs$folds[[paste(task$fold,"O20",sep="/")]]$train$cohort_id)
    expect_identical(ocean_hashes(file.path(output_dir,"data","freeze_manifest.json")),
      ocean_hashes(file.path(data,"freeze_manifest.json")))
    for(p in ocean_task_files(task,settings$profile,settings$action)) {
      dest<-ocean_path(output_dir,p)
      if(grepl("[.]rds$",p))mixedgp_atomic_rds(list(chains=1:4),dest)
      else mixedgp_atomic_csv(data.frame(rhat=1,ess_bulk=500,ess_tail=500),dest)
    }
  }
  cfg<-ocean_settings(profile="pilot",action="fit",folds=1L,cells="O20",modules="eiv",cores=1L)
  execute(repo,data,root,cfg);execute(repo,data,root,cfg)
  expect_equal(calls,1L)
  expect_true(file.exists(file.path(root,"eiv_fold01_O20","fit.rds")))
  expect_true(file.exists(file.path(root,"_meta/tasks/fit_eiv_fold01_O20.done.rds")))
  expect_false(dir.exists(file.path(root,"full5_class4")))
  cfg$action<-"report";execute(repo,data,root,cfg)
  expect_equal(calls,1L)
  expect_false(dir.exists(file.path(root,".workflow.lock")))
})

test_that("completed tasks import across layouts in both directions", {
  parent<-tempfile();dir.create(parent);on.exit(unlink(parent,recursive=TRUE),add=TRUE)
  old<-file.path(parent,"old");flat<-file.path(parent,"flat");back<-file.path(parent,"back")
  dir.create(old);ocean_storage_init(flat);dir.create(back)
  identity<-list(profile="paper",code="fixture")
  for(root in c(old,flat,back))mixedgp_atomic_rds(identity,ocean_path(root,"workflow_identity.rds"))
  rel<-"full5_class4/fold01/O20/fit.rds"
  mixedgp_atomic_rds(list(chains=1:4),file.path(old,rel))
  mixedgp_atomic_rds(list(identity=identity,files=rel,hashes=ocean_hashes(file.path(old,rel))),file.path(old,"workflow_tasks/task.done.rds"))
  ocean_combine_results(old,flat,identity)
  expect_true(ocean_verify_receipt(flat,"task",identity))
  expect_equal(readRDS(ocean_path(flat,"workflow_tasks/task.done.rds"))$files,"eiv_fold01_O20/fit.rds")
  ocean_combine_results(flat,back,identity)
  expect_true(ocean_verify_receipt(back,"task",identity))
  expect_true(file.exists(file.path(back,rel)))
})

test_that("standalone CLI selects its library before resolving storage paths", {
  repo<-normalizePath(testthat::test_path("../.."))
  skip_if_not(dir.exists(file.path(repo,"data/ocean_silicate/n200")))
  log<-system2(file.path(R.home("bin"),"Rscript"),
    c(shQuote(file.path(repo,"tools/ocean.R")),"--action=plan","--profile=paper","--dataset=ocean200","--cores=1"),
    env=paste0("EIVGP_R_LIBRARY=",shQuote(dirname(find.package("eivGP")))),stdout=TRUE,stderr=TRUE)
  expect_null(attr(log,"status"))
  expect_true(any(grepl("Total shared core budget:",log,fixed=TRUE)))
})

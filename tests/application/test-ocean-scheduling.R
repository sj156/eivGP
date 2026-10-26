repo<-Sys.getenv('EIVGP_TEST_REPO');stopifnot(nzchar(repo))
source(file.path(repo,'tools/ocean-tools.R'),local=TRUE)
test_that("Ocean resource plans accept different machines and configurable chain counts", {
  for(cores in c(1L,3L,12L,16L,56L)) {
    p<-mixedgp_core_plan(cores,15L,4L)
    expect_equal(p$active_cores,min(cores,60L))
    expect_true(all(p$task_cores>=1 & p$task_cores<=4))
    expect_equal(sum(p$task_cores),p$active_cores)
  }
  expect_identical(vapply(c(12L,16L,56L),function(x)mixedgp_core_plan(x,15L)$concurrent_fits,integer(1)),c(3L,4L,14L))
  expect_equal(mixedgp_core_plan(56,5,4)$active_cores,20)
  expect_identical(mixedgp_core_plan(56,5,12)$task_cores,c(12L,12L,12L,12L,8L))
  expect_error(ocean_settings(profile="pilot",chains=8),"chains=min")
  expect_error(ocean_settings(chains=0),"chains")
  expect_error(ocean_settings(chains=1.5),"chains")
  expect_identical(ocean_settings(profile="quick",cores=2)$chains,2L)
  expect_identical(mixedgp_core_plan(56,0,4)$task_cores,integer())
})

test_that("the shared queue bounds concurrent reservations and runs all tasks", {
  skip_if(.Platform$OS.type!="unix")
  starts<-list();ends<-list();indices<-seq_len(5L);caps<-c(4L,4L,1L,1L,1L)
  mixedgp_schedule(indices,caps,6L,
    function(i,n) {Sys.sleep(.1);list(i=i,cores=n,pid=Sys.getpid())},
    function(i,n,pid,active)starts[[as.character(i)]]<<-list(n=n,pid=pid,active=active),
    function(i,result,active)ends[[as.character(i)]]<<-list(result=result,active=active))
  expect_equal(length(starts),5L);expect_equal(length(ends),5L)
  expect_identical(starts[["1"]]$n,4L);expect_identical(starts[["2"]]$n,2L)
  expect_equal(max(vapply(starts,`[[`,integer(1),"active")),6L)
  for(i in indices) {
    key<-as.character(i)
    expect_lte(starts[[key]]$n,caps[i]);expect_lte(starts[[key]]$active,6)
    expect_identical(ends[[key]]$result$i,i)
    expect_identical(ends[[key]]$result$cores,starts[[key]]$n)
    expect_identical(ends[[key]]$result$pid,starts[[key]]$pid)
    expect_false(starts[[key]]$pid==Sys.getpid())
  }
})

test_that("single-core analysis stays in the coordinator process", {
  result<-list();active<-integer()
  mixedgp_schedule(1:2,c(4L,4L),1L,function(i,n)list(pid=Sys.getpid(),cores=n),
    function(i,n,pid,a)active<<-c(active,a),
    function(i,x,a)result[[i]]<<-x)
  expect_identical(active,c(1L,1L))
  expect_true(all(vapply(result,function(x)x$pid==Sys.getpid() && x$cores==1L,logical(1))))
})

test_that("imports retain both fitting and analysis receipts", {
  root<-tempfile();source<-file.path(root,"source");output<-file.path(root,"output")
  dir.create(source,recursive=TRUE);dir.create(output)
  on.exit(unlink(root,recursive=TRUE),add=TRUE)
  identity<-list(version="test");mixedgp_atomic_rds(identity,file.path(source,"workflow_identity.rds"))
  task<-data.frame(task_id="eiv_fold01_O20")
  for(id in c("fit_eiv_fold01_O20","eiv_fold01_O20")) {
    payload<-paste0(id,".rds");mixedgp_atomic_rds(id,file.path(source,payload))
    mixedgp_atomic_rds(list(identity=identity,task=task,files=payload,
      hashes=ocean_hashes(file.path(source,payload))),file.path(source,"workflow_tasks",paste0(id,".done.rds")))
  }
  ocean_combine_results(source,output,identity)
  expect_true(ocean_verify_receipt(output,"fit_eiv_fold01_O20",identity))
  expect_true(ocean_verify_receipt(output,"eiv_fold01_O20",identity))
})


test_that("ADNI and Ocean use identical cores-only resident allocations", {
  adni<-new.env(parent=globalenv())
  sys.source(file.path(repo,'tools/adni-tools.R'),adni)
  for(cores in 1:16)for(chains in c(1L,4L,8L))for(tasks in c(1L,3L,15L)) {
    expected<-as.integer(min(tasks,max(1L,floor(cores/chains))))
    for(stage in c('fit','continue','analyze','functionals','predict'))
      expect_identical(mixedgp_stage_workers(cores,stage,tasks,chains),expected)
    for(phase in c('fit','evaluate'))
      expect_identical(adni$adni_resident_workers(cores,phase,tasks,chains),expected)
  }
})


test_that("core budgets select actual chains and leave competitors single-core", {
  adni<-new.env(parent=globalenv());sys.source(file.path(repo,'tools/adni-tools.R'),adni)
  for(n in c(1L,2L,3L,4L,6L,8L,10L,12L)) {
    for(profile in c('quick','pilot','paper')) {
      expect_identical(adni$adni_settings(profile,cores=n)$chains,min(n,4L))
      expect_identical(ocean_settings(profile,cores=n)$chains,min(n,4L))
      expect_identical(ocean_settings(profile,cores=n)$cores,n)
    }
    expect_equal(mixedgp_application_allocation(n)$jobs,n %/% min(n,4L))
    expect_equal(mixedgp_stage_workers(n,'fit',30L,1L),n)
  }
  expect_false(ocean_cc_controls('paper')$LVGP$parallel)
  expect_error(ocean_settings('paper',cores=2,ess_min=100),'four chains')
  expect_error(adni$adni_settings('paper',cores=2,ess_min=100),'four chains')
  expect_identical(ocean_settings('paper',action='report',cores=2,chains=4)$chains,4L)
  expect_error(ocean_settings('paper',action='continue',cores=2,chains=4),'chains=min')
})

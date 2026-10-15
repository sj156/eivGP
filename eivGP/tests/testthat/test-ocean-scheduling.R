test_that("Ocean resource plans accept different machines and configurable chain counts", {
  for(cores in c(1L,3L,12L,16L,56L)) {
    p<-ocean_core_plan(cores,15L,4L)
    expect_equal(p$active_cores,min(cores,60L))
    expect_true(all(p$task_cores>=1 & p$task_cores<=4))
    expect_equal(sum(p$task_cores),p$active_cores)
  }
  expect_identical(vapply(c(12L,16L,56L),function(x)ocean_core_plan(x,15L)$concurrent_fits,integer(1)),c(3L,4L,14L))
  expect_equal(ocean_core_plan(56,5,4)$active_cores,20)
  expect_identical(ocean_core_plan(56,5,12)$task_cores,c(12L,12L,12L,12L,8L))
  expect_identical(ocean_settings(profile="pilot",chains=8)$chains,8L)
  expect_error(ocean_settings(chains=0),"chains")
  expect_error(ocean_settings(chains=1.5),"chains")
  expect_error(ocean_settings(profile="quick",chains=4),"quick")
  expect_identical(ocean_core_plan(56,0,4)$task_cores,integer())
})

test_that("the shared queue bounds concurrent reservations and runs all tasks", {
  skip_if(.Platform$OS.type!="unix")
  starts<-list();ends<-list();indices<-seq_len(5L);caps<-c(4L,4L,1L,1L,1L)
  ocean_schedule(indices,caps,6L,
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
  ocean_schedule(1:2,c(4L,4L),1L,function(i,n)list(pid=Sys.getpid(),cores=n),
    function(i,n,pid,a)active<<-c(active,a),
    function(i,x,a)result[[i]]<<-x)
  expect_identical(active,c(1L,1L))
  expect_true(all(vapply(result,function(x)x$pid==Sys.getpid() && x$cores==1L,logical(1))))
})

test_that("imports retain both fitting and analysis receipts", {
  root<-tempfile();source<-file.path(root,"source");output<-file.path(root,"output")
  dir.create(source,recursive=TRUE);dir.create(output)
  on.exit(unlink(root,recursive=TRUE),add=TRUE)
  identity<-list(version="test");ocean_atomic_rds(identity,file.path(source,"workflow_identity.rds"))
  task<-data.frame(task_id="eiv_fold01_O20")
  for(id in c("fit_eiv_fold01_O20","eiv_fold01_O20")) {
    payload<-paste0(id,".rds");ocean_atomic_rds(id,file.path(source,payload))
    ocean_atomic_rds(list(identity=identity,task=task,files=payload,
      hashes=ocean_hashes(file.path(source,payload))),file.path(source,"workflow_tasks",paste0(id,".done.rds")))
  }
  ocean_combine_results(source,output,identity)
  expect_true(ocean_verify_receipt(output,"fit_eiv_fold01_O20",identity))
  expect_true(ocean_verify_receipt(output,"eiv_fold01_O20",identity))
})

source('codes/simulation_helpers.R');source('codes/publication_recovery.R')
e<-mixedgp_simulation_engine('codes')
a<-mixedgp_recovery_core_settings(56L)
stopifnot(a$dataset_workers==14L,a$chain_workers==4L,a$competitor_workers==56L,
 mixedgp_recovery_core_settings(56L,dataset_workers=7L)$dataset_workers==7L)
# Migration permits only the exact known runner and otherwise unchanged identity.
cur<-list(code_md5=c('/a/publication_recovery.R'='new','/a/sampler.R'='fixed'),R_version='same')
old<-cur;old$code_md5[1]<-'2841978643ece443c3a17525459aed6b'
stopifnot(mixedgp_recovery_scheduler_upgrade(old,cur))
bad<-old;bad$code_md5[2]<-'changed';stopifnot(!mixedgp_recovery_scheduler_upgrade(bad,cur))
bad<-old;bad$R_version<-'changed';stopifnot(!mixedgp_recovery_scheduler_upgrade(bad,cur))
if(.Platform$OS.type!='windows'){
  events<-tempfile();dir.create(events)
  result<-list();order<-integer()
  worker<-function(i){
    saveRDS(as.numeric(Sys.time()),file.path(events,paste0('start',i)))
    if(i==4L)stop('deliberate worker failure')
    # Actual nested forked chains under a two-dataset x two-chain allocation.
    chains<-e$mixedgp_parallel_lapply(as.list(1:2),function(j){
      Sys.sleep(if(i==1L)1.2 else .05)
      list(pid=Sys.getpid(),draw=runif(3))
    },n_cores=2L,seeds=as.integer(100*i+1:2))
    saveRDS(as.numeric(Sys.time()),file.path(events,paste0('end',i)))
    list(dataset_pid=Sys.getpid(),chains=chains)
  }
  mixedgp_recovery_dispatch(as.list(1:4),worker,function(ans,i){
    result[[as.character(i)]]<<-ans;order<<-c(order,i)
  },workers=2L)
  stopifnot(length(result)==4L,!is.null(result[['4']]$worker_error),
    readRDS(file.path(events,'start3'))<readRDS(file.path(events,'end1')))
  for(i in 1:3){
    z<-result[[as.character(i)]]
    stopifnot(z$dataset_pid!=Sys.getpid(),length(unique(vapply(z$chains,`[[`,integer(1),'pid')))==2L)
    for(j in 1:2){set.seed(100*i+j);stopifnot(identical(z$chains[[j]]$draw,runif(3)))}
  }
  message('Dynamic refill, nested forked chains, deterministic seeds, and worker-error isolation passed.')
}

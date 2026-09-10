source("codes/simulation_helpers.R");source("codes/publication_recovery.R")
real_engine <- mixedgp_simulation_engine
e <- real_engine("codes")
root<-tempfile("recovery-integration-");dir.create(root)
archive<-file.path(root,"archive");dir.create(archive)
data_root<-file.path(root,"data")
run<-file.path(archive,"study2-publication-fixture");dir.create(file.path(run,"combined"),recursive=TRUE)
cfg<-study2_simulation_config("publication",code_dir="codes",core_budget=1L,data_root=file.path(data_root,"study2"))
cfg$cells<-cfg$cells[1L];cell<-cfg$cells[[1]]
manifest<-mixedgp_generate_cell_data(cfg,cell,e);manifest$cell_id<-cell$id
saveRDS(list(config=cfg,data_manifest=manifest),file.path(run,"run_summary.rds"))
fit_dir<-file.path(run,"cells",cell$id,"results");dir.create(fit_dir,recursive=TRUE)
saveRDS(list(fit=list(sampler_version="0.3.1",priors=e$mixedgp_v030_priors(p=2,d=2,m_vec=rep(4L,cell$q)),
 control=e$mixedgp_v030_control(n_mis=50L),kernel=list(name="se",matern_nu=2.5),
 data=list(q=cell$q,d=2L,ident="lower_triangular"))),file.path(fit_dir,"flagged_fit_primary_rep001_cal050.rds"))

base<-data.frame(cell_id=cell$id,rep=1:50,n_calib=50L,method="EIV-GP",evaluation_stratum="overall",
 RMSE=.2,MAE=.1,CRPS=.1,NLPD=0,Coverage95=.95,Width95=1,IntervalScore95=1)
saveRDS(list(predictive_metrics=base[-1,]),file.path(run,"combined","all_raw_outputs.rds"))
comp<-do.call(rbind,lapply(c("UC-GP","LVGP","EzGP"),function(m){z<-base;z$method<-m;z$n_calib<-NA_integer_;z}))
missing<-with(comp,(method=="UC-GP"&rep==2)|(method=="LVGP"&rep==3)|(method=="EzGP"&rep==4))
comp<-comp[!missing,]
report<-file.path(archive,"competitor-reports","study2","publication");dir.create(report,recursive=TRUE)
write.csv(comp,file.path(report,"metrics.csv"),row.names=FALSE)
status<-comp[,c("cell_id","rep","method")];status$status<-"success";status$message<-""
write.csv(status,file.path(report,"statuses.csv"),row.names=FALSE)
pred<-do.call(rbind,lapply(seq_len(nrow(comp)),function(i){
 d<-mixedgp_read_cell_replication(cfg,cell,comp$rep[i])$data
 data.frame(cell_id=cell$id,rep=comp$rep[i],method=comp$method[i],test_id=seq_along(d$test$y),y=d$test$y,mean=0,variance=1)
}))
write.csv(pred,file.path(report,"predictions.csv"),row.names=FALSE)
expected_chain_workers<-4L
fit_calls<-competitor_calls<-0L
mixedgp_recovery_source_cell <- function(path,env){
 fit_calls<<-fit_calls+1L
 stopifnot(env$STUDY2_PARALLEL_LEVEL=="hybrid",env$STUDY2_CHAIN_WORKERS==expected_chain_workers,env$STUDY2_DATASET_WORKERS==1L,env$STUDY2_REP_IDS==1L,identical(env$STUDY2_RECOVERY_MODEL$priors$signal_shape,c(13,3)))
 env$raw_outputs<-list(predictive_metrics=base[1,,drop=FALSE])
}
e$mixedgp_cached_competitors<-function(X_train,y_train,C_train,X_test,C_test,n_draw,seed,m_vec,methods,controls,cache_root,allow_fit,retry_failed){
 competitor_calls<<-competitor_calls+1L;method<-methods[1]
 draw<-e$attach_predictive_normal_components(matrix(rep(qnorm((seq_len(n_draw)-.5)/n_draw),nrow(X_test)),n_draw,nrow(X_test)),rep(0,nrow(X_test)),rep(1,nrow(X_test)))
 list(status=e$mixedgp_competitor_status(method,status="success"),draws=setNames(list(draw),method),
 predictive_means=setNames(list(rep(0,nrow(X_test))),method),predictive_variances=setNames(list(rep(1,nrow(X_test))),method))
}
mixedgp_simulation_engine<-function(code_dir)e
out<-file.path(root,"recovered")
plan<-mixedgp_recover_publication(getwd(),archive,data_root,out,"plan",studies="study2",competitor_workers=1L)
stopifnot(fit_calls==0L,competitor_calls==0L,sum(!plan$complete)==4L)
first<-mixedgp_recover_publication(getwd(),archive,data_root,out,"run",studies="study2",competitor_workers=1L)
stopifnot(all(first$complete),fit_calls==1L,competitor_calls==3L)
second<-mixedgp_recover_publication(getwd(),archive,data_root,out,"run",studies="study2",competitor_workers=1L)
stopifnot(all(second$complete),fit_calls==1L,competitor_calls==3L)
state<-readRDS(file.path(out,"state.rds"))
stopifnot(nrow(state$study2$metrics)==150L,nrow(state$study2$raw$predictive_metrics)==50L)
# Existing successful metrics retain their exact values and row IDs.
key<-function(d)paste(d$cell_id,d$rep,d$method)
z<-state$study2$metrics[match(key(comp),key(state$study2$metrics)),names(comp)]
stopifnot(isTRUE(all.equal(z,comp,check.attributes=FALSE)))
# Original archive/input collections are unchanged, and mismatched test data stop
# before any fit: changing a retained response in a recovered state is detected.
state$study2$predictions$y[1]<-state$study2$predictions$y[1]+1
saveRDS(state,file.path(out,"state.rds"))
err<-try(mixedgp_recover_publication(getwd(),archive,data_root,out,"check",studies="study2",competitor_workers=1L),silent=TRUE)
stopifnot(inherits(err,"try-error"),grepl("do not match",err),fit_calls==1L,competitor_calls==3L)
message("Full recovery orchestration passed: plan/no fits, original IDs, only missing tasks, resume, retained scores, input mismatch refusal.")

# Certification requires every Table 2 row and explicit provenance validation.
fullcfg<-study2_simulation_config("publication",code_dir="codes",core_budget=1L)
all_eiv<-all_comp<-list()
for(cell in fullcfg$cells){
 for(cal in cell$calibration_grid){z<-base;z$cell_id<-cell$id;z$n_calib<-cal;all_eiv[[length(all_eiv)+1L]]<-z}
 for(method in c("UC-GP","LVGP","EzGP")){z<-base;z$cell_id<-cell$id;z$method<-method;z$n_calib<-NA_integer_;all_comp[[length(all_comp)+1L]]<-z}
}
state2<-list(study2=list(raw=list(predictive_metrics=do.call(rbind,all_eiv)),metrics=do.call(rbind,all_comp)))
check_out<-file.path(root,"certify");dir.create(check_out)
a<-mixedgp_recovery_tables(state2,list(study2=fullcfg),check_out)
stopifnot(all(a$complete),grepl("INCOMPLETE",readLines(file.path(check_out,"table2_complete.tex"))))
a<-mixedgp_recovery_tables(state2,list(study2=fullcfg),check_out,certify=TRUE)
tex<-readLines(file.path(check_out,"table2_complete.tex"))
stopifnot(any(grepl("\\begin{table}[H]",tex,fixed=TRUE)),any(grepl("same 50 datasets",tex,fixed=TRUE)),
 !any(grepl("Method & R &",tex,fixed=TRUE)))
state2$study2$metrics$CRPS[1]<-NA_real_
a<-mixedgp_recovery_tables(state2,list(study2=fullcfg),check_out,certify=TRUE)
stopifnot(grepl("INCOMPLETE",readLines(file.path(check_out,"table2_complete.tex"))))
message("Table certification and stale-complete-table invalidation passed.")

# Real forked orchestration: distinct processes, identical seeded metrics, resume.
if(.Platform$OS.type!="windows"){
  calls_dir<-file.path(root,"parallel-calls");dir.create(calls_dir)
  serial_mock<-e$mixedgp_cached_competitors
  e$mixedgp_cached_competitors<-function(...){
    args<-list(...)
    writeLines(as.character(Sys.getpid()),file.path(calls_dir,paste0(args$methods,"-",args$seed)))
    do.call(serial_mock,args)
  }
  parallel_out<-file.path(root,"parallel-recovered")
  expected_chain_workers<-2L
  par_audit<-mixedgp_recover_publication(getwd(),archive,data_root,parallel_out,"run",studies="study2",competitor_workers=8L,core_budget=2L)
  call_files<-list.files(calls_dir,full.names=TRUE)
  pids<-vapply(call_files,function(f)readLines(f)[1],character(1))
  stopifnot(all(par_audit$complete),length(call_files)==3L,length(unique(pids))>=2L,sum(pids!=as.character(Sys.getpid()))>=2L)
  parallel_state<-readRDS(file.path(parallel_out,"state.rds"))
  stopifnot(identical(parallel_state$study2$metrics,state$study2$metrics))
  before<-tools::md5sum(call_files)
  mixedgp_recover_publication(getwd(),archive,data_root,parallel_out,"run",studies="study2",competitor_workers=2L)
  stopifnot(identical(before,tools::md5sum(call_files)))
  message("Multicore orchestration passed: distinct worker PIDs, serial-equivalent scores, no refitting on resume.")
}

# CPU allocation is validated and caps both phases, without increasing defaults.
stopifnot(identical(mixedgp_recovery_core_settings(2L,4L,8L),
  list(core_budget=2L,chain_workers=2L,competitor_workers=2L)),
  mixedgp_recovery_core_settings(8L)$competitor_workers==4L,
  mixedgp_recovery_core_settings(8L,competitor_workers=8L)$competitor_workers==8L)
for(bad in list(0,-1,1.5,NA_real_,Inf,c(2L,4L)))
  stopifnot(inherits(try(mixedgp_recovery_core_settings(bad),silent=TRUE),"try-error"))
message("CPU budget validation and phase caps passed.")

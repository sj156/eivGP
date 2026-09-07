source("codes/simulation_helpers.R")
source("codes/combine_experiment_results.R")
cells <- list(list(id="a",n_rep=3L))
own <- data.frame(cell_id="a",rep=1:2,method="EIV-GP",n_calib=20,RMSE=c(1,2))
comp <- data.frame(cell_id="a",rep=1:2,method="UC-GP",n_calib=NA_integer_,RMSE=c(2,4))
status <- data.frame(cell_id="a",rep=1:2,method="UC-GP",status="success",dataset_md5=c("x","y"))
ids <- status[c("cell_id","rep","dataset_md5")]
ans <- mixedgp_combine_predictions(own,comp,status,ids,cells)
stopifnot(nrow(ans$summary)==2,all(ans$summary$n_expected==3),
  ans$paired$n_pairs==2,ans$paired$mean_eiv_minus_competitor == -1.5)
bad <- status; bad$dataset_md5[1]<-"wrong"
stopifnot(inherits(tryCatch(mixedgp_combine_predictions(own,comp,bad,ids,cells),error=identity),"error"))
stopifnot(inherits(tryCatch(mixedgp_combine_predictions(own,rbind(comp,comp),status,ids,cells),error=identity),"error"))
missing <- mixedgp_combine_predictions(own,data.frame(),data.frame(),ids,cells)
stopifnot(nrow(missing$summary)==1,nrow(missing$paired)==0,!any(missing$availability$competitor_metrics_available))
## Portable on-disk fixture for each study; no frozen data/model cache required.
for(study in c("study1","study2")) {
  root<-tempfile(); dir.create(root)
  run<-file.path(root,"mcmc"); cp<-file.path(root,"competitors")
  dir.create(file.path(run,"config"),recursive=TRUE); dir.create(cp)
  dir.create(file.path(run,"cells","a"),recursive=TRUE)
  config<-list(study=study,cells=cells)
  saveRDS(config,file.path(run,"config","resolved_config.rds"))
  context<-list(mc_results=own)
  saveRDS(list(study=study,context=context),file.path(run,"cells","a","report_inputs.rds"))
  write.csv(ids,file.path(run,"cells","a","dataset_identity.csv"),row.names=FALSE)
  saveRDS(list(config=config),file.path(cp,"provenance.rds"))
  write.csv(comp,file.path(cp,"metrics.csv"),row.names=FALSE)
  write.csv(status,file.path(cp,"statuses.csv"),row.names=FALSE)
  out<-file.path(root,"paper")
  mixedgp_combine_saved_runs(run,cp,out)
  stopifnot(file.exists(file.path(out,"predictive_summary.tex")))
}
message("Independent reporting: both studies, missing methods, matched pairs, duplicate and checksum rejection passed.")

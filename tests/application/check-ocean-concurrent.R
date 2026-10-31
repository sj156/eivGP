# Author: Xinye Chen. Independent Ocean processes and saved-result equivalence.
# Reference roots contain five-fold Ocean200 pilot results.
script<-sub('^--file=','',grep('^--file=',commandArgs(FALSE),value=TRUE)[1])
repo<-normalizePath(file.path(dirname(normalizePath(script)),'../..'))
source(file.path(repo,'tools/ocean-tools.R'))
p<-list()
for(arg in commandArgs(TRUE)) {
  if(!grepl('^--(output|reference-eiv|reference-competitors|cores)=.+$',arg))stop('Use --output=NEW --reference-eiv=ROOT --reference-competitors=ROOT --cores=4')
  key<-sub('^--([^=]+)=.*$','\\1',arg);if(key%in%names(p))stop('Repeated option: ',key)
  p[[key]]<-sub('^[^=]+=','',arg)
}
stopifnot(length(p[['reference-eiv']])>0L,length(p[['reference-competitors']])>0L)
reference<-lapply(p[c('reference-eiv','reference-competitors')],ocean_resolve_path,repository_dir=repo)
out<-ocean_resolve_path(if(is.null(p$output))'results/ocean_concurrent_check'else p$output,repo)
stopifnot(!dir.exists(out));dir.create(out,recursive=TRUE)
workers<-if(is.null(p$cores))4L else as.integer(p$cores);stopifnot(workers>=4L,.Platform$OS.type=='unix')
Sys.setenv(R_LIBS=paste(.libPaths(),collapse=.Platform$path.sep))
result<-file.path(out,'run')
common<-c('--profile=pilot','--dataset=ocean200',paste0('--cores=',workers),paste0('--output=',result))
run<-function(group,extra=character(),label=group)system2(file.path(R.home('bin'),'Rscript'),
  shQuote(c(file.path(repo,'tools',paste0('ocean-',group,'.R')),common,extra)),
  stdout=file.path(out,paste0(label,'.log')),stderr=file.path(out,paste0(label,'.log')))
jobs<-setNames(lapply(c('competitors','eivGP'),function(group)
  parallel::mcparallel(run(group),mc.set.seed=FALSE)),c('competitors','eivGP'))
both_owned<-FALSE;finished<-list();start<-proc.time()[['elapsed']]
while(length(jobs)) {
  if(all(file.exists(file.path(result,'_meta/locks',paste0(c('eivGP','competitors'),'.lock/owner.json')))))both_owned<-TRUE
  got<-parallel::mccollect(jobs,wait=FALSE)
  if(length(got)){finished<-c(finished,got);jobs<-Filter(function(x)!as.character(x$pid)%in%names(got),jobs)}
  if(proc.time()[['elapsed']]-start>1800)stop('Concurrent check timed out; stop its workers and inspect ',out)
  if(length(jobs))Sys.sleep(.1)
}
stopifnot(both_owned,all(vapply(finished,function(x)identical(x,0L),logical(1))))
for(group in c('eivGP','competitors'))
  stopifnot(readRDS(ocean_path(result,paste0('execution_settings_',group,'_fit.rds')))$cores==workers)
compared<-character()
for(group in names(reference)) {
  root<-reference[[group]]
  files<-list.files(root,'[.]csv$',recursive=TRUE,full.names=TRUE)
  relative<-substring(files,nchar(root)+2L)
  keep<-if(group=='reference-eiv')grepl('^eiv_fold[0-9]+_O[0-9]+/[^/]+[.]csv$',relative)else
    grepl('^(baseline_fold[0-9]+|cc_fold[0-9]+_O[0-9]+_[^/]+)/[^/]+[.]csv$',relative)
  keep<-keep&grepl('/(predictive_metrics|metrics|predictions|mixed_predictions|baseline_metrics|baseline_predictions|U_posterior_intervals|imputations|imputations_model_scale|convergence_parameters|convergence_summary)[.]csv$',relative)
  stopifnot(any(keep))
  for(i in which(keep)) {
    stopifnot(isTRUE(all.equal(read.csv(files[i]),read.csv(file.path(result,relative[i])),tolerance=1e-12)))
    compared<-c(compared,relative[i])
  }
}
before<-tools::md5sum(list.files(result,'^(fit|package_output|baseline_package_outputs)[.]rds$',recursive=TRUE,full.names=TRUE))
stopifnot(length(before)==50L)
for(group in c('eivGP','competitors')) {
  stopifnot(file.exists(file.path(result,'reports',group,'index.html')),
    run(group,'--reuse-only=true',paste0(group,'_reuse'))==0L)
}
stopifnot(identical(before,tools::md5sum(names(before))),run('report')==0L,
  file.exists(file.path(result,'ocean.html')),
  !any(dir.exists(file.path(result,'_meta/locks',paste0(c('eivGP','competitors'),'.lock')))))
ocean_atomic_json(list(status='PASS',profile='pilot',dataset='ocean200',folds=1:5,
  simultaneously_owned_group_locks=both_owned,cores_per_process=workers,
  metric_and_prediction_tables_compared=length(compared),comparison_tolerance=1e-12,
  fitted_object_hashes_unchanged=TRUE,combined_report=TRUE),file.path(out,'CHECK.json'))
cat('PASS: concurrent Ocean CLIs, five-fold metric equivalence, unchanged saved fits and combined HTML.\n')

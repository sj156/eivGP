# Author: Xinye Chen. Two real OAI CLI processes share one synthetic result root.
script<-sub('^--file=','',grep('^--file=',commandArgs(FALSE),value=TRUE)[1])
repo<-normalizePath(file.path(dirname(normalizePath(script)),'../..'))
source(file.path(repo,'tools/oai.R'));api<-oai_load_api(repo)
p<-list()
for(arg in commandArgs(TRUE)) {
  if(!grepl('^--(output|reference|cores)=.+$',arg))stop('Options: --output=NEW --reference=COMPLETED_SYNTHETIC_CHECK --cores=4')
  key<-sub('^--([^=]+)=.*$','\\1',arg)
  if(key%in%names(p))stop('Repeated option: ',key)
  p[[key]]<-sub('^[^=]+=','',arg)
}
if(is.null(p$reference))stop('Supply the completed output root from tests/application/check-oai.R as --reference.')
reference<-api$oai_resolve_path(p$reference,repo)
check<-jsonlite::fromJSON(file.path(reference,'CHECK.json'));stopifnot(isTRUE(check$synthetic_only),check$status=='PASS')
out<-api$oai_resolve_path(if(is.null(p$output))'results/oai_concurrent_check'else p$output,repo)
stopifnot(!dir.exists(out));dir.create(out,recursive=TRUE)
workers<-if(is.null(p$cores))4L else api$oai_integer(p$cores,'cores')
stopifnot(workers>=4L,.Platform$OS.type=='unix')
Sys.setenv(R_LIBS=paste(.libPaths(),collapse=.Platform$path.sep))
result<-file.path(out,'run');input<-file.path(reference,'inputs')
common<-c('--profile=paper',paste0('--data-root=',input),paste0('--output=',result),
  paste0('--cores=',workers),'--folds=1','--cells=O0,O20,O50',
  '--iter=40','--burn=10','--stages=40','--first-checkpoint=40','--checkpoint-step=40',
  '--draws-per-chain=5','--baseline-profile=quick','--functional-draws-per-chain=4','--n-latent=16','--appendix=false')
run<-function(group,extra=character(),label=group) {
  log<-file.path(out,paste0(label,'.log'))
  system2(file.path(R.home('bin'),'Rscript'),
    shQuote(c(file.path(repo,'tools',paste0('oai-',group,'.R')),common,extra)),stdout=log,stderr=log)
}
jobs<-setNames(lapply(c('eivGP','competitors'),function(group)
  parallel::mcparallel(run(group),mc.set.seed=FALSE)),c('eivGP','competitors'))
both_owned<-FALSE;finished<-list();start<-proc.time()[['elapsed']]
while(length(jobs)) {
  locks<-file.path(result,'_meta/locks',paste0(c('eivGP','competitors'),'.lock/owner.json'))
  if(all(file.exists(locks)))both_owned<-TRUE
  got<-parallel::mccollect(jobs,wait=FALSE)
  if(length(got)) {
    finished<-c(finished,got)
    jobs<-Filter(function(x)!as.character(x$pid)%in%names(got),jobs)
  }
  if(proc.time()[['elapsed']]-start>900)stop('Concurrent CLI check timed out; inspect the test output and stop its workers.')
  if(length(jobs))Sys.sleep(.1)
}
stopifnot(both_owned,all(vapply(finished,function(x)identical(x,0L),logical(1))))
for(group in c('eivGP','competitors')) {
  text<-readLines(file.path(out,paste0(group,'.log')))
  stopifnot(any(grepl(paste0('; ',workers,' cores;'),text,fixed=TRUE)))
  stopifnot(any(grepl('START fit',text,fixed=TRUE)),any(grepl('DONE fit',text,fixed=TRUE)),
    any(grepl('COMPLETE | report:',text,fixed=TRUE)),!any(grepl('ADNI',text,fixed=TRUE)),
    file.exists(file.path(result,paste0('PROGRESS_',group,'.txt'))),
    file.exists(file.path(result,'reports',group,'index.html')))
}
expected<-file.path(reference,'run')
metrics<-list.files(expected,'^(prediction_metrics|metrics)[.]csv$',recursive=TRUE,full.names=TRUE)
metrics<-metrics[grepl('/(eiv_fold|baseline_fold|cc_fold)',metrics)]
stopifnot(length(metrics)==12L)
for(path in metrics) {
  relative<-substring(path,nchar(expected)+2L)
  stopifnot(isTRUE(all.equal(read.csv(path),read.csv(file.path(result,relative)),tolerance=1e-12)))
}
before<-tools::md5sum(list.files(result,'^fit[.]rds$',recursive=TRUE,full.names=TRUE))
stopifnot(run('eivGP',label='eivGP_reuse')==0L,run('competitors',label='competitors_reuse')==0L)
stopifnot(identical(tools::md5sum(names(before)),before),run('report')==0L)
api$atomic_json(list(schema='oai-concurrent-cli-check-v1',status='PASS',synthetic_only=TRUE,
  simultaneously_owned_group_locks=both_owned,independent_processes=2L,cores_per_process=workers,
  comparison_with_sequential_run='12 metric tables agree to 1e-12',fit_reuse_unchanged=TRUE,
  terminal_start_done_and_complete=TRUE,combined_report=TRUE),file.path(out,'CHECK.json'))
cat('PASS: two concurrent CLI processes, sequential metric equivalence, terminal progress, fit reuse and combined HTML.\n')

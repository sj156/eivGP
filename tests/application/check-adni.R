# Author: Xinye Chen. Portable synthetic check of the three ADNI CLI entrypoints.
script<-sub('^--file=','',grep('^--file=',commandArgs(FALSE),value=TRUE)[1])
repo<-normalizePath(file.path(dirname(normalizePath(script)),'../..'))
source(file.path(repo,'tools/adni.R'))
source(file.path(repo,'tools/adni-tools.R'))
params<-adni_cli_params(commandArgs(TRUE))
if(any(!names(params)%in%c('cores','output')))stop('Check options: --cores=N --output=results/NEW_DIRECTORY')
cores<-adni_integer(if(is.null(params$cores))2L else params$cores,'cores')
api<-adni_load_api(repo)
out<-adni_resolve_path(if(is.null(params$output))file.path('results',paste0('adni_quick_check_',format(Sys.time(),'%Y%m%d_%H%M%S')))else params$output,repo)
api$application$adni_output_check(out,repo)
if(dir.exists(out))stop('Use a new check directory: ',out)
dir.create(out,recursive=TRUE)
set.seed(20261023L);n<-40L
d<-data.frame(RID=seq_len(n),X1=60+runif(n)*20,X2=10+runif(n)*12,U1=1.2+rnorm(n,sd=.15),U2=4.5+rnorm(n,sd=.6))
d$Y<-1.5*d$U1+.2*d$U2+.02*(d$X2-16)+rnorm(n,sd=.4)
for(j in 1:3)d[[paste0('C',j)]]<-sample(rep(1:4,length.out=n))
meta<-list(dataset_id='synthetic_check',population='Synthetic workflow check; no ADNI participant records',
  Y='synthetic_response',X_order=c('age','education'),U_order=c('FDG','HIP1000_ICV'),
  C_order=c('MEMORY3','MEMORY8','DIVATT3'),m_vec=rep(4L,3),fold_seed=9L,calibration_seed=10L,test_seed=11L,model_seed=12L)
write.csv(d,file.path(out,'synthetic.csv'),row.names=FALSE)
jsonlite::write_json(meta,file.path(out,'synthetic.json'),auto_unbox=TRUE,pretty=TRUE)
input<-file.path(out,'inputs');result<-file.path(out,'run')
adni_prepare_folds(file.path(out,'synthetic.csv'),file.path(out,'synthetic.json'),input)
relative<-function(p)if(startsWith(p,paste0(repo,'/')))substring(p,nchar(repo)+2L)else p
Sys.setenv(R_LIBS=paste(.libPaths(),collapse=.Platform$path.sep))
common<-c('--profile=quick',paste0('--data-root=',relative(input)),paste0('--output=',relative(result)),
  '--folds=1','--cells=O50','--functional-draws-per-chain=8','--n-latent=16','--appendix=false')
run_cli<-function(group,label,workers=cores,extra=character()) {
  log<-file.path(out,paste0(label,'.log'))
  old<-getwd();setwd(tempdir());on.exit(setwd(old))
  code<-system2(file.path(R.home('bin'),'Rscript'),
    shQuote(c(file.path(repo,'tools',paste0('adni-',group,'.R')),common,paste0('--cores=',workers),extra)),
    stdout=log,stderr=log)
  if(code!=0L)stop('CLI failed: ',label,'; inspect ',log)
}
run_cli('eivGP','01_eivGP');run_cli('competitors','02_competitors');run_cli('report','03_report')
fit<-api$adni_plan_fit_path(result,list(fold=1L,cell='O50'));fit_before<-tools::md5sum(fit)
scores<-api$adni_path(result,'S80/analysis/fold01/O50/prediction_metrics.csv')
before<-tools::md5sum(scores)
combined<-file.path(result,'adni.html');combined_before<-tools::md5sum(combined)
run_cli('eivGP','04_core_only_reuse',if(cores>=4L)4L else cores)
stopifnot(identical(tools::md5sum(combined),combined_before))
run_cli('eivGP','05_saved_only',1L,'--reuse-only=true')
stopifnot(identical(tools::md5sum(combined),combined_before))
stopifnot(identical(tools::md5sum(fit),fit_before),identical(tools::md5sum(scores),before))
for(method in c('UC-GP','LVGP','EzGP'))for(cell in c('OI','O50')) {
  path<-if(cell=='OI')api$adni_path(result,'baselines/fold01',method)else api$adni_path(result,'CC/O50/baselines/fold01',method)
  score<-read.csv(file.path(path,'metrics.csv'));stopifnot(all(is.finite(score$RMSE_SD)))
}
for(p in c(file.path(result,'adni.html'),file.path(result,'reports/eivGP/index.html'),file.path(result,'reports/competitors/index.html')))
  stopifnot(file.exists(p),file.info(p)$size>1000)
stopifnot(api$adni_storage_flat(result),
  file.exists(file.path(result,'reports/paper/artifact_manifest.csv')),
  file.exists(file.path(result,'data/fold01_O50_train.csv')),
  !any(dir.exists(file.path(result,c('fits','S80','CC','workflow','receipts','report-assets')))))
api$atomic_json(list(schema='adni-cli-check-v1',synthetic_only=TRUE,status='PASS',
  entrypoints=c('eivGP','competitors','report'),cores=cores,core_only_reuse=TRUE,
  relative_paths_tested=!identical(relative(input),input),caller_outside_repository=TRUE,
  native_comparator_scores=6L,output=relative(out),convergence_required=FALSE),file.path(out,'CHECK.json'))
cat('PASS: three ADNI CLIs, six comparator scores, saved HTML and unchanged fit after core-only reuse.\n',
  'Check record: ',file.path(out,'CHECK.json'),'\n',sep='')

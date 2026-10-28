# Author: Xinye Chen. Synthetic execution check; no OAI participant data.
script<-sub('^--file=','',grep('^--file=',commandArgs(FALSE),value=TRUE)[1])
repo<-normalizePath(file.path(dirname(normalizePath(script)),'../..'))
source(file.path(repo,'tools/oai.R'))
api<-oai_load_api(repo)
params<-api$adni_cli_params(commandArgs(TRUE))
if(any(!names(params)%in%c('cores','output')))stop('Check options: --cores=N --output=results/NEW_DIRECTORY')
cores<-api$adni_integer(if(is.null(params$cores))4L else params$cores,'cores')
out<-api$adni_resolve_path(if(is.null(params$output))file.path('results',paste0('oai_quick_check_',format(Sys.time(),'%Y%m%d_%H%M%S')))else params$output,repo)
api$application$adni_output_check(out,repo)
if(dir.exists(out))stop('Use a new check directory: ',out)
dir.create(out,recursive=TRUE)
set.seed(20261028L);n<-60L
d<-data.frame(RID=seq_len(n),X1=50+runif(n)*25,X2=20+runif(n)*12,X3=1500+runif(n)*200)
d$U1<-3+.4*rnorm(n)
d$Y<-.8*d$U1+.01*(d$X1-60)+rnorm(n,sd=.2)
grade<-rep(1:3,length.out=n)
for(j in 1:3)d[[paste0('C',j)]]<-grade
meta<-list(dataset_id='oai',population='Synthetic workflow check; no participant records',
  Y='synthetic_response',Y_label='Synthetic joint-space response',X_order=c('age','BMI','height'),
  U_order='cartilage',U_labels=list(cartilage='Synthetic cartilage reference'),
  C_order=paste0('grade',1:3),C_labels=as.list(setNames(paste('Synthetic native grade',1:3),paste0('grade',1:3))),
  grade_encoding=setNames(rep(list(as.list(c('0 -> 1','1 -> 2','2 -> 3'))),3),paste0('C',1:3)),
  m_vec=rep(3L,3),stratification_columns=c('C1','C2'),
  fold_seed=9L,calibration_seed=10L,test_seed=11L,model_seed=12L,
  imputation_display_fold=4L,test_U_selection=list(outcome_selected=FALSE))
write.csv(d,file.path(out,'synthetic.csv'),row.names=FALSE)
jsonlite::write_json(meta,file.path(out,'synthetic.json'),auto_unbox=TRUE,pretty=TRUE)
input<-file.path(out,'inputs');result<-file.path(out,'run')
api$oai_prepare_folds(file.path(out,'synthetic.csv'),file.path(out,'synthetic.json'),input)
relative<-function(p)if(startsWith(p,paste0(repo,'/')))substring(p,nchar(repo)+2L)else p
Sys.setenv(R_LIBS=paste(.libPaths(),collapse=.Platform$path.sep))
common<-c('--profile=quick',paste0('--data-root=',relative(input)),paste0('--output=',relative(result)),
  '--folds=1','--cells=O0,O20,O50','--iter=40','--burn=10','--stages=40',
  '--first-checkpoint=40','--checkpoint-step=40','--draws-per-chain=5',
  '--functional-draws-per-chain=4','--n-latent=16','--appendix=false')
run_cli<-function(group,label,workers=cores,extra=character()) {
  log<-file.path(out,paste0(label,'.log'))
  old<-getwd();setwd(tempdir());on.exit(setwd(old))
  code<-system2(file.path(R.home('bin'),'Rscript'),
    shQuote(c(file.path(repo,'tools',paste0('oai-',group,'.R')),common,paste0('--cores=',workers),extra)),
    stdout=log,stderr=log)
  if(code!=0L)stop('CLI failed: ',label,'; inspect ',log)
}
run_cli('eivGP','01_eivGP');run_cli('competitors','02_competitors');run_cli('report','03_report')
fits<-vapply(c('O0','O20','O50'),function(cell)api$adni_plan_fit_path(result,list(fold=1L,cell=cell)),character(1))
fit_before<-tools::md5sum(fits)
scores<-api$adni_path(result,'S40/analysis/fold01/O50/prediction_metrics.csv');before<-tools::md5sum(scores)
combined<-file.path(result,'oai.html');combined_before<-tools::md5sum(combined)
run_cli('eivGP','04_core_only_reuse',if(cores>=4L)4L else cores)
stopifnot(identical(tools::md5sum(combined),combined_before))
run_cli('eivGP','05_saved_only',1L,'--reuse-only=true')
stopifnot(identical(tools::md5sum(combined),combined_before),
  identical(tools::md5sum(fits),fit_before),identical(tools::md5sum(scores),before))
for(method in c('UC-GP','LVGP','EzGP'))for(cell in c('OI','O20','O50')) {
  path<-if(cell=='OI')api$adni_path(result,'baselines/fold01',method)else api$adni_path(result,paste0('CC/',cell,'/baselines/fold01'),method)
  score<-read.csv(file.path(path,'metrics.csv'));stopifnot(all(is.finite(score$RMSE_SD)))
}
for(p in c(combined,file.path(result,'reports/eivGP/index.html'),file.path(result,'reports/competitors/index.html')))
  stopifnot(file.exists(p),file.info(p)$size>1000)
stopifnot(api$adni_storage_flat(result),file.exists(file.path(result,'reports/paper/artifact_manifest.csv')),
  file.exists(file.path(result,'data/fold01_O50_train.csv')),
  !any(dir.exists(file.path(result,c('fits','S40','CC','workflow','receipts','report-assets')))))
api$atomic_json(list(schema='oai-cli-check-v1',synthetic_only=TRUE,status='PASS',
  entrypoints=c('eivGP','competitors','report'),cores=cores,core_only_reuse=TRUE,
  relative_paths_tested=!identical(relative(input),input),caller_outside_repository=TRUE,
  EIV_cells=c('O0','O20','O50'),native_comparator_scores=9L,output=relative(out),convergence_required=FALSE),file.path(out,'CHECK.json'))
cat('PASS: three OAI CLIs, three EIV cells, nine comparator scores, HTML and saved-fit reuse.\n',
  'Check record: ',file.path(out,'CHECK.json'),'\n',sep='')

# Author: Xinye Chen. Exact sampler recovery, including an interrupted warm-up.
repo<-Sys.getenv('EIVGP_TEST_REPO');stopifnot(nzchar(repo))
source(file.path(repo,'tools/oai.R'))

testthat::test_that('HMC and ESS recover exact chains after warm-up and full-fit interruptions',{
  root<-tempfile();dir.create(root);on.exit(unlink(root,recursive=TRUE))
  e<-oai_load_api(repo);set.seed(71L);n<-30L
  d<-data.frame(RID=1:n,Y=rnorm(n),X1=runif(n),X2=runif(n),X3=runif(n),U1=rnorm(n))
  for(j in 1:8)d[[paste0('C',j)]]<-rep(1:3,length.out=n)
  meta<-list(dataset_id='oai',Y='synthetic',X_order=c('age','BMI','height'),U_order='cartilage',
    C_order=paste0('grade',1:8),m_vec=rep(3L,8),fold_seed=9L,calibration_seed=10L,test_seed=11L,model_seed=12L)
  write.csv(d,file.path(root,'cohort.csv'),row.names=FALSE)
  jsonlite::write_json(meta,file.path(root,'meta.json'),auto_unbox=TRUE)
  e$oai_prepare_folds(file.path(root,'cohort.csv'),file.path(root,'meta.json'),file.path(root,'inputs'))
  input<-e$oai_plan_inputs(file.path(root,'inputs'))
  s<-e$oai_screen_settings(list(profile='quick',smoke=TRUE,chains=4L,burn=8L,stages=32L,
    first_checkpoint=16L,checkpoint_step=8L,folds=1L,cells=c('O20','O50'),trial_id='oai'),32L,1L)
  s$final_stage<-FALSE
  real_write<-getFromNamespace('mixedgp_checkpoint_write','eivGP')
  for(cell in c('O20','O50')) {
    task<-list(fold=1L,cell=cell);out<-file.path(root,cell)
    e$oai_storage_init(out,32L)
    # Fail immediately after publishing a genuine chain checkpoint at burn.
    testthat::expect_error(suppressWarnings(testthat::with_mocked_bindings(
      e$oai_plan_fit_task(task,input,out,s),
      mixedgp_checkpoint_write=function(object,path) {
        real_write(object,path)
        if(grepl('chain-[0-9]+[.]rds$',path)&&object$iteration==8L)stop('test warm-up interruption')
        invisible(path)
      },.package='eivGP')))
    path<-e$oai_plan_fit_path(out,task)
    testthat::expect_false(file.exists(path))
    warm<-list.files(dirname(path),'^chain-[0-9]+[.]rds$',recursive=TRUE,full.names=TRUE)
    testthat::expect_true(length(warm)>=1L)
    testthat::expect_identical(readRDS(warm[1])$iteration,8L)
    save_fit<-e$oai_plan_save_fit
    e$oai_plan_save_fit<-function(fit,path) {
      save_fit(fit,path)
      if(fit$checkpoint$iteration==24L)stop('test published-fit interruption')
      invisible(path)
    }
    testthat::expect_error(suppressWarnings(e$oai_plan_fit_task(task,input,out,s)),'published-fit interruption')
    testthat::expect_identical(e$oai_plan_read_fit(path)$checkpoint$iteration,24L)
    e$oai_plan_save_fit<-save_fit;s$cores<-2L
    suppressWarnings(e$oai_plan_fit_task(task,input,out,s))
    resumed<-e$oai_plan_read_fit(path)
    direct<-suppressWarnings(do.call(eivGP::fit_eivgp,c(e$oai_plan_fit_args(input$cells[[paste(1,cell,sep='/')]],s),
      list(n_iter=32L,burn=8L))))
    testthat::expect_identical(resumed$mcmc$samples_by_chain,direct$mcmc$samples_by_chain)
    testthat::expect_identical(resumed$checkpoint$states,direct$checkpoint$states)
    testthat::expect_true(!is.null(resumed$mcmc$samples_by_chain$S))
    z<-input$cells[[paste(1,cell,sep='/')]]
    predict<-function(fit) {
      C<-as.matrix(z$test[,z$c]);colnames(C)<-z$protocol$C_order
      eivGP::predict_eivgp(fit,new_X=as.matrix(z$test[,z$x]),new_C=C,
        target='response',draw_ids=1:4,seed=z$protocol$seeds$prediction)
    }
    testthat::expect_identical(predict(resumed),predict(direct))
    progress<-read.csv(file.path(dirname(path),'mcmc_progress.csv'))
    testthat::expect_identical(progress$iteration,32L)
    testthat::expect_identical(progress$retained_per_chain,24L)
    testthat::expect_identical(progress$event,'FIT_COMPLETE')
    testthat::expect_length(list.files(dirname(path),'^chain-[0-9]+[.]rds$',recursive=TRUE),0L)
    s$cores<-1L
  }
  e$oai_render_report<-function(...)invisible(NULL)
  args<-list(profile='quick',cores=4L,data_root=file.path(root,'inputs'),output=file.path(root,'split'),
    iter=24L,burn=8L,first_checkpoint=16L,checkpoint_step=8L,folds='1:2',cells='O20',
    draws_per_chain=5L,functional_draws_per_chain=4L,n_latent=8L,appendix=FALSE,run_folds='1')
  suppressWarnings(e$oai_method_run(repo,args,'eivGP'))
  identity<-tools::md5sum(e$oai_path(args$output,'workflow_identity.rds'))
  args$run_folds<-'2';suppressWarnings(e$oai_method_run(repo,args,'eivGP'))
  testthat::expect_identical(tools::md5sum(e$oai_path(args$output,'workflow_identity.rds')),identity)
  split<-args$output;args$output<-file.path(root,'batch');args$run_folds<-NULL
  suppressWarnings(e$oai_method_run(repo,args,'eivGP'))
  for(fold in 1:2) {
    task<-list(fold=fold,cell='O20')
    one<-e$oai_plan_read_fit(e$oai_plan_fit_path(split,task))
    batch<-e$oai_plan_read_fit(e$oai_plan_fit_path(args$output,task))
    testthat::expect_identical(one$mcmc$samples_by_chain,batch$mcmc$samples_by_chain)
    testthat::expect_identical(one$checkpoint$states,batch$checkpoint$states)
    metric<-function(out)e$oai_path(out,'S24/analysis',sprintf('fold%02d',fold),'O20/prediction_metrics.csv')
    testthat::expect_identical(read.csv(metric(split)),read.csv(metric(args$output)))
  }
})

testthat::test_that('execution-only fold selection leaves the full frozen manifest and identity unchanged',{
  root<-tempfile();dir.create(root);on.exit(unlink(root,recursive=TRUE))
  e<-oai_load_api(repo)
  input<-file.path(repo,'data/oai')
  testthat::skip_if_not(file.exists(file.path(input,'DATASET.json')),'Private fixture not supplied')
  out<-file.path(root,'result');seen<-list()
  e$oai_screen_run<-function(repo,p,group){seen[[length(seen)+1L]]<<-p;invisible(NULL)}
  e$oai_render_report<-function(...)invisible(NULL)
  args<-list(profile='paper',cores=4L,data_root=input,output=out,run_folds='1')
  e$oai_method_run(repo,args,'eivGP')
  before<-tools::md5sum(e$oai_path(out,'application_manifest.json'))
  args$run_folds<-'2';e$oai_method_run(repo,args,'eivGP')
  args$run_folds<-NULL;e$oai_method_run(repo,args,'eivGP')
  testthat::expect_identical(tools::md5sum(e$oai_path(out,'application_manifest.json')),before)
  m<-jsonlite::fromJSON(e$oai_path(out,'application_manifest.json'))
  testthat::expect_identical(m$folds,1:5)
  testthat::expect_identical(m$checkpoint_step,5000L)
  testthat::expect_identical(m$first_checkpoint,10000L)
  testthat::expect_identical(vapply(seen[1:2],function(x)x$run_folds,1L),c(1L,2L))
  testthat::expect_null(seen[[3]]$run_folds)
})

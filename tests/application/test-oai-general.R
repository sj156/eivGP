# Author: Xinye Chen. OAI dimensions, information boundaries and privacy contracts.
repo<-Sys.getenv('EIVGP_TEST_REPO');stopifnot(nzchar(repo))
source(file.path(repo,'tools/oai.R'))
testthat::test_that('OAI has exactly the same six flat entry/tool files as ADNI',{
  expected<-paste0('oai',c('.R','-tools.R','-eivGP.R','-competitors.R','-report.R','-prepare-folds.R'))
  testthat::expect_setequal(basename(list.files(file.path(repo,'tools'),'^oai.*[.]R$',full.names=TRUE)),expected)
  api<-oai_load_api(repo)
  testthat::expect_true(all(file.exists(api$oai_source_files(repo))))
  testthat::expect_true(all(expected%in%basename(api$oai_source_files(repo))))
  testthat::expect_false('adni-tools.R'%in%basename(api$oai_source_files(repo)))
  testthat::expect_false(any(grepl('^adni_',ls(api))))
  testthat::expect_false(any(grepl('^(oai_|OAI_)',ls(asNamespace('eivGP'),all.names=TRUE))))
  testthat::expect_equal(api$oai_default_input('oai',repo),file.path(repo,'data/oai'))
  testthat::expect_equal(api$oai_default_input('other',repo),file.path(repo,'data/oai_other'))
  testthat::expect_equal(api$oai_default_output('oai','paper',repo),file.path(repo,'results/oai'))
})
testthat::test_that('the OAI API loads in a checkout containing no ADNI tools',{
  root<-tempfile();dir.create(root);on.exit(unlink(root,recursive=TRUE))
  dir.create(file.path(root,'tools'));dir.create(file.path(root,'paper/applications'),recursive=TRUE)
  dir.create(file.path(root,'data/oai'),recursive=TRUE)
  for(path in oai_load_api(repo)$oai_source_files(repo)) {
    relative<-substring(path,nchar(repo)+2L)
    testthat::expect_true(file.copy(path,file.path(root,relative)))
  }
  testthat::expect_true(file.copy(file.path(repo,'data/oai/DATASET.example.json'),file.path(root,'data/oai/DATASET.example.json')))
  env<-new.env(parent=globalenv());sys.source(file.path(root,'tools/oai.R'),envir=env)
  api<-env$oai_load_api(root)
  testthat::expect_identical(api$oai_settings('paper',4L)$iter,50000L)
  testthat::expect_false(file.exists(file.path(root,'tools/adni-tools.R')))
  testthat::expect_true(all(file.exists(api$oai_source_files(root))))
})
testthat::test_that('OAI pilot preserves total1500/burn500 and inherits current paper resources',{
  api<-oai_load_api(repo)
  for(cores in c(1L,2L,4L,8L)) {
    s<-api$oai_settings('pilot',cores)
    testthat::expect_equal(s$chains,min(cores,4L));testthat::expect_equal(s$stages,1500L)
    testthat::expect_equal(s$iter,1500L);testthat::expect_equal(s$burn,500L)
    testthat::expect_equal(s$checkpoint_step,500L);testthat::expect_equal(s$first_checkpoint,750L)
    p<-api$oai_settings('paper',cores)
    testthat::expect_equal(p$chains,min(cores,4L));testthat::expect_equal(p$iter,50000L);testthat::expect_equal(p$burn,5000L)
  }
  s<-api$oai_settings('pilot',cores=4L,iter=2000L,burn=600L)
  testthat::expect_equal(s$stages,2000L);testthat::expect_equal(s$burn,600L)
  testthat::expect_error(api$oai_settings('pilot',cores=2L,chains=4L),'chains must equal')
  s<-api$oai_settings('pilot',cores=4L)
  testthat::expect_equal(c(s$sampler_o0,s$sampler_o20,s$sampler_o50),c('HMC','HMC','ESS'))
})
testthat::test_that('variable C counts use the shared C-only folds and mask hidden U',{
  api<-oai_load_api(repo);root<-tempfile();dir.create(root);on.exit(unlink(root,recursive=TRUE))
  set.seed(991L);n<-90L;grade<-rep(1:3,length.out=n)
  d<-data.frame(RID=seq_len(n),Y=rnorm(n),X1=50+runif(n)*30,X2=20+runif(n)*15,X3=1500+runif(n)*250,U1=rnorm(n))
  for(j in 1:8)d[[paste0('C',j)]]<-grade
  for(q in c(3L,8L)) {
    meta<-list(dataset_id='oai',Y='Synthetic response',X_order=c('age','BMI','height'),U_order='cartilage',C_order=paste0('grade',seq_len(q)),
      m_vec=rep(3L,q),fold_seed=7L,calibration_seed=8L,test_seed=9L,model_seed=10L,stratification_columns=c('C1','C2'))
    csv<-file.path(root,paste0('q',q,'.csv'));json<-file.path(root,paste0('q',q,'.json'));out<-file.path(root,paste0('q',q))
    write.csv(d[,c('RID','Y','X1','X2','X3','U1',paste0('C',seq_len(q)))],csv,row.names=FALSE)
    jsonlite::write_json(meta,json,auto_unbox=TRUE)
    api$oai_prepare_folds(csv,json,out);x<-api$oai_plan_inputs(out)
    testthat::expect_length(list.files(out),46L);testthat::expect_equal(x$N,n)
    for(f in 1:5) {
      a<-x$cells[[paste(f,'O20',sep='/')]];b<-x$cells[[paste(f,'O50',sep='/')]]
      testthat::expect_true(all(a$train$calibration<=b$train$calibration));testthat::expect_true(all(!a$available|b$available))
      args<-api$oai_plan_fit_args(a,api$oai_settings('pilot',4L))
      testthat::expect_equal(ncol(args$C),q);testthat::expect_true(all(is.na(args$U_obs[a$train$calibration==0,])))
      changed<-a;changed$train$U1[changed$train$calibration==0]<-1e8;changed$test$Y<--1e8;changed$test$U1<-1e8
      testthat::expect_identical(args,api$oai_plan_fit_args(changed,api$oai_settings('pilot',4L)))
    }
  }
})
testthat::test_that('selected-mask metadata is disclosed but never drives mask generation',{
  api<-oai_load_api(repo)
  testthat::expect_identical(api$oai_mask_notice(list(metadata=list(test_U_selection=list(outcome_selected=FALSE)))), '')
  notice<-api$oai_mask_notice(list(metadata=list(test_U_selection=list(outcome_selected=TRUE))))
  testthat::expect_match(notice,'exploratory');testthat::expect_match(notice,'no performance-based mask search')
  testthat::expect_match(api$oai_mask_notice(list(metadata=list(test_U_design='Post-hoc outcome-selected random mask'))),'exploratory')
})
